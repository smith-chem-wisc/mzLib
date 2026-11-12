using System.Collections.Concurrent;
using Chromatography.RetentionTimeCalibration;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using MassSpectrometry;
using MzLibUtil;
using Readers.ExternalResults.ResultFiles;

namespace Readers.DeNovoMapping;

/// <summary>
/// Gives each de novo prediction an observed and a predicted retention time.
/// Observed times come from the mzTab when Casanovo wrote them; otherwise the spectra file is read.
/// Predicted times come from Chronologer, in its own units (% acetonitrile), and are mapped onto each file's
/// minutes by a calibration fitted once per file (<see cref="FitCalibrations"/>).
/// </summary>
public class DeNovoRetentionTimeMapper : IDisposable
{
    private readonly RetentionTimePredictor _predictor;
    private readonly bool _ownsPredictor;
    private readonly ConcurrentDictionary<string, double?> _predictionCache = new(StringComparer.Ordinal);
    private readonly Dictionary<string, IrtCalibrationModel> _calibrations = new(StringComparer.Ordinal);

    /// <param name="predictor">
    /// The predictor to use. When null, a Chronologer predictor is created and disposed with this mapper;
    /// a predictor passed in stays the caller's to dispose.
    /// </param>
    public DeNovoRetentionTimeMapper(RetentionTimePredictor? predictor = null)
    {
        _ownsPredictor = predictor is null;
        _predictor = predictor ?? new ChronologerRetentionTimePredictor();
    }

    public List<string> Warnings { get; } = new();

    /// <summary>
    /// The calibration from predicted to observed retention time, per mzTab file (<see cref="DeNovoPrediction.FileKey"/>).
    /// Files that could not be fitted are absent.
    /// </summary>
    public IReadOnlyDictionary<string, IrtCalibrationModel> Calibrations => _calibrations;

    /// <summary>
    /// Fills in the observed retention time of every prediction in <paramref name="predictions"/>, which
    /// must all come from <paramref name="mzTab"/>. The spectra file is only read if any time is missing.
    /// </summary>
    public void MapObservedRetentionTimes(CasanovoMzTabFile mzTab, IReadOnlyList<DeNovoPrediction> predictions,
        IReadOnlyList<string> spectraFilePaths, int maxThreads)
    {
        if (predictions.All(p => p.ObservedRetentionTime.HasValue))
            return;

        // Casanovo normally writes one mzTab per spectra file, so this is usually a single ms_run.
        foreach (var group in predictions.GroupBy(p => SplitSpectraRef(p.Record.SpectraRef).MsRun))
        {
            if (!mzTab.MsRunDictionary.TryGetValue(group.Key, out var recordedPath))
            {
                Warnings.Add($"{Path.GetFileName(mzTab.FilePath)}: no ms_run location for {group.Key}; retention times left empty.");
                continue;
            }

            string? spectraPath = ResolveSpectraFile(recordedPath, spectraFilePaths);
            if (spectraPath is null)
            {
                Warnings.Add($"{Path.GetFileName(mzTab.FilePath)}: spectra file {Path.GetFileName(recordedPath)} not found; retention times left empty.");
                continue;
            }

            MsDataFile dataFile = MsDataFileReader.GetDataFile(spectraPath).LoadAllStaticData(maxThreads: maxThreads);

            // Casanovo only reads MS2 spectra, and "index=N" counts those in file order.
            var ms2Scans = dataFile.GetAllScansList()
                .Where(s => s.MsnOrder > 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();

            foreach (var prediction in group)
            {
                if (prediction.ObservedRetentionTime.HasValue)
                    continue;

                var scan = FindScan(prediction.Record.SpectraRef, ms2Scans, dataFile);
                if (scan is not null)
                    prediction.ObservedRetentionTime = scan.RetentionTime;
            }
        }
    }

    /// <summary>
    /// Predicts a retention time for every prediction. Sequences Chronologer cannot handle are left null.
    /// </summary>
    public void PredictRetentionTimes(IReadOnlyList<DeNovoPrediction> predictions, int maxThreads)
    {
        var toPredict = predictions
            .Where(p => !_predictionCache.ContainsKey(p.FullSequence))
            .GroupBy(p => p.FullSequence, StringComparer.Ordinal)
            .Select(g => (IRetentionPredictable)g.First())
            .ToList();

        foreach (var (value, peptide, _) in _predictor.PredictRetentionTimeEquivalents(toPredict, maxThreads))
            _predictionCache[peptide.FullSequence] = value;

        foreach (var prediction in predictions)
            prediction.PredictedRetentionTime = _predictionCache.GetValueOrDefault(prediction.FullSequence);
    }

    /// <summary>
    /// Fits, for each mzTab file, the predicted retention time against the observed one with mzLib's robust LOWESS
    /// (<see cref="IrtCalibration.Fit"/>), then sets <see cref="DeNovoPrediction.CalibratedRetentionTime"/> on every
    /// prediction of that file. Anchors are the file's predictions scoring at least <paramref name="minimumAnchorScore"/>,
    /// one per peptide (the best-scoring one). They do not depend on any database, so every database is scored against
    /// the same calibration and the errors are comparable between databases. A file that cannot be fitted, for example
    /// with fewer than 20 anchors, is warned about and its predictions get no calibrated time.
    /// </summary>
    public void FitCalibrations(IReadOnlyList<DeNovoPrediction> predictions, double minimumAnchorScore)
    {
        foreach (var file in predictions.GroupBy(p => p.FileKey, StringComparer.Ordinal))
        {
            var anchors = file
                .Where(p => p.Score >= minimumAnchorScore
                            && p.ObservedRetentionTime is { } observed && double.IsFinite(observed)
                            && p.PredictedRetentionTime is { } predicted && double.IsFinite(predicted))
                .GroupBy(p => p.MappingSequence, StringComparer.Ordinal)
                .Select(g => g.MaxBy(p => p.Score)!)
                .Select(p => (Rt: new RtMinutes(p.ObservedRetentionTime!.Value), Irt: new Irt(p.PredictedRetentionTime!.Value)))
                .ToList();

            IrtCalibrationModel calibration;
            try
            {
                calibration = IrtCalibration.Fit(anchors);
            }
            catch (ArgumentException e)
            {
                Warnings.Add($"{file.Key}: predicted retention times could not be calibrated to this file's retention times ({e.Message}); " +
                             "no retention time errors are reported for it.");
                continue;
            }

            _calibrations[file.Key] = calibration;
            foreach (var prediction in file)
            {
                if (prediction.PredictedRetentionTime is { } predicted)
                    prediction.CalibratedRetentionTime = calibration.ToRtMinutes(new Irt(predicted)).Value;
            }
        }
    }

    /// <summary>
    /// The factor that turns an mzTab's PSM retention_time into minutes. mzTab 1.0 reports it in seconds unless
    /// a colunit-psm entry gives another unit; only minutes (UO:0000031) is recognized as such.
    /// </summary>
    public static double MzTabRetentionTimeToMinutes(CasanovoMzTabFile mzTab)
    {
        foreach (var line in mzTab.HeaderLines)
        {
            var fields = line.Split('	');
            if (fields.Length < 3 || fields[1] != "colunit-psm"
                || !fields[2].StartsWith("retention_time=", StringComparison.OrdinalIgnoreCase))
                continue;

            string unit = fields[2];
            if (unit.Contains("UO:0000031", StringComparison.OrdinalIgnoreCase)
                || unit.Contains("minute", StringComparison.OrdinalIgnoreCase))
                return 1;
        }

        return 1 / 60.0;
    }

    public void Dispose()
    {
        if (_ownsPredictor)
            _predictor.Dispose();
        GC.SuppressFinalize(this);
    }

    /// <summary>
    /// Picks the spectra file for an mzTab ms_run: a supplied file with the same name (without extension),
    /// otherwise the recorded location if it still exists.
    /// </summary>
    public static string? ResolveSpectraFile(string recordedPath, IReadOnlyList<string> spectraFilePaths)
    {
        string name = Path.GetFileNameWithoutExtension(recordedPath);
        var match = spectraFilePaths.FirstOrDefault(p =>
            string.Equals(Path.GetFileNameWithoutExtension(p), name, StringComparison.OrdinalIgnoreCase));
        if (match is not null)
            return match;

        return File.Exists(recordedPath) ? recordedPath : null;
    }

    /// <summary>
    /// Resolves a Casanovo spectra_ref. "index=N" is the zero-based position among the MS2 spectra in the
    /// file (for an MGF this is not the SCANS number); "scan=N" is a scan number.
    /// </summary>
    internal static MsDataScan? FindScan(string spectraRef, MsDataScan[] ms2ScansInFileOrder, MsDataFile dataFile)
    {
        string reference = SplitSpectraRef(spectraRef).Reference;

        if (reference.StartsWith("index=") && int.TryParse(reference.AsSpan(6), out int index))
            return index >= 0 && index < ms2ScansInFileOrder.Length ? ms2ScansInFileOrder[index] : null;

        foreach (var part in reference.Split(' '))
        {
            if (part.StartsWith("scan=") && int.TryParse(part.AsSpan(5), out int scanNumber))
            {
                try
                {
                    return dataFile.GetOneBasedScan(scanNumber);
                }
                catch (Exception e) when (e is IndexOutOfRangeException or ArgumentOutOfRangeException)
                {
                    return null;
                }
            }
        }

        return null;
    }

    private static (string MsRun, string Reference) SplitSpectraRef(string spectraRef)
    {
        int colon = spectraRef.IndexOf(':');
        return colon < 0 ? (string.Empty, spectraRef) : (spectraRef[..colon], spectraRef[(colon + 1)..]);
    }
}
