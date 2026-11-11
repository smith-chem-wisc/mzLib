using System.Collections.Concurrent;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using MassSpectrometry;
using Readers.ExternalResults.ResultFiles;

namespace Readers.DeNovoMapping;

/// <summary>
/// Gives each de novo prediction an observed and a predicted retention time.
/// Observed times come from the mzTab when Casanovo wrote them; otherwise the spectra file is read.
/// Predicted times come from Chronologer.
/// </summary>
public class DeNovoRetentionTimeMapper
{
    private readonly RetentionTimePredictor _predictor;
    private readonly ConcurrentDictionary<string, double?> _predictionCache = new(StringComparer.Ordinal);

    public DeNovoRetentionTimeMapper(RetentionTimePredictor? predictor = null)
    {
        _predictor = predictor ?? new ChronologerRetentionTimePredictor();
    }

    public List<string> Warnings { get; } = new();

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
