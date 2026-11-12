using System.Collections.Concurrent;
using Chromatography.RetentionTimeCalibration;
using Omics;
using Omics.Modifications;
using Proteomics;
using Readers.ExternalResults.ResultFiles;
using UsefulProteomicsDatabases;

namespace Readers.DeNovoMapping;

/// <summary>
/// Maps Casanovo de novo predictions onto protein databases. A prediction maps to a database when its
/// peptide (I/L-normalized, modifications ignored) is produced by digesting a protein in it. Each database
/// gets one <see cref="DeNovoMappingResult"/>: how many predictions hit its targets and its decoys, their
/// scores, and the retention-time error of the target hits in minutes. A decoy peptide that equals a target
/// peptide of the same database is not a decoy, as in <c>MslLibraryBuilder</c>.
/// </summary>
/// <example>
/// <code>
/// var mapper = new DeNovoMapper(mzTabPaths, parameters);
/// var results = mapper.MapDatabases(fastaPaths);
/// DeNovoMapper.WriteSummary(results, "MappingSummary.tsv");
/// </code>
/// </example>
public class DeNovoMapper
{
    private readonly DeNovoMappingParameters _parameters;
    private readonly Dictionary<string, int[]> _predictionsBySequence;

    /// <summary>
    /// Reads the mzTab files and, when <see cref="DeNovoMappingParameters.MapRetentionTimes"/> is set,
    /// assigns observed and predicted retention times and calibrates the predictions to each file's minutes.
    /// </summary>
    public DeNovoMapper(IEnumerable<string> mzTabPaths, DeNovoMappingParameters? parameters = null)
        : this(mzTabPaths.Select(p => new CasanovoMzTabFile(p)), parameters)
    {
    }

    public DeNovoMapper(IEnumerable<CasanovoMzTabFile> mzTabFiles, DeNovoMappingParameters? parameters = null)
    {
        _parameters = parameters ?? new DeNovoMappingParameters();
        var mzTabList = mzTabFiles.ToList();
        var predictions = new List<DeNovoPrediction>();
        var predictionsByMzTab = new List<(CasanovoMzTabFile MzTab, List<DeNovoPrediction> Predictions)>();

        foreach (var mzTab in mzTabList)
        {
            string fileKey = Path.GetFileNameWithoutExtension(mzTab.FilePath);
            double toMinutes = DeNovoRetentionTimeMapper.MzTabRetentionTimeToMinutes(mzTab);
            var filePredictions = mzTab.Results
                .Select(r => new DeNovoPrediction(r, fileKey, _parameters.ReplaceIWithL, toMinutes))
                .Where(p => p.MappingSequence.Length > 0)
                .ToList();

            predictionsByMzTab.Add((mzTab, filePredictions));
            predictions.AddRange(filePredictions);
        }

        Predictions = predictions;

        // All retention-time work happens here, so the mapper (and a predictor it creates) lives only this long.
        if (_parameters.MapRetentionTimes)
        {
            using var retentionTimeMapper = new DeNovoRetentionTimeMapper(_parameters.RetentionTimePredictor);
            foreach (var (mzTab, filePredictions) in predictionsByMzTab)
                retentionTimeMapper.MapObservedRetentionTimes(mzTab, filePredictions, _parameters.SpectraFilePaths, _parameters.MaxThreads);

            retentionTimeMapper.PredictRetentionTimes(predictions, _parameters.MaxThreads);
            retentionTimeMapper.FitCalibrations(predictions, _parameters.CalibrationScoreThreshold);

            RetentionTimeCalibrations = retentionTimeMapper.Calibrations.ToDictionary(p => p.Key, p => p.Value, StringComparer.Ordinal);
            Warnings.AddRange(retentionTimeMapper.Warnings);
        }

        _predictionsBySequence = predictions
            .Select((p, i) => (p.MappingSequence, i))
            .GroupBy(t => t.MappingSequence, StringComparer.Ordinal)
            .ToDictionary(g => g.Key, g => g.Select(t => t.i).ToArray(), StringComparer.Ordinal);
    }

    public IReadOnlyList<DeNovoPrediction> Predictions { get; }

    /// <summary>
    /// Each mzTab file's calibration from predicted retention time to minutes, keyed by
    /// <see cref="DeNovoPrediction.FileKey"/>. Empty when retention times are not mapped; files that could not be
    /// fitted are absent.
    /// </summary>
    public IReadOnlyDictionary<string, IrtCalibrationModel> RetentionTimeCalibrations { get; } =
        new Dictionary<string, IrtCalibrationModel>();

    /// <summary>
    /// Problems that did not stop the mapping, such as a spectra file that could not be found.
    /// </summary>
    public List<string> Warnings { get; } = new();

    /// <summary>
    /// Reads each database (FASTA or UniProt XML), adds decoys of <see cref="DeNovoMappingParameters.DecoyType"/>
    /// and maps the predictions onto it. The database identifier is the file name without extension.
    /// </summary>
    public List<DeNovoMappingResult> MapDatabases(IEnumerable<string> databasePaths, CancellationToken cancellationToken = default)
    {
        var results = new List<DeNovoMappingResult>();
        foreach (var databasePath in databasePaths)
        {
            cancellationToken.ThrowIfCancellationRequested();
            results.Add(MapDatabase(Path.GetFileNameWithoutExtension(databasePath), LoadProteins(databasePath), cancellationToken));
        }
        return results;
    }

    /// <summary>
    /// Maps the predictions onto proteins that are already loaded. Decoy proteins
    /// (<see cref="IBioPolymer.IsDecoy"/>) count as decoys; supply them to get decoy counts. Decoy peptides
    /// that equal one of these targets' peptides are dropped, so a prediction is never both target and decoy.
    /// </summary>
    public DeNovoMappingResult MapDatabase(string databaseIdentifier, IEnumerable<IBioPolymer> proteins,
        CancellationToken cancellationToken = default)
    {
        var proteinList = proteins as IReadOnlyList<IBioPolymer> ?? proteins.ToList();
        var targetHit = new bool[Predictions.Count];
        var decoyHit = new bool[Predictions.Count];
        var databaseTargetPeptides = new ConcurrentDictionary<string, byte>(StringComparer.Ordinal);
        var mappedTargetPeptides = new ConcurrentDictionary<string, byte>(StringComparer.Ordinal);
        var mappedTargetProteins = new ConcurrentDictionary<string, byte>(StringComparer.Ordinal);

        var options = new ParallelOptions
        {
            MaxDegreeOfParallelism = Math.Max(1, _parameters.MaxThreads),
            CancellationToken = cancellationToken,
        };

        // Digest first: a decoy peptide is only a decoy if no target in this database produces it too.
        var digests = new HashSet<string>[proteinList.Count];
        Parallel.For(0, proteinList.Count, options, i =>
        {
            digests[i] = DigestToMappingSequences(proteinList[i]);
            if (!proteinList[i].IsDecoy)
            {
                foreach (var peptide in digests[i])
                    databaseTargetPeptides.TryAdd(peptide, 0);
            }
        });

        Parallel.For(0, proteinList.Count, options, i =>
        {
            bool isDecoy = proteinList[i].IsDecoy;
            bool proteinMatched = false;

            foreach (var peptide in digests[i])
            {
                if (isDecoy && databaseTargetPeptides.ContainsKey(peptide))
                    continue;

                if (!_predictionsBySequence.TryGetValue(peptide, out var predictionIndices))
                    continue;

                // Several threads may set the same flag; they all write true, so no lock is needed.
                var hits = isDecoy ? decoyHit : targetHit;
                foreach (int p in predictionIndices)
                    hits[p] = true;

                if (!isDecoy)
                {
                    mappedTargetPeptides.TryAdd(peptide, 0);
                    proteinMatched = true;
                }
            }

            if (proteinMatched)
                mappedTargetProteins.TryAdd(proteinList[i].Accession, 0);
        });

        var result = new DeNovoMappingResult
        {
            DatabaseIdentifier = databaseIdentifier,
            DatabaseProteinCount = proteinList.Count(p => !p.IsDecoy),
            DatabasePeptideCount = databaseTargetPeptides.Count,
            UniquePeptidesMapped = mappedTargetPeptides.Count,
            UniqueProteinsMapped = mappedTargetProteins.Count,
        };

        // Each prediction counts once, as a target or a decoy, however many proteins share its peptide.
        for (int i = 0; i < Predictions.Count; i++)
        {
            var prediction = Predictions[i];
            if (targetHit[i])
            {
                result.TargetPredictions++;
                result.TargetPredictionScores.Add(prediction.Score);
                if (prediction is { ObservedRetentionTime: { } observed, CalibratedRetentionTime: { } calibrated })
                    result.RetentionTimeErrors.Add(observed - calibrated);
            }

            if (decoyHit[i])
            {
                result.DecoyPredictions++;
                result.DecoyPredictionScores.Add(prediction.Score);
            }

            if (targetHit[i] || decoyHit[i])
                result.PredictionCountsByFile[prediction.FileKey] = result.PredictionCountsByFile.GetValueOrDefault(prediction.FileKey) + 1;
        }

        result.FinalizeValues();
        return result;
    }

    /// <summary>
    /// Writes the mapping summary as a tab-separated file.
    /// </summary>
    public static void WriteSummary(IEnumerable<DeNovoMappingResult> results, string outputPath)
    {
        new DeNovoMappingResultFile { Results = results.ToList() }.WriteResults(outputPath);
    }

    /// <summary>
    /// Distinct unmodified peptides of the protein, normalized the same way as the predictions.
    /// </summary>
    private HashSet<string> DigestToMappingSequences(IBioPolymer protein)
    {
        var peptides = new HashSet<string>(StringComparer.Ordinal);
        foreach (var peptide in protein.Digest(_parameters.DigestionParams, new List<Modification>(), new List<Modification>()))
            peptides.Add(DeNovoPrediction.NormalizeSequence(peptide.BaseSequence, _parameters.ReplaceIWithL));
        return peptides;
    }

    private List<IBioPolymer> LoadProteins(string databasePath)
    {
        bool isXml = databasePath.EndsWith(".xml", StringComparison.OrdinalIgnoreCase)
                     || databasePath.EndsWith(".xml.gz", StringComparison.OrdinalIgnoreCase);

        List<Protein> proteins = isXml
            ? ProteinDbLoader.LoadProteinXML(databasePath, true, _parameters.DecoyType, new List<Modification>(), false,
                new List<string>(), out _, _parameters.MaxThreads)
            : ProteinDbLoader.LoadProteinFasta(databasePath, true, _parameters.DecoyType, false, out _,
                maxThreads: _parameters.MaxThreads);

        return proteins.Cast<IBioPolymer>().ToList();
    }
}
