using System.Collections.Concurrent;
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
/// scores, and the retention-time error of the target hits.
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
    /// assigns observed and predicted retention times.
    /// </summary>
    public DeNovoMapper(IEnumerable<string> mzTabPaths, DeNovoMappingParameters? parameters = null)
        : this(mzTabPaths.Select(p => new CasanovoMzTabFile(p)), parameters)
    {
    }

    public DeNovoMapper(IEnumerable<CasanovoMzTabFile> mzTabFiles, DeNovoMappingParameters? parameters = null)
    {
        _parameters = parameters ?? new DeNovoMappingParameters();
        var retentionTimeMapper = new DeNovoRetentionTimeMapper(_parameters.RetentionTimePredictor);
        var predictions = new List<DeNovoPrediction>();

        foreach (var mzTab in mzTabFiles)
        {
            string fileKey = Path.GetFileNameWithoutExtension(mzTab.FilePath);
            var filePredictions = mzTab.Results
                .Select(r => new DeNovoPrediction(r, fileKey, _parameters.ReplaceIWithL))
                .Where(p => p.MappingSequence.Length > 0)
                .ToList();

            if (_parameters.MapRetentionTimes)
                retentionTimeMapper.MapObservedRetentionTimes(mzTab, filePredictions, _parameters.SpectraFilePaths, _parameters.MaxThreads);

            predictions.AddRange(filePredictions);
        }

        if (_parameters.MapRetentionTimes)
            retentionTimeMapper.PredictRetentionTimes(predictions, _parameters.MaxThreads);

        Predictions = predictions;
        Warnings.AddRange(retentionTimeMapper.Warnings);
        _predictionsBySequence = predictions
            .Select((p, i) => (p.MappingSequence, i))
            .GroupBy(t => t.MappingSequence, StringComparer.Ordinal)
            .ToDictionary(g => g.Key, g => g.Select(t => t.i).ToArray(), StringComparer.Ordinal);
    }

    public IReadOnlyList<DeNovoPrediction> Predictions { get; }

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
    /// (<see cref="IBioPolymer.IsDecoy"/>) count as decoys; supply them to get decoy counts.
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

        Parallel.ForEach(proteinList, options, protein =>
        {
            bool isDecoy = protein.IsDecoy;
            bool proteinMatched = false;

            foreach (var peptide in DigestToMappingSequences(protein))
            {
                if (!isDecoy)
                    databaseTargetPeptides.TryAdd(peptide, 0);

                if (!_predictionsBySequence.TryGetValue(peptide, out var predictionIndices))
                    continue;

                // Several threads may set the same flag; they all write true, so no lock is needed.
                var hits = isDecoy ? decoyHit : targetHit;
                foreach (int i in predictionIndices)
                    hits[i] = true;

                if (!isDecoy)
                {
                    mappedTargetPeptides.TryAdd(peptide, 0);
                    proteinMatched = true;
                }
            }

            if (proteinMatched)
                mappedTargetProteins.TryAdd(protein.Accession, 0);
        });

        var result = new DeNovoMappingResult
        {
            DatabaseIdentifier = databaseIdentifier,
            DatabaseProteinCount = proteinList.Count(p => !p.IsDecoy),
            DatabasePeptideCount = databaseTargetPeptides.Count,
            UniquePeptidesMapped = mappedTargetPeptides.Count,
            UniqueProteinsMapped = mappedTargetProteins.Count,
        };

        // Each prediction counts once as a target and/or once as a decoy, however many proteins share its peptide.
        for (int i = 0; i < Predictions.Count; i++)
        {
            var prediction = Predictions[i];
            if (targetHit[i])
            {
                result.TargetPredictions++;
                result.TargetPredictionScores.Add(prediction.Score);
                if (prediction is { ObservedRetentionTime: not null, PredictedRetentionTime: not null })
                    result.RetentionTimeErrors.Add(prediction.ObservedRetentionTime.Value - prediction.PredictedRetentionTime.Value);
            }

            if (decoyHit[i])
            {
                result.DecoyPredictions++;
                result.DecoyPredictionScores.Add(prediction.Score);
            }

            if (targetHit[i] || decoyHit[i])
                result.PredictionCountsByFile[prediction.FileKey] = result.PredictionCountsByFile.GetValueOrDefault(prediction.FileKey) + (targetHit[i] ? 1 : 0) + (decoyHit[i] ? 1 : 0);
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
