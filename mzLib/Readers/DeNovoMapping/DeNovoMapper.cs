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
/// With <see cref="DeNovoMappingParameters.HostProteins"/>, each result also separates the hits whose peptide the host
/// produces too (host-ambiguous) from those only the database explains (unambiguous).
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
    /// <see cref="_predictionsBySequence"/> looked up by a slice of a protein, so digestion makes no string per peptide.
    /// </summary>
    private readonly Dictionary<string, int[]>.AlternateLookup<ReadOnlySpan<char>> _sequenceLookup;

    private readonly MappingSequenceDigester _digester;

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
        _sequenceLookup = _predictionsBySequence.GetAlternateLookup<ReadOnlySpan<char>>();
        _digester = new MappingSequenceDigester(_parameters.DigestionParams, _parameters.ReplaceIWithL);

        // After calibration, which keeps every prediction: a host-shared peptide is as good an anchor as any.
        if (_parameters.HostProteins is { } hostProteins)
            MarkHostSharedPredictions(hostProteins);
    }

    public IReadOnlyList<DeNovoPrediction> Predictions { get; }

    /// <summary>
    /// Target proteins of <see cref="DeNovoMappingParameters.HostProteins"/> that were digested, or -1 when no host was given.
    /// </summary>
    public int HostProteinCount { get; private set; } = -1;

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
        var database = RentAccumulator();

        // One thread needs no partitions to merge, and most databases are mapped many at a time, one thread each.
        if (_parameters.MaxThreads <= 1)
        {
            foreach (var protein in proteinList)
            {
                cancellationToken.ThrowIfCancellationRequested();
                database.Add(protein);
            }
        }
        else
        {
            var options = new ParallelOptions { MaxDegreeOfParallelism = _parameters.MaxThreads, CancellationToken = cancellationToken };
            Parallel.ForEach(proteinList, options,
                RentAccumulator,
                (protein, _, local) => { local.Add(protein); return local; },
                local =>
                {
                    lock (database)
                        database.Merge(local);
                    ReturnAccumulator(local);
                });
        }

        // A decoy peptide is only a decoy if no target in this database produces it too. A hit here is always a
        // prediction sequence's, so taking away the target hits takes away every target peptide that matters.
        database.DecoyHits.ExceptWith(database.TargetHits);

        var result = new DeNovoMappingResult
        {
            DatabaseIdentifier = databaseIdentifier,
            DatabaseProteinCount = proteinList.Count(p => !p.IsDecoy),
            DatabasePeptideCount = database.TargetPeptides.Count,
            UniquePeptidesMapped = database.TargetHits.Count,
            UniqueProteinsMapped = database.MatchedProteins.Count,
            HostProteinCount = HostProteinCount,
            UniqueUnambiguousPeptidesMapped = HasHost ? database.TargetHits.Count(hit => !IsHostShared(hit)) : 0,
            UniqueUnambiguousProteinsMapped = database.UnambiguouslyMatchedProteins.Count,
            UniqueDecoyPeptidesMapped = database.DecoyHits.Count,
            UniqueUnambiguousDecoyPeptidesMapped = HasHost ? database.DecoyHits.Count(hit => !IsHostShared(hit)) : -1,
        };

        AddPeptideValues(result, database.TargetHits);

        // Hits in prediction order, so the score and retention-time lists come out as if every prediction were visited.
        var hits = new List<(int Prediction, bool IsTarget)>();
        foreach (var predictionIndices in database.TargetHits)
        foreach (int p in predictionIndices)
            hits.Add((p, true));
        foreach (var predictionIndices in database.DecoyHits)
        foreach (int p in predictionIndices)
            hits.Add((p, false));
        hits.Sort((a, b) => a.Prediction.CompareTo(b.Prediction));

        // Each prediction counts once, as a target or a decoy, however many proteins share its peptide. With a host,
        // it is also unambiguous or ambiguous by the same rule for targets and decoys: whether the host produces its
        // peptide. The totals keep counting every hit.
        foreach (var (i, isTarget) in hits)
        {
            var prediction = Predictions[i];
            bool unambiguous = HasHost && !prediction.IsHostShared;
            if (isTarget)
            {
                result.TargetPredictions++;
                result.TargetPredictionScores.Add(prediction.Score);
                if (unambiguous)
                {
                    result.UnambiguousTargetPredictions++;
                    result.UnambiguousTargetPredictionScores.Add(prediction.Score);
                }

                if (prediction is { ObservedRetentionTime: { } observed, CalibratedRetentionTime: { } calibrated })
                {
                    result.RetentionTimeErrors.Add(observed - calibrated);
                    if (unambiguous)
                        result.UnambiguousRetentionTimeErrors.Add(observed - calibrated);
                }
            }
            else
            {
                result.DecoyPredictions++;
                result.DecoyPredictionScores.Add(prediction.Score);

                if (unambiguous)
                {
                    result.UnambiguousDecoyPredictions++;
                    result.UnambiguousDecoyPredictionScores.Add(prediction.Score);
                }
            }

            result.PredictionCountsByFile[prediction.FileKey] = result.PredictionCountsByFile.GetValueOrDefault(prediction.FileKey) + 1;
        }

        ReturnAccumulator(database);

        if (HasHost)
            result.AmbiguousTargetPredictions = result.TargetPredictions - result.UnambiguousTargetPredictions;
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
    /// Digests the host's target proteins once and flags every prediction whose peptide one of them produces. Only
    /// host peptides that some prediction has are kept, so a whole proteome as host costs one digestion and little memory.
    /// </summary>
    private void MarkHostSharedPredictions(IEnumerable<IBioPolymer> hostProteins)
    {
        var hostTargets = hostProteins.Where(p => !p.IsDecoy).ToList();
        var sharedSequences = new HashSet<string>(StringComparer.Ordinal);
        var options = new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, _parameters.MaxThreads) };

        Parallel.ForEach(hostTargets, options,
            () => new PredictedSequenceCollector(this),
            (protein, _, local) => { _digester.Digest(protein, local); return local; },
            local => { lock (sharedSequences) sharedSequences.UnionWith(local.Sequences); });

        foreach (var sequence in sharedSequences)
        foreach (int p in _predictionsBySequence[sequence])
            Predictions[p].IsHostShared = true;

        HostProteinCount = hostTargets.Count;
    }

    private bool HasHost => HostProteinCount >= 0;

    /// <summary>
    /// Accumulators kept between databases, so mapping thousands of small databases does not grow fresh sets (whose
    /// target-peptide set reaches the large object heap) for each one.
    /// </summary>
    private readonly ConcurrentBag<DatabaseAccumulator> _accumulatorPool = new();

    private DatabaseAccumulator RentAccumulator() => _accumulatorPool.TryTake(out var accumulator) ? accumulator : new DatabaseAccumulator(this);

    private void ReturnAccumulator(DatabaseAccumulator accumulator)
    {
        accumulator.Clear();
        _accumulatorPool.Add(accumulator);
    }

    /// <summary>
    /// One value per target peptide: the median retention-time error and the median score of its predictions, so a
    /// peptide predicted many times counts once. Peptides come in the order of their first prediction.
    /// </summary>
    private void AddPeptideValues(DeNovoMappingResult result, HashSet<int[]> targetHits)
    {
        var errors = new List<double>();
        var scores = new List<double>();
        foreach (var predictionIndices in targetHits.OrderBy(hit => hit[0]))
        {
            errors.Clear();
            scores.Clear();
            foreach (int p in predictionIndices)
            {
                var prediction = Predictions[p];
                scores.Add(prediction.Score);
                if (prediction is { ObservedRetentionTime: { } observed, CalibratedRetentionTime: { } calibrated })
                    errors.Add(observed - calibrated);
            }

            bool unambiguous = HasHost && !IsHostShared(predictionIndices);
            double score = DeNovoMappingResult.Median(scores);
            result.PeptidePredictionScores.Add(score);
            if (unambiguous)
                result.UnambiguousPeptidePredictionScores.Add(score);

            // A peptide none of whose predictions has both retention times has no error to report.
            if (errors.Count == 0)
                continue;
            double error = DeNovoMappingResult.Median(errors);
            result.PeptideRetentionTimeErrors.Add(error);
            if (unambiguous)
                result.UnambiguousPeptideRetentionTimeErrors.Add(error);
        }
    }

    /// <summary>
    /// True when the host produces the sequence of these predictions, which all share it and so its host status.
    /// </summary>
    private bool IsHostShared(int[] predictionIndices) => Predictions[predictionIndices[0]].IsHostShared;

    /// <summary>
    /// Collects the prediction sequences a digestion produces.
    /// </summary>
    private sealed class PredictedSequenceCollector(DeNovoMapper mapper) : IMappingWindowSink
    {
        public HashSet<string> Sequences { get; } = new(StringComparer.Ordinal);

        public void OnWindow(string normalizedSource, int start, int length)
        {
            if (mapper._sequenceLookup.TryGetValue(normalizedSource.AsSpan(start, length), out string? sequence, out _))
                Sequences.Add(sequence);
        }
    }

    /// <summary>
    /// What one database's digestion matched, for the proteins one thread saw; merged across threads. A matched prediction
    /// sequence is kept as its array of prediction indices (one array per sequence, compared by reference), so the work
    /// after digestion grows with the hits, not with the number of predictions.
    /// </summary>
    private sealed class DatabaseAccumulator(DeNovoMapper mapper) : IMappingWindowSink
    {
        private bool _isDecoy;
        private bool _proteinMatched;
        private bool _proteinMatchedUnambiguously;

        /// <summary>Prediction sequences a target protein produces.</summary>
        public HashSet<int[]> TargetHits { get; } = new(ReferenceEqualityComparer.Instance);

        /// <summary>Prediction sequences a decoy protein produces, before the targets' are taken away.</summary>
        public HashSet<int[]> DecoyHits { get; } = new(ReferenceEqualityComparer.Instance);

        /// <summary>Every distinct target peptide, predicted or not, as a window of its normalized protein.</summary>
        public HashSet<SequenceWindow> TargetPeptides { get; } = new(SequenceWindow.Comparer);

        public HashSet<string> MatchedProteins { get; } = new(StringComparer.Ordinal);
        public HashSet<string> UnambiguouslyMatchedProteins { get; } = new(StringComparer.Ordinal);

        public void Add(IBioPolymer protein)
        {
            _isDecoy = protein.IsDecoy;
            _proteinMatched = false;
            _proteinMatchedUnambiguously = false;

            mapper._digester.Digest(protein, this);

            if (_proteinMatched)
                MatchedProteins.Add(protein.Accession);
            if (_proteinMatchedUnambiguously)
                UnambiguouslyMatchedProteins.Add(protein.Accession);
        }

        /// <summary>
        /// Empties the sets for the next database, keeping their capacity.
        /// </summary>
        public void Clear()
        {
            TargetHits.Clear();
            DecoyHits.Clear();
            TargetPeptides.Clear();
            MatchedProteins.Clear();
            UnambiguouslyMatchedProteins.Clear();
        }

        public void OnWindow(string normalizedSource, int start, int length)
        {
            if (!_isDecoy)
                TargetPeptides.Add(new SequenceWindow(normalizedSource, start, length));

            if (!mapper._sequenceLookup.TryGetValue(normalizedSource.AsSpan(start, length), out _, out int[]? predictionIndices))
                return;

            if (_isDecoy)
            {
                DecoyHits.Add(predictionIndices);
                return;
            }

            TargetHits.Add(predictionIndices);
            _proteinMatched = true;
            if (mapper.HasHost && !mapper.IsHostShared(predictionIndices))
                _proteinMatchedUnambiguously = true;
        }

        public void Merge(DatabaseAccumulator other)
        {
            TargetHits.UnionWith(other.TargetHits);
            DecoyHits.UnionWith(other.DecoyHits);
            TargetPeptides.UnionWith(other.TargetPeptides);
            MatchedProteins.UnionWith(other.MatchedProteins);
            UnambiguouslyMatchedProteins.UnionWith(other.UnambiguouslyMatchedProteins);
        }
    }

    /// <summary>
    /// A peptide as a slice of its normalized protein, compared by its residues, so counting distinct peptides needs no
    /// string per peptide.
    /// </summary>
    private readonly record struct SequenceWindow(string Source, int Start, int Length)
    {
        public static readonly IEqualityComparer<SequenceWindow> Comparer = new WindowComparer();

        private sealed class WindowComparer : IEqualityComparer<SequenceWindow>
        {
            public bool Equals(SequenceWindow x, SequenceWindow y) =>
                x.Source.AsSpan(x.Start, x.Length).SequenceEqual(y.Source.AsSpan(y.Start, y.Length));

            public int GetHashCode(SequenceWindow window) => string.GetHashCode(window.Source.AsSpan(window.Start, window.Length));
        }
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
