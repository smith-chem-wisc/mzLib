using System;
using System.Collections.Generic;
using System.Linq;

namespace StatisticalModels
{
    /// <summary>How a target-decoy rescoring ended.</summary>
    public enum RescoreStatus
    {
        /// <summary>Every fold trained a discriminant.</summary>
        Rescored,

        /// <summary>There were no decoys, so nothing could be scored; every score is NaN.</summary>
        NoDecoys,

        /// <summary>There were no targets, so nothing could be scored; every score is NaN.</summary>
        NoTargets,

        /// <summary>
        /// At least one fold's training rows lacked positives or decoys to train on. That fold was scored by the best single
        /// feature instead.
        /// </summary>
        FoldStarved,
    }

    /// <summary>The model the rescorer fits within each fold.</summary>
    public enum RescoreModel
    {
        /// <summary>A ridge Fisher linear discriminant (<see cref="LinearDiscriminant"/>), refit each iteration.</summary>
        LinearDiscriminant,

        /// <summary>
        /// The linear iterations choose the training rows, then an averaged ensemble of small tanh networks
        /// (<see cref="MultilayerPerceptron"/>) is trained on them, all targets against all decoys (as DIA-NN does), and scores
        /// by the logit of its mean probability.
        /// </summary>
        NeuralNetworkEnsemble,
    }

    /// <summary>The combined score of each candidate, the fold that scored it, and how the rescoring ended.</summary>
    public sealed class RescoreResult
    {
        internal RescoreResult(double[] scores, int[] folds, RescoreStatus status)
        {
            Scores = scores;
            Folds = folds;
            Status = status;
        }

        /// <summary>
        /// One score per candidate, in input order; higher is better. Each fold's scores are normalized on its own training
        /// rows (0 at that fold's q-value cutoff, −1 at its median decoy), so scores from different folds are comparable.
        /// </summary>
        public double[] Scores { get; }

        /// <summary>The fold each candidate belonged to; that fold's model never trained on it.</summary>
        public int[] Folds { get; }

        public RescoreStatus Status { get; }
    }

    /// <summary>
    /// Semi-supervised target-decoy rescoring in the style of mProphet and Percolator: combine many features into one
    /// score for target-decoy q-values, trained only on the search's own targets and decoys.
    /// <para>
    /// Candidates are split into folds by group (for example one peptide's charge states), and a group never straddles
    /// folds. Fold assignment uses only the group keys, never the labels. Each fold is scored by a model trained on the
    /// other folds, so no candidate is scored by a model that saw it. Within the training rows, the first iteration ranks by
    /// the best single feature (either sign). Each later iteration fits a <see cref="LinearDiscriminant"/> of targets
    /// passing the q-value cutoff against all decoys, then re-ranks.
    /// </para>
    /// </summary>
    /// <summary>Which rows train the network when a fold has more than the cap.</summary>
    public enum NetworkTrainingSample
    {
        /// <summary>A random sample of the fold's training rows.</summary>
        Random,

        /// <summary>Half from the targets the linear model ranks highest, half from the decoys it ranks highest.</summary>
        Confident,

        /// <summary>
        /// The rows the linear model ranks highest, targets and decoys together: one score threshold, so no score band
        /// trains on decoys alone.
        /// </summary>
        ConfidentPooled,
    }

    public static class TargetDecoyRescorer
    {
        /// <param name="networkPositiveQValue">
        /// Network model, confident sample only: positives are the targets passing this q-value in the linear ranking (with as
        /// many top decoys) rather than the top half-cap of targets whatever their q. Null keeps the latter.
        /// </param>
        /// <param name="networkTargetFraction">
        /// Network model, confident sample only: the share of the cap taken from the top targets, the rest from the top decoys.
        /// One half by default; strictly between 0 and 1.
        /// </param>
        /// <param name="networkLayers">
        /// Network model only: units in each hidden layer, input side first. DIA-NN 2020's 25-20-15-10-5 when null.
        /// </param>
        /// <param name="networkPasses">
        /// Network model only: training passes. The first trains on each candidate group's top row by the linear model;
        /// each later pass re-picks the top rows with the previous network and trains a new one, as DIA-NN trains twice.
        /// 1 by default.
        /// </param>
        /// <param name="networkTrainingSample">
        /// Network model only: which rows train it when a fold has more than <paramref name="maxNetworkTrainingRows"/>. Random by
        /// default; Confident takes the highest-ranked targets and decoys, half each, as DIA-NN trains after removing
        /// low-confidence identifications.
        /// </param>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">Lengths disagree, rows are ragged, or a value is not finite.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A count or the q-value cutoff is out of range.</exception>
        public static RescoreResult Score(IReadOnlyList<double[]> features, IReadOnlyList<bool> isDecoy, IReadOnlyList<string> groupKeys,
            int folds = 3, int iterations = 3, double positiveQValue = 0.01, IReadOnlyList<int>? candidateGroups = null,
            RescoreModel model = RescoreModel.LinearDiscriminant, int? maxNetworkTrainingRows = null, int randomSeed = 0,
            int networkMembers = NetworkMembers, int networkEpochs = NetworkEpochs, int networkPasses = 1,
            NetworkTrainingSample networkTrainingSample = NetworkTrainingSample.Random, int? normalizationGroups = null,
            IReadOnlyList<int>? networkLayers = null, double? networkPositiveQValue = null, double networkTargetFraction = 0.5)
            => ScoreRows(features, isDecoy, groupKeys, folds, iterations, positiveQValue, candidateGroups, model, maxNetworkTrainingRows, randomSeed,
                networkMembers, networkEpochs, networkPasses, networkTrainingSample, normalizationGroups, networkLayers, networkPositiveQValue, networkTargetFraction);

        /// <summary>
        /// <see cref="Score(IReadOnlyList{double[]}, IReadOnlyList{bool}, IReadOnlyList{string}, int, int, double, IReadOnlyList{int}?, RescoreModel, int?, int, int, int, int, NetworkTrainingSample, int?, IReadOnlyList{int}?, double?, double)"/>
        /// on float features. Every value is widened exactly to a double before any arithmetic, so the scores are bit for bit
        /// those of the same values in double rows, at half the memory: a DIA search keeps tens of millions of candidate rows.
        /// A name of its own, not an overload, so that a call passing <c>null</c> or <c>[]</c> stays unambiguous.
        /// </summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">Lengths disagree, rows are ragged, or a value is not finite.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A count or the q-value cutoff is out of range.</exception>
        public static RescoreResult ScoreSinglePrecision(IReadOnlyList<float[]> features, IReadOnlyList<bool> isDecoy, IReadOnlyList<string> groupKeys,
            int folds = 3, int iterations = 3, double positiveQValue = 0.01, IReadOnlyList<int>? candidateGroups = null,
            RescoreModel model = RescoreModel.LinearDiscriminant, int? maxNetworkTrainingRows = null, int randomSeed = 0,
            int networkMembers = NetworkMembers, int networkEpochs = NetworkEpochs, int networkPasses = 1,
            NetworkTrainingSample networkTrainingSample = NetworkTrainingSample.Random, int? normalizationGroups = null,
            IReadOnlyList<int>? networkLayers = null, double? networkPositiveQValue = null, double networkTargetFraction = 0.5)
            => ScoreRows(features, isDecoy, groupKeys, folds, iterations, positiveQValue, candidateGroups, model, maxNetworkTrainingRows, randomSeed,
                networkMembers, networkEpochs, networkPasses, networkTrainingSample, normalizationGroups, networkLayers, networkPositiveQValue, networkTargetFraction);

        private static RescoreResult ScoreRows<T>(IReadOnlyList<T[]> features, IReadOnlyList<bool> isDecoy, IReadOnlyList<string> groupKeys,
            int folds, int iterations, double positiveQValue, IReadOnlyList<int>? candidateGroups,
            RescoreModel model, int? maxNetworkTrainingRows, int randomSeed,
            int networkMembers, int networkEpochs, int networkPasses,
            NetworkTrainingSample networkTrainingSample, int? normalizationGroups,
            IReadOnlyList<int>? networkLayers, double? networkPositiveQValue, double networkTargetFraction)
            where T : unmanaged, System.Numerics.IBinaryFloatingPointIeee754<T>
        {
            if (!(networkTargetFraction > 0 && networkTargetFraction < 1))
                throw new ArgumentOutOfRangeException(nameof(networkTargetFraction), "The target share must lie strictly between 0 and 1.");
            ArgumentNullException.ThrowIfNull(features);
            ArgumentNullException.ThrowIfNull(isDecoy);
            ArgumentNullException.ThrowIfNull(groupKeys);
            if (isDecoy.Count != features.Count)
                throw new ArgumentException($"There are {features.Count} candidates but {isDecoy.Count} decoy labels.", nameof(isDecoy));
            if (groupKeys.Count != features.Count)
                throw new ArgumentException($"There are {features.Count} candidates but {groupKeys.Count} group keys.", nameof(groupKeys));
            if (candidateGroups is not null && candidateGroups.Count != features.Count)
                throw new ArgumentException($"There are {features.Count} candidates but {candidateGroups.Count} candidate groups.", nameof(candidateGroups));
            if (folds < 2)
                throw new ArgumentOutOfRangeException(nameof(folds), folds, "There must be at least 2 folds.");
            if (iterations < 1)
                throw new ArgumentOutOfRangeException(nameof(iterations), iterations, "There must be at least 1 iteration.");
            if (networkMembers < 1)
                throw new ArgumentOutOfRangeException(nameof(networkMembers), networkMembers, "The ensemble needs at least one network.");
            if (networkEpochs < 1)
                throw new ArgumentOutOfRangeException(nameof(networkEpochs), networkEpochs, "Training needs at least one epoch.");
            if (networkPasses < 1)
                throw new ArgumentOutOfRangeException(nameof(networkPasses), networkPasses, "The network needs at least one training pass.");
            if (networkLayers is not null && (networkLayers.Count == 0 || networkLayers.Any(units => units < 1)))
                throw new ArgumentOutOfRangeException(nameof(networkLayers), "The network needs at least one hidden layer, each of at least one unit.");
            if (normalizationGroups is < 1)
                throw new ArgumentOutOfRangeException(nameof(normalizationGroups), normalizationGroups, "Normalisation needs at least one group.");
            if (maxNetworkTrainingRows is < 2)
                throw new ArgumentOutOfRangeException(nameof(maxNetworkTrainingRows), maxNetworkTrainingRows, "The network needs at least 2 training rows.");
            if (!(positiveQValue > 0 && positiveQValue < 1))
                throw new ArgumentOutOfRangeException(nameof(positiveQValue), positiveQValue, "The q-value cutoff must be in (0, 1).");
            int n = features.Count;
            int p = n == 0 ? 0 : features[0].Length;
            foreach (var row in features)
                if (row is null || row.Length != p || !FeatureRows.AllFinite(row))
                    throw new ArgumentException("Every candidate needs the same number of finite features.", nameof(features));

            int[] fold = AssignFolds(groupKeys, folds);
            var scores = Enumerable.Repeat(double.NaN, n).ToArray();
            if (!isDecoy.Any(d => d))
                return new RescoreResult(scores, fold, RescoreStatus.NoDecoys);
            if (isDecoy.All(d => d))
                return new RescoreResult(scores, fold, RescoreStatus.NoTargets);

            // Folds are independent (each writes only its own held-out rows, with its own seed), so they run in parallel
            var starved = new bool[folds];
            System.Threading.Tasks.Parallel.For(0, folds, f =>
            {
                int[] train = Enumerable.Range(0, n).Where(i => fold[i] != f).ToArray();
                int[] test = Enumerable.Range(0, n).Where(i => fold[i] == f).ToArray();
                if (test.Length == 0)
                    return;

                var trainGroups = RowGroups.Of(train, candidateGroups);
                // The folds run at once, so each seed step takes its share of the cores (and of the workers' buffers)
                var seed = BestSingleFeature(features, isDecoy, trainGroups, positiveQValue,
                    maxParallelism: Math.Max(1, Environment.ProcessorCount / Math.Max(1, folds)));
                Func<T[], double> scorer = x => seed.Sign * double.CreateTruncating(x[seed.Feature]);
                MultilayerPerceptronEnsemble? network = null; // the scorer, when it is the network
                bool trained = false;
                for (int iteration = 0; iteration < iterations; iteration++)
                {
                    // Only each candidate group's top row under the current model trains it (pyProphet)
                    int[] active = TopPerGroup(trainGroups, i => scorer(features[i]));
                    double[] trainScores = active.Select(i => scorer(features[i])).ToArray();
                    double[] q = QValues(trainScores, active.Select(i => isDecoy[i]).ToArray());
                    double cutoff = TrainingCutoff(q, active.Select(i => isDecoy[i]).ToArray(), positiveQValue);
                    var rows = new List<T[]>();
                    var positive = new List<bool>();
                    for (int t = 0; t < active.Length; t++)
                    {
                        bool decoy = isDecoy[active[t]];
                        if (decoy || q[t] <= cutoff)
                        {
                            rows.Add(features[active[t]]);
                            positive.Add(!decoy);
                        }
                    }
                    if (positive.Count(b => b) < 2 || positive.All(b => b))
                        break;
                    var fit = LinearDiscriminant.FitRows(rows, positive);
                    scorer = x => fit.ScoreOf<T>(x);
                    trained = true;
                }
                if (!trained)
                    starved[f] = true;

                if (model == RescoreModel.NeuralNetworkEnsemble && trained)
                {
                    for (int pass = 0; pass < networkPasses; pass++)
                    {
                        // Each group's top row trains the network: by the line on the first pass, and on each later pass by
                        // the previous network, which can tell the real candidate where the line cannot (DIA-NN trains twice)
                        Func<int, double> pick;
                        if (pass == 0)
                            pick = i => scorer(features[i]);
                        else
                        {
                            var previous = scorer;
                            var current = new double[n];
                            System.Threading.Tasks.Parallel.ForEach(train, i => current[i] = previous(features[i]));
                            pick = i => current[i];
                        }
                        int[] rows = TopPerGroup(trainGroups, pick);
                        if (maxNetworkTrainingRows is int pooledCap && rows.Length > pooledCap && networkTrainingSample == NetworkTrainingSample.ConfidentPooled)
                        {
                            // One threshold for both labels: equal halves leave a band where only decoys train
                            rows = PooledTrainingRows(rows.OrderByDescending(pick).ThenBy(i => i).ToArray(), isDecoy, pooledCap);
                        }
                        else if (maxNetworkTrainingRows is int cap && rows.Length > cap && networkTrainingSample == NetworkTrainingSample.Confident)
                        {
                            // As DIA-NN removes low-confidence identifications before training: the targets and the decoys
                            // ranked highest, half the cap each, so real targets are a large share of the positives
                            var ranked = rows.OrderByDescending(pick).ThenBy(i => i).ToArray();
                            rows = ConfidentTrainingRows(ranked, ranked.Select(pick).ToArray(), isDecoy, cap, networkPositiveQValue, networkTargetFraction);
                        }
                        else if (maxNetworkTrainingRows is int cap2 && rows.Length > cap2)
                        {
                            // A random subsample of the fold's own training rows, seeded by the fold. Not the top rows by the
                            // linear score: where the line misses the signal, its top rows are the wrong ones.
                            rows = RandomTrainingRows(rows, isDecoy, cap2, new Random(31 + f + 1000 * randomSeed + 100_000 * pass));
                        }
                        // The training rows as doubles (a float row widened exactly; a double row as it is): a capped sample
                        var ensemble = MultilayerPerceptron.TrainEnsemble(rows.Select(i => FeatureRows.Widened(features[i])).ToList(), rows.Select(i => !isDecoy[i]).ToList(),
                            networkMembers, networkLayers ?? NetworkLayers, networkEpochs, seed: 17 + f + 1000 * randomSeed + 100_000 * pass);
                        scorer = x => ensemble.PredictLogitOf(x);
                        network = ensemble;
                    }
                }

                // Every row scored once, in parallel: the folds alone use three cores, and a network ensemble over millions
                // of rows otherwise costs more than training it
                var final = new double[n];
                var normalizers = trainGroups;
                if (normalizationGroups is int sampleSize && trainGroups.Count > sampleSize)
                {
                    // A random sample of training groups, seeded by the fold, sets the normalisation: scoring every training
                    // row only for that was twice the work of scoring the fold's own rows
                    int[] sampledRows = SampleGroupRows(trainGroups, sampleSize, new Random(53 + f + 1000 * randomSeed));
                    normalizers = RowGroups.Of(sampledRows, candidateGroups);
                    int[] needed = sampledRows.Concat(test).ToArray();
                    ScoreRows(needed);
                }
                else
                    ScoreRows(null);

                // The network scores a block of rows at a time (bit-identical to one at a time, several times cheaper)
                void ScoreRows(int[]? which)
                {
                    int count = which?.Length ?? n;
                    if (network is null)
                    {
                        System.Threading.Tasks.Parallel.For(0, count, k => { int i = which?[k] ?? k; final[i] = scorer(features[i]); });
                        return;
                    }
                    System.Threading.Tasks.Parallel.ForEach(System.Collections.Concurrent.Partitioner.Create(0, count, 4096), range =>
                    {
                        var block = new int[range.Item2 - range.Item1];
                        for (int k = 0; k < block.Length; k++)
                            block[k] = which?[range.Item1 + k] ?? range.Item1 + k;
                        network.PredictLogitsOf(features, block, final);
                    });
                }

                // Normalize on the training rows so that folds are comparable when pooled
                int[] finalActive = TopPerGroup(normalizers, i => final[i]);
                double[] finalTrain = finalActive.Select(i => final[i]).ToArray();
                bool[] trainDecoy = finalActive.Select(i => isDecoy[i]).ToArray();
                double[] finalQ = QValues(finalTrain, trainDecoy);
                double[] decoyScores = finalTrain.Where((_, t) => trainDecoy[t]).Order().ToArray();
                double medianDecoy = decoyScores.Length == 0 ? 0 : decoyScores[decoyScores.Length / 2];
                double[] passing = finalTrain.Where((_, t) => !trainDecoy[t] && finalQ[t] <= positiveQValue).ToArray();
                double threshold = passing.Length > 0 ? passing.Min() : (decoyScores.Length > 0 ? decoyScores[^1] : medianDecoy + 1);
                double scale = threshold > medianDecoy ? threshold - medianDecoy : 1;
                foreach (int i in test)
                    scores[i] = (final[i] - threshold) / scale;
            });
            return new RescoreResult(scores, fold, starved.Any(s => s) ? RescoreStatus.FoldStarved : RescoreStatus.Rescored);
        }

        /// <summary>
        /// Groups sorted by key (ordinal) and dealt round-robin, so a group stays whole and the assignment depends only on the
        /// keys, never on the labels or the row order.
        /// </summary>
        internal static int[] AssignFolds(IReadOnlyList<string> groupKeys, int folds)
        {
            var foldOfGroup = groupKeys.Distinct().Order(StringComparer.Ordinal)
                .Select((key, index) => (key, index))
                .ToDictionary(g => g.key, g => g.index % folds, StringComparer.Ordinal);
            return groupKeys.Select(key => foldOfGroup[key]).ToArray();
        }

        public const int NetworkMembers = 5;
        public const int NetworkEpochs = 10;
        private static readonly int[] NetworkLayers = [25, 20, 15, 10, 5]; // DIA-NN 2020's architecture

        /// <summary>Fewer training positives than this and the training cutoff is relaxed.</summary>
        internal const int MinimumPositives = 10;

        private static readonly double[] RelaxedCutoffs = [0.05, 0.10, 0.25, 0.50];

        /// <summary>
        /// The q-value cutoff for choosing training positives: the requested one, relaxed step by step only while fewer than
        /// <see cref="MinimumPositives"/> targets pass. Weak features rarely pass targets at 1% on the first iteration, so this
        /// seeds a first model, as Percolator and mokapot do. It affects only which training rows are positives; the reported
        /// q-values are the caller's.
        /// </summary>
        internal static double TrainingCutoff(double[] q, bool[] isDecoy, double requested)
        {
            int Passing(double cutoff) => Enumerable.Range(0, q.Length).Count(t => !isDecoy[t] && q[t] <= cutoff);
            if (Passing(requested) >= MinimumPositives)
                return requested;
            foreach (double relaxed in RelaxedCutoffs.Where(c => c > requested))
                if (Passing(relaxed) >= MinimumPositives)
                    return relaxed;
            return RelaxedCutoffs[^1];
        }

        /// <summary>
        /// The single feature, and sign, that passes the most training targets. It uses the requested cutoff, or the first
        /// relaxed one at which some feature passes <see cref="MinimumPositives"/>. Ties, including the case where nothing
        /// passes anywhere, go to the larger standardized target-minus-decoy mean difference. That difference always carries the
        /// right sign, so a direction is never chosen by the order the features happen to be in.
        /// </summary>
        internal static (int Feature, int Sign) BestSingleFeature(IReadOnlyList<double[]> features, IReadOnlyList<bool> isDecoy, int[] train, double cutoff,
            IReadOnlyList<int>? candidateGroups = null, int maxParallelism = -1)
            => BestSingleFeature(features, isDecoy, RowGroups.Of(train, candidateGroups), cutoff, maxParallelism);

        /// <summary><see cref="BestSingleFeature(IReadOnlyList{double[]}, IReadOnlyList{bool}, int[], double, IReadOnlyList{int}?, int)"/> on rows already grouped.</summary>
        internal static (int Feature, int Sign) BestSingleFeature<T>(IReadOnlyList<T[]> features, IReadOnlyList<bool> isDecoy, RowGroups groups, double cutoff,
            int maxParallelism) where T : unmanaged, System.Numerics.IBinaryFloatingPointIeee754<T>
        {
            int p = features.Count == 0 ? 0 : features[0].Length;
            int rowCount = groups.Count; // the rows each feature and sign judges
            // The requested cutoff, then the relaxed ones above it: the choice below needs only each slot's passing count at
            // these, so a slot keeps its counts, not its q-values
            double[] cutoffs = new[] { cutoff }.Concat(RelaxedCutoffs.Where(c => c > cutoff)).ToArray();
            // Every feature and sign is judged independently, so in parallel, each into its own slot (same order as before).
            // Each worker ranks into a set of buffers from a pool, reused for every slot it takes and handed back when it
            // stops (the loop may retire a worker and start another): fresh arrays per slot, millions of rows each for 152
            // slots in each of three folds, were most of a DIA search's peak memory.
            var slots = new (int Feature, int Sign, double Separation, int[] Passing)[2 * p];
            var pool = new System.Collections.Concurrent.ConcurrentBag<SeedBuffers>();
            System.Threading.Tasks.Parallel.For(0, 2 * p, new System.Threading.Tasks.ParallelOptions { MaxDegreeOfParallelism = maxParallelism },
                () => pool.TryTake(out var reused) ? reused : new SeedBuffers(rowCount),
                (k, _, buffers) =>
                {
                    int j = k / 2, sign = k % 2 == 0 ? 1 : -1;
                    // With candidate groups, each feature and sign judges only the group rows it ranks top
                    int count = TopPerGroupInto(groups, i => sign * double.CreateTruncating(features[i][j]), buffers.Rows);
                    for (int t = 0; t < count; t++)
                    {
                        int row = buffers.Rows[t];
                        buffers.Decoy[t] = isDecoy[row];
                        buffers.Values[t] = double.CreateTruncating(features[row][j]);
                    }
                    double separation = StandardizedMeanDifference(buffers.Values, buffers.Decoy, count);
                    // The signed values in place (sign * v, exactly): the unsigned ones are not needed again
                    if (sign < 0)
                        for (int t = 0; t < count; t++)
                            buffers.Values[t] = sign * buffers.Values[t];
                    // Each target's q-value is read at its rank, the order the ranking leaves them in: the counts need no
                    // q-values in row order
                    TargetDecoyQValues.Rank(buffers.Values, buffers.Decoy, count, buffers.Order, buffers.Ranked);
                    var passing = new int[cutoffs.Length];
                    for (int r = 0; r < count; r++)
                        if (!buffers.Decoy[buffers.Order[r]])
                            for (int c = 0; c < cutoffs.Length; c++)
                                if (buffers.Ranked[r] <= cutoffs[c])
                                    passing[c]++;
                    slots[k] = (j, sign, sign * separation, passing);
                    return buffers;
                },
                pool.Add);

            // The first cutoff at which some slot passes enough targets; the requested one when none does
            int chosen = 0;
            for (int c = 0; c < cutoffs.Length; c++)
                if (slots.Any(s => s.Passing[c] >= MinimumPositives))
                {
                    chosen = c;
                    break;
                }
            var best = slots
                .OrderByDescending(s => s.Passing[chosen])
                .ThenByDescending(s => s.Separation)
                .ThenBy(s => s.Feature).ThenByDescending(s => s.Sign)
                .First();
            return (best.Feature, best.Sign);
        }

        /// <summary>
        /// The rows that train: all of them, or, with candidate groups, each group's top-scoring row (ties go to the earlier row),
        /// in input order.
        /// </summary>
        private static int[] TopPerGroup(RowGroups groups, Func<int, double> score)
        {
            if (groups.MultipleCount == 0)
                return groups.Singles;
            var top = new int[groups.Count];
            TopPerGroupInto(groups, score, top);
            return top;
        }

        /// <summary>
        /// <see cref="TopPerGroup"/> written into <paramref name="into"/> (at least one entry per group), in ascending row order;
        /// returns how many rows it wrote.
        /// </summary>
        private static int TopPerGroupInto(RowGroups groups, Func<int, double> score, int[] into)
        {
            groups.Singles.CopyTo(into, 0);
            int k = groups.Singles.Length;
            int[] rows = groups.MultipleRows, starts = groups.MultipleStarts;
            for (int g = 0; g < groups.MultipleCount; g++)
            {
                // Highest score; ties go to the earlier row (a group's rows are in ascending order)
                int best = rows[starts[g]];
                double bestScore = score(best);
                for (int m = starts[g] + 1; m < starts[g + 1]; m++)
                {
                    double s = score(rows[m]);
                    if (s > bestScore)
                    {
                        best = rows[m];
                        bestScore = s;
                    }
                }
                into[k++] = best;
            }
            // Distinct rows, so any sort gives the one ascending order (the singles alone are ascending already)
            if (groups.MultipleCount > 0)
                into.AsSpan(0, k).Sort();
            return k;
        }

        /// <summary>
        /// The rows of <paramref name="sampleSize"/> groups drawn at random, ascending. Each group, the singles first and then
        /// the groups of two or more, takes the next draw of <paramref name="sampler"/> in turn; the lowest draws win, and a tie
        /// goes to the earlier group. That is <c>OrderBy(_ =&gt; sampler.Next()).Take(sampleSize)</c> over the groups, which
        /// this replaces, without a list per group.
        /// </summary>
        internal static int[] SampleGroupRows(RowGroups groups, int sampleSize, Random sampler)
        {
            int count = groups.Count;
            var draws = new int[count];
            var order = new int[count];
            for (int g = 0; g < count; g++)
            {
                draws[g] = sampler.Next();
                order[g] = g;
            }
            order.AsSpan().Sort((a, b) =>
            {
                int c = draws[a].CompareTo(draws[b]);
                return c != 0 ? c : a.CompareTo(b);
            });
            int take = Math.Min(sampleSize, count), singles = groups.Singles.Length, total = 0;
            for (int k = 0; k < take; k++)
                total += order[k] < singles ? 1 : groups.Multiple(order[k] - singles).Length;
            var rows = new int[total];
            int at = 0;
            for (int k = 0; k < take; k++)
            {
                if (order[k] < singles)
                    rows[at++] = groups.Singles[order[k]];
                else
                {
                    var members = groups.Multiple(order[k] - singles);
                    members.CopyTo(rows.AsSpan(at));
                    at += members.Length;
                }
            }
            rows.AsSpan().Sort();
            return rows;
        }

        /// <summary>One worker's buffers for the seed step (BestSingleFeature), each holding one feature and sign's rows.</summary>
        private sealed class SeedBuffers(int rows)
        {
            public int[] Rows { get; } = new int[rows];
            public bool[] Decoy { get; } = new bool[rows];
            public double[] Values { get; } = new double[rows];
            public int[] Order { get; } = new int[rows];
            public double[] Ranked { get; } = new double[rows];
        }

        /// <summary>
        /// Rows split by candidate group once, so that picking each group's top row under a new score is one pass with no
        /// regrouping. Without candidate groups every row is its own group. The groups of two or more rows are held flat, one
        /// array of their rows and one of where each starts: a list per group was millions of small objects per fold.
        /// </summary>
        internal sealed class RowGroups
        {
            private RowGroups(int[] singles, int[] multipleRows, int[] multipleStarts)
            {
                Singles = singles;
                MultipleRows = multipleRows;
                MultipleStarts = multipleStarts;
            }

            /// <summary>Rows alone in their group, in ascending order.</summary>
            public int[] Singles { get; }

            /// <summary>
            /// The rows of every group of two or more, group after group in the order each group first appears in the rows
            /// given, each group's rows ascending.
            /// </summary>
            public int[] MultipleRows { get; }

            /// <summary>Where each group of two or more starts in <see cref="MultipleRows"/>, and one past the last.</summary>
            public int[] MultipleStarts { get; }

            /// <summary>The groups of two or more rows.</summary>
            public int MultipleCount => MultipleStarts.Length - 1;

            /// <summary>Every group: the singles, then the groups of two or more.</summary>
            public int Count => Singles.Length + MultipleCount;

            /// <summary>The rows of group <paramref name="g"/> of two or more (0 is the first of them).</summary>
            public ReadOnlySpan<int> Multiple(int g) => MultipleRows.AsSpan(MultipleStarts[g], MultipleStarts[g + 1] - MultipleStarts[g]);

            public static RowGroups Of(int[] rows, IReadOnlyList<int>? candidateGroups)
            {
                if (candidateGroups is null)
                    return new RowGroups(rows, [], [0]);

                // Each group's number, in the order it first appears; an array keyed by group when the groups are dense
                // (precursor indices are), a dictionary otherwise
                var number = new int[rows.Length];
                var sizes = new List<int>();
                int low = int.MaxValue, high = int.MinValue;
                foreach (int i in rows)
                {
                    low = Math.Min(low, candidateGroups[i]);
                    high = Math.Max(high, candidateGroups[i]);
                }
                if (rows.Length > 0 && (long)high - low < 4L * rows.Length + 1024)
                {
                    var numberOf = new int[high - low + 1];
                    Array.Fill(numberOf, -1);
                    for (int t = 0; t < rows.Length; t++)
                    {
                        ref int g = ref numberOf[candidateGroups[rows[t]] - low];
                        if (g < 0)
                        {
                            g = sizes.Count;
                            sizes.Add(0);
                        }
                        number[t] = g;
                        sizes[g]++;
                    }
                }
                else
                {
                    var numberOf = new Dictionary<int, int>();
                    for (int t = 0; t < rows.Length; t++)
                    {
                        if (!numberOf.TryGetValue(candidateGroups[rows[t]], out int g))
                        {
                            numberOf[candidateGroups[rows[t]]] = g = sizes.Count;
                            sizes.Add(0);
                        }
                        number[t] = g;
                        sizes[g]++;
                    }
                }

                // The groups of two or more keep their order of first appearance; each one's rows, written in the order
                // given, are then sorted
                var multipleIndex = new int[sizes.Count];
                int multiples = 0, singleCount = 0;
                for (int g = 0; g < sizes.Count; g++)
                {
                    if (sizes[g] == 1)
                    {
                        multipleIndex[g] = -1;
                        singleCount++;
                    }
                    else
                        multipleIndex[g] = multiples++;
                }
                var starts = new int[multiples + 1];
                for (int g = 0; g < sizes.Count; g++)
                    if (multipleIndex[g] >= 0)
                        starts[multipleIndex[g] + 1] = sizes[g];
                for (int m = 0; m < multiples; m++)
                    starts[m + 1] += starts[m];
                var singles = new int[singleCount];
                var multipleRows = new int[starts[multiples]];
                var cursor = starts[..^1];
                int s = 0;
                for (int t = 0; t < rows.Length; t++)
                {
                    int m = multipleIndex[number[t]];
                    if (m < 0)
                        singles[s++] = rows[t];
                    else
                        multipleRows[cursor[m]++] = rows[t];
                }
                singles.AsSpan().Sort();
                for (int m = 0; m < multiples; m++)
                    multipleRows.AsSpan(starts[m], starts[m + 1] - starts[m]).Sort();
                return new RowGroups(singles, multipleRows, starts);
            }
        }

        /// <summary>(mean of targets − mean of decoys) / pooled SD; 0 when either class is missing or the feature is constant.</summary>
        internal static double StandardizedMeanDifference(double[] values, bool[] isDecoy)
            => StandardizedMeanDifference(values, isDecoy, values.Length);

        /// <summary>
        /// <see cref="StandardizedMeanDifference(double[], bool[])"/> over the first <paramref name="count"/> entries. Every mean
        /// and sum is taken as LINQ's Average and Sum take it over an array, one addition at a time in row order, so the result
        /// is bit for bit the same, without the class arrays: this runs for every feature and sign of every fold.
        /// </summary>
        internal static double StandardizedMeanDifference(double[] values, bool[] isDecoy, int count)
        {
            double targetSum = 0, decoySum = 0, total = 0;
            int targets = 0, decoys = 0;
            for (int t = 0; t < count; t++)
            {
                total += values[t];
                if (isDecoy[t])
                {
                    decoySum += values[t];
                    decoys++;
                }
                else
                {
                    targetSum += values[t];
                    targets++;
                }
            }
            if (targets == 0 || decoys == 0)
                return 0;
            double mean = total / count;
            double squares = 0;
            for (int t = 0; t < count; t++)
                squares += (values[t] - mean) * (values[t] - mean);
            double sd = Math.Sqrt(squares / Math.Max(1, count - 1));
            return sd > LinearDiscriminant.ConstantTolerance ? (targetSum / targets - decoySum / decoys) / sd : 0;
        }

        /// <summary>Target-decoy q-values in input order, from the shared <see cref="TargetDecoyQValues"/>.</summary>
        internal static double[] QValues(double[] scores, bool[] isDecoy) => TargetDecoyQValues.Compute(scores, isDecoy);

        /// <summary>
        /// The pooled confident sample: the first <paramref name="cap"/> of <paramref name="ranked"/> (best first), in ascending
        /// row order. If those are all one label, the network cannot train on them, so the top targets and decoys, half each.
        /// </summary>
        internal static int[] PooledTrainingRows(int[] ranked, IReadOnlyList<bool> isDecoy, int cap)
        {
            int[] top = ranked.Take(cap).ToArray();
            if (top.Any(i => isDecoy[i]) && top.Any(i => !isDecoy[i]))
                return top.Order().ToArray();
            return ranked.Where(i => !isDecoy[i]).Take(cap / 2).Concat(ranked.Where(i => isDecoy[i]).Take(cap - cap / 2)).Order().ToArray();
        }

        /// <summary>
        /// A random <paramref name="cap"/> of <paramref name="rows"/>, in ascending row order. The network cannot train on one
        /// class, so a draw that holds one label while the rows hold both gives up its last row to the first row of the other
        /// label in the same random order. A draw that holds both is unchanged.
        /// </summary>
        internal static int[] RandomTrainingRows(int[] rows, IReadOnlyList<bool> isDecoy, int cap, Random sampler)
        {
            int[] shuffled = rows.OrderBy(_ => sampler.Next()).ToArray();
            int[] sample = shuffled.Take(cap).ToArray();
            bool hasTarget = sample.Any(i => !isDecoy[i]), hasDecoy = sample.Any(i => isDecoy[i]);
            if (hasTarget != hasDecoy && sample.Length < shuffled.Length)
            {
                // The missing label is the decoys when only targets were drawn, and the other way round
                int other = Array.FindIndex(shuffled, sample.Length, i => isDecoy[i] == hasTarget);
                if (other >= 0)
                    sample[^1] = shuffled[other];
            }
            return sample.Order().ToArray();
        }

        /// <summary>
        /// The confident training sample from rows <paramref name="ranked"/> best first (their scores in
        /// <paramref name="rankedScores"/>): the top half-cap of targets and of decoys. With <paramref name="positiveQValue"/>,
        /// only targets whose q-value in that ranking passes it are positives (at most half the cap), with as many of the top
        /// decoys; if none passes, the plain sample. Returned in ascending row order.
        /// </summary>
        internal static int[] ConfidentTrainingRows(int[] ranked, double[] rankedScores, IReadOnlyList<bool> isDecoy, int cap, double? positiveQValue,
            double targetFraction = 0.5)
        {
            int targetCount = (int)Math.Floor(cap * targetFraction + 1e-9); // cap / 2 at one half, as before
            // The network cannot train on one class: a small cap or share that rounds one away keeps a row of each (cap >= 2)
            if (ranked.Any(i => !isDecoy[i]) && ranked.Any(i => isDecoy[i]))
                targetCount = Math.Clamp(targetCount, 1, cap - 1);
            if (positiveQValue is double cutoff)
            {
                double[] q = QValues(rankedScores, ranked.Select(i => isDecoy[i]).ToArray());
                int[] targets = ranked.Where((i, r) => !isDecoy[i] && q[r] <= cutoff).Take(cap / 2).ToArray();
                if (targets.Length > 0)
                    return targets.Concat(ranked.Where(i => isDecoy[i]).Take(Math.Min(cap - cap / 2, targets.Length))).Order().ToArray();
            }
            return ranked.Where(i => !isDecoy[i]).Take(targetCount).Concat(ranked.Where(i => isDecoy[i]).Take(cap - targetCount)).Order().ToArray();
        }
    }
}
