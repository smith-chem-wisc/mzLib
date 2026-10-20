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
                if (row is null || row.Length != p || !row.All(double.IsFinite))
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
                var seed = BestSingleFeature(features, isDecoy, train, positiveQValue, candidateGroups);
                Func<double[], double> scorer = x => seed.Sign * x[seed.Feature];
                MultilayerPerceptronEnsemble? network = null; // the scorer, when it is the network
                bool trained = false;
                for (int iteration = 0; iteration < iterations; iteration++)
                {
                    // Only each candidate group's top row under the current model trains it (pyProphet)
                    int[] active = TopPerGroup(trainGroups, i => scorer(features[i]));
                    double[] trainScores = active.Select(i => scorer(features[i])).ToArray();
                    double[] q = QValues(trainScores, active.Select(i => isDecoy[i]).ToArray());
                    double cutoff = TrainingCutoff(q, active.Select(i => isDecoy[i]).ToArray(), positiveQValue);
                    var rows = new List<double[]>();
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
                    var fit = LinearDiscriminant.Fit(rows, positive);
                    scorer = x => fit.Score(x);
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
                            var sampler = new Random(31 + f + 1000 * randomSeed + 100_000 * pass);
                            rows = rows.OrderBy(_ => sampler.Next()).Take(cap2).Order().ToArray();
                        }
                        var ensemble = MultilayerPerceptron.TrainEnsemble(rows.Select(i => features[i]).ToList(), rows.Select(i => !isDecoy[i]).ToList(),
                            networkMembers, networkLayers ?? NetworkLayers, networkEpochs, seed: 17 + f + 1000 * randomSeed + 100_000 * pass);
                        scorer = x => ensemble.PredictLogit(x);
                        network = ensemble;
                    }
                }

                // Every row scored once, in parallel: the folds alone use three cores, and a network ensemble over millions
                // of rows otherwise costs more than training it
                var final = new double[n];
                var normalizers = trainGroups;
                if (normalizationGroups is int sampleSize && trainGroups.Singles.Length + trainGroups.Multiple.Length > sampleSize)
                {
                    // A random sample of training groups, seeded by the fold, sets the normalisation: scoring every training
                    // row only for that was twice the work of scoring the fold's own rows
                    var sampler = new Random(53 + f + 1000 * randomSeed);
                    int[] sampledRows = trainGroups.Singles.Select(i => new[] { i }).Concat(trainGroups.Multiple)
                        .OrderBy(_ => sampler.Next()).Take(sampleSize).SelectMany(g => g).Order().ToArray();
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
                        network.PredictLogits(features, block, final);
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
            IReadOnlyList<int>? candidateGroups = null)
        {
            int p = features.Count == 0 ? 0 : features[0].Length;
            var groups = RowGroups.Of(train, candidateGroups);
            // Every feature and sign is judged independently, so in parallel, each into its own slot (same order as before)
            var slots = new (int Feature, int Sign, double Separation, double[] Q, bool[] Decoy)[2 * p];
            System.Threading.Tasks.Parallel.For(0, 2 * p, k =>
            {
                int j = k / 2, sign = k % 2 == 0 ? 1 : -1;
                // With candidate groups, each feature and sign judges only the group rows it ranks top
                int[] rows = TopPerGroup(groups, i => sign * features[i][j]);
                bool[] rowDecoy = rows.Select(i => isDecoy[i]).ToArray();
                double[] values = rows.Select(i => features[i][j]).ToArray();
                double separation = StandardizedMeanDifference(values, rowDecoy);
                slots[k] = (j, sign, sign * separation, QValues(values.Select(v => sign * v).ToArray(), rowDecoy), rowDecoy);
            });
            var candidates = slots.ToList();

            int PassingOf((int Feature, int Sign, double Separation, double[] Q, bool[] Decoy) k, double c) =>
                Enumerable.Range(0, k.Q.Length).Count(t => !k.Decoy[t] && k.Q[t] <= c);
            double chosenCutoff = new[] { cutoff }.Concat(RelaxedCutoffs.Where(c => c > cutoff))
                .FirstOrDefault(c => candidates.Any(k => PassingOf(k, c) >= MinimumPositives), cutoff);
            var best = candidates
                .OrderByDescending(k => PassingOf(k, chosenCutoff))
                .ThenByDescending(k => k.Separation)
                .ThenBy(k => k.Feature).ThenByDescending(k => k.Sign)
                .First();
            return (best.Feature, best.Sign);
        }

        /// <summary>
        /// The rows that train: all of them, or, with candidate groups, each group's top-scoring row (ties go to the earlier row),
        /// in input order.
        /// </summary>
        private static int[] TopPerGroup(RowGroups groups, Func<int, double> score)
        {
            if (groups.Multiple.Length == 0)
                return groups.Singles;
            var top = new int[groups.Singles.Length + groups.Multiple.Length];
            groups.Singles.CopyTo(top, 0);
            int k = groups.Singles.Length;
            foreach (int[] members in groups.Multiple)
            {
                // Highest score; ties go to the earlier row (members are in ascending row order)
                int best = members[0];
                double bestScore = score(best);
                for (int m = 1; m < members.Length; m++)
                {
                    double s = score(members[m]);
                    if (s > bestScore)
                    {
                        best = members[m];
                        bestScore = s;
                    }
                }
                top[k++] = best;
            }
            Array.Sort(top);
            return top;
        }

        /// <summary>
        /// Rows split by candidate group once, so that picking each group's top row under a new score is one pass with no
        /// regrouping. Without candidate groups every row is its own group.
        /// </summary>
        private sealed class RowGroups
        {
            private RowGroups(int[] singles, int[][] multiple)
            {
                Singles = singles;
                Multiple = multiple;
            }

            /// <summary>Rows alone in their group, in ascending order.</summary>
            public int[] Singles { get; }

            /// <summary>Groups of two or more rows, each in ascending row order.</summary>
            public int[][] Multiple { get; }

            public static RowGroups Of(int[] rows, IReadOnlyList<int>? candidateGroups)
            {
                if (candidateGroups is null)
                    return new RowGroups(rows, []);
                var byGroup = new Dictionary<int, List<int>>();
                foreach (int i in rows)
                {
                    if (!byGroup.TryGetValue(candidateGroups[i], out var list))
                        byGroup[candidateGroups[i]] = list = [];
                    list.Add(i);
                }
                var singles = byGroup.Values.Where(g => g.Count == 1).Select(g => g[0]).Order().ToArray();
                var multiple = byGroup.Values.Where(g => g.Count > 1).Select(g => g.Order().ToArray()).ToArray();
                return new RowGroups(singles, multiple);
            }
        }

        /// <summary>(mean of targets − mean of decoys) / pooled SD; 0 when either class is missing or the feature is constant.</summary>
        internal static double StandardizedMeanDifference(double[] values, bool[] isDecoy)
        {
            double[] targets = values.Where((_, t) => !isDecoy[t]).ToArray();
            double[] decoys = values.Where((_, t) => isDecoy[t]).ToArray();
            if (targets.Length == 0 || decoys.Length == 0)
                return 0;
            double mean = values.Average();
            double sd = Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / Math.Max(1, values.Length - 1));
            return sd > LinearDiscriminant.ConstantTolerance ? (targets.Average() - decoys.Average()) / sd : 0;
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
        /// The confident training sample from rows <paramref name="ranked"/> best first (their scores in
        /// <paramref name="rankedScores"/>): the top half-cap of targets and of decoys. With <paramref name="positiveQValue"/>,
        /// only targets whose q-value in that ranking passes it are positives (at most half the cap), with as many of the top
        /// decoys; if none passes, the plain sample. Returned in ascending row order.
        /// </summary>
        internal static int[] ConfidentTrainingRows(int[] ranked, double[] rankedScores, IReadOnlyList<bool> isDecoy, int cap, double? positiveQValue,
            double targetFraction = 0.5)
        {
            int targetCount = (int)Math.Floor(cap * targetFraction + 1e-9); // cap / 2 at one half, as before
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
