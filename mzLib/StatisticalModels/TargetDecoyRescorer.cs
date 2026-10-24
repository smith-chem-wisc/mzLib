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
    public static class TargetDecoyRescorer
    {
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">Lengths disagree, rows are ragged, or a value is not finite.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A count or the q-value cutoff is out of range.</exception>
        public static RescoreResult Score(IReadOnlyList<double[]> features, IReadOnlyList<bool> isDecoy, IReadOnlyList<string> groupKeys,
            int folds = 3, int iterations = 3, double positiveQValue = 0.01)
        {
            ArgumentNullException.ThrowIfNull(features);
            ArgumentNullException.ThrowIfNull(isDecoy);
            ArgumentNullException.ThrowIfNull(groupKeys);
            if (isDecoy.Count != features.Count)
                throw new ArgumentException($"There are {features.Count} candidates but {isDecoy.Count} decoy labels.", nameof(isDecoy));
            if (groupKeys.Count != features.Count)
                throw new ArgumentException($"There are {features.Count} candidates but {groupKeys.Count} group keys.", nameof(groupKeys));
            if (folds < 2)
                throw new ArgumentOutOfRangeException(nameof(folds), folds, "There must be at least 2 folds.");
            if (iterations < 1)
                throw new ArgumentOutOfRangeException(nameof(iterations), iterations, "There must be at least 1 iteration.");
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

            var status = RescoreStatus.Rescored;
            for (int f = 0; f < folds; f++)
            {
                int[] train = Enumerable.Range(0, n).Where(i => fold[i] != f).ToArray();
                int[] test = Enumerable.Range(0, n).Where(i => fold[i] == f).ToArray();
                if (test.Length == 0)
                    continue;

                var seed = BestSingleFeature(features, isDecoy, train, positiveQValue);
                Func<double[], double> scorer = x => seed.Sign * x[seed.Feature];
                bool trained = false;
                for (int iteration = 0; iteration < iterations; iteration++)
                {
                    double[] trainScores = train.Select(i => scorer(features[i])).ToArray();
                    double[] q = QValues(trainScores, train.Select(i => isDecoy[i]).ToArray());
                    double cutoff = TrainingCutoff(q, train.Select(i => isDecoy[i]).ToArray(), positiveQValue);
                    var rows = new List<double[]>();
                    var positive = new List<bool>();
                    for (int t = 0; t < train.Length; t++)
                    {
                        bool decoy = isDecoy[train[t]];
                        if (decoy || q[t] <= cutoff)
                        {
                            rows.Add(features[train[t]]);
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
                    status = RescoreStatus.FoldStarved;

                // Normalize on the training rows so that folds are comparable when pooled
                double[] finalTrain = train.Select(i => scorer(features[i])).ToArray();
                bool[] trainDecoy = train.Select(i => isDecoy[i]).ToArray();
                double[] finalQ = QValues(finalTrain, trainDecoy);
                double[] decoyScores = finalTrain.Where((_, t) => trainDecoy[t]).Order().ToArray();
                double medianDecoy = decoyScores.Length == 0 ? 0 : decoyScores[decoyScores.Length / 2];
                double[] passing = finalTrain.Where((_, t) => !trainDecoy[t] && finalQ[t] <= positiveQValue).ToArray();
                double threshold = passing.Length > 0 ? passing.Min() : (decoyScores.Length > 0 ? decoyScores[^1] : medianDecoy + 1);
                double scale = threshold > medianDecoy ? threshold - medianDecoy : 1;
                foreach (int i in test)
                    scores[i] = (scorer(features[i]) - threshold) / scale;
            }
            return new RescoreResult(scores, fold, status);
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
        internal static (int Feature, int Sign) BestSingleFeature(IReadOnlyList<double[]> features, IReadOnlyList<bool> isDecoy, int[] train, double cutoff)
        {
            bool[] decoy = train.Select(i => isDecoy[i]).ToArray();
            int p = features.Count == 0 ? 0 : features[0].Length;
            var candidates = new List<(int Feature, int Sign, double Separation, double[] Q)>();
            for (int j = 0; j < p; j++)
            {
                double[] values = train.Select(i => features[i][j]).ToArray();
                double separation = StandardizedMeanDifference(values, decoy);
                foreach (int sign in new[] { 1, -1 })
                    candidates.Add((j, sign, sign * separation, QValues(values.Select(v => sign * v).ToArray(), decoy)));
            }

            int Passing(double[] q, double c) => Enumerable.Range(0, q.Length).Count(t => !decoy[t] && q[t] <= c);
            double chosenCutoff = new[] { cutoff }.Concat(RelaxedCutoffs.Where(c => c > cutoff))
                .FirstOrDefault(c => candidates.Any(k => Passing(k.Q, c) >= MinimumPositives), cutoff);
            var best = candidates
                .OrderByDescending(k => Passing(k.Q, chosenCutoff))
                .ThenByDescending(k => k.Separation)
                .ThenBy(k => k.Feature).ThenByDescending(k => k.Sign)
                .First();
            return (best.Feature, best.Sign);
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

        /// <summary>Target-decoy q-values, (D + 1) / T and monotone, in input order. Decoys receive the value at their rank.</summary>
        internal static double[] QValues(double[] scores, bool[] isDecoy)
        {
            int[] order = Enumerable.Range(0, scores.Length).OrderByDescending(i => scores[i]).ThenBy(i => i).ToArray();
            var ranked = new double[order.Length];
            int decoys = 0, targets = 0;
            for (int k = 0; k < order.Length; k++)
            {
                if (isDecoy[order[k]]) decoys++; else targets++;
                ranked[k] = targets == 0 ? 1 : Math.Min(1, (decoys + 1.0) / targets);
            }
            for (int k = order.Length - 2; k >= 0; k--)
                ranked[k] = Math.Min(ranked[k], ranked[k + 1]);
            var q = new double[scores.Length];
            for (int k = 0; k < order.Length; k++)
                q[order[k]] = ranked[k];
            return q;
        }
    }
}
