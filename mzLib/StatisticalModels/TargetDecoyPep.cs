using System;
using System.Collections.Generic;
using System.Linq;

namespace StatisticalModels
{
    /// <summary>How a target-decoy PEP estimate ended.</summary>
    public enum TargetDecoyPepStatus
    {
        /// <summary>Every candidate has a PEP.</summary>
        Estimated,

        /// <summary>There were no decoys, so nothing could be estimated; every PEP is NaN.</summary>
        NoDecoys,

        /// <summary>There were no targets, so nothing could be estimated; every PEP is NaN.</summary>
        NoTargets,
    }

    /// <summary>The estimated posterior error probability of each candidate, and how the estimate ended.</summary>
    public sealed class TargetDecoyPepResult
    {
        internal TargetDecoyPepResult(double[] peps, TargetDecoyPepStatus status)
        {
            Peps = peps;
            Status = status;
        }

        /// <summary>
        /// One PEP per candidate, in input order, in [0, 1]. A decoy gets the value at its score. Every value is NaN when
        /// <see cref="Status"/> is not <see cref="TargetDecoyPepStatus.Estimated"/>.
        /// </summary>
        public double[] Peps { get; }

        /// <summary>Whether the estimate could be made.</summary>
        public TargetDecoyPepStatus Status { get; }
    }

    /// <summary>
    /// Posterior error probabilities read off the target-decoy distribution of a score.
    /// <para>
    /// A classifier's probability is calibrated to the class balance of the rows it was trained on, which is rarely the
    /// population it is then applied to, so it is not the probability that a match is wrong. This estimate keeps the
    /// score's ranking and reads the error probability off the targets and decoys that share it, as Percolator's qvality
    /// and mokapot do.
    /// </para>
    /// <list type="number">
    /// <item>Rank the candidates by score, best first. Exactly tied scores form one block and get one value.</item>
    /// <item>Fit the decoy fraction <c>r</c> along that ranking by weighted isotonic regression (pool adjacent violators),
    /// constrained not to fall as the score falls.</item>
    /// <item><c>PEP = min(1, r / (1 − r))</c>: the local decoy-to-target ratio. Under target-decoy competition each decoy
    /// stands for one false target, the same assumption behind the q-value (D + 1) / T.</item>
    /// </list>
    /// <para>
    /// There is no bandwidth, bin width or smoothing constant. Two consequences follow from the data, not from a choice:
    /// candidates scoring above every decoy get PEP 0, and wherever decoys are at least as frequent as targets the PEP is 1.
    /// The PEP never rises as the score rises.
    /// </para>
    /// </summary>
    public static class TargetDecoyPep
    {
        /// <summary>Estimates a PEP for each candidate from its score and whether it is a decoy.</summary>
        /// <param name="scores">One finite score per candidate; higher is better.</param>
        /// <param name="isDecoy">Whether each candidate is a decoy.</param>
        /// <returns>A PEP per candidate in input order, and the status.</returns>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">The lengths disagree, or a score is not finite.</exception>
        public static TargetDecoyPepResult Estimate(IReadOnlyList<double> scores, IReadOnlyList<bool> isDecoy)
        {
            ArgumentNullException.ThrowIfNull(scores);
            ArgumentNullException.ThrowIfNull(isDecoy);
            return Estimate(scores, isDecoy.Select(d => d ? 1.0 : 0.0).ToArray());
        }

        /// <summary>
        /// Estimates a PEP for each candidate from its score and the share of it that is a decoy. A candidate ambiguous
        /// between one target and one decoy hypothesis, for example, is half a decoy (0.5), and is counted that way.
        /// </summary>
        /// <param name="scores">One finite score per candidate; higher is better.</param>
        /// <param name="decoyFraction">The decoy share of each candidate, in [0, 1]: 0 is a target, 1 a decoy.</param>
        /// <returns>A PEP per candidate in input order, and the status.</returns>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">
        /// The lengths disagree, a score is not finite, or a decoy fraction is not in [0, 1].
        /// </exception>
        public static TargetDecoyPepResult Estimate(IReadOnlyList<double> scores, IReadOnlyList<double> decoyFraction)
        {
            ArgumentNullException.ThrowIfNull(scores);
            ArgumentNullException.ThrowIfNull(decoyFraction);
            if (scores.Count != decoyFraction.Count)
                throw new ArgumentException($"There are {scores.Count} scores but {decoyFraction.Count} decoy fractions.", nameof(decoyFraction));
            double decoys = 0, targets = 0;
            for (int i = 0; i < scores.Count; i++)
            {
                if (!double.IsFinite(scores[i]))
                    throw new ArgumentException($"Score {i} is {scores[i]}; every score must be finite.", nameof(scores));
                double f = decoyFraction[i];
                if (!(f >= 0 && f <= 1))
                    throw new ArgumentException($"Decoy fraction {i} is {f}; it must be in [0, 1].", nameof(decoyFraction));
                decoys += f;
                targets += 1 - f;
            }

            var peps = new double[scores.Count];
            if (decoys == 0 || targets == 0)
            {
                Array.Fill(peps, double.NaN);
                return new TargetDecoyPepResult(peps, decoys == 0 ? TargetDecoyPepStatus.NoDecoys : TargetDecoyPepStatus.NoTargets);
            }

            // Best first; the index makes the order total, so the result never depends on how the sort breaks ties.
            int[] order = Enumerable.Range(0, scores.Count).OrderByDescending(i => scores[i]).ThenBy(i => i).ToArray();

            // One block per distinct score: its weight is its size, its value the decoy share within it.
            var blockStart = new List<int>();
            var blockDecoys = new List<double>();
            var blockWeight = new List<double>();
            for (int k = 0; k < order.Length;)
            {
                int end = k;
                double d = 0;
                do
                {
                    d += decoyFraction[order[end]];
                    end++;
                }
                while (end < order.Length && scores[order[end]] == scores[order[k]]);
                blockStart.Add(k);
                blockDecoys.Add(d);
                blockWeight.Add(end - k);
                k = end;
            }

            double[] rate = IsotonicNonDecreasing(
                blockDecoys.Zip(blockWeight, (d, w) => d / w).ToArray(), blockWeight.ToArray());

            for (int b = 0; b < rate.Length; b++)
            {
                double pep = rate[b] >= 0.5 ? 1.0 : rate[b] / (1.0 - rate[b]);
                int end = b + 1 < rate.Length ? blockStart[b + 1] : order.Length;
                for (int k = blockStart[b]; k < end; k++)
                    peps[order[k]] = pep;
            }
            return new TargetDecoyPepResult(peps, TargetDecoyPepStatus.Estimated);
        }

        /// <summary>
        /// Weighted isotonic regression by pool adjacent violators: the non-decreasing sequence closest to
        /// <paramref name="values"/> in weighted least squares. Adjacent values that violate the order are pooled into their
        /// weighted mean until none does, so the fit is a step function with flat plateaus.
        /// </summary>
        /// <param name="values">The values, in order.</param>
        /// <param name="weights">A positive weight per value.</param>
        /// <returns>The fitted value at each position.</returns>
        internal static double[] IsotonicNonDecreasing(double[] values, double[] weights)
        {
            // Each pooled block keeps its weighted sum and total weight, so its mean is formed once, never updated in place.
            var sum = new double[values.Length];
            var weight = new double[values.Length];
            var size = new int[values.Length];
            int top = -1;
            for (int i = 0; i < values.Length; i++)
            {
                top++;
                sum[top] = values[i] * weights[i];
                weight[top] = weights[i];
                size[top] = 1;
                while (top > 0 && sum[top - 1] / weight[top - 1] > sum[top] / weight[top])
                {
                    sum[top - 1] += sum[top];
                    weight[top - 1] += weight[top];
                    size[top - 1] += size[top];
                    top--;
                }
            }
            var fitted = new double[values.Length];
            int position = 0;
            for (int b = 0; b <= top; b++)
            {
                double mean = sum[b] / weight[b];
                for (int j = 0; j < size[b]; j++)
                    fitted[position++] = mean;
            }
            return fitted;
        }
    }
}
