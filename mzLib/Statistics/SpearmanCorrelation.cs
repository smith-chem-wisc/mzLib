using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;

namespace Statistics
{
    /// <summary>How a Spearman p-value was computed.</summary>
    public enum SpearmanPValueMethod
    {
        /// <summary>
        /// Exact: the permutation distribution of S = Σ dᵢ², enumerated in full. Used only for
        /// n ≤ <see cref="SpearmanCorrelation.ExactMaxN"/> with no ties, where it equals R's
        /// cor.test(method = "spearman", exact = TRUE).
        /// </summary>
        Exact,
        /// <summary>
        /// Asymptotic: t = ρ √((n − 2) / (1 − ρ²)) referred to Student's t with n − 2 df, as R's
        /// cor.test(method = "spearman", exact = FALSE). Used for larger n and whenever there are ties.
        /// For 10 ≤ n &lt; 1290 without ties R's default is instead an Edgeworth series (AS 89), so
        /// p-values there differ slightly from R's default output.
        /// </summary>
        Asymptotic,
        /// <summary>No p-value: fewer than 3 complete pairs, or a constant variable.</summary>
        NotEstimable,
    }

    /// <summary>Spearman's rank correlation of one pair of variables.</summary>
    public sealed class SpearmanResult
    {
        internal SpearmanResult() { }

        /// <summary>Spearman's ρ: Pearson's correlation of the ranks, ties given their average rank. NaN if not estimable.</summary>
        public double Rho { get; internal set; }
        /// <summary>Two-sided p-value against ρ = 0. NaN if not estimable.</summary>
        public double PValue { get; internal set; }
        /// <summary>Number of complete pairs (both values finite) the correlation was computed on.</summary>
        public int N { get; internal set; }
        /// <summary>Whether either variable had tied values among the complete pairs.</summary>
        public bool HasTies { get; internal set; }
        /// <summary>How <see cref="PValue"/> was computed.</summary>
        public SpearmanPValueMethod Method { get; internal set; }
    }

    /// <summary>Spearman's rank correlation with a two-sided test of ρ = 0.</summary>
    /// <remarks>
    /// Pairs with a non-finite value on either side are omitted (pairwise-complete), never imputed. The
    /// estimator of the p-value is named on every result (<see cref="SpearmanResult.Method"/>).
    /// </remarks>
    public static class SpearmanCorrelation
    {
        /// <summary>Largest n for which the exact permutation p-value is enumerated (9! = 362,880 orderings).</summary>
        public const int ExactMaxN = 9;

        private static readonly ConcurrentDictionary<int, long[]> ExactCounts = new();

        /// <summary>Spearman's ρ of <paramref name="x"/> and <paramref name="y"/>, and its two-sided p-value.</summary>
        public static SpearmanResult Correlate(IReadOnlyList<double> x, IReadOnlyList<double> y)
        {
            ArgumentNullException.ThrowIfNull(x);
            ArgumentNullException.ThrowIfNull(y);
            if (x.Count != y.Count)
                throw new ArgumentException($"x has {x.Count} values and y has {y.Count}.", nameof(y));

            var xs = new List<double>(x.Count);
            var ys = new List<double>(x.Count);
            for (int i = 0; i < x.Count; i++)
                if (double.IsFinite(x[i]) && double.IsFinite(y[i])) { xs.Add(x[i]); ys.Add(y[i]); }
            int n = xs.Count;
            var result = new SpearmanResult { N = n, Rho = double.NaN, PValue = double.NaN, Method = SpearmanPValueMethod.NotEstimable };
            if (n < 3) return result;

            var rx = Ranks(xs, out bool tiesX);
            var ry = Ranks(ys, out bool tiesY);
            result.HasTies = tiesX || tiesY;
            double rho = Pearson(rx, ry);
            if (double.IsNaN(rho)) return result;
            result.Rho = rho;

            if (!result.HasTies && n <= ExactMaxN)
            {
                result.Method = SpearmanPValueMethod.Exact;
                result.PValue = ExactPValue(rx, ry, n);
            }
            else
            {
                result.Method = SpearmanPValueMethod.Asymptotic;
                double r2 = Math.Min(rho * rho, 1);
                result.PValue = r2 >= 1 ? 0 : 2 * StudentT.CDF(0, 1, n - 2, -Math.Abs(rho) * Math.Sqrt((n - 2) / (1 - r2)));
            }
            return result;
        }

        /// <summary>Average ranks, 1-based; ties share the mean of the ranks they span.</summary>
        internal static double[] Ranks(IReadOnlyList<double> v, out bool ties)
        {
            int n = v.Count;
            var order = Enumerable.Range(0, n).OrderBy(i => v[i]).ThenBy(i => i).ToArray();
            var r = new double[n];
            ties = false;
            for (int i = 0; i < n;)
            {
                int j = i;
                while (j + 1 < n && v[order[j + 1]] == v[order[i]]) j++;
                if (j > i) ties = true;
                double avg = (i + j) / 2.0 + 1;
                for (int k = i; k <= j; k++) r[order[k]] = avg;
                i = j + 1;
            }
            return r;
        }

        private static double Pearson(double[] a, double[] b)
        {
            double ma = a.Average(), mb = b.Average(), sab = 0, saa = 0, sbb = 0;
            for (int i = 0; i < a.Length; i++)
            {
                double da = a[i] - ma, db = b[i] - mb;
                sab += da * db; saa += da * da; sbb += db * db;
            }
            if (saa == 0 || sbb == 0) return double.NaN;
            return Math.Clamp(sab / Math.Sqrt(saa * sbb), -1, 1);
        }

        /// <summary>
        /// Two-sided exact p-value as R computes it: with S = Σ dᵢ² and its mean (n³ − n)/6, the tail on the
        /// side S falls, doubled and capped at 1.
        /// </summary>
        private static double ExactPValue(double[] rx, double[] ry, int n)
        {
            long s = 0;
            for (int i = 0; i < n; i++) { long d = (long)(rx[i] - ry[i]); s += d * d; }
            var counts = ExactCounts.GetOrAdd(n, Enumerate);
            double total = counts.Sum();
            double mean = (Math.Pow(n, 3) - n) / 6;
            double tail = 0;
            if (s > mean) { for (long k = s; k < counts.Length; k++) tail += counts[k]; }
            else { for (long k = 0; k <= s; k++) tail += counts[k]; }
            return Math.Min(1, 2 * tail / total);
        }

        /// <summary>Counts of S = Σ (i − π(i))² over all n! permutations π, indexed by S.</summary>
        private static long[] Enumerate(int n)
        {
            int max = 0;
            for (int i = 0; i < n; i++) max += (n - 1 - 2 * i) * (n - 1 - 2 * i);
            var counts = new long[max + 1];
            var used = new bool[n];
            void Recurse(int pos, int partial)
            {
                if (pos == n) { counts[partial]++; return; }
                for (int v = 0; v < n; v++)
                {
                    if (used[v]) continue;
                    used[v] = true;
                    Recurse(pos + 1, partial + (pos - v) * (pos - v));
                    used[v] = false;
                }
            }
            Recurse(0, 0);
            return counts;
        }
    }
}
