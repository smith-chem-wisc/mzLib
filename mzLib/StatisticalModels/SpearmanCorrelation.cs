using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;

namespace StatisticalModels
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

        /// <summary>
        /// Partial Spearman correlation of <paramref name="x"/> and <paramref name="y"/> given one covariate
        /// <paramref name="z"/>: <see cref="PartialCorrelate(IReadOnlyList{double}, IReadOnlyList{double}, IReadOnlyList{IReadOnlyList{double}})"/>
        /// with one covariate, where ρ_xy·z = (ρ_xy − ρ_xz ρ_yz) / √((1 − ρ_xz²)(1 − ρ_yz²)) and the test has n − 3 df.
        /// </summary>
        public static SpearmanResult PartialCorrelate(IReadOnlyList<double> x, IReadOnlyList<double> y, IReadOnlyList<double> z)
        {
            ArgumentNullException.ThrowIfNull(z);
            return PartialCorrelate(x, y, new[] { z });
        }

        /// <summary>
        /// Partial Spearman correlation of <paramref name="x"/> and <paramref name="y"/> given the covariates
        /// <paramref name="covariates"/>, and its two-sided p-value. Every variable is ranked over the complete cases
        /// (all values finite). ρ_xy·z is the correlation of the ranks of x and y after each is regressed on the ranks of
        /// the covariates; with R the rank correlation matrix of (x, y, z₁ … z_k) and P = R⁻¹, it is −P₀₁ / √(P₀₀ P₁₁).
        /// p: t = ρ √((n − 2 − k) / (1 − ρ²)) referred to Student's t with n − 2 − k df, as R's
        /// ppcor::pcor.test(method = "spearman"). If x and y are left perfectly correlated, ρ = ±1 and p = 0, as
        /// <see cref="Correlate"/> gives. A covariate whose ranks the others carry entirely (two level series that rank
        /// the runs identically) adds nothing and is dropped, and k counts the covariates kept.
        /// <see cref="SpearmanResult.Method"/> is <see cref="SpearmanPValueMethod.Asymptotic"/>, or NotEstimable with
        /// fewer than k + 3 complete cases (k the covariates kept, so a duplicated covariate needs no extra case), a
        /// constant variable, or covariates that carry x or y entirely (nothing of it is left to correlate).
        /// </summary>
        public static SpearmanResult PartialCorrelate(IReadOnlyList<double> x, IReadOnlyList<double> y, IReadOnlyList<IReadOnlyList<double>> covariates)
        {
            ArgumentNullException.ThrowIfNull(x);
            ArgumentNullException.ThrowIfNull(y);
            ArgumentNullException.ThrowIfNull(covariates);
            if (covariates.Count == 0)
                throw new ArgumentException("At least one covariate is needed; use Correlate for none.", nameof(covariates));
            if (x.Count != y.Count || covariates.Any(z => z is null || z.Count != x.Count))
                throw new ArgumentException($"x has {x.Count} values; y and every covariate must have as many.", nameof(covariates));

            int k = covariates.Count;
            var variables = new IReadOnlyList<double>[k + 2];
            variables[0] = x;
            variables[1] = y;
            for (int j = 0; j < k; j++) variables[j + 2] = covariates[j];
            var complete = Enumerable.Range(0, x.Count).Where(i => variables.All(v => double.IsFinite(v[i]))).ToArray();
            int n = complete.Length;
            var result = new SpearmanResult { N = n, Rho = double.NaN, PValue = double.NaN, Method = SpearmanPValueMethod.NotEstimable };
            // At least one covariate is always kept, so fewer than 4 cases leave no df. The full check, on the
            // covariates kept, follows the basis: a covariate the others carry does not cost a case.
            if (n < 4) return result;

            var ranks = new double[k + 2][];
            bool ties = false;
            for (int j = 0; j < k + 2; j++)
            {
                ranks[j] = Ranks(complete.Select(i => variables[j][i]).ToArray(), out bool t);
                ties |= t;
            }
            result.HasTies = ties;
            // Each rank vector is centred and scaled to unit length, and the covariates are orthonormalized in turn
            // (modified Gram-Schmidt). The test is on what is left of x and y after projecting them off that basis, so
            // only nothing left of x or of y refuses it: x and y left identical in rank give ρ = ±1, not a refusal.
            var unit = new double[k + 2][];
            for (int j = 0; j < k + 2; j++)
            {
                double mean = ranks[j].Average();
                unit[j] = ranks[j].Select(v => v - mean).ToArray();
                double norm = Math.Sqrt(Dot(unit[j], unit[j]));
                if (norm == 0) return result;
                for (int i = 0; i < n; i++) unit[j][i] /= norm;
            }
            var basis = new List<double[]>(k);
            for (int j = 2; j < k + 2; j++)
            {
                var v = Residual(unit[j], basis);
                double left = Dot(v, v);
                if (left < CarriedTolerance) continue;
                double norm = Math.Sqrt(left);
                basis.Add(v.Select(e => e / norm).ToArray());
            }
            if (n < basis.Count + 3) return result;
            var ex = Residual(unit[0], basis);
            var ey = Residual(unit[1], basis);
            double sxx = Dot(ex, ex), syy = Dot(ey, ey);
            if (sxx < CarriedTolerance || syy < CarriedTolerance) return result;
            double rho = Math.Clamp(Dot(ex, ey) / Math.Sqrt(sxx * syy), -1, 1);
            result.Rho = rho;
            result.Method = SpearmanPValueMethod.Asymptotic;
            double r2 = rho * rho;
            int df = n - 2 - basis.Count;
            result.PValue = r2 >= 1 ? 0 : 2 * StudentT.CDF(0, 1, df, -Math.Abs(rho) * Math.Sqrt(df / (1 - r2)));
            return result;
        }

        /// <summary>
        /// Share of a rank vector's variance, left after the covariates, below which
        /// <see cref="PartialCorrelate(IReadOnlyList{double}, IReadOnlyList{double}, IReadOnlyList{IReadOnlyList{double}})"/>
        /// takes the covariates to carry it entirely. One adjacent swap in n ranks leaves about 24 / n³, so distinct rank
        /// orders stay above it for n up to about 28,000.
        /// </summary>
        private const double CarriedTolerance = 1e-12;

        /// <summary><paramref name="v"/> minus its projection on each vector of the orthonormal <paramref name="basis"/>, in turn.</summary>
        private static double[] Residual(double[] v, List<double[]> basis)
        {
            var e = (double[])v.Clone();
            foreach (var q in basis)
            {
                double c = Dot(e, q);
                for (int i = 0; i < e.Length; i++) e[i] -= c * q[i];
            }
            return e;
        }

        private static double Dot(double[] a, double[] b)
        {
            double s = 0;
            for (int i = 0; i < a.Length; i++) s += a[i] * b[i];
            return s;
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
