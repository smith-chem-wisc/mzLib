using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;

namespace StatisticalModels
{
    /// <summary>Whether <see cref="EmpiricalNull.Fit"/> could fit the empirical null.</summary>
    public enum EmpiricalNullStatus
    {
        /// <summary>The empirical null was fitted.</summary>
        Estimated,

        /// <summary>Fewer than <see cref="EmpiricalNull.MinimumScores"/> scores in all, or in a window. Every estimate is
        /// NaN.</summary>
        TooFewScores,

        /// <summary>No normal fits the scores in a window: the likelihood rises without bound as the spread grows, so the
        /// window holds scores at least as spread out as a uniform over it. Every estimate is NaN.</summary>
        NoMaximum,

        /// <summary>Newton's method did not settle within <see cref="EmpiricalNull.MaxIterations"/> steps. Every
        /// estimate is NaN.</summary>
        NotConverged,
    }

    /// <summary>An empirical null N(δ0, σ0²) fitted to a set of z-values, with the windows it was fitted on.</summary>
    public sealed class EmpiricalNullResult
    {
        internal EmpiricalNullResult() { }

        /// <summary>Whether the empirical null was fitted.</summary>
        public EmpiricalNullStatus Status { get; internal set; }
        /// <summary>Scores given.</summary>
        public int Scores { get; internal set; }
        /// <summary>δ0, the centre of the null scores. NaN when not estimated.</summary>
        public double Delta { get; internal set; }
        /// <summary>σ0, the spread of the null scores. NaN when not estimated.</summary>
        public double Sigma { get; internal set; }
        /// <summary>p0, the share of scores that are null: the share inside the final window divided by the fitted
        /// normal's probability of that window (Efron 2007, eq. 4.13). It can exceed 1. NaN when not estimated.</summary>
        public double P0 { get; internal set; }
        /// <summary>The first pass's δ0, fitted on the window around the median. NaN when that pass failed.</summary>
        public double FirstDelta { get; internal set; }
        /// <summary>The first pass's σ0. NaN when that pass failed.</summary>
        public double FirstSigma { get; internal set; }
        /// <summary>The final window's lower bound. NaN when no final window was formed.</summary>
        public double WindowLow { get; internal set; }
        /// <summary>The final window's upper bound. NaN when no final window was formed.</summary>
        public double WindowHigh { get; internal set; }
        /// <summary>Scores inside the final window (0 when none was formed).</summary>
        public int ScoresInWindow { get; internal set; }
    }

    /// <summary>
    /// Efron's empirical null: the normal N(δ0, σ0²) that the bulk of a large set of z-values follows, fitted by maximum
    /// likelihood to the scores in a central window (Efron 2004, J. Am. Stat. Assoc. 99:96; Efron 2007, Ann. Stat.
    /// 35:1351, section 4).
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>The likelihood</b> is Efron's (2007, eq. 4.12): the scores inside a window [A, B] follow the normal truncated to
    /// it, and how many fall inside is binomial. δ0 and σ0 maximise the truncated-normal part; p0 then follows from the
    /// binomial part (eq. 4.13). The truncated normal on a bounded window is an exponential family in (δ0/σ0²,
    /// −1/(2σ0²)), whose log-likelihood is concave, so its maximum is unique when it exists. It is found by Newton's
    /// method with step halving, on scores rescaled to the window [−1, 1], from their mean and standard deviation; the
    /// truncated normal's moments come from the standard recurrence (Johnson, Kotz and Balakrishnan 1994, ch. 13).
    /// </para>
    /// <para>
    /// <b>The window</b> follows the default of R's locfdr 1.1.8 (Efron, Turnbull and Narasimhan), which the papers do not
    /// state: the first pass is centred at the median, with half-width b × (Q3 − Q1) / (2 Φ⁻¹(0.75)), where
    /// b = 4.3 exp(−0.26 log10 N) for N scores (about 1.8 at N = 2,000) and the quartiles are R's default (type 7). The
    /// second pass refits on the window δ0 ± b σ0 of the first. This implementation is written from those descriptions, not from
    /// locfdr's source, and is checked against locfdr's results (Test/StatisticalModels/ReferenceData).
    /// </para>
    /// </remarks>
    public static class EmpiricalNull
    {
        /// <summary>The fewest scores, in all and in each window, for which a null is fitted.</summary>
        public const int MinimumScores = 3;

        /// <summary>The most scores accepted: above this locfdr widens its window, which is not implemented.</summary>
        public const int MaximumScores = 500_000;

        /// <summary>The most Newton steps in each pass.</summary>
        public const int MaxIterations = 200;

        /// <summary>A pass has settled when no parameter moves by more than this, relative to 1 + its size.</summary>
        public const double ConvergenceTolerance = 1e-12;

        /// <summary>Beyond this σ0, in units of the window's half-width, the likelihood is treated as rising without
        /// bound (<see cref="EmpiricalNullStatus.NoMaximum"/>).</summary>
        private const double SpreadBound = 1e4;

        /// <summary>
        /// Fits the empirical null to <paramref name="z"/>.
        /// </summary>
        /// <param name="z">The z-values, all finite.</param>
        /// <returns>The fit; when it fails, its status says why and every estimate is NaN.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="z"/> is null.</exception>
        /// <exception cref="ArgumentException">A z-value is not finite.</exception>
        /// <exception cref="ArgumentOutOfRangeException">More than <see cref="MaximumScores"/> z-values.</exception>
        public static EmpiricalNullResult Fit(IReadOnlyList<double> z)
        {
            ArgumentNullException.ThrowIfNull(z);
            if (z.Count > MaximumScores)
                throw new ArgumentOutOfRangeException(nameof(z), z.Count, $"At most {MaximumScores} z-values are accepted.");
            for (int i = 0; i < z.Count; i++)
                if (!double.IsFinite(z[i]))
                    throw new ArgumentException($"z-value {i} is not finite.", nameof(z));

            var result = new EmpiricalNullResult
            {
                Scores = z.Count, Delta = double.NaN, Sigma = double.NaN, P0 = double.NaN,
                FirstDelta = double.NaN, FirstSigma = double.NaN, WindowLow = double.NaN, WindowHigh = double.NaN,
            };
            if (z.Count < MinimumScores)
            {
                result.Status = EmpiricalNullStatus.TooFewScores;
                return result;
            }

            double[] sorted = z.ToArray();
            Array.Sort(sorted);
            int n = sorted.Length;
            double b = 4.3 * Math.Exp(-0.26 * Math.Log10(n));
            double median = SortedArrayStatistics.QuantileCustom(sorted, 0.5, QuantileDefinition.R7);
            double q1 = SortedArrayStatistics.QuantileCustom(sorted, 0.25, QuantileDefinition.R7);
            double q3 = SortedArrayStatistics.QuantileCustom(sorted, 0.75, QuantileDefinition.R7);
            double scale = (q3 - q1) / (2 * Normal.InvCDF(0, 1, 0.75));

            var (status, d1, s1, _) = FitWindow(sorted, median, b * scale);
            if (status != EmpiricalNullStatus.Estimated)
            {
                result.Status = status;
                return result;
            }
            result.FirstDelta = d1;
            result.FirstSigma = s1;

            double half = b * s1;
            var (status2, d2, s2, inside) = FitWindow(sorted, d1, half);
            result.WindowLow = d1 - half;
            result.WindowHigh = d1 + half;
            result.ScoresInWindow = inside;
            if (status2 != EmpiricalNullStatus.Estimated)
            {
                result.Status = status2;
                return result;
            }

            double h0 = WindowProbability((result.WindowLow - d2) / s2, (result.WindowHigh - d2) / s2);
            result.Status = EmpiricalNullStatus.Estimated;
            result.Delta = d2;
            result.Sigma = s2;
            result.P0 = (double)inside / n / h0;
            return result;
        }

        /// <summary>The truncated-normal maximum-likelihood fit to the sorted scores in [centre − half, centre + half],
        /// bounds included.</summary>
        private static (EmpiricalNullStatus Status, double Delta, double Sigma, int Inside) FitWindow(double[] sorted,
            double centre, double half)
        {
            double low = centre - half, high = centre + half;
            var inside = sorted.Where(v => v >= low && v <= high).ToArray();
            // A window of no width (scores whose quartiles tie) holds no normal with a positive spread.
            if (!(half > 0) || !double.IsFinite(half))
                return (EmpiricalNullStatus.NoMaximum, double.NaN, double.NaN, inside.Length);
            if (inside.Length < MinimumScores)
                return (EmpiricalNullStatus.TooFewScores, double.NaN, double.NaN, inside.Length);

            // On u = (z - centre) / half the window is [-1, 1]; the fit is equivariant, so it is mapped back at the end.
            double m1 = 0, m2 = 0;
            foreach (double v in inside)
            {
                double u = (v - centre) / half;
                m1 += u;
                m2 += u * u;
            }
            m1 /= inside.Length;
            m2 /= inside.Length;
            double variance = m2 - m1 * m1;
            if (!(variance > 0) || !HasNormalMaximum(m1, m2))
                return (EmpiricalNullStatus.NoMaximum, double.NaN, double.NaN, inside.Length);

            // Natural parameters of N(mu, sigma^2) restricted to [-1, 1]: e1 = mu / sigma^2, e2 = -1 / (2 sigma^2).
            double e1 = m1 / variance, e2 = -0.5 / variance;
            double ll = LogLikelihood(e1, e2, m1, m2);
            for (int iteration = 0; iteration < MaxIterations; iteration++)
            {
                var (t1, t2, v11, v12, v22) = Moments(e1, e2);
                double g1 = m1 - t1, g2 = m2 - t2;
                double det = v11 * v22 - v12 * v12;
                if (!(det > 0) || !double.IsFinite(det))
                    return (EmpiricalNullStatus.NotConverged, double.NaN, double.NaN, inside.Length);
                double step1 = (v22 * g1 - v12 * g2) / det, step2 = (v11 * g2 - v12 * g1) / det;

                double t = 1;
                double n1 = e1, n2 = e2, nll = double.NaN;
                for (int halving = 0; halving < 60; halving++, t /= 2)
                {
                    n1 = e1 + t * step1;
                    n2 = e2 + t * step2;
                    if (!(n2 < 0)) continue;
                    nll = LogLikelihood(n1, n2, m1, m2);
                    if (double.IsFinite(nll) && nll >= ll - 1e-15 * Math.Abs(ll)) break;
                }
                if (!(n2 < 0) || !double.IsFinite(nll))
                    return (EmpiricalNullStatus.NotConverged, double.NaN, double.NaN, inside.Length);

                bool settled = Math.Abs(n1 - e1) <= ConvergenceTolerance * (1 + Math.Abs(e1))
                    && Math.Abs(n2 - e2) <= ConvergenceTolerance * (1 + Math.Abs(e2));
                e1 = n1;
                e2 = n2;
                ll = nll;
                double sigmaU = 1 / Math.Sqrt(-2 * e2);
                if (sigmaU > SpreadBound)
                    return (EmpiricalNullStatus.NoMaximum, double.NaN, double.NaN, inside.Length);
                if (settled)
                {
                    double muU = e1 * sigmaU * sigmaU;
                    return (EmpiricalNullStatus.Estimated, centre + half * muU, half * sigmaU, inside.Length);
                }
            }
            return (EmpiricalNullStatus.NotConverged, double.NaN, double.NaN, inside.Length);
        }

        /// <summary>
        /// Whether a normal truncated to [-1, 1] maximises the likelihood of scores with mean <paramref name="m1"/> and
        /// mean square <paramref name="m2"/> there.
        /// </summary>
        /// <remarks>
        /// The family exp(e1 u + e2 u²) on [-1, 1] has a concave log-likelihood over every (e1, e2), and its maximum
        /// exists for scores that are not all equal. It is a normal only when e2 &lt; 0. At e2 = 0 the family is the
        /// exponential exp(λ u), whose best λ matches the mean, E[u] = coth λ − 1/λ = m1. If the log-likelihood still rises
        /// in e2 there (m2 ≥ E[u²] = 1 − 2 E[u] / λ), the maximum has e2 ≥ 0 and no normal attains it: the likelihood over
        /// normals rises without bound as the spread grows.
        /// </remarks>
        private static bool HasNormalMaximum(double m1, double m2)
        {
            if (!(Math.Abs(m1) < 1)) return false;
            // E[u] is increasing in λ, from -1 to 1; bracket the root, then bisect to rounding.
            double lo = -1, hi = 1;
            while (ExponentialMean(lo) > m1) lo *= 2;
            while (ExponentialMean(hi) < m1) hi *= 2;
            for (int i = 0; i < 200 && hi - lo > 1e-15 * (1 + Math.Abs(lo)); i++)
            {
                double mid = 0.5 * (lo + hi);
                if (ExponentialMean(mid) < m1) lo = mid; else hi = mid;
            }
            double lambda = 0.5 * (lo + hi);
            return m2 < ExponentialMeanSquare(lambda);
        }

        /// <summary>E[u] under exp(λ u) on [-1, 1]: coth λ − 1/λ, by its series near 0.</summary>
        private static double ExponentialMean(double lambda)
        {
            if (Math.Abs(lambda) < 1e-3)
            {
                double l2 = lambda * lambda;
                return lambda / 3 * (1 - l2 / 15 * (1 - 2 * l2 / 21));
            }
            return 1 / Math.Tanh(lambda) - 1 / lambda;
        }

        /// <summary>E[u²] under exp(λ u) on [-1, 1]: 1 − 2 E[u] / λ, and 1/3 at λ = 0, by its series near 0.</summary>
        private static double ExponentialMeanSquare(double lambda)
        {
            if (Math.Abs(lambda) < 1e-3)
            {
                double l2 = lambda * lambda;
                return 1.0 / 3 + 2 * l2 / 45 * (1 - 2 * l2 / 21);
            }
            return 1 - 2 * ExponentialMean(lambda) / lambda;
        }

        /// <summary>The mean log-likelihood per score in the window [-1, 1]: e1 m1 + e2 m2 − log ∫ exp(e1 u + e2 u²) du.</summary>
        private static double LogLikelihood(double e1, double e2, double m1, double m2)
        {
            double sigma = 1 / Math.Sqrt(-2 * e2), mu = e1 * sigma * sigma;
            double h0 = WindowProbability((-1 - mu) / sigma, (1 - mu) / sigma);
            // log ∫ exp(e1 u + e2 u²) du over [-1, 1] = mu² / (2 sigma²) + log(sigma √(2π)) + log h0
            double logNormaliser = mu * mu / (2 * sigma * sigma) + Math.Log(sigma * Math.Sqrt(2 * Math.PI)) + Math.Log(h0);
            return e1 * m1 + e2 * m2 - logNormaliser;
        }

        /// <summary>E[u], E[u²] and the covariance of (u, u²) under N(mu, sigma²) truncated to [-1, 1].</summary>
        private static (double T1, double T2, double V11, double V12, double V22) Moments(double e1, double e2)
        {
            double sigma = 1 / Math.Sqrt(-2 * e2), mu = e1 * sigma * sigma;
            double a = (-1 - mu) / sigma, c = (1 - mu) / sigma;
            double h0 = WindowProbability(a, c);
            double fa = Normal.PDF(0, 1, a), fc = Normal.PDF(0, 1, c);
            // Moments of the standard normal truncated to [a, c]: s0 = 1, s1 = (fa - fc) / h0,
            // s_k = (k - 1) s_{k-2} + (a^{k-1} fa - c^{k-1} fc) / h0.
            var s = new double[5];
            s[0] = 1;
            s[1] = (fa - fc) / h0;
            for (int k = 2; k <= 4; k++)
                s[k] = (k - 1) * s[k - 2] + (Math.Pow(a, k - 1) * fa - Math.Pow(c, k - 1) * fc) / h0;
            // Raw moments of u = mu + sigma y.
            var e = new double[5];
            for (int j = 1; j <= 4; j++)
            {
                double sum = 0;
                for (int i = 0; i <= j; i++)
                    sum += Binomial(j, i) * Math.Pow(mu, j - i) * Math.Pow(sigma, i) * s[i];
                e[j] = sum;
            }
            return (e[1], e[2], e[2] - e[1] * e[1], e[3] - e[1] * e[2], e[4] - e[2] * e[2]);
        }

        /// <summary>Φ(c) − Φ(a), taken in the upper tail when both are positive so it keeps its digits.</summary>
        private static double WindowProbability(double a, double c) =>
            a > 0 ? Normal.CDF(0, 1, -a) - Normal.CDF(0, 1, -c) : Normal.CDF(0, 1, c) - Normal.CDF(0, 1, a);

        private static double Binomial(int n, int k) => k switch
        {
            0 => 1,
            _ when k == n => 1,
            _ => Binomial(n - 1, k - 1) + Binomial(n - 1, k),
        };
    }
}
