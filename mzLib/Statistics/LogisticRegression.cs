using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;

namespace Statistics
{
    /// <summary>
    /// Per-feature logistic regressions of one design against many binary features. Estimator: maximum
    /// likelihood by iteratively reweighted least squares (IRLS); test: Wald z. A value that could not be
    /// computed is NaN, and the reason is in <see cref="Status"/>.
    /// </summary>
    public sealed class LogisticModelFit
    {
        internal LogisticModelFit(int features, IReadOnlyList<string> coefficientNames)
        {
            CoefficientNames = coefficientNames;
            int p = coefficientNames.Count;
            CoefficientMatrix = Filled(features, p);
            StandardErrorMatrix = Filled(features, p);
            DevianceValues = Enumerable.Repeat(double.NaN, features).ToArray();
            IterationValues = new int[features];
            ObservedValues = new int[features];
            EventValues = new int[features];
            StatusValues = new FeatureFitStatus[features];
        }

        internal double[,] CoefficientMatrix { get; }
        internal double[,] StandardErrorMatrix { get; }
        internal double[] DevianceValues { get; }
        internal int[] IterationValues { get; }
        internal int[] ObservedValues { get; }
        internal int[] EventValues { get; }
        internal FeatureFitStatus[] StatusValues { get; }

        /// <summary>Names of the design columns, in order.</summary>
        public IReadOnlyList<string> CoefficientNames { get; }

        /// <summary>Number of features (rows of the response matrix).</summary>
        public int FeatureCount => DevianceValues.Length;

        /// <summary>Maximum-likelihood coefficient on the log-odds scale.</summary>
        public double Coefficient(int feature, int coefficient) => CoefficientMatrix[feature, coefficient];

        /// <summary>Standard error from the inverse Fisher information at the estimate.</summary>
        public double StandardError(int feature, int coefficient) => StandardErrorMatrix[feature, coefficient];

        /// <summary>Wald statistic, coefficient / standard error.</summary>
        public double ZValue(int feature, int coefficient) => CoefficientMatrix[feature, coefficient] / StandardErrorMatrix[feature, coefficient];

        /// <summary>Two-sided Wald p-value against a coefficient of zero (standard normal reference).</summary>
        public double PValue(int feature, int coefficient)
        {
            double z = ZValue(feature, coefficient);
            return double.IsNaN(z) ? double.NaN : 2 * Normal.CDF(0, 1, -Math.Abs(z));
        }

        /// <summary>Residual deviance, −2 × log-likelihood, per feature.</summary>
        public IReadOnlyList<double> Deviance => DevianceValues;

        /// <summary>IRLS iterations used, per feature (0 when no fit was attempted).</summary>
        public IReadOnlyList<int> Iterations => IterationValues;

        /// <summary>Number of samples with an observed (0 or 1) response, per feature.</summary>
        public IReadOnlyList<int> Observed => ObservedValues;

        /// <summary>Number of observed responses equal to 1, per feature.</summary>
        public IReadOnlyList<int> Events => EventValues;

        /// <summary>Whether each feature was fitted, and if not, why.</summary>
        public IReadOnlyList<FeatureFitStatus> Status => StatusValues;

        /// <summary>Index of a coefficient by name, or an exception naming the ones that exist.</summary>
        public int IndexOf(string coefficientName)
        {
            for (int i = 0; i < CoefficientNames.Count; i++)
                if (CoefficientNames[i] == coefficientName) return i;
            throw new ArgumentException(
                $"No coefficient named '{coefficientName}'. The design has: {string.Join(", ", CoefficientNames)}.",
                nameof(coefficientName));
        }

        private static double[,] Filled(int rows, int cols)
        {
            var a = new double[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    a[i, j] = double.NaN;
            return a;
        }
    }

    /// <summary>
    /// Fits logit P(y = 1) = Xβ separately for every feature, omitting that feature's missing samples.
    /// </summary>
    /// <remarks>
    /// Estimator: maximum likelihood by IRLS (Fisher scoring, which for the canonical logit link is
    /// Newton-Raphson), started as R's glm does and iterated until the relative change in deviance is
    /// below <c>tolerance</c>. Standard errors come from the inverse Fisher information and tests are
    /// Wald z. When the outcomes are separated by the design the estimate does not exist; the feature is
    /// reported <see cref="FeatureFitStatus.Separated"/> with NaN values rather than the arbitrarily large
    /// coefficients an unchecked IRLS returns. Output does not depend on the thread count.
    /// </remarks>
    public static class LogisticRegression
    {
        /// <summary>Fitted probabilities closer than this to 0 or 1 mark a separated fit.</summary>
        public const double SeparationThreshold = 1e-10;

        /// <summary>Fits one design matrix against every feature (row) of <paramref name="responses"/>.</summary>
        /// <param name="responses">
        /// [feature, sample]. Each value is 0, 1, or non-finite (missing, omitted). Any other value is an
        /// error: a detection indicator or other binary outcome is expected, not a proportion.
        /// </param>
        /// <param name="design">[sample, coefficient]. Must be finite and of full column rank.</param>
        /// <param name="coefficientNames">One name per design column; defaults to "c0", "c1", ….</param>
        /// <param name="maxThreads">-1 (default) uses all cores but one; below 1 means one. Results are identical for every value.</param>
        /// <param name="maxIterations">IRLS iteration limit per feature.</param>
        /// <param name="tolerance">Convergence when |dev − dev_prev| / (|dev| + 0.1) is below this (R glm's rule).</param>
        public static LogisticModelFit Fit(double[,] responses, double[,] design,
            IReadOnlyList<string>? coefficientNames = null, int maxThreads = -1,
            int maxIterations = 100, double tolerance = 1e-12)
        {
            ArgumentNullException.ThrowIfNull(responses);
            ArgumentNullException.ThrowIfNull(design);
            int features = responses.GetLength(0), samples = responses.GetLength(1);
            int n = design.GetLength(0), p = design.GetLength(1);
            if (n != samples)
                throw new ArgumentException($"The design has {n} rows but the response matrix has {samples} samples.", nameof(design));
            if (p == 0)
                throw new ArgumentException("The design has no columns.", nameof(design));
            for (int i = 0; i < n; i++)
                for (int j = 0; j < p; j++)
                    if (!double.IsFinite(design[i, j]))
                        throw new ArgumentException($"The design is not finite at sample {i}, column {j}.", nameof(design));
            for (int f = 0; f < features; f++)
                for (int s = 0; s < samples; s++)
                {
                    double v = responses[f, s];
                    if (double.IsFinite(v) && v != 0 && v != 1)
                        throw new ArgumentException($"Response at feature {f}, sample {s} is {v}; expected 0, 1 or missing.", nameof(responses));
                }
            coefficientNames ??= Enumerable.Range(0, p).Select(j => $"c{j}").ToList();
            if (coefficientNames.Count != p)
                throw new ArgumentException($"{coefficientNames.Count} coefficient names for {p} design columns.", nameof(coefficientNames));
            if (!LinearModel.IsFullRank(Matrix<double>.Build.DenseOfArray(design)))
                throw new ArgumentException("The design is not of full column rank over all samples; a column is redundant.", nameof(design));
            if (maxIterations < 1) throw new ArgumentOutOfRangeException(nameof(maxIterations));
            if (!(tolerance > 0)) throw new ArgumentOutOfRangeException(nameof(tolerance));

            var fit = new LogisticModelFit(features, coefficientNames);
            var options = new ParallelOptions { MaxDegreeOfParallelism = LinearModel.ResolveThreads(maxThreads) };
            Parallel.ForEach(Partitioner.Create(0, Math.Max(features, 1)), options, range =>
            {
                for (int f = range.Item1; f < Math.Min(range.Item2, features); f++)
                    FitOne(responses, design, f, fit, maxIterations, tolerance);
            });
            return fit;
        }

        private static void FitOne(double[,] responses, double[,] design, int f, LogisticModelFit fit,
            int maxIterations, double tolerance)
        {
            int samples = responses.GetLength(1), p = design.GetLength(1);
            var rows = new List<int>(samples);
            int events = 0;
            for (int s = 0; s < samples; s++)
            {
                double v = responses[f, s];
                if (double.IsFinite(v)) { rows.Add(s); if (v == 1) events++; }
            }
            int m = rows.Count;
            fit.ObservedValues[f] = m;
            fit.EventValues[f] = events;
            if (m <= p)
            {
                fit.StatusValues[f] = FeatureFitStatus.TooFewObservations;
                return;
            }
            var x = Matrix<double>.Build.Dense(m, p, (i, j) => design[rows[i], j]);
            if (!LinearModel.IsFullRank(x))
            {
                fit.StatusValues[f] = FeatureFitStatus.RankDeficient;
                return;
            }
            if (events == 0 || events == m)
            {
                // Every outcome the same: the intercept runs to ±∞, which is separation in its simplest form.
                fit.StatusValues[f] = FeatureFitStatus.Separated;
                return;
            }

            var y = Vector<double>.Build.Dense(m, i => responses[f, rows[i]]);
            var mu = y.Map(v => (v + 0.5) / 2);
            var eta = mu.Map(Logit);
            var beta = Vector<double>.Build.Dense(p);
            double dev = Deviance(y, mu), devOld;
            bool converged = false;
            int iter = 0;
            while (iter < maxIterations)
            {
                iter++;
                var w = mu.Map(u => u * (1 - u));
                var z = eta + (y - mu).PointwiseDivide(w);
                var sw = w.PointwiseSqrt();
                var xw = Matrix<double>.Build.Dense(m, p, (i, j) => x[i, j] * sw[i]);
                var zw = z.PointwiseMultiply(sw);
                beta = xw.QR(MathNet.Numerics.LinearAlgebra.Factorization.QRMethod.Thin).Solve(zw);
                eta = x * beta;
                mu = eta.Map(Inverse);
                devOld = dev;
                dev = Deviance(y, mu);
                if (Math.Abs(dev - devOld) / (Math.Abs(dev) + 0.1) < tolerance) { converged = true; break; }
            }
            fit.IterationValues[f] = iter;

            bool extreme = mu.Any(u => u < SeparationThreshold || u > 1 - SeparationThreshold);
            if (extreme)
            {
                fit.StatusValues[f] = FeatureFitStatus.Separated;
                return;
            }
            if (!converged)
            {
                fit.StatusValues[f] = FeatureFitStatus.NotConverged;
                return;
            }

            // Var(β̂) = (XᵀWX)⁻¹ at the estimate = R⁻¹R⁻ᵀ from the QR of W^½X.
            var sqrtW = mu.Map(u => Math.Sqrt(u * (1 - u)));
            var xwFinal = Matrix<double>.Build.Dense(m, p, (i, j) => x[i, j] * sqrtW[i]);
            var rInv = xwFinal.QR(MathNet.Numerics.LinearAlgebra.Factorization.QRMethod.Thin).R.Inverse();
            for (int j = 0; j < p; j++)
            {
                fit.CoefficientMatrix[f, j] = beta[j];
                double d = 0;
                for (int k = 0; k < p; k++) d += rInv[j, k] * rInv[j, k];
                fit.StandardErrorMatrix[f, j] = Math.Sqrt(d);
            }
            fit.DevianceValues[f] = dev;
            fit.StatusValues[f] = FeatureFitStatus.Fitted;
        }

        private static double Logit(double u) => Math.Log(u / (1 - u));

        /// <summary>Inverse logit, clamped as R's binomial()$linkinv is so that 0 &lt; μ &lt; 1 in floating point.</summary>
        private static double Inverse(double eta)
        {
            const double thresh = 30, eps = 2.220446049250313e-16;
            double e = eta < -thresh ? eps : eta > thresh ? 1 / eps : Math.Exp(eta);
            return e / (1 + e);
        }

        private static double Deviance(Vector<double> y, Vector<double> mu)
        {
            double d = 0;
            for (int i = 0; i < y.Count; i++)
                d -= 2 * (y[i] == 1 ? Math.Log(mu[i]) : Math.Log(1 - mu[i]));
            return d;
        }
    }
}
