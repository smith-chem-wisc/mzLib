using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace Statistics
{
    /// <summary>Why a feature does or does not carry a fitted model.</summary>
    public enum FeatureFitStatus
    {
        /// <summary>The model was fitted and has at least one residual degree of freedom.</summary>
        Fitted,
        /// <summary>Fewer observed values than coefficients plus one, so no residual variance exists.</summary>
        TooFewObservations,
        /// <summary>
        /// The design restricted to this feature's observed samples is not of full rank, e.g. every
        /// observed sample shares one sex so the sex coefficient cannot be estimated.
        /// </summary>
        RankDeficient,
        /// <summary>
        /// Logistic regression only: the outcomes are (quasi-)completely separated by the design, so the
        /// maximum-likelihood estimate does not exist (a coefficient runs to ±∞).
        /// </summary>
        Separated,
        /// <summary>An iterative fit did not converge within its iteration limit.</summary>
        NotConverged,
        /// <summary>
        /// Mixed model only: fewer than two groups, or no group observed twice, so the between-group variance
        /// cannot be told apart from the residual variance.
        /// </summary>
        TooFewGroups,
    }

    /// <summary>
    /// Per-feature ordinary least-squares fits of one design against many features. A value that could not
    /// be computed is NaN, and the reason is in <see cref="Status"/>. The result is owned by the fit and is
    /// read-only to callers.
    /// </summary>
    public sealed class LinearModelFit
    {
        internal LinearModelFit(int features, IReadOnlyList<string> coefficientNames)
        {
            CoefficientNames = coefficientNames;
            int p = coefficientNames.Count;
            CoefficientMatrix = Filled(features, p);
            StdevUnscaledMatrix = Filled(features, p);
            SigmaValues = Enumerable.Repeat(double.NaN, features).ToArray();
            DfResidualValues = new int[features];
            ObservedValues = new int[features];
            AverageResponseValues = Enumerable.Repeat(double.NaN, features).ToArray();
            StatusValues = new FeatureFitStatus[features];
        }

        internal double[,] CoefficientMatrix { get; }
        internal double[,] StdevUnscaledMatrix { get; }
        internal double[] SigmaValues { get; }
        internal int[] DfResidualValues { get; }
        internal int[] ObservedValues { get; }
        internal double[] AverageResponseValues { get; }
        internal FeatureFitStatus[] StatusValues { get; }

        /// <summary>Names of the design columns, in order.</summary>
        public IReadOnlyList<string> CoefficientNames { get; }

        /// <summary>Number of features (rows of the response matrix).</summary>
        public int FeatureCount => SigmaValues.Length;

        /// <summary>Least-squares coefficient <paramref name="coefficient"/> of <paramref name="feature"/>.</summary>
        public double Coefficient(int feature, int coefficient) => CoefficientMatrix[feature, coefficient];

        /// <summary>
        /// sqrt of the diagonal element of (XᵀX)⁻¹ for the feature's observed samples: the coefficient's
        /// standard error divided by the residual standard deviation. Moderation replaces the residual
        /// standard deviation and keeps this factor.
        /// </summary>
        public double StdevUnscaled(int feature, int coefficient) => StdevUnscaledMatrix[feature, coefficient];

        /// <summary>Residual standard deviation, sqrt(RSS / residual df), per feature.</summary>
        public IReadOnlyList<double> Sigma => SigmaValues;

        /// <summary>Residual degrees of freedom per feature: observed samples minus coefficients.</summary>
        public IReadOnlyList<int> DfResidual => DfResidualValues;

        /// <summary>Number of samples with a finite response, per feature.</summary>
        public IReadOnlyList<int> Observed => ObservedValues;

        /// <summary>Mean of each feature's observed responses; the covariate an intensity trend is fitted on.</summary>
        public IReadOnlyList<double> AverageResponse => AverageResponseValues;

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
    /// Fits y = Xβ + ε separately for every feature, omitting that feature's missing samples.
    /// </summary>
    /// <remarks>
    /// Missing values are NaN and are OMITTED per feature, never imputed and never read as zero: a
    /// feature observed in 7 of 10 samples is fitted on those 7, against the 7 matching rows of the
    /// design. The fit is by QR decomposition. Output does not depend on the thread count: every
    /// feature writes only its own row.
    /// </remarks>
    public static class LinearModel
    {
        /// <summary>
        /// Relative tolerance below which a diagonal element of R marks the restricted design as rank
        /// deficient.
        /// </summary>
        public const double RankTolerance = 1e-7;

        /// <summary>
        /// Fits one design matrix against every feature (row) of <paramref name="responses"/> by least squares.
        /// </summary>
        /// <param name="responses">
        /// [feature, sample], e.g. log2 intensities. Any non-finite value (NaN, ±Infinity) is treated as
        /// missing and omitted. A 0 is an observation: pass NaN for "not measured", or log-transform first
        /// so that an unmeasured 0 becomes -Infinity. Quantification elsewhere writes missing as 0.
        /// </param>
        /// <param name="design">[sample, coefficient]. Must be finite and of full column rank.</param>
        /// <param name="coefficientNames">One name per design column; defaults to "c0", "c1", ….</param>
        /// <param name="maxThreads">
        /// Degree of parallelism over features. -1 (the default) uses all cores but one, as in FlashLFQ;
        /// any other value below 1 means one thread. Results are identical for every value.
        /// </param>
        public static LinearModelFit Fit(double[,] responses, double[,] design,
            IReadOnlyList<string>? coefficientNames = null, int maxThreads = -1)
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
            coefficientNames ??= Enumerable.Range(0, p).Select(j => $"c{j}").ToList();
            if (coefficientNames.Count != p)
                throw new ArgumentException($"{coefficientNames.Count} coefficient names for {p} design columns.", nameof(coefficientNames));
            if (!IsFullRank(Matrix<double>.Build.DenseOfArray(design)))
                throw new ArgumentException("The design is not of full column rank over all samples; a column is redundant.", nameof(design));

            var fit = new LinearModelFit(features, coefficientNames);
            var options = new ParallelOptions { MaxDegreeOfParallelism = ResolveThreads(maxThreads) };
            Parallel.ForEach(Partitioner.Create(0, Math.Max(features, 1)), options, range =>
            {
                for (int f = range.Item1; f < Math.Min(range.Item2, features); f++)
                    FitOne(responses, design, f, fit);
            });
            return fit;
        }

        /// <summary>FlashLFQ's convention: -1, or at least the core count, means all cores but one; anything else below 1 means one.</summary>
        internal static int ResolveThreads(int maxThreads)
        {
            if (maxThreads == -1 || maxThreads >= Environment.ProcessorCount)
                return Math.Max(1, Environment.ProcessorCount - 1);
            return Math.Max(1, maxThreads);
        }

        private static void FitOne(double[,] responses, double[,] design, int f, LinearModelFit fit)
        {
            int samples = responses.GetLength(1), p = design.GetLength(1);
            var rows = new List<int>(samples);
            double sum = 0;
            for (int s = 0; s < samples; s++)
            {
                double v = responses[f, s];
                if (double.IsFinite(v)) { rows.Add(s); sum += v; }
            }
            int m = rows.Count;
            fit.ObservedValues[f] = m;
            if (m > 0) fit.AverageResponseValues[f] = sum / m;
            if (m <= p)
            {
                fit.StatusValues[f] = FeatureFitStatus.TooFewObservations;
                return;
            }

            var x = Matrix<double>.Build.Dense(m, p, (i, j) => design[rows[i], j]);
            var y = Vector<double>.Build.Dense(m, i => responses[f, rows[i]]);
            var qr = x.QR(MathNet.Numerics.LinearAlgebra.Factorization.QRMethod.Thin);
            if (!IsFullRank(qr.R))
            {
                fit.StatusValues[f] = FeatureFitStatus.RankDeficient;
                return;
            }

            var beta = qr.Solve(y);
            var residual = y - x * beta;
            int df = m - p;
            // (XᵀX)⁻¹ = R⁻¹ R⁻ᵀ, so its diagonal is the row sums of squares of R⁻¹.
            var rInv = qr.R.Inverse();
            for (int j = 0; j < p; j++)
            {
                fit.CoefficientMatrix[f, j] = beta[j];
                double d = 0;
                for (int k = 0; k < p; k++) d += rInv[j, k] * rInv[j, k];
                fit.StdevUnscaledMatrix[f, j] = Math.Sqrt(d);
            }
            fit.SigmaValues[f] = Math.Sqrt(residual.DotProduct(residual) / df);
            fit.DfResidualValues[f] = df;
            fit.StatusValues[f] = FeatureFitStatus.Fitted;
        }

        internal static bool IsFullRank(Matrix<double> m)
        {
            var r = m.RowCount == m.ColumnCount && IsUpperTriangular(m) ? m : m.QR(MathNet.Numerics.LinearAlgebra.Factorization.QRMethod.Thin).R;
            double max = 0;
            for (int j = 0; j < r.ColumnCount; j++) max = Math.Max(max, Math.Abs(r[j, j]));
            if (max == 0) return false;
            for (int j = 0; j < r.ColumnCount; j++)
                if (Math.Abs(r[j, j]) <= RankTolerance * max) return false;
            return true;
        }

        private static bool IsUpperTriangular(Matrix<double> m)
        {
            for (int i = 1; i < m.RowCount; i++)
                for (int j = 0; j < i; j++)
                    if (m[i, j] != 0) return false;
            return true;
        }
    }
}
