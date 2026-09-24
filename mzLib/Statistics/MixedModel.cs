using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;

namespace Statistics
{
    /// <summary>How the variance components of a mixed model are estimated.</summary>
    public enum VarianceEstimator
    {
        /// <summary>
        /// Restricted maximum likelihood: variances estimated from the likelihood of the residuals after the
        /// fixed effects are removed, so they are not biased downward by the number of coefficients. The default.
        /// </summary>
        Reml,
        /// <summary>
        /// Maximum likelihood. Biased variances in small samples, but its likelihoods can be compared between
        /// models with different fixed effects (likelihood-ratio tests), which REML likelihoods cannot.
        /// </summary>
        MaximumLikelihood,
    }

    /// <summary>
    /// Per-feature random-intercept linear mixed models, y = Xβ + u_group + ε. The variance estimator is
    /// named in <see cref="Estimator"/>. Fixed-effect tests are Wald t with degrees of freedom by nlme's
    /// "containment" rule (see <see cref="MixedModel"/>). A value that could not be computed is NaN, and the
    /// reason is in <see cref="Status"/>.
    /// </summary>
    public sealed class MixedModelFit
    {
        internal MixedModelFit(int features, IReadOnlyList<string> coefficientNames, VarianceEstimator estimator)
        {
            CoefficientNames = coefficientNames;
            Estimator = estimator;
            int p = coefficientNames.Count;
            CoefficientMatrix = Filled(features, p);
            StandardErrorMatrix = Filled(features, p);
            DfMatrix = Filled(features, p);
            ResidualVarianceValues = Enumerable.Repeat(double.NaN, features).ToArray();
            GroupVarianceValues = Enumerable.Repeat(double.NaN, features).ToArray();
            LogLikelihoodValues = Enumerable.Repeat(double.NaN, features).ToArray();
            ObservedValues = new int[features];
            GroupValues = new int[features];
            StatusValues = new FeatureFitStatus[features];
        }

        internal double[,] CoefficientMatrix { get; }
        internal double[,] StandardErrorMatrix { get; }
        internal double[,] DfMatrix { get; }
        internal double[] ResidualVarianceValues { get; }
        internal double[] GroupVarianceValues { get; }
        internal double[] LogLikelihoodValues { get; }
        internal int[] ObservedValues { get; }
        internal int[] GroupValues { get; }
        internal FeatureFitStatus[] StatusValues { get; }

        /// <summary>Names of the design columns, in order.</summary>
        public IReadOnlyList<string> CoefficientNames { get; }

        /// <summary>The variance-component estimator every feature was fitted with.</summary>
        public VarianceEstimator Estimator { get; }

        /// <summary>Number of features (rows of the response matrix).</summary>
        public int FeatureCount => ResidualVarianceValues.Length;

        /// <summary>Generalized-least-squares fixed-effect coefficient at the estimated variances.</summary>
        public double Coefficient(int feature, int coefficient) => CoefficientMatrix[feature, coefficient];

        /// <summary>Standard error, √ of the diagonal of σ̂²(XᵀV⁻¹X)⁻¹.</summary>
        public double StandardError(int feature, int coefficient) => StandardErrorMatrix[feature, coefficient];

        /// <summary>Denominator degrees of freedom of the coefficient's t test (containment rule). NaN when not positive.</summary>
        public double DegreesOfFreedom(int feature, int coefficient) => DfMatrix[feature, coefficient];

        /// <summary>Wald statistic, coefficient / standard error.</summary>
        public double TValue(int feature, int coefficient) => CoefficientMatrix[feature, coefficient] / StandardErrorMatrix[feature, coefficient];

        /// <summary>Two-sided p-value of the t test against a coefficient of zero.</summary>
        public double PValue(int feature, int coefficient)
        {
            double t = TValue(feature, coefficient), df = DfMatrix[feature, coefficient];
            return double.IsNaN(t) || double.IsNaN(df) ? double.NaN : 2 * StudentT.CDF(0, 1, df, -Math.Abs(t));
        }

        /// <summary>Residual (within-group) variance σ², per feature.</summary>
        public IReadOnlyList<double> ResidualVariance => ResidualVarianceValues;

        /// <summary>Between-group (random-intercept) variance τ², per feature. Zero when the estimate is on the boundary.</summary>
        public IReadOnlyList<double> GroupVariance => GroupVarianceValues;

        /// <summary>Maximized log-likelihood (restricted when <see cref="Estimator"/> is REML), per feature.</summary>
        public IReadOnlyList<double> LogLikelihood => LogLikelihoodValues;

        /// <summary>Number of samples with a finite response, per feature.</summary>
        public IReadOnlyList<int> Observed => ObservedValues;

        /// <summary>Number of distinct groups among the observed samples, per feature.</summary>
        public IReadOnlyList<int> Groups => GroupValues;

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
    /// Fits a random-intercept linear mixed model, y = Xβ + u_g + ε with u_g ~ N(0, τ²) and ε ~ N(0, σ²),
    /// separately for every feature, omitting that feature's missing samples.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The typical use is one group per dataset (or per donor), so that a covariate's effect is estimated
    /// within groups while groups keep their own baseline.
    /// </para>
    /// <para>
    /// Estimation: σ² is profiled out and the likelihood (REML or ML, <see cref="VarianceEstimator"/>) is
    /// maximized over the variance ratio θ = τ²/σ² ≥ 0, a one-dimensional search, using the closed-form
    /// inverse and determinant of each group's compound-symmetric block. θ = 0 (no between-group variance)
    /// is admitted as a boundary estimate. Fixed effects are then generalized least squares at θ̂.
    /// </para>
    /// <para>
    /// Degrees of freedom follow nlme's containment rule for one grouping level: a coefficient whose column
    /// varies within some group is tested with m − G − p_within df; one that is constant within every group
    /// (a group-level covariate) with G − 1 − p_between df; a column constant over all observed samples
    /// (the intercept) takes the within-group df. Here m is the feature's observed samples and G its
    /// observed groups. The rule is exact for balanced designs and an approximation otherwise (Satterthwaite
    /// or Kenward-Roger are the refinements, not implemented).
    /// </para>
    /// Output does not depend on the thread count.
    /// </remarks>
    public static class MixedModel
    {
        /// <summary>Fits one design and one grouping against every feature (row) of <paramref name="responses"/>.</summary>
        /// <param name="responses">[feature, sample]. Any non-finite value is missing and omitted.</param>
        /// <param name="design">[sample, coefficient]. Must be finite and of full column rank.</param>
        /// <param name="groups">One group label per sample (e.g. the dataset accession). Must be non-empty strings.</param>
        /// <param name="coefficientNames">One name per design column; defaults to "c0", "c1", ….</param>
        /// <param name="estimator">REML (default) or ML.</param>
        /// <param name="maxThreads">-1 (default) uses all cores but one; below 1 means one. Results are identical for every value.</param>
        public static MixedModelFit Fit(double[,] responses, double[,] design, IReadOnlyList<string> groups,
            IReadOnlyList<string>? coefficientNames = null, VarianceEstimator estimator = VarianceEstimator.Reml,
            int maxThreads = -1)
        {
            ArgumentNullException.ThrowIfNull(responses);
            ArgumentNullException.ThrowIfNull(design);
            ArgumentNullException.ThrowIfNull(groups);
            int features = responses.GetLength(0), samples = responses.GetLength(1);
            int n = design.GetLength(0), p = design.GetLength(1);
            if (n != samples)
                throw new ArgumentException($"The design has {n} rows but the response matrix has {samples} samples.", nameof(design));
            if (groups.Count != samples)
                throw new ArgumentException($"{groups.Count} group labels for {samples} samples.", nameof(groups));
            for (int s = 0; s < samples; s++)
                if (string.IsNullOrEmpty(groups[s]))
                    throw new ArgumentException($"Sample {s} has no group label.", nameof(groups));
            if (p == 0)
                throw new ArgumentException("The design has no columns.", nameof(design));
            for (int i = 0; i < n; i++)
                for (int j = 0; j < p; j++)
                    if (!double.IsFinite(design[i, j]))
                        throw new ArgumentException($"The design is not finite at sample {i}, column {j}.", nameof(design));
            coefficientNames ??= Enumerable.Range(0, p).Select(j => $"c{j}").ToList();
            if (coefficientNames.Count != p)
                throw new ArgumentException($"{coefficientNames.Count} coefficient names for {p} design columns.", nameof(coefficientNames));
            if (!LinearModel.IsFullRank(Matrix<double>.Build.DenseOfArray(design)))
                throw new ArgumentException("The design is not of full column rank over all samples; a column is redundant.", nameof(design));

            // Group labels to dense indices in order of first appearance, so the result never depends on hashing.
            var index = new Dictionary<string, int>(StringComparer.Ordinal);
            var groupOf = new int[samples];
            for (int s = 0; s < samples; s++)
            {
                if (!index.TryGetValue(groups[s], out int g)) { g = index.Count; index[groups[s]] = g; }
                groupOf[s] = g;
            }

            var fit = new MixedModelFit(features, coefficientNames, estimator);
            var options = new ParallelOptions { MaxDegreeOfParallelism = LinearModel.ResolveThreads(maxThreads) };
            Parallel.ForEach(Partitioner.Create(0, Math.Max(features, 1)), options, range =>
            {
                for (int f = range.Item1; f < Math.Min(range.Item2, features); f++)
                    FitOne(responses, design, groupOf, index.Count, f, fit, estimator);
            });
            return fit;
        }

        /// <summary>Per-group sufficient statistics; everything the likelihood needs for any θ.</summary>
        private sealed class Blocks
        {
            public required int[] Sizes;                    // n_g
            public required Vector<double>[] ColumnSums;    // s_g = X_gᵀ1
            public required double[] ResponseSums;          // t_g = Σ y_g
            public required Matrix<double> XtX;
            public required Vector<double> Xty;
            public required double Yty;
        }

        private readonly record struct Profile(double Objective, Vector<double> Beta, Matrix<double> AInverse, double Rss);

        private static void FitOne(double[,] responses, double[,] design, int[] groupOf, int groupCount,
            int f, MixedModelFit fit, VarianceEstimator estimator)
        {
            int samples = responses.GetLength(1), p = design.GetLength(1);
            var rows = new List<int>(samples);
            for (int s = 0; s < samples; s++)
                if (double.IsFinite(responses[f, s])) rows.Add(s);
            int m = rows.Count;
            fit.ObservedValues[f] = m;

            // Observed groups, renumbered densely in order of first appearance.
            var local = new int[groupCount];
            Array.Fill(local, -1);
            var rowGroup = new int[m];
            int G = 0;
            for (int i = 0; i < m; i++)
            {
                int g = groupOf[rows[i]];
                if (local[g] < 0) local[g] = G++;
                rowGroup[i] = local[g];
            }
            fit.GroupValues[f] = G;
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
            var sizes = new int[G];
            foreach (int g in rowGroup) sizes[g]++;
            if (G < 2 || sizes.Max() < 2)
            {
                fit.StatusValues[f] = FeatureFitStatus.TooFewGroups;
                return;
            }

            var y = Vector<double>.Build.Dense(m, i => responses[f, rows[i]]);
            var colSums = Enumerable.Range(0, G).Select(_ => Vector<double>.Build.Dense(p)).ToArray();
            var ySums = new double[G];
            for (int i = 0; i < m; i++)
            {
                int g = rowGroup[i];
                for (int j = 0; j < p; j++) colSums[g][j] += x[i, j];
                ySums[g] += y[i];
            }
            var blocks = new Blocks
            {
                Sizes = sizes, ColumnSums = colSums, ResponseSums = ySums,
                XtX = x.TransposeThisAndMultiply(x), Xty = x.TransposeThisAndMultiply(y), Yty = y.DotProduct(y),
            };
            bool reml = estimator == VarianceEstimator.Reml;
            int dfResid = reml ? m - p : m;

            double theta = MaximizeOverTheta(blocks, m, p, reml);
            var best = Evaluate(blocks, theta, m, p, reml);
            if (!double.IsFinite(best.Objective) || !(best.Rss > 0))
            {
                fit.StatusValues[f] = FeatureFitStatus.NotConverged;
                return;
            }

            double sigma2 = best.Rss / dfResid;
            fit.ResidualVarianceValues[f] = sigma2;
            fit.GroupVarianceValues[f] = theta * sigma2;
            fit.LogLikelihoodValues[f] = -0.5 * best.Objective;

            var (dfWithin, dfBetween, isBetween) = ContainmentDf(x, rowGroup, G, m);
            for (int j = 0; j < p; j++)
            {
                fit.CoefficientMatrix[f, j] = best.Beta[j];
                fit.StandardErrorMatrix[f, j] = Math.Sqrt(sigma2 * best.AInverse[j, j]);
                double df = isBetween[j] ? dfBetween : dfWithin;
                fit.DfMatrix[f, j] = df > 0 ? df : double.NaN;
            }
            fit.StatusValues[f] = FeatureFitStatus.Fitted;
        }

        /// <summary>
        /// −2 × profiled log-likelihood at variance ratio θ, with the GLS solution. With V_g = I + θJ,
        /// V_g⁻¹ = I − c_g J where c_g = θ / (1 + n_g θ), and log|V_g| = log(1 + n_g θ).
        /// </summary>
        private static Profile Evaluate(Blocks b, double theta, int m, int p, bool reml)
        {
            var a = b.XtX.Clone();
            var r = b.Xty.Clone();
            double q = b.Yty, logDetV = 0;
            for (int g = 0; g < b.Sizes.Length; g++)
            {
                double c = theta / (1 + b.Sizes[g] * theta);
                logDetV += Math.Log(1 + b.Sizes[g] * theta);
                if (c == 0) continue;
                var s = b.ColumnSums[g];
                a -= c * s.OuterProduct(s);
                r -= c * b.ResponseSums[g] * s;
                q -= c * b.ResponseSums[g] * b.ResponseSums[g];
            }
            var chol = a.Cholesky();
            var beta = chol.Solve(r);
            double rss = q - r.DotProduct(beta);
            int dfResid = reml ? m - p : m;
            if (!(rss > 0)) return new Profile(double.NaN, beta, chol.Solve(Matrix<double>.Build.DenseIdentity(p)), rss);
            double objective = dfResid * Math.Log(2 * Math.PI * rss / dfResid) + logDetV + dfResid;
            if (reml) objective += chol.DeterminantLn;
            return new Profile(objective, beta, chol.Solve(Matrix<double>.Build.DenseIdentity(p)), rss);
        }

        /// <summary>
        /// θ̂ ≥ 0 minimizing the profiled objective: a grid over log θ brackets the minimum, Brent's method
        /// refines it, and the boundary θ = 0 is kept when it is at least as good.
        /// </summary>
        private static double MaximizeOverTheta(Blocks b, int m, int p, bool reml)
        {
            double Obj(double u) => Evaluate(b, Math.Exp(u), m, p, reml).Objective is var o && double.IsFinite(o) ? o : double.MaxValue;

            const double lo = -20, hi = 12, step = 0.5;
            double bestU = lo, bestObj = double.MaxValue;
            for (double u = lo; u <= hi + 1e-12; u += step)
            {
                double o = Obj(u);
                if (o < bestObj) { bestObj = o; bestU = u; }
            }
            double left = Math.Max(lo, bestU - step), right = Math.Min(hi, bestU + step);
            double uStar = Brent(Obj, left, right, 1e-12);
            double thetaStar = Math.Exp(uStar);
            double atZero = Evaluate(b, 0, m, p, reml).Objective;
            double atStar = Evaluate(b, thetaStar, m, p, reml).Objective;
            return double.IsFinite(atZero) && atZero <= atStar ? 0 : thetaStar;
        }

        /// <summary>Brent's (1973) derivative-free minimizer on [a, b].</summary>
        internal static double Brent(Func<double, double> f, double a, double b, double tol)
        {
            const double golden = 0.3819660112501051;
            double x = a + golden * (b - a), w = x, v = x, fx = f(x), fw = fx, fv = fx, d = 0, e = 0;
            for (int iter = 0; iter < 200; iter++)
            {
                double mid = 0.5 * (a + b), tol1 = tol * Math.Abs(x) + 1e-15, tol2 = 2 * tol1;
                if (Math.Abs(x - mid) <= tol2 - 0.5 * (b - a)) break;
                bool golden_step = true;
                if (Math.Abs(e) > tol1)
                {
                    double r = (x - w) * (fx - fv), q = (x - v) * (fx - fw), pp = (x - v) * q - (x - w) * r;
                    q = 2 * (q - r);
                    if (q > 0) pp = -pp; else q = -q;
                    double eTemp = e;
                    e = d;
                    if (Math.Abs(pp) < Math.Abs(0.5 * q * eTemp) && pp > q * (a - x) && pp < q * (b - x))
                    {
                        d = pp / q;
                        double u0 = x + d;
                        if (u0 - a < tol2 || b - u0 < tol2) d = mid >= x ? tol1 : -tol1;
                        golden_step = false;
                    }
                }
                if (golden_step) { e = (x >= mid ? a : b) - x; d = golden * e; }
                double u = Math.Abs(d) >= tol1 ? x + d : x + (d > 0 ? tol1 : -tol1);
                double fu = f(u);
                if (fu <= fx)
                {
                    if (u >= x) a = x; else b = x;
                    v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
                }
                else
                {
                    if (u < x) a = u; else b = u;
                    if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
                    else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
                }
            }
            return x;
        }

        /// <summary>nlme's containment degrees of freedom for one grouping level (see the class remarks).</summary>
        private static (double within, double between, bool[] isBetween) ContainmentDf(Matrix<double> x, int[] rowGroup, int G, int m)
        {
            int p = x.ColumnCount;
            var isBetween = new bool[p];
            int pWithin = 0, pBetween = 0;
            for (int j = 0; j < p; j++)
            {
                bool constantOverall = true, variesWithin = false;
                var first = new double?[G];
                for (int i = 0; i < m; i++)
                {
                    if (x[i, j] != x[0, j]) constantOverall = false;
                    int g = rowGroup[i];
                    if (first[g] is double v0) { if (x[i, j] != v0) variesWithin = true; }
                    else first[g] = x[i, j];
                }
                if (constantOverall) continue;          // intercept: tested at the within-group df, counted at neither level
                if (variesWithin) pWithin++;
                else { pBetween++; isBetween[j] = true; }
            }
            return (m - G - pWithin, G - 1 - pBetween, isBetween);
        }
    }
}
