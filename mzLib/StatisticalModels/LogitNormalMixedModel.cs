using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Optimization;

namespace StatisticalModels
{
    /// <summary>A count cell: of <see cref="Trials"/> chances, <see cref="Successes"/> succeeded.</summary>
    public readonly record struct BinomialCount(int Successes, int Trials);

    /// <summary>
    /// How the latent value of a <see cref="LogitNormalMixedModel"/> sets the chance g that one trial of a count cell
    /// succeeds: logit g = <see cref="Intercept"/> + <see cref="Slope"/> × latent + <see cref="DepthSlope"/> × depth,
    /// on natural log-odds. Trials within one cell are correlated by <see cref="Overdispersion"/> ρ (a beta-binomial
    /// with mean g and intra-cell correlation ρ); ρ = 0 is the binomial.
    /// </summary>
    /// <remarks>
    /// The identity link (0, 1) makes g the latent's own inverse logit: a count cell then samples the same
    /// probability a measured cell reads. A slope below 1 says the chance of a success rises more slowly than the
    /// latent, e.g. when detecting a form depends on how much of it there is rather than on its share.
    /// </remarks>
    public sealed record DetectionLink(double Intercept, double Slope, double DepthSlope = 0, double Overdispersion = 0)
    {
        /// <summary>g = the latent's inverse logit, binomial trials.</summary>
        public static DetectionLink Identity { get; } = new(0, 1);
    }

    /// <summary>
    /// One cell for <see cref="LogitNormalMixedModel.FitLink"/>: its count, and the latent value it is taken to sit at,
    /// as a normal with mean <see cref="LatentMean"/> and standard deviation <see cref="LatentSd"/> (e.g. the mean and
    /// SD of the feature's measured cells).
    /// </summary>
    public readonly record struct LinkCell(double LatentMean, double LatentSd, int Successes, int Trials, double Depth = 0);

    /// <summary>
    /// Per-feature fits of <see cref="LogitNormalMixedModel"/>. Fixed-effect tests are Wald t with nlme's containment
    /// degrees of freedom over every sample with a response of either kind (see <see cref="MixedModel"/>). A value that
    /// could not be computed is NaN, and the reason is in <see cref="Status"/>.
    /// </summary>
    public sealed class LogitNormalMixedModelFit
    {
        internal LogitNormalMixedModelFit(int features, IReadOnlyList<string> coefficientNames, VarianceEstimator estimator, DetectionLink link)
        {
            CoefficientNames = coefficientNames;
            Estimator = estimator;
            Link = link;
            int p = coefficientNames.Count;
            CoefficientMatrix = Filled(features, p);
            StandardErrorMatrix = Filled(features, p);
            DfMatrix = Filled(features, p);
            ResidualVarianceValues = Enumerable.Repeat(double.NaN, features).ToArray();
            GroupVarianceValues = Enumerable.Repeat(double.NaN, features).ToArray();
            LogLikelihoodValues = Enumerable.Repeat(double.NaN, features).ToArray();
            MeasuredValues = new int[features];
            CountedValues = new int[features];
            GroupValues = new int[features];
            ResidualVarianceGivenValues = new bool[features];
            StatusValues = new FeatureFitStatus[features];
        }

        internal double[,] CoefficientMatrix { get; }
        internal double[,] StandardErrorMatrix { get; }
        internal double[,] DfMatrix { get; }
        internal double[] ResidualVarianceValues { get; }
        internal double[] GroupVarianceValues { get; }
        internal double[] LogLikelihoodValues { get; }
        internal int[] MeasuredValues { get; }
        internal int[] CountedValues { get; }
        internal int[] GroupValues { get; }
        internal bool[] ResidualVarianceGivenValues { get; }
        internal FeatureFitStatus[] StatusValues { get; }

        /// <summary>Names of the design columns, in order.</summary>
        public IReadOnlyList<string> CoefficientNames { get; }

        /// <summary>The variance-component estimator every feature was fitted with.</summary>
        public VarianceEstimator Estimator { get; }

        /// <summary>The link every count cell was read through.</summary>
        public DetectionLink Link { get; }

        /// <summary>Number of features (rows of the response matrix).</summary>
        public int FeatureCount => ResidualVarianceValues.Length;

        /// <summary>Fixed-effect coefficient on the latent scale, at the estimated variances.</summary>
        public double Coefficient(int feature, int coefficient) => CoefficientMatrix[feature, coefficient];

        /// <summary>
        /// Standard error: √ of the diagonal of the inverse of −∂²ℓ/∂β² at the estimate, ℓ the log-likelihood with the
        /// random effects integrated out. With no count cells this is (XᵀV⁻¹X)⁻¹, as <see cref="MixedModel"/> reports
        /// under REML.
        /// </summary>
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

        /// <summary>Residual (run-to-run) variance σ² of the latent, per feature: estimated, or the value given.</summary>
        public IReadOnlyList<double> ResidualVariance => ResidualVarianceValues;

        /// <summary>Whether σ² was given by the caller rather than estimated, per feature.</summary>
        public IReadOnlyList<bool> ResidualVarianceGiven => ResidualVarianceGivenValues;

        /// <summary>Between-group (random-intercept) variance τ², per feature.</summary>
        public IReadOnlyList<double> GroupVariance => GroupVarianceValues;

        /// <summary>
        /// Maximized log-likelihood, per feature: under REML the Laplace-approximated restricted likelihood (the fixed
        /// effects integrated out), exactly the restricted likelihood when there are no count cells.
        /// </summary>
        public IReadOnlyList<double> LogLikelihood => LogLikelihoodValues;

        /// <summary>Number of samples with a measured (Gaussian) response, per feature.</summary>
        public IReadOnlyList<int> Measured => MeasuredValues;

        /// <summary>Number of samples with a count, per feature.</summary>
        public IReadOnlyList<int> Counted => CountedValues;

        /// <summary>Number of distinct groups among the samples with a response of either kind, per feature.</summary>
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
    /// A random-intercept mixed model on a latent scale (e.g. log-odds of a site's occupancy) where each sample is
    /// either a measurement of the latent or a count governed by it, fitted separately for every feature.
    /// </summary>
    /// <remarks>
    /// <para>
    /// For sample s in group g: the latent is η_s = x_sᵀβ + u_g + e_s, with u_g ~ N(0, τ²) and e_s ~ N(0, σ²). A measured
    /// sample reads y_s = η_s. A counted sample has k_s of n_s trials succeed, each with chance g_s, where
    /// logit g_s = a + b η_s + c d_s (<see cref="DetectionLink"/>; d_s the sample's depth), binomially or, with
    /// overdispersion ρ, beta-binomially. A sample with neither is missing. With no counted samples the model is
    /// <see cref="MixedModel"/>'s.
    /// </para>
    /// <para>
    /// Likelihood: per group, the measured samples are integrated against u exactly (they are Gaussian), which leaves
    /// u normal with the measured samples' posterior mean and variance. The counted samples are then integrated over
    /// that normal by 40-node adaptive Gauss-Hermite quadrature (as glmer's nAGQ), each one's e_s by 12 nodes. The u rule
    /// sits at the mode of the group's integrand over (u, e) and the curvature of u's Laplace marginal there; each e rule,
    /// per u node, at the mode and curvature of that count's integrand. A rule centred on the prior instead fails when
    /// a group's counts say much more about u than τ does, or a count says more about its e than σ does. The rules are
    /// held fixed while β is optimized, so ℓ's derivatives are exact for them, then re-centred at the optimum until
    /// they stop moving. A count's integrand over e is log-concave and close to normal, so 12 nodes suffice there; u keeps
    /// 40 for a group whose integrand is one-sided (every count 0, or every count full) at a large τ. Against lme4's exact
    /// likelihood the fits agree to about 1e-6 (ReferenceComparisonLogitNormalTests).
    /// </para>
    /// <para>
    /// Estimation: for given variances, β maximizes the log-likelihood ℓ (Newton, with ℓ's gradient and Hessian exact
    /// for the quadrature). REML (default) maximizes over the variances
    /// ℓ(β̂) − ½ log|H| + (p/2) log 2π, H = −∂²ℓ/∂β² at β̂: the fixed effects integrated out by a Laplace approximation,
    /// as glmmTMB's REML does. It is exact when no sample is counted, so it is REML then. ML maximizes ℓ(β̂). The
    /// variances are searched on log τ and log σ by Nelder-Mead, or on log τ alone by Brent's method when σ² is given.
    /// Standard errors are √diag H⁻¹ at the estimate. A coefficient beyond ±<see cref="SeparationLimit"/>, or an H that
    /// is not positive definite, is <see cref="FeatureFitStatus.Separated"/>: the counts are (quasi-)completely
    /// separated by the design and the estimate runs to infinity.
    /// </para>
    /// Output does not depend on the thread count.
    /// </remarks>
    public static class LogitNormalMixedModel
    {
        /// <summary>A coefficient beyond this on the latent scale is taken as running to infinity (separation).</summary>
        public const double SeparationLimit = 30;

        private const int UNodes = 40, ENodes = 12;
        private const double LogSdMin = -8, LogSdMax = 5;
        private static readonly (double[] X, double[] W) UQuadrature = GaussHermite(UNodes);
        private static readonly (double[] X, double[] W) EQuadrature = GaussHermite(ENodes);

        /// <summary>log w + x², the log weights of a rule for ∫ f(x) dx rather than ∫ f(x) e^(−x²) dx.</summary>
        private static readonly double[] ULogWeights = UQuadrature.X.Select((x, i) => Math.Log(UQuadrature.W[i]) + x * x).ToArray();
        private static readonly double[] ELogWeights = EQuadrature.X.Select((x, i) => Math.Log(EQuadrature.W[i]) + x * x).ToArray();

        /// <summary>Fits one design, grouping and link against every feature (row).</summary>
        /// <param name="responses">[feature, sample] measured latent values; any non-finite value is not measured.</param>
        /// <param name="counts">[feature, sample] counts, or null where not counted. A sample is measured or counted, not both.</param>
        /// <param name="design">[sample, coefficient]. Must be finite and of full column rank.</param>
        /// <param name="groups">One group label per sample. Must be non-empty strings.</param>
        /// <param name="link">How the latent sets the chance per trial.</param>
        /// <param name="depth">One value per sample for the link's depth term. Needed when <see cref="DetectionLink.DepthSlope"/> is not 0.</param>
        /// <param name="residualVariance">
        /// Per feature, a σ² to use instead of estimating it (e.g. borrowed when a feature has too few measured samples),
        /// or NaN to estimate. Null (default): estimate every one.
        /// </param>
        /// <param name="coefficientNames">One name per design column; defaults to "c0", "c1", ….</param>
        /// <param name="estimator">REML (default) or ML.</param>
        /// <param name="maxThreads">-1 (default) uses all cores but one; below 1 means one. Results are identical for every value.</param>
        public static LogitNormalMixedModelFit Fit(double[,] responses, BinomialCount?[,] counts, double[,] design,
            IReadOnlyList<string> groups, DetectionLink link, IReadOnlyList<double>? depth = null,
            IReadOnlyList<double>? residualVariance = null, IReadOnlyList<string>? coefficientNames = null,
            VarianceEstimator estimator = VarianceEstimator.Reml, int maxThreads = -1)
        {
            ArgumentNullException.ThrowIfNull(responses);
            ArgumentNullException.ThrowIfNull(counts);
            ArgumentNullException.ThrowIfNull(design);
            ArgumentNullException.ThrowIfNull(groups);
            ArgumentNullException.ThrowIfNull(link);
            int features = responses.GetLength(0), samples = responses.GetLength(1);
            int n = design.GetLength(0), p = design.GetLength(1);
            if (counts.GetLength(0) != features || counts.GetLength(1) != samples)
                throw new ArgumentException($"The counts are {counts.GetLength(0)} × {counts.GetLength(1)}; the responses are {features} × {samples}.", nameof(counts));
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
            if (!double.IsFinite(link.Intercept) || !double.IsFinite(link.Slope) || link.Slope == 0 || !double.IsFinite(link.DepthSlope)
                || !(link.Overdispersion >= 0 && link.Overdispersion < 1))
                throw new ArgumentException("The link needs a finite intercept and depth slope, a finite non-zero slope, and overdispersion in [0, 1).", nameof(link));
            if (link.DepthSlope != 0 && depth == null)
                throw new ArgumentException("The link has a depth term; supply one depth per sample.", nameof(depth));
            if (depth != null && (depth.Count != samples || depth.Any(d => !double.IsFinite(d))))
                throw new ArgumentException($"Supply one finite depth per sample ({samples}).", nameof(depth));
            if (residualVariance != null && (residualVariance.Count != features || residualVariance.Any(v => !double.IsNaN(v) && !(v > 0 && double.IsFinite(v)))))
                throw new ArgumentException($"Supply one σ² per feature ({features}), positive and finite, or NaN to estimate it.", nameof(residualVariance));
            for (int f = 0; f < features; f++)
                for (int s = 0; s < samples; s++)
                {
                    if (counts[f, s] is not BinomialCount c) continue;
                    if (double.IsFinite(responses[f, s]))
                        throw new ArgumentException($"Feature {f}, sample {s} is both measured and counted.", nameof(counts));
                    if (c.Trials < 1 || c.Successes < 0 || c.Successes > c.Trials)
                        throw new ArgumentException($"Feature {f}, sample {s}: {c.Successes} of {c.Trials} is not a count.", nameof(counts));
                }
            coefficientNames ??= Enumerable.Range(0, p).Select(j => $"c{j}").ToList();
            if (coefficientNames.Count != p)
                throw new ArgumentException($"{coefficientNames.Count} coefficient names for {p} design columns.", nameof(coefficientNames));
            if (!LinearModel.IsFullRank(Matrix<double>.Build.DenseOfArray(design)))
                throw new ArgumentException("The design is not of full column rank over all samples; a column is redundant.", nameof(design));

            var index = new Dictionary<string, int>(StringComparer.Ordinal);
            var groupOf = new int[samples];
            for (int s = 0; s < samples; s++)
            {
                if (!index.TryGetValue(groups[s], out int g)) { g = index.Count; index[groups[s]] = g; }
                groupOf[s] = g;
            }

            var fit = new LogitNormalMixedModelFit(features, coefficientNames, estimator, link);
            var options = new ParallelOptions { MaxDegreeOfParallelism = LinearModel.ResolveThreads(maxThreads) };
            Parallel.ForEach(Partitioner.Create(0, Math.Max(features, 1)), options, range =>
            {
                for (int f = range.Item1; f < Math.Min(range.Item2, features); f++)
                    FitOne(responses, counts, design, groupOf, index.Count, link, depth, residualVariance?[f] ?? double.NaN, f, fit, estimator);
            });
            return fit;
        }

        /// <summary>
        /// Fits a <see cref="DetectionLink"/> by maximum likelihood over <paramref name="cells"/>, each cell's probability
        /// integrated over its latent normal (12-node adaptive Gauss-Hermite), summed over cells.
        /// </summary>
        /// <param name="cells">Covered cells with their latent mean and SD. Each needs Trials ≥ 1, finite values, SD ≥ 0.</param>
        /// <param name="withDepth">Fit the depth slope (else 0).</param>
        /// <param name="withOverdispersion">Fit the beta-binomial ρ (else 0, binomial).</param>
        /// <param name="detectionOnly">
        /// True (default): fit on whether any trial succeeded, P(k = 0) against P(k ≥ 1), the event that makes a
        /// covered-zero. False: fit on the full count, P(k of n). On PTM data the two differ (results/g46_trials_v1 in
        /// the ptmQtl project): among cells where the modified form was seen, how many of the covering peptidoforms carry
        /// it follows occupancy more weakly than whether it was seen at all, and a link fitted on the full counts
        /// predicts covered-zeros worse.
        /// </param>
        /// <returns>The link, or null when the search did not converge.</returns>
        public static DetectionLink? FitLink(IReadOnlyList<LinkCell> cells, bool withDepth = false, bool withOverdispersion = false,
            bool detectionOnly = true)
        {
            ArgumentNullException.ThrowIfNull(cells);
            if (cells.Count == 0) throw new ArgumentException("No cells.", nameof(cells));
            foreach (var c in cells)
                if (!double.IsFinite(c.LatentMean) || !(c.LatentSd >= 0) || !double.IsFinite(c.LatentSd) || !double.IsFinite(c.Depth)
                    || c.Trials < 1 || c.Successes < 0 || c.Successes > c.Trials)
                    throw new ArgumentException("Every cell needs a finite latent mean, SD ≥ 0, finite depth, and 0 ≤ successes ≤ trials ≥ 1.", nameof(cells));

            DetectionLink LinkOf(Vector<double> v)
            {
                int i = 2;
                double c = withDepth ? v[i++] : 0;
                double rho = withOverdispersion ? 1 / (1 + Math.Exp(-v[i])) : 0;
                return new DetectionLink(v[0], v[1], c, rho);
            }
            double Objective(Vector<double> v)
            {
                var link = LinkOf(v);
                double sum = 0;
                foreach (var c in cells)
                {
                    // With detectionOnly, the probability of no success; its complement is used when there was one.
                    int k = detectionOnly ? 0 : c.Successes;
                    double logChoose = LogChoose(c.Trials, k);
                    double offset = link.Intercept + link.Slope * c.LatentMean + link.DepthSlope * c.Depth;
                    double lse;
                    if (c.LatentSd == 0)
                        lse = LogCountProbability(k, c.Trials, offset, link.Overdispersion);
                    else
                    {
                        var (centre, scale) = CellRule(k, c.Trials, logChoose, offset, link.Slope, link.Overdispersion, c.LatentSd, 0);
                        lse = LogCellIntegral(k, c.Trials, logChoose, offset, link.Slope, link.Overdispersion, c.LatentSd, centre, scale);
                    }
                    if (detectionOnly && c.Successes > 0) lse = Math.Log(-double.ExpM1(lse));
                    sum += lse;
                }
                return double.IsFinite(sum) ? -sum : double.MaxValue;
            }
            int dims = 2 + (withDepth ? 1 : 0) + (withOverdispersion ? 1 : 0);
            var start = Vector<double>.Build.Dense(dims);
            start[1] = 1;
            if (withOverdispersion) start[dims - 1] = -2;
            var step = Vector<double>.Build.Dense(dims, 0.3);
            try
            {
                var result = NelderMeadSimplex.Minimum(ObjectiveFunction.Value(Objective), start, step, 1e-10, 20000);
                return LinkOf(result.MinimizingPoint);
            }
            catch (MaximumIterationsException)
            {
                return null;
            }
        }

        /// <summary>
        /// Gauss-Hermite nodes and weights for ∫ f(x) e^(−x²) dx, by Golub-Welsch: the nodes are the eigenvalues of the
        /// Jacobi matrix with off-diagonal √(k/2), the weights √π times the squared first eigenvector components.
        /// </summary>
        internal static (double[] X, double[] W) GaussHermite(int n)
        {
            var j = Matrix<double>.Build.Dense(n, n);
            for (int k = 1; k < n; k++) j[k, k - 1] = j[k - 1, k] = Math.Sqrt(k / 2.0);
            var evd = j.Evd(Symmetricity.Symmetric);
            var x = new double[n];
            var w = new double[n];
            var order = Enumerable.Range(0, n).OrderBy(i => evd.EigenValues[i].Real).ToArray();
            for (int i = 0; i < n; i++)
            {
                int o = order[i];
                x[i] = evd.EigenValues[o].Real;
                double v0 = evd.EigenVectors[0, o];
                w[i] = Math.Sqrt(Math.PI) * v0 * v0;
            }
            // Symmetrize: the rule is exactly symmetric, and averaging removes the eigensolver's last-digit asymmetry.
            for (int i = 0; i < n / 2; i++)
            {
                double xs = 0.5 * (x[n - 1 - i] - x[i]), ws = 0.5 * (w[i] + w[n - 1 - i]);
                x[i] = -xs; x[n - 1 - i] = xs; w[i] = w[n - 1 - i] = ws;
            }
            if (n % 2 == 1) x[n / 2] = 0;
            return (x, w);
        }

        /// <summary>One feature's samples, in the order of the observed rows.</summary>
        private sealed class Data
        {
            public required Matrix<double> X;          // [rows, p]
            public required double[][] Rows;           // X's rows
            public required double[] LogChoose;        // log C(N, K) per counted row
            public required int[] Group;               // dense, 0..G-1
            public required int G;
            public required bool[] IsMeasured;
            public required double[] Y;                // measured value, or NaN
            public required int[] K, N;                // count, or 0/0
            public required double[] Depth;
            public required DetectionLink Link;
            public required List<int>[] MeasuredRows, CountedRows;   // per group
        }

        private static void FitOne(double[,] responses, BinomialCount?[,] counts, double[,] design, int[] groupOf, int groupCount,
            DetectionLink link, IReadOnlyList<double>? depth, double givenSigma2, int f, LogitNormalMixedModelFit fit, VarianceEstimator estimator)
        {
            int samples = responses.GetLength(1), p = design.GetLength(1);
            var rows = new List<int>(samples);
            for (int s = 0; s < samples; s++)
                if (double.IsFinite(responses[f, s]) || counts[f, s] is not null) rows.Add(s);
            int m = rows.Count;
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
            var measured = rows.Select(s => double.IsFinite(responses[f, s])).ToArray();
            fit.MeasuredValues[f] = measured.Count(b => b);
            fit.CountedValues[f] = m - fit.MeasuredValues[f];
            fit.GroupValues[f] = G;
            bool sigmaGiven = !double.IsNaN(givenSigma2);
            fit.ResidualVarianceGivenValues[f] = sigmaGiven;
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

            var data = new Data
            {
                X = x, Group = rowGroup, G = G, IsMeasured = measured, Link = link,
                Rows = Enumerable.Range(0, m).Select(i => x.Row(i).ToArray()).ToArray(),
                LogChoose = rows.Select(s => counts[f, s] is BinomialCount c ? LogChoose(c.Trials, c.Successes) : 0).ToArray(),
                Y = rows.Select(s => responses[f, s]).ToArray(),
                K = rows.Select(s => counts[f, s]?.Successes ?? 0).ToArray(),
                N = rows.Select(s => counts[f, s]?.Trials ?? 0).ToArray(),
                Depth = rows.Select(s => depth?[s] ?? 0).ToArray(),
                MeasuredRows = Enumerable.Range(0, G).Select(_ => new List<int>()).ToArray(),
                CountedRows = Enumerable.Range(0, G).Select(_ => new List<int>()).ToArray(),
            };
            for (int i = 0; i < m; i++) (measured[i] ? data.MeasuredRows : data.CountedRows)[rowGroup[i]].Add(i);

            bool reml = estimator == VarianceEstimator.Reml;
            double logSigmaGiven = sigmaGiven ? 0.5 * Math.Log(givenSigma2) : double.NaN;
            var (logTau0, logSigma0) = StartingVariances(data);
            var beta = StartingBeta(data);

            // The outer objective: −(profiled) log-likelihood at the variances, β maximized inside (warm-started).
            double Outer(double logTau, double logSigma)
            {
                if (!(logTau >= LogSdMin && logTau <= LogSdMax && logSigma >= LogSdMin && logSigma <= LogSdMax)) return double.MaxValue;
                var inner = MaximizeBeta(data, beta, logTau, logSigma);
                if (inner is not (Vector<double> b, double ll, Matrix<double> h)) return double.MaxValue;
                beta = b;
                double value = ll;
                if (reml)
                {
                    double logDet;
                    try { logDet = h.Cholesky().DeterminantLn; }
                    catch (ArgumentException) { return double.MaxValue; }
                    value += -0.5 * logDet + 0.5 * p * Math.Log(2 * Math.PI);
                }
                return double.IsFinite(value) ? -value : double.MaxValue;
            }

            double logTauHat, logSigmaHat;
            if (sigmaGiven)
            {
                logSigmaHat = logSigmaGiven;
                logTauHat = MinimizeOneDimension(u => Outer(u, logSigmaGiven), logTau0);
            }
            else
            {
                try
                {
                    var r = NelderMeadSimplex.Minimum(ObjectiveFunction.Value(v => Outer(v[0], v[1])),
                        Vector<double>.Build.DenseOfArray(new[] { logTau0, logSigma0 }), Vector<double>.Build.DenseOfArray(new[] { 0.5, 0.5 }),
                        1e-10, 4000);
                    (logTauHat, logSigmaHat) = Polish(Outer, r.MinimizingPoint[0], r.MinimizingPoint[1]);
                }
                catch (MaximumIterationsException)
                {
                    fit.StatusValues[f] = FeatureFitStatus.NotConverged;
                    return;
                }
            }

            double objective = Outer(logTauHat, logSigmaHat);
            // The boundary τ = 0, as MixedModel keeps θ = 0 when it is at least as good (here τ = e^LogSdMin, negligible).
            double atZero = Outer(LogSdMin, logSigmaHat);
            if (atZero <= objective) (logTauHat, objective) = (LogSdMin, atZero);
            var final = MaximizeBeta(data, beta, logTauHat, logSigmaHat);
            if (!(objective < double.MaxValue) || final is not (Vector<double> betaHat, double _, Matrix<double> hessian))
            {
                fit.StatusValues[f] = FeatureFitStatus.NotConverged;
                return;
            }
            if (betaHat.Any(v => Math.Abs(v) > SeparationLimit))
            {
                fit.StatusValues[f] = FeatureFitStatus.Separated;
                return;
            }
            Matrix<double> covariance;
            try { covariance = hessian.Cholesky().Solve(Matrix<double>.Build.DenseIdentity(p)); }
            catch (ArgumentException)
            {
                fit.StatusValues[f] = FeatureFitStatus.Separated;
                return;
            }

            fit.ResidualVarianceValues[f] = Math.Exp(2 * logSigmaHat);
            fit.GroupVarianceValues[f] = logTauHat <= LogSdMin ? 0 : Math.Exp(2 * logTauHat);
            fit.LogLikelihoodValues[f] = -objective;
            var (dfWithin, dfBetween, isBetween) = MixedModel.ContainmentDf(x, rowGroup, G, m);
            for (int j = 0; j < p; j++)
            {
                fit.CoefficientMatrix[f, j] = betaHat[j];
                fit.StandardErrorMatrix[f, j] = Math.Sqrt(covariance[j, j]);
                double df = isBetween[j] ? dfBetween : dfWithin;
                fit.DfMatrix[f, j] = df > 0 ? df : double.NaN;
            }
            fit.StatusValues[f] = FeatureFitStatus.Fitted;
        }

        /// <summary>Starting log SDs: from the measured samples (pooled within-group and between-group spread), else 0.</summary>
        private static (double logTau, double logSigma) StartingVariances(Data d)
        {
            double within = 0;
            int dfWithin = 0;
            var means = new List<double>();
            for (int g = 0; g < d.G; g++)
            {
                var ys = d.MeasuredRows[g].Select(i => d.Y[i]).ToList();
                if (ys.Count == 0) continue;
                double mean = ys.Average();
                means.Add(mean);
                within += ys.Sum(y => (y - mean) * (y - mean));
                dfWithin += ys.Count - 1;
            }
            double sigma = dfWithin > 0 && within > 0 ? Math.Sqrt(within / dfWithin) : 1;
            double between = means.Count > 1 ? means.Select(v => (v - means.Average()) * (v - means.Average())).Sum() / (means.Count - 1) : 0;
            double tau = Math.Sqrt(Math.Max(between, 0.01 * sigma * sigma));
            return (Math.Clamp(Math.Log(tau), LogSdMin + 1, LogSdMax - 1), Math.Clamp(Math.Log(sigma), LogSdMin + 1, LogSdMax - 1));
        }

        /// <summary>Starting β: least squares on the measured samples when they identify it, else 0.</summary>
        private static Vector<double> StartingBeta(Data d)
        {
            int p = d.X.ColumnCount;
            var measuredRows = Enumerable.Range(0, d.Y.Length).Where(i => d.IsMeasured[i]).ToArray();
            if (measuredRows.Length > p)
            {
                var xm = Matrix<double>.Build.Dense(measuredRows.Length, p, (i, j) => d.X[measuredRows[i], j]);
                if (LinearModel.IsFullRank(xm))
                    return xm.QR().Solve(Vector<double>.Build.Dense(measuredRows.Length, i => d.Y[measuredRows[i]]));
            }
            return Vector<double>.Build.Dense(p);
        }

        /// <summary>
        /// Newton steps on (log τ, log σ) from Nelder-Mead's point, with central-difference derivatives (step 1e-4); a
        /// step is kept only if it lowers the objective. Nelder-Mead stops with the variances good to about 1e-5; this
        /// takes them to the objective's own precision.
        /// </summary>
        private static (double, double) Polish(Func<double, double, double> objective, double u, double v)
        {
            const double h = 1e-4;
            double f0 = objective(u, v);
            for (int iter = 0; iter < 5; iter++)
            {
                if (u - 2 * h < LogSdMin || u + 2 * h > LogSdMax || v - 2 * h < LogSdMin || v + 2 * h > LogSdMax) break;
                double fu1 = objective(u + h, v), fu0 = objective(u - h, v), fv1 = objective(u, v + h), fv0 = objective(u, v - h);
                double fpp = objective(u + h, v + h), fpm = objective(u + h, v - h), fmp = objective(u - h, v + h), fmm = objective(u - h, v - h);
                double gu = (fu1 - fu0) / (2 * h), gv = (fv1 - fv0) / (2 * h);
                double huu = (fu1 - 2 * f0 + fu0) / (h * h), hvv = (fv1 - 2 * f0 + fv0) / (h * h), huv = (fpp - fpm - fmp + fmm) / (4 * h * h);
                double det = huu * hvv - huv * huv;
                if (!(huu > 0 && det > 0) || !double.IsFinite(det)) break;
                double du = -(hvv * gu - huv * gv) / det, dv = -(huu * gv - huv * gu) / det;
                double f1 = objective(u + du, v + dv);
                if (!(f1 < f0)) break;
                (u, v, f0) = (u + du, v + dv, f1);
                if (Math.Max(Math.Abs(du), Math.Abs(dv)) < 1e-10) break;
            }
            return (u, v);
        }

        /// <summary>Brent's method over [LogSdMin, LogSdMax], bracketed by a coarse grid around the start.</summary>
        private static double MinimizeOneDimension(Func<double, double> f, double start)
        {
            const double step = 0.5;
            double bestU = start, best = f(start);
            for (double u = LogSdMin; u <= LogSdMax + 1e-12; u += step)
            {
                double o = f(u);
                if (o < best) { best = o; bestU = u; }
            }
            return MixedModel.Brent(f, Math.Max(LogSdMin, bestU - step), Math.Min(LogSdMax, bestU + step), 1e-10);
        }

        /// <summary>
        /// β maximizing ℓ at the given variances, by Newton with step halving, from <paramref name="start"/>; with ℓ and
        /// H = −∂²ℓ/∂β² there. The quadrature rules are held fixed during the Newton steps and re-centred at the optimum,
        /// until re-centring no longer moves them. Null when ℓ cannot be evaluated.
        /// </summary>
        private static (Vector<double> Beta, double LogLik, Matrix<double> H)? MaximizeBeta(Data d, Vector<double> start, double logTau, double logSigma)
        {
            int p = start.Count;
            var beta = start.ToArray();
            var grad = new double[p];
            var h = new double[p, p];
            var rules = Centre(d, beta, logTau, logSigma);
            double ll = Evaluate(d, beta, logTau, logSigma, rules, grad, h);
            if (!double.IsFinite(ll))
            {
                beta = new double[p];
                rules = Centre(d, beta, logTau, logSigma);
                ll = Evaluate(d, beta, logTau, logSigma, rules, grad, h);
                if (!double.IsFinite(ll)) return null;
            }
            for (int round = 0; round < 10; round++)
            {
                for (int iter = 0; iter < 100; iter++)
                {
                    var hm = Matrix<double>.Build.DenseOfArray(h);
                    var g = Vector<double>.Build.DenseOfArray(grad);
                    Vector<double> step;
                    try { step = hm.Cholesky().Solve(g); }
                    catch (ArgumentException) { step = g / Math.Max(1, hm.Diagonal().AbsoluteMaximum()); }   // not concave here: a gradient step
                    bool moved = false, small = false;
                    double t = 1;
                    for (int half = 0; half < 40; half++, t /= 2)
                    {
                        var candidate = beta.Select((b, j) => b + t * step[j]).ToArray();
                        if (candidate.Any(v => Math.Abs(v) > 4 * SeparationLimit)) continue;
                        double c = Evaluate(d, candidate, logTau, logSigma, rules, null, null);
                        if (!(double.IsFinite(c) && c >= ll - 1e-12)) continue;
                        small = t * step.AbsoluteMaximum() < 1e-10 * Math.Max(1, beta.Max(Math.Abs));
                        moved = c > ll;
                        beta = candidate;
                        break;
                    }
                    ll = Evaluate(d, beta, logTau, logSigma, rules, grad, h);
                    if (!moved || small) break;
                }
                if (d.CountedRows.All(r => r.Count == 0)) break;
                var next = Centre(d, beta, logTau, logSigma);
                bool settled = next.IsCloseTo(rules);
                rules = next;
                ll = Evaluate(d, beta, logTau, logSigma, rules, grad, h);
                if (settled || !double.IsFinite(ll)) break;
            }
            if (!double.IsFinite(ll)) return null;
            return (Vector<double>.Build.DenseOfArray(beta), ll, Matrix<double>.Build.DenseOfArray(h));
        }

        /// <summary>
        /// One feature's quadrature rules at given β and variances, per group with counted samples: the u rule's centre and
        /// scale, and per u node and counted sample (index node × counts + count), the e rule's.
        /// </summary>
        private sealed class Rules
        {
            public required double[] UCentre, UScale;
            public required double[][] ECentre, EScale;

            /// <summary>Whether every centre moved by less than 1e-8 of its scale, and every scale by less than 1e-8 of itself.</summary>
            public bool IsCloseTo(Rules other)
            {
                const double tolerance = 1e-8;
                static bool Close(double c0, double s0, double c1, double s1) =>
                    Math.Abs(c0 - c1) <= tolerance * s0 && Math.Abs(s0 - s1) <= tolerance * s0;
                for (int g = 0; g < UCentre.Length; g++)
                {
                    if (ECentre[g] == null) continue;
                    if (!Close(UCentre[g], UScale[g], other.UCentre[g], other.UScale[g])) return false;
                    for (int i = 0; i < ECentre[g].Length; i++)
                        if (!Close(ECentre[g][i], EScale[g][i], other.ECentre[g][i], other.EScale[g][i])) return false;
                }
                return true;
            }
        }

        /// <summary>
        /// The adaptive rules at β (see the class remarks). The u rule: the joint mode of the group's integrand over
        /// (u, e_1..e_m), found by Newton on its arrow-shaped Hessian, and 1/√ of u's marginal precision there,
        /// 1/v + Σ w_c (1/σ²) / (1/σ² + w_c), w_c the count's curvature in its logit times b². Each e rule: the mode and
        /// curvature of the count's integrand at that u node, from the mode's linear shift with u.
        /// </summary>
        private static Rules Centre(Data d, double[] beta, double logTau, double logSigma)
        {
            double s2 = Math.Exp(2 * logSigma), t2 = Math.Exp(2 * logTau), sigma = Math.Exp(logSigma);
            var link = d.Link;
            double b = link.Slope, rho = link.Overdispersion;
            var rules = new Rules { UCentre = new double[d.G], UScale = new double[d.G], ECentre = new double[d.G][], EScale = new double[d.G][] };
            for (int g = 0; g < d.G; g++)
            {
                var counted = d.CountedRows[g];
                int nc = counted.Count;
                if (nc == 0) continue;
                int n = d.MeasuredRows[g].Count;
                double sr = 0;
                foreach (int i in d.MeasuredRows[g]) sr += d.Y[i] - Dot(d.Rows[i], beta);
                double mean = n > 0 ? t2 * sr / (s2 + n * t2) : 0, variance = n > 0 ? t2 * s2 / (s2 + n * t2) : t2;
                var offset = new double[nc];
                for (int ci = 0; ci < nc; ci++)
                {
                    int i = counted[ci];
                    offset[ci] = link.Intercept + b * Dot(d.Rows[i], beta) + link.DepthSlope * d.Depth[i];
                }

                double Joint(double u, double[] e)
                {
                    double v = -(u - mean) * (u - mean) / (2 * variance);
                    for (int ci = 0; ci < nc; ci++)
                    {
                        int i = counted[ci];
                        v += -e[ci] * e[ci] / (2 * s2) + LogCountProbability(d.K[i], d.N[i], offset[ci] + b * (u + e[ci]), rho, d.LogChoose[i], out _, out _);
                    }
                    return v;
                }
                // Curvatures at (u, e): w_c = b² × the count's curvature in its logit, floored at 0 so the system stays positive definite.
                var w = new double[nc];
                var ge = new double[nc];
                double Gradient(double u, double[] e)
                {
                    double gu = -(u - mean) / variance;
                    for (int ci = 0; ci < nc; ci++)
                    {
                        int i = counted[ci];
                        LogCountProbability(d.K[i], d.N[i], offset[ci] + b * (u + e[ci]), rho, d.LogChoose[i], out double d1, out double d2);
                        w[ci] = Math.Max(-b * b * d2, 0);
                        gu += b * d1;
                        ge[ci] = -e[ci] / s2 + b * d1;
                    }
                    return gu;
                }

                double uHat = mean;
                var eHat = new double[nc];
                double objective = Joint(uHat, eHat);
                var de = new double[nc];
                var candidate = new double[nc];
                for (int iter = 0; iter < 100; iter++)
                {
                    double gu = Gradient(uHat, eHat);
                    double schur = 1 / variance, rhs = gu;
                    for (int ci = 0; ci < nc; ci++)
                    {
                        double dc = 1 / s2 + w[ci];
                        schur += w[ci] - w[ci] * w[ci] / dc;
                        rhs -= w[ci] * ge[ci] / dc;
                    }
                    double du = rhs / schur, largest = Math.Abs(du);
                    for (int ci = 0; ci < nc; ci++)
                    {
                        de[ci] = (ge[ci] - w[ci] * du) / (1 / s2 + w[ci]);
                        largest = Math.Max(largest, Math.Abs(de[ci]));
                    }
                    double t = 1;
                    int half = 0;
                    for (; half < 40; half++, t /= 2)
                    {
                        for (int ci = 0; ci < nc; ci++) candidate[ci] = eHat[ci] + t * de[ci];
                        double c = Joint(uHat + t * du, candidate);
                        if (c >= objective - 1e-12)
                        {
                            uHat += t * du;
                            Array.Copy(candidate, eHat, nc);
                            objective = c;
                            break;
                        }
                    }
                    if (half == 40 || t * largest < 1e-10) break;
                }
                Gradient(uHat, eHat);
                double precision = 1 / variance;
                for (int ci = 0; ci < nc; ci++) precision += w[ci] * (1 / s2) / (1 / s2 + w[ci]);
                double uScale = 1 / Math.Sqrt(precision);
                rules.UCentre[g] = uHat;
                rules.UScale[g] = uScale;

                var centres = new double[UNodes * nc];
                var scales = new double[UNodes * nc];
                for (int a = 0; a < UNodes; a++)
                {
                    double u = uHat + Math.Sqrt(2) * uScale * UQuadrature.X[a];
                    for (int ci = 0; ci < nc; ci++)
                    {
                        int i = counted[ci];
                        double start = eHat[ci] - (u - uHat) * w[ci] / (1 / s2 + w[ci]);
                        (centres[a * nc + ci], scales[a * nc + ci]) =
                            CellRule(d.K[i], d.N[i], d.LogChoose[i], offset[ci] + b * u, b, rho, sigma, double.IsFinite(start) ? start : 0);
                    }
                }
                rules.ECentre[g] = centres;
                rules.EScale[g] = scales;
            }
            return rules;
        }

        /// <summary>
        /// The adaptive rule for log ∫ N(z; 0, sd²) P(k of n | logit = offset + slope·z) dz: the integrand's mode, by
        /// Newton with step halving from <paramref name="start"/>, and 1/√ of its curvature there.
        /// </summary>
        private static (double Centre, double Scale) CellRule(int k, int n, double logChoose, double offset, double slope, double rho, double sd, double start)
        {
            double precision0 = 1 / (sd * sd);
            double Value(double z) => -0.5 * z * z * precision0 + LogCountProbability(k, n, offset + slope * z, rho, logChoose, out _, out _);
            double zHat = start, value = Value(zHat);
            for (int iter = 0; iter < 100; iter++)
            {
                LogCountProbability(k, n, offset + slope * zHat, rho, logChoose, out double d1, out double d2);
                double curvature = precision0 + Math.Max(-slope * slope * d2, 0);
                double step = (-zHat * precision0 + slope * d1) / curvature;
                double t = 1;
                int half = 0;
                for (; half < 40; half++, t /= 2)
                {
                    double c = Value(zHat + t * step);
                    if (c >= value - 1e-14)
                    {
                        zHat += t * step;
                        value = c;
                        break;
                    }
                }
                if (half == 40 || Math.Abs(t * step) * Math.Sqrt(curvature) < 1e-10) break;
            }
            LogCountProbability(k, n, offset + slope * zHat, rho, logChoose, out _, out double d2Hat);
            return (zHat, 1 / Math.Sqrt(precision0 + Math.Max(-slope * slope * d2Hat, 0)));
        }

        /// <summary>log ∫ N(z; 0, sd²) P(k of n | logit = offset + slope·z) dz, by the e rule at (centre, scale).</summary>
        private static double LogCellIntegral(int k, int n, double logChoose, double offset, double slope, double rho, double sd, double centre, double scale)
        {
            double lse = double.NegativeInfinity, logNorm = -0.5 * (LogTwoPi + 2 * Math.Log(sd)) + Math.Log(Math.Sqrt(2) * scale);
            for (int j = 0; j < ENodes; j++)
            {
                double z = centre + Math.Sqrt(2) * scale * EQuadrature.X[j];
                lse = LogAdd(lse, ELogWeights[j] + logNorm - z * z / (2 * sd * sd)
                    + LogCountProbability(k, n, offset + slope * z, rho, logChoose, out _, out _));
            }
            return lse;
        }

        private static readonly double LogTwoPi = Math.Log(2 * Math.PI);

        /// <summary>
        /// ℓ(β) at the given variances, by the given quadrature rules; with <paramref name="grad"/> and <paramref name="h"/>,
        /// also its gradient and H = −∂²ℓ/∂β², exact for those rules.
        /// </summary>
        /// <remarks>
        /// Per group: the measured samples give N(Xβ, σ²I + τ²J), and leave u normal with mean m = κ Σ r and variance
        /// v = τ²σ² / (σ² + nτ²), κ = τ² / (σ² + nτ²), r the measured residuals. On node a of the u rule, F_a is the node's
        /// weight times N(u_a; m, v) times Π_c q_c(u_a), q_c(u) a counted sample's probability integrated over e by its
        /// rule; L = log Σ_a F_a. The nodes are fixed, so β enters through m (∂m/∂β = −κ Σ x, the measured rows' sum) and
        /// each count's latent location xᵀβ + u_a:
        /// ∂ log F_a = −κ Σx (u_a − m)/v + Σ_c b r₁ x_c, ∂² log F_a = −κ² Σx Σxᵀ / v + Σ_c b² (r₂ − r₁²) x_c x_cᵀ, where r₁ and
        /// r₂ are the e-weighted means of ∂ log P/∂ℓ and (∂ log P/∂ℓ)² + ∂² log P/∂ℓ² (ℓ the per-trial logit, b the link's
        /// slope); then ∂L = Σ π_a ∂ log F_a and ∂²L = Σ π_a (∂² log F_a + ∂ log F_a ∂ log F_aᵀ) − ∂L ∂Lᵀ, π_a = F_a / Σ F.
        /// </remarks>
        private static double Evaluate(Data d, double[] beta, double logTau, double logSigma, Rules rules, double[]? grad, double[,]? h)
        {
            int p = beta.Length;
            bool deriv = grad != null;
            if (deriv) { Array.Clear(grad!); Array.Clear(h!); }
            double s2 = Math.Exp(2 * logSigma), t2 = Math.Exp(2 * logTau), total = 0;
            double logNormE = -0.5 * (LogTwoPi + 2 * logSigma);
            var link = d.Link;
            double slope = link.Slope;
            var sx = new double[p];
            var logF = new double[UNodes];
            var gF = deriv ? new double[UNodes, p] : null;
            var hF = deriv ? new double[UNodes, p, p] : null;
            var gL = new double[p];
            var lw = new double[ENodes];
            var s1 = new double[ENodes];
            var s2c = new double[ENodes];
            for (int g = 0; g < d.G; g++)
            {
                var measured = d.MeasuredRows[g];
                int n = measured.Count;
                Array.Clear(sx);
                double sr = 0, sr2 = 0;
                foreach (int i in measured)
                {
                    var xi = d.Rows[i];
                    double r = d.Y[i] - Dot(xi, beta);
                    sr += r;
                    sr2 += r * r;
                    for (int j = 0; j < p; j++)
                    {
                        sx[j] += xi[j];
                        if (!deriv) continue;
                        grad![j] += xi[j] * r / s2;
                        for (int k = 0; k < p; k++) h![j, k] += xi[j] * xi[k] / s2;
                    }
                }
                double kappa = 0, mean = 0, variance = t2;
                if (n > 0)
                {
                    double denom = s2 + n * t2;
                    total += -0.5 * (n * LogTwoPi + (n - 1) * Math.Log(s2) + Math.Log(denom) + sr2 / s2 - t2 * sr * sr / (s2 * denom));
                    if (deriv)
                    {
                        double c = t2 / (s2 * denom);
                        for (int j = 0; j < p; j++)
                        {
                            grad![j] -= c * sx[j] * sr;
                            for (int k = 0; k < p; k++) h![j, k] -= c * sx[j] * sx[k];
                        }
                    }
                    kappa = t2 / denom;
                    mean = kappa * sr;
                    variance = t2 * s2 / denom;
                }

                var counted = d.CountedRows[g];
                int nc = counted.Count;
                if (nc == 0) continue;
                double uCentre = rules.UCentre[g], uScale = rules.UScale[g];
                var eCentre = rules.ECentre[g];
                var eScale = rules.EScale[g];
                double logNormU = -0.5 * (LogTwoPi + Math.Log(variance)) + Math.Log(Math.Sqrt(2) * uScale);
                for (int a = 0; a < UNodes; a++)
                {
                    double u = uCentre + Math.Sqrt(2) * uScale * UQuadrature.X[a];
                    double du = u - mean;
                    logF[a] = ULogWeights[a] + logNormU - du * du / (2 * variance);
                    if (deriv)
                        for (int j = 0; j < p; j++)
                        {
                            gF![a, j] = -kappa * sx[j] * du / variance;
                            for (int k = 0; k < p; k++) hF![a, j, k] = -kappa * kappa * sx[j] * sx[k] / variance;
                        }
                    for (int ci = 0; ci < nc; ci++)
                    {
                        int i = counted[ci];
                        var xi = d.Rows[i];
                        double location = link.Intercept + slope * (Dot(xi, beta) + u) + link.DepthSlope * d.Depth[i];
                        double centre = eCentre[a * nc + ci], scale = eScale[a * nc + ci];
                        double logNorm = logNormE + Math.Log(Math.Sqrt(2) * scale);
                        double top = double.NegativeInfinity;
                        for (int b = 0; b < ENodes; b++)
                        {
                            double e = centre + Math.Sqrt(2) * scale * EQuadrature.X[b];
                            lw[b] = ELogWeights[b] + logNorm - e * e / (2 * s2)
                                + LogCountProbability(d.K[i], d.N[i], location + slope * e, link.Overdispersion, d.LogChoose[i], out s1[b], out s2c[b]);
                            top = Math.Max(top, lw[b]);
                        }
                        // The node terms relative to the largest (one exponential each), then log q = top + log Σ.
                        double sum = 0;
                        for (int b = 0; b < ENodes; b++) sum += lw[b] = Math.Exp(lw[b] - top);
                        double logq = top + Math.Log(sum);
                        logF[a] += logq;
                        if (!deriv) continue;
                        double r1 = 0, r2 = 0;
                        for (int b = 0; b < ENodes; b++)
                        {
                            double w = lw[b] / sum;
                            r1 += w * s1[b];
                            r2 += w * (s1[b] * s1[b] + s2c[b]);
                        }
                        double curvature = slope * slope * (r2 - r1 * r1);
                        for (int j = 0; j < p; j++)
                        {
                            gF![a, j] += slope * r1 * xi[j];
                            for (int k = 0; k < p; k++) hF![a, j, k] += curvature * xi[j] * xi[k];
                        }
                    }
                }
                double lse = double.NegativeInfinity;
                for (int a = 0; a < UNodes; a++) lse = LogAdd(lse, logF[a]);
                total += lse;
                if (!deriv) continue;
                Array.Clear(gL);
                for (int a = 0; a < UNodes; a++)
                {
                    double pi = Math.Exp(logF[a] - lse);
                    for (int j = 0; j < p; j++) gL[j] += pi * gF![a, j];
                }
                for (int a = 0; a < UNodes; a++)
                {
                    double pi = Math.Exp(logF[a] - lse);
                    if (pi == 0) continue;
                    for (int j = 0; j < p; j++)
                        for (int k = 0; k < p; k++)
                            h![j, k] -= pi * (hF![a, j, k] + gF![a, j] * gF[a, k]);
                }
                for (int j = 0; j < p; j++)
                {
                    grad![j] += gL[j];
                    for (int k = 0; k < p; k++) h![j, k] += gL[j] * gL[k];
                }
            }
            return total;
        }

        private static double Dot(double[] x, double[] beta)
        {
            double s = 0;
            for (int j = 0; j < x.Length; j++) s += x[j] * beta[j];
            return s;
        }

        /// <summary>
        /// log P(k of n) for a per-trial chance g = 1 / (1 + e^(−logit)): binomial when ρ = 0, else beta-binomial with
        /// mean g and intra-cell correlation ρ, written as products so no gamma function is needed.
        /// </summary>
        internal static double LogCountProbability(int k, int n, double logit, double rho) =>
            LogCountProbability(k, n, logit, rho, LogChoose(n, k), out _, out _);

        /// <summary>
        /// <see cref="LogCountProbability(int, int, double, double)"/> with log C(n, k) given, and its first and second
        /// derivatives in the logit: binomial k − ng and −ng(1 − g); beta-binomial through g' = g(1 − g), g'' = g'(1 − 2g).
        /// </summary>
        internal static double LogCountProbability(int k, int n, double logit, double rho, double logChoose, out double d1, out double d2)
        {
            double logG = -Log1pExp(-logit), log1mG = -Log1pExp(logit);
            double g = Math.Exp(logG), h = Math.Exp(log1mG);
            if (rho == 0)
            {
                d1 = k - n * g;
                d2 = -n * g * h;
                return logChoose + k * logG + (n - k) * log1mG;
            }
            double s = (1 - rho) / rho, sum = logChoose, t1 = 0, t2 = 0;
            for (int j = 0; j < k; j++)
            {
                double a = g * s + j;
                sum += Math.Log(a);
                t1 += s / a;
                t2 -= s * s / (a * a);
            }
            for (int j = 0; j < n - k; j++)
            {
                double b = h * s + j;
                sum += Math.Log(b);
                t1 -= s / b;
                t2 -= s * s / (b * b);
            }
            for (int j = 0; j < n; j++) sum -= Math.Log(s + j);
            double g1 = g * h, g2 = g1 * (h - g);
            d1 = t1 * g1;
            d2 = t2 * g1 * g1 + t1 * g2;
            return sum;
        }

        private static double LogChoose(int n, int k)
        {
            double sum = 0;
            for (int j = 1; j <= Math.Min(k, n - k); j++) sum += Math.Log((double)(n - Math.Min(k, n - k) + j) / j);
            return sum;
        }

        /// <summary>log(1 + e^x) without overflow.</summary>
        private static double Log1pExp(double x) => x > 0 ? x + Math.Log(1 + Math.Exp(-x)) : Math.Log(1 + Math.Exp(x));

        private static double LogAdd(double a, double b)
        {
            if (double.IsNegativeInfinity(a)) return b;
            if (double.IsNegativeInfinity(b)) return a;
            return Math.Max(a, b) + Math.Log(1 + Math.Exp(-Math.Abs(a - b)));
        }
    }
}
