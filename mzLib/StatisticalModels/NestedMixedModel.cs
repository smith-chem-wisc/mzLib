using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics.Optimization;

namespace StatisticalModels
{
    /// <summary>The random intercepts a feature's model carries.</summary>
    public enum RandomStructure
    {
        /// <summary>One random intercept per sample: <c>(1 | sample)</c>. No individual was given.</summary>
        Sample,
        /// <summary>
        /// One random intercept per individual: <c>(1 | individual)</c>. Individuals were given but each has a single
        /// sample among the feature's values, so a sample-within-individual intercept would be the same grouping.
        /// </summary>
        Individual,
        /// <summary>
        /// Individual, and sample within individual: <c>(1 | individual) + (1 | individual:sample)</c> (QuantProject
        /// GR-6), so comparisons within and between individuals are both estimable.
        /// </summary>
        IndividualAndSample,
    }

    /// <summary>
    /// One feature's data for <see cref="NestedMixedModel"/>: a response per observation (e.g. a log2 intensity per
    /// peptide and sample), the design columns every feature shares, the feature's own nuisance columns (e.g. its
    /// peptides, treatment-coded), and the sample and individual each observation came from. The arrays are held, not
    /// copied.
    /// </summary>
    public sealed class NestedMixedModelProblem
    {
        /// <param name="response">One value per observation. A non-finite value is missing and omitted.</param>
        /// <param name="common">[observation, column] in the order of the fit's common names; finite.</param>
        /// <param name="nuisance">[observation, column] fitted but not reported (e.g. peptide effects); finite; or null.</param>
        /// <param name="sample">The sample of each observation; non-empty.</param>
        /// <param name="individual">The individual of each observation, all non-empty; or null when none is known.</param>
        public NestedMixedModelProblem(IReadOnlyList<double> response, double[,] common, double[,]? nuisance,
            IReadOnlyList<string> sample, IReadOnlyList<string?>? individual)
        {
            ArgumentNullException.ThrowIfNull(response);
            ArgumentNullException.ThrowIfNull(common);
            ArgumentNullException.ThrowIfNull(sample);
            int n = response.Count;
            if (common.GetLength(0) != n)
                throw new ArgumentException($"The common design has {common.GetLength(0)} rows but there are {n} responses.", nameof(common));
            if (nuisance is not null && nuisance.GetLength(0) != n)
                throw new ArgumentException($"The nuisance design has {nuisance.GetLength(0)} rows but there are {n} responses.", nameof(nuisance));
            if (sample.Count != n)
                throw new ArgumentException($"There are {sample.Count} sample labels for {n} responses.", nameof(sample));
            if (individual is not null && individual.Count != n)
                throw new ArgumentException($"There are {individual.Count} individual labels for {n} responses.", nameof(individual));
            CheckFinite(common, nameof(common));
            if (nuisance is not null) CheckFinite(nuisance, nameof(nuisance));
            for (int i = 0; i < n; i++)
            {
                if (string.IsNullOrEmpty(sample[i]))
                    throw new ArgumentException($"Observation {i} has no sample.", nameof(sample));
                if (individual is not null && string.IsNullOrEmpty(individual[i]))
                    throw new ArgumentException($"Observation {i} has no individual; give one for every observation, or none.", nameof(individual));
            }
            Response = response;
            Common = common;
            Nuisance = nuisance;
            Sample = sample;
            Individual = individual;
        }

        internal IReadOnlyList<double> Response { get; }
        internal double[,] Common { get; }
        internal double[,]? Nuisance { get; }
        internal IReadOnlyList<string> Sample { get; }
        internal IReadOnlyList<string?>? Individual { get; }

        private static void CheckFinite(double[,] m, string name)
        {
            for (int i = 0; i < m.GetLength(0); i++)
                for (int j = 0; j < m.GetLength(1); j++)
                    if (!double.IsFinite(m[i, j]))
                        throw new ArgumentException($"The design is not finite at observation {i}, column {j}.", name);
        }
    }

    /// <summary>
    /// Per-feature linear mixed model fits (<see cref="NestedMixedModel.Fit"/>). A value that could not be computed is
    /// NaN and the reason is in <see cref="Status"/>. Read-only to callers.
    /// </summary>
    public sealed class NestedMixedModelFit
    {
        internal NestedMixedModelFit(int features, IReadOnlyList<string> commonNames)
        {
            CommonNames = commonNames;
            int p = commonNames.Count;
            double[] Nan(int n) => Enumerable.Repeat(double.NaN, n).ToArray();
            CoefficientValues = Nan(features * p);
            UnscaledValues = Nan(features * p * p);
            ResidualVarianceValues = Nan(features);
            IndividualVarianceValues = Nan(features);
            SampleVarianceValues = Nan(features);
            RemlValues = Nan(features);
            ObservedValues = new int[features];
            StatusValues = new FeatureFitStatus[features];
            StructureValues = new RandomStructure?[features];
            Jacobians = new double[features][];
            VarParCovariances = new double[features][];
        }

        internal double[] CoefficientValues { get; }
        internal double[] UnscaledValues { get; }
        internal double[] ResidualVarianceValues { get; }
        internal double[] IndividualVarianceValues { get; }
        internal double[] SampleVarianceValues { get; }
        internal double[] RemlValues { get; }
        internal int[] ObservedValues { get; }
        internal FeatureFitStatus[] StatusValues { get; }
        internal RandomStructure?[] StructureValues { get; }
        /// <summary>Per feature: ∂Cov(β_common)/∂φ_k for each variance parameter φ = (θ…, σ), k-major, p×p each.</summary>
        internal double[]?[] Jacobians { get; }
        /// <summary>Per feature: the asymptotic covariance of φ, 2·H⁺ (H = Hessian of the REML deviance), k×k.</summary>
        internal double[]?[] VarParCovariances { get; }

        /// <summary>Names of the common design columns, in order.</summary>
        public IReadOnlyList<string> CommonNames { get; }

        /// <summary>Number of features (problems).</summary>
        public int FeatureCount => ObservedValues.Length;

        /// <summary>REML estimate of common coefficient <paramref name="coefficient"/> of <paramref name="feature"/>.</summary>
        public double Coefficient(int feature, int coefficient) => CoefficientValues[feature * CommonNames.Count + coefficient];

        /// <summary>
        /// Element (i, j) of (XᵀV⁻¹X)⁻¹ for the common columns, with V the observations' covariance divided by the residual
        /// variance: Cov(β) / σ². lme4's <c>unsc()</c>, and what moderation rescales (GR-21).
        /// </summary>
        public double UnscaledCovariance(int feature, int i, int j)
        {
            int p = CommonNames.Count;
            return UnscaledValues[(feature * p + i) * p + j];
        }

        /// <summary>REML residual variance σ² per feature.</summary>
        public IReadOnlyList<double> ResidualVariance => ResidualVarianceValues;
        /// <summary>Variance of the individual intercept; NaN when the structure has none.</summary>
        public IReadOnlyList<double> IndividualVariance => IndividualVarianceValues;
        /// <summary>Variance of the sample (or sample-within-individual) intercept; NaN when the structure has none.</summary>
        public IReadOnlyList<double> SampleVariance => SampleVarianceValues;
        /// <summary>The REML criterion (−2 × restricted log-likelihood), as lme4's <c>REMLcrit</c>.</summary>
        public IReadOnlyList<double> RemlCriterion => RemlValues;
        /// <summary>Observations with a finite response, per feature.</summary>
        public IReadOnlyList<int> Observed => ObservedValues;
        /// <summary>Whether each feature was fitted, and if not, why.</summary>
        public IReadOnlyList<FeatureFitStatus> Status => StatusValues;
        /// <summary>The random intercepts each fitted feature carried; null when it was not fitted.</summary>
        public IReadOnlyList<RandomStructure?> Structure => StructureValues;

        /// <summary>
        /// Satterthwaite degrees of freedom of the contrast c′β (lmerTest's <c>contest1D</c>):
        /// 2·(c′Cc)² / (g′Ag), with C = Cov(β), g_k = c′(∂C/∂φ_k)c, A the asymptotic covariance of the variance
        /// parameters. Not moderated. NaN for a feature that was not fitted.
        /// </summary>
        public double SatterthwaiteDf(int feature, IReadOnlyList<double> contrast)
        {
            Check(contrast, nameof(contrast));
            if (StatusValues[feature] != FeatureFitStatus.Fitted) return double.NaN;
            return Satterthwaite(feature, contrast.ToArray(), out _);
        }

        /// <summary>
        /// Denominator degrees of freedom of the joint F-test of several contrasts (lmerTest's <c>contestMD</c>): the
        /// contrasts are rotated to the eigenvectors of L C L′ (eigenvalues above √ε of the largest), each direction gets
        /// its Satterthwaite ν, and the F df is 2E/(E − q) with E = Σ ν/(ν − 2), or 2 if any ν ≤ 2. NaN when not fitted.
        /// </summary>
        public double JointSatterthwaiteDf(int feature, IReadOnlyList<IReadOnlyList<double>> rows)
        {
            ArgumentNullException.ThrowIfNull(rows);
            if (rows.Count == 0) throw new ArgumentException("No contrast rows were given.", nameof(rows));
            foreach (var r in rows) Check(r, nameof(rows));
            if (StatusValues[feature] != FeatureFitStatus.Fitted) return double.NaN;
            int p = CommonNames.Count, q = rows.Count;
            var cov = CovarianceMatrix(feature);
            var l = Matrix<double>.Build.Dense(q, p, (i, j) => rows[i][j]);
            var vll = l * cov * l.Transpose();
            var evd = vll.Evd(Symmetricity.Symmetric);
            var values = evd.EigenValues.Select(z => z.Real).ToArray();
            double largest = values.Max();
            double tol = Math.Sqrt(2.220446049250313e-16) * largest;
            var nu = new List<double>();
            for (int k = 0; k < values.Length; k++)
            {
                if (!(values[k] > Math.Max(tol, 0))) continue;
                var direction = evd.EigenVectors.Column(k) * l; // a row of the rotated contrast matrix
                nu.Add(Satterthwaite(feature, direction.ToArray(), out _));
            }
            if (nu.Count == 0) return double.NaN;
            if (nu.Count == 1) return nu[0];
            if (nu.Max() - nu.Min() < 1e-8) return nu.Average();
            if (nu.Any(v => v <= 2)) return 2;
            double e = nu.Sum(v => v / (v - 2));
            return 2 * e / (e - nu.Count);
        }

        internal Matrix<double> CovarianceMatrix(int feature)
        {
            int p = CommonNames.Count;
            double s2 = ResidualVarianceValues[feature];
            return Matrix<double>.Build.Dense(p, p, (i, j) => s2 * UnscaledCovariance(feature, i, j));
        }

        internal double Satterthwaite(int feature, double[] c, out double variance)
        {
            int p = CommonNames.Count;
            var jac = Jacobians[feature]!;
            var a = VarParCovariances[feature]!;
            int k = a.Length == 0 ? 0 : (int)Math.Round(Math.Sqrt(a.Length));
            var cov = CovarianceMatrix(feature);
            var cv = Vector<double>.Build.DenseOfArray(c);
            variance = cv * cov * cv;
            var g = new double[k];
            for (int t = 0; t < k; t++)
            {
                double s = 0;
                for (int i = 0; i < p; i++)
                    for (int j = 0; j < p; j++) s += c[i] * c[j] * jac[t * p * p + i * p + j];
                g[t] = s;
            }
            double denominator = 0;
            for (int i = 0; i < k; i++)
                for (int j = 0; j < k; j++) denominator += g[i] * a[i * k + j] * g[j];
            return denominator > 0 ? 2 * variance * variance / denominator : double.NaN;
        }

        private void Check(IReadOnlyList<double> contrast, string name)
        {
            ArgumentNullException.ThrowIfNull(contrast, name);
            if (contrast.Count != CommonNames.Count)
                throw new ArgumentException(
                    $"The contrast has {contrast.Count} weights but there are {CommonNames.Count} common coefficients " +
                    $"({string.Join(", ", CommonNames)}).", name);
            if (contrast.Any(w => !double.IsFinite(w)))
                throw new ArgumentException("A contrast weight is not finite.", name);
        }
    }

    /// <summary>
    /// A linear mixed model fitted separately for each feature, each with its own design: the response (e.g. every
    /// peptide's log2 intensity in every sample) on common fixed effects (the experimental design) plus the feature's own
    /// nuisance fixed effects (e.g. its peptides), with random intercepts for the sample, or the individual, or both
    /// nested (<see cref="RandomStructure"/>). Fitted by REML; degrees of freedom by Satterthwaite's approximation; variance
    /// moderated across features by <see cref="ModerateContrasts"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Fit.</b> lme4's formulation (Bates, Mächler, Bolker &amp; Walker 2015, J. Stat. Softw. 67:1): with relative
    /// standard deviations θ (one per random term) and Λ = diag(θ), the REML criterion is
    /// log|ΛZ′ZΛ + I| + log|X′V⁻¹X| + (m − p)(1 + log(2π r²/(m − p))), V = I + ZΛΛZ′, r² the penalised residual sum of
    /// squares, σ² = r²/(m − p). It is minimised over θ ≥ 0 by a grid, MathNet's Nelder–Mead (two random terms) or Brent
    /// (one), and Newton polishing, with every boundary (θ = 0) tried explicitly. It reproduces
    /// <c>lme4::lmer(REML = TRUE)</c> to 1e-6.
    /// </para>
    /// <para>
    /// <b>Degrees of freedom.</b> lmerTest's Satterthwaite method (Kuznetsova, Brockhoff &amp; Christensen 2017,
    /// J. Stat. Softw. 82:13): the variance parameters are φ = (θ…, σ); A = 2·H⁺, H the Hessian of the REML deviance in φ
    /// (eigenvalues ≤ 1e-8 dropped, as lmerTest); H and ∂Cov(β)/∂φ by Richardson-extrapolated central differences with
    /// numDeriv's defaults, which lmerTest uses (relative step 0.1 for the Hessian, 1e-4 for the Jacobian; 4 halvings),
    /// so the df carry the same numerical error and match lmerTest. r² is computed from the residuals, as lme4 does:
    /// the df are sensitive to the optimum near a small variance component. Never the containment rule (QuantProject GR-7).
    /// </para>
    /// <para>
    /// <b>Not fittable.</b> Too few observations for the fixed effects (<see cref="FeatureFitStatus.TooFewObservations"/>);
    /// fixed effects not of full rank on the observed values (<see cref="FeatureFitStatus.RankDeficient"/>); a random
    /// term with fewer than 2 levels or as many levels as observations, e.g. a protein with one peptide
    /// (<see cref="FeatureFitStatus.TooFewGroups"/>; lme4 refuses the same); a criterion that cannot be evaluated, or a
    /// simplex search that does not converge (<see cref="FeatureFitStatus.NotConverged"/>). Output does not depend on the
    /// thread count.
    /// </para>
    /// </remarks>
    public static class NestedMixedModel
    {
        /// <summary>Fits every problem.</summary>
        /// <param name="problems">One per feature.</param>
        /// <param name="commonNames">One name per common design column; every problem has these columns.</param>
        /// <param name="maxThreads">As <see cref="LinearModel.Fit"/>: -1 uses all cores but one. Results are identical for every value.</param>
        public static NestedMixedModelFit Fit(IReadOnlyList<NestedMixedModelProblem> problems, IReadOnlyList<string> commonNames,
            int maxThreads = -1)
        {
            ArgumentNullException.ThrowIfNull(problems);
            ArgumentNullException.ThrowIfNull(commonNames);
            if (commonNames.Count == 0) throw new ArgumentException("There are no common columns.", nameof(commonNames));
            for (int f = 0; f < problems.Count; f++)
            {
                ArgumentNullException.ThrowIfNull(problems[f], nameof(problems));
                if (problems[f].Common.GetLength(1) != commonNames.Count)
                    throw new ArgumentException(
                        $"Problem {f} has {problems[f].Common.GetLength(1)} common columns but {commonNames.Count} names were given.",
                        nameof(commonNames));
            }
            var fit = new NestedMixedModelFit(problems.Count, commonNames);
            var options = new ParallelOptions { MaxDegreeOfParallelism = LinearModel.ResolveThreads(maxThreads) };
            Parallel.ForEach(Partitioner.Create(0, Math.Max(problems.Count, 1)), options, range =>
            {
                for (int f = range.Item1; f < Math.Min(range.Item2, problems.Count); f++)
                    new FeatureModel(problems[f], commonNames.Count).Fit(f, fit);
            });
            return fit;
        }

        /// <summary>
        /// Moderated contrasts with confidence intervals, as MSstatsTMT 2.20.0 combines moderation with Satterthwaite df,
        /// extended to several condition terms (QuantProject GR-21). Per fitted feature: s² = σ²; its df = the joint
        /// Satterthwaite F df of <paramref name="varianceDfRows"/> (the condition terms); the prior is fitted from (s², df)
        /// across features as <see cref="EmpiricalBayes.FitPrior"/> with <see cref="VariancePriorEstimator.MomentsLegacy"/>
        /// (limma <c>squeezeVar(legacy = TRUE)</c>). For each contrast c: SE = s̃ · sqrt(c′Uc) with U the unscaled
        /// covariance; df = the contrast's unmoderated Satterthwaite df + d0, capped at the sum of the features' variance df;
        /// t, p and a t-based interval on that df; Benjamini–Hochberg over the features with a result. Everything after the
        /// prior is the moderation <see cref="EmpiricalBayes.ModerateContrasts"/> applies, from the same code; only the
        /// variance's df and the contrast's df differ.
        /// </summary>
        /// <param name="fit">Output of <see cref="Fit"/>.</param>
        /// <param name="contrasts">At least one; names unique; one finite weight per common coefficient, not all zero.</param>
        /// <param name="varianceDfRows">The condition terms whose joint F-test gives each feature's variance df.</param>
        /// <param name="confidenceLevel">Two-sided interval level, strictly between 0 and 1.</param>
        public static IReadOnlyList<ModeratedContrast> ModerateContrasts(NestedMixedModelFit fit,
            IReadOnlyList<ContrastWeights> contrasts, IReadOnlyList<IReadOnlyList<double>> varianceDfRows,
            double confidenceLevel = 0.95)
        {
            ArgumentNullException.ThrowIfNull(fit);
            ArgumentNullException.ThrowIfNull(contrasts);
            ArgumentNullException.ThrowIfNull(varianceDfRows);
            int p = fit.CommonNames.Count;
            EmpiricalBayes.CheckContrasts(contrasts, confidenceLevel, c =>
            {
                if (c.Weights is null || c.Weights.Count != p || c.Weights.Any(w => !double.IsFinite(w)))
                    throw new ArgumentException($"Contrast '{c.Name}' needs {p} finite weights.", nameof(contrasts));
            });
            if (varianceDfRows.Count == 0) throw new ArgumentException("No variance-df rows were given.", nameof(varianceDfRows));
            foreach (var r in varianceDfRows)
                if (r is null || r.Count != p || r.Any(w => !double.IsFinite(w)))
                    throw new ArgumentException($"Each variance-df row needs {p} finite weights.", nameof(varianceDfRows));

            int n = fit.FeatureCount;
            var s2 = new double[n];
            var df = new double[n];
            var usable = new List<int>();
            for (int f = 0; f < n; f++)
            {
                s2[f] = double.NaN;
                if (fit.Status[f] != FeatureFitStatus.Fitted) continue;
                double d = fit.JointSatterthwaiteDf(f, varianceDfRows);
                if (!(double.IsFinite(d) && d > 0)) continue;
                s2[f] = fit.ResidualVariance[f];
                df[f] = d;
                usable.Add(f);
            }
            if (usable.Count < 2)
                throw new ArgumentException($"Moderation needs at least 2 fitted features with a variance df; this fit has {usable.Count}.", nameof(fit));
            var prior = EmpiricalBayes.FitPrior(s2, df, null, null, VariancePriorEstimator.MomentsLegacy);
            return EmpiricalBayes.ModerateContrastsCore(contrasts, prior, s2, df, usable, fit.Status, confidenceLevel, (f, w) =>
            {
                double unscaled = 0, estimate = 0;
                for (int i = 0; i < p; i++)
                {
                    estimate += w[i] * fit.Coefficient(f, i);
                    for (int j = 0; j < p; j++) unscaled += w[i] * w[j] * fit.UnscaledCovariance(f, i, j);
                }
                return new FeatureContrast(estimate, Math.Sqrt(unscaled), fit.Satterthwaite(f, w, out _));
            });
        }

        /// <summary>One feature's fit: data reduced to cross-products, the REML criterion, its optimum and derivatives.</summary>
        private sealed class FeatureModel
        {
            private readonly NestedMixedModelProblem _problem;
            private readonly int _common;
            private int _m, _p, _q;
            private int[] _termOfColumn = Array.Empty<int>();
            private Matrix<double> _ztz = null!, _ztx = null!, _xtx = null!, _x = null!;
            private Vector<double> _zty = null!, _xty = null!, _y = null!;
            private List<int[]> _columns = new();

            public FeatureModel(NestedMixedModelProblem problem, int common)
            {
                _problem = problem;
                _common = common;
            }

            public void Fit(int f, NestedMixedModelFit fit)
            {
                var rows = Enumerable.Range(0, _problem.Response.Count).Where(i => double.IsFinite(_problem.Response[i])).ToArray();
                _m = rows.Length;
                fit.ObservedValues[f] = _m;
                int nuisance = _problem.Nuisance?.GetLength(1) ?? 0;
                _p = _common + nuisance;
                if (_m <= _p) { fit.StatusValues[f] = FeatureFitStatus.TooFewObservations; return; }

                var x = Matrix<double>.Build.Dense(_m, _p, (i, j) =>
                    j < _common ? _problem.Common[rows[i], j] : _problem.Nuisance![rows[i], j - _common]);
                if (!LinearModel.IsFullRank(x)) { fit.StatusValues[f] = FeatureFitStatus.RankDeficient; return; }
                var y = Vector<double>.Build.Dense(_m, i => _problem.Response[rows[i]]);

                // Random terms: individual and/or sample (within individual), levels in ordinal order.
                var sample = rows.Select(i => _problem.Sample[i]).ToArray();
                var terms = new List<string[]>();
                RandomStructure structure;
                if (_problem.Individual is null)
                {
                    terms.Add(sample);
                    structure = RandomStructure.Sample;
                }
                else
                {
                    var individual = rows.Select(i => _problem.Individual[i]!).ToArray();
                    bool oneSampleEach = individual.Zip(sample).GroupBy(t => t.First, StringComparer.Ordinal)
                        .All(g => g.Select(t => t.Second).Distinct(StringComparer.Ordinal).Count() == 1);
                    terms.Add(individual);
                    structure = RandomStructure.Individual;
                    if (!oneSampleEach)
                    {
                        terms.Add(individual.Zip(sample, (a, b) => a + "\u001f" + b).ToArray());
                        structure = RandomStructure.IndividualAndSample;
                    }
                }
                var columns = new List<int[]>(); // per term: column index (within Z) of each observation
                var termOfColumn = new List<int>();
                int q = 0;
                for (int t = 0; t < terms.Count; t++)
                {
                    var levels = terms[t].Distinct(StringComparer.Ordinal).OrderBy(s => s, StringComparer.Ordinal).ToList();
                    if (levels.Count < 2 || levels.Count >= _m) { fit.StatusValues[f] = FeatureFitStatus.TooFewGroups; return; }
                    var index = levels.Select((l, i) => (l, i)).ToDictionary(z => z.l, z => q + z.i, StringComparer.Ordinal);
                    columns.Add(terms[t].Select(l => index[l]).ToArray());
                    termOfColumn.AddRange(Enumerable.Repeat(t, levels.Count));
                    q += levels.Count;
                }
                _q = q;
                _termOfColumn = termOfColumn.ToArray();
                _ztz = Matrix<double>.Build.Dense(q, q);
                _ztx = Matrix<double>.Build.Dense(q, _p);
                _zty = Vector<double>.Build.Dense(q);
                for (int i = 0; i < _m; i++)
                    foreach (var a in columns)
                    {
                        int ca = a[i];
                        _zty[ca] += y[i];
                        for (int j = 0; j < _p; j++) _ztx[ca, j] += x[i, j];
                        foreach (var b in columns) _ztz[ca, b[i]] += 1;
                    }
                _xtx = x.TransposeThisAndMultiply(x);
                _xty = x.TransposeThisAndMultiply(y);
                _x = x;
                _y = y;
                _columns = columns;

                var theta = Minimise(terms.Count);
                if (theta is null) { fit.StatusValues[f] = FeatureFitStatus.NotConverged; return; }
                var best = Evaluate(theta);
                if (!best.Ok) { fit.StatusValues[f] = FeatureFitStatus.NotConverged; return; }
                double sigma2 = best.R2 / (_m - _p);
                int p = _common;
                for (int i = 0; i < p; i++)
                {
                    fit.CoefficientValues[f * p + i] = best.Beta![i];
                    for (int j = 0; j < p; j++) fit.UnscaledValues[(f * p + i) * p + j] = best.AInverse![i, j];
                }
                fit.ResidualVarianceValues[f] = sigma2;
                fit.RemlValues[f] = best.Criterion;
                if (structure == RandomStructure.Sample) fit.SampleVarianceValues[f] = theta[0] * theta[0] * sigma2;
                else
                {
                    fit.IndividualVarianceValues[f] = theta[0] * theta[0] * sigma2;
                    if (structure == RandomStructure.IndividualAndSample) fit.SampleVarianceValues[f] = theta[1] * theta[1] * sigma2;
                }

                // Satterthwaite pieces in φ = (θ…, σ), as lmerTest.
                int k = theta.Length + 1;
                var phi = theta.Append(Math.Sqrt(sigma2)).ToArray();
                var h = Richardson.Hessian(DevianceInPhi, phi);
                var evd = Matrix<double>.Build.DenseOfArray(h).Evd(Symmetricity.Symmetric);
                var a2 = new double[k * k];
                for (int e = 0; e < k; e++)
                {
                    double lambda = evd.EigenValues[e].Real;
                    if (!(lambda > 1e-8)) continue;
                    var v = evd.EigenVectors.Column(e);
                    for (int i = 0; i < k; i++)
                        for (int j = 0; j < k; j++) a2[i * k + j] += 2 * v[i] * v[j] / lambda;
                }
                var jac = new double[k * p * p];
                for (int t = 0; t < k; t++)
                {
                    var d = Richardson.Derivative(CovarianceInPhi, phi, t);
                    Array.Copy(d, 0, jac, t * p * p, p * p);
                }
                if (jac.Any(z => !double.IsFinite(z)) || a2.Any(z => !double.IsFinite(z)))
                {
                    fit.StatusValues[f] = FeatureFitStatus.NotConverged;
                    return;
                }
                fit.Jacobians[f] = jac;
                fit.VarParCovariances[f] = a2;
                fit.StructureValues[f] = structure;
                fit.StatusValues[f] = FeatureFitStatus.Fitted;
            }

            private readonly record struct Result(bool Ok, double Criterion, double LogDetL, double LogDetRx, double R2,
                Vector<double>? Beta, Matrix<double>? AInverse);

            private Result Evaluate(double[] theta)
            {
                try
                {
                    var lambda = Vector<double>.Build.Dense(_q, j => theta[_termOfColumn[j]]);
                    var m = Matrix<double>.Build.Dense(_q, _q, (i, j) => lambda[i] * _ztz[i, j] * lambda[j] + (i == j ? 1 : 0));
                    var chol = m.Cholesky();
                    var lztx = Matrix<double>.Build.Dense(_q, _p, (i, j) => lambda[i] * _ztx[i, j]);
                    var lzty = Vector<double>.Build.Dense(_q, i => lambda[i] * _zty[i]);
                    var s = chol.Solve(lztx);
                    var sv = chol.Solve(lzty);
                    var a = _xtx - lztx.TransposeThisAndMultiply(s);
                    var b = _xty - lztx.TransposeThisAndMultiply(sv);
                    var cholA = a.Cholesky();
                    var beta = cholA.Solve(b);
                    // The penalised residual sum of squares from the residuals themselves, as lme4 computes it:
                    // r² = |y − Xβ − ZΛu|² + |u|², u = (ΛZ′ZΛ + I)⁻¹ΛZ′(y − Xβ). Subtracting cross-products instead loses
                    // ~1e-10 of the criterion to cancellation (log2 intensities near 22), enough to move the optimum.
                    var e = _y - _x * beta;
                    var lzte = Vector<double>.Build.Dense(_q);
                    for (int i = 0; i < _m; i++) foreach (var c in _columns) lzte[c[i]] += e[i];
                    lzte = lzte.PointwiseMultiply(lambda);
                    var u = chol.Solve(lzte);
                    var lu = u.PointwiseMultiply(lambda);
                    double r2 = u.DotProduct(u);
                    for (int i = 0; i < _m; i++)
                    {
                        double ri = e[i];
                        foreach (var c in _columns) ri -= lu[c[i]];
                        r2 += ri * ri;
                    }
                    if (!(r2 > 0)) return default;
                    int dfr = _m - _p;
                    double criterion = chol.DeterminantLn + cholA.DeterminantLn + dfr * (1 + Math.Log(2 * Math.PI * r2 / dfr));
                    if (!double.IsFinite(criterion)) return default;
                    return new Result(true, criterion, chol.DeterminantLn, cholA.DeterminantLn, r2, beta,
                        cholA.Solve(Matrix<double>.Build.DenseIdentity(_p)));
                }
                catch (Exception e) when (e is ArgumentException or InvalidOperationException)
                {
                    return default; // not positive definite: never an exception out of the parallel fit
                }
            }

            private double Criterion(double[] theta)
            {
                var r = Evaluate(theta);
                return r.Ok ? r.Criterion : double.PositiveInfinity;
            }

            /// <summary>lmerTest's <c>devfun_vp</c>: the REML deviance with σ as a free parameter.</summary>
            private double DevianceInPhi(double[] phi)
            {
                var r = Evaluate(phi[..^1]);
                double s2 = phi[^1] * phi[^1];
                return r.LogDetL + r.LogDetRx + r.R2 / s2 + (_m - _p) * Math.Log(2 * Math.PI * s2);
            }

            /// <summary>lmerTest's <c>get_covbeta</c>, common block only, row-major.</summary>
            private double[] CovarianceInPhi(double[] phi)
            {
                var r = Evaluate(phi[..^1]);
                double s2 = phi[^1] * phi[^1];
                var c = new double[_common * _common];
                for (int i = 0; i < _common; i++)
                    for (int j = 0; j < _common; j++)
                        c[i * _common + j] = r.Ok ? s2 * r.AInverse![i, j] : double.NaN;
                return c;
            }

            /// <summary>Whether <paramref name="candidate"/> is no worse than <paramref name="best"/> beyond rounding (1e-10 relative).</summary>
            private static bool AtLeastAsGood(double candidate, double best) =>
                candidate <= best + 1e-10 * Math.Max(1, Math.Abs(best));

            private static readonly double[] Grid = { 0, 1e-3, 1e-2, 3e-2, 0.1, 0.2, 0.4, 0.7, 1, 1.5, 2.5, 4, 7, 12, 25, 60, 150 };

            /// <summary>θ minimising the REML criterion; null when the simplex search did not converge.</summary>
            private double[]? Minimise(int k)
            {
                if (k == 1)
                {
                    int bestI = 0;
                    double bestF = double.PositiveInfinity;
                    for (int i = 0; i < Grid.Length; i++)
                    {
                        double v = Criterion(new[] { Grid[i] });
                        if (v < bestF) { bestF = v; bestI = i; }
                    }
                    double lo = bestI == 0 ? 0 : Grid[bestI - 1], hi = bestI == Grid.Length - 1 ? Grid[^1] * 4 : Grid[bestI + 1];
                    var t = Polish(new[] { MixedModel.Brent(z => Criterion(new[] { z }), lo, hi, 1e-12) });
                    return AtLeastAsGood(Criterion(new[] { 0.0 }), Criterion(t)) ? new[] { 0.0 } : t;
                }

                var start = new[] { 1.0, 1.0 };
                double startF = double.PositiveInfinity;
                foreach (double a in Grid)
                    foreach (double b in Grid)
                    {
                        double v = Criterion(new[] { a, b });
                        if (v < startF) { startF = v; start = new[] { a, b }; }
                    }
                var simplexMinimum = NelderMead(start);
                if (simplexMinimum is null) return null;
                var candidates = new List<double[]> { Polish(simplexMinimum) };
                // Each boundary, explicitly: one component at zero, the other optimised in one dimension.
                for (int fixedAt = 0; fixedAt < 2; fixedAt++)
                {
                    int free = 1 - fixedAt;
                    double Along(double z) { var v = new double[2]; v[free] = z; return Criterion(v); }
                    double t = MixedModel.Brent(Along, 0, Grid[^1] * 4, 1e-12);
                    var c = new double[2]; c[free] = AtLeastAsGood(Along(0), Along(t)) ? 0 : t;
                    candidates.Add(c);
                }
                candidates.Add(new[] { 0.0, 0.0 });
                // Fewer non-zero components wins a tie: within rounding, a boundary is kept (as MixedModel keeps θ = 0).
                var ranked = candidates.Select(c => (c, f: Criterion(c), zeros: c.Count(z => z == 0))).ToList();
                double lowest = ranked.Min(r => r.f);
                return ranked.Where(r => AtLeastAsGood(r.f, lowest)).OrderByDescending(r => r.zeros).ThenBy(r => r.f).First().c;
            }

            /// <summary>
            /// The simplex stops when the criterion at its best and worst vertices agrees to this relative difference, twice
            /// running (MathNet's rule). It must sit well above the criterion's rounding noise, about 1e-14 relative: at 1e-15
            /// a simplex collapsed to adjacent floating-point θ cycled until the evaluation limit. <see cref="Polish"/>
            /// finishes from the simplex's minimum.
            /// </summary>
            private const double SimplexTolerance = 1e-10;
            /// <summary>
            /// The most criterion evaluations the simplex may take (the reference fixtures need at most 64); a search still
            /// going then is not converged.
            /// </summary>
            private const int SimplexEvaluations = 4000;

            /// <summary>
            /// MathNet's Nelder–Mead from <paramref name="start"/>, first steps max(0.1, 0.2·θᵢ); θ is clamped at 0 inside the
            /// criterion and a θ the criterion cannot be evaluated at counts as the worst value. Null when the search has not
            /// converged within <see cref="SimplexEvaluations"/>.
            /// </summary>
            private double[]? NelderMead(double[] start)
            {
                double F(Vector<double> v)
                {
                    double c = Criterion(v.Select(z => Math.Max(0, z)).ToArray());
                    return double.IsFinite(c) ? c : double.MaxValue;
                }
                var step = Vector<double>.Build.Dense(start.Length, i => Math.Max(0.1, 0.2 * start[i]));
                try
                {
                    var result = NelderMeadSimplex.Minimum(ObjectiveFunction.Value(F), Vector<double>.Build.DenseOfArray(start),
                        step, SimplexTolerance, SimplexEvaluations);
                    return result.MinimizingPoint.Select(z => Math.Max(0, z)).ToArray();
                }
                catch (MaximumIterationsException)
                {
                    return null;
                }
            }

            /// <summary>Newton steps on the interior components (Richardson derivatives), each accepted only if it lowers the criterion.</summary>
            private double[] Polish(double[] theta)
            {
                var current = (double[])theta.Clone();
                double fc = Criterion(current);
                for (int iteration = 0; iteration < 30; iteration++)
                {
                    var free = Enumerable.Range(0, current.Length).Where(i => current[i] > 1e-8).ToArray();
                    if (free.Length == 0) break;
                    double Sub(double[] z) { var v = (double[])current.Clone(); for (int i = 0; i < free.Length; i++) v[free[i]] = z[i]; return Criterion(v); }
                    var z0 = free.Select(i => current[i]).ToArray();
                    var g = Enumerable.Range(0, free.Length).Select(i => Richardson.Derivative(v => new[] { Sub(v) }, z0, i)[0]).ToArray();
                    var hm = Matrix<double>.Build.DenseOfArray(Richardson.Hessian(Sub, z0));
                    Vector<double> step;
                    try { step = hm.Cholesky().Solve(Vector<double>.Build.DenseOfArray(g)); }
                    catch (Exception e) when (e is ArgumentException or InvalidOperationException) { break; }
                    bool improved = false;
                    for (double s = 1; s > 1e-6; s /= 2)
                    {
                        var trial = (double[])current.Clone();
                        for (int i = 0; i < free.Length; i++) trial[free[i]] = Math.Max(0, current[free[i]] - s * step[i]);
                        double ft = Criterion(trial);
                        if (ft < fc) { current = trial; improved = fc - ft > 1e-15 * Math.Abs(fc); fc = ft; break; }
                    }
                    if (!improved || step.AbsoluteMaximum() < 1e-12) break;
                }
                return current;
            }
        }
    }

    /// <summary>
    /// Richardson-extrapolated central differences with numDeriv's defaults (4 halvings; relative step 1e-4 for a
    /// derivative and 0.1 for a Hessian, plus 1e-4 where |x| is near zero), the derivatives lmerTest's Satterthwaite
    /// method takes.
    /// </summary>
    internal static class Richardson
    {
        private const double RelativeStep = 1e-4, HessianRelativeStep = 0.1, AbsoluteStep = 1e-4, ZeroTolerance = 1.7818416954376e-5;
        private const int Levels = 4;

        private static double Step(double x) => Math.Abs(RelativeStep * x) + (Math.Abs(x) < ZeroTolerance ? AbsoluteStep : 0);

        /// <summary>∂f/∂x_k, elementwise for a vector-valued f.</summary>
        public static double[] Derivative(Func<double[], double[]> f, double[] x, int k)
        {
            double h0 = Step(x[k]);
            var table = new double[Levels][];
            for (int j = 0; j < Levels; j++)
            {
                double h = h0 / Math.Pow(2, j);
                var up = (double[])x.Clone(); up[k] += h;
                var down = (double[])x.Clone(); down[k] -= h;
                var fu = f(up); var fd = f(down);
                table[j] = fu.Zip(fd, (a, b) => (a - b) / (2 * h)).ToArray();
            }
            return Extrapolate(table);
        }

        /// <summary>
        /// The Hessian of a scalar f exactly as numDeriv's <c>hessian</c> (<c>genD</c>) computes it, because lmerTest's
        /// Satterthwaite df inherit its truncation error: relative step 0.1 (not the 1e-4 of the gradient), 4 halvings,
        /// diagonals from second differences, off-diagonals from (f(x+hᵢ+hⱼ) − 2f(x) + f(x−hᵢ−hⱼ) − Hᵢᵢhᵢ² − Hⱼⱼhⱼ²)/(2hᵢhⱼ)
        /// with the extrapolated diagonals, each Richardson-extrapolated.
        /// </summary>
        public static double[,] Hessian(Func<double[], double> f, double[] x)
        {
            int n = x.Length;
            double f0 = f(x);
            var h0 = x.Select(v => Math.Abs(HessianRelativeStep * v) + (Math.Abs(v) < ZeroTolerance ? AbsoluteStep : 0)).ToArray();
            var diagonal = new double[n];
            for (int i = 0; i < n; i++)
            {
                var table = new double[Levels][];
                double h = h0[i];
                for (int k = 0; k < Levels; k++)
                {
                    var up = (double[])x.Clone(); up[i] += h;
                    var down = (double[])x.Clone(); down[i] -= h;
                    table[k] = new[] { (f(up) - 2 * f0 + f(down)) / (h * h) };
                    h /= 2;
                }
                diagonal[i] = Extrapolate(table)[0];
            }
            var hessian = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                hessian[i, i] = diagonal[i];
                for (int j = 0; j < i; j++)
                {
                    var table = new double[Levels][];
                    double hi = h0[i], hj = h0[j];
                    for (int k = 0; k < Levels; k++)
                    {
                        var up = (double[])x.Clone(); up[i] += hi; up[j] += hj;
                        var down = (double[])x.Clone(); down[i] -= hi; down[j] -= hj;
                        table[k] = new[] { (f(up) - 2 * f0 + f(down) - diagonal[i] * hi * hi - diagonal[j] * hj * hj) / (2 * hi * hj) };
                        hi /= 2; hj /= 2;
                    }
                    hessian[i, j] = hessian[j, i] = Extrapolate(table)[0];
                }
            }
            return hessian;
        }
        private static double[] Extrapolate(double[][] table)
        {
            var t = table.Select(r => (double[])r.Clone()).ToArray();
            for (int m = 1; m < t.Length; m++)
            {
                double factor = Math.Pow(4, m);
                for (int j = t.Length - 1; j >= m; j--)
                    for (int e = 0; e < t[j].Length; e++)
                        t[j][e] = (factor * t[j][e] - t[j - 1][e]) / (factor - 1);
            }
            return t[^1];
        }
    }
}
