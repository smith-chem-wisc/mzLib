#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Fitting;

public sealed record GlobalAbundanceRefitOptions(
    int MaxIterations = 8,
    double ConvergenceTolerance = 1e-3,
    double MinimumAbundance = 0.0,
    bool Verbose = false,
    long MaxBasisCacheBytes = 1L << 30,
    double BasisSparsityThreshold = 1e-10);

/// <summary>
/// How much of the observed signal the fitted abundances fail to account for, and how much signal
/// they invent, reported separately.
/// </summary>
/// <param name="UnexplainedFraction">
/// Σ(observed − predicted)² / Σ observed², over <b>distinct experimental peaks</b>. Overlapping
/// proteoforms are matched to the same physical peak once each; counting that peak once per
/// claimant inflates both sums by an amount that depends on how crowded the region is, making the
/// number incomparable between runs.
/// </param>
/// <param name="OverpredictedFraction">
/// Σ predicted² over samples where no experimental peak matched at all, on the same denominator.
/// This is signal the model puts somewhere the instrument saw nothing — the opposite failure to
/// <paramref name="UnexplainedFraction"/>, with the opposite fix, so the two are not folded into
/// one ratio.
/// <para>
/// <b>Unlike <paramref name="UnexplainedFraction"/>, this numerator is summed over samples, not
/// over distinct positions.</b> A sample that matched nothing has no peak index to key on, so
/// there is nothing to deduplicate it by; two proteoforms probing the same empty neighbourhood
/// each contribute. It therefore still scales with proteoform crowding and should be read as a
/// within-run diagnostic rather than compared across runs of differing density.
/// </para>
/// </param>
/// <param name="UnmatchedSamples">Count of samples, not of distinct empty positions.</param>
public sealed record ResidualEnergySummary(
    double UnexplainedFraction,
    double OverpredictedFraction,
    int DistinctPeaks,
    int UnmatchedSamples);

public sealed record GlobalAbundanceRefitResult(
    FittedProteoform[] FittedProteoforms,
    int IterationsCompleted,
    bool Converged,
    double InitialResidualFraction,
    double FinalResidualFraction,
    IReadOnlyList<double> ResidualFractionByIteration,
    ResidualEnergySummary? InitialResiduals = null,
    ResidualEnergySummary? FinalResiduals = null);

/// <summary>
/// Refines per-proteoform abundance values jointly over potentially overlapping
/// signal using a non-negative coordinate-descent least-squares update.
/// Shape terms (mass, RT profile, charge distribution, sigma) are kept fixed.
/// </summary>
/// <remarks>
/// The sweep is Gauss-Seidel: updating proteoform <c>p</c> immediately publishes its new abundance
/// into the predicted totals, so proteoform <c>p + 1</c> reasons about a residual that already
/// accounts for it. Coordinate descent is only coordinate descent if that holds — with totals
/// frozen for the length of a sweep the update is Jacobi instead, which converges more slowly and
/// can oscillate on strongly coupled coordinates. Co-eluting, mass-adjacent proteoforms are exactly
/// that case.
/// </remarks>
public sealed class GlobalAbundanceRefitter
{
    private const double BasisInclusionThreshold = 1e-12;

    /// <summary>Marks a sample where the extractor found no experimental peak at all.</summary>
    private const long UnmatchedSourceKey = long.MinValue;

    /// <summary>
    /// High bit reserved for the synthetic keys handed to truths that carry no peak identity, so
    /// they can never collide with a real (scan, peak) key and nothing is deduplicated by accident.
    /// </summary>
    private const long SyntheticSourceKeyBase = 1L << 62;

    /// <param name="SourceKeys">
    /// Identity of the experimental peak behind each sample. Samples sharing a key across sample
    /// sets are the same physical peak claimed by different proteoforms.
    /// </param>
    /// <param name="SourceDistances">
    /// |observed m/z − this claimant's theoretical centroid|. Each claimant recorded the sample at
    /// its <i>own</i> centroid, which is up to a tolerance width away, so this is what says whose
    /// sample sits closest to the real peak.
    /// </param>
    private sealed record SampleSet(
        double[] Times,
        double[] Mzs,
        double[] Observed,
        double[] Basis,
        long[] SourceKeys,
        double[] SourceDistances);

    /// <summary>
    /// Compressed-row storage of the (sample x model) basis matrix for one sample set:
    /// entries for sample <c>k</c> are <c>[RowStart[k], RowStart[k + 1])</c>, model indices
    /// ascending so summation order matches the dense formulation.
    /// </summary>
    private sealed record SparseBasis(int[] RowStart, int[] ModelIndices, double[] Values);

    /// <summary>
    /// The transpose of <see cref="SparseBasis"/>: for each model, every (sample, basis value) pair
    /// it contributes to, over all sample sets at once. Samples are addressed by a flat index —
    /// <c>SampleSetOffsets[p] + k</c> is sample <c>k</c> of set <c>p</c> — so one array of totals
    /// covers the whole problem.
    /// </summary>
    /// <remarks>
    /// Built from the entries of the forward matrices rather than re-derived, so the two cannot
    /// disagree about which terms survived <see cref="GlobalAbundanceRefitOptions.BasisSparsityThreshold"/>.
    /// This is what lets a coordinate update publish its change in time proportional to the number
    /// of samples that model actually touches, instead of a full recompute per coordinate.
    /// </remarks>
    private sealed record BasisTranspose(
        int[] SampleSetOffsets,
        int[] ModelRowStart,
        int[] SampleIndices,
        double[] Values);

    private readonly GlobalAbundanceRefitOptions _options;

    public GlobalAbundanceRefitter(GlobalAbundanceRefitOptions? options = null)
    {
        _options = options ?? new GlobalAbundanceRefitOptions();
        if (_options.MaxIterations < 1)
            throw new ArgumentOutOfRangeException(nameof(options), "MaxIterations must be at least 1.");
        if (_options.ConvergenceTolerance <= 0)
            throw new ArgumentOutOfRangeException(nameof(options), "ConvergenceTolerance must be positive.");
    }

    public GlobalAbundanceRefitResult Refit(
        IReadOnlyList<FittedProteoform> fitted,
        IReadOnlyList<ProteoformGroundTruth> truths,
        int minCharge,
        int maxCharge,
        double sigmaMz) =>
        Refit(fitted, truths, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz));

    public GlobalAbundanceRefitResult Refit(
        IReadOnlyList<FittedProteoform> fitted,
        IReadOnlyList<ProteoformGroundTruth> truths,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel)
    {
        if (fitted.Count != truths.Count)
            throw new ArgumentException("fitted and truths must be parallel arrays with identical length.");
        if (minCharge < 1 || maxCharge < minCharge)
            throw new ArgumentException("Charge range must satisfy 1 <= minCharge <= maxCharge.");
        if (widthModel is null)
            throw new ArgumentNullException(nameof(widthModel));

        if (fitted.Count == 0)
        {
            return new GlobalAbundanceRefitResult(
                Array.Empty<FittedProteoform>(),
                IterationsCompleted: 0,
                Converged: true,
                InitialResidualFraction: 0,
                FinalResidualFraction: 0,
                ResidualFractionByIteration: Array.Empty<double>(),
                InitialResiduals: new ResidualEnergySummary(0, 0, 0, 0),
                FinalResiduals: new ResidualEnergySummary(0, 0, 0, 0));
        }

        var models = fitted.Select(f => f.Model).ToArray();
        var kernels = models.Select(m => new IsotopeEnvelopeKernel(m.MonoisotopicMass)).ToArray();

        var sampleSets = BuildSampleSets(models, truths, kernels, minCharge, maxCharge, widthModel);
        var sampleSetOffsets = BuildSampleSetOffsets(sampleSets);

        // The basis value of model m at sample k depends only on shape terms — mass, RT profile,
        // charge distribution, sigma — none of which this refit touches. Only Abundance changes, so
        // the whole (sample x model) matrix is invariant and is computed once instead of on every
        // one of the MaxIterations + 2 passes below. Cached values are the same numbers the
        // uncached path would recompute, so caching alone does not change results.
        var basisMatrices = TryBuildBasisMatrices(
            models, kernels, sampleSets, minCharge, maxCharge, widthModel,
            _options.MaxBasisCacheBytes, _options.BasisSparsityThreshold);

        // Model-major view of the same entries. Null only when the basis itself could not be cached,
        // in which case the sweep falls back to recomputing one sample set's totals at a time.
        var transpose = basisMatrices is null
            ? null
            : BuildTranspose(basisMatrices, sampleSets, models.Length);
        double[]? liveTotals = transpose is null
            ? null
            : new double[transpose.SampleSetOffsets[sampleSets.Length]];

        var initialPredictedTotals = ComputePredictedTotals(models, kernels, sampleSets, basisMatrices, minCharge, maxCharge, widthModel, _options.BasisSparsityThreshold);
        var initialResiduals = ComputeResidualEnergy(sampleSets, Flatten(initialPredictedTotals, sampleSetOffsets), sampleSetOffsets);
        var residualByIteration = new List<double>(_options.MaxIterations);
        int completedIterations = 0;
        bool converged = false;

        for (int iteration = 0; iteration < _options.MaxIterations; iteration++)
        {
            // Rebuilt from the current abundances rather than carried over, so rounding error from
            // the previous sweep's incremental updates cannot accumulate across sweeps.
            if (transpose is not null)
                RefreshLiveTotals(models, transpose, liveTotals!);

            double maxRelativeChange = 0;

            for (int p = 0; p < models.Length; p++)
            {
                double oldAbundance = models[p].Abundance;

                double[] totals;
                int totalsOffset;
                if (transpose is not null)
                {
                    totals = liveTotals!;
                    totalsOffset = transpose.SampleSetOffsets[p];
                }
                else
                {
                    // Recomputed immediately before use, which is what makes this branch
                    // Gauss-Seidel too. Over a whole sweep it touches each sample set exactly once,
                    // so it costs no more than the single up-front pass it replaces.
                    totals = ComputeSampleSetTotals(models, kernels, sampleSets[p], minCharge, maxCharge, widthModel, _options.BasisSparsityThreshold);
                    totalsOffset = 0;
                }

                double updatedAbundance = SolveCoordinateUpdate(
                    sampleSet: sampleSets[p],
                    predictedTotals: totals,
                    totalsOffset: totalsOffset,
                    oldAbundance,
                    _options.MinimumAbundance);

                if (oldAbundance > 0)
                    maxRelativeChange = Math.Max(maxRelativeChange, Math.Abs(updatedAbundance - oldAbundance) / oldAbundance);
                else if (updatedAbundance > 0)
                    maxRelativeChange = 1;

                models[p] = models[p] with { Abundance = updatedAbundance };

                if (transpose is not null)
                    ApplyAbundanceDelta(transpose, p, updatedAbundance - oldAbundance, liveTotals!);
            }

            completedIterations++;

            double[] sweepTotals = transpose is not null
                ? liveTotals!
                : Flatten(
                    ComputePredictedTotals(models, kernels, sampleSets, basisMatrices, minCharge, maxCharge, widthModel, _options.BasisSparsityThreshold),
                    sampleSetOffsets);
            double sweepResidual = ComputeResidualEnergy(sampleSets, sweepTotals, sampleSetOffsets).UnexplainedFraction;
            residualByIteration.Add(sweepResidual);

            if (_options.Verbose)
                Console.WriteLine($"Global abundance refit iteration {completedIterations}/{_options.MaxIterations}: max relative change = {maxRelativeChange:G4}, residual fraction = {sweepResidual:G6}");

            if (maxRelativeChange <= _options.ConvergenceTolerance)
            {
                converged = true;
                break;
            }
        }

        var finalPredictedTotals = ComputePredictedTotals(models, kernels, sampleSets, basisMatrices, minCharge, maxCharge, widthModel, _options.BasisSparsityThreshold);
        var finalResiduals = ComputeResidualEnergy(sampleSets, Flatten(finalPredictedTotals, sampleSetOffsets), sampleSetOffsets);
        var updatedFits = new FittedProteoform[fitted.Count];
        for (int i = 0; i < fitted.Count; i++)
        {
            updatedFits[i] = fitted[i] with { Model = models[i] };
        }

        return new GlobalAbundanceRefitResult(
            updatedFits,
            completedIterations,
            converged,
            initialResiduals.UnexplainedFraction,
            finalResiduals.UnexplainedFraction,
            residualByIteration,
            initialResiduals,
            finalResiduals);
    }

    private static int[] BuildSampleSetOffsets(IReadOnlyList<SampleSet> sampleSets)
    {
        var offsets = new int[sampleSets.Count + 1];
        for (int p = 0; p < sampleSets.Count; p++)
            offsets[p + 1] = offsets[p] + sampleSets[p].Observed.Length;

        return offsets;
    }

    private static double[] Flatten(IReadOnlyList<double[]> jagged, int[] offsets)
    {
        var flat = new double[offsets[^1]];
        for (int p = 0; p < jagged.Count; p++)
            Array.Copy(jagged[p], 0, flat, offsets[p], jagged[p].Length);

        return flat;
    }

    /// <summary>
    /// Groups the entries of <paramref name="matrices"/> by model instead of by sample. Each
    /// model's row comes out ordered by ascending flat sample index.
    /// </summary>
    /// <remarks>
    /// Agreement with the forward pass's summation order does <b>not</b> come from this method — it
    /// comes from <see cref="RefreshLiveTotals"/> accumulating models in ascending order, which is
    /// what makes each sample's contributions land in the same sequence the sample-major pass adds
    /// them in. Reordering or parallelising that loop would break the bit-agreement the tests lean
    /// on, however this method is written.
    /// </remarks>
    private static BasisTranspose BuildTranspose(
        SparseBasis[] matrices,
        IReadOnlyList<SampleSet> sampleSets,
        int modelCount)
    {
        var offsets = new int[sampleSets.Count + 1];
        for (int p = 0; p < sampleSets.Count; p++)
            offsets[p + 1] = offsets[p] + sampleSets[p].Observed.Length;

        var rowStart = new int[modelCount + 1];
        foreach (var matrix in matrices)
        {
            foreach (int m in matrix.ModelIndices)
                rowStart[m + 1]++;
        }

        for (int m = 0; m < modelCount; m++)
            rowStart[m + 1] += rowStart[m];

        var sampleIndices = new int[rowStart[modelCount]];
        var values = new double[rowStart[modelCount]];
        var cursor = (int[])rowStart.Clone();

        for (int p = 0; p < matrices.Length; p++)
        {
            var matrix = matrices[p];
            int setOffset = offsets[p];
            int sampleCount = sampleSets[p].Observed.Length;
            for (int k = 0; k < sampleCount; k++)
            {
                int end = matrix.RowStart[k + 1];
                for (int e = matrix.RowStart[k]; e < end; e++)
                {
                    int slot = cursor[matrix.ModelIndices[e]]++;
                    sampleIndices[slot] = setOffset + k;
                    values[slot] = matrix.Values[e];
                }
            }
        }

        return new BasisTranspose(offsets, rowStart, sampleIndices, values);
    }

    private static void RefreshLiveTotals(
        IReadOnlyList<ProteoformModel> models,
        BasisTranspose transpose,
        double[] liveTotals)
    {
        Array.Clear(liveTotals);
        for (int m = 0; m < models.Count; m++)
        {
            double abundance = models[m].Abundance;
            int end = transpose.ModelRowStart[m + 1];
            for (int e = transpose.ModelRowStart[m]; e < end; e++)
                liveTotals[transpose.SampleIndices[e]] += abundance * transpose.Values[e];
        }
    }

    private static void ApplyAbundanceDelta(
        BasisTranspose transpose,
        int model,
        double delta,
        double[] liveTotals)
    {
        if (delta == 0)
            return;

        int end = transpose.ModelRowStart[model + 1];
        for (int e = transpose.ModelRowStart[model]; e < end; e++)
            liveTotals[transpose.SampleIndices[e]] += delta * transpose.Values[e];
    }

    private static double SolveCoordinateUpdate(
        SampleSet sampleSet,
        double[] predictedTotals,
        int totalsOffset,
        double oldAbundance,
        double minimumAbundance)
    {
        if (sampleSet.Basis.Length == 0)
            return Math.Max(minimumAbundance, oldAbundance);

        double numerator = 0;
        double denominator = 0;

        for (int k = 0; k < sampleSet.Basis.Length; k++)
        {
            double basis = sampleSet.Basis[k];
            double observed = sampleSet.Observed[k];
            double residualExcludingCurrent = observed - (predictedTotals[totalsOffset + k] - oldAbundance * basis);
            numerator += basis * residualExcludingCurrent;
            denominator += basis * basis;
        }

        if (denominator <= 0)
            return Math.Max(minimumAbundance, oldAbundance);

        double updated = numerator / denominator;
        if (double.IsNaN(updated) || double.IsInfinity(updated))
            return Math.Max(minimumAbundance, oldAbundance);

        return Math.Max(minimumAbundance, updated);
    }

    private static SampleSet[] BuildSampleSets(
        IReadOnlyList<ProteoformModel> models,
        IReadOnlyList<ProteoformGroundTruth> truths,
        IReadOnlyList<IsotopeEnvelopeKernel> kernels,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel)
    {
        var sampleSets = new SampleSet[truths.Count];

        for (int p = 0; p < truths.Count; p++)
        {
            var truth = truths[p];
            var model = models[p];
            var kernel = kernels[p];
            bool hasIdentity = truth.HasSourcePeakIdentity;

            var times = new List<double>();
            var mzs = new List<double>();
            var observed = new List<double>();
            var basis = new List<double>();
            var sourceKeys = new List<long>();
            var sourceDistances = new List<double>();

            for (int c = 0; c < truth.ChargeCount; c++)
            {
                var chargeMzs = truth.CentroidMzs[c];
                for (int i = 0; i < chargeMzs.Length; i++)
                {
                    double mz = chargeMzs[i];
                    for (int s = 0; s < truth.ScanCount; s++)
                    {
                        double b = EvaluateUnitContribution(model, kernel, truth.ScanTimes[s], mz, minCharge, maxCharge, widthModel);
                        if (b <= BasisInclusionThreshold)
                            continue;

                        times.Add(truth.ScanTimes[s]);
                        mzs.Add(mz);
                        observed.Add(truth.IsotopologueIntensities[c][i][s]);
                        basis.Add(b);

                        if (hasIdentity)
                        {
                            int peakIndex = truth.SourcePeakIndices[c][i][s];
                            sourceKeys.Add(peakIndex < 0
                                ? UnmatchedSourceKey
                                : ((long)truth.ZeroBasedScanIndices[s] << 32) | (uint)peakIndex);
                            sourceDistances.Add(truth.SourcePeakDeltas[c][i][s]);
                        }
                        else
                        {
                            // No identity to dedup on, so give each sample a key of its own and let
                            // the accounting behave as it did before peak identity existed.
                            sourceKeys.Add(SyntheticSourceKeyBase + ((long)p << 32) + sourceKeys.Count);
                            sourceDistances.Add(double.NaN);
                        }
                    }
                }
            }

            sampleSets[p] = new SampleSet(
                times.ToArray(), mzs.ToArray(), observed.ToArray(), basis.ToArray(),
                sourceKeys.ToArray(), sourceDistances.ToArray());
        }

        return sampleSets;
    }

    /// <summary>
    /// Precomputes the (sample x model) basis matrix for every sample set in compressed-row form.
    /// Returns null when it would exceed <paramref name="maxBytes"/>, in which case callers fall
    /// back to recomputing basis values on demand — slower, but the same arithmetic, because the
    /// fallback applies the identical sparsity cutoff.
    /// </summary>
    private static SparseBasis[]? TryBuildBasisMatrices(
        IReadOnlyList<ProteoformModel> models,
        IReadOnlyList<IsotopeEnvelopeKernel> kernels,
        IReadOnlyList<SampleSet> sampleSets,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        long maxBytes,
        double sparsityThreshold)
    {
        const int bytesPerEntry = sizeof(int) + sizeof(double);
        int modelCount = models.Count;
        long budgetEntries = maxBytes / bytesPerEntry;
        long usedEntries = 0;

        var matrices = new SparseBasis[sampleSets.Count];
        var modelIndices = new List<int>();
        var values = new List<double>();

        for (int p = 0; p < sampleSets.Count; p++)
        {
            var set = sampleSets[p];
            int sampleCount = set.Observed.Length;
            var rowStart = new int[sampleCount + 1];

            modelIndices.Clear();
            values.Clear();

            for (int k = 0; k < sampleCount; k++)
            {
                rowStart[k] = modelIndices.Count;
                double time = set.Times[k];
                double mz = set.Mzs[k];
                double cutoff = SparsityCutoff(set.Basis[k], sparsityThreshold);

                for (int m = 0; m < modelCount; m++)
                {
                    double basis = EvaluateUnitContribution(models[m], kernels[m], time, mz, minCharge, maxCharge, widthModel);
                    if (basis <= 0 || basis < cutoff)
                        continue;

                    modelIndices.Add(m);
                    values.Add(basis);
                }

                if (usedEntries + modelIndices.Count > budgetEntries)
                    return null;
            }

            rowStart[sampleCount] = modelIndices.Count;
            usedEntries += modelIndices.Count;
            matrices[p] = new SparseBasis(rowStart, modelIndices.ToArray(), values.ToArray());
        }

        return matrices;
    }

    /// <summary>
    /// Terms below this fraction of the sample's own-model basis are dropped. Scaling by the
    /// own-model basis makes the cutoff relative to the signal actually present at that sample
    /// rather than an absolute intensity, so it behaves the same across abundance ranges.
    /// A threshold of zero keeps every positive term.
    /// </summary>
    private static double SparsityCutoff(double ownBasis, double sparsityThreshold) =>
        sparsityThreshold > 0 ? sparsityThreshold * ownBasis : 0;

    private static double[][] ComputePredictedTotals(
        IReadOnlyList<ProteoformModel> models,
        IReadOnlyList<IsotopeEnvelopeKernel> kernels,
        IReadOnlyList<SampleSet> sampleSets,
        SparseBasis[]? basisMatrices,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        double sparsityThreshold = 0)
    {
        var predicted = new double[sampleSets.Count][];

        for (int p = 0; p < sampleSets.Count; p++)
        {
            var sparse = basisMatrices?[p];
            predicted[p] = sparse is not null
                ? ComputeSampleSetTotalsFromSparse(models, sampleSets[p], sparse)
                : ComputeSampleSetTotals(models, kernels, sampleSets[p], minCharge, maxCharge, widthModel, sparsityThreshold);
        }

        return predicted;
    }

    private static double[] ComputeSampleSetTotalsFromSparse(
        IReadOnlyList<ProteoformModel> models,
        SampleSet set,
        SparseBasis sparse)
    {
        var totals = new double[set.Observed.Length];
        for (int k = 0; k < totals.Length; k++)
        {
            double sum = 0;
            int end = sparse.RowStart[k + 1];
            for (int e = sparse.RowStart[k]; e < end; e++)
                sum += models[sparse.ModelIndices[e]].Abundance * sparse.Values[e];

            totals[k] = sum;
        }

        return totals;
    }

    /// <summary>
    /// Predicted totals for one sample set with basis values evaluated on demand. Applies the same
    /// sparsity cutoff as <see cref="TryBuildBasisMatrices"/>, so it is the same arithmetic the
    /// cached path would have produced — only slower.
    /// </summary>
    private static double[] ComputeSampleSetTotals(
        IReadOnlyList<ProteoformModel> models,
        IReadOnlyList<IsotopeEnvelopeKernel> kernels,
        SampleSet set,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        double sparsityThreshold)
    {
        int modelCount = models.Count;
        var totals = new double[set.Observed.Length];

        for (int k = 0; k < totals.Length; k++)
        {
            double time = set.Times[k];
            double mz = set.Mzs[k];
            double cutoff = SparsityCutoff(set.Basis[k], sparsityThreshold);

            double sum = 0;
            for (int m = 0; m < modelCount; m++)
            {
                double basis = EvaluateUnitContribution(models[m], kernels[m], time, mz, minCharge, maxCharge, widthModel);
                if (basis <= 0 || basis < cutoff)
                    continue;

                sum += models[m].Abundance * basis;
            }

            totals[k] = sum;
        }

        return totals;
    }

    private static double EvaluateUnitContribution(
        ProteoformModel model,
        IsotopeEnvelopeKernel kernel,
        double time,
        double mz,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel)
    {
        double rt = model.RtProfile.Evaluate(time);
        if (rt <= 0)
            return 0;

        double chargeSum = 0;
        for (int z = minCharge; z <= maxCharge; z++)
        {
            double fz = model.ChargeDistribution.Evaluate(z);
            if (fz <= 0)
                continue;

            chargeSum += fz * kernel.Evaluate(mz, z, widthModel);
        }

        return chargeSum > 0 ? rt * chargeSum : 0;
    }

    /// <summary>
    /// Energy accounting over <b>distinct experimental peaks</b> rather than over samples.
    /// </summary>
    /// <remarks>
    /// Sample sets are built one per proteoform from that proteoform's own centroid grid, so when
    /// two proteoforms overlap in m/z and RT the extractor matches the same experimental peak into
    /// both. Summing over sample sets independently counts that peak's intensity once per claimant,
    /// inflating the denominator by an amount that depends on how crowded the neighbourhood is. The
    /// result is a weighted average over per-proteoform neighbourhoods, not a fraction of the
    /// file's energy, and it is not comparable across runs with different degrees of overlap —
    /// which is the only comparison anyone would want to make with it.
    /// <para>
    /// Deduplication cannot key on the recorded m/z: each claimant stored the sample at its own
    /// theoretical centroid, and those differ by Δmass/z even though the source peak is identical.
    /// It keys on (scan, peak index) instead, and where several claimants share a peak it keeps the
    /// one whose theoretical centroid sits closest to where the peak actually is, that being the
    /// best available estimate of what the model predicts at the real peak position.
    /// </para>
    /// </remarks>
    private static ResidualEnergySummary ComputeResidualEnergy(
        IReadOnlyList<SampleSet> sampleSets,
        double[] flatTotals,
        int[] sampleSetOffsets)
    {
        var bestBySourcePeak = new Dictionary<long, (double Observed, double Predicted, double Distance)>();
        double overpredictedEnergy = 0;
        int unmatchedSamples = 0;

        for (int p = 0; p < sampleSets.Count; p++)
        {
            var set = sampleSets[p];
            int offset = sampleSetOffsets[p];
            for (int k = 0; k < set.Observed.Length; k++)
            {
                double predicted = flatTotals[offset + k];
                long key = set.SourceKeys[k];

                if (key == UnmatchedSourceKey)
                {
                    // Nothing was there to explain. Predicting signal here is a different failure
                    // from failing to explain signal that is there, so it is reported separately.
                    overpredictedEnergy += predicted * predicted;
                    unmatchedSamples++;
                    continue;
                }

                double distance = set.SourceDistances[k];
                if (!bestBySourcePeak.TryGetValue(key, out var existing) || distance < existing.Distance)
                    bestBySourcePeak[key] = (set.Observed[k], predicted, distance);
            }
        }

        double observedEnergy = 0;
        double residualEnergy = 0;
        foreach (var entry in bestBySourcePeak.Values)
        {
            double residual = entry.Observed - entry.Predicted;
            observedEnergy += entry.Observed * entry.Observed;
            residualEnergy += residual * residual;
        }

        return new ResidualEnergySummary(
            UnexplainedFraction: observedEnergy > 0 ? residualEnergy / observedEnergy : 0,
            OverpredictedFraction: observedEnergy > 0 ? overpredictedEnergy / observedEnergy : 0,
            DistinctPeaks: bestBySourcePeak.Count,
            UnmatchedSamples: unmatchedSamples);
    }
}
