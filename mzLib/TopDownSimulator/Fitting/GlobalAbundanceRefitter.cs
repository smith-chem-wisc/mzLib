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

public sealed record GlobalAbundanceRefitResult(
    FittedProteoform[] FittedProteoforms,
    int IterationsCompleted,
    bool Converged,
    double InitialResidualFraction,
    double FinalResidualFraction,
    IReadOnlyList<double> ResidualFractionByIteration);

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

    private sealed record SampleSet(double[] Times, double[] Mzs, double[] Observed, double[] Basis);

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
        double sigmaMz)
    {
        if (fitted.Count != truths.Count)
            throw new ArgumentException("fitted and truths must be parallel arrays with identical length.");
        if (minCharge < 1 || maxCharge < minCharge)
            throw new ArgumentException("Charge range must satisfy 1 <= minCharge <= maxCharge.");
        if (sigmaMz <= 0)
            throw new ArgumentOutOfRangeException(nameof(sigmaMz));

        if (fitted.Count == 0)
        {
            return new GlobalAbundanceRefitResult(
                Array.Empty<FittedProteoform>(),
                IterationsCompleted: 0,
                Converged: true,
                InitialResidualFraction: 0,
                FinalResidualFraction: 0,
                ResidualFractionByIteration: Array.Empty<double>());
        }

        var models = fitted.Select(f => f.Model).ToArray();
        var kernels = models.Select(m => new IsotopeEnvelopeKernel(m.MonoisotopicMass)).ToArray();

        var sampleSets = BuildSampleSets(models, truths, kernels, minCharge, maxCharge, sigmaMz);

        // The basis value of model m at sample k depends only on shape terms — mass, RT profile,
        // charge distribution, sigma — none of which this refit touches. Only Abundance changes, so
        // the whole (sample x model) matrix is invariant and is computed once instead of on every
        // one of the MaxIterations + 2 passes below. Cached values are the same numbers the
        // uncached path would recompute, so caching alone does not change results.
        var basisMatrices = TryBuildBasisMatrices(
            models, kernels, sampleSets, minCharge, maxCharge, sigmaMz,
            _options.MaxBasisCacheBytes, _options.BasisSparsityThreshold);

        // Model-major view of the same entries. Null only when the basis itself could not be cached,
        // in which case the sweep falls back to recomputing one sample set's totals at a time.
        var transpose = basisMatrices is null
            ? null
            : BuildTranspose(basisMatrices, sampleSets, models.Length);
        double[]? liveTotals = transpose is null
            ? null
            : new double[transpose.SampleSetOffsets[sampleSets.Length]];

        var initialPredictedTotals = ComputePredictedTotals(models, kernels, sampleSets, basisMatrices, minCharge, maxCharge, sigmaMz, _options.BasisSparsityThreshold);
        double initialResidual = ComputeResidualEnergyFraction(sampleSets, initialPredictedTotals);
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
                    totals = ComputeSampleSetTotals(models, kernels, sampleSets[p], minCharge, maxCharge, sigmaMz, _options.BasisSparsityThreshold);
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

            double sweepResidual = transpose is not null
                ? ComputeResidualEnergyFraction(sampleSets, liveTotals!, transpose.SampleSetOffsets)
                : ComputeResidualEnergyFraction(
                    sampleSets,
                    ComputePredictedTotals(models, kernels, sampleSets, basisMatrices, minCharge, maxCharge, sigmaMz, _options.BasisSparsityThreshold));
            residualByIteration.Add(sweepResidual);

            if (_options.Verbose)
                Console.WriteLine($"Global abundance refit iteration {completedIterations}/{_options.MaxIterations}: max relative change = {maxRelativeChange:G4}, residual fraction = {sweepResidual:G6}");

            if (maxRelativeChange <= _options.ConvergenceTolerance)
            {
                converged = true;
                break;
            }
        }

        var finalPredictedTotals = ComputePredictedTotals(models, kernels, sampleSets, basisMatrices, minCharge, maxCharge, sigmaMz, _options.BasisSparsityThreshold);
        double finalResidual = ComputeResidualEnergyFraction(sampleSets, finalPredictedTotals);
        var updatedFits = new FittedProteoform[fitted.Count];
        for (int i = 0; i < fitted.Count; i++)
        {
            updatedFits[i] = fitted[i] with { Model = models[i] };
        }

        return new GlobalAbundanceRefitResult(
            updatedFits,
            completedIterations,
            converged,
            initialResidual,
            finalResidual,
            residualByIteration);
    }

    /// <summary>
    /// Groups the entries of <paramref name="matrices"/> by model instead of by sample. Iterating
    /// models in ascending order while filling means each sample's contributions land in ascending
    /// model order, matching the forward pass's summation order exactly.
    /// </summary>
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
        double sigmaMz)
    {
        var sampleSets = new SampleSet[truths.Count];

        for (int p = 0; p < truths.Count; p++)
        {
            var truth = truths[p];
            var model = models[p];
            var kernel = kernels[p];

            var times = new List<double>();
            var mzs = new List<double>();
            var observed = new List<double>();
            var basis = new List<double>();

            for (int c = 0; c < truth.ChargeCount; c++)
            {
                var chargeMzs = truth.CentroidMzs[c];
                for (int i = 0; i < chargeMzs.Length; i++)
                {
                    double mz = chargeMzs[i];
                    for (int s = 0; s < truth.ScanCount; s++)
                    {
                        double b = EvaluateUnitContribution(model, kernel, truth.ScanTimes[s], mz, minCharge, maxCharge, sigmaMz);
                        if (b <= BasisInclusionThreshold)
                            continue;

                        times.Add(truth.ScanTimes[s]);
                        mzs.Add(mz);
                        observed.Add(truth.IsotopologueIntensities[c][i][s]);
                        basis.Add(b);
                    }
                }
            }

            sampleSets[p] = new SampleSet(times.ToArray(), mzs.ToArray(), observed.ToArray(), basis.ToArray());
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
        double sigmaMz,
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
                    double basis = EvaluateUnitContribution(models[m], kernels[m], time, mz, minCharge, maxCharge, sigmaMz);
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
        double sigmaMz,
        double sparsityThreshold = 0)
    {
        var predicted = new double[sampleSets.Count][];

        for (int p = 0; p < sampleSets.Count; p++)
        {
            var sparse = basisMatrices?[p];
            predicted[p] = sparse is not null
                ? ComputeSampleSetTotalsFromSparse(models, sampleSets[p], sparse)
                : ComputeSampleSetTotals(models, kernels, sampleSets[p], minCharge, maxCharge, sigmaMz, sparsityThreshold);
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
        double sigmaMz,
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
                double basis = EvaluateUnitContribution(models[m], kernels[m], time, mz, minCharge, maxCharge, sigmaMz);
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
        double sigmaMz)
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

            chargeSum += fz * kernel.Evaluate(mz, z, sigmaMz);
        }

        return chargeSum > 0 ? rt * chargeSum : 0;
    }

    private static double ComputeResidualEnergyFraction(
        IReadOnlyList<SampleSet> sampleSets,
        IReadOnlyList<double[]> predictedTotals)
    {
        double observedEnergy = 0;
        double residualEnergy = 0;

        for (int p = 0; p < sampleSets.Count; p++)
        {
            var set = sampleSets[p];
            var totals = predictedTotals[p];
            for (int k = 0; k < set.Observed.Length; k++)
            {
                double observed = set.Observed[k];
                double residual = observed - totals[k];
                observedEnergy += observed * observed;
                residualEnergy += residual * residual;
            }
        }

        return observedEnergy > 0 ? residualEnergy / observedEnergy : 0;
    }

    /// <summary>
    /// Same quantity as the jagged overload, read off the flat live totals so a sweep can report
    /// its residual without a second pass over the basis.
    /// </summary>
    private static double ComputeResidualEnergyFraction(
        IReadOnlyList<SampleSet> sampleSets,
        double[] liveTotals,
        int[] sampleSetOffsets)
    {
        double observedEnergy = 0;
        double residualEnergy = 0;

        for (int p = 0; p < sampleSets.Count; p++)
        {
            var set = sampleSets[p];
            int offset = sampleSetOffsets[p];
            for (int k = 0; k < set.Observed.Length; k++)
            {
                double observed = set.Observed[k];
                double residual = observed - liveTotals[offset + k];
                observedEnergy += observed * observed;
                residualEnergy += residual * residual;
            }
        }

        return observedEnergy > 0 ? residualEnergy / observedEnergy : 0;
    }
}
