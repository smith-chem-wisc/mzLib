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
    double FinalResidualFraction);

/// <summary>
/// Refines per-proteoform abundance values jointly over potentially overlapping
/// signal using a non-negative coordinate-descent least-squares update.
/// Shape terms (mass, RT profile, charge distribution, sigma) are kept fixed.
/// </summary>
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
                FinalResidualFraction: 0);
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

        var initialPredictedTotals = ComputePredictedTotals(models, kernels, sampleSets, basisMatrices, minCharge, maxCharge, sigmaMz, _options.BasisSparsityThreshold);
        double initialResidual = ComputeResidualEnergyFraction(sampleSets, initialPredictedTotals);
        int completedIterations = 0;
        bool converged = false;

        for (int iteration = 0; iteration < _options.MaxIterations; iteration++)
        {
            var predictedTotals = ComputePredictedTotals(models, kernels, sampleSets, basisMatrices, minCharge, maxCharge, sigmaMz, _options.BasisSparsityThreshold);
            double maxRelativeChange = 0;

            for (int p = 0; p < models.Length; p++)
            {
                double oldAbundance = models[p].Abundance;
                double updatedAbundance = SolveCoordinateUpdate(
                    sampleSet: sampleSets[p],
                    predictedTotals: predictedTotals[p],
                    oldAbundance,
                    _options.MinimumAbundance);

                if (oldAbundance > 0)
                    maxRelativeChange = Math.Max(maxRelativeChange, Math.Abs(updatedAbundance - oldAbundance) / oldAbundance);
                else if (updatedAbundance > 0)
                    maxRelativeChange = 1;

                models[p] = models[p] with { Abundance = updatedAbundance };
            }

            completedIterations++;
            if (_options.Verbose)
                Console.WriteLine($"Global abundance refit iteration {completedIterations}/{_options.MaxIterations}: max relative change = {maxRelativeChange:G4}");

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
            finalResidual);
    }

    private static double SolveCoordinateUpdate(
        SampleSet sampleSet,
        double[] predictedTotals,
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
            double residualExcludingCurrent = observed - (predictedTotals[k] - oldAbundance * basis);
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
        int modelCount = models.Count;

        for (int p = 0; p < sampleSets.Count; p++)
        {
            var set = sampleSets[p];
            var totals = new double[set.Observed.Length];
            var sparse = basisMatrices?[p];

            for (int k = 0; k < totals.Length; k++)
            {
                double sum = 0;

                if (sparse is not null)
                {
                    int end = sparse.RowStart[k + 1];
                    for (int e = sparse.RowStart[k]; e < end; e++)
                        sum += models[sparse.ModelIndices[e]].Abundance * sparse.Values[e];
                }
                else
                {
                    double time = set.Times[k];
                    double mz = set.Mzs[k];
                    double cutoff = SparsityCutoff(set.Basis[k], sparsityThreshold);

                    for (int m = 0; m < modelCount; m++)
                    {
                        double basis = EvaluateUnitContribution(models[m], kernels[m], time, mz, minCharge, maxCharge, sigmaMz);
                        if (basis <= 0 || basis < cutoff)
                            continue;

                        sum += models[m].Abundance * basis;
                    }
                }

                totals[k] = sum;
            }

            predicted[p] = totals;
        }

        return predicted;
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
}
