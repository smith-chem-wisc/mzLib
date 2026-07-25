#nullable enable
using System;
using System.Collections.Generic;
using Chemistry;
using Omics;
using Omics.Fragmentation;
using TopDownSimulator.Model;

namespace TopDownSimulator.Ms2;

/// <summary>
/// The per-bond arithmetic of the cascade, kept separate from ion assignment so it can be inspected
/// and tested on its own.
/// </summary>
/// <param name="BaseSequence">Sequence walked, N-terminus to C-terminus.</param>
/// <param name="StartingIonCount">R_0.</param>
/// <param name="CleavageProbabilities">p_i for each bond, length <c>BaseSequence.Length - 1</c>.</param>
/// <param name="RemainingBeforeBond">
/// R_i, the intensity still intact when bond i is reached. <c>RemainingBeforeBond[0] == R_0</c>.
/// </param>
/// <param name="CleavedIntensity">p_i * R_i, the intensity consumed by bond i.</param>
/// <param name="SurvivingPrecursorIntensity">
/// R_n after the last bond, i.e. <c>R_0 * Π (1 - p_i)</c>.
/// </param>
public sealed record BondCleavageCascade(
    string BaseSequence,
    double StartingIonCount,
    double[] CleavageProbabilities,
    double[] RemainingBeforeBond,
    double[] CleavedIntensity,
    double SurvivingPrecursorIntensity)
{
    /// <summary>Number of backbone bonds walked.</summary>
    public int BondCount => CleavageProbabilities.Length;
}

/// <summary>
/// Basic MS2 model: a sequential cleavage cascade along the backbone from the N-terminus, driven by
/// an <see cref="IBondCleavageModel"/>.
/// </summary>
/// <remarks>
/// <para>
/// Start with a fixed ion count <c>R_0</c> (<see cref="PropensityCascadeOptions.StartingIonCount"/>).
/// Walk the backbone from the N-terminus. At bond <c>i</c> with cleavage probability <c>p_i</c>:
/// </para>
/// <code>
/// cleaved_i = p_i * R_i
/// R_(i+1)   = R_i * (1 - p_i)
/// </code>
/// <para>
/// So an 11% bond takes 11% of whatever intensity is still intact when it is reached, and the
/// remaining 89% is carried forward to the next bond. Whatever survives the last bond stays as
/// intact precursor, which is why
/// <c>Σ cleaved_i + SurvivingPrecursorIntensity == R_0</c> exactly.
/// </para>
/// <para>
/// <c>cleaved_i</c> is the intensity of a <em>cleavage event</em>. It is then split between the
/// N-terminal ion and its complementary C-terminal ion by
/// <see cref="PropensityCascadeOptions.IonAssignment"/> /
/// <see cref="PropensityCascadeOptions.NTerminalIonFraction"/>, and each of those is split evenly
/// across <see cref="PropensityCascadeOptions.FragmentChargeStates"/>. Every split is a partition,
/// so conservation survives all of it.
/// </para>
/// <para>
/// Fragment masses and m/z are not computed here. The precursor's own
/// <c>IBioPolymerWithSetMods.Fragment</c> supplies mzLib <see cref="Product"/> objects (correct
/// terminal caps, product-type mass shifts and localized modifications), and
/// <c>Chemistry.ClassExtensions.ToMz</c> converts them.
/// </para>
/// </remarks>
public sealed class PropensityCascadeFragmentationModel : IMs2FragmentationModel
{
    public PropensityCascadeFragmentationModel(
        IBondCleavageModel? bondCleavageModel = null,
        PropensityCascadeOptions? options = null)
    {
        BondCleavageModel = bondCleavageModel ?? new ResiduePairCleavagePropensityModel();
        Options = options ?? new PropensityCascadeOptions();
        Options.Validate();
    }

    /// <summary>Where the per-bond cleavage probabilities come from.</summary>
    public IBondCleavageModel BondCleavageModel { get; }

    public PropensityCascadeOptions Options { get; }

    /// <summary>
    /// Runs the intensity cascade over a sequence without touching fragment chemistry. Exposed so
    /// the recursion can be checked directly.
    /// </summary>
    public BondCleavageCascade ComputeCascade(string baseSequence)
    {
        if (string.IsNullOrEmpty(baseSequence))
            throw new ArgumentException("A base sequence is required.", nameof(baseSequence));

        int bondCount = Math.Max(baseSequence.Length - 1, 0);
        var probabilities = new double[bondCount];
        var remaining = new double[bondCount];
        var cleaved = new double[bondCount];

        double r = Options.StartingIonCount;
        for (int i = 0; i < bondCount; i++)
        {
            double p = BondCleavageModel.GetCleavageProbability(baseSequence, i);
            if (double.IsNaN(p) || p < 0 || p > 1)
                throw new InvalidOperationException(
                    $"Bond cleavage model returned {p} for bond {i} of '{baseSequence}'; a probability in [0, 1] is required.");

            probabilities[i] = p;
            remaining[i] = r;
            cleaved[i] = p * r;
            r -= cleaved[i];
        }

        return new BondCleavageCascade(baseSequence, Options.StartingIonCount, probabilities, remaining, cleaved, r);
    }

    public Ms2FragmentationResult Fragment(IBioPolymerWithSetMods precursor)
    {
        if (precursor is null)
            throw new ArgumentNullException(nameof(precursor));

        string baseSequence = precursor.BaseSequence;
        var cascade = ComputeCascade(baseSequence);
        var (nType, cType) = Options.ResolveProductTypes();

        // mzLib owns all fragment mass chemistry.
        var products = new List<Product>();
        precursor.Fragment(Options.DissociationType, FragmentationTerminus.Both, products);

        var nTerminalByFragmentNumber = new Dictionary<int, Product>();
        var cTerminalByFragmentNumber = new Dictionary<int, Product>();
        foreach (var product in products)
        {
            if (product.IsInternalFragment)
                continue;

            if (product.Terminus == FragmentationTerminus.N && product.ProductType == nType)
                nTerminalByFragmentNumber[product.FragmentNumber] = product;
            else if (product.Terminus == FragmentationTerminus.C && product.ProductType == cType)
                cTerminalByFragmentNumber[product.FragmentNumber] = product;
        }

        double nFraction = Options.EffectiveNTerminalIonFraction;
        var fragments = new List<SimulatedFragmentIon>();
        double unassigned = 0;

        for (int bond = 0; bond < cascade.BondCount; bond++)
        {
            int oneBasedBondNumber = bond + 1;
            double cleaved = cascade.CleavedIntensity[bond];
            if (cleaved <= 0)
                continue;

            double nShare = cleaved * nFraction;
            double cShare = cleaved - nShare;

            if (nShare > 0)
            {
                if (nTerminalByFragmentNumber.TryGetValue(oneBasedBondNumber, out var nProduct))
                    AddPeaks(fragments, nProduct, oneBasedBondNumber, nShare);
                else
                    unassigned += nShare;
            }

            if (cShare > 0)
            {
                int cFragmentNumber = baseSequence.Length - oneBasedBondNumber;
                if (cTerminalByFragmentNumber.TryGetValue(cFragmentNumber, out var cProduct))
                    AddPeaks(fragments, cProduct, oneBasedBondNumber, cShare);
                else
                    unassigned += cShare;
            }
        }

        // Intensity that had no product ion to land on (e.g. a terminus mzLib does not emit for the
        // requested series) is returned to the surviving precursor so the books still balance.
        return new Ms2FragmentationResult(
            baseSequence,
            cascade.StartingIonCount,
            fragments,
            cascade.SurvivingPrecursorIntensity + unassigned);
    }

    private void AddPeaks(List<SimulatedFragmentIon> sink, Product product, int oneBasedBondNumber, double intensity)
    {
        var charges = Options.FragmentChargeStates;
        double perCharge = intensity / charges.Count;

        IsotopeEnvelopeKernel? kernel = Options.UseIsotopeEnvelopes
            ? new IsotopeEnvelopeKernel(product.NeutralMass)
            : null;

        foreach (int charge in charges)
        {
            if (kernel is null)
            {
                sink.Add(new SimulatedFragmentIon(
                    product, oneBasedBondNumber, charge, 0, product.NeutralMass.ToMz(charge), perCharge));
                continue;
            }

            double[] centroids = kernel.CentroidMzs(charge);
            for (int i = 0; i < centroids.Length; i++)
            {
                double weight = kernel.Intensity(i);
                if (weight <= 0)
                    continue;

                sink.Add(new SimulatedFragmentIon(
                    product, oneBasedBondNumber, charge, i, centroids[i], perCharge * weight));
            }
        }
    }
}
