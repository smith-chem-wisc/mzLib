#nullable enable
using System;
using System.Collections.Generic;
using MassSpectrometry;
using Omics.Fragmentation;

namespace TopDownSimulator.Ms2;

/// <summary>
/// How a cleavage event is turned into observable ions.
/// </summary>
/// <remarks>
/// The cascade assigns intensity per <em>cleavage event</em>, not per ion. A single backbone
/// cleavage physically yields two complementary pieces — an N-terminal fragment and a C-terminal
/// fragment — but only the one that retains the charge is detected. This enum makes the choice
/// explicit rather than accidental.
/// </remarks>
public enum FragmentIonAssignment
{
    /// <summary>All of the cleaved intensity goes to the N-terminal series (b or c ions).</summary>
    NTerminalOnly,

    /// <summary>All of the cleaved intensity goes to the C-terminal series (y or z ions).</summary>
    CTerminalOnly,

    /// <summary>
    /// The cleaved intensity is split between the N-terminal ion and its complementary C-terminal
    /// ion according to <see cref="PropensityCascadeOptions.NTerminalIonFraction"/>. Conserves total
    /// intensity, since the two shares sum to the cleavage event's intensity.
    /// </summary>
    Both
}

/// <summary>
/// Knobs for <see cref="PropensityCascadeFragmentationModel"/>.
/// </summary>
public sealed record PropensityCascadeOptions
{
    /// <summary>
    /// R_0: ion count present in the isolated precursor before any bond has cleaved. Fragment
    /// intensities are reported in these units, so this scales the whole spectrum.
    /// </summary>
    public double StartingIonCount { get; init; } = 1e6;

    /// <summary>
    /// Activation method, used only to pick the default ion series (see
    /// <see cref="NTerminalProductType"/>) and to stamp the written mzML. Defaults to
    /// <see cref="DissociationType.ETD"/>: electron-transfer activation is the workhorse of
    /// top-down proteomics because it preserves labile modifications and cleaves the N-Cα bond
    /// broadly along a large backbone, and the propensity table this model is built on was measured
    /// on top-down ETD data.
    /// </summary>
    public DissociationType DissociationType { get; init; } = DissociationType.ETD;

    /// <summary>
    /// N-terminal ion series. Null (the default) means "derive from
    /// <see cref="DissociationType"/>": c for ETD/ECD/EThcD, b for CID/HCD/IRMPD.
    /// </summary>
    public ProductType? NTerminalProductType { get; init; }

    /// <summary>
    /// C-terminal ion series. Null (the default) means "derive from
    /// <see cref="DissociationType"/>": z-dot for ETD/ECD/EThcD, y for CID/HCD/IRMPD.
    /// </summary>
    public ProductType? CTerminalProductType { get; init; }

    /// <summary>Which side of a cleavage is detected. Defaults to <see cref="FragmentIonAssignment.Both"/>.</summary>
    public FragmentIonAssignment IonAssignment { get; init; } = FragmentIonAssignment.Both;

    /// <summary>
    /// Share of a cleavage event's intensity given to the N-terminal ion when
    /// <see cref="IonAssignment"/> is <see cref="FragmentIonAssignment.Both"/>. The C-terminal ion
    /// gets the remainder. Defaults to 0.5 — an even split, which is the neutral assumption in the
    /// absence of a charge-retention model.
    /// </summary>
    public double NTerminalIonFraction { get; init; } = 0.5;

    /// <summary>
    /// Charge states at which each fragment is observed. A fragment's intensity is split evenly
    /// across them. Defaults to singly charged only, which keeps the basic model's peak list one
    /// peak per fragment ion.
    /// </summary>
    public IReadOnlyList<int> FragmentChargeStates { get; init; } = new[] { 1 };

    /// <summary>
    /// When true, each fragment's intensity is spread over an averagine isotope envelope via
    /// <see cref="Model.IsotopeEnvelopeKernel"/> instead of landing entirely on the monoisotopic
    /// peak. Defaults to false: a single peak per fragment is the honest "basic" model and keeps the
    /// intensity bookkeeping trivially checkable, at the cost of spectra that look unrealistically
    /// clean for fragments above a few kDa, where the monoisotopic peak is not the tallest one.
    /// </summary>
    public bool UseIsotopeEnvelopes { get; init; }

    /// <summary>
    /// Default N-terminal / C-terminal series for an activation method. Kept separate from mzLib's
    /// <c>DissociationTypeCollection.ProductsFromDissociationType</c> because that returns every
    /// series mzLib will search for (ETD lists c, y <em>and</em> z-dot), whereas the cascade needs
    /// exactly one series per terminus.
    /// </summary>
    public static (ProductType NTerminal, ProductType CTerminal) DefaultProductTypes(DissociationType dissociationType) =>
        dissociationType switch
        {
            DissociationType.ETD or DissociationType.ECD or DissociationType.EThcD
                => (ProductType.c, ProductType.zDot),
            DissociationType.CID or DissociationType.LowCID or DissociationType.HCD
                or DissociationType.IRMPD or DissociationType.AnyActivationType
                => (ProductType.b, ProductType.y),
            _ => throw new ArgumentOutOfRangeException(nameof(dissociationType), dissociationType,
                "No default ion series for this dissociation type. Set NTerminalProductType and " +
                "CTerminalProductType explicitly.")
        };

    /// <summary>Resolves the ion series actually used, honouring any explicit overrides.</summary>
    public (ProductType NTerminal, ProductType CTerminal) ResolveProductTypes()
    {
        if (NTerminalProductType.HasValue && CTerminalProductType.HasValue)
            return (NTerminalProductType.Value, CTerminalProductType.Value);

        var defaults = DefaultProductTypes(DissociationType);
        return (NTerminalProductType ?? defaults.NTerminal, CTerminalProductType ?? defaults.CTerminal);
    }

    /// <summary>Share of a cleavage event given to the N-terminal ion under the current assignment.</summary>
    public double EffectiveNTerminalIonFraction => IonAssignment switch
    {
        FragmentIonAssignment.NTerminalOnly => 1.0,
        FragmentIonAssignment.CTerminalOnly => 0.0,
        FragmentIonAssignment.Both => NTerminalIonFraction,
        _ => throw new ArgumentOutOfRangeException(nameof(IonAssignment), IonAssignment, null)
    };

    internal void Validate()
    {
        if (StartingIonCount <= 0)
            throw new ArgumentOutOfRangeException(nameof(StartingIonCount), "Starting ion count must be positive.");
        if (NTerminalIonFraction < 0 || NTerminalIonFraction > 1)
            throw new ArgumentOutOfRangeException(nameof(NTerminalIonFraction), "Must be in [0, 1].");
        if (FragmentChargeStates is null || FragmentChargeStates.Count == 0)
            throw new ArgumentException("At least one fragment charge state is required.", nameof(FragmentChargeStates));

        foreach (int z in FragmentChargeStates)
        {
            if (z < 1)
                throw new ArgumentOutOfRangeException(nameof(FragmentChargeStates), z, "Charge states must be >= 1.");
        }

        // Throws with a clear message if the dissociation type has no default series and none was supplied.
        ResolveProductTypes();
    }
}
