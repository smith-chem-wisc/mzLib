using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using Chemistry;
using MassSpectrometry;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;

namespace Test.Omics;

[ExcludeFromCodeCoverage]
public class MockBioPolymerWithSetMods : IBioPolymerWithSetMods
{
    public string BaseSequence { get; }
    public string FullSequence { get; }
    public double MostAbundantMonoisotopicMass { get; } = 0;
    public double MonoisotopicMass { get; } = 0;
    public string SequenceWithChemicalFormulas => BaseSequence;
    public int OneBasedStartResidue { get; }
    public int OneBasedEndResidue { get; }
    public int MissedCleavages => 0;
    public string Description => "";
    public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; } = CleavageSpecificity.Full;
    public char PreviousResidue => '-';
    public char NextResidue => '-';
    public IDigestionParams DigestionParams => null!;
    public Dictionary<int, Modification> AllModsOneIsNterminus { get; }
    public int NumMods => AllModsOneIsNterminus?.Count ?? 0;
    public int NumFixedMods => 0;
    public int NumVariableMods => NumMods;
    public int Length => BaseSequence.Length;
    public IBioPolymer Parent { get; }
    public ChemicalFormula ThisChemicalFormula => new();
    public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

    public MockBioPolymerWithSetMods(string baseSequence, string fullSequence, IBioPolymer parent = null,
        int startResidue = 1, int endResidue = 0, Dictionary<int, Modification> mods = null)
    {
        BaseSequence = baseSequence;
        FullSequence = fullSequence;
        Parent = parent;
        OneBasedStartResidue = startResidue;
        OneBasedEndResidue = endResidue > 0 ? endResidue : startResidue + baseSequence.Length - 1;
        AllModsOneIsNterminus = mods ?? new Dictionary<int, Modification>();
        MostAbundantMonoisotopicMass = MonoisotopicMass = 500.0;
    }

    public MockBioPolymerWithSetMods(IBioPolymer parent, int start, int end, Dictionary<int, Modification>? mods = null)
    {
        Parent = parent;
        OneBasedStartResidue = start;
        OneBasedEndResidue = end;
        BaseSequence = parent.BaseSequence.Substring(start - 1, end - start + 1);
        AllModsOneIsNterminus = mods ?? new();
        FullSequence = BaseSequence; // Simplified
    }

    public void Fragment(DissociationType d, FragmentationTerminus t, List<Product> p, FragmentationParams? f = null) { }
    public void FragmentInternally(DissociationType d, int m, List<Product> p, FragmentationParams? f = null) { }
    public IBioPolymerWithSetMods Localize(int i, double m) => this;
    public bool Equals(IBioPolymerWithSetMods? other) => other != null && BaseSequence == other.BaseSequence;
    public override bool Equals(object? obj) => obj is IBioPolymerWithSetMods other && Equals(other);
    public override int GetHashCode() => BaseSequence.GetHashCode();
}