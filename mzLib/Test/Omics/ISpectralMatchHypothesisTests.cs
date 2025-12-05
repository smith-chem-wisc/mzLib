using System;
using System.Collections.Generic;
using NUnit.Framework;
using Omics;
using Chemistry;
using Omics.Digestion;
using Omics.Modifications;
using MassSpectrometry;
using Omics.Fragmentation; // required for IHasChemicalFormula / ChemicalFormula

namespace Test.Omics
{
    // Minimal test stub implementing IBioPolymerWithSetMods.
    // Only members required by the unit test are implemented with simple defaults.
    internal class TestBioPolymerWithSetMods : IBioPolymerWithSetMods
    {
        public string BaseSequence { get; }
        public string FullSequence { get; }

        // Name difference in production code: IBioPolymerWithSetMods exposes MostAbundantMonoisotopicMass,
        // while IHasMass (via IHasChemicalFormula) requires MonoisotopicMass. Provide both.
        public double MostAbundantMonoisotopicMass { get; }
        public double MonoisotopicMass => MostAbundantMonoisotopicMass;

        // Minimal ChemicalFormula implementation for tests
        public ChemicalFormula ThisChemicalFormula => new ChemicalFormula();

        public string SequenceWithChemicalFormulas => FullSequence;
        public int OneBasedStartResidue { get; } = 1;
        public int OneBasedEndResidue { get; } = 1;
        public int MissedCleavages { get; } = 0;
        public string Description { get; } = string.Empty;
        public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; } = CleavageSpecificity.Full;
        public char PreviousResidue { get; } = '-';
        public char NextResidue { get; } = '-';
        public IDigestionParams DigestionParams => throw new NotImplementedException();
        public Dictionary<int, Modification> AllModsOneIsNterminus => new();
        public int NumMods => 0;
        public int NumFixedMods => 0;
        public int NumVariableMods => 0;
        public int Length => BaseSequence?.Length ?? 0;

        // Indexer required by IBioPolymerWithSetMods
        public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

        public IBioPolymer Parent => throw new NotImplementedException();

        public TestBioPolymerWithSetMods(string baseSeq, string fullSeq, double mass = 0)
        {
            BaseSequence = baseSeq;
            FullSequence = fullSeq;
            MostAbundantMonoisotopicMass = mass;
        }

        public void Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus, List<Product> products, FragmentationParams? fragmentationParams = null)
            => throw new NotImplementedException();

        public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments, List<Product> products, FragmentationParams? fragmentationParams = null)
            => throw new NotImplementedException();

        public IBioPolymerWithSetMods Localize(int indexOfMass, double massToLocalize) => this;

        public bool Equals(IBioPolymerWithSetMods? other)
        {
            if (other == null) return false;
            return string.Equals(BaseSequence, other.BaseSequence, StringComparison.Ordinal)
                && string.Equals(FullSequence, other.FullSequence, StringComparison.Ordinal);
        }

        public override bool Equals(object? obj) => Equals(obj as IBioPolymerWithSetMods);

        public override int GetHashCode()
            => HashCode.Combine(StringComparer.Ordinal.GetHashCode(BaseSequence ?? string.Empty),
                                StringComparer.Ordinal.GetHashCode(FullSequence ?? string.Empty));
    }

}
