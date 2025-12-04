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
        public string BaseSequence { get; init; }
        public string FullSequence { get; init; }

        // Name difference in production code: IBioPolymerWithSetMods exposes MostAbundantMonoisotopicMass,
        // while IHasMass (via IHasChemicalFormula) requires MonoisotopicMass. Provide both.
        public double MostAbundantMonoisotopicMass { get; init; }
        public double MonoisotopicMass => MostAbundantMonoisotopicMass;

        // Minimal ChemicalFormula implementation for tests
        public ChemicalFormula ThisChemicalFormula => new ChemicalFormula();

        public string SequenceWithChemicalFormulas => FullSequence;
        public int OneBasedStartResidue { get; init; } = 1;
        public int OneBasedEndResidue { get; init; } = 1;
        public int MissedCleavages { get; init; } = 0;
        public string Description { get; init; } = string.Empty;
        public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; } = CleavageSpecificity.Full;
        public char PreviousResidue { get; init; } = '-';
        public char NextResidue { get; init; } = '-';
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

    // Concrete, test-only implementation of ISpectralMatchHypothesis.
    internal class TestSpectralMatchHypothesis : ISpectralMatchHypothesis
    {
        public double Score { get; init; }
        public IBioPolymerWithSetMods SpecificBioPolymer { get; init; }
        public bool IsDecoy { get; init; }
        public string FullSequence { get; init; }
        public double? QValueNotch { get; init; }

        public TestSpectralMatchHypothesis(string fullSequence, IBioPolymerWithSetMods specific, double score = 0, bool isDecoy = false, double? q = null)
        {
            FullSequence = fullSequence;
            SpecificBioPolymer = specific;
            Score = score;
            IsDecoy = isDecoy;
            QValueNotch = q;
        }

        // Value-based equality: FullSequence + SpecificBioPolymer.FullSequence + IsDecoy.
        public bool Equals(ISpectralMatchHypothesis? other)
        {
            if (other == null) return false;
            if (ReferenceEquals(this, other)) return true;

            bool bioEqual = (SpecificBioPolymer == null && other.SpecificBioPolymer == null)
                || (SpecificBioPolymer != null && SpecificBioPolymer.Equals(other.SpecificBioPolymer));

            return string.Equals(FullSequence, other.FullSequence, StringComparison.Ordinal)
                && bioEqual
                && IsDecoy == other.IsDecoy;
        }

        public override bool Equals(object? obj) => Equals(obj as ISpectralMatchHypothesis);

        public override int GetHashCode()
        {
            int bioHash = SpecificBioPolymer?.GetHashCode() ?? 0;
            return HashCode.Combine(
                StringComparer.Ordinal.GetHashCode(FullSequence ?? string.Empty),
                bioHash,
                IsDecoy);
        }
    }

    [TestFixture]
    public class ISpectralMatchHypothesisTests
    {
        [Test]
        public void Equality_SameValues_AreEqual()
        {
            var bio1 = new TestBioPolymerWithSetMods("PEPTIDE", "PEPTIDE");
            var a = new TestSpectralMatchHypothesis("PEPTIDE", bio1, score: 100, isDecoy: false);
            var b = new TestSpectralMatchHypothesis("PEPTIDE", new TestBioPolymerWithSetMods("PEPTIDE", "PEPTIDE"), score: 90, isDecoy: false);

            Assert.That(a.Equals(b), Is.True, "Hypotheses with same sequence and bio polymer should be equal");
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()), "Equal hypotheses must produce equal hash codes");
        }

        [Test]
        public void Inequality_DifferentDecoyFlag_NotEqual()
        {
            var bio = new TestBioPolymerWithSetMods("PEPTIDE", "PEPTIDE");
            var target = new TestSpectralMatchHypothesis("PEPTIDE", bio, isDecoy: false);
            var decoy = new TestSpectralMatchHypothesis("PEPTIDE", bio, isDecoy: true);

            Assert.That(target.Equals(decoy), Is.False);
        }

        [Test]
        public void Inequality_DifferentSequence_NotEqual()
        {
            var bio = new TestBioPolymerWithSetMods("PEPTIDE", "PEPTIDE");
            var a = new TestSpectralMatchHypothesis("PEPTIDE", bio);
            var b = new TestSpectralMatchHypothesis("PEPTIDEX", bio);

            Assert.That(a.Equals(b), Is.False);
        }
    }
}
