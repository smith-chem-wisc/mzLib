using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using System;
using System.Collections.Generic;

namespace Test.Omics
{
    // Minimal test-only implementation of ISpectralMatchHypothesis used by the ISpectralMatch test.
    internal class SimpleSpectralMatchHypothesis : ISpectralMatchHypothesis
    {
        public double Score { get; init; }
        public IBioPolymerWithSetMods SpecificBioPolymer { get; init; }
        public bool IsDecoy { get; init; }
        public string FullSequence { get; init; }
        public double? QValueNotch { get; init; }

        public SimpleSpectralMatchHypothesis(string fullSequence, IBioPolymerWithSetMods bio, double score = 0, bool isDecoy = false, double? q = null)
        {
            FullSequence = fullSequence;
            SpecificBioPolymer = bio;
            Score = score;
            IsDecoy = isDecoy;
            QValueNotch = q;
        }

        public bool Equals(ISpectralMatchHypothesis? other)
        {
            if (other is null) return false;
            return string.Equals(FullSequence, other.FullSequence, StringComparison.Ordinal)
                   && (SpecificBioPolymer?.Equals(other.SpecificBioPolymer) ?? other.SpecificBioPolymer is null)
                   && IsDecoy == other.IsDecoy;
        }

        public override bool Equals(object? obj) => Equals(obj as ISpectralMatchHypothesis);

        public override int GetHashCode()
            => HashCode.Combine(StringComparer.Ordinal.GetHashCode(FullSequence ?? string.Empty),
                                SpecificBioPolymer?.GetHashCode() ?? 0,
                                IsDecoy);
    }

    // Minimal test-only implementation of IBioPolymerWithSetMods required by the hypotheses above.
    internal class SimpleBioPolymerWithSetMods : IBioPolymerWithSetMods
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

        public SimpleBioPolymerWithSetMods(string baseSeq, string fullSeq, double mass = 0)
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

    // Concrete test implementation of ISpectralMatch
    internal class TestSpectralMatch : ISpectralMatch
    {
        public string FullFilePath { get; init; }
        public string FullSequence { get; init; }
        public string BaseSequence { get; init; }
        public List<ISpectralMatchHypothesis> BestMatchingBioPolymersWithSetMods { get; init; } = new();
        public double Score { get; init; }

        public TestSpectralMatch(string filePath, string fullSequence, string baseSequence, double score = 0)
        {
            FullFilePath = filePath ?? string.Empty;
            FullSequence = fullSequence ?? string.Empty;
            BaseSequence = baseSequence ?? string.Empty;
            Score = score;
        }

        // IComparable<ISpectralMatch>.CompareTo implementation
        public int CompareTo(ISpectralMatch? other)
        {
            if (other is null) return 1;
            // Order by Score descending (higher is better), then by FullSequence, then file path to make ordering stable.
            int scoreCmp = Score.CompareTo(other.Score);
            if (scoreCmp != 0) return scoreCmp > 0 ? 1 : -1;
            int seqCmp = string.Compare(FullSequence, other.FullSequence, StringComparison.Ordinal);
            if (seqCmp != 0) return seqCmp;
            return string.Compare(FullFilePath, other.FullFilePath, StringComparison.Ordinal);
        }

        public override bool Equals(object? obj)
        {
            var o = obj as ISpectralMatch;
            if (o == null) return false;
            return string.Equals(FullFilePath, o.FullFilePath, StringComparison.Ordinal)
                && string.Equals(FullSequence, o.FullSequence, StringComparison.Ordinal)
                && Score.Equals(o.Score);
        }

        public override int GetHashCode()
            => HashCode.Combine(StringComparer.Ordinal.GetHashCode(FullFilePath ?? string.Empty),
                                StringComparer.Ordinal.GetHashCode(FullSequence ?? string.Empty),
                                Score);

    }

    [TestFixture]
    internal class ISpectralMatchTests
    {
        [Test]
        public void CompareByScore_HigherScoreIsGreater()
        {
            var a = new TestSpectralMatch("file1", "PEPTIDE", "PEPTIDE", score: 100);
            var b = new TestSpectralMatch("file1", "PEPTIDE", "PEPTIDE", score: 50);

            Assert.That(a.CompareTo(b), Is.GreaterThan(0));
            Assert.That(b.CompareTo(a), Is.LessThan(0));
        }

        [Test]
        public void TieBreakBySequence_UsesFullSequence()
        {
            var a = new TestSpectralMatch("file1", "AAA", "AAA", score: 100);
            var b = new TestSpectralMatch("file1", "BBB", "BBB", score: 100);

            Assert.That(string.Compare(a.FullSequence, b.FullSequence, StringComparison.Ordinal), Is.LessThan(0));
            Assert.That(a.CompareTo(b), Is.LessThan(0));
        }

        [Test]
        public void HashAndEquals_EqualProperties_ProduceEqualHash()
        {
            var a = new TestSpectralMatch("path", "PEPTIDE", "PEPTIDE", score: 42);
            var b = new TestSpectralMatch("path", "PEPTIDE", "PEPTIDE", score: 42);

            Assert.That(a.Equals(b), Is.True);
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void BestMatchingHypotheses_CanBePopulatedAndRead()
        {
            var polymer = new SimpleBioPolymerWithSetMods("PEP", "PEP");
            var hyp = new SimpleSpectralMatchHypothesis("PEP", polymer, score: 10);
            var match = new TestSpectralMatch("f", "PEP", "PEP", score: 10)
            {
                BestMatchingBioPolymersWithSetMods = new List<ISpectralMatchHypothesis> { hyp }
            };

            Assert.That(match.BestMatchingBioPolymersWithSetMods.Count, Is.EqualTo(1));
            Assert.That(match.BestMatchingBioPolymersWithSetMods[0].FullSequence, Is.EqualTo("PEP"));
        }
    }
}