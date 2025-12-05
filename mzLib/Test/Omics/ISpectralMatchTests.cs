using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.Omics
{
    // Minimal test-only implementation of IBioPolymerWithSetMods required by the tests.
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

    /// <summary>
    /// Concrete test implementation of <see cref="ISpectralMatch"/>.
    /// Stores an internal, defensive copy of identified biopolymers and exposes them
    /// via <see cref="GetIdentifiedBioPolymersWithSetMods"/> as a read-only collection.
    /// </summary>
    internal class TestSpectralMatch : ISpectralMatch
    {
        private readonly List<IBioPolymerWithSetMods> _identified;

        public string FullFilePath { get; init; }
        public string FullSequence { get; init; }
        public string BaseSequence { get; init; }
        public double Score { get; init; }

        /// <summary>
        /// Construct a test spectral match.
        /// The identified collection is defensively copied; passing null creates an empty set.
        /// </summary>
        public TestSpectralMatch(string filePath, string fullSequence, string baseSequence, double score = 0, IEnumerable<IBioPolymerWithSetMods>? identified = null)
        {
            FullFilePath = filePath ?? string.Empty;
            FullSequence = fullSequence ?? string.Empty;
            BaseSequence = baseSequence ?? string.Empty;
            Score = score;
            // defensive copy to prevent external mutation
            _identified = identified?.ToList() ?? new List<IBioPolymerWithSetMods>();
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

        /// <summary>
        /// Return the identified biopolymer objects for this match.
        /// Returns a read-only snapshot; callers cannot modify the internal collection.
        /// May return zero items; null entries in the original input are preserved.
        /// </summary>
        public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods()
            => _identified.AsReadOnly();

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
        /// <summary>
        /// Higher score should sort as greater than lower score.
        /// </summary>
        [Test]
        public void CompareByScore_HigherScoreIsGreater()
        {
            var a = new TestSpectralMatch("file1", "PEPTIDE", "PEPTIDE", score: 100);
            var b = new TestSpectralMatch("file1", "PEPTIDE", "PEPTIDE", score: 50);

            Assert.That(a.CompareTo(b), Is.GreaterThan(0));
            Assert.That(b.CompareTo(a), Is.LessThan(0));
        }

        /// <summary>
        /// When scores tie, FullSequence is used as a deterministic tiebreaker.
        /// </summary>
        [Test]
        public void TieBreakBySequence_UsesFullSequence()
        {
            var a = new TestSpectralMatch("file1", "AAA", "AAA", score: 100);
            var b = new TestSpectralMatch("file1", "BBB", "BBB", score: 100);

            Assert.That(string.Compare(a.FullSequence, b.FullSequence, StringComparison.Ordinal), Is.LessThan(0));
            Assert.That(a.CompareTo(b), Is.LessThan(0));
        }

        /// <summary>
        /// Equals and GetHashCode behave consistently for two objects that share the same visible identity.
        /// </summary>
        [Test]
        public void HashAndEquals_EqualProperties_ProduceEqualHash()
        {
            var a = new TestSpectralMatch("path", "PEPTIDE", "PEPTIDE", score: 42);
            var b = new TestSpectralMatch("path", "PEPTIDE", "PEPTIDE", score: 42);

            Assert.That(a.Equals(b), Is.True);
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
        }

        /// <summary>
        /// Verify that GetIdentifiedBioPolymersWithSetMods returns the biopolymers provided at construction.
        /// Ensures order is preserved.
        /// </summary>
        [Test]
        public void GetIdentifiedBioPolymers_ReturnsProvidedBiopolymers()
        {
            var polymer1 = new SimpleBioPolymerWithSetMods("PEP1", "PEP1");
            var polymer2 = new SimpleBioPolymerWithSetMods("PEP2", "PEP2");
            var match = new TestSpectralMatch("f", "PEP1", "PEP1", score: 10, identified: new[] { polymer1, polymer2 });

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(2));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("PEP1"));
            Assert.That(identified[1].BaseSequence, Is.EqualTo("PEP2"));
        }

        /// <summary>
        /// When no identified biopolymers were provided, method returns an empty enumeration (not null).
        /// </summary>
        [Test]
        public void GetIdentifiedBioPolymers_EmptyWhenNone()
        {
            var match = new TestSpectralMatch("f", "X", "X", score: 0, identified: null);
            var identified = match.GetIdentifiedBioPolymersWithSetMods();
            Assert.That(identified, Is.Not.Null);
            Assert.That(identified, Is.Empty);
        }

        /// <summary>
        /// The implementation defensively copies the input list: mutating the original source after
        /// construction does not affect the match's returned collection.
        /// </summary>
        [Test]
        public void GetIdentifiedBioPolymers_DefensiveCopy_OriginalMutated()
        {
            var source = new List<IBioPolymerWithSetMods>
            {
                new SimpleBioPolymerWithSetMods("A","A")
            };
            var match = new TestSpectralMatch("f", "A", "A", score: 1, identified: source);

            // mutate source after construction
            source.Add(new SimpleBioPolymerWithSetMods("B", "B"));

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            // match should have only the original snapshot (one element)
            Assert.That(identified.Count, Is.EqualTo(1));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("A"));
        }

        /// <summary>
        /// Returned collection is read-only: attempts to modify it through the IList interface throw.
        /// </summary>
        [Test]
        public void GetIdentifiedBioPolymers_ReturnsReadOnlyCollection()
        {
            var polymer1 = new SimpleBioPolymerWithSetMods("P1", "P1");
            var match = new TestSpectralMatch("f", "P1", "P1", score: 1, identified: new[] { polymer1 });

            var coll = match.GetIdentifiedBioPolymersWithSetMods() as IList<IBioPolymerWithSetMods>;
            Assert.That(coll, Is.Not.Null, "Implementation should expose an IList wrapper (read-only).");
            Assert.Throws<NotSupportedException>(() => coll.Add(new SimpleBioPolymerWithSetMods("X", "X")));
        }

        /// <summary>
        /// Null entries provided in the identified list are preserved in the returned enumeration.
        /// This test documents the current behavior; callers should avoid nulls in production.
        /// </summary>
        [Test]
        public void GetIdentifiedBioPolymers_PreservesNullEntries()
        {
            var source = new List<IBioPolymerWithSetMods?>
            {
                null,
                new SimpleBioPolymerWithSetMods("Z","Z")
            };
            // cast to non-nullable interface type for constructor (we intentionally include a null)
            var match = new TestSpectralMatch("f", "Z", "Z", score: 1, identified: source.Cast<IBioPolymerWithSetMods?>().Select(x => x as IBioPolymerWithSetMods));

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(2));
            Assert.That(identified[0], Is.Null);
            Assert.That(identified[1]?.BaseSequence, Is.EqualTo("Z"));
        }
    }
}