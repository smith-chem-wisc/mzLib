using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectralMatch;
using System;
using System.Collections.Generic;
using System.Linq;
using static Test.Omics.BioPolymerGroupSequenceCoverageTests;

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

    [TestFixture]
    internal class ISpectralMatchTests
    {
        /// <summary>
        /// Higher score should sort as greater than lower score.
        /// CompareTo returns a positive value when this.Score &gt; other.Score (higher score is preferred).
        /// </summary>
        [Test]
        public void CompareByScore_HigherScoreIsGreater()
        {
            ISpectralMatch a = new CoverageSpectralMatch("file1", "PEPTIDE", "PEPTIDE", score: 100, scanNumber: 1);
            ISpectralMatch b = new CoverageSpectralMatch("file1", "PEPTIDE", "PEPTIDE", score: 50, scanNumber: 1);

            Assert.That(a.CompareTo(b), Is.GreaterThan(0));
            Assert.That(b.CompareTo(a), Is.LessThan(0));
        }

        /// <summary>
        /// When all above tie, ScanNumber is used as final tie-breaker (ascending).
        /// </summary>
        [Test]
        public void TieBreakByScanNumber_UsesScanNumber()
        {
            ISpectralMatch a = new CoverageSpectralMatch("file", "SEQ", "BASE", score: 100, scanNumber: 5);
            ISpectralMatch b = new CoverageSpectralMatch("file", "SEQ", "BASE", score: 100, scanNumber: 10);

            Assert.That(a.CompareTo(b), Is.GreaterThan(0));
            Assert.That(b.CompareTo(a), Is.LessThan(0));
        }

        /// <summary>
        /// CompareTo with null should return positive (this &gt; null).
        /// </summary>
        [Test]
        public void CompareTo_Null_ReturnsNegative()
        {
            ISpectralMatch a = new CoverageSpectralMatch("file", "SEQ", "BASE", score: 1, scanNumber: 1);
            Assert.That(a.CompareTo(null), Is.LessThan(0));
        }

        /// <summary>
        /// Two objects with identical ordering-relevant fields compare equal (CompareTo == 0).
        /// </summary>
        [Test]
        public void CompareTo_EqualObjects_ReturnsZero()
        {
            ISpectralMatch a = new CoverageSpectralMatch("file", "SEQ", "BASE", score: 10, scanNumber: 2);
            ISpectralMatch b = new CoverageSpectralMatch("file", "SEQ", "BASE", score: 10, scanNumber: 2);
            Assert.That(a.CompareTo(b), Is.EqualTo(0));
        }

        /// <summary>
        /// Equals and GetHashCode behave consistently for two objects that share the same visible identity.
        /// ScanNumber is part of equality semantics in the test implementation.
        /// </summary>
        [Test]
        public void HashAndEquals_EqualProperties_ProduceEqualHash()
        {
            ISpectralMatch a = new CoverageSpectralMatch("path", "PEPTIDE", "PEPTIDE", score: 42, scanNumber: 7);
            ISpectralMatch b = new CoverageSpectralMatch("path", "PEPTIDE", "PEPTIDE", score: 42, scanNumber: 7);

            Assert.That(a.Equals(b), Is.True);
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
        }

        /// <summary>
        /// Additional equality check: differing scan numbers make objects not equal.
        /// </summary>
        [Test]
        public void DifferentScanNumbers_AreNotEqual()
        {
            ISpectralMatch a = new CoverageSpectralMatch("path", "PEP", "PEP", score: 5, scanNumber: 1);
            ISpectralMatch b = new CoverageSpectralMatch("path", "PEP", "PEP", score: 5, scanNumber: 2);

            Assert.That(a.Equals(b), Is.False);
            Assert.That(a.GetHashCode(), Is.Not.EqualTo(b.GetHashCode()));
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
            ISpectralMatch match = new CoverageSpectralMatch("f", "PEP1", "PEP1", score: 10, scanNumber: 1, identified: new[] { polymer1, polymer2 });

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(2));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("PEP1"));
            Assert.That(identified[1].BaseSequence, Is.EqualTo("PEP2"));
        }

        /// <summary>
        /// When no identified biopolymers were provided, method returns an empty enumeration (not null).
        /// </summary>
        [Test]
        public void GetIdentifiedBioPolymersWithSetMods_EmptyWhenNone()
        {
            ISpectralMatch match = new CoverageSpectralMatch("f", "X", "X", score: 0, scanNumber: 1, identified: null);
            var identified = match.GetIdentifiedBioPolymersWithSetMods();
            Assert.That(identified, Is.Not.Null);
            Assert.That(identified, Is.Empty);
        }

        /// <summary>
        /// The implementation defensively copies the input list: mutating the original source after
        /// construction does not affect the match's returned collection.
        /// </summary>
        [Test]
        public void GetIdentifiedBioPolymersWithSetMods_DefensiveCopy_OriginalMutated()
        {
            var source = new List<IBioPolymerWithSetMods>
            {
                new SimpleBioPolymerWithSetMods("A","A")
            };
            ISpectralMatch match = new CoverageSpectralMatch("f", "A", "A", score: 1, scanNumber: 1, identified: source);

            // mutate source after construction
            source.Add(new SimpleBioPolymerWithSetMods("B", "B"));

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            // match should have only the original snapshot (one element)
            Assert.That(identified.Count, Is.EqualTo(1));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("A"));
        }

        /// <summary>
        /// Null entries provided in the identified list are preserved in the returned enumeration.
        /// This test verifies that behavior by asserting the first returned element is null
        /// and that a subsequent element matches the expected polymer.
        /// </summary>
        [Test]
        public void GetIdentifiedBioPolymersWithSetMods_PreservesNullEntries()
        {
            // Arrange: create a source list containing a null entry and a real biopolymer
            var polymer = new SimpleBioPolymerWithSetMods("Z", "Z");
            var identifiedList = new List<IBioPolymerWithSetMods?> { null, polymer };

            // Create match with the test list (constructor makes a defensive copy)
            var match = new CoverageSpectralMatch("f", "Z", "Z", score: 1, scanNumber: 1, identified: identifiedList);

            // Act - call the method under test and materialize the results
            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();

            // Assert - use NUnit collection constraints and Matches predicate:
            // - exactly one null element
            // - exactly one element matching the predicate (BaseSequence == "Z")
            Assert.That(identified, Has.Exactly(1).Null, "Expected exactly one null entry in the identified collection.");
            Assert.That(identified, Has.Exactly(1).Matches<IBioPolymerWithSetMods>(p => p != null && p.BaseSequence == "Z"),
                "Expected exactly one identified biopolymer with BaseSequence == \"Z\".");

            // Positional verification: first element must be null and second must be the polymer
            Assert.That(identified.Count, Is.EqualTo(2));
            Assert.That(identified[0], Is.Null, "The first returned element should be the preserved null entry.");
            Assert.That(identified[1]?.BaseSequence, Is.EqualTo("Z"));
        }

        #region GetSequenceCoverage Tests

        /// <summary>
        /// GetSequenceCoverage returns empty list when no fragment ions are set.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_NoFragments_ReturnsEmpty()
        {
            var match = new CoverageSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            match.MatchedFragmentIons = new List<MatchedFragmentIon>();
            match.GetSequenceCoverage();
            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Empty);
        }

        /// <summary>
        /// GetSequenceCoverage returns empty list for empty base sequence.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_EmptySequence_ReturnsEmpty()
        {
            var match = new CoverageSpectralMatch("f", "", "", score: 1, scanNumber: 1);
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 200.0, 10.0, 1),
            };
            match.GetSequenceCoverage();
            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Empty);
        }

        /// <summary>
        /// Sequential N-terminal fragments cover the second position.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_SequentialNTermFragments_CoversSecondPosition()
        {
            var match = new CoverageSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 200.0, 10.0, 1),
            };
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(1)); // First position covered
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(2)); // Sequential coverage
        }

        /// <summary>
        /// Sequential C-terminal fragments cover the first position.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_SequentialCTermFragments_CoversFirstPosition()
        {
            var match = new CoverageSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 2, 2, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 3, 3, 0), 200.0, 10.0, 1),
            };
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(1)); // Position 2 covers first residue
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(2)); // Sequential coverage
        }

        /// <summary>
        /// C-terminal fragment at position 2 covers the first residue.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_CTermAtPosition2_CoversFirstResidue()
        {
            var match = new CoverageSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 2, 2, 0), 100.0, 10.0, 1),
            };
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(1));
        }

        /// <summary>
        /// N-terminal fragment at the last position covers the last residue.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_NTermAtLastPosition_CoversLastResidue()
        {
            var match = new CoverageSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            // BaseSequence.Length is 7, so last N-term position is 6 (Length - 1)
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 6, 6, 0), 100.0, 10.0, 1),
            };
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(7)); // Last residue covered
        }

        /// <summary>
        /// C-terminal fragment at the last position covers the last residue.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_CTermAtLastPosition_CoversLastResidue()
        {
            var match = new CoverageSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 7, 7, 0), 100.0, 10.0, 1),
            };
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(7));
        }

        /// <summary>
        /// Coverage from both N and C terminal fragments at the same position.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_BothTerminiAtSamePosition_CoversPosition()
        {
            var match = new CoverageSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 3, 3, 0), 200.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 3, 3, 0), 300.0, 10.0, 1),
            };
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(3));
        }

        #endregion
    }
}