using NUnit.Framework;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Omics;
using static Test.Omics.BioPolymerGroupSequenceCoverageTests;

namespace Test.Omics
{
    /// <summary>
    /// Minimal tests for ISpectralMatch interface behavior.
    /// Note: Most tests use TestSpectralMatch helper class defined in test files that need it.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class ISpectralMatchTests
    {
        /// <summary>
        /// Verifies ISpectralMatch.CompareTo orders by Score descending (higher scores first).
        /// Critical: Determines PSM ranking for FDR calculation and best-match selection.
        /// </summary>
        [Test]
        public void CompareTo_HigherScoreComesFirst()
        {
            ISpectralMatch a = new CoverageSpectralMatch("file1", "PEPTIDE", "PEPTIDE", score: 100, scanNumber: 1);
            ISpectralMatch b = new CoverageSpectralMatch("file1", "PEPTIDE", "PEPTIDE", score: 50, scanNumber: 1);

            // Higher score should compare as greater (comes first when sorted descending)
            Assert.That(highScore.CompareTo(lowScore), Is.GreaterThan(0));
            Assert.That(lowScore.CompareTo(highScore), Is.LessThan(0));
        }

        /// <summary>
        /// Verifies GetIdentifiedBioPolymersWithSetMods returns empty enumeration (not null) when none provided.
        /// Critical: Prevents null reference exceptions in protein grouping code.
        /// </summary>
        [Test]
        public void GetIdentifiedBioPolymersWithSetMods_ReturnsEmptyNotNull()
        {
            var match = new TestSpectralMatch("file", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);

            var identified = match.GetIdentifiedBioPolymersWithSetMods();

            Assert.That(identified, Is.Not.Null);
            Assert.That(identified, Is.Empty);
        }

        #region Test Helper

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
        public void GetIdentifiedBioPolymersWithSetMods_ReturnsReadOnlyCollection()
        {
            var polymer1 = new SimpleBioPolymerWithSetMods("P1", "P1");
            var match = new TestSpectralMatch("f", "P1", "P1", score: 1, scanNumber: 1, identified: new[] { polymer1 });

            var coll = match.GetIdentifiedBioPolymersWithSetMods() as IList<IBioPolymerWithSetMods>;
            Assert.That(coll, Is.Not.Null, "Implementation should expose an IList wrapper (read-only).");
            Assert.Throws<NotSupportedException>(() => coll.Add(new SimpleBioPolymerWithSetMods("X", "X")));
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
            var match = new TestSpectralMatch("f", "Z", "Z", score: 1, scanNumber: 1, identified: identifiedList);

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
        /// GetSequenceCoverage returns null when no fragment positions are set.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_NoFragments_ReturnsNull()
        {
            var match = new TestSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            match.GetSequenceCoverage();
            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Null);
        }

        /// <summary>
        /// GetSequenceCoverage returns null for empty base sequence.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_EmptySequence_ReturnsNull()
        {
            var match = new TestSpectralMatch("f", "", "", score: 1, scanNumber: 1);
            match.NTerminalFragmentPositions = new List<int> { 1, 2 };
            match.GetSequenceCoverage();
            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Null);
        }

        /// <summary>
        /// Sequential N-terminal fragments cover the second position.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_SequentialNTermFragments_CoversSecondPosition()
        {
            var match = new TestSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            match.NTerminalFragmentPositions = new List<int> { 1, 2 };
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
            var match = new TestSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            match.CTerminalFragmentPositions = new List<int> { 2, 3 };
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
            var match = new TestSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            match.CTerminalFragmentPositions = new List<int> { 2 };
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
            var match = new TestSpectralMatch("f", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);
            // BaseSequence.Length is 7, so last N-term position is 6 (Length - 1)
            match.NTerminalFragmentPositions = new List<int> { 6 };
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