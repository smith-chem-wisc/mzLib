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
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using static Test.Omics.BioPolymerGroupSequenceCoverageTests;

namespace Test.Omics
{
    /// <summary>
    /// Tests for BaseSpectralMatch class.
    /// Validates PSM storage, scoring, fragment coverage, and collection behavior.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class SpectralMatchTests
    {
        #region Constructor Tests

        /// <summary>
        /// Verifies constructor sets all properties and handles null inputs gracefully.
        /// Critical: PSM data must be correctly stored for downstream analysis.
        /// </summary>
        [Test]
        public void Constructor_SetsAllPropertiesAndHandlesNulls()
        {
            MockSpectralMatch match = new MockSpectralMatch("path/to/file.mzML", "PEP[Phospho]TIDE", "PEPTIDE", 100.5, 42);

            Assert.Multiple(() =>
            {
                Assert.That(match.FullFilePath, Is.EqualTo("path/to/file.mzML"));
                Assert.That(match.OneBasedScanNumber, Is.EqualTo(42));
                Assert.That(match.Score, Is.EqualTo(100.5));
                Assert.That(match.FullSequence, Is.EqualTo("PEP[Phospho]TIDE"));
                Assert.That(match.BaseSequence, Is.EqualTo("PEPTIDE"));
            });

            // Null handling
            var nullMatch = new MockSpectralMatch(null!, null!, null!, 10.0, 1);
            Assert.That(nullMatch.FullFilePath, Is.EqualTo(string.Empty));
            Assert.That(nullMatch.FullSequence, Is.EqualTo(string.Empty));
            Assert.That(nullMatch.BaseSequence, Is.EqualTo(string.Empty));
        }

        #endregion

        #region Identified BioPolymers Tests

        /// <summary>
        /// Verifies biopolymers can be added and retrieved, with defensive copy protection.
        /// Critical: Protein inference requires correct peptide-to-protein mapping.
        /// </summary>
        [Test]
        public void IdentifiedBioPolymers_AddAndRetrieveWithDefensiveCopy()
        {
            var polymer1 = new MockBioPolymerWithSetMods("ABC", "ABC");
            var polymer2 = new MockBioPolymerWithSetMods("DEF", "DEF");

            // Constructor with biopolymers
            var match = new MockSpectralMatch("file", "ABC", "ABC", 10.0, 1, new[] { polymer1, polymer2 });

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(2));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("ABC"));
            Assert.That(identified[1].BaseSequence, Is.EqualTo("DEF"));

            // Verify defensive copy - original list mutation doesn't affect match
            var source = new List<IBioPolymerWithSetMods> { new MockBioPolymerWithSetMods("X", "X") };
            var match2 = new MockSpectralMatch("file", "X", "X", 10.0, 1, source);
            source.Add(new MockBioPolymerWithSetMods("Y", "Y"));
            Assert.That(match2.GetIdentifiedBioPolymersWithSetMods().Count(), Is.EqualTo(1));
        }

        #endregion

        #region GetSequenceCoverage Tests

        /// <summary>
        /// Verifies fragment coverage calculation from N-terminal and C-terminal ions.
        /// Critical: Sequence coverage is used for protein-level coverage reporting.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_CalculatesFromBothTermini()
        {
            var match = new MockSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            // Sequential N-term fragments cover positions
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 200.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 3, 3, 0), 300.0, 10.0, 1),
            };
            match.GetSequenceCoverage();
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(1));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(2));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(3));

            // C-term fragments also work
            var match2 = new MockSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);
            match2.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 2, 2, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 3, 3, 0), 200.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 4, 4, 0), 300.0, 10.0, 1),
            };
            match2.GetSequenceCoverage();
            Assert.That(match2.FragmentCoveragePositionInPeptide, Contains.Item(1)); // y2 covers first

            // Both termini combined
            var match3 = new MockSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);
            match3.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 200.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 6, 6, 0), 300.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 7, 7, 0), 400.0, 10.0, 1),
            };
            match3.GetSequenceCoverage();
            Assert.That(match3.FragmentCoveragePositionInPeptide, Contains.Item(1));
            Assert.That(match3.FragmentCoveragePositionInPeptide, Contains.Item(7));
        }

        /// <summary>
        /// Verifies edge cases: empty sequence and no fragments return empty coverage.
        /// Critical: Prevents crashes with invalid input.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_EdgeCases_ReturnEmpty()
        {
            var emptySeq = new MockSpectralMatch("file", "", "", 10.0, 1);
            emptySeq.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
            };
            emptySeq.GetSequenceCoverage();
            Assert.That(emptySeq.FragmentCoveragePositionInPeptide, Is.Empty);

            var noFragments = new MockSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);
            noFragments.MatchedFragmentIons = new List<MatchedFragmentIon>();
            noFragments.GetSequenceCoverage();
            Assert.That(noFragments.FragmentCoveragePositionInPeptide, Is.Empty);
        }

        #endregion

        #region GetSequenceCoverage Advanced Tests

        /// <summary>
        /// Verifies that when the final N-terminal fragment (Length-1) is present, the last residue is covered.
        /// Critical: Ensures complete sequence coverage is reported when b(n-1) ion is matched.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_FinalNTerminalFragment_CoversLastResidue()
        {
            // For "PEPTIDE" (length 7), b6 would be position 6 (Length-1)
            var match = new MockSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            // Position 6 = BaseSequence.Length - 1 = 7 - 1 = 6
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 6, 6, 0), 100.0, 10.0, 1),
            };
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(7),
                "Last residue (position 7) should be covered when b6 is present");
        }

        /// <summary>
        /// Verifies residues are covered when both N-terminal and C-terminal fragments match at the same position.
        /// Critical: Bidirectional fragment evidence provides strongest sequence confirmation.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_BothDirectionsInclusive_CoversPosition()
        {
            var match = new MockSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            // N-term positions 2,3 and C-term position 3 - position 3 should be covered (inclusive)
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 3, 3, 0), 200.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 3, 3, 0), 300.0, 10.0, 1),
            };
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(3),
                "Position 3 should be covered when both b3 and y3 are present (inclusive coverage)");
        }

        /// <summary>
        /// Verifies residues are covered via exclusive bidirectional evidence (n-term + c-term+2).
        /// Critical: Handles the offset relationship between b and y ion numbering for internal residues.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_BothDirectionsExclusive_CoversPosition()
        {
            var match = new MockSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            // N-term positions 2,3 and C-term position 5 (3+2) should cover position 4
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 3, 3, 0), 200.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 5, 5, 0), 300.0, 10.0, 1),
            };
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(4),
                "Position 4 should be covered when b3 is present and y5 (3+2) provides exclusive coverage");
        }

        #endregion

        #region CompareTo IHasSequenceCoverageFromFragments Tests

        /// <summary>
        /// Verifies CompareTo(IHasSequenceCoverageFromFragments) delegates to CompareTo(ISpectralMatch) when applicable.
        /// Critical: Ensures consistent ordering when BaseSpectralMatch is used in sorted collections via either interface.
        /// </summary>
        [Test]
        public void CompareTo_IHasSequenceCoverageFromFragments_DelegatesToISpectralMatch()
        {
            var high = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var low = new MockSpectralMatch("file", "SEQ", "SEQ", 50.0, 1);

            // Cast to IHasSequenceCoverageFromFragments
            IHasSequenceCoverageFromFragments highCoverage = high;
            IHasSequenceCoverageFromFragments lowCoverage = low;

            // Should delegate to ISpectralMatch comparison (score ascending)
            Assert.That(high.CompareTo(lowCoverage), Is.GreaterThan(0),
                "Higher score should come first when comparing via IHasSequenceCoverageFromFragments");
            Assert.That(low.CompareTo(highCoverage), Is.LessThan(0));
        }

        /// <summary>
        /// Verifies CompareTo(IHasSequenceCoverageFromFragments) returns 0 for non-ISpectralMatch implementations.
        /// Critical: Prevents crashes when mixed types implement only IHasSequenceCoverageFromFragments.
        /// </summary>
        [Test]
        public void CompareTo_IHasSequenceCoverageFromFragments_NonSpectralMatch_ReturnsZero()
        {
            var match = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);

            // Create a mock that implements IHasSequenceCoverageFromFragments but not ISpectralMatch
            var nonSpectralMatch = new NonSpectralMatchCoverageProvider();

            Assert.That(match.CompareTo(nonSpectralMatch), Is.EqualTo(0),
                "Should return 0 for IHasSequenceCoverageFromFragments that is not ISpectralMatch");
        }

        /// <summary>
        /// Verifies CompareTo(IHasSequenceCoverageFromFragments) handles null correctly.
        /// Critical: Null safety for interface-based comparisons.
        /// </summary>
        [Test]
        public void CompareTo_IHasSequenceCoverageFromFragments_Null_ReturnsNegative()
        {
            var match = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);

            Assert.That(match.CompareTo((IHasSequenceCoverageFromFragments?)null), Is.LessThan(0));
        }

        // Helper class for testing non-ISpectralMatch comparison
        public class NonSpectralMatchCoverageProvider : IHasSequenceCoverageFromFragments
        {
            public string BaseSequence { get; } = "PEPTIDE";
            public List<MatchedFragmentIon> MatchedFragmentIons { get; set; } = new();
            public List<int> FragmentCoveragePositionInPeptide { get; set; }
            public void GetSequenceCoverage() { }
        }

        #endregion

        #region Equals(object) Tests

        /// <summary>
        /// Verifies Equals(object) correctly handles BaseSpectralMatch type.
        /// Critical: Required for correct behavior in non-generic collections like ArrayList.
        /// </summary>
        [Test]
        public void Equals_Object_BaseSpectralMatchType_ComparesCorrectly()
        {
            var match1 = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var match2 = new MockSpectralMatch("file", "SEQ", "SEQ", 50.0, 1); // Different score
            object boxedMatch2 = match2;

            Assert.That(match1.Equals(boxedMatch2), Is.True,
                "Should be equal when boxed as object - score is not part of equality");
        }

        /// <summary>
        /// Verifies Equals(object) returns false for non-BaseSpectralMatch types.
        /// Critical: Prevents incorrect equality with unrelated types.
        /// </summary>
        [Test]
        public void Equals_Object_UnrelatedType_ReturnsFalse()
        {
            var match = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);

            Assert.Multiple(() =>
            {
                Assert.That(match.Equals("string"), Is.False);
                Assert.That(match.Equals(123), Is.False);
                Assert.That(match.Equals(new List<string>()), Is.False);
                Assert.That(match.Equals((object?)null), Is.False);
            });
        }

        #endregion

        #region Equality Operator Tests

        /// <summary>
        /// Verifies == operator correctly compares two BaseSpectralMatch instances.
        /// Critical: Operators must be consistent with Equals() for predictable behavior.
        /// </summary>
        [Test]
        public void EqualityOperator_ComparesCorrectly()
        {
            var match1 = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var match2 = new MockSpectralMatch("file", "SEQ", "SEQ", 50.0, 1); // Different score
            var match3 = new MockSpectralMatch("other", "SEQ", "SEQ", 100.0, 1); // Different file

            Assert.Multiple(() =>
            {
                Assert.That(match1 == match2, Is.True, "Same file/scan/sequence should be equal");
                Assert.That(match1 == match3, Is.False, "Different file should not be equal");
            });
        }

        /// <summary>
        /// Verifies != operator correctly compares two BaseSpectralMatch instances.
        /// Critical: Inequality operator must be logically opposite of equality operator.
        /// </summary>
        [Test]
        public void InequalityOperator_ComparesCorrectly()
        {
            var match1 = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var match2 = new MockSpectralMatch("file", "SEQ", "SEQ", 50.0, 1); // Different score
            var match3 = new MockSpectralMatch("other", "SEQ", "SEQ", 100.0, 1); // Different file

            Assert.Multiple(() =>
            {
                Assert.That(match1 != match2, Is.False, "Same file/scan/sequence should be equal");
                Assert.That(match1 != match3, Is.True, "Different file should not be equal");
            });
        }

        /// <summary>
        /// Verifies == and != operators handle null correctly on both sides.
        /// Critical: Null-safe operators prevent NullReferenceException in comparisons.
        /// </summary>
        [Test]
        public void EqualityOperators_HandleNullCorrectly()
        {
            var match = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            MockSpectralMatch? nullMatch = null;

            Assert.Multiple(() =>
            {
                // Left null
                Assert.That(nullMatch == match, Is.False);
                Assert.That(nullMatch != match, Is.True);

                // Right null
                Assert.That(match == nullMatch, Is.False);
                Assert.That(match != nullMatch, Is.True);

                // Both null
                Assert.That(nullMatch == null, Is.True);
                Assert.That(nullMatch != null, Is.False);
            });
        }

        /// <summary>
        /// Verifies reference equality is handled correctly by operators.
        /// Critical: Same instance should always be equal to itself.
        /// </summary>
        [Test]
        public void EqualityOperators_SameReference_AreEqual()
        {
            var match = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var sameRef = match;

            Assert.Multiple(() =>
            {
                Assert.That(match == sameRef, Is.True);
                Assert.That(match != sameRef, Is.False);
            });
        }

        #endregion

        #region CompareTo Tests

        /// <summary>
        /// Verifies CompareTo orders by Score descending, then FilePath, then ScanNumber.
        /// Critical: Determines PSM ranking for FDR calculation and best-match selection.
        /// </summary>
        [Test]
        public void CompareTo_OrdersByScoreDescendingThenFilePathThenScan()
        {
            var highScore = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var lowScore = new MockSpectralMatch("file", "SEQ", "SEQ", 50.0, 1);

            // Higher score comes last (negative comparison = comes before)
            Assert.That(highScore.CompareTo((ISpectralMatch)lowScore), Is.GreaterThan(0));
            Assert.That(lowScore.CompareTo((ISpectralMatch)highScore), Is.LessThan(0));

            // Equal score - are equal
            var fileA = new MockSpectralMatch("a_file", "SEQ", "SEQ", 100.0, 1);
            var fileB = new MockSpectralMatch("b_file", "SEQ", "SEQ", 100.0, 1);
            Assert.That(fileA.CompareTo((ISpectralMatch)fileB), Is.EqualTo(0));

            // Equal score - compare by scan
            var scan5 = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 5);
            var scan10 = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 10);
            Assert.That(scan5.CompareTo((ISpectralMatch)scan10), Is.GreaterThan(0));

            // Null handling
            Assert.That(highScore.CompareTo((ISpectralMatch?)null), Is.GreaterThan(0));

            // Equal matches
            var equal1 = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var equal2 = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            Assert.That(equal1.CompareTo((ISpectralMatch)equal2), Is.EqualTo(0));
        }

        #endregion

        #region Equality and Collection Tests

        /// <summary>
        /// Verifies equality based on FilePath, ScanNumber, FullSequence (not Score).
        /// Critical: Required for HashSet/Dictionary operations in result aggregation.
        /// </summary>
        [Test]
        public void Equality_BasedOnPathScanSequence_NotScore()
        {
            var a = new MockSpectralMatch("file", "PEP[Mod]TIDE", "PEPTIDE", 100.0, 42);
            var b = new MockSpectralMatch("file", "PEP[Mod]TIDE", "PEPTIDE", 50.0, 42); // Different score

            Assert.That(a.Equals(b), Is.True, "Score should not affect equality");
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));

            // Different file = not equal
            var c = new MockSpectralMatch("other", "PEP[Mod]TIDE", "PEPTIDE", 100.0, 42);
            Assert.That(a.Equals(c), Is.False);

            // Null handling
            Assert.That(a.Equals(null), Is.False);
        }

        /// <summary>
        /// Verifies SpectralMatch works correctly in HashSet and sorting.
        /// Critical: Result aggregation uses HashSets; FDR uses sorted lists.
        /// </summary>
        [Test]
        public void CollectionBehavior_HashSetAndSorting()
        {
            var match1 = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var duplicate = new MockSpectralMatch("file", "SEQ", "SEQ", 50.0, 1); // Same identity, different score
            var different = new MockSpectralMatch("file", "SEQ", "SEQ", 100.0, 2);

            // HashSet deduplication
            var set = new HashSet<MockSpectralMatch> { match1, duplicate, different };
            Assert.That(set.Count, Is.EqualTo(2));

            // Sorting by score descending
            var matches = new List<MockSpectralMatch>
            {
                new("file", "SEQ", "SEQ", 50.0, 1),
                new("file", "SEQ", "SEQ", 100.0, 2),
                new("file", "SEQ", "SEQ", 75.0, 3)
            };
            matches.Sort();

            Assert.That(matches[2].Score, Is.EqualTo(100.0));
            Assert.That(matches[1].Score, Is.EqualTo(75.0));
            Assert.That(matches[0].Score, Is.EqualTo(50.0));
        }

        #endregion

        #region ToString Test

        /// <summary>
        /// Verifies ToString returns expected format for debugging/logging.
        /// </summary>
        [Test]
        public void ToString_ReturnsFormattedString()
        {
            var match = new MockSpectralMatch("file", "PEP[Mod]TIDE", "PEPTIDE", 123.456, 42);

            var result = match.ToString();

            Assert.That(result, Is.EqualTo("Scan 42: PEP[Mod]TIDE (Score: 123.46)"));
        }

        #endregion
    }
}