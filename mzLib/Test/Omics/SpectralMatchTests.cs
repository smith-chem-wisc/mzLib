using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

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
        #region Helper Classes

        private class TestBioPolymerWithSetMods : IBioPolymerWithSetMods
        {
            public string BaseSequence { get; init; }
            public string FullSequence { get; init; }
            public double MostAbundantMonoisotopicMass { get; init; }
            public double MonoisotopicMass => MostAbundantMonoisotopicMass;
            public ChemicalFormula ThisChemicalFormula => new();
            public string SequenceWithChemicalFormulas => FullSequence;
            public int OneBasedStartResidue { get; init; } = 1;
            public int OneBasedEndResidue { get; init; } = 1;
            public int MissedCleavages => 0;
            public string Description => "";
            public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; } = CleavageSpecificity.Full;
            public char PreviousResidue => '-';
            public char NextResidue => '-';
            public IDigestionParams DigestionParams => null!;
            public Dictionary<int, Modification> AllModsOneIsNterminus => new();
            public int NumMods => 0;
            public int NumFixedMods => 0;
            public int NumVariableMods => 0;
            public int Length => BaseSequence?.Length ?? 0;
            public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];
            public IBioPolymer Parent { get; init; }

            public TestBioPolymerWithSetMods(string baseSeq, string fullSeq)
            {
                BaseSequence = baseSeq;
                FullSequence = fullSeq;
            }

            public void Fragment(DissociationType d, FragmentationTerminus t, List<Product> p, FragmentationParams? f = null) { }
            public void FragmentInternally(DissociationType d, int m, List<Product> p, FragmentationParams? f = null) { }
            public IBioPolymerWithSetMods Localize(int i, double m) => this;
            public bool Equals(IBioPolymerWithSetMods? other) => other != null && BaseSequence == other.BaseSequence;
            public override bool Equals(object? obj) => Equals(obj as IBioPolymerWithSetMods);
            public override int GetHashCode() => BaseSequence?.GetHashCode() ?? 0;
        }

        #endregion

        #region Constructor Tests

        /// <summary>
        /// Verifies constructor sets all properties and handles null inputs gracefully.
        /// Critical: PSM data must be correctly stored for downstream analysis.
        /// </summary>
        [Test]
        public void Constructor_SetsAllPropertiesAndHandlesNulls()
        {
            var match = new BaseSpectralMatch("path/file.mzML", 42, 100.5, "PEP[Phospho]TIDE", "PEPTIDE");

            Assert.Multiple(() =>
            {
                Assert.That(match.FullFilePath, Is.EqualTo("path/file.mzML"));
                Assert.That(match.OneBasedScanNumber, Is.EqualTo(42));
                Assert.That(match.Score, Is.EqualTo(100.5));
                Assert.That(match.FullSequence, Is.EqualTo("PEP[Phospho]TIDE"));
                Assert.That(match.BaseSequence, Is.EqualTo("PEPTIDE"));
            });

            // Null handling
            var nullMatch = new BaseSpectralMatch(null!, 1, 10.0, null!, null!);
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
            var polymer1 = new TestBioPolymerWithSetMods("ABC", "ABC");
            var polymer2 = new TestBioPolymerWithSetMods("DEF", "DEF");

            // Constructor with biopolymers
            var match = new BaseSpectralMatch("file", 1, 10.0, "ABC", "ABC", new[] { polymer1 });

            // Add more
            match.AddIdentifiedBioPolymer(polymer2);

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(2));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("ABC"));
            Assert.That(identified[1].BaseSequence, Is.EqualTo("DEF"));

            // Verify defensive copy - original list mutation doesn't affect match
            var source = new List<IBioPolymerWithSetMods> { new TestBioPolymerWithSetMods("X", "X") };
            var match2 = new BaseSpectralMatch("file", 1, 10.0, "X", "X", source);
            source.Add(new TestBioPolymerWithSetMods("Y", "Y"));
            Assert.That(match2.GetIdentifiedBioPolymersWithSetMods().Count(), Is.EqualTo(1));
        }

        /// <summary>
        /// Verifies null biopolymers are filtered out when adding.
        /// Critical: Prevents null reference exceptions in protein grouping.
        /// </summary>
        [Test]
        public void AddIdentifiedBioPolymers_FiltersNulls()
        {
            var match = new BaseSpectralMatch("file", 1, 10.0, "SEQ", "SEQ");

            match.AddIdentifiedBioPolymers(new IBioPolymerWithSetMods?[]
            {
                new TestBioPolymerWithSetMods("A", "A"),
                null,
                new TestBioPolymerWithSetMods("B", "B")
            }!);

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(2));
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
            var match = new BaseSpectralMatch("file", 1, 10.0, "PEPTIDE", "PEPTIDE");

            // Sequential N-term fragments cover positions
            match.GetSequenceCoverage(new List<int> { 1, 2, 3 }, null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(1));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(2));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(3));

            // C-term fragments also work
            var match2 = new BaseSpectralMatch("file", 1, 10.0, "PEPTIDE", "PEPTIDE");
            match2.GetSequenceCoverage(null, new List<int> { 2, 3, 4 });
            Assert.That(match2.FragmentCoveragePositionInPeptide, Contains.Item(1)); // y2 covers first

            // Both termini combined
            var match3 = new BaseSpectralMatch("file", 1, 10.0, "PEPTIDE", "PEPTIDE");
            match3.GetSequenceCoverage(new List<int> { 1, 2 }, new List<int> { 6, 7 });
            Assert.That(match3.FragmentCoveragePositionInPeptide, Contains.Item(1));
            Assert.That(match3.FragmentCoveragePositionInPeptide, Contains.Item(7));
        }

        /// <summary>
        /// Verifies edge cases: empty sequence and no fragments return null coverage.
        /// Critical: Prevents crashes with invalid input.
        /// </summary>
        [Test]
        public void GetSequenceCoverage_EdgeCases_ReturnNull()
        {
            var emptySeq = new BaseSpectralMatch("file", 1, 10.0, "", "");
            emptySeq.GetSequenceCoverage(new List<int> { 1 }, null);
            Assert.That(emptySeq.FragmentCoveragePositionInPeptide, Is.Null);

            var noFragments = new BaseSpectralMatch("file", 1, 10.0, "PEPTIDE", "PEPTIDE");
            noFragments.GetSequenceCoverage(null, null);
            Assert.That(noFragments.FragmentCoveragePositionInPeptide, Is.Null);
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
            var match = new BaseSpectralMatch("file", 1, 10.0, "PEPTIDE", "PEPTIDE");

            // Position 6 = BaseSequence.Length - 1 = 7 - 1 = 6
            match.GetSequenceCoverage(new List<int> { 6 }, null);

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
            var match = new BaseSpectralMatch("file", 1, 10.0, "PEPTIDE", "PEPTIDE");

            // N-term positions 2,3 and C-term position 3 - position 3 should be covered (inclusive)
            match.GetSequenceCoverage(new List<int> { 2, 3 }, new List<int> { 3 });

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
            var match = new BaseSpectralMatch("file", 1, 10.0, "PEPTIDE", "PEPTIDE");

            // N-term positions 2,3 and C-term position 5 (3+2) should cover position 4
            match.GetSequenceCoverage(new List<int> { 2, 3 }, new List<int> { 5 });

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
            var high = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");
            var low = new BaseSpectralMatch("file", 1, 50.0, "SEQ", "SEQ");

            // Cast to IHasSequenceCoverageFromFragments
            IHasSequenceCoverageFromFragments highCoverage = high;
            IHasSequenceCoverageFromFragments lowCoverage = low;

            // Should delegate to ISpectralMatch comparison (score descending)
            Assert.That(high.CompareTo(lowCoverage), Is.LessThan(0),
                "Higher score should come first when comparing via IHasSequenceCoverageFromFragments");
            Assert.That(low.CompareTo(highCoverage), Is.GreaterThan(0));
        }

        /// <summary>
        /// Verifies CompareTo(IHasSequenceCoverageFromFragments) returns 0 for non-ISpectralMatch implementations.
        /// Critical: Prevents crashes when mixed types implement only IHasSequenceCoverageFromFragments.
        /// </summary>
        [Test]
        public void CompareTo_IHasSequenceCoverageFromFragments_NonSpectralMatch_ReturnsZero()
        {
            var match = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");

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
            var match = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");

            Assert.That(match.CompareTo((IHasSequenceCoverageFromFragments?)null), Is.LessThan(0));
        }

        // Helper class for testing non-ISpectralMatch comparison
        private class NonSpectralMatchCoverageProvider : IHasSequenceCoverageFromFragments
        {
            public HashSet<int>? FragmentCoveragePositionInPeptide { get; set; }
            public void GetSequenceCoverage() { }
            public int CompareTo(IHasSequenceCoverageFromFragments? other) => 0;
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
            var match1 = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");
            var match2 = new BaseSpectralMatch("file", 1, 50.0, "SEQ", "SEQ"); // Different score
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
            var match = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");

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
            var match1 = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");
            var match2 = new BaseSpectralMatch("file", 1, 50.0, "SEQ", "SEQ"); // Different score
            var match3 = new BaseSpectralMatch("other", 1, 100.0, "SEQ", "SEQ"); // Different file

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
            var match1 = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");
            var match2 = new BaseSpectralMatch("file", 1, 50.0, "SEQ", "SEQ"); // Different score
            var match3 = new BaseSpectralMatch("other", 1, 100.0, "SEQ", "SEQ"); // Different file

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
            var match = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");
            BaseSpectralMatch? nullMatch = null;

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
            var match = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");
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
            var highScore = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");
            var lowScore = new BaseSpectralMatch("file", 1, 50.0, "SEQ", "SEQ");

            // Higher score comes first (negative comparison = comes before)
            Assert.That(highScore.CompareTo((ISpectralMatch)lowScore), Is.LessThan(0));
            Assert.That(lowScore.CompareTo((ISpectralMatch)highScore), Is.GreaterThan(0));

            // Equal score - compare by file path
            var fileA = new BaseSpectralMatch("a_file", 1, 100.0, "SEQ", "SEQ");
            var fileB = new BaseSpectralMatch("b_file", 1, 100.0, "SEQ", "SEQ");
            Assert.That(fileA.CompareTo((ISpectralMatch)fileB), Is.LessThan(0));

            // Equal score and path - compare by scan
            var scan5 = new BaseSpectralMatch("file", 5, 100.0, "SEQ", "SEQ");
            var scan10 = new BaseSpectralMatch("file", 10, 100.0, "SEQ", "SEQ");
            Assert.That(scan5.CompareTo((ISpectralMatch)scan10), Is.LessThan(0));

            // Null handling
            Assert.That(highScore.CompareTo((ISpectralMatch?)null), Is.LessThan(0));

            // Equal matches
            var equal1 = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");
            var equal2 = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");
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
            var a = new BaseSpectralMatch("file", 42, 100.0, "PEP[Mod]TIDE", "PEPTIDE");
            var b = new BaseSpectralMatch("file", 42, 50.0, "PEP[Mod]TIDE", "PEPTIDE"); // Different score

            Assert.That(a.Equals(b), Is.True, "Score should not affect equality");
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));

            // Different file = not equal
            var c = new BaseSpectralMatch("other", 42, 100.0, "PEP[Mod]TIDE", "PEPTIDE");
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
            var match1 = new BaseSpectralMatch("file", 1, 100.0, "SEQ", "SEQ");
            var duplicate = new BaseSpectralMatch("file", 1, 50.0, "SEQ", "SEQ"); // Same identity, different score
            var different = new BaseSpectralMatch("file", 2, 100.0, "SEQ", "SEQ");

            // HashSet deduplication
            var set = new HashSet<BaseSpectralMatch> { match1, duplicate, different };
            Assert.That(set.Count, Is.EqualTo(2));

            // Sorting by score descending
            var matches = new List<BaseSpectralMatch>
            {
                new("file", 1, 50.0, "SEQ", "SEQ"),
                new("file", 2, 100.0, "SEQ", "SEQ"),
                new("file", 3, 75.0, "SEQ", "SEQ")
            };
            matches.Sort();

            Assert.That(matches[0].Score, Is.EqualTo(100.0));
            Assert.That(matches[1].Score, Is.EqualTo(75.0));
            Assert.That(matches[2].Score, Is.EqualTo(50.0));
        }

        #endregion

        #region ToString Test

        /// <summary>
        /// Verifies ToString returns expected format for debugging/logging.
        /// </summary>
        [Test]
        public void ToString_ReturnsFormattedString()
        {
            var match = new BaseSpectralMatch("file", 42, 123.456, "PEP[Mod]TIDE", "PEPTIDE");

            var result = match.ToString();

            Assert.That(result, Is.EqualTo("Scan 42: PEP[Mod]TIDE (Score: 123.46)"));
        }

        #endregion
    }
}