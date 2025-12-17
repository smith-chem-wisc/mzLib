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