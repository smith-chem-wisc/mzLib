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
    [TestFixture]
    internal class SpectralMatchTests
    {
        #region Helper Classes

        /// <summary>
        /// Minimal test implementation of IBioPolymerWithSetMods for SpectralMatch tests.
        /// </summary>
        private class TestBioPolymerWithSetMods : IBioPolymerWithSetMods
        {
            public string BaseSequence { get; init; }
            public string FullSequence { get; init; }
            public double MostAbundantMonoisotopicMass { get; init; }
            public double MonoisotopicMass => MostAbundantMonoisotopicMass;
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
            public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];
            public IBioPolymer Parent { get; init; }

            public TestBioPolymerWithSetMods(string baseSeq, string fullSeq, double mass = 0)
            {
                BaseSequence = baseSeq;
                FullSequence = fullSeq;
                MostAbundantMonoisotopicMass = mass;
            }

            public void Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus,
                List<Product> products, FragmentationParams? fragmentationParams = null)
                => throw new NotImplementedException();

            public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments,
                List<Product> products, FragmentationParams? fragmentationParams = null)
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
                => HashCode.Combine(
                    StringComparer.Ordinal.GetHashCode(BaseSequence ?? string.Empty),
                    StringComparer.Ordinal.GetHashCode(FullSequence ?? string.Empty));
        }

        /// <summary>
        /// Derived SpectralMatch for testing protected members.
        /// </summary>
        private class TestableSpectralMatch : SpectralMatchWithSequenceCoverage
        {
            public TestableSpectralMatch(
                string fullFilePath,
                int oneBasedScanNumber,
                double score,
                string fullSequence,
                string baseSequence,
                IEnumerable<IBioPolymerWithSetMods>? identifiedBioPolymers = null)
                : base(fullFilePath, oneBasedScanNumber, score, fullSequence, baseSequence, identifiedBioPolymers)
            {
            }

            public void SetNTerminalPositions(List<int>? positions)
            {
                // Access protected setter via reflection or expose through derived class
                SetFragmentPositions(positions, null);
            }

            public void SetCTerminalPositions(List<int>? positions)
            {
                SetFragmentPositions(null, positions);
            }

            public void SetBothTerminalPositions(List<int>? nTermPositions, List<int>? cTermPositions)
            {
                SetFragmentPositions(nTermPositions, cTermPositions);
            }
        }

        /// <summary>
        /// Non-SpectralMatch implementation of IHasSequenceCoverageFromFragments for testing interface comparison.
        /// </summary>
        private class NonSpectralMatchCoverageProvider : IHasSequenceCoverageFromFragments
        {
            public HashSet<int>? FragmentCoveragePositionInPeptide { get; private set; }

            public void GetSequenceCoverage() { }

            public int CompareTo(IHasSequenceCoverageFromFragments? other) => 0;
        }

        #endregion

        #region Constructor Tests

        [Test]
        public void Constructor_WithValidParameters_SetsAllProperties()
        {
            var match = new SpectralMatchWithSequenceCoverage("path/to/file.mzML", 42, 100.5, "PEP[Phospho]TIDE", "PEPTIDE");

            Assert.That(match.FullFilePath, Is.EqualTo("path/to/file.mzML"));
            Assert.That(match.OneBasedScanNumber, Is.EqualTo(42));
            Assert.That(match.Score, Is.EqualTo(100.5));
            Assert.That(match.FullSequence, Is.EqualTo("PEP[Phospho]TIDE"));
            Assert.That(match.BaseSequence, Is.EqualTo("PEPTIDE"));
        }

        [Test]
        public void Constructor_WithNullFilePath_SetsEmptyString()
        {
            var match = new SpectralMatchWithSequenceCoverage(null!, 1, 10.0, "SEQ", "SEQ");

            Assert.That(match.FullFilePath, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithNullFullSequence_SetsEmptyString()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, null!, "SEQ");

            Assert.That(match.FullSequence, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithNullBaseSequence_SetsEmptyString()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", null!);

            Assert.That(match.BaseSequence, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithIdentifiedBioPolymers_StoresThem()
        {
            var polymer1 = new TestBioPolymerWithSetMods("ABC", "ABC");
            var polymer2 = new TestBioPolymerWithSetMods("DEF", "DEF");
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "ABC", "ABC", new[] { polymer1, polymer2 });

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();

            Assert.That(identified.Count, Is.EqualTo(2));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("ABC"));
            Assert.That(identified[1].BaseSequence, Is.EqualTo("DEF"));
        }

        [Test]
        public void Constructor_WithNullBioPolymers_CreatesEmptyList()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ", null);

            var identified = match.GetIdentifiedBioPolymersWithSetMods();

            Assert.That(identified, Is.Not.Null);
            Assert.That(identified, Is.Empty);
        }

        [Test]
        public void Constructor_DefensiveCopy_OriginalListMutationDoesNotAffectMatch()
        {
            var source = new List<IBioPolymerWithSetMods>
            {
                new TestBioPolymerWithSetMods("A", "A")
            };
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "A", "A", source);

            // Mutate original list after construction
            source.Add(new TestBioPolymerWithSetMods("B", "B"));

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(1));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("A"));
        }

        #endregion

        #region GetIdentifiedBioPolymersWithSetMods Tests

        [Test]
        public void GetIdentifiedBioPolymersWithSetMods_ReturnsAllProvided()
        {
            var polymers = Enumerable.Range(1, 5)
                .Select(i => new TestBioPolymerWithSetMods($"SEQ{i}", $"SEQ{i}"))
                .ToList();
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ1", "SEQ1", polymers);

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();

            Assert.That(identified.Count, Is.EqualTo(5));
            for (int i = 0; i < 5; i++)
            {
                Assert.That(identified[i].BaseSequence, Is.EqualTo($"SEQ{i + 1}"));
            }
        }

        [Test]
        public void GetIdentifiedBioPolymersWithSetMods_PreservesOrder()
        {
            var polymers = new[]
            {
                new TestBioPolymerWithSetMods("ZEBRA", "ZEBRA"),
                new TestBioPolymerWithSetMods("APPLE", "APPLE"),
                new TestBioPolymerWithSetMods("MANGO", "MANGO")
            };
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ", polymers);

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();

            Assert.That(identified[0].BaseSequence, Is.EqualTo("ZEBRA"));
            Assert.That(identified[1].BaseSequence, Is.EqualTo("APPLE"));
            Assert.That(identified[2].BaseSequence, Is.EqualTo("MANGO"));
        }

        #endregion

        #region AddIdentifiedBioPolymer Tests

        [Test]
        public void AddIdentifiedBioPolymer_AddsToCollection()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");
            var polymer = new TestBioPolymerWithSetMods("NEW", "NEW");

            match.AddIdentifiedBioPolymer(polymer);

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(1));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("NEW"));
        }

        [Test]
        public void AddIdentifiedBioPolymer_WithNull_DoesNotAdd()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");

            match.AddIdentifiedBioPolymer(null!);

            var identified = match.GetIdentifiedBioPolymersWithSetMods();
            Assert.That(identified, Is.Empty);
        }

        [Test]
        public void AddIdentifiedBioPolymer_MultipleCalls_AddsAll()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");

            match.AddIdentifiedBioPolymer(new TestBioPolymerWithSetMods("A", "A"));
            match.AddIdentifiedBioPolymer(new TestBioPolymerWithSetMods("B", "B"));
            match.AddIdentifiedBioPolymer(new TestBioPolymerWithSetMods("C", "C"));

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(3));
        }

        #endregion

        #region AddIdentifiedBioPolymers Tests

        [Test]
        public void AddIdentifiedBioPolymers_AddsAllToCollection()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");
            var polymers = new[]
            {
                new TestBioPolymerWithSetMods("A", "A"),
                new TestBioPolymerWithSetMods("B", "B")
            };

            match.AddIdentifiedBioPolymers(polymers);

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(2));
        }

        [Test]
        public void AddIdentifiedBioPolymers_WithNull_DoesNothing()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");

            match.AddIdentifiedBioPolymers(null!);

            var identified = match.GetIdentifiedBioPolymersWithSetMods();
            Assert.That(identified, Is.Empty);
        }

        [Test]
        public void AddIdentifiedBioPolymers_FiltersOutNullEntries()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");
            var polymers = new IBioPolymerWithSetMods?[]
            {
                new TestBioPolymerWithSetMods("A", "A"),
                null,
                new TestBioPolymerWithSetMods("B", "B")
            };

            match.AddIdentifiedBioPolymers(polymers!);

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(2));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("A"));
            Assert.That(identified[1].BaseSequence, Is.EqualTo("B"));
        }

        #endregion

        #region SetFragmentPositions Tests

        [Test]
        public void SetFragmentPositions_SetsPositionsForCoverageCalculation()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");

            match.SetFragmentPositions(new List<int> { 1, 2 }, new List<int> { 6, 7 });
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
        }

        [Test]
        public void SetFragmentPositions_WithNulls_ClearsPositions()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");

            match.SetFragmentPositions(null, null);
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Null);
        }

        #endregion

        #region GetSequenceCoverage Tests

        [Test]
        public void GetSequenceCoverage_EmptyBaseSequence_ReturnsNull()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "", "");
            match.SetFragmentPositions(new List<int> { 1, 2 }, null);

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Null);
        }

        [Test]
        public void GetSequenceCoverage_NoFragmentPositions_ReturnsNull()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Null);
        }

        [Test]
        public void GetSequenceCoverage_SequentialNTermFragments_CoversResidues()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");
            match.SetFragmentPositions(new List<int> { 1, 2, 3 }, null);

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(1));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(2));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(3));
        }

        [Test]
        public void GetSequenceCoverage_SequentialCTermFragments_CoversResidues()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");
            match.SetFragmentPositions(null, new List<int> { 2, 3, 4 });

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(1)); // y2 covers first
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(2));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(3));
        }

        [Test]
        public void GetSequenceCoverage_NTermAtLastPosition_CoversLastResidue()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");
            // PEPTIDE has 7 characters, so last N-term position is 6 (Length - 1)
            match.SetFragmentPositions(new List<int> { 6 }, null);

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(7));
        }

        [Test]
        public void GetSequenceCoverage_CTermAtLastPosition_CoversLastResidue()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");
            match.SetFragmentPositions(null, new List<int> { 7 });

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(7));
        }

        [Test]
        public void GetSequenceCoverage_BothTerminiCoverage_CombinesResults()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");
            match.SetFragmentPositions(new List<int> { 1, 2 }, new List<int> { 6, 7 });

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(1));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(2));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(7));
        }

        [Test]
        public void GetSequenceCoverage_OverlappingTermini_CoversPosition()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");
            match.SetFragmentPositions(new List<int> { 1, 3 }, new List<int> { 3 });

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(3));
        }

        [Test]
        public void GetSequenceCoverage_CalledMultipleTimes_UpdatesResult()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");

            // First call with some positions
            match.SetFragmentPositions(new List<int> { 1 }, null);
            match.GetSequenceCoverage();
            var firstResult = match.FragmentCoveragePositionInPeptide?.ToHashSet();

            // Second call with different positions
            match.SetFragmentPositions(new List<int> { 1, 2, 3 }, null);
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide.Count, Is.GreaterThan(firstResult?.Count ?? 0));
        }

        #endregion

        #region CompareTo ISpectralMatch Tests

        [Test]
        public void CompareTo_ISpectralMatch_HigherScore_ReturnsNegative()
        {
            var higher = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var lower = new SpectralMatchWithSequenceCoverage("file", 1, 50.0, "SEQ", "SEQ");

            // Higher score should come first (descending), so higher.CompareTo(lower) < 0
            Assert.That(higher.CompareTo((ISpectralMatch)lower), Is.LessThan(0));
            Assert.That(lower.CompareTo((ISpectralMatch)higher), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_ISpectralMatch_EqualScore_ComparesFilePath()
        {
            var a = new SpectralMatchWithSequenceCoverage("a_file", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("b_file", 1, 100.0, "SEQ", "SEQ");

            Assert.That(a.CompareTo((ISpectralMatch)b), Is.LessThan(0));
            Assert.That(b.CompareTo((ISpectralMatch)a), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_ISpectralMatch_EqualScoreAndPath_ComparesScanNumber()
        {
            var scan5 = new SpectralMatchWithSequenceCoverage("file", 5, 100.0, "SEQ", "SEQ");
            var scan10 = new SpectralMatchWithSequenceCoverage("file", 10, 100.0, "SEQ", "SEQ");

            Assert.That(scan5.CompareTo((ISpectralMatch)scan10), Is.LessThan(0));
            Assert.That(scan10.CompareTo((ISpectralMatch)scan5), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_ISpectralMatch_Null_ReturnsNegative()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");

            Assert.That(match.CompareTo((ISpectralMatch?)null), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ISpectralMatch_EqualMatches_ReturnsZero()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");

            Assert.That(a.CompareTo((ISpectralMatch)b), Is.EqualTo(0));
        }

        #endregion

        #region CompareTo IHasSequenceCoverageFromFragments Tests

        [Test]
        public void CompareTo_IHasSequenceCoverageFromFragments_WithSpectralMatch_DelegatesToISpectralMatch()
        {
            var higher = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var lower = new SpectralMatchWithSequenceCoverage("file", 1, 50.0, "SEQ", "SEQ");

            // When comparing SpectralMatch objects via the interface, should use ISpectralMatch comparison
            Assert.That(higher.CompareTo((IHasSequenceCoverageFromFragments)lower), Is.LessThan(0));
            Assert.That(lower.CompareTo((IHasSequenceCoverageFromFragments)higher), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_IHasSequenceCoverageFromFragments_Null_ReturnsNegative()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");

            Assert.That(match.CompareTo((IHasSequenceCoverageFromFragments?)null), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_IHasSequenceCoverageFromFragments_NonSpectralMatch_ReturnsZero()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");
            var nonSpectralMatch = new NonSpectralMatchCoverageProvider();

            // When comparing with a non-ISpectralMatch implementation, should return 0
            Assert.That(match.CompareTo(nonSpectralMatch), Is.EqualTo(0));
        }

        [Test]
        public void CompareTo_IHasSequenceCoverageFromFragments_EqualSpectralMatches_ReturnsZero()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");

            Assert.That(a.CompareTo((IHasSequenceCoverageFromFragments)b), Is.EqualTo(0));
        }

        #endregion

        #region IHasSequenceCoverageFromFragments Interface Tests

        [Test]
        public void SpectralMatch_ImplementsIHasSequenceCoverageFromFragments()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");

            Assert.That(match, Is.InstanceOf<IHasSequenceCoverageFromFragments>());
        }

        [Test]
        public void IHasSequenceCoverageFromFragments_GetSequenceCoverage_CanBeCalledViaInterface()
        {
            IHasSequenceCoverageFromFragments match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");

            // Should not throw
            Assert.DoesNotThrow(() => match.GetSequenceCoverage());
        }

        [Test]
        public void IHasSequenceCoverageFromFragments_Sorting_WorksCorrectly()
        {
            var matches = new List<IHasSequenceCoverageFromFragments>
                {
                    new SpectralMatchWithSequenceCoverage("file", 1, 50.0, "SEQ", "SEQ"),
                    new SpectralMatchWithSequenceCoverage("file", 2, 100.0, "SEQ", "SEQ"),
                    new SpectralMatchWithSequenceCoverage("file", 3, 75.0, "SEQ", "SEQ")
                };

            // Sort using the ISpectralMatch comparison since all items are SpectralMatch
            matches.Sort((x, y) => ((ISpectralMatch)x).CompareTo((ISpectralMatch)y));

            // Higher scores should come first (as SpectralMatch implements descending score order)
            Assert.That(((SpectralMatchWithSequenceCoverage)matches[0]).Score, Is.EqualTo(100.0));
            Assert.That(((SpectralMatchWithSequenceCoverage)matches[1]).Score, Is.EqualTo(75.0));
            Assert.That(((SpectralMatchWithSequenceCoverage)matches[2]).Score, Is.EqualTo(50.0));
        }


        #endregion

        #region Equals Tests

        [Test]
        public void Equals_SameProperties_ReturnsTrue()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 42, 100.0, "PEP[Mod]TIDE", "PEPTIDE");
            var b = new SpectralMatchWithSequenceCoverage("file", 42, 100.0, "PEP[Mod]TIDE", "PEPTIDE");

            Assert.That(a.Equals(b), Is.True);
            Assert.That(b.Equals(a), Is.True);
        }

        [Test]
        public void Equals_DifferentFilePath_ReturnsFalse()
        {
            var a = new SpectralMatchWithSequenceCoverage("file1", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file2", 1, 100.0, "SEQ", "SEQ");

            Assert.That(a.Equals(b), Is.False);
        }

        [Test]
        public void Equals_DifferentScanNumber_ReturnsFalse()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file", 2, 100.0, "SEQ", "SEQ");

            Assert.That(a.Equals(b), Is.False);
        }

        [Test]
        public void Equals_DifferentFullSequence_ReturnsFalse()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ1", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ2", "SEQ");

            Assert.That(a.Equals(b), Is.False);
        }

        [Test]
        public void Equals_DifferentScore_StillReturnsTrue()
        {
            // Score is not part of equality
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file", 1, 50.0, "SEQ", "SEQ");

            Assert.That(a.Equals(b), Is.True);
        }

        [Test]
        public void Equals_Null_ReturnsFalse()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");

            Assert.That(match.Equals(null), Is.False);
        }

        [Test]
        public void Equals_SameReference_ReturnsTrue()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");

            Assert.That(match.Equals(match), Is.True);
        }

        [Test]
        public void Equals_Object_WorksCorrectly()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            object boxed = b;

            Assert.That(a.Equals(boxed), Is.True);
        }

        [Test]
        public void Equals_NonSpectralMatchObject_ReturnsFalse()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "SEQ", "SEQ");

            Assert.That(match.Equals("not a spectral match"), Is.False);
            Assert.That(match.Equals(42), Is.False);
        }

        #endregion

        #region GetHashCode Tests

        [Test]
        public void GetHashCode_EqualObjects_ProduceSameHash()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");

            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFilePath_ProducesDifferentHash()
        {
            var a = new SpectralMatchWithSequenceCoverage("file1", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file2", 1, 100.0, "SEQ", "SEQ");

            Assert.That(a.GetHashCode(), Is.Not.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentScanNumber_ProducesDifferentHash()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file", 2, 100.0, "SEQ", "SEQ");

            Assert.That(a.GetHashCode(), Is.Not.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFullSequence_ProducesDifferentHash()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ1", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ2", "SEQ");

            Assert.That(a.GetHashCode(), Is.Not.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void GetHashCode_ConsistentAcrossCalls()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");

            var hash1 = match.GetHashCode();
            var hash2 = match.GetHashCode();

            Assert.That(hash1, Is.EqualTo(hash2));
        }

        #endregion

        #region Operator Tests

        [Test]
        public void EqualityOperator_EqualObjects_ReturnsTrue()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");

            Assert.That(a == b, Is.True);
        }

        [Test]
        public void EqualityOperator_DifferentObjects_ReturnsFalse()
        {
            var a = new SpectralMatchWithSequenceCoverage("file1", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file2", 1, 100.0, "SEQ", "SEQ");

            Assert.That(a == b, Is.False);
        }

        [Test]
        public void EqualityOperator_BothNull_ReturnsTrue()
        {
            SpectralMatchWithSequenceCoverage? a = null;
            SpectralMatchWithSequenceCoverage? b = null;

            Assert.That(a == b, Is.True);
        }

        [Test]
        public void EqualityOperator_OneNull_ReturnsFalse()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            SpectralMatchWithSequenceCoverage? b = null;

            Assert.That(a == b, Is.False);
            Assert.That(b == a, Is.False);
        }

        [Test]
        public void InequalityOperator_DifferentObjects_ReturnsTrue()
        {
            var a = new SpectralMatchWithSequenceCoverage("file1", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file2", 1, 100.0, "SEQ", "SEQ");

            Assert.That(a != b, Is.True);
        }

        [Test]
        public void InequalityOperator_EqualObjects_ReturnsFalse()
        {
            var a = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var b = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");

            Assert.That(a != b, Is.False);
        }

        #endregion

        #region ToString Tests

        [Test]
        public void ToString_ReturnsFormattedString()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 42, 123.456, "PEP[Mod]TIDE", "PEPTIDE");

            var result = match.ToString();

            Assert.That(result, Does.Contain("42"));
            Assert.That(result, Does.Contain("PEP[Mod]TIDE"));
            Assert.That(result, Does.Contain("123.46")); // Score formatted to 2 decimal places
        }

        [Test]
        public void ToString_Format_MatchesExpected()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 100.00, "SEQ", "SEQ");

            var result = match.ToString();

            Assert.That(result, Is.EqualTo("Scan 1: SEQ (Score: 100.00)"));
        }

        #endregion

        #region FragmentCoveragePositionInPeptide Initial State Tests

        [Test]
        public void FragmentCoveragePositionInPeptide_InitiallyNull()
        {
            var match = new SpectralMatchWithSequenceCoverage("file", 1, 10.0, "PEPTIDE", "PEPTIDE");

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Null);
        }

        #endregion

        #region Integration Tests

        [Test]
        public void FullWorkflow_CreateAddCoverageCompare()
        {
            // Create match
            var match = new SpectralMatchWithSequenceCoverage("test.mzML", 100, 85.5, "PEP[Phospho]TIDE", "PEPTIDE");

            // Add biopolymers
            match.AddIdentifiedBioPolymer(new TestBioPolymerWithSetMods("PEPTIDE", "PEP[Phospho]TIDE", 800.0));

            // Set fragment positions and calculate coverage
            match.SetFragmentPositions(new List<int> { 1, 2, 3 }, new List<int> { 5, 6, 7 });
            match.GetSequenceCoverage();

            // Verify all components
            Assert.That(match.FullFilePath, Is.EqualTo("test.mzML"));
            Assert.That(match.OneBasedScanNumber, Is.EqualTo(100));
            Assert.That(match.Score, Is.EqualTo(85.5));
            Assert.That(match.GetIdentifiedBioPolymersWithSetMods().Count(), Is.EqualTo(1));
            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide.Count, Is.GreaterThan(0));
        }

        [Test]
        public void HashSetBehavior_WorksCorrectly()
        {
            var match1 = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ");
            var match2 = new SpectralMatchWithSequenceCoverage("file", 1, 100.0, "SEQ", "SEQ"); // Equal to match1
            var match3 = new SpectralMatchWithSequenceCoverage("file", 2, 100.0, "SEQ", "SEQ"); // Different scan

            var hashSet = new HashSet<SpectralMatchWithSequenceCoverage> { match1, match2, match3 };

            Assert.That(hashSet.Count, Is.EqualTo(2)); // match1 and match2 are equal
            Assert.That(hashSet.Contains(match1), Is.True);
            Assert.That(hashSet.Contains(match2), Is.True);
            Assert.That(hashSet.Contains(match3), Is.True);
        }

        [Test]
        public void SortingBehavior_OrdersByScoreDescending()
        {
            var matches = new List<SpectralMatchWithSequenceCoverage>
            {
                new SpectralMatchWithSequenceCoverage("file", 1, 50.0, "SEQ", "SEQ"),
                new SpectralMatchWithSequenceCoverage("file", 2, 100.0, "SEQ", "SEQ"),
                new SpectralMatchWithSequenceCoverage("file", 3, 75.0, "SEQ", "SEQ")
            };

            matches.Sort();

            // Higher scores should come first
            Assert.That(matches[0].Score, Is.EqualTo(100.0));
            Assert.That(matches[1].Score, Is.EqualTo(75.0));
            Assert.That(matches[2].Score, Is.EqualTo(50.0));
        }

        #endregion
    }
}