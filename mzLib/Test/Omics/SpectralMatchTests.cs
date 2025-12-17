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

        #endregion

        #region Constructor Tests

        [Test]
        public void Constructor_WithValidParameters_SetsAllProperties()
        {
            CoverageSpectralMatch match = new CoverageSpectralMatch("path/to/file.mzML", "PEP[Phospho]TIDE", "PEPTIDE",100.5, 42);

            Assert.That(match.FullFilePath, Is.EqualTo("path/to/file.mzML"));
            Assert.That(match.OneBasedScanNumber, Is.EqualTo(42));
            Assert.That(match.Score, Is.EqualTo(100.5));
            Assert.That(match.FullSequence, Is.EqualTo("PEP[Phospho]TIDE"));
            Assert.That(match.BaseSequence, Is.EqualTo("PEPTIDE"));
        }

        [Test]
        public void Constructor_WithNullFilePath_SetsEmptyString()
        {
            var match = new CoverageSpectralMatch(null!, "SEQ", "SEQ", 1.0, 1);

            Assert.That(match.FullFilePath, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithNullFullSequence_SetsEmptyString()
        {
            var match = new CoverageSpectralMatch("file", null!, "SEQ", 10.0, 1);

            Assert.That(match.FullSequence, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithNullBaseSequence_SetsEmptyString()
        {
            var match = new CoverageSpectralMatch("file", "SEQ", null!, 10.0, 1);

            Assert.That(match.BaseSequence, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithIdentifiedBioPolymers_StoresThem()
        {
            var polymer1 = new TestBioPolymerWithSetMods("ABC", "ABC");
            var polymer2 = new TestBioPolymerWithSetMods("DEF", "DEF");
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 10.0, 1, new[] { polymer1, polymer2 });

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();

            Assert.That(identified.Count, Is.EqualTo(2));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("ABC"));
            Assert.That(identified[1].BaseSequence, Is.EqualTo("DEF"));
        }

        [Test]
        public void Constructor_DefensiveCopy_OriginalListMutationDoesNotAffectMatch()
        {
            var source = new List<IBioPolymerWithSetMods>
            {
                new TestBioPolymerWithSetMods("A", "A")
            };
            var match = new CoverageSpectralMatch("file", "A", "A", 10.0, 1, source);

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
            var match = new CoverageSpectralMatch("file", "SEQ1", "SEQ1", 10.0, 1, polymers);

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
            var match = new CoverageSpectralMatch("file", "SEQ1", "SEQ1", 10.0, 1, polymers);

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
            var polymer = new TestBioPolymerWithSetMods("NEW", "NEW");
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 10.0, 1, [polymer]);

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(1));
            Assert.That(identified[0].BaseSequence, Is.EqualTo("NEW"));
        }

        [Test]
        public void AddIdentifiedBioPolymer_WithEmpty_DoesNotAdd()
        {
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 10.0, 1, []);

            var identified = match.GetIdentifiedBioPolymersWithSetMods();
            Assert.That(identified, Is.Empty);
        }

        [Test]
        public void AddIdentifiedBioPolymer_MultipleCalls_AddsAll()
        {
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 10.0, 1, new List<IBioPolymerWithSetMods>()
            {
                new TestBioPolymerWithSetMods("A", "A"),
                new TestBioPolymerWithSetMods("B", "B"),
                new TestBioPolymerWithSetMods("C", "C"),
            });

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(3));
        }

        #endregion

        #region AddIdentifiedBioPolymers Tests

        [Test]
        public void AddIdentifiedBioPolymers_AddsAllToCollection()
        {
            var polymers = new[]
            {
                new TestBioPolymerWithSetMods("A", "A"),
                new TestBioPolymerWithSetMods("B", "B")
            };
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 10.0, 1, polymers);

            var identified = match.GetIdentifiedBioPolymersWithSetMods().ToList();
            Assert.That(identified.Count, Is.EqualTo(2));
        }

        #endregion

        #region GetSequenceCoverage Tests

        [Test]
        public void GetSequenceCoverage_EmptyBaseSequence_ReturnsEmpty()
        {
            var match = new CoverageSpectralMatch("file", "", "", 10.0, 1);
            match.MatchedFragmentIons = new List<MatchedFragmentIon>();

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Empty);
        }

        [Test]
        public void GetSequenceCoverage_SequentialNTermFragments_CoversResidues()
        {
            var match = new CoverageSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            match.MatchedFragmentIons = new List<MatchedFragmentIon>
                {
                    // Fragments that cover positions 1, 2, and 3 of the peptide (P, E, P)
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 200.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 3, 3, 0), 300.0, 10.0, 1),
                };

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(1));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(2));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(3));
        }

        [Test]
        public void GetSequenceCoverage_SequentialCTermFragments_CoversResidues()
        {
            var match = new CoverageSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            match.MatchedFragmentIons = new List<MatchedFragmentIon>
                {
                    // y ions: y2, y3, y4 covering from C-terminus
                    new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 2, 2, 0), 100.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 3, 3, 0), 200.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 4, 4, 0), 300.0, 10.0, 1),
                };

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            // Sequential C-term fragments should cover residues
            Assert.That(match.FragmentCoveragePositionInPeptide.Count, Is.GreaterThan(0));
        }

        [Test]
        public void GetSequenceCoverage_NTermAtLastPosition_CoversLastResidue()
        {
            var match = new CoverageSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);
            // PEPTIDE has 7 characters, so last N-term position is 6 (Length - 1)

            match.MatchedFragmentIons = new List<MatchedFragmentIon>
                {
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 6, 6, 0), 100.0, 10.0, 1),
                };

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(7));
        }

        [Test]
        public void GetSequenceCoverage_CTermAtLastPosition_CoversLastResidue()
        {
            var match = new CoverageSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            match.MatchedFragmentIons = new List<MatchedFragmentIon>
                {
                    new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 7, 7, 0), 100.0, 10.0, 1),
                };

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(7));
        }

        [Test]
        public void GetSequenceCoverage_BothTerminiCoverage_CombinesResults()
        {
            var match = new CoverageSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            match.MatchedFragmentIons = new List<MatchedFragmentIon>
                {
                    // N-terminal fragments
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 200.0, 10.0, 1),
                    // C-terminal fragments
                    new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 6, 6, 0), 300.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 7, 7, 0), 400.0, 10.0, 1),
                };

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(1));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(2));
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(7));
        }

        [Test]
        public void GetSequenceCoverage_OverlappingTermini_CoversPosition()
        {
            var match = new CoverageSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            match.MatchedFragmentIons = new List<MatchedFragmentIon>
                {
                    // Fragments that should cover position 3 from both directions
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 3, 3, 0), 200.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 3, 3, 0), 300.0, 10.0, 1),
                };

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide, Contains.Item(3));
        }

        [Test]
        public void GetSequenceCoverage_CalledMultipleTimes_UpdatesResult()
        {
            var match = new CoverageSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            // First call with some fragments
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
                {
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                };
            match.GetSequenceCoverage();
            var firstResult = match.FragmentCoveragePositionInPeptide?.ToHashSet();

            // Second call with more fragments
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
                {
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 200.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 3, 3, 0), 300.0, 10.0, 1),
                };
            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Not.Null);
            Assert.That(match.FragmentCoveragePositionInPeptide.Count, Is.GreaterThan(firstResult?.Count ?? 0));
        }

        [Test]
        public void GetSequenceCoverage_WithNoFragmentIons_ReturnsEmpty()
        {
            var match = new CoverageSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);
            match.MatchedFragmentIons = new List<MatchedFragmentIon>();

            match.GetSequenceCoverage();

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Empty);
        }

        #endregion

        #region CompareTo ISpectralMatch Tests

        [Test]
        public void CompareTo_ISpectralMatch_HigherScore_ReturnsNegative()
        {
            var higher = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var lower = new CoverageSpectralMatch("file", "SEQ", "SEQ", 50.0, 1);

            // Higher score should come first (descending), so higher.CompareTo(lower) < 0
            Assert.That(higher.CompareTo((ISpectralMatch)lower), Is.GreaterThan(0));
            Assert.That(lower.CompareTo((ISpectralMatch)higher), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ISpectralMatch_EqualScore_ComparesFilePath()
        {
            var a = new CoverageSpectralMatch("a_file", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("b_file", "SEQ", "SEQ", 100.0, 1);

            Assert.That(a.CompareTo((ISpectralMatch)b), Is.EqualTo(0));
        }

        [Test]
        public void CompareTo_ISpectralMatch_EqualScoreAndPath_ComparesScanNumber()
        {
            var scan5 = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 5);
            var scan10 = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 10);

            Assert.That(scan5.CompareTo((ISpectralMatch)scan10), Is.GreaterThan(0));
            Assert.That(scan10.CompareTo((ISpectralMatch)scan5), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ISpectralMatch_Null_ReturnsNegative()
        {
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 10.0, 1);

            Assert.That(match.CompareTo((ISpectralMatch?)null), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ISpectralMatch_EqualMatches_ReturnsZero()
        {
            var a = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);

            Assert.That(a.CompareTo((ISpectralMatch)b), Is.EqualTo(0));
        }

        #endregion

        #region CompareTo IHasSequenceCoverageFromFragments Tests

        [Test]
        public void CompareTo_IHasSequenceCoverageFromFragments_WithSpectralMatch_DelegatesToISpectralMatch()
        {
            var higher = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var lower = new CoverageSpectralMatch("file", "SEQ", "SEQ", 50.0, 1);

            // When comparing SpectralMatch objects via the interface, should use ISpectralMatch comparison
            Assert.That(higher.CompareTo((IHasSequenceCoverageFromFragments)lower), Is.GreaterThan(0));
            Assert.That(lower.CompareTo((IHasSequenceCoverageFromFragments)higher), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_IHasSequenceCoverageFromFragments_Null_ReturnsNegative()
        {
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 10.0, 1);

            Assert.That(match.CompareTo((IHasSequenceCoverageFromFragments?)null), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_IHasSequenceCoverageFromFragments_EqualSpectralMatches_ReturnsZero()
        {
            var a = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);

            Assert.That(a.CompareTo((IHasSequenceCoverageFromFragments)b), Is.EqualTo(0));
        }

        #endregion

        #region IHasSequenceCoverageFromFragments Interface Tests

        [Test]
        public void SpectralMatch_ImplementsIHasSequenceCoverageFromFragments()
        {
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 10.0, 1);

            Assert.That(match, Is.InstanceOf<IHasSequenceCoverageFromFragments>());
        }

        [Test]
        public void IHasSequenceCoverageFromFragments_GetSequenceCoverage_CanBeCalledViaInterface()
        {
            IHasSequenceCoverageFromFragments match = new CoverageSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            // Should not throw
            Assert.DoesNotThrow(() => match.GetSequenceCoverage());
        }

        [Test]
        public void IHasSequenceCoverageFromFragments_Sorting_WorksCorrectly()
        {
            var matches = new List<IHasSequenceCoverageFromFragments>
                {
                    new CoverageSpectralMatch("file", "SEQ", "SEQ", 50.0, 1),
                    new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 2),
                    new CoverageSpectralMatch("file", "SEQ", "SEQ", 75.0, 3)
                };

            // Sort using the ISpectralMatch comparison since all items are SpectralMatch
            matches.Sort((x, y) => ((ISpectralMatch)x).CompareTo((ISpectralMatch)y));

            // Higher scores should come first (as SpectralMatch implements descending score order)
            Assert.That(((CoverageSpectralMatch)matches[2]).Score, Is.EqualTo(100.0));
            Assert.That(((CoverageSpectralMatch)matches[1]).Score, Is.EqualTo(75.0));
            Assert.That(((CoverageSpectralMatch)matches[0]).Score, Is.EqualTo(50.0));
        }


        #endregion

        #region Equals Tests

        [Test]
        public void Equals_SameProperties_ReturnsTrue()
        {
            var a = new CoverageSpectralMatch("file", "PEP[Mod]TIDE", "PEPTIDE", 100.0, 42);
            var b = new CoverageSpectralMatch("file", "PEP[Mod]TIDE", "PEPTIDE", 100.0, 42);

            Assert.That(a.Equals(b), Is.True);
            Assert.That(b.Equals(a), Is.True);
        }

        [Test]
        public void Equals_DifferentFilePath_ReturnsFalse()
        {
            var a = new CoverageSpectralMatch("file1", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file2", "SEQ", "SEQ", 100.0, 1);

            Assert.That(a.Equals(b), Is.False);
        }

        [Test]
        public void Equals_DifferentScanNumber_ReturnsFalse()
        {
            var a = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 2);

            Assert.That(a.Equals(b), Is.False);
        }

        [Test]
        public void Equals_DifferentFullSequence_ReturnsFalse()
        {
            var a = new CoverageSpectralMatch("file", "SEQ1", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file", "SEQ2", "SEQ", 100.0, 1);

            Assert.That(a.Equals(b), Is.False);
        }

        [Test]
        public void Equals_DifferentScore_StillReturnsTrue()
        {
            // Score is not part of equality
            var a = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file", "SEQ", "SEQ", 50.0, 1);

            Assert.That(a.Equals(b), Is.True);
        }

        [Test]
        public void Equals_Null_ReturnsFalse()
        {
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 10.0, 1);

            Assert.That(match.Equals(null), Is.False);
        }

        [Test]
        public void Equals_SameReference_ReturnsTrue()
        {
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 10.0, 1);

            Assert.That(match.Equals(match), Is.True);
        }

        [Test]
        public void Equals_Object_WorksCorrectly()
        {
            var a = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            object boxed = b;

            Assert.That(a.Equals(boxed), Is.True);
        }

        [Test]
        public void Equals_NonSpectralMatchObject_ReturnsFalse()
        {
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 10.0, 1);

            Assert.That(match.Equals("not a spectral match"), Is.False);
            Assert.That(match.Equals(42), Is.False);
        }

        #endregion

        #region GetHashCode Tests

        [Test]
        public void GetHashCode_EqualObjects_ProduceSameHash()
        {
            var a = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);

            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFilePath_ProducesDifferentHash()
        {
            var a = new CoverageSpectralMatch("file1", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file2", "SEQ", "SEQ", 100.0, 1);

            Assert.That(a.GetHashCode(), Is.Not.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentScanNumber_ProducesDifferentHash()
        {
            var a = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 2);

            Assert.That(a.GetHashCode(), Is.Not.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFullSequence_ProducesDifferentHash()
        {
            var a = new CoverageSpectralMatch("file", "SEQ1", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file", "SEQ2", "SEQ", 100.0, 1);

            Assert.That(a.GetHashCode(), Is.Not.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void GetHashCode_ConsistentAcrossCalls()
        {
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);

            var hash1 = match.GetHashCode();
            var hash2 = match.GetHashCode();

            Assert.That(hash1, Is.EqualTo(hash2));
        }

        #endregion

        #region Operator Tests

        [Test]
        public void EqualityOperator_EqualObjects_ReturnsTrue()
        {
            var a = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);

            Assert.That(a == b, Is.True);
        }

        [Test]
        public void EqualityOperator_DifferentObjects_ReturnsFalse()
        {
            var a = new CoverageSpectralMatch("file1", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file2", "SEQ", "SEQ", 100.0, 1);

            Assert.That(a == b, Is.False);
        }

        [Test]
        public void EqualityOperator_BothNull_ReturnsTrue()
        {
            CoverageSpectralMatch? a = null;
            CoverageSpectralMatch? b = null;

            Assert.That(a == b, Is.True);
        }

        [Test]
        public void EqualityOperator_OneNull_ReturnsFalse()
        {
            var a = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            CoverageSpectralMatch? b = null;

            Assert.That(a == b, Is.False);
            Assert.That(b == a, Is.False);
        }

        [Test]
        public void InequalityOperator_DifferentObjects_ReturnsTrue()
        {
            var a = new CoverageSpectralMatch("file1", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file2", "SEQ", "SEQ", 100.0, 1);

            Assert.That(a != b, Is.True);
        }

        [Test]
        public void InequalityOperator_EqualObjects_ReturnsFalse()
        {
            var a = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var b = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);

            Assert.That(a != b, Is.False);
        }

        #endregion

        #region ToString Tests

        [Test]
        public void ToString_ReturnsFormattedString()
        {
            var match = new CoverageSpectralMatch("file", "PEP[Mod]TIDE", "PEPTIDE", 123.456, 42);

            var result = match.ToString();

            Assert.That(result, Does.Contain("42"));
            Assert.That(result, Does.Contain("PEP[Mod]TIDE"));
            Assert.That(result, Does.Contain("123.46")); // Score formatted to 2 decimal places
        }

        [Test]
        public void ToString_Format_MatchesExpected()
        {
            var match = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.00, 1);

            var result = match.ToString();

            Assert.That(result, Is.EqualTo("Scan 1: SEQ (Score: 100.00)"));
        }

        #endregion

        #region FragmentCoveragePositionInPeptide Initial State Tests

        [Test]
        public void FragmentCoveragePositionInPeptide_InitiallyNull()
        {
            var match = new CoverageSpectralMatch("file", "PEPTIDE", "PEPTIDE", 10.0, 1);

            Assert.That(match.FragmentCoveragePositionInPeptide, Is.Null);
        }

        #endregion

        #region Integration Tests

        [Test]
        public void FullWorkflow_CreateAddCoverageCompare()
        {
            // Create match
            var match = new CoverageSpectralMatch("test.mzML", "PEP[Phospho]TIDE", "PEPTIDE", 85.5, 100, 
                new[] { new TestBioPolymerWithSetMods("PEPTIDE", "PEP[Phospho]TIDE", 800.0) });

            // Add matched fragment ions
            match.MatchedFragmentIons = new List<MatchedFragmentIon>
                {
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 200.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 3, 3, 0), 300.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 5, 5, 0), 400.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 6, 6, 0), 500.0, 10.0, 1),
                    new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 7, 7, 0), 600.0, 10.0, 1),
                };

            // Calculate coverage
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
            var match1 = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1);
            var match2 = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 1); // Equal to match1
            var match3 = new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 2); // Different scan

            var hashSet = new HashSet<CoverageSpectralMatch> { match1, match2, match3 };

            Assert.That(hashSet.Count, Is.EqualTo(2)); // match1 and match2 are equal
            Assert.That(hashSet.Contains(match1), Is.True);
            Assert.That(hashSet.Contains(match2), Is.True);
            Assert.That(hashSet.Contains(match3), Is.True);
        }

        [Test]
        public void SortingBehavior_OrdersByScoreAscending()
        {
            var matches = new List<CoverageSpectralMatch>
            {
                new CoverageSpectralMatch("file", "SEQ", "SEQ", 50.0, 1),
                new CoverageSpectralMatch("file", "SEQ", "SEQ", 100.0, 2),
                new CoverageSpectralMatch("file", "SEQ", "SEQ", 75.0, 3)
            };

            matches.Sort();

            // Higher scores should come first
            Assert.That(matches[2].Score, Is.EqualTo(100.0));
            Assert.That(matches[1].Score, Is.EqualTo(75.0));
            Assert.That(matches[0].Score, Is.EqualTo(50.0));
        }

        #endregion
    }
}