using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymer;
using Omics.BioPolymerGroup;
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
    /// Tests for BioPolymerGroup.SequenceCoverageResult nested class.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SequenceCoverageResultTests
    {
        #region Helper Classes

        private class TestBioPolymer : IBioPolymer
        {
            public string BaseSequence { get; }
            public string Accession { get; }
            public string Organism { get; } = "";
            public string Name { get; } = "";
            public string FullName { get; } = "";
            public List<Tuple<string, string>> GeneNames { get; } = new();
            public bool IsDecoy { get; }
            public bool IsContaminant { get; }
            public string DatabaseFilePath { get; } = "";
            public int Length => BaseSequence.Length;
            public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; } = new Dictionary<int, List<Modification>>();
            public string SampleNameForVariants { get; set; } = "";
            public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; } = new Dictionary<int, List<Modification>>();
            public IBioPolymer ConsensusVariant => this;
            public List<SequenceVariation> AppliedSequenceVariations { get; } = new();
            public List<SequenceVariation> SequenceVariations { get; } = new();
            public List<TruncationProduct> TruncationProducts { get; } = new();

            public TestBioPolymer(string sequence, string accession, bool isDecoy = false, bool isContaminant = false)
            {
                BaseSequence = sequence;
                Accession = accession;
                IsDecoy = isDecoy;
                IsContaminant = isContaminant;
            }

            public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams,
                List<Modification> allKnownFixedModifications, List<Modification> variableModifications,
                List<SilacLabel>? silacLabels = null, (SilacLabel, SilacLabel)? turnoverLabels = null,
                bool topDownTruncationSearch = false) => Enumerable.Empty<IBioPolymerWithSetMods>();

            public IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence, IDictionary<int, List<Modification>>? newMods)
                => new TestBioPolymer(newBaseSequence, Accession, IsDecoy, IsContaminant);

            public TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original,
                IEnumerable<SequenceVariation> appliedSequenceVariants, IEnumerable<TruncationProduct> applicableProteolysisProducts,
                IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
                where TBioPolymerType : IHasSequenceVariants => original;

            public bool Equals(IBioPolymer? other) => other != null && Accession == other.Accession && BaseSequence == other.BaseSequence;
            public override bool Equals(object? obj) => obj is IBioPolymer other && Equals(other);
            public override int GetHashCode() => HashCode.Combine(Accession, BaseSequence);
        }

        private class TestBioPolymerWithSetMods : IBioPolymerWithSetMods
        {
            public string BaseSequence { get; }
            public string FullSequence { get; }
            public double MostAbundantMonoisotopicMass { get; } = 0;
            public double MonoisotopicMass { get; } = 0;
            public string SequenceWithChemicalFormulas => BaseSequence;
            public int OneBasedStartResidue { get; }
            public int OneBasedEndResidue { get; }
            public int MissedCleavages => 0;
            public string Description => "";
            public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; } = CleavageSpecificity.Full;
            public char PreviousResidue => '-';
            public char NextResidue => '-';
            public IDigestionParams DigestionParams => null!;
            public Dictionary<int, Modification> AllModsOneIsNterminus { get; } = new();
            public int NumMods => 0;
            public int NumFixedMods => 0;
            public int NumVariableMods => 0;
            public int Length => BaseSequence.Length;
            public IBioPolymer Parent { get; }
            public ChemicalFormula ThisChemicalFormula => new();
            public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

            public TestBioPolymerWithSetMods(string baseSequence, string fullSequence, IBioPolymer parent,
                int startResidue = 1, int endResidue = 0)
            {
                BaseSequence = baseSequence;
                FullSequence = fullSequence;
                Parent = parent;
                OneBasedStartResidue = startResidue;
                OneBasedEndResidue = endResidue > 0 ? endResidue : startResidue + baseSequence.Length - 1;
            }

            public void Fragment(DissociationType d, FragmentationTerminus t, List<Product> p, FragmentationParams? f = null) { }
            public void FragmentInternally(DissociationType d, int m, List<Product> p, FragmentationParams? f = null) { }
            public IBioPolymerWithSetMods Localize(int i, double m) => this;
            public bool Equals(IBioPolymerWithSetMods? other) => other != null && BaseSequence == other.BaseSequence;
            public override bool Equals(object? obj) => obj is IBioPolymerWithSetMods other && Equals(other);
            public override int GetHashCode() => BaseSequence.GetHashCode();
        }

        private class TestSpectralMatch : ISpectralMatch
        {
            private readonly List<IBioPolymerWithSetMods> _identified;
            public string FullFilePath { get; }
            public string FullSequence { get; }
            public string BaseSequence { get; }
            public double Score { get; }
            public int OneBasedScanNumber { get; }

            public TestSpectralMatch(string filePath, string baseSequence, string fullSequence, double score, int scanNumber,
                IEnumerable<IBioPolymerWithSetMods> identified = null)
            {
                FullFilePath = filePath;
                BaseSequence = baseSequence;
                FullSequence = fullSequence;
                Score = score;
                OneBasedScanNumber = scanNumber;
                _identified = identified?.ToList() ?? new List<IBioPolymerWithSetMods>();
            }

            public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => _identified;
            public int CompareTo(ISpectralMatch? other) => other is null ? 1 : other.Score.CompareTo(Score);
        }

        #endregion

        #region Constructor Tests

        /// <summary>
        /// Verifies constructor initializes all lists as empty, non-null collections.
        /// Critical: Prevents null reference exceptions when populating coverage data.
        /// </summary>
        [Test]
        public void Constructor_InitializesEmptyLists()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            Assert.Multiple(() =>
            {
                Assert.That(result.SequenceCoverageFraction, Is.Not.Null.And.Empty);
                Assert.That(result.SequenceCoverageDisplayList, Is.Not.Null.And.Empty);
                Assert.That(result.SequenceCoverageDisplayListWithMods, Is.Not.Null.And.Empty);
                Assert.That(result.FragmentSequenceCoverageDisplayList, Is.Not.Null.And.Empty);
                Assert.That(result.ModsInfo, Is.Not.Null.And.Empty);
            });
        }

        /// <summary>
        /// Verifies each instance has independent lists.
        /// Critical: Prevents cross-contamination between protein groups.
        /// </summary>
        [Test]
        public void Constructor_ListsAreIndependent()
        {
            var result1 = new BioPolymerGroup.SequenceCoverageResult();
            var result2 = new BioPolymerGroup.SequenceCoverageResult();

            result1.SequenceCoverageFraction.Add(0.5);
            result1.SequenceCoverageDisplayList.Add("TEST");

            Assert.That(result2.SequenceCoverageFraction, Is.Empty);
            Assert.That(result2.SequenceCoverageDisplayList, Is.Empty);
        }

        #endregion

        #region List Functionality Tests

        /// <summary>
        /// Verifies all lists support standard List operations.
        /// Critical: Lists must be mutable for CalculateSequenceCoverage to populate them.
        /// </summary>
        [Test]
        public void Lists_SupportStandardOperations()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            // Add
            result.SequenceCoverageFraction.Add(0.5);
            result.SequenceCoverageDisplayList.Add("ACDEFGHIK");
            result.SequenceCoverageDisplayListWithMods.Add("[Acetyl]-ACDEF");
            result.FragmentSequenceCoverageDisplayList.Add("ACDefghik");
            result.ModsInfo.Add("#aa3[Phospho,info:occupancy=0.50(1/2)]");

            Assert.Multiple(() =>
            {
                Assert.That(result.SequenceCoverageFraction.Count, Is.EqualTo(1));
                Assert.That(result.SequenceCoverageDisplayList[0], Is.EqualTo("ACDEFGHIK"));
                Assert.That(result.SequenceCoverageDisplayListWithMods[0], Does.Contain("[Acetyl]"));
                Assert.That(result.FragmentSequenceCoverageDisplayList[0], Is.EqualTo("ACDefghik"));
                Assert.That(result.ModsInfo[0], Does.Contain("occupancy"));
            });
        }

        #endregion

        #region Integration Tests

        /// <summary>
        /// Verifies SequenceCoverageResult is correctly populated by CalculateSequenceCoverage.
        /// Critical: Ensures coverage calculation produces valid output for ToString.
        /// </summary>
        [Test]
        public void IntegrationWithBioPolymerGroup_PopulatesCorrectly()
        {
            var bioPolymer = new TestBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new TestBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer, 1, 5);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new TestSpectralMatch(@"C:\test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();

            // Coverage fraction should be 5/9 ≈ 0.556
            Assert.That(output, Does.Contain("0.5"));
            // Coverage display should show uppercase for covered residues
            Assert.That(output, Does.Contain("ACDEF"));
        }

        /// <summary>
        /// Verifies ToString works before CalculateSequenceCoverage is called.
        /// Critical: Prevents crashes when coverage hasn't been calculated yet.
        /// </summary>
        [Test]
        public void IntegrationWithBioPolymerGroup_ToStringWorksBeforeCalculation()
        {
            var bioPolymer = new TestBioPolymer("ACDEFGHIK", "P00001");

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            Assert.DoesNotThrow(() => group.ToString());
        }

        /// <summary>
        /// Verifies repeated CalculateSequenceCoverage calls replace rather than append results.
        /// Critical: Prevents duplicate data accumulation.
        /// </summary>
        [Test]
        public void IntegrationWithBioPolymerGroup_ReplacesOnRecalculation()
        {
            var bioPolymer = new TestBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new TestBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer, 1, 5);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new TestSpectralMatch(@"C:\test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();
            var output1 = group.ToString();

            group.CalculateSequenceCoverage();
            var output2 = group.ToString();

            Assert.That(output1, Is.EqualTo(output2));
        }

        #endregion
    }
}