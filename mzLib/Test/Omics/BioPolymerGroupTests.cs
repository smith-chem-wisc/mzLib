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
    /// Unit tests for BioPolymerGroup class.
    /// Tests core functionality for protein/gene grouping in mass spectrometry analysis.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class BioPolymerGroupTests
    {
        #region Test Data and Setup

        private TestBioPolymer _bioPolymer1;
        private TestBioPolymer _bioPolymer2;
        private TestBioPolymer _decoyBioPolymer;
        private TestBioPolymer _contaminantBioPolymer;
        private TestBioPolymerWithSetMods _sequence1;
        private TestBioPolymerWithSetMods _uniqueSequence;
        private HashSet<IBioPolymer> _bioPolymers;
        private HashSet<IBioPolymerWithSetMods> _allSequences;
        private HashSet<IBioPolymerWithSetMods> _uniqueSequences;
        private BioPolymerGroup _bioPolymerGroup;

        [SetUp]
        public void Setup()
        {
            _bioPolymer1 = new TestBioPolymer("ACGTACGT", "BP12345",
                organism: "Homo sapiens",
                geneNames: new List<Tuple<string, string>> { new("primary", "GENE1") });

            _bioPolymer2 = new TestBioPolymer("TGCATGCA", "BP67890",
                organism: "Homo sapiens",
                geneNames: new List<Tuple<string, string>> { new("primary", "GENE2") });

            _decoyBioPolymer = new TestBioPolymer("DECOYSEQ", "DECOY_BP12345", isDecoy: true);
            _contaminantBioPolymer = new TestBioPolymer("CONTAMINANT", "CONT_BP99999", isContaminant: true);

            _sequence1 = new TestBioPolymerWithSetMods("ACGT", "ACGT");
            _uniqueSequence = new TestBioPolymerWithSetMods("UNIQUE", "UNIQUE");

            _bioPolymers = new HashSet<IBioPolymer> { _bioPolymer1, _bioPolymer2 };
            _allSequences = new HashSet<IBioPolymerWithSetMods> { _sequence1, _uniqueSequence };
            _uniqueSequences = new HashSet<IBioPolymerWithSetMods> { _uniqueSequence };

            _bioPolymerGroup = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
        }

        #endregion

        #region Helper Classes

        private class TestBioPolymer : IBioPolymer
        {
            public string BaseSequence { get; }
            public string Accession { get; }
            public string Organism { get; }
            public string Name { get; }
            public string FullName { get; }
            public List<Tuple<string, string>> GeneNames { get; }
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

            public TestBioPolymer(string sequence, string accession,
                string organism = "", string name = "", string fullName = "",
                List<Tuple<string, string>> geneNames = null,
                bool isDecoy = false, bool isContaminant = false)
            {
                BaseSequence = sequence;
                Accession = accession;
                Organism = organism;
                Name = name;
                FullName = fullName;
                GeneNames = geneNames ?? new List<Tuple<string, string>>();
                IsDecoy = isDecoy;
                IsContaminant = isContaminant;
            }

            public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams,
                List<Modification> allKnownFixedModifications, List<Modification> variableModifications,
                List<SilacLabel>? silacLabels = null, (SilacLabel, SilacLabel)? turnoverLabels = null,
                bool topDownTruncationSearch = false) => Enumerable.Empty<IBioPolymerWithSetMods>();

            public IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence, IDictionary<int, List<Modification>>? newMods)
                => new TestBioPolymer(newBaseSequence, Accession, Organism, Name, FullName, GeneNames, IsDecoy, IsContaminant);

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

            public TestBioPolymerWithSetMods(string baseSequence, string fullSequence, IBioPolymer parent = null)
            {
                BaseSequence = baseSequence;
                FullSequence = fullSequence;
                Parent = parent;
                OneBasedStartResidue = 1;
                OneBasedEndResidue = baseSequence.Length;
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
            private readonly List<IBioPolymerWithSetMods> _identified = new();
            public string FullFilePath { get; }
            public string FullSequence { get; }
            public string BaseSequence { get; }
            public double Score { get; }
            public int OneBasedScanNumber { get; }

            public TestSpectralMatch(string filePath, string baseSequence, string fullSequence, double score, int scanNumber)
            {
                FullFilePath = filePath;
                BaseSequence = baseSequence;
                FullSequence = fullSequence;
                Score = score;
                OneBasedScanNumber = scanNumber;
            }

            public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => _identified;
            public int CompareTo(ISpectralMatch? other) => other is null ? 1 : other.Score.CompareTo(Score);
        }

        #endregion

        #region Constructor Tests

        /// <summary>
        /// Verifies constructor properly initializes all properties and collections.
        /// Critical: Ensures object is in valid state for subsequent operations.
        /// </summary>
        [Test]
        public void Constructor_InitializesAllPropertiesCorrectly()
        {
            var bg = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);

            Assert.Multiple(() =>
            {
                Assert.That(bg.BioPolymers, Is.EqualTo(_bioPolymers));
                Assert.That(bg.AllBioPolymersWithSetMods, Is.EqualTo(_allSequences));
                Assert.That(bg.UniqueBioPolymersWithSetMods, Is.EqualTo(_uniqueSequences));
                Assert.That(bg.BioPolymerGroupScore, Is.EqualTo(0));
                Assert.That(bg.IsDecoy, Is.False);
                Assert.That(bg.IsContaminant, Is.False);
                Assert.That(bg.AllPsmsBelowOnePercentFDR, Is.Empty);
                Assert.That(bg.BioPolymerGroupName, Is.EqualTo("BP12345|BP67890"));
            });
        }

        /// <summary>
        /// Verifies decoy and contaminant flags are set based on member biopolymers.
        /// Critical: These flags determine FDR calculation and result filtering.
        /// </summary>
        [Test]
        [TestCase(true, false, true, false, Description = "Decoy only")]
        [TestCase(false, true, false, true, Description = "Contaminant only")]
        [TestCase(true, true, true, true, Description = "Both decoy and contaminant")]
        public void Constructor_SetsDecoyAndContaminantFlags(bool includeDecoy, bool includeContaminant,
            bool expectedIsDecoy, bool expectedIsContaminant)
        {
            var bioPolymers = new HashSet<IBioPolymer>();
            if (includeDecoy) bioPolymers.Add(_decoyBioPolymer);
            if (includeContaminant) bioPolymers.Add(_contaminantBioPolymer);
            if (!includeDecoy && !includeContaminant) bioPolymers.Add(_bioPolymer1);

            var bg = new BioPolymerGroup(bioPolymers, _allSequences, _uniqueSequences);

            Assert.That(bg.IsDecoy, Is.EqualTo(expectedIsDecoy));
            Assert.That(bg.IsContaminant, Is.EqualTo(expectedIsContaminant));
        }

        /// <summary>
        /// Verifies biopolymers are ordered alphabetically by accession.
        /// Critical: Ensures deterministic, reproducible output across runs.
        /// </summary>
        [Test]
        public void Constructor_OrdersBioPolymersByAccession()
        {
            var orderedList = _bioPolymerGroup.ListOfBioPolymersOrderedByAccession;

            Assert.That(orderedList[0].Accession, Is.EqualTo("BP12345"));
            Assert.That(orderedList[1].Accession, Is.EqualTo("BP67890"));
        }

        #endregion

        #region Score Tests

        /// <summary>
        /// Verifies scoring sums best score per unique base sequence.
        /// Critical: This scoring algorithm determines protein group ranking for FDR.
        /// </summary>
        [Test]
        public void Score_SumsBestScorePerUniqueBaseSequence()
        {
            var psm1 = new TestSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 100, scanNumber: 1);
            var psm2 = new TestSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 150, scanNumber: 2);
            var psm3 = new TestSpectralMatch(@"C:\test.raw", "TGCA", "TGCA", score: 200, scanNumber: 3);

            _bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3 };
            _bioPolymerGroup.Score();

            Assert.That(_bioPolymerGroup.BioPolymerGroupScore, Is.EqualTo(350));
        }

        /// <summary>
        /// Verifies empty PSM collection results in zero score.
        /// Critical: Prevents null reference exceptions and ensures valid state.
        /// </summary>
        [Test]
        public void Score_WithNoPsms_ReturnsZero()
        {
            _bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>();
            _bioPolymerGroup.Score();

            Assert.That(_bioPolymerGroup.BioPolymerGroupScore, Is.EqualTo(0));
        }

        #endregion

        #region MergeWith Tests

        /// <summary>
        /// Verifies merge combines biopolymers, sequences, PSMs and updates group name.
        /// Critical: Protein grouping relies on correct merging during parsimony.
        /// </summary>
        [Test]
        public void MergeWith_CombinesAllCollectionsAndUpdatesName()
        {
            var otherBioPolymer = new TestBioPolymer("MERGESEQ", "A00001");
            var otherSequence = new TestBioPolymerWithSetMods("MERGED", "MERGED");
            var otherPsm = new TestSpectralMatch(@"C:\test.raw", "SEQ", "SEQ", score: 50, scanNumber: 1);

            var otherGroup = new BioPolymerGroup(
                new HashSet<IBioPolymer> { otherBioPolymer },
                new HashSet<IBioPolymerWithSetMods> { otherSequence },
                new HashSet<IBioPolymerWithSetMods> { otherSequence });
            otherGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { otherPsm };
            otherGroup.BioPolymerGroupScore = 100;

            _bioPolymerGroup.MergeWith(otherGroup);

            Assert.Multiple(() =>
            {
                Assert.That(_bioPolymerGroup.BioPolymers.Count, Is.EqualTo(3));
                Assert.That(_bioPolymerGroup.AllBioPolymersWithSetMods.Count, Is.EqualTo(3));
                Assert.That(_bioPolymerGroup.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(1));
                Assert.That(_bioPolymerGroup.BioPolymerGroupName, Does.StartWith("A00001"));
                Assert.That(otherGroup.BioPolymerGroupScore, Is.EqualTo(0), "Merged group score should reset");
            });
        }

        #endregion

        #region Equality Tests

        /// <summary>
        /// Verifies equality is based on BioPolymerGroupName.
        /// Critical: Required for HashSet/Dictionary operations during protein inference.
        /// </summary>
        [Test]
        public void Equality_BasedOnGroupName()
        {
            var bg1 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
            var bg2 = new BioPolymerGroup(_bioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            var bg3 = new BioPolymerGroup(new HashSet<IBioPolymer> { _bioPolymer1 }, _allSequences, _uniqueSequences);

            Assert.Multiple(() =>
            {
                Assert.That(bg1.Equals(bg2), Is.True, "Same biopolymers = same name = equal");
                Assert.That(bg1.Equals(bg3), Is.False, "Different biopolymers = different name = not equal");
                Assert.That(bg1.Equals((BioPolymerGroup)null), Is.False);
                Assert.That(bg1.GetHashCode(), Is.EqualTo(bg2.GetHashCode()));
            });
        }

        /// <summary>
        /// Verifies BioPolymerGroup works correctly in HashSet for deduplication.
        /// Critical: Protein inference uses HashSets extensively.
        /// </summary>
        [Test]
        public void Equality_WorksInHashSet()
        {
            var bg1 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
            var bg2 = new BioPolymerGroup(_bioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            var set = new HashSet<BioPolymerGroup> { bg1, bg2 };

            Assert.That(set.Count, Is.EqualTo(1), "Duplicate groups should be deduplicated");
        }

        #endregion

        #region ConstructSubsetBioPolymerGroup Tests

        /// <summary>
        /// Verifies subset construction filters PSMs and samples to specific file.
        /// Critical: Per-file analysis requires correct data partitioning.
        /// </summary>
        [Test]
        public void ConstructSubsetBioPolymerGroup_FiltersByFile()
        {
            var file1 = new SpectraFileInfo(@"C:\test1.raw", "Control", 1, 1, 0);
            var file2 = new SpectraFileInfo(@"C:\test2.raw", "Control", 1, 1, 0);
            var psm1 = new TestSpectralMatch(@"C:\test1.raw", "SEQ", "SEQ", score: 100, scanNumber: 1);
            var psm2 = new TestSpectralMatch(@"C:\test2.raw", "SEQ", "SEQ", score: 100, scanNumber: 2);

            _bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { file1, file2 };
            _bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { file1, 1000.0 },
                { file2, 2000.0 }
            };
            _bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            var subset = _bioPolymerGroup.ConstructSubsetBioPolymerGroup(@"C:\test1.raw");

            Assert.Multiple(() =>
            {
                Assert.That(subset.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(1));
                Assert.That(subset.SamplesForQuantification.Count, Is.EqualTo(1));
                Assert.That(subset.IntensitiesBySample.Count, Is.EqualTo(1));
                Assert.That(subset.IntensitiesBySample.First().Value, Is.EqualTo(1000.0));
            });
        }

        #endregion

        #region ToString and Header Tests

        /// <summary>
        /// Verifies ToString and GetTabSeparatedHeader have matching column counts.
        /// Critical: Mismatched columns corrupt output files and break downstream tools.
        /// </summary>
        [Test]
        public void ToStringAndHeader_HaveMatchingColumnCounts()
        {
            var header = _bioPolymerGroup.GetTabSeparatedHeader();
            var row = _bioPolymerGroup.ToString();

            var headerColumns = header.Split('\t').Length;
            var rowColumns = row.Split('\t').Length;

            Assert.That(headerColumns, Is.EqualTo(rowColumns));
        }

        /// <summary>
        /// Verifies decoy/contaminant/target markers are correctly included in output.
        /// Critical: These markers determine how results are filtered and reported.
        /// </summary>
        [Test]
        [TestCase(false, false, "\tT\t", Description = "Target")]
        [TestCase(true, false, "\tD\t", Description = "Decoy")]
        [TestCase(false, true, "\tC\t", Description = "Contaminant")]
        public void ToString_ContainsCorrectTypeMarker(bool isDecoy, bool isContaminant, string expectedMarker)
        {
            var bioPolymer = isDecoy ? _decoyBioPolymer : (isContaminant ? _contaminantBioPolymer : _bioPolymer1);
            var bg = new BioPolymerGroup(new HashSet<IBioPolymer> { bioPolymer }, _allSequences, _uniqueSequences);

            var result = bg.ToString();

            Assert.That(result, Does.Contain(expectedMarker));
        }

        #endregion

        #region Quantification Tests

        /// <summary>
        /// Verifies both SpectraFileInfo and IsobaricQuantSampleInfo can be stored for quantification.
        /// Critical: Supports both label-free and TMT/iTRAQ quantification workflows.
        /// </summary>
        [Test]
        public void Quantification_SupportsBothSampleTypes()
        {
            var spectraFile = new SpectraFileInfo(@"C:\test1.raw", "Control", 1, 1, 0);
            var isobaricSample = new IsobaricQuantSampleInfo(@"C:\test2.raw", "Treatment", 1, 1, 0, 1, "126", 126.0, false);

            _bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile, isobaricSample };
            _bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { spectraFile, 1000.0 },
                { isobaricSample, 2000.0 }
            };

            // Verify samples are stored correctly
            Assert.Multiple(() =>
            {
                Assert.That(_bioPolymerGroup.SamplesForQuantification.Count, Is.EqualTo(2));
                Assert.That(_bioPolymerGroup.SamplesForQuantification[0], Is.InstanceOf<SpectraFileInfo>());
                Assert.That(_bioPolymerGroup.SamplesForQuantification[1], Is.InstanceOf<IsobaricQuantSampleInfo>());
                Assert.That(_bioPolymerGroup.IntensitiesBySample[spectraFile], Is.EqualTo(1000.0));
                Assert.That(_bioPolymerGroup.IntensitiesBySample[isobaricSample], Is.EqualTo(2000.0));
            });

            // Verify ToString doesn't throw with mixed types
            Assert.DoesNotThrow(() => _bioPolymerGroup.ToString());
        }

        /// <summary>
        /// Verifies isobaric channels are ordered correctly in output.
        /// Critical: Channel order must match header for correct data alignment.
        /// </summary>
        [Test]
        public void Quantification_IsobaricChannelsOrderedCorrectly()
        {
            var sample126 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var sample127 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 2, "127N", 127.0, false);

            _bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { sample127, sample126 };
            _bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { sample126, 1111.0 },
                { sample127, 2222.0 }
            };

            var result = _bioPolymerGroup.ToString();
            var index1111 = result.IndexOf("1111");
            var index2222 = result.IndexOf("2222");

            Assert.That(index1111, Is.LessThan(index2222), "126 should appear before 127N");
        }

        #endregion

        #region Edge Cases

        /// <summary>
        /// Verifies long strings are truncated to prevent Excel compatibility issues.
        /// Critical: Excel has 32,767 character cell limit; exceeding corrupts files.
        /// </summary>
        [Test]
        public void ToString_TruncatesLongStrings()
        {
            var longName = new string('A', 50000);
            var bioPolymer = new TestBioPolymer("SEQ", "BP00001", fullName: longName);
            var bg = new BioPolymerGroup(new HashSet<IBioPolymer> { bioPolymer }, _allSequences, _uniqueSequences);

            var result = bg.ToString();

            Assert.That(result.Length, Is.LessThan(longName.Length));
        }

        /// <summary>
        /// Verifies null/empty collections are handled without throwing.
        /// Critical: Prevents crashes during edge-case processing.
        /// </summary>
        [Test]
        public void HandlesNullAndEmptyCollections()
        {
            _bioPolymerGroup.SamplesForQuantification = null;
            _bioPolymerGroup.IntensitiesBySample = null;

            Assert.DoesNotThrow(() => _bioPolymerGroup.ToString());
            Assert.DoesNotThrow(() => _bioPolymerGroup.GetTabSeparatedHeader());
            Assert.DoesNotThrow(() => _bioPolymerGroup.ConstructSubsetBioPolymerGroup(@"C:\test.raw"));
        }

        #endregion
    }
}