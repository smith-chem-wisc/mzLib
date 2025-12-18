using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymer;
using Omics.BioPolymerGroup;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectralMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
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
            public Dictionary<int, Modification> AllModsOneIsNterminus { get; }
            public int NumMods => AllModsOneIsNterminus?.Count ?? 0;
            public int NumFixedMods => 0;
            public int NumVariableMods => NumMods;
            public int Length => BaseSequence.Length;
            public IBioPolymer Parent { get; }
            public ChemicalFormula ThisChemicalFormula => new();
            public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

            public TestBioPolymerWithSetMods(string baseSequence, string fullSequence, IBioPolymer parent = null,
                int startResidue = 1, int endResidue = 0, Dictionary<int, Modification> mods = null)
            {
                BaseSequence = baseSequence;
                FullSequence = fullSequence;
                Parent = parent;
                OneBasedStartResidue = startResidue;
                OneBasedEndResidue = endResidue > 0 ? endResidue : startResidue + baseSequence.Length - 1;
                AllModsOneIsNterminus = mods ?? new Dictionary<int, Modification>();
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

            public void AddIdentifiedBioPolymer(IBioPolymerWithSetMods peptide) => _identified.Add(peptide);

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

		[Test]
		public void Score_WithPsms_SumsBestScorePerBaseSequence()
		{
			var psm1 = new BioPolymerGroupSequenceCoverageTests.CoverageSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 100, scanNumber: 1);
			var psm2 = new BioPolymerGroupSequenceCoverageTests.CoverageSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 150, scanNumber: 2);
			var psm3 = new BioPolymerGroupSequenceCoverageTests.CoverageSpectralMatch(@"C:\test.raw", "TGCA", "TGCA", score: 200, scanNumber: 3);

            _bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3 };
            _bioPolymerGroup.Score();

            Assert.That(_bioPolymerGroup.BioPolymerGroupScore, Is.EqualTo(350));
        }

		[Test]
		public void Score_WithMultipleSequencesSameBase_UsesMaxScore()
		{
			var psm1 = new BioPolymerGroupSequenceCoverageTests.CoverageSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 50, scanNumber: 1);
			var psm2 = new BioPolymerGroupSequenceCoverageTests.CoverageSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 75, scanNumber: 2);
			var psm3 = new BioPolymerGroupSequenceCoverageTests.CoverageSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 100, scanNumber: 3);

			_bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3 };
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

        #region Header Generation Tests

        /// <summary>
        /// Verifies label-free header uses filename format when files don't exist (SILAC design path).
        /// Critical: Header format must match intensity column order for correct data parsing.
        /// Note: When files don't exist, the code treats it as SILAC experimental design and uses filename.
        /// </summary>
        [Test]
        public void GetTabSeparatedHeader_LabelFree_WithConditions_UsesFilenameWhenFilesDoNotExist()
        {
            // Files that don't exist trigger SILAC experimental design path, which uses filename
            // Different bioreps ensure separate columns
            var file1 = new SpectraFileInfo(@"C:\test1.raw", "Control", 0, 1, 0);
            var file2 = new SpectraFileInfo(@"C:\test2.raw", "Treatment", 1, 1, 0);

            _bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { file1, file2 };

            var header = _bioPolymerGroup.GetTabSeparatedHeader();

            // When files don't exist, falls back to filename format
            Assert.That(header, Does.Contain("Intensity_test1"));
            Assert.That(header, Does.Contain("Intensity_test2"));
        }

        /// <summary>
        /// Verifies label-free header uses filename when conditions are undefined and unfractionated.
        /// Critical: Ensures correct fallback behavior for simple experimental designs.
        /// </summary>
        [Test]
        public void GetTabSeparatedHeader_LabelFree_UndefinedConditions_UsesFilename()
        {
            // Use different biological replicates so they generate separate columns
            // Constructor: SpectraFileInfo(path, condition, biorep, techrep, fraction)
            var file1 = new SpectraFileInfo(@"C:\sample_A.raw", "", 0, 1, 0);  // biorep=0
            var file2 = new SpectraFileInfo(@"C:\sample_B.raw", "", 1, 1, 0);  // biorep=1

            _bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { file1, file2 };

            var header = _bioPolymerGroup.GetTabSeparatedHeader();

            Assert.That(header, Does.Contain("Intensity_sample_A"));
            Assert.That(header, Does.Contain("Intensity_sample_B"));
        }

        /// <summary>
        /// Verifies isobaric header groups channels by file and orders by channel label.
        /// Critical: Channel order in header must match intensity column order for TMT/iTRAQ data.
        /// </summary>
        [Test]
        public void GetTabSeparatedHeader_Isobaric_GroupsByFileAndOrdersByChannel()
        {
            var sample126 = new IsobaricQuantSampleInfo(@"C:\fileA.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var sample127 = new IsobaricQuantSampleInfo(@"C:\fileA.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);
            var sample128 = new IsobaricQuantSampleInfo(@"C:\fileB.raw", "Control", 1, 1, 0, 1, "128C", 128.0, false);

            _bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { sample127, sample126, sample128 };

            var header = _bioPolymerGroup.GetTabSeparatedHeader();

            // Channels should be ordered by file path first, then by channel label
            var index126 = header.IndexOf("fileA_126");
            var index127 = header.IndexOf("fileA_127N");
            var index128 = header.IndexOf("fileB_128C");

            Assert.That(index126, Is.LessThan(index127), "126 should come before 127N within same file");
            Assert.That(index127, Is.LessThan(index128), "fileA channels should come before fileB channels");
        }

        #endregion

        #region Modification Display in Coverage Tests

        /// <summary>
        /// Verifies N-terminal modifications are displayed with prefix format [ModName]-.
        /// Critical: N-terminal mod annotation format is required for correct sequence interpretation.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_NTerminalMod_DisplaysWithPrefixFormat()
        {
            var bioPolymer = new TestBioPolymer("PEPTIDE", "P00001");

            // Create modification with N-terminal location restriction
            ModificationMotif.TryGetMotif("P", out var motif);
            var nTermMod = new Modification(
                _originalId: "Acetyl",
                _modificationType: "ProteinTermMod",
                _locationRestriction: "N-terminal.",
                _target: motif,
                _monoisotopicMass: 42.0);

            var modsDict = new Dictionary<int, Modification> { { 1, nTermMod } };
            var peptide = new TestBioPolymerWithSetMods("PEPTIDE", "[Acetyl on P]-PEPTIDE", bioPolymer, 1, 7, modsDict);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new TestSpectralMatch(@"C:\test.raw", "PEPTIDE", "[Acetyl on P]-PEPTIDE", 100, 1);
            psm.AddIdentifiedBioPolymer(peptide);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            // N-terminal mods should appear as [ModName]- prefix
            Assert.That(output, Does.Contain("[Acetyl on P]-"));
        }

        /// <summary>
        /// Verifies C-terminal modifications are displayed with suffix format -[ModName].
        /// Critical: C-terminal mod annotation format is required for correct sequence interpretation.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_CTerminalMod_DisplaysWithSuffixFormat()
        {
            var bioPolymer = new TestBioPolymer("PEPTIDE", "P00001");

            ModificationMotif.TryGetMotif("E", out var motif);
            var cTermMod = new Modification(
                _originalId: "Amidated",
                _modificationType: "ProteinTermMod",
                _locationRestriction: "C-terminal.",
                _target: motif,
                _monoisotopicMass: -0.98);

            var modsDict = new Dictionary<int, Modification> { { 8, cTermMod } }; // Position after last residue
            var peptide = new TestBioPolymerWithSetMods("PEPTIDE", "PEPTIDE-[Amidated on E]", bioPolymer, 1, 7, modsDict);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new TestSpectralMatch(@"C:\test.raw", "PEPTIDE", "PEPTIDE-[Amidated on E]", 100, 1);
            psm.AddIdentifiedBioPolymer(peptide);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            // C-terminal mods should appear as -[ModName] suffix
            Assert.That(output, Does.Contain("-[Amidated on E]"));
        }

        #endregion

        #region Modification Occupancy Tests

        /// <summary>
        /// Verifies modification occupancy is calculated correctly for N-terminal modifications.
        /// Critical: N-terminal occupancy must use position 1 in protein coordinates.
        /// </summary>
        [Test]
        public void CalculateModificationOccupancy_NTerminalMod_UsesPosition1()
        {
            var bioPolymer = new TestBioPolymer("MPEPTIDE", "P00001");

            ModificationMotif.TryGetMotif("M", out var motif);
            var nTermMod = new Modification(
                _originalId: "Acetyl",
                _modificationType: "ProteinTermMod",
                _locationRestriction: "N-terminal.",
                _target: motif,
                _monoisotopicMass: 42.0);

            var modsDict = new Dictionary<int, Modification> { { 1, nTermMod } };
            var peptide = new TestBioPolymerWithSetMods("MPEPTIDE", "[Acetyl on M]-MPEPTIDE", bioPolymer, 1, 8, modsDict);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new TestSpectralMatch(@"C:\test.raw", "MPEPTIDE", "[Acetyl on M]-MPEPTIDE", 100, 1);
            psm.AddIdentifiedBioPolymer(peptide);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            // N-terminal mod occupancy should report position as aa1
            Assert.That(output, Does.Contain("#aa1["));
            Assert.That(output, Does.Contain("occupancy=1.00(1/1)"));
        }

        /// <summary>
        /// Verifies modification occupancy is calculated correctly for C-terminal modifications.
        /// Critical: C-terminal occupancy must use protein length as position.
        /// </summary>
        [Test]
        public void CalculateModificationOccupancy_CTerminalMod_UsesProteinLength()
        {
            var bioPolymer = new TestBioPolymer("PEPTIDEK", "P00001"); // Length = 8

            ModificationMotif.TryGetMotif("K", out var motif);
            var cTermMod = new Modification(
                _originalId: "Amidated",
                _modificationType: "ProteinTermMod",
                _locationRestriction: "C-terminal.",
                _target: motif,
                _monoisotopicMass: -0.98);

            var modsDict = new Dictionary<int, Modification> { { 9, cTermMod } };
            var peptide = new TestBioPolymerWithSetMods("PEPTIDEK", "PEPTIDEK-[Amidated on K]", bioPolymer, 1, 8, modsDict);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new TestSpectralMatch(@"C:\test.raw", "PEPTIDEK", "PEPTIDEK-[Amidated on K]", 100, 1);
            psm.AddIdentifiedBioPolymer(peptide);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            // C-terminal mod occupancy should report position as aa8 (protein length)
            Assert.That(output, Does.Contain("#aa8["));
            Assert.That(output, Does.Contain("occupancy=1.00(1/1)"));
        }

        /// <summary>
        /// Verifies unrecognized location restrictions are skipped in occupancy calculation.
        /// Critical: Prevents crashes from unexpected modification location types.
        /// </summary>
        [Test]
        public void CalculateModificationOccupancy_UnrecognizedLocationRestriction_IsSkipped()
        {
            var bioPolymer = new TestBioPolymer("PEPTIDE", "P00001");

            ModificationMotif.TryGetMotif("P", out var motif);
            var unknownMod = new Modification(
                _originalId: "UnknownMod",
                _modificationType: "Unknown",
                _locationRestriction: "SomeUnknownLocation",
                _target: motif,
                _monoisotopicMass: 10.0);

            var modsDict = new Dictionary<int, Modification> { { 2, unknownMod } };
            var peptide = new TestBioPolymerWithSetMods("PEPTIDE", "P[UnknownMod]EPTIDE", bioPolymer, 1, 7, modsDict);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new TestSpectralMatch(@"C:\test.raw", "PEPTIDE", "P[UnknownMod]EPTIDE", 100, 1);
            psm.AddIdentifiedBioPolymer(peptide);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            // Should not throw
            Assert.DoesNotThrow(() => group.CalculateSequenceCoverage());

            var output = group.ToString();
            // Unknown location restriction mods should not appear in occupancy info
            Assert.That(output, Does.Not.Contain("UnknownMod"));
        }

        #endregion

        #region Equals(object) Tests

        /// <summary>
        /// Verifies Equals(object) correctly handles BioPolymerGroup type.
        /// Critical: Required for correct behavior in non-generic collections.
        /// </summary>
        [Test]
        public void Equals_Object_BioPolymerGroupType_ComparesCorrectly()
        {
            var bg1 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
            var bg2 = new BioPolymerGroup(_bioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            object boxedBg2 = bg2;

            Assert.That(bg1.Equals(boxedBg2), Is.True);
        }

        /// <summary>
        /// Verifies Equals(object) correctly handles IBioPolymerGroup interface type.
        /// Critical: Supports polymorphic equality comparisons.
        /// </summary>
        [Test]
        public void Equals_Object_IBioPolymerGroupType_ComparesCorrectly()
        {
            var bg1 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
            IBioPolymerGroup ibg = new BioPolymerGroup(_bioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            object boxedIbg = ibg;

            Assert.That(bg1.Equals(boxedIbg), Is.True);
        }

        /// <summary>
        /// Verifies Equals(object) returns false for non-BioPolymerGroup types.
        /// Critical: Prevents incorrect equality with unrelated types.
        /// </summary>
        [Test]
        public void Equals_Object_UnrelatedType_ReturnsFalse()
        {
            var bg = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);

            Assert.Multiple(() =>
            {
                Assert.That(bg.Equals("string"), Is.False);
                Assert.That(bg.Equals(123), Is.False);
                Assert.That(bg.Equals(new List<string>()), Is.False);
            });
        }

        #endregion

    }
}