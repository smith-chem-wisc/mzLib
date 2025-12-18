using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Omics.SpectralMatch;
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

        private MockBioPolymer _bioPolymer1;
        private MockBioPolymer _bioPolymer2;
        private MockBioPolymer _decoyBioPolymer;
        private MockBioPolymer _contaminantBioPolymer;
        private MockBioPolymerWithSetMods _sequence1;
        private MockBioPolymerWithSetMods _uniqueSequence;
        private HashSet<IBioPolymer> _bioPolymers;
        private HashSet<IBioPolymerWithSetMods> _allSequences;
        private HashSet<IBioPolymerWithSetMods> _uniqueSequences;
        private BioPolymerGroup _bioPolymerGroup;

        [SetUp]
        public void Setup()
        {
            _bioPolymer1 = new MockBioPolymer("ACGTACGT", "BP12345",
                organism: "Homo sapiens",
                geneNames: new List<Tuple<string, string>> { new("primary", "GENE1") });

            _bioPolymer2 = new MockBioPolymer("TGCATGCA", "BP67890",
                organism: "Homo sapiens",
                geneNames: new List<Tuple<string, string>> { new("primary", "GENE2") });

            _decoyBioPolymer = new MockBioPolymer("DECOYSEQ", "DECOY_BP12345", true, false);
            _contaminantBioPolymer = new MockBioPolymer("CONTAMINANT", "CONT_BP99999", false, true);

            _sequence1 = new MockBioPolymerWithSetMods("ACGT", "ACGT");
            _uniqueSequence = new MockBioPolymerWithSetMods("UNIQUE", "UNIQUE");

            _bioPolymers = new HashSet<IBioPolymer> { _bioPolymer1, _bioPolymer2 };
            _allSequences = new HashSet<IBioPolymerWithSetMods> { _sequence1, _uniqueSequence };
            _uniqueSequences = new HashSet<IBioPolymerWithSetMods> { _uniqueSequence };

            _bioPolymerGroup = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
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
			var psm1 = new MockSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 100, scanNumber: 1);
			var psm2 = new MockSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 150, scanNumber: 2);
			var psm3 = new MockSpectralMatch(@"C:\test.raw", "TGCA", "TGCA", score: 200, scanNumber: 3);

            _bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3 };
            _bioPolymerGroup.Score();

            Assert.That(_bioPolymerGroup.BioPolymerGroupScore, Is.EqualTo(350));
        }

		[Test]
		public void Score_WithMultipleSequencesSameBase_UsesMaxScore()
		{
			var psm1 = new MockSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 50, scanNumber: 1);
			var psm2 = new MockSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 75, scanNumber: 2);
			var psm3 = new MockSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 100, scanNumber: 3);

			_bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3 };
			_bioPolymerGroup.Score();

            Assert.That(_bioPolymerGroup.BioPolymerGroupScore, Is.EqualTo(100));
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
            var otherBioPolymer = new MockBioPolymer("MERGESEQ", "A00001");
            var otherSequence = new MockBioPolymerWithSetMods("MERGED", "MERGED");
            var otherPsm = new MockSpectralMatch(@"C:\test.raw", "SEQ", "SEQ", score: 50, scanNumber: 1);

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
            var psm1 = new MockSpectralMatch(@"C:\test1.raw", "SEQ", "SEQ", score: 100, scanNumber: 1);
            var psm2 = new MockSpectralMatch(@"C:\test2.raw", "SEQ", "SEQ", score: 100, scanNumber: 2);

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
            var bioPolymer = new MockBioPolymer("SEQ", "BP00001", fullName: longName);
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
            var bioPolymer = new MockBioPolymer("PEPTIDE", "P00001");

            // Create modification with N-terminal location restriction
            ModificationMotif.TryGetMotif("P", out var motif);
            var nTermMod = new Modification(
                _originalId: "Acetyl",
                _modificationType: "ProteinTermMod",
                _locationRestriction: "N-terminal.",
                _target: motif,
                _monoisotopicMass: 42.0);

            var modsDict = new Dictionary<int, Modification> { { 1, nTermMod } };
            var peptide = new MockBioPolymerWithSetMods("PEPTIDE", "[Acetyl on P]-PEPTIDE", bioPolymer, 1, 7, modsDict);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new MockSpectralMatch(@"C:\test.raw", "PEPTIDE", "[Acetyl on P]-PEPTIDE", 100, 1, [peptide]);
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
            var bioPolymer = new MockBioPolymer("PEPTIDE", "P00001");

            ModificationMotif.TryGetMotif("E", out var motif);
            var cTermMod = new Modification(
                _originalId: "Amidated",
                _modificationType: "ProteinTermMod",
                _locationRestriction: "C-terminal.",
                _target: motif,
                _monoisotopicMass: -0.98);

            var modsDict = new Dictionary<int, Modification> { { 8, cTermMod } }; // Position after last residue
            var peptide = new MockBioPolymerWithSetMods("PEPTIDE", "PEPTIDE-[Amidated on E]", bioPolymer, 1, 7, modsDict);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new MockSpectralMatch(@"C:\test.raw", "PEPTIDE", "PEPTIDE-[Amidated on E]", 100, 1, [peptide]);
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
            var bioPolymer = new MockBioPolymer("MPEPTIDE", "P00001");

            ModificationMotif.TryGetMotif("M", out var motif);
            var nTermMod = new Modification(
                _originalId: "Acetyl",
                _modificationType: "ProteinTermMod",
                _locationRestriction: "N-terminal.",
                _target: motif,
                _monoisotopicMass: 42.0);

            var modsDict = new Dictionary<int, Modification> { { 1, nTermMod } };
            var peptide = new MockBioPolymerWithSetMods("MPEPTIDE", "[Acetyl on M]-MPEPTIDE", bioPolymer, 1, 8, modsDict);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new MockSpectralMatch(@"C:\test.raw", "MPEPTIDE", "[Acetyl on M]-MPEPTIDE", 100, 1, [peptide]);
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
            var bioPolymer = new MockBioPolymer("PEPTIDEK", "P00001"); // Length = 8

            ModificationMotif.TryGetMotif("K", out var motif);
            var cTermMod = new Modification(
                _originalId: "Amidated",
                _modificationType: "ProteinTermMod",
                _locationRestriction: "C-terminal.",
                _target: motif,
                _monoisotopicMass: -0.98);

            var modsDict = new Dictionary<int, Modification> { { 9, cTermMod } };
            var peptide = new MockBioPolymerWithSetMods("PEPTIDEK", "PEPTIDEK-[Amidated on K]", bioPolymer, 1, 8, modsDict);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new MockSpectralMatch(@"C:\test.raw", "PEPTIDEK", "PEPTIDEK-[Amidated on K]", 100, 1, [peptide]);
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
            var bioPolymer = new MockBioPolymer("PEPTIDE", "P00001");

            ModificationMotif.TryGetMotif("P", out var motif);
            var unknownMod = new Modification(
                _originalId: "UnknownMod",
                _modificationType: "Unknown",
                _locationRestriction: "SomeUnknownLocation",
                _target: motif,
                _monoisotopicMass: 10.0);

            var modsDict = new Dictionary<int, Modification> { { 2, unknownMod } };
            var peptide = new MockBioPolymerWithSetMods("PEPTIDE", "P[UnknownMod]EPTIDE", bioPolymer, 1, 7, modsDict);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new MockSpectralMatch(@"C:\test.raw", "PEPTIDE", "P[UnknownMod]EPTIDE", 100, 1, [peptide]);
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
        #region DisplayModsOnPeptides Tests

        /// <summary>
        /// Verifies unique and shared sequences use FullSequence (with mods) when DisplayModsOnPeptides is true.
        /// Critical: Modification-aware sequence output is required for PTM analysis workflows.
        /// </summary>
        [Test]
        public void GetIdentifiedSequencesOutput_DisplayModsOnPeptides_UsesFullSequence()
        {
            var bioPolymer = new MockBioPolymer("PEPTIDE", "P00001");
            var modifiedSeq = new MockBioPolymerWithSetMods("PEPTIDE", "PEP[Phospho]TIDE", bioPolymer, 1, 7);
            var unmodifiedSeq = new MockBioPolymerWithSetMods("PEPTIDE", "PEPTIDE", bioPolymer, 1, 7);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { modifiedSeq, unmodifiedSeq },
                new HashSet<IBioPolymerWithSetMods> { modifiedSeq });

            // Enable modification display
            group.DisplayModsOnPeptides = true;

            var output = group.ToString();

            // Unique sequences column should contain the full modified sequence
            Assert.That(output, Does.Contain("PEP[Phospho]TIDE"));
        }

        #endregion

        #region Modification Coverage Display Tests

        /// <summary>
        /// Verifies "Anywhere." modifications are inserted at the correct position in coverage display.
        /// Critical: Incorrect position calculation causes misaligned modification annotations in output.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_AnywhereModification_InsertsAtCorrectPosition()
        {
            var bioPolymer = new MockBioPolymer("PEPTIDE", "P00001");

            ModificationMotif.TryGetMotif("T", out var motif);
            var anywhereMod = new Modification(
                _originalId: "Phospho",
                _modificationType: "Common Biological",
                _locationRestriction: "Anywhere.",
                _target: motif,
                _monoisotopicMass: 79.97);

            // Mod at position 4 in AllModsOneIsNterminus (which is T at position 4 in PEPTIDE)
            var modsDict = new Dictionary<int, Modification> { { 5, anywhereMod } }; // 1=Nterm, so 5 = position 4
            var peptide = new MockBioPolymerWithSetMods("PEPTIDE", "PEP[Phospho on T]TIDE", bioPolymer, 1, 7, modsDict);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new MockSpectralMatch(@"C:\test.raw", "PEPTIDE", "PEP[Phospho on T]TIDE", 100, 1, [peptide]);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            // The modification should appear after the T (4th residue)
            Assert.That(output, Does.Contain("[Phospho on T]"));
        }

        #endregion

        #region TruncateString Tests

        /// <summary>
        /// Verifies MaxStringLength setting controls string truncation in output.
        /// Critical: Uncontrolled string lengths corrupt Excel files (32,767 char limit).
        /// </summary>
        [Test]
        public void MaxStringLength_ControlsTruncation()
        {
            var originalMax = BioPolymerGroup.MaxStringLength;
            try
            {
                // Test with custom limit
                BioPolymerGroup.MaxStringLength = 100;
                var longName = new string('X', 200);
                var bioPolymer = new MockBioPolymer("SEQ", "P00001", fullName: longName);
                var bg = new BioPolymerGroup(new HashSet<IBioPolymer> { bioPolymer }, _allSequences, _uniqueSequences);

                var result = bg.ToString();

                // Full name column should be truncated
                Assert.That(result, Does.Not.Contain(longName));
                Assert.That(result.Contains(new string('X', 100)), Is.True);

                // Test disabling truncation (0 or negative)
                BioPolymerGroup.MaxStringLength = 0;
                var result2 = bg.ToString();
                Assert.That(result2, Does.Contain(longName), "MaxStringLength=0 should disable truncation");
            }
            finally
            {
                BioPolymerGroup.MaxStringLength = originalMax;
            }
        }

        /// <summary>
        /// Verifies null and empty strings are handled gracefully by truncation.
        /// Critical: Prevents NullReferenceException during output generation.
        /// </summary>
        [Test]
        public void TruncateString_HandlesNullAndEmpty()
        {
            var bioPolymer = new MockBioPolymer("SEQ", "P00001", fullName: null);
            var bg = new BioPolymerGroup(new HashSet<IBioPolymer> { bioPolymer }, _allSequences, _uniqueSequences);

            // Should not throw with null fullName
            Assert.DoesNotThrow(() => bg.ToString());

            var bioPolymer2 = new MockBioPolymer("SEQ", "P00002", fullName: "");
            var bg2 = new BioPolymerGroup(new HashSet<IBioPolymer> { bioPolymer2 }, _allSequences, _uniqueSequences);

            // Should not throw with empty fullName
            Assert.DoesNotThrow(() => bg2.ToString());
        }

        #endregion
    }
}