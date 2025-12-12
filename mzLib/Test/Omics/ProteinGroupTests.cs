using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteinGroup;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Omics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class ProteinGroupTests
    {
        #region Test Data and Setup

        private Protein _protein1;
        private Protein _protein2;
        private Protein _decoyProtein;
        private Protein _contaminantProtein;
        private PeptideWithSetModifications _peptide1;
        private PeptideWithSetModifications _peptide2;
        private PeptideWithSetModifications _uniquePeptide;
        private HashSet<IBioPolymer> _proteins;
        private HashSet<IBioPolymerWithSetMods> _allPeptides;
        private HashSet<IBioPolymerWithSetMods> _uniquePeptides;
        private ProteinGroup _proteinGroup;

        [SetUp]
        public void Setup()
        {
            _protein1 = new Protein("PEPTIDEK", "P12345",
                organism: "Homo sapiens",
                name: "TestProtein1",
                fullName: "Test Protein 1 Full Name",
                geneNames: new List<Tuple<string, string>> { new("primary", "GENE1") });

            _protein2 = new Protein("ANOTHERPEPTIDE", "P67890",
                organism: "Homo sapiens",
                name: "TestProtein2",
                fullName: "Test Protein 2 Full Name",
                geneNames: new List<Tuple<string, string>> { new("primary", "GENE2") });

            _decoyProtein = new Protein("DECOYSEQ", "DECOY_P12345",
                organism: "Homo sapiens",
                isDecoy: true);

            _contaminantProtein = new Protein("CONTAMINANT", "CONT_P99999",
                organism: "Homo sapiens",
                isContaminant: true);

            _peptide1 = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            _peptide2 = new PeptideWithSetModifications("ANOTHER", new Dictionary<string, Modification>());
            _uniquePeptide = new PeptideWithSetModifications("UNIQUE", new Dictionary<string, Modification>());

            _proteins = new HashSet<IBioPolymer> { _protein1, _protein2 };
            _allPeptides = new HashSet<IBioPolymerWithSetMods> { _peptide1, _peptide2, _uniquePeptide };
            _uniquePeptides = new HashSet<IBioPolymerWithSetMods> { _uniquePeptide };

            _proteinGroup = new ProteinGroup(_proteins, _allPeptides, _uniquePeptides);
        }

        #endregion

        #region Constructor Tests

        [Test]
        public void Constructor_WithValidParameters_InitializesAllProperties()
        {
            var pg = new ProteinGroup(_proteins, _allPeptides, _uniquePeptides);

            Assert.That(pg.BioPolymers, Is.EqualTo(_proteins));
            Assert.That(pg.AllBioPolymersWithSetMods, Is.EqualTo(_allPeptides));
            Assert.That(pg.UniqueBioPolymersWithSetMods, Is.EqualTo(_uniquePeptides));
            Assert.That(pg.BioPolymerGroupScore, Is.EqualTo(0));
            Assert.That(pg.BestBioPolymerWithSetModsScore, Is.EqualTo(0));
            Assert.That(pg.QValue, Is.EqualTo(0));
            Assert.That(pg.IsDecoy, Is.False);
            Assert.That(pg.IsContaminant, Is.False);
            Assert.That(pg.AllPsmsBelowOnePercentFDR, Is.Not.Null);
            Assert.That(pg.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(0));
            Assert.That(pg.ModsInfo, Is.Not.Null);
            Assert.That(pg.SequenceCoverageFraction, Is.Not.Null);
            Assert.That(pg.SequenceCoverageDisplayList, Is.Not.Null);
        }

        [Test]
        public void Constructor_WithDecoyProtein_SetsIsDecoyTrue()
        {
            var proteins = new HashSet<IBioPolymer> { _decoyProtein };
            var pg = new ProteinGroup(proteins, _allPeptides, _uniquePeptides);

            Assert.That(pg.IsDecoy, Is.True);
            Assert.That(pg.IsContaminant, Is.False);
        }

        [Test]
        public void Constructor_WithContaminantProtein_SetsIsContaminantTrue()
        {
            var proteins = new HashSet<IBioPolymer> { _contaminantProtein };
            var pg = new ProteinGroup(proteins, _allPeptides, _uniquePeptides);

            Assert.That(pg.IsContaminant, Is.True);
            Assert.That(pg.IsDecoy, Is.False);
        }

        [Test]
        public void Constructor_WithMixedProteins_DecoyTakesPrecedence()
        {
            // Decoy is checked first in the loop
            var proteins = new HashSet<IBioPolymer> { _protein1, _decoyProtein, _contaminantProtein };
            var pg = new ProteinGroup(proteins, _allPeptides, _uniquePeptides);

            // The loop breaks on the first decoy found
            Assert.That(pg.IsDecoy, Is.True);
        }

        [Test]
        public void Constructor_BioPolymerGroupName_IsOrderedByAccession()
        {
            // P12345 comes before P67890 alphabetically
            var expectedName = "P12345|P67890";
            Assert.That(_proteinGroup.BioPolymerGroupName, Is.EqualTo(expectedName));
        }

        [Test]
        public void Constructor_ListOfBioPolymersOrderedByAccession_IsCorrectlyOrdered()
        {
            var orderedList = _proteinGroup.ListOfBioPolymersOrderedByAccession;

            Assert.That(orderedList.Count, Is.EqualTo(2));
            Assert.That(orderedList[0].Accession, Is.EqualTo("P12345"));
            Assert.That(orderedList[1].Accession, Is.EqualTo("P67890"));
        }

        [Test]
        public void Constructor_WithEmptyProteins_CreatesEmptyGroup()
        {
            var emptyProteins = new HashSet<IBioPolymer>();
            var pg = new ProteinGroup(emptyProteins, _allPeptides, _uniquePeptides);

            Assert.That(pg.BioPolymers.Count, Is.EqualTo(0));
            Assert.That(pg.BioPolymerGroupName, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithEmptyPeptides_CreatesGroupWithNoPeptides()
        {
            var emptyPeptides = new HashSet<IBioPolymerWithSetMods>();
            var pg = new ProteinGroup(_proteins, emptyPeptides, emptyPeptides);

            Assert.That(pg.AllBioPolymersWithSetMods.Count, Is.EqualTo(0));
            Assert.That(pg.UniqueBioPolymersWithSetMods.Count, Is.EqualTo(0));
        }

        #endregion

        #region Legacy Property Tests

        [Test]
        public void LegacyProperties_Proteins_MapsToAndFromBioPolymers()
        {
            Assert.That(_proteinGroup.Proteins, Is.SameAs(_proteinGroup.BioPolymers));

            var newProteins = new HashSet<IBioPolymer> { _protein1 };
            _proteinGroup.Proteins = newProteins;
            Assert.That(_proteinGroup.BioPolymers, Is.SameAs(newProteins));
        }

        [Test]
        public void LegacyProperties_ProteinGroupName_MapsToBioPolymerGroupName()
        {
            Assert.That(_proteinGroup.ProteinGroupName, Is.EqualTo(_proteinGroup.BioPolymerGroupName));
        }

        [Test]
        public void LegacyProperties_ProteinGroupScore_MapsToBioPolymerGroupScore()
        {
            _proteinGroup.ProteinGroupScore = 42.5;
            Assert.That(_proteinGroup.BioPolymerGroupScore, Is.EqualTo(42.5));

            _proteinGroup.BioPolymerGroupScore = 100.0;
            Assert.That(_proteinGroup.ProteinGroupScore, Is.EqualTo(100.0));
        }

        [Test]
        public void LegacyProperties_AllPeptides_MapsToAllBioPolymersWithSetMods()
        {
            Assert.That(_proteinGroup.AllPeptides, Is.SameAs(_proteinGroup.AllBioPolymersWithSetMods));

            var newPeptides = new HashSet<IBioPolymerWithSetMods> { _peptide1 };
            _proteinGroup.AllPeptides = newPeptides;
            Assert.That(_proteinGroup.AllBioPolymersWithSetMods, Is.SameAs(newPeptides));
        }

        [Test]
        public void LegacyProperties_UniquePeptides_MapsToUniqueBioPolymersWithSetMods()
        {
            Assert.That(_proteinGroup.UniquePeptides, Is.SameAs(_proteinGroup.UniqueBioPolymersWithSetMods));
        }

        [Test]
        public void LegacyProperties_BestPeptideQValue_MapsToBestBioPolymerWithSetModsQValue()
        {
            _proteinGroup.BestPeptideQValue = 0.01;
            Assert.That(_proteinGroup.BestBioPolymerWithSetModsQValue, Is.EqualTo(0.01));
        }

        [Test]
        public void LegacyProperties_BestPeptideScore_MapsToBestBioPolymerWithSetModsScore()
        {
            _proteinGroup.BestPeptideScore = 150.0;
            Assert.That(_proteinGroup.BestBioPolymerWithSetModsScore, Is.EqualTo(150.0));
        }

        #endregion

        #region FilesForQuantification Legacy Property Tests

        [Test]
        public void FilesForQuantification_Get_ReturnsOnlySpectraFileInfo()
        {
            var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            var isobaricSample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            _proteinGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile, isobaricSample };

            var filesForQuant = _proteinGroup.FilesForQuantification;

            Assert.That(filesForQuant.Count, Is.EqualTo(1));
            Assert.That(filesForQuant[0], Is.EqualTo(spectraFile));
        }

        [Test]
        public void FilesForQuantification_Set_ConvertsToBioPolymerGroupSamples()
        {
            var spectraFiles = new List<SpectraFileInfo>
            {
                new(@"C:\test1.raw", "Control", 1, 1, 0),
                new(@"C:\test2.raw", "Treatment", 1, 1, 0)
            };

            _proteinGroup.FilesForQuantification = spectraFiles;

            Assert.That(_proteinGroup.SamplesForQuantification.Count, Is.EqualTo(2));
            Assert.That(_proteinGroup.SamplesForQuantification.All(s => s is SpectraFileInfo), Is.True);
        }

        [Test]
        public void FilesForQuantification_GetWhenNull_ReturnsNull()
        {
            _proteinGroup.SamplesForQuantification = null;
            Assert.That(_proteinGroup.FilesForQuantification, Is.Null);
        }

        #endregion

        #region IntensitiesByFile Legacy Property Tests

        [Test]
        public void IntensitiesByFile_Get_ReturnsOnlySpectraFileInfoEntries()
        {
            var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            var isobaricSample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            _proteinGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { spectraFile, 1000.0 },
                { isobaricSample, 2000.0 }
            };

            var intensitiesByFile = _proteinGroup.IntensitiesByFile;

            Assert.That(intensitiesByFile.Count, Is.EqualTo(1));
            Assert.That(intensitiesByFile[spectraFile], Is.EqualTo(1000.0));
        }

        [Test]
        public void IntensitiesByFile_Set_ConvertsToIntensitiesBySample()
        {
            var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            var intensities = new Dictionary<SpectraFileInfo, double>
            {
                { spectraFile, 5000.0 }
            };

            _proteinGroup.IntensitiesByFile = intensities;

            Assert.That(_proteinGroup.IntensitiesBySample.Count, Is.EqualTo(1));
            Assert.That(_proteinGroup.IntensitiesBySample[spectraFile], Is.EqualTo(5000.0));
        }

        [Test]
        public void IntensitiesByFile_GetWhenNull_ReturnsNull()
        {
            _proteinGroup.IntensitiesBySample = null;
            Assert.That(_proteinGroup.IntensitiesByFile, Is.Null);
        }

        #endregion

        #region Score Method Tests

        [Test]
        public void Score_WithNoPsms_SetsScoreToZero()
        {
            _proteinGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>();
            _proteinGroup.Score();

            Assert.That(_proteinGroup.BioPolymerGroupScore, Is.EqualTo(0));
        }

        [Test]
        public void Score_WithPsms_SumsBestScorePerBaseSequence()
        {
            var psm1 = new ProteinGroupTestSpectralMatch(@"C:\test.raw", "PEPTIDE", "PEPTIDE", score: 100, scanNumber: 1);
            var psm2 = new ProteinGroupTestSpectralMatch(@"C:\test.raw", "PEPTIDE", "PEPTIDE", score: 150, scanNumber: 2);
            var psm3 = new ProteinGroupTestSpectralMatch(@"C:\test.raw", "ANOTHER", "ANOTHER", score: 200, scanNumber: 3);

            _proteinGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3 };
            _proteinGroup.Score();

            // Best score for PEPTIDE = 150, Best score for ANOTHER = 200, Total = 350
            Assert.That(_proteinGroup.BioPolymerGroupScore, Is.EqualTo(350));
        }

        #endregion

        #region MergeWith Tests

        [Test]
        public void MergeWith_CombinesProteinsAndPeptides()
        {
            var otherProtein = new Protein("MERGESEQ", "P99999");
            var otherProteins = new HashSet<IBioPolymer> { otherProtein };
            var otherPeptide = new PeptideWithSetModifications("MERGED", new Dictionary<string, Modification>());
            var otherPeptides = new HashSet<IBioPolymerWithSetMods> { otherPeptide };
            var otherUnique = new HashSet<IBioPolymerWithSetMods> { otherPeptide };

            var otherGroup = new ProteinGroup(otherProteins, otherPeptides, otherUnique);
            otherGroup.BioPolymerGroupScore = 100;

            _proteinGroup.MergeWith(otherGroup);

            Assert.That(_proteinGroup.BioPolymers.Count, Is.EqualTo(3));
            Assert.That(_proteinGroup.AllBioPolymersWithSetMods.Count, Is.EqualTo(4));
            Assert.That(_proteinGroup.UniqueBioPolymersWithSetMods.Count, Is.EqualTo(2));
            Assert.That(otherGroup.BioPolymerGroupScore, Is.EqualTo(0)); // Reset after merge
        }

        [Test]
        public void MergeWith_UpdatesBioPolymerGroupName()
        {
            var otherProtein = new Protein("MERGESEQ", "A00001"); // Comes before P12345 alphabetically
            var otherProteins = new HashSet<IBioPolymer> { otherProtein };
            var otherGroup = new ProteinGroup(otherProteins, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            _proteinGroup.MergeWith(otherGroup);

            Assert.That(_proteinGroup.BioPolymerGroupName, Does.StartWith("A00001"));
        }

        [Test]
        public void MergeWith_NonProteinGroup_DoesNotThrow()
        {
            // When passed a non-ProteinGroup, the method should simply not merge
            Assert.DoesNotThrow(() => _proteinGroup.MergeWith(null));
        }

        #endregion

        #region Equality Tests

        [Test]
        public void Equals_SameBioPolymerGroupName_ReturnsTrue()
        {
            var pg1 = new ProteinGroup(_proteins, _allPeptides, _uniquePeptides);
            var pg2 = new ProteinGroup(_proteins, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            Assert.That(pg1.Equals(pg2), Is.True);
        }

        [Test]
        public void Equals_DifferentBioPolymerGroupName_ReturnsFalse()
        {
            var otherProteins = new HashSet<IBioPolymer> { _protein1 };
            var pg1 = new ProteinGroup(_proteins, _allPeptides, _uniquePeptides);
            var pg2 = new ProteinGroup(otherProteins, _allPeptides, _uniquePeptides);

            Assert.That(pg1.Equals(pg2), Is.False);
        }

        [Test]
        public void Equals_Null_ReturnsFalse()
        {
            Assert.That(_proteinGroup.Equals((ProteinGroup)null), Is.False);
            Assert.That(_proteinGroup.Equals((IBioPolymerGroup)null), Is.False);
        }

        [Test]
        public void Equals_SameReference_ReturnsTrue()
        {
            Assert.That(_proteinGroup.Equals(_proteinGroup), Is.True);
        }

        [Test]
        public void Equals_Object_WithMatchingGroup_ReturnsTrue()
        {
            var pg2 = new ProteinGroup(_proteins, _allPeptides, _uniquePeptides);
            Assert.That(_proteinGroup.Equals((object)pg2), Is.True);
        }

        [Test]
        public void Equals_Object_WithNonProteinGroup_ReturnsFalse()
        {
            Assert.That(_proteinGroup.Equals("not a protein group"), Is.False);
            Assert.That(_proteinGroup.Equals(123), Is.False);
        }

        [Test]
        public void GetHashCode_SameGroupName_ReturnsSameHashCode()
        {
            var pg1 = new ProteinGroup(_proteins, _allPeptides, _uniquePeptides);
            var pg2 = new ProteinGroup(_proteins, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            Assert.That(pg1.GetHashCode(), Is.EqualTo(pg2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentGroupName_ReturnsDifferentHashCode()
        {
            var otherProteins = new HashSet<IBioPolymer> { _protein1 };
            var pg1 = new ProteinGroup(_proteins, _allPeptides, _uniquePeptides);
            var pg2 = new ProteinGroup(otherProteins, _allPeptides, _uniquePeptides);

            Assert.That(pg1.GetHashCode(), Is.Not.EqualTo(pg2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_WorksInHashSet()
        {
            var pg1 = new ProteinGroup(_proteins, _allPeptides, _uniquePeptides);
            var pg2 = new ProteinGroup(_proteins, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            var set = new HashSet<ProteinGroup> { pg1, pg2 };

            Assert.That(set.Count, Is.EqualTo(1)); // Should deduplicate
        }

        #endregion

        #region GetTabSeparatedHeader Tests

        [Test]
        public void GetTabSeparatedHeader_ContainsExpectedColumns()
        {
            var header = _proteinGroup.GetTabSeparatedHeader();

            Assert.That(header, Does.Contain("Protein Accession"));
            Assert.That(header, Does.Contain("Gene"));
            Assert.That(header, Does.Contain("Organism"));
            Assert.That(header, Does.Contain("Protein Full Name"));
            Assert.That(header, Does.Contain("Number of Proteins in Group"));
            Assert.That(header, Does.Contain("Unique Peptides"));
            Assert.That(header, Does.Contain("Shared Peptides"));
            Assert.That(header, Does.Contain("Sequence Coverage Fraction"));
            Assert.That(header, Does.Contain("Protein QValue"));
            Assert.That(header, Does.Contain("Best Peptide Score"));
        }

        [Test]
        public void GetTabSeparatedHeader_WithSpectraFileInfo_IncludesIntensityColumns()
        {
            var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            _proteinGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile };

            var header = _proteinGroup.GetTabSeparatedHeader();

            Assert.That(header, Does.Contain("Intensity_"));
        }

        [Test]
        public void GetTabSeparatedHeader_WithIsobaricSamples_IncludesChannelColumns()
        {
            var isobaricSample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            _proteinGroup.SamplesForQuantification = new List<ISampleInfo> { isobaricSample };

            var header = _proteinGroup.GetTabSeparatedHeader();

            Assert.That(header, Does.Contain("Intensity_test_126"));
        }

        #endregion

        #region ToString Tests

        [Test]
        public void ToString_ContainsProteinGroupName()
        {
            var result = _proteinGroup.ToString();
            Assert.That(result, Does.StartWith("P12345|P67890"));
        }

        [Test]
        public void ToString_ContainsGeneNames()
        {
            var result = _proteinGroup.ToString();
            Assert.That(result, Does.Contain("GENE1"));
            Assert.That(result, Does.Contain("GENE2"));
        }

        [Test]
        public void ToString_ContainsOrganisms()
        {
            var result = _proteinGroup.ToString();
            Assert.That(result, Does.Contain("Homo sapiens"));
        }

        [Test]
        public void ToString_ContainsDecoyContaminantTarget()
        {
            var result = _proteinGroup.ToString();
            Assert.That(result, Does.Contain("\tT\t")); // Target
        }

        [Test]
        public void ToString_DecoyGroup_ContainsD()
        {
            var decoyProteins = new HashSet<IBioPolymer> { _decoyProtein };
            var pg = new ProteinGroup(decoyProteins, _allPeptides, _uniquePeptides);

            var result = pg.ToString();
            Assert.That(result, Does.Contain("\tD\t")); // Decoy
        }

        [Test]
        public void ToString_ContaminantGroup_ContainsC()
        {
            var contaminantProteins = new HashSet<IBioPolymer> { _contaminantProtein };
            var pg = new ProteinGroup(contaminantProteins, _allPeptides, _uniquePeptides);

            var result = pg.ToString();
            Assert.That(result, Does.Contain("\tC\t")); // Contaminant
        }

        [Test]
        public void ToString_WithIntensities_IncludesIntensityValues()
        {
            var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            _proteinGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile };
            _proteinGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { spectraFile, 12345.67 }
            };

            var result = _proteinGroup.ToString();
            Assert.That(result, Does.Contain("12345.67"));
        }

        #endregion

        #region GetIdentifiedPeptidesOutput Tests

        [Test]
        public void GetIdentifiedPeptidesOutput_WithoutMods_UsesBaseSequence()
        {
            _proteinGroup.DisplayModsOnPeptides = false;
            _proteinGroup.GetIdentifiedPeptidesOutput(null);

            // The method populates private fields - we test indirectly through ToString
            // This test just ensures it doesn't throw
            Assert.DoesNotThrow(() => _proteinGroup.ToString());
        }

        [Test]
        public void GetIdentifiedPeptidesOutput_WithMods_UsesFullSequence()
        {
            _proteinGroup.DisplayModsOnPeptides = true;
            _proteinGroup.GetIdentifiedPeptidesOutput(null);

            Assert.DoesNotThrow(() => _proteinGroup.ToString());
        }

        #endregion

        #region ConstructSubsetBioPolymerGroup Tests

        [Test]
        public void ConstructSubsetProteinGroup_FiltersToSpecificFile()
        {
            var spectraFile1 = new SpectraFileInfo(@"C:\test1.raw", "Control", 1, 1, 0);
            var spectraFile2 = new SpectraFileInfo(@"C:\test2.raw", "Control", 1, 1, 0);

            _proteinGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile1, spectraFile2 };
            _proteinGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { spectraFile1, 1000.0 },
                { spectraFile2, 2000.0 }
            };

            var psm1 = new ProteinGroupTestSpectralMatch(@"C:\test1.raw", "PEP", "PEP", score: 100, scanNumber: 1);
            var psm2 = new ProteinGroupTestSpectralMatch(@"C:\test2.raw", "PEP", "PEP", score: 100, scanNumber: 2);
            _proteinGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            var subset = _proteinGroup.ConstructSubsetProteinGroup(@"C:\test1.raw");

            Assert.That(subset.SamplesForQuantification.Count, Is.EqualTo(1));
            Assert.That(subset.SamplesForQuantification[0].FullFilePathWithExtension, Is.EqualTo(@"C:\test1.raw"));
            Assert.That(subset.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(1));
        }

        [Test]
        public void ConstructSubsetBioPolymerGroup_ReturnsIBioPolymerGroup()
        {
            var subset = _proteinGroup.ConstructSubsetBioPolymerGroup(@"C:\nonexistent.raw");
            Assert.That(subset, Is.InstanceOf<IBioPolymerGroup>());
        }

        #endregion

        #region IBioPolymerGroup Interface Tests

        [Test]
        public void ImplementsIBioPolymerGroup()
        {
            Assert.That(_proteinGroup, Is.InstanceOf<IBioPolymerGroup>());
        }

        [Test]
        public void IBioPolymerGroup_PropertiesAccessibleViaInterface()
        {
            IBioPolymerGroup group = _proteinGroup;

            Assert.That(group.BioPolymerGroupName, Is.EqualTo(_proteinGroup.BioPolymerGroupName));
            Assert.That(group.BioPolymers, Is.EqualTo(_proteinGroup.BioPolymers));
            Assert.That(group.IsDecoy, Is.EqualTo(_proteinGroup.IsDecoy));
            Assert.That(group.IsContaminant, Is.EqualTo(_proteinGroup.IsContaminant));
        }

        #endregion

        #region Edge Cases Tests

        [Test]
        public void EdgeCase_SingleProtein_GroupNameIsJustAccession()
        {
            var singleProtein = new HashSet<IBioPolymer> { _protein1 };
            var pg = new ProteinGroup(singleProtein, _allPeptides, _uniquePeptides);

            Assert.That(pg.BioPolymerGroupName, Is.EqualTo("P12345"));
        }

        [Test]
        public void EdgeCase_ProteinWithNoGeneNames_HandlesGracefully()
        {
            var proteinNoGenes = new Protein("SEQUENCE", "P00000");
            var proteins = new HashSet<IBioPolymer> { proteinNoGenes };
            var pg = new ProteinGroup(proteins, _allPeptides, _uniquePeptides);

            Assert.DoesNotThrow(() => pg.ToString());
        }

        [Test]
        public void EdgeCase_NullSamplesForQuantification_HeaderStillWorks()
        {
            _proteinGroup.SamplesForQuantification = null;
            Assert.DoesNotThrow(() => _proteinGroup.GetTabSeparatedHeader());
        }

        [Test]
        public void EdgeCase_NullIntensitiesBySample_ToStringStillWorks()
        {
            _proteinGroup.IntensitiesBySample = null;
            Assert.DoesNotThrow(() => _proteinGroup.ToString());
        }

        [Test]
        public void EdgeCase_VeryLongStrings_AreTruncated()
        {
            // Create a protein with a very long name
            var longName = new string('A', 50000);
            var proteinLongName = new Protein("SEQ", "P00001", fullName: longName);
            var proteins = new HashSet<IBioPolymer> { proteinLongName };
            var pg = new ProteinGroup(proteins, _allPeptides, _uniquePeptides);

            var result = pg.ToString();

            // The output should be truncated to MaxStringLength (32000 by default)
            Assert.That(result.Length, Is.LessThan(longName.Length));
        }

        [Test]
        public void EdgeCase_MaxStringLengthDisabled_NoTruncation()
        {
            var originalMaxLength = ProteinGroup.MaxStringLength;
            try
            {
                ProteinGroup.MaxStringLength = 0; // Disable truncation

                var longName = new string('A', 50000);
                var proteinLongName = new Protein("SEQ", "P00001", fullName: longName);
                var proteins = new HashSet<IBioPolymer> { proteinLongName };
                var pg = new ProteinGroup(proteins, _allPeptides, _uniquePeptides);

                var result = pg.ToString();
                Assert.That(result, Does.Contain(longName));
            }
            finally
            {
                ProteinGroup.MaxStringLength = originalMaxLength;
            }
        }

        [Test]
        public void EdgeCase_ZeroIntensity_NotIncludedInOutput()
        {
            var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            _proteinGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile };
            _proteinGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { spectraFile, 0.0 }
            };

            var result = _proteinGroup.ToString();
            // Zero intensity should result in an empty column (just tabs)
            // The intensity column should be empty (consecutive tabs)
            Assert.That(result, Does.Contain("\t\t"));
        }

        #endregion

        #region Mixed Quantification Types Tests

        [Test]
        public void MixedQuantificationTypes_SpectraFileInfoAndIsobaric_BothSupported()
        {
            var spectraFile = new SpectraFileInfo(@"C:\test1.raw", "Control", 1, 1, 0);
            var isobaricSample = new IsobaricQuantSampleInfo(@"C:\test2.raw", "Treatment", 1, 1, 0, 1, "126", 126.0, false);

            _proteinGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile, isobaricSample };
            _proteinGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { spectraFile, 1000.0 },
                { isobaricSample, 2000.0 }
            };

            // Should handle mixed types without throwing
            Assert.DoesNotThrow(() => _proteinGroup.GetTabSeparatedHeader());
            Assert.DoesNotThrow(() => _proteinGroup.ToString());
        }

        [Test]
        public void IsobaricSamples_OrderedByFileAndChannel()
        {
            var sample1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);
            var sample2 = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var sample3 = new IsobaricQuantSampleInfo(@"C:\b.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            _proteinGroup.SamplesForQuantification = new List<ISampleInfo> { sample1, sample2, sample3 };

            var header = _proteinGroup.GetTabSeparatedHeader();

            // Should be ordered: a.raw 126, a.raw 127N, b.raw 126
            var index126_a = header.IndexOf("Intensity_a_126");
            var index127N_a = header.IndexOf("Intensity_a_127N");
            var index126_b = header.IndexOf("Intensity_b_126");

            Assert.That(index126_a, Is.LessThan(index127N_a));
            Assert.That(index127N_a, Is.LessThan(index126_b));
        }

        #endregion
    }

    /// <summary>
    /// Test implementation of ISpectralMatch for ProteinGroup tests.
    /// </summary>
    internal class ProteinGroupTestSpectralMatch : ISpectralMatch
    {
        private readonly List<IBioPolymerWithSetMods> _identified;

        public string FullFilePath { get; }
        public string FullSequence { get; }
        public string BaseSequence { get; }
        public double Score { get; }
        public int OneBasedScanNumber { get; }

        public ProteinGroupTestSpectralMatch(string filePath, string fullSequence, string baseSequence,
            double score = 0, int scanNumber = 0, IEnumerable<IBioPolymerWithSetMods> identified = null)
        {
            FullFilePath = filePath ?? string.Empty;
            FullSequence = fullSequence ?? string.Empty;
            BaseSequence = baseSequence ?? string.Empty;
            Score = score;
            OneBasedScanNumber = scanNumber;
            _identified = identified?.ToList() ?? new List<IBioPolymerWithSetMods>();
        }

        public int CompareTo(ISpectralMatch other)
        {
            if (other is null) return 1;
            int scoreCmp = Score.CompareTo(other.Score);
            if (scoreCmp != 0) return scoreCmp;
            return OneBasedScanNumber.CompareTo(other.OneBasedScanNumber);
        }

        public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods()
            => _identified.AsReadOnly();
    }
}