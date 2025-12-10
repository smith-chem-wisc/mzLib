using Chemistry;
using Easy.Common.Extensions;
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
using System.IO;
using System.Linq;

namespace Test.Omics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class IBioPolymerGroupTests
    {
        private TestBioPolymer _protein1;
        private TestBioPolymer _protein2;
        private TestBioPolymer _decoyProtein;
        private TestBioPolymer _contaminantProtein;
        private TestBioPolymerWithSetMods _peptide1;
        private TestBioPolymerWithSetMods _peptide2;
        private TestBioPolymerWithSetMods _uniquePeptide;
        private SpectraFileInfo _file1;
        private SpectraFileInfo _file2;

        [SetUp]
        public void Setup()
        {
            _protein1 = new TestBioPolymer
            {
                Accession = "P12345",
                BaseSequence = "ACDEFGHIK",
                IsDecoy = false,
                IsContaminant = false,
                Organism = "Homo sapiens",
                FullName = "Test Protein 1",
                GeneNames = new List<Tuple<string, string>> { new("primary", "GENE1") }
            };

            _protein2 = new TestBioPolymer
            {
                Accession = "P67890",
                BaseSequence = "KLMNPQRST",
                IsDecoy = false,
                IsContaminant = false,
                Organism = "Homo sapiens",
                FullName = "Test Protein 2",
                GeneNames = new List<Tuple<string, string>> { new("primary", "GENE2") }
            };

            _decoyProtein = new TestBioPolymer
            {
                Accession = "DECOY_P12345",
                BaseSequence = "KIHGFEDCA",
                IsDecoy = true,
                IsContaminant = false,
                Organism = "Homo sapiens",
                FullName = "Decoy Protein"
            };

            _contaminantProtein = new TestBioPolymer
            {
                Accession = "CONT_P99999",
                BaseSequence = "CONTAMINANT",
                IsDecoy = false,
                IsContaminant = true,
                Organism = "Contaminant",
                FullName = "Contaminant Protein"
            };

            _peptide1 = new TestBioPolymerWithSetMods
            {
                BaseSequence = "ACDEFGHIK",
                FullSequence = "ACDEFGHIK",
                MostAbundantMonoisotopicMass = 1000.0
            };

            _peptide2 = new TestBioPolymerWithSetMods
            {
                BaseSequence = "KLMNPQR",
                FullSequence = "KLMNPQR",
                MostAbundantMonoisotopicMass = 800.0
            };

            _uniquePeptide = new TestBioPolymerWithSetMods
            {
                BaseSequence = "UNIQUE",
                FullSequence = "UNIQUE",
                MostAbundantMonoisotopicMass = 600.0
            };

            _file1 = new SpectraFileInfo("path/to/file1.mzML", "Condition1", 1, 1, 1);
            _file2 = new SpectraFileInfo("path/to/file2.mzML", "Condition2", 1, 1, 1);
        }

        #region Constructor Tests

        [Test]
        public void Constructor_WithValidParameters_CreatesGroup()
        {
            var bioPolymers = new HashSet<IBioPolymer> { _protein1 };
            var allBpwsm = new HashSet<IBioPolymerWithSetMods> { _peptide1 };
            var uniqueBpwsm = new HashSet<IBioPolymerWithSetMods> { _uniquePeptide };

            var group = new TestBioPolymerGroup<SpectraFileInfo>(bioPolymers, allBpwsm, uniqueBpwsm);

            Assert.That(group.BioPolymers, Is.EqualTo(bioPolymers));
            Assert.That(group.AllBioPolymerWithSetMods, Is.EqualTo(allBpwsm));
            Assert.That(group.UniqueBioPolymerWithSetMods, Is.EqualTo(uniqueBpwsm));
            Assert.That(group.IsDecoy, Is.False);
            Assert.That(group.IsContaminant, Is.False);
        }

        [Test]
        public void Constructor_WithNullBioPolymers_ThrowsArgumentNullException()
        {
            Assert.Throws<ArgumentNullException>(() =>
                new TestBioPolymerGroup<SpectraFileInfo>(null!, new HashSet<IBioPolymerWithSetMods>(),
                    new HashSet<IBioPolymerWithSetMods>()));
        }

        [Test]
        public void Constructor_WithNullAllBpwsm_CreatesEmptySet()
        {
            var bioPolymers = new HashSet<IBioPolymer> { _protein1 };

            var group = new TestBioPolymerGroup<SpectraFileInfo>(bioPolymers, null!, null!);

            Assert.That(group.AllBioPolymerWithSetMods, Is.Not.Null);
            Assert.That(group.AllBioPolymerWithSetMods, Is.Empty);
            Assert.That(group.UniqueBioPolymerWithSetMods, Is.Not.Null);
            Assert.That(group.UniqueBioPolymerWithSetMods, Is.Empty);
        }

        [Test]
        public void Constructor_WithDecoyProteins_SetsIsDecoyTrue()
        {
            var bioPolymers = new HashSet<IBioPolymer> { _decoyProtein };

            var group = new TestBioPolymerGroup<SpectraFileInfo>(bioPolymers,
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group.IsDecoy, Is.True);
        }

        [Test]
        public void Constructor_WithMixedDecoyAndTarget_SetsIsDecoyFalse()
        {
            var bioPolymers = new HashSet<IBioPolymer> { _protein1, _decoyProtein };

            var group = new TestBioPolymerGroup<SpectraFileInfo>(bioPolymers,
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group.IsDecoy, Is.False);
        }

        [Test]
        public void Constructor_WithContaminantProtein_SetsIsContaminantTrue()
        {
            var bioPolymers = new HashSet<IBioPolymer> { _contaminantProtein };

            var group = new TestBioPolymerGroup<SpectraFileInfo>(bioPolymers,
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group.IsContaminant, Is.True);
        }

        [Test]
        public void Constructor_InitializesDefaultValues()
        {
            var bioPolymers = new HashSet<IBioPolymer> { _protein1 };

            var group = new TestBioPolymerGroup<SpectraFileInfo>(bioPolymers,
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group.QValue, Is.EqualTo(double.MaxValue));
            Assert.That(group.BestBiopolymerWithSetsModQValue, Is.EqualTo(double.MaxValue));
            Assert.That(group.BestBioPolymerWithSetsModScore, Is.EqualTo(0));
            Assert.That(group.BioPolymerGroupScore, Is.EqualTo(0));
            Assert.That(group.FilesForQuantification, Is.Empty);
            Assert.That(group.IntensitiesByFile, Is.Empty);
            Assert.That(group.ModsInfo, Is.Empty);
            Assert.That(group.AllPsmsBelowOnePercentFDR, Is.Empty);
        }

        #endregion

        #region BioPolymerGroupName Tests

        [Test]
        public void BioPolymerGroupName_SingleProtein_ReturnsAccession()
        {
            var bioPolymers = new HashSet<IBioPolymer> { _protein1 };
            var group = new TestBioPolymerGroup<SpectraFileInfo>(bioPolymers,
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group.BioPolymerGroupName, Is.EqualTo("P12345"));
        }

        [Test]
        public void BioPolymerGroupName_MultipleProteins_ReturnsPipeDelimitedOrderedAccessions()
        {
            var bioPolymers = new HashSet<IBioPolymer> { _protein2, _protein1 }; // Add in reverse order
            var group = new TestBioPolymerGroup<SpectraFileInfo>(bioPolymers,
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            // Should be ordered alphabetically
            Assert.That(group.BioPolymerGroupName, Is.EqualTo("P12345|P67890"));
        }

        [Test]
        public void ListOfBioPolymersOrderedByAccession_IsCached()
        {
            var bioPolymers = new HashSet<IBioPolymer> { _protein2, _protein1 };
            var group = new TestBioPolymerGroup<SpectraFileInfo>(bioPolymers,
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            var list1 = group.ListOfBioPolymersOrderedByAccession;
            var list2 = group.ListOfBioPolymersOrderedByAccession;

            Assert.That(list1, Is.SameAs(list2)); // Same reference = cached
        }

        #endregion

        #region Score Tests

        [Test]
        public void Score_WithUniquePeptides_CalculatesCorrectScore()
        {
            var bioPolymers = new HashSet<IBioPolymer> { _protein1 };
            var allBpwsm = new HashSet<IBioPolymerWithSetMods> { _peptide1, _peptide2 };
            var uniqueBpwsm = new HashSet<IBioPolymerWithSetMods> { _uniquePeptide };

            var group = new TestBioPolymerGroup<SpectraFileInfo>(bioPolymers, allBpwsm, uniqueBpwsm);
            group.Score();

            // Score = uniqueCount + (allCount - uniqueCount) * 0.5 = 1 + (2 - 1) * 0.5 = 1.5
            Assert.That(group.BioPolymerGroupScore, Is.EqualTo(1.5));
        }

        [Test]
        public void Score_WithNoUniquePeptides_CalculatesCorrectScore()
        {
            var bioPolymers = new HashSet<IBioPolymer> { _protein1 };
            var allBpwsm = new HashSet<IBioPolymerWithSetMods> { _peptide1, _peptide2 };
            var uniqueBpwsm = new HashSet<IBioPolymerWithSetMods>();

            var group = new TestBioPolymerGroup<SpectraFileInfo>(bioPolymers, allBpwsm, uniqueBpwsm);
            group.Score();

            // Score = 0 + 2 * 0.5 = 1.0
            Assert.That(group.BioPolymerGroupScore, Is.EqualTo(1.0));
        }

        #endregion

        #region MergeWith Tests

        [Test]
        public void MergeWith_NullGroup_ThrowsArgumentNullException()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            Assert.Throws<ArgumentNullException>(() => group.MergeWith(null!));
        }

        [Test]
        public void MergeWith_MergesBioPolymers()
        {
            var group1 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods> { _peptide1 },
                new HashSet<IBioPolymerWithSetMods> { _uniquePeptide });

            var group2 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein2 },
                new HashSet<IBioPolymerWithSetMods> { _peptide2 },
                new HashSet<IBioPolymerWithSetMods>());

            group1.MergeWith(group2);

            Assert.That(group1.BioPolymers, Has.Count.EqualTo(2));
            Assert.That(group1.BioPolymers, Contains.Item(_protein1));
            Assert.That(group1.BioPolymers, Contains.Item(_protein2));
        }

        [Test]
        public void MergeWith_MergesPeptides()
        {
            var group1 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods> { _peptide1 },
                new HashSet<IBioPolymerWithSetMods> { _uniquePeptide });

            var group2 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein2 },
                new HashSet<IBioPolymerWithSetMods> { _peptide2 },
                new HashSet<IBioPolymerWithSetMods>());

            group1.MergeWith(group2);

            Assert.That(group1.AllBioPolymerWithSetMods, Has.Count.EqualTo(2));
        }

        [Test]
        public void MergeWith_SumsIntensities()
        {
            var group1 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());
            group1.IntensitiesByFile[_file1] = 1000.0;
            group1.FilesForQuantification.Add(_file1);

            var group2 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein2 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());
            group2.IntensitiesByFile[_file1] = 500.0;
            group2.FilesForQuantification.Add(_file1);

            group1.MergeWith(group2);

            Assert.That(group1.IntensitiesByFile[_file1], Is.EqualTo(1500.0));
        }

        [Test]
        public void MergeWith_UpdatesBestScores()
        {
            var group1 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>())
            {
                BestBioPolymerWithSetsModScore = 50.0,
                BestBiopolymerWithSetsModQValue = 0.05
            };

            var group2 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein2 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>())
            {
                BestBioPolymerWithSetsModScore = 100.0,
                BestBiopolymerWithSetsModQValue = 0.01
            };

            group1.MergeWith(group2);

            Assert.That(group1.BestBioPolymerWithSetsModScore, Is.EqualTo(100.0));
            Assert.That(group1.BestBiopolymerWithSetsModQValue, Is.EqualTo(0.01));
        }

        [Test]
        public void MergeWith_InvalidatesCachedList()
        {
            var group1 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein2 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            // Access to cache the list
            var originalName = group1.BioPolymerGroupName;
            Assert.That(originalName, Is.EqualTo("P67890"));

            var group2 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            group1.MergeWith(group2);

            // Name should now include both proteins
            Assert.That(group1.BioPolymerGroupName, Is.EqualTo("P12345|P67890"));
        }

        #endregion

        #region Equality Tests

        [Test]
        public void Equals_SameGroupName_ReturnsTrue()
        {
            var group1 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            var group2 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods> { _peptide1 }, // Different peptides
                new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group1.Equals(group2), Is.True);
        }

        [Test]
        public void Equals_DifferentGroupName_ReturnsFalse()
        {
            var group1 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            var group2 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein2 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group1.Equals(group2), Is.False);
        }

        [Test]
        public void Equals_Null_ReturnsFalse()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group.Equals(null), Is.False);
        }

        [Test]
        public void Equals_SameReference_ReturnsTrue()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group.Equals(group), Is.True);
        }

        [Test]
        public void GetHashCode_SameGroupName_ReturnsSameHashCode()
        {
            var group1 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            var group2 = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group1.GetHashCode(), Is.EqualTo(group2.GetHashCode()));
        }

        #endregion

        #region GetTabSeparatedHeader Tests

        [Test]
        public void GetTabSeparatedHeader_ContainsRequiredColumns()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            var header = group.GetTabSeparatedHeader();

            Assert.That(header, Does.Contain("Accession"));
            Assert.That(header, Does.Contain("Gene"));
            Assert.That(header, Does.Contain("Organism"));
            Assert.That(header, Does.Contain("IsDecoy"));
            Assert.That(header, Does.Contain("IsContaminant"));
            Assert.That(header, Does.Contain("Score"));
            Assert.That(header, Does.Contain("QValue"));
        }

        [Test]
        public void GetTabSeparatedHeader_IncludesIntensityColumns()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());
            group.FilesForQuantification.Add(_file1);
            group.FilesForQuantification.Add(_file2);

            var header = group.GetTabSeparatedHeader();

            Assert.That(header, Does.Contain("Intensity_"));
        }

        #endregion

        #region ToString Tests

        [Test]
        public void ToString_ContainsGroupName()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            var result = group.ToString();

            Assert.That(result, Does.Contain("P12345"));
        }

        [Test]
        public void ToString_ContainsOrganism()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            var result = group.ToString();

            Assert.That(result, Does.Contain("Homo sapiens"));
        }

        [Test]
        public void ToString_IncludesIntensities()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());
            group.FilesForQuantification.Add(_file1);
            group.IntensitiesByFile[_file1] = 12345.67;

            var result = group.ToString();

            Assert.That(result, Does.Contain("12346")); // Formatted as F0
        }

        #endregion

        #region ConstructSubsetBioPolymerGroup Tests

        [Test]
        public void ConstructSubsetBioPolymerGroup_NullFilePath_ThrowsArgumentNullException()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            Assert.Throws<ArgumentNullException>(() => group.ConstructSubsetBioPolymerGroup(null!));
        }

        [Test]
        public void ConstructSubsetBioPolymerGroup_EmptyFilePath_ThrowsArgumentNullException()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            Assert.Throws<ArgumentNullException>(() => group.ConstructSubsetBioPolymerGroup(string.Empty));
        }

        [Test]
        public void ConstructSubsetBioPolymerGroup_PreservesBioPolymers()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1, _protein2 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            var subset = group.ConstructSubsetBioPolymerGroup("path/to/file.mzML");

            Assert.That(subset.BioPolymers, Has.Count.EqualTo(2));
        }

        [Test]
        public void ConstructSubsetBioPolymerGroup_PreservesQValues()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>())
            {
                QValue = 0.01,
                BestBiopolymerWithSetsModQValue = 0.005,
                BestBioPolymerWithSetsModScore = 100.0
            };

            var subset = group.ConstructSubsetBioPolymerGroup("path/to/file.mzML");

            Assert.That(subset.QValue, Is.EqualTo(0.01));
            Assert.That(subset.BestBiopolymerWithSetsModQValue, Is.EqualTo(0.005));
            Assert.That(subset.BestBioPolymerWithSetsModScore, Is.EqualTo(100.0));
        }

        #endregion

        #region Edge Cases

        [Test]
        public void EmptyBioPolymers_GroupNameIsEmpty()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer>(),
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group.BioPolymerGroupName, Is.EqualTo(string.Empty));
        }

        [Test]
        public void EmptyBioPolymers_IsDecoyFalse()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer>(),
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            Assert.That(group.IsDecoy, Is.False);
        }

        [Test]
        public void IntensitiesByFile_MissingFile_ReturnsZeroInToString()
        {
            var group = new TestBioPolymerGroup<SpectraFileInfo>(
                new HashSet<IBioPolymer> { _protein1 },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());
            group.FilesForQuantification.Add(_file1);
            // Don't add intensity for _file1

            var result = group.ToString();

            Assert.That(result, Does.Contain("\t0\t").Or.EndWith("\t0"));
        }

        #endregion
    }

    #region Test Helper Classes

    /// <summary>
    /// Test implementation of IBioPolymerGroup for unit testing purposes.
    /// </summary>
    internal class TestBioPolymerGroup<TSampleInfo> : IBioPolymerGroup<TSampleInfo>
        where TSampleInfo : notnull
    {
        private List<IBioPolymer>? _listOfBioPolymersOrderedByAccession;

        public TestBioPolymerGroup(
            HashSet<IBioPolymer> bioPolymers,
            HashSet<IBioPolymerWithSetMods> allBioPolymerWithSetMods,
            HashSet<IBioPolymerWithSetMods> uniqueBioPolymerWithSetMods)
        {
            BioPolymers = bioPolymers ?? throw new ArgumentNullException(nameof(bioPolymers));
            AllBioPolymerWithSetMods = allBioPolymerWithSetMods ?? new HashSet<IBioPolymerWithSetMods>();
            UniqueBioPolymerWithSetMods = uniqueBioPolymerWithSetMods ?? new HashSet<IBioPolymerWithSetMods>();
            AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>();
            FilesForQuantification = new List<TSampleInfo>();
            IntensitiesByFile = new Dictionary<TSampleInfo, double>();
            QValue = double.MaxValue;
            BestBiopolymerWithSetsModQValue = double.MaxValue;
            BestBioPolymerWithSetsModScore = 0;
            BioPolymerGroupScore = 0;
        }

        public bool IsDecoy => BioPolymers.Count > 0 && BioPolymers.All(p => p.IsDecoy);
        public bool IsContaminant => BioPolymers.Any(p => p.IsContaminant);
        public List<TSampleInfo> FilesForQuantification { get; set; }
        public HashSet<IBioPolymer> BioPolymers { get; set; }
        public string BioPolymerGroupName => string.Join("|", ListOfBioPolymersOrderedByAccession.Select(p => p.Accession));
        public double BioPolymerGroupScore { get; set; }
        public HashSet<IBioPolymerWithSetMods> AllBioPolymerWithSetMods { get; set; }
        public HashSet<IBioPolymerWithSetMods> UniqueBioPolymerWithSetMods { get; set; }
        public HashSet<ISpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }
        public double QValue { get; set; }
        public double BestBiopolymerWithSetsModQValue { get; set; }
        public double BestBioPolymerWithSetsModScore { get; set; }

        public List<string> ModsInfo
        {
            get
            {
                var modsInfo = new List<string>();
                var allMods = AllBioPolymerWithSetMods
                    .SelectMany(p => p.AllModsOneIsNterminus.Values)
                    .GroupBy(m => m.IdWithMotif)
                    .OrderByDescending(g => g.Count());

                foreach (var modGroup in allMods)
                {
                    modsInfo.Add($"{modGroup.Key}({modGroup.Count()})");
                }

                return modsInfo;
            }
        }

        public Dictionary<TSampleInfo, double> IntensitiesByFile { get; set; }

        public List<IBioPolymer> ListOfBioPolymersOrderedByAccession
        {
            get
            {
                _listOfBioPolymersOrderedByAccession ??= BioPolymers
                    .OrderBy(p => p.Accession, StringComparer.Ordinal)
                    .ToList();
                return _listOfBioPolymersOrderedByAccession;
            }
        }

        public string GetTabSeparatedHeader()
        {
            var headers = new List<string>
            {
                "Accession", "Gene", "Organism", "IsDecoy", "IsContaminant",
                "UniqueSequences", "SharedSequences", "PSMs", "Score", "QValue",
                "BestPeptideScore", "BestPeptideQValue", "Modifications"
            };

            foreach (var file in FilesForQuantification)
            {
                headers.Add($"Intensity_{file}");
            }

            return string.Join("\t", headers);
        }

        public override string ToString()
        {
            var values = new List<string>
            {
                BioPolymerGroupName,
                string.Join(";", BioPolymers.SelectMany(p => p.GeneNames).Select(g => g.Item2).Distinct()),
                string.Join(";", BioPolymers.Select(p => p.Organism).Distinct()),
                IsDecoy.ToString(),
                IsContaminant.ToString(),
                UniqueBioPolymerWithSetMods.Count.ToString(),
                (AllBioPolymerWithSetMods.Count - UniqueBioPolymerWithSetMods.Count).ToString(),
                AllPsmsBelowOnePercentFDR.Count.ToString(),
                BioPolymerGroupScore.ToString("F4"),
                QValue.ToString("F4"),
                BestBioPolymerWithSetsModScore.ToString("F4"),
                BestBiopolymerWithSetsModQValue.ToString("F4"),
                string.Join(";", ModsInfo)
            };

            foreach (var file in FilesForQuantification)
            {
                values.Add(IntensitiesByFile.TryGetValue(file, out var intensity) ? intensity.ToString("F0") : "0");
            }

            return string.Join("\t", values);
        }

        public void Score()
        {
            BioPolymerGroupScore = UniqueBioPolymerWithSetMods.Count +
                (AllBioPolymerWithSetMods.Count - UniqueBioPolymerWithSetMods.Count) * 0.5;
        }

        public void MergeWith(IBioPolymerGroup<TSampleInfo> otherBioPolymerGroup)
        {
            if (otherBioPolymerGroup == null)
                throw new ArgumentNullException(nameof(otherBioPolymerGroup));

            foreach (var bp in otherBioPolymerGroup.BioPolymers) BioPolymers.Add(bp);
            foreach (var bpwsm in otherBioPolymerGroup.AllBioPolymerWithSetMods) AllBioPolymerWithSetMods.Add(bpwsm);
            foreach (var bpwsm in otherBioPolymerGroup.UniqueBioPolymerWithSetMods) UniqueBioPolymerWithSetMods.Add(bpwsm);
            foreach (var psm in otherBioPolymerGroup.AllPsmsBelowOnePercentFDR) AllPsmsBelowOnePercentFDR.Add(psm);

            foreach (var kvp in otherBioPolymerGroup.IntensitiesByFile)
            {
                IntensitiesByFile[kvp.Key] = IntensitiesByFile.TryGetValue(kvp.Key, out var existing)
                    ? existing + kvp.Value : kvp.Value;
            }

            foreach (var file in otherBioPolymerGroup.FilesForQuantification)
            {
                if (!FilesForQuantification.Contains(file)) FilesForQuantification.Add(file);
            }

            if (otherBioPolymerGroup.BestBioPolymerWithSetsModScore > BestBioPolymerWithSetsModScore)
                BestBioPolymerWithSetsModScore = otherBioPolymerGroup.BestBioPolymerWithSetsModScore;

            if (otherBioPolymerGroup.BestBiopolymerWithSetsModQValue < BestBiopolymerWithSetsModQValue)
                BestBiopolymerWithSetsModQValue = otherBioPolymerGroup.BestBiopolymerWithSetsModQValue;

            _listOfBioPolymersOrderedByAccession = null;
            Score();
        }

        public IBioPolymerGroup<TSampleInfo> ConstructSubsetBioPolymerGroup(string fullFilePath, List<SilacLabel>? silacLabels = null)
        {
            if (string.IsNullOrEmpty(fullFilePath))
                throw new ArgumentNullException(nameof(fullFilePath));

            var subsetPsms = AllPsmsBelowOnePercentFDR
                .Where(psm => psm.FullFilePath == fullFilePath ||
                              Path.GetFileNameWithoutExtension(psm.FullFilePath) == Path.GetFileNameWithoutExtension(fullFilePath))
                .ToHashSet();

            var subsetSequences = subsetPsms.Select(psm => psm.FullSequence).ToHashSet();

            var subsetBpwsm = AllBioPolymerWithSetMods
                .Where(p => subsetSequences.Contains(p.FullSequence))
                .ToHashSet();

            var subsetUniqueBpwsm = UniqueBioPolymerWithSetMods
                .Where(p => subsetSequences.Contains(p.FullSequence))
                .ToHashSet();

            var subsetGroup = new TestBioPolymerGroup<TSampleInfo>(
                new HashSet<IBioPolymer>(BioPolymers),
                subsetBpwsm,
                subsetUniqueBpwsm)
            {
                QValue = QValue,
                BestBiopolymerWithSetsModQValue = BestBiopolymerWithSetsModQValue,
                BestBioPolymerWithSetsModScore = BestBioPolymerWithSetsModScore,
                AllPsmsBelowOnePercentFDR = subsetPsms
            };

            var matchingFile = FilesForQuantification.FirstOrDefault(f => f?.ToString() == fullFilePath);
            if (matchingFile != null && IntensitiesByFile.TryGetValue(matchingFile, out var intensity))
            {
                subsetGroup.IntensitiesByFile[matchingFile] = intensity;
                subsetGroup.FilesForQuantification.Add(matchingFile);
            }

            subsetGroup.Score();
            return subsetGroup;
        }

        public bool Equals(IBioPolymerGroup<TSampleInfo>? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            return BioPolymerGroupName == other.BioPolymerGroupName;
        }

        public override bool Equals(object? obj) => Equals(obj as IBioPolymerGroup<TSampleInfo>);
        public override int GetHashCode() => StringComparer.Ordinal.GetHashCode(BioPolymerGroupName ?? string.Empty);
    }

    /// <summary>
    /// Minimal test implementation of IBioPolymer for unit tests.
    /// </summary>
    internal class TestBioPolymer : IBioPolymer
    {
        public string Accession { get; init; } = string.Empty;
        public string BaseSequence { get; init; } = string.Empty;
        public bool IsDecoy { get; init; }
        public bool IsContaminant { get; init; }
        public string Organism { get; init; } = string.Empty;
        public string FullName { get; init; } = string.Empty;
        public string Name { get; init; } = string.Empty;
        public string DatabaseFilePath { get; init; } = string.Empty;
        public int Length => BaseSequence?.Length ?? 0;
        public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];
        public List<Tuple<string, string>> GeneNames { get; init; } = new();
        public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; init; } = new Dictionary<int, List<Modification>>();
        public string SampleNameForVariants { get; init; } = string.Empty;
        public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; } = new Dictionary<int, List<Modification>>();
        public IBioPolymer ConsensusVariant => this;
        public List<SequenceVariation> AppliedSequenceVariations { get; init; } = new();
        public List<SequenceVariation> SequenceVariations { get; init; } = new();
        public List<TruncationProduct> TruncationProducts { get; init; } = new();

        public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams, List<Modification> allKnownFixedModifications,
            List<Modification> variableModifications, List<SilacLabel>? silacLabels = null,
            (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null, bool topDownTruncationSearch = false)
            => Enumerable.Empty<IBioPolymerWithSetMods>();

        public IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence, IDictionary<int, List<Modification>>? newMods)
            => new TestBioPolymer
            {
                Accession = Accession,
                BaseSequence = newBaseSequence,
                IsDecoy = IsDecoy,
                IsContaminant = IsContaminant,
                Organism = Organism,
                FullName = FullName,
                Name = Name,
                DatabaseFilePath = DatabaseFilePath,
                GeneNames = GeneNames,
                OneBasedPossibleLocalizedModifications = newMods ?? new Dictionary<int, List<Modification>>()
            };

        public TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original,
            IEnumerable<SequenceVariation> appliedSequenceVariants, IEnumerable<TruncationProduct> applicableProteolysisProducts,
            IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants) where TBioPolymerType : IHasSequenceVariants
            => (TBioPolymerType)(object)new TestBioPolymer
            {
                Accession = Accession,
                BaseSequence = variantBaseSequence,
                GeneNames = GeneNames,
                OneBasedPossibleLocalizedModifications = oneBasedModifications,
                AppliedSequenceVariations = appliedSequenceVariants.ToList()
            };

        public bool Equals(IBioPolymer? other) => other != null && Accession == other.Accession && BaseSequence == other.BaseSequence;
        public override bool Equals(object? obj) => Equals(obj as IBioPolymer);
        public override int GetHashCode() => HashCode.Combine(Accession, BaseSequence);
    }

    /// <summary>
    /// Minimal test implementation of IBioPolymerWithSetMods for unit tests.
    /// </summary>
    internal class TestBioPolymerWithSetMods : IBioPolymerWithSetMods
    {
        public string BaseSequence { get; init; } = string.Empty;
        public string FullSequence { get; init; } = string.Empty;
        public double MostAbundantMonoisotopicMass { get; init; }
        public double MonoisotopicMass => MostAbundantMonoisotopicMass;
        public ChemicalFormula ThisChemicalFormula => new();
        public string SequenceWithChemicalFormulas => FullSequence;
        public int OneBasedStartResidue { get; init; } = 1;
        public int OneBasedEndResidue { get; init; } = 1;
        public int MissedCleavages { get; init; }
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
        public IBioPolymer Parent => throw new NotImplementedException();

        public void Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus,
            List<Product> products, FragmentationParams? fragmentationParams = null) => throw new NotImplementedException();
        public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments,
            List<Product> products, FragmentationParams? fragmentationParams = null) => throw new NotImplementedException();
        public IBioPolymerWithSetMods Localize(int indexOfMass, double massToLocalize) => this;

        public bool Equals(IBioPolymerWithSetMods? other) => other != null && FullSequence == other.FullSequence;
        public override bool Equals(object? obj) => Equals(obj as IBioPolymerWithSetMods);
        public override int GetHashCode() => StringComparer.Ordinal.GetHashCode(FullSequence ?? string.Empty);
    }

    /// <summary>
    /// Minimal test implementation of ISpectralMatch for unit tests.
    /// </summary>
    internal class TestOmicsSpectralMatch : ISpectralMatch
    {
        public string FullFilePath { get; init; } = string.Empty;
        public string FullSequence { get; init; } = string.Empty;
        public string BaseSequence { get; init; } = string.Empty;
        public double Score { get; init; }
        public int OneBasedScanNumber { get; init; }

        public int CompareTo(ISpectralMatch? other)
        {
            if (other is null) return 1;
            int cmp = Score.CompareTo(other.Score);
            return cmp != 0 ? cmp : OneBasedScanNumber.CompareTo(other.OneBasedScanNumber);
        }

        public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => new List<IBioPolymerWithSetMods>();
    }

    #endregion
}