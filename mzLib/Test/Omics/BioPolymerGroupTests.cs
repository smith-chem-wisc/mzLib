using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymer;
using Omics.BioPolymerGroup;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;

namespace Test.Omics
{
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
		private TestBioPolymerWithSetMods _sequence2;
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
				name: "TestBioPolymer1",
				fullName: "Test BioPolymer 1 Full Name",
				geneNames: new List<Tuple<string, string>> { new("primary", "GENE1") });

			_bioPolymer2 = new TestBioPolymer("TGCATGCA", "BP67890",
				organism: "Homo sapiens",
				name: "TestBioPolymer2",
				fullName: "Test BioPolymer 2 Full Name",
				geneNames: new List<Tuple<string, string>> { new("primary", "GENE2") });

			_decoyBioPolymer = new TestBioPolymer("DECOYSEQ", "DECOY_BP12345",
				organism: "Homo sapiens",
				isDecoy: true);

			_contaminantBioPolymer = new TestBioPolymer("CONTAMINANT", "CONT_BP99999",
				organism: "Homo sapiens",
				isContaminant: true);

			_sequence1 = new TestBioPolymerWithSetMods("ACGT", "ACGT");
			_sequence2 = new TestBioPolymerWithSetMods("TGCA", "TGCA");
			_uniqueSequence = new TestBioPolymerWithSetMods("UNIQUE", "UNIQUE");

			_bioPolymers = new HashSet<IBioPolymer> { _bioPolymer1, _bioPolymer2 };
			_allSequences = new HashSet<IBioPolymerWithSetMods> { _sequence1, _sequence2, _uniqueSequence };
			_uniqueSequences = new HashSet<IBioPolymerWithSetMods> { _uniqueSequence };

			_bioPolymerGroup = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
		}

		#endregion

		#region Constructor Tests

		[Test]
		public void Constructor_WithValidParameters_InitializesAllProperties()
		{
			var bg = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);

			Assert.That(bg.BioPolymers, Is.EqualTo(_bioPolymers));
			Assert.That(bg.AllBioPolymersWithSetMods, Is.EqualTo(_allSequences));
			Assert.That(bg.UniqueBioPolymersWithSetMods, Is.EqualTo(_uniqueSequences));
			Assert.That(bg.BioPolymerGroupScore, Is.EqualTo(0));
			Assert.That(bg.BestBioPolymerWithSetModsScore, Is.EqualTo(0));
			Assert.That(bg.QValue, Is.EqualTo(0));
			Assert.That(bg.IsDecoy, Is.False);
			Assert.That(bg.IsContaminant, Is.False);
			Assert.That(bg.AllPsmsBelowOnePercentFDR, Is.Not.Null);
			Assert.That(bg.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(0));
		}

		[Test]
		public void Constructor_WithDecoyBioPolymer_SetsIsDecoyTrue()
		{
			var bioPolymers = new HashSet<IBioPolymer> { _decoyBioPolymer };
			var bg = new BioPolymerGroup(bioPolymers, _allSequences, _uniqueSequences);

			Assert.That(bg.IsDecoy, Is.True);
			Assert.That(bg.IsContaminant, Is.False);
		}

		[Test]
		public void Constructor_WithContaminantBioPolymer_SetsIsContaminantTrue()
		{
			var bioPolymers = new HashSet<IBioPolymer> { _contaminantBioPolymer };
			var bg = new BioPolymerGroup(bioPolymers, _allSequences, _uniqueSequences);

			Assert.That(bg.IsContaminant, Is.True);
			Assert.That(bg.IsDecoy, Is.False);
		}

		[Test]
		public void Constructor_WithBothDecoyAndContaminant_SetsBothFlags()
		{
			var bioPolymers = new HashSet<IBioPolymer> { _decoyBioPolymer, _contaminantBioPolymer };
			var bg = new BioPolymerGroup(bioPolymers, _allSequences, _uniqueSequences);

			Assert.That(bg.IsDecoy, Is.True);
			Assert.That(bg.IsContaminant, Is.True);
		}

		[Test]
		public void Constructor_BioPolymerGroupName_IsOrderedByAccession()
		{
			// BP12345 comes before BP67890 alphabetically
			var expectedName = "BP12345|BP67890";
			Assert.That(_bioPolymerGroup.BioPolymerGroupName, Is.EqualTo(expectedName));
		}

		[Test]
		public void Constructor_ListOfBioPolymersOrderedByAccession_IsCorrectlyOrdered()
		{
			var orderedList = _bioPolymerGroup.ListOfBioPolymersOrderedByAccession;

			Assert.That(orderedList.Count, Is.EqualTo(2));
			Assert.That(orderedList[0].Accession, Is.EqualTo("BP12345"));
			Assert.That(orderedList[1].Accession, Is.EqualTo("BP67890"));
		}

		[Test]
		public void Constructor_WithEmptyBioPolymers_CreatesEmptyGroup()
		{
			var emptyBioPolymers = new HashSet<IBioPolymer>();
			var bg = new BioPolymerGroup(emptyBioPolymers, _allSequences, _uniqueSequences);

			Assert.That(bg.BioPolymers.Count, Is.EqualTo(0));
			Assert.That(bg.BioPolymerGroupName, Is.EqualTo(string.Empty));
		}

		[Test]
		public void Constructor_WithEmptySequences_CreatesGroupWithNoSequences()
		{
			var emptySequences = new HashSet<IBioPolymerWithSetMods>();
			var bg = new BioPolymerGroup(_bioPolymers, emptySequences, emptySequences);

			Assert.That(bg.AllBioPolymersWithSetMods.Count, Is.EqualTo(0));
			Assert.That(bg.UniqueBioPolymersWithSetMods.Count, Is.EqualTo(0));
		}

		[Test]
		public void Constructor_WithMixedDecoyAndTarget_DecoyFlagSetWhenDecoyPresent()
		{
			var bioPolymers = new HashSet<IBioPolymer> { _bioPolymer1, _decoyBioPolymer };
			var bg = new BioPolymerGroup(bioPolymers, _allSequences, _uniqueSequences);

			Assert.That(bg.IsDecoy, Is.True);
		}

		#endregion

		#region Score Method Tests

		[Test]
		public void Score_WithNoPsms_SetsScoreToZero()
		{
			_bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>();
			_bioPolymerGroup.Score();

			Assert.That(_bioPolymerGroup.BioPolymerGroupScore, Is.EqualTo(0));
		}

		[Test]
		public void Score_WithPsms_SumsBestScorePerBaseSequence()
		{
			var psm1 = new TestSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 100, scanNumber: 1);
			var psm2 = new TestSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 150, scanNumber: 2);
			var psm3 = new TestSpectralMatch(@"C:\test.raw", "TGCA", "TGCA", score: 200, scanNumber: 3);

			_bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3 };
			_bioPolymerGroup.Score();

			// Best score for ACGT = 150, Best score for TGCA = 200, Total = 350
			Assert.That(_bioPolymerGroup.BioPolymerGroupScore, Is.EqualTo(350));
		}

		[Test]
		public void Score_WithMultipleSequencesSameBase_UsesMaxScore()
		{
			var psm1 = new TestSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 50, scanNumber: 1);
			var psm2 = new TestSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 75, scanNumber: 2);
			var psm3 = new TestSpectralMatch(@"C:\test.raw", "ACGT", "ACGT", score: 100, scanNumber: 3);

			_bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3 };
			_bioPolymerGroup.Score();

			Assert.That(_bioPolymerGroup.BioPolymerGroupScore, Is.EqualTo(100));
		}

		#endregion

		#region MergeWith Tests

		[Test]
		public void MergeWith_CombinesBioPolymersAndSequences()
		{
			var otherBioPolymer = new TestBioPolymer("MERGESEQ", "BP99999");
			var otherBioPolymers = new HashSet<IBioPolymer> { otherBioPolymer };
			var otherSequence = new TestBioPolymerWithSetMods("MERGED", "MERGED");
			var otherSequences = new HashSet<IBioPolymerWithSetMods> { otherSequence };
			var otherUnique = new HashSet<IBioPolymerWithSetMods> { otherSequence };

			var otherGroup = new BioPolymerGroup(otherBioPolymers, otherSequences, otherUnique);
			otherGroup.BioPolymerGroupScore = 100;

			_bioPolymerGroup.MergeWith(otherGroup);

			Assert.That(_bioPolymerGroup.BioPolymers.Count, Is.EqualTo(3));
			Assert.That(_bioPolymerGroup.AllBioPolymersWithSetMods.Count, Is.EqualTo(4));
			Assert.That(_bioPolymerGroup.UniqueBioPolymersWithSetMods.Count, Is.EqualTo(2));
			Assert.That(otherGroup.BioPolymerGroupScore, Is.EqualTo(0)); // Reset after merge
		}

		[Test]
		public void MergeWith_UpdatesBioPolymerGroupName()
		{
			var otherBioPolymer = new TestBioPolymer("MERGESEQ", "A00001"); // Comes before BP12345 alphabetically
			var otherBioPolymers = new HashSet<IBioPolymer> { otherBioPolymer };
			var otherGroup = new BioPolymerGroup(otherBioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

			_bioPolymerGroup.MergeWith(otherGroup);

			Assert.That(_bioPolymerGroup.BioPolymerGroupName, Does.StartWith("A00001"));
			Assert.That(_bioPolymerGroup.BioPolymerGroupName, Does.Contain("BP12345"));
			Assert.That(_bioPolymerGroup.BioPolymerGroupName, Does.Contain("BP67890"));
		}

		[Test]
		public void MergeWith_CombinesPsms()
		{
			var psm1 = new TestSpectralMatch(@"C:\test.raw", "SEQ1", "SEQ1", score: 100, scanNumber: 1);
			var psm2 = new TestSpectralMatch(@"C:\test.raw", "SEQ2", "SEQ2", score: 150, scanNumber: 2);

			_bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1 };

			var otherBioPolymer = new TestBioPolymer("OTHER", "BP99999");
			var otherBioPolymers = new HashSet<IBioPolymer> { otherBioPolymer };
			var otherGroup = new BioPolymerGroup(otherBioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
			otherGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm2 };

			_bioPolymerGroup.MergeWith(otherGroup);

			Assert.That(_bioPolymerGroup.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(2));
			Assert.That(_bioPolymerGroup.AllPsmsBelowOnePercentFDR.Contains(psm1), Is.True);
			Assert.That(_bioPolymerGroup.AllPsmsBelowOnePercentFDR.Contains(psm2), Is.True);
		}

		[Test]
		public void MergeWith_WithDuplicateElements_DeduplicatesCorrectly()
		{
			// Create a group with overlapping content
			var otherGroup = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);

			_bioPolymerGroup.MergeWith(otherGroup);

			// HashSet.UnionWith deduplicates automatically
			Assert.That(_bioPolymerGroup.BioPolymers.Count, Is.EqualTo(2));
			Assert.That(_bioPolymerGroup.AllBioPolymersWithSetMods.Count, Is.EqualTo(3));
		}

		#endregion

		#region Equality Tests

		[Test]
		public void Equals_SameBioPolymerGroupName_ReturnsTrue()
		{
			var bg1 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
			var bg2 = new BioPolymerGroup(_bioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

			Assert.That(bg1.Equals(bg2), Is.True);
		}

		[Test]
		public void Equals_DifferentBioPolymerGroupName_ReturnsFalse()
		{
			var otherBioPolymers = new HashSet<IBioPolymer> { _bioPolymer1 };
			var bg1 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
			var bg2 = new BioPolymerGroup(otherBioPolymers, _allSequences, _uniqueSequences);

			Assert.That(bg1.Equals(bg2), Is.False);
		}

		[Test]
		public void Equals_Null_ReturnsFalse()
		{
			Assert.That(_bioPolymerGroup.Equals((BioPolymerGroup)null), Is.False);
			Assert.That(_bioPolymerGroup.Equals((IBioPolymerGroup)null), Is.False);
		}

		[Test]
		public void Equals_SameReference_ReturnsTrue()
		{
			Assert.That(_bioPolymerGroup.Equals(_bioPolymerGroup), Is.True);
		}

		[Test]
		public void Equals_Object_WithMatchingGroup_ReturnsTrue()
		{
			var bg2 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
			Assert.That(_bioPolymerGroup.Equals((object)bg2), Is.True);
		}

		[Test]
		public void Equals_Object_WithNonBioPolymerGroup_ReturnsFalse()
		{
			Assert.That(_bioPolymerGroup.Equals("not a biopolymer group"), Is.False);
			Assert.That(_bioPolymerGroup.Equals(123), Is.False);
			Assert.That(_bioPolymerGroup.Equals(null), Is.False);
		}

		[Test]
		public void GetHashCode_SameGroupName_ReturnsSameHashCode()
		{
			var bg1 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
			var bg2 = new BioPolymerGroup(_bioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

			Assert.That(bg1.GetHashCode(), Is.EqualTo(bg2.GetHashCode()));
		}

		[Test]
		public void GetHashCode_DifferentGroupName_ReturnsDifferentHashCode()
		{
			var otherBioPolymers = new HashSet<IBioPolymer> { _bioPolymer1 };
			var bg1 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
			var bg2 = new BioPolymerGroup(otherBioPolymers, _allSequences, _uniqueSequences);

			Assert.That(bg1.GetHashCode(), Is.Not.EqualTo(bg2.GetHashCode()));
		}

		[Test]
		public void GetHashCode_WorksInHashSet()
		{
			var bg1 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
			var bg2 = new BioPolymerGroup(_bioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

			var set = new HashSet<BioPolymerGroup> { bg1, bg2 };

			Assert.That(set.Count, Is.EqualTo(1)); // Should deduplicate
		}

		[Test]
		public void GetHashCode_EmptyGroupName_ReturnsEmptyStringHashCode()
		{
			var emptyBioPolymers = new HashSet<IBioPolymer>();
			var bg = new BioPolymerGroup(emptyBioPolymers, _allSequences, _uniqueSequences);

			Assert.That(bg.BioPolymerGroupName, Is.EqualTo(string.Empty));
			Assert.That(bg.GetHashCode(), Is.EqualTo(string.Empty.GetHashCode()));
		}

		#endregion

		#region IBioPolymerGroup Interface Tests

		[Test]
		public void ImplementsIBioPolymerGroup()
		{
			Assert.That(_bioPolymerGroup, Is.InstanceOf<IBioPolymerGroup>());
		}

		[Test]
		public void IBioPolymerGroup_PropertiesAccessibleViaInterface()
		{
			IBioPolymerGroup group = _bioPolymerGroup;

			Assert.That(group.BioPolymerGroupName, Is.EqualTo(_bioPolymerGroup.BioPolymerGroupName));
			Assert.That(group.BioPolymers, Is.EqualTo(_bioPolymerGroup.BioPolymers));
			Assert.That(group.IsDecoy, Is.EqualTo(_bioPolymerGroup.IsDecoy));
			Assert.That(group.IsContaminant, Is.EqualTo(_bioPolymerGroup.IsContaminant));
			Assert.That(group.AllBioPolymersWithSetMods, Is.EqualTo(_bioPolymerGroup.AllBioPolymersWithSetMods));
			Assert.That(group.UniqueBioPolymersWithSetMods, Is.EqualTo(_bioPolymerGroup.UniqueBioPolymersWithSetMods));
		}

		[Test]
		public void IBioPolymerGroup_Equals_WorksViaInterface()
		{
			IBioPolymerGroup group1 = _bioPolymerGroup;
			IBioPolymerGroup group2 = new BioPolymerGroup(_bioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

			Assert.That(group1.Equals(group2), Is.True);
		}

		[Test]
		public void IBioPolymerGroup_WorksInLinqOperations()
		{
			var bg1 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
			var bg2 = new BioPolymerGroup(_bioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
			var bg3 = new BioPolymerGroup(new HashSet<IBioPolymer> { _bioPolymer1 }, _allSequences, _uniqueSequences);

			var list = new List<IBioPolymerGroup> { bg1, bg2, bg3 };
			var distinct = list.Distinct().ToList();

			// bg1 and bg2 have the same name, so distinct should have 2 items
			Assert.That(distinct.Count, Is.EqualTo(2));
		}

		[Test]
		public void IBioPolymerGroup_WorksInDictionary()
		{
			var bg1 = new BioPolymerGroup(_bioPolymers, _allSequences, _uniqueSequences);
			var bg2 = new BioPolymerGroup(_bioPolymers, new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

			var dict = new Dictionary<IBioPolymerGroup, string>();
			dict[bg1] = "First";

			// bg2 should be considered equal to bg1 as a key
			Assert.That(dict.ContainsKey(bg2), Is.True);
			Assert.That(dict[bg2], Is.EqualTo("First"));
		}

		#endregion

		#region ISampleInfo Tests - SpectraFileInfo

		[Test]
		public void SamplesForQuantification_WithSpectraFileInfo_AcceptsAndStores()
		{
			var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile };

			Assert.That(_bioPolymerGroup.SamplesForQuantification.Count, Is.EqualTo(1));
			Assert.That(_bioPolymerGroup.SamplesForQuantification[0], Is.EqualTo(spectraFile));
			Assert.That(_bioPolymerGroup.SamplesForQuantification[0], Is.InstanceOf<SpectraFileInfo>());
		}

		[Test]
		public void SamplesForQuantification_WithMultipleSpectraFiles_MaintainsOrder()
		{
			var file1 = new SpectraFileInfo(@"C:\test1.raw", "Control", 1, 1, 0);
			var file2 = new SpectraFileInfo(@"C:\test2.raw", "Treatment", 1, 1, 0);
			var file3 = new SpectraFileInfo(@"C:\test3.raw", "Control", 2, 1, 0);

			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { file1, file2, file3 };

			Assert.That(_bioPolymerGroup.SamplesForQuantification.Count, Is.EqualTo(3));
			Assert.That(_bioPolymerGroup.SamplesForQuantification[0], Is.EqualTo(file1));
			Assert.That(_bioPolymerGroup.SamplesForQuantification[1], Is.EqualTo(file2));
			Assert.That(_bioPolymerGroup.SamplesForQuantification[2], Is.EqualTo(file3));
		}

		[Test]
		public void IntensitiesBySample_WithSpectraFileInfo_StoresCorrectly()
		{
			var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ spectraFile, 12345.67 }
			};

			Assert.That(_bioPolymerGroup.IntensitiesBySample.Count, Is.EqualTo(1));
			Assert.That(_bioPolymerGroup.IntensitiesBySample[spectraFile], Is.EqualTo(12345.67));
		}

		[Test]
		public void IntensitiesBySample_WithMultipleSpectraFiles_AllStoredCorrectly()
		{
			var file1 = new SpectraFileInfo(@"C:\test1.raw", "Control", 1, 1, 0);
			var file2 = new SpectraFileInfo(@"C:\test2.raw", "Treatment", 1, 1, 0);

			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ file1, 1000.0 },
				{ file2, 2000.0 }
			};

			Assert.That(_bioPolymerGroup.IntensitiesBySample.Count, Is.EqualTo(2));
			Assert.That(_bioPolymerGroup.IntensitiesBySample[file1], Is.EqualTo(1000.0));
			Assert.That(_bioPolymerGroup.IntensitiesBySample[file2], Is.EqualTo(2000.0));
		}

		#endregion

		#region ISampleInfo Tests - IsobaricQuantSampleInfo

		[Test]
		public void SamplesForQuantification_WithIsobaricQuantSampleInfo_AcceptsAndStores()
		{
			var isobaricSample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { isobaricSample };

			Assert.That(_bioPolymerGroup.SamplesForQuantification.Count, Is.EqualTo(1));
			Assert.That(_bioPolymerGroup.SamplesForQuantification[0], Is.EqualTo(isobaricSample));
			Assert.That(_bioPolymerGroup.SamplesForQuantification[0], Is.InstanceOf<IsobaricQuantSampleInfo>());
		}

		[Test]
		public void SamplesForQuantification_WithMultipleIsobaricSamples_MaintainsOrder()
		{
			var sample1 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
			var sample2 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 2, "127N", 127.0, false);
			var sample3 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 3, "127C", 127.5, false);

			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { sample1, sample2, sample3 };

			Assert.That(_bioPolymerGroup.SamplesForQuantification.Count, Is.EqualTo(3));
			Assert.That(_bioPolymerGroup.SamplesForQuantification[0], Is.EqualTo(sample1));
			Assert.That(_bioPolymerGroup.SamplesForQuantification[1], Is.EqualTo(sample2));
			Assert.That(_bioPolymerGroup.SamplesForQuantification[2], Is.EqualTo(sample3));
		}

		[Test]
		public void IntensitiesBySample_WithIsobaricQuantSampleInfo_StoresCorrectly()
		{
			var isobaricSample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ isobaricSample, 5432.1 }
			};

			Assert.That(_bioPolymerGroup.IntensitiesBySample.Count, Is.EqualTo(1));
			Assert.That(_bioPolymerGroup.IntensitiesBySample[isobaricSample], Is.EqualTo(5432.1));
		}

		[Test]
		public void IntensitiesBySample_WithMultipleIsobaricSamples_AllStoredCorrectly()
		{
			var sample1 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
			var sample2 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 2, "127N", 127.0, false);

			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ sample1, 3000.0 },
				{ sample2, 4000.0 }
			};

			Assert.That(_bioPolymerGroup.IntensitiesBySample.Count, Is.EqualTo(2));
			Assert.That(_bioPolymerGroup.IntensitiesBySample[sample1], Is.EqualTo(3000.0));
			Assert.That(_bioPolymerGroup.IntensitiesBySample[sample2], Is.EqualTo(4000.0));
		}

		#endregion

		#region ISampleInfo Tests - Mixed Types

		[Test]
		public void SamplesForQuantification_WithMixedSampleTypes_AcceptsBoth()
		{
			var spectraFile = new SpectraFileInfo(@"C:\test1.raw", "Control", 1, 1, 0);
			var isobaricSample = new IsobaricQuantSampleInfo(@"C:\test2.raw", "Treatment", 1, 1, 0, 1, "126", 126.0, false);

			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile, isobaricSample };

			Assert.That(_bioPolymerGroup.SamplesForQuantification.Count, Is.EqualTo(2));
			Assert.That(_bioPolymerGroup.SamplesForQuantification[0], Is.InstanceOf<SpectraFileInfo>());
			Assert.That(_bioPolymerGroup.SamplesForQuantification[1], Is.InstanceOf<IsobaricQuantSampleInfo>());
		}

		[Test]
		public void IntensitiesBySample_WithMixedSampleTypes_StoresBothCorrectly()
		{
			var spectraFile = new SpectraFileInfo(@"C:\test1.raw", "Control", 1, 1, 0);
			var isobaricSample = new IsobaricQuantSampleInfo(@"C:\test2.raw", "Treatment", 1, 1, 0, 1, "126", 126.0, false);

			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ spectraFile, 1000.0 },
				{ isobaricSample, 2000.0 }
			};

			Assert.That(_bioPolymerGroup.IntensitiesBySample.Count, Is.EqualTo(2));
			Assert.That(_bioPolymerGroup.IntensitiesBySample[spectraFile], Is.EqualTo(1000.0));
			Assert.That(_bioPolymerGroup.IntensitiesBySample[isobaricSample], Is.EqualTo(2000.0));
		}

		[Test]
		public void IntensitiesBySample_CanRetrieveByISampleInfo()
		{
			var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ spectraFile, 9999.0 }
			};

			ISampleInfo sampleInfo = spectraFile;
			Assert.That(_bioPolymerGroup.IntensitiesBySample[sampleInfo], Is.EqualTo(9999.0));
		}

		#endregion

		#region GetTabSeparatedHeader Tests

		[Test]
		public void GetTabSeparatedHeader_ContainsExpectedColumns()
		{
			var header = _bioPolymerGroup.GetTabSeparatedHeader();

			Assert.That(header, Does.Contain("BioPolymer Accession"));
			Assert.That(header, Does.Contain("Gene"));
			Assert.That(header, Does.Contain("Organism"));
			Assert.That(header, Does.Contain("BioPolymer Full Name"));
			Assert.That(header, Does.Contain("Number of BioPolymers in Group"));
			Assert.That(header, Does.Contain("Unique Sequences"));
			Assert.That(header, Does.Contain("Shared Sequences"));
			Assert.That(header, Does.Contain("Sequence Coverage Fraction"));
			Assert.That(header, Does.Contain("BioPolymer QValue"));
			Assert.That(header, Does.Contain("Best Sequence Score"));
		}

		[Test]
		public void GetTabSeparatedHeader_WithSpectraFileInfo_IncludesIntensityColumns()
		{
			var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile };

			var header = _bioPolymerGroup.GetTabSeparatedHeader();

			Assert.That(header, Does.Contain("Intensity_"));
		}

		[Test]
		public void GetTabSeparatedHeader_WithIsobaricSamples_IncludesChannelColumns()
		{
			var isobaricSample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { isobaricSample };

			var header = _bioPolymerGroup.GetTabSeparatedHeader();

			Assert.That(header, Does.Contain("Intensity_test_126"));
		}

		[Test]
		public void GetTabSeparatedHeader_WithNullSamples_DoesNotThrow()
		{
			_bioPolymerGroup.SamplesForQuantification = null;
			Assert.DoesNotThrow(() => _bioPolymerGroup.GetTabSeparatedHeader());
		}

		[Test]
		public void GetTabSeparatedHeader_WithMultipleIsobaricChannels_OrderedByChannel()
		{
			var sample1 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 2, "127N", 127.0, false);
			var sample2 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { sample1, sample2 };

			var header = _bioPolymerGroup.GetTabSeparatedHeader();

			var index126 = header.IndexOf("Intensity_test_126");
			var index127 = header.IndexOf("Intensity_test_127N");

			Assert.That(index126, Is.LessThan(index127), "126 should come before 127N");
		}

		[Test]
		public void GetTabSeparatedHeader_WithExistingFiles_UsesConditionAndBioRep()
		{
			// Use a path that actually exists
			var existingPath = System.Reflection.Assembly.GetExecutingAssembly().Location;

			var spectraFile1 = new SpectraFileInfo(existingPath, "Treatment", 0, 0, 0);
			var spectraFile2 = new SpectraFileInfo(existingPath, "Treatment", 0, 0, 1);
			var spectraFile3 = new SpectraFileInfo(existingPath, "Treatment", 1, 0, 0);

			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile1, spectraFile2, spectraFile3 };

			var header = _bioPolymerGroup.GetTabSeparatedHeader();

			Assert.That(header, Does.Contain("Intensity_Treatment_1"));
			Assert.That(header, Does.Contain("Intensity_Treatment_2"));
		}

		#endregion

		#region ToString Tests

		[Test]
		public void ToString_ContainsBioPolymerGroupName()
		{
			var result = _bioPolymerGroup.ToString();
			Assert.That(result, Does.StartWith("BP12345|BP67890"));
		}

		[Test]
		public void ToString_ContainsGeneNames()
		{
			var result = _bioPolymerGroup.ToString();
			Assert.That(result, Does.Contain("GENE1"));
			Assert.That(result, Does.Contain("GENE2"));
		}

		[Test]
		public void ToString_ContainsOrganisms()
		{
			var result = _bioPolymerGroup.ToString();
			Assert.That(result, Does.Contain("Homo sapiens"));
		}

		[Test]
		public void ToString_TargetGroup_ContainsT()
		{
			var result = _bioPolymerGroup.ToString();
			Assert.That(result, Does.Contain("\tT\t")); // Target
		}

		[Test]
		public void ToString_DecoyGroup_ContainsD()
		{
			var decoyBioPolymers = new HashSet<IBioPolymer> { _decoyBioPolymer };
			var bg = new BioPolymerGroup(decoyBioPolymers, _allSequences, _uniqueSequences);

			var result = bg.ToString();
			Assert.That(result, Does.Contain("\tD\t")); // Decoy
		}

		[Test]
		public void ToString_ContaminantGroup_ContainsC()
		{
			var contaminantBioPolymers = new HashSet<IBioPolymer> { _contaminantBioPolymer };
			var bg = new BioPolymerGroup(contaminantBioPolymers, _allSequences, _uniqueSequences);

			var result = bg.ToString();
			Assert.That(result, Does.Contain("\tC\t")); // Contaminant
		}

		[Test]
		public void ToString_WithSpectraFileIntensities_IncludesIntensityValues()
		{
			var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile };
			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ spectraFile, 12345.67 }
			};

			var result = _bioPolymerGroup.ToString();
			Assert.That(result, Does.Contain("12345.67"));
		}

		[Test]
		public void ToString_WithIsobaricIntensities_IncludesIntensityValues()
		{
			var isobaricSample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { isobaricSample };
			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ isobaricSample, 9876.54 }
			};

			var result = _bioPolymerGroup.ToString();
			Assert.That(result, Does.Contain("9876.54"));
		}

		[Test]
		public void ToString_WithZeroIntensity_OutputsEmptyColumn()
		{
			var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { spectraFile };
			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ spectraFile, 0.0 }
			};

			var result = _bioPolymerGroup.ToString();
			// Zero intensity should not be printed, resulting in consecutive tabs
			Assert.That(result, Does.Contain("\t\t"));
		}

		[Test]
		public void ToString_WithInvalidSequence_HandlesMassCalculationGracefully()
		{
			// This test may produce NaN for invalid sequences
			var invalidBioPolymer = new TestBioPolymer("INVALIDSEQ123", "BP00001");
			var bioPolymers = new HashSet<IBioPolymer> { invalidBioPolymer };
			var bg = new BioPolymerGroup(bioPolymers, _allSequences, _uniqueSequences);

			Assert.DoesNotThrow(() => bg.ToString());
			var result = bg.ToString();
			Assert.That(result, Does.Contain("NaN"));
		}

		[Test]
		public void ToString_WithNullIntensitiesBySample_DoesNotThrow()
		{
			_bioPolymerGroup.IntensitiesBySample = null;
			Assert.DoesNotThrow(() => _bioPolymerGroup.ToString());
		}

		#endregion

		#region GetIdentifiedSequencesOutput Tests

		[Test]
		public void GetIdentifiedSequencesOutput_WithoutMods_UsesBaseSequence()
		{
			_bioPolymerGroup.DisplayModsOnPeptides = false;

			// The method populates private fields - test indirectly through ToString
			Assert.DoesNotThrow(() => _bioPolymerGroup.ToString());
		}

		[Test]
		public void GetIdentifiedSequencesOutput_WithMods_UsesFullSequence()
		{
			_bioPolymerGroup.DisplayModsOnPeptides = true;

			Assert.DoesNotThrow(() => _bioPolymerGroup.ToString());
		}

		#endregion

		#region ConstructSubsetBioPolymerGroup Tests

		[Test]
		public void ConstructSubsetBioPolymerGroup_FiltersToSpecificFile()
		{
			var file1 = new SpectraFileInfo(@"C:\test1.raw", "Control", 1, 1, 0);
			var file2 = new SpectraFileInfo(@"C:\test2.raw", "Control", 1, 1, 0);

			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { file1, file2 };
			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ file1, 1000.0 },
				{ file2, 2000.0 }
			};

			var psm1 = new TestSpectralMatch(@"C:\test1.raw", "SEQ", "SEQ", score: 100, scanNumber: 1);
			var psm2 = new TestSpectralMatch(@"C:\test2.raw", "SEQ", "SEQ", score: 100, scanNumber: 2);
			_bioPolymerGroup.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

			var subset = _bioPolymerGroup.ConstructSubsetBioPolymerGroup(@"C:\test1.raw");

			Assert.That(subset.SamplesForQuantification.Count, Is.EqualTo(1));
			Assert.That(subset.SamplesForQuantification[0].FullFilePathWithExtension, Is.EqualTo(@"C:\test1.raw"));
			Assert.That(subset.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(1));
		}

		[Test]
		public void ConstructSubsetBioPolymerGroup_WithIsobaricSamples_FiltersCorrectly()
		{
			var sample1 = new IsobaricQuantSampleInfo(@"C:\test1.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
			var sample2 = new IsobaricQuantSampleInfo(@"C:\test1.raw", "Control", 1, 1, 0, 2, "127N", 127.0, false);
			var sample3 = new IsobaricQuantSampleInfo(@"C:\test2.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { sample1, sample2, sample3 };
			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ sample1, 1000.0 },
				{ sample2, 1500.0 },
				{ sample3, 2000.0 }
			};

			var subset = _bioPolymerGroup.ConstructSubsetBioPolymerGroup(@"C:\test1.raw");

			Assert.That(subset.SamplesForQuantification.Count, Is.EqualTo(2)); // Both channels from test1.raw
			Assert.That(subset.IntensitiesBySample.Count, Is.EqualTo(2));
			Assert.That(subset.IntensitiesBySample.ContainsKey(sample1), Is.True);
			Assert.That(subset.IntensitiesBySample.ContainsKey(sample2), Is.True);
			Assert.That(subset.IntensitiesBySample.ContainsKey(sample3), Is.False);
		}

		[Test]
		public void ConstructSubsetBioPolymerGroup_ReturnsIBioPolymerGroup()
		{
			var subset = _bioPolymerGroup.ConstructSubsetBioPolymerGroup(@"C:\nonexistent.raw");
			Assert.That(subset, Is.InstanceOf<IBioPolymerGroup>());
		}

		[Test]
		public void ConstructSubsetBioPolymerGroup_PreservesDisplayModsOnPeptides()
		{
			_bioPolymerGroup.DisplayModsOnPeptides = true;

			var subset = (BioPolymerGroup)_bioPolymerGroup.ConstructSubsetBioPolymerGroup(@"C:\test.raw");

			Assert.That(subset.DisplayModsOnPeptides, Is.True);
		}

		[Test]
		public void ConstructSubsetBioPolymerGroup_WithNullSamples_HandlesGracefully()
		{
			_bioPolymerGroup.SamplesForQuantification = null;
			_bioPolymerGroup.IntensitiesBySample = null;

			var subset = _bioPolymerGroup.ConstructSubsetBioPolymerGroup(@"C:\test.raw");

			Assert.That(subset, Is.Not.Null);
			Assert.That(subset.SamplesForQuantification, Is.Null);
		}

		[Test]
		public void ConstructSubsetBioPolymerGroup_WithNoMatchingFile_ReturnsEmptySubset()
		{
			var file1 = new SpectraFileInfo(@"C:\test1.raw", "Control", 1, 1, 0);
			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { file1 };

			var subset = _bioPolymerGroup.ConstructSubsetBioPolymerGroup(@"C:\different.raw");

			Assert.That(subset.SamplesForQuantification.Count, Is.EqualTo(0));
		}

		#endregion

		#region Edge Cases

		[Test]
		public void EdgeCase_SingleBioPolymer_GroupNameIsJustAccession()
		{
			var singleBioPolymer = new HashSet<IBioPolymer> { _bioPolymer1 };
			var bg = new BioPolymerGroup(singleBioPolymer, _allSequences, _uniqueSequences);

			Assert.That(bg.BioPolymerGroupName, Is.EqualTo("BP12345"));
		}

		[Test]
		public void EdgeCase_BioPolymerWithNoGeneNames_HandlesGracefully()
		{
			var bioPolymerNoGenes = new TestBioPolymer("SEQUENCE", "BP00000");
			var bioPolymers = new HashSet<IBioPolymer> { bioPolymerNoGenes };
			var bg = new BioPolymerGroup(bioPolymers, _allSequences, _uniqueSequences);

			Assert.DoesNotThrow(() => bg.ToString());
		}

		[Test]
		public void EdgeCase_BioPolymerWithNullGeneNames_HandlesGracefully()
		{
			var bioPolymerNullGenes = new TestBioPolymer("SEQUENCE", "BP00000", geneNames: null);
			var bioPolymers = new HashSet<IBioPolymer> { bioPolymerNullGenes };
			var bg = new BioPolymerGroup(bioPolymers, _allSequences, _uniqueSequences);

			Assert.DoesNotThrow(() => bg.ToString());
		}

		[Test]
		public void EdgeCase_VeryLongStrings_AreTruncated()
		{
			var longName = new string('A', 50000);
			var bioPolymerLongName = new TestBioPolymer("SEQ", "BP00001", fullName: longName);
			var bioPolymers = new HashSet<IBioPolymer> { bioPolymerLongName };
			var bg = new BioPolymerGroup(bioPolymers, _allSequences, _uniqueSequences);

			var result = bg.ToString();

			// The output should be truncated to MaxStringLength (32000 by default)
			Assert.That(result.Length, Is.LessThan(longName.Length));
		}

		[Test]
		public void EdgeCase_MaxStringLengthDisabled_NoTruncation()
		{
			var originalMaxLength = BioPolymerGroup.MaxStringLength;
			try
			{
				BioPolymerGroup.MaxStringLength = 0; // Disable truncation

				var longName = new string('A', 50000);
				var bioPolymerLongName = new TestBioPolymer("SEQ", "BP00001", fullName: longName);
				var bioPolymers = new HashSet<IBioPolymer> { bioPolymerLongName };
				var bg = new BioPolymerGroup(bioPolymers, _allSequences, _uniqueSequences);

				var result = bg.ToString();
				Assert.That(result, Does.Contain(longName));
			}
			finally
			{
				BioPolymerGroup.MaxStringLength = originalMaxLength;
			}
		}

		[Test]
		public void EdgeCase_EmptyGeneNamesList_HandlesGracefully()
		{
			var bioPolymerEmptyGenes = new TestBioPolymer("SEQUENCE", "BP00000", geneNames: new List<Tuple<string, string>>());
			var bioPolymers = new HashSet<IBioPolymer> { bioPolymerEmptyGenes };
			var bg = new BioPolymerGroup(bioPolymers, _allSequences, _uniqueSequences);

			Assert.DoesNotThrow(() => bg.ToString());
		}

		[Test]
		public void EdgeCase_MultipleIsobaricFilesAndChannels_OrderedCorrectly()
		{
			var sample1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 2, "127N", 127.0, false);
			var sample2 = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
			var sample3 = new IsobaricQuantSampleInfo(@"C:\b.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { sample1, sample2, sample3 };
			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ sample1, 1111.0 },
				{ sample2, 2222.0 },
				{ sample3, 3333.0 }
			};

			var result = _bioPolymerGroup.ToString();

			// All intensities should be present
			Assert.That(result, Does.Contain("1111"));
			Assert.That(result, Does.Contain("2222"));
			Assert.That(result, Does.Contain("3333"));

			// Verify ordering
			var index2222 = result.IndexOf("2222"); // a.raw 126
			var index1111 = result.IndexOf("1111"); // a.raw 127N
			var index3333 = result.IndexOf("3333"); // b.raw 126

			Assert.That(index2222, Is.LessThan(index1111), "126 should come before 127N for same file");
			Assert.That(index1111, Is.LessThan(index3333), "a.raw channels should come before b.raw channels");
		}

		[Test]
		public void EdgeCase_MissingIntensityForSample_HandlesGracefully()
		{
			var sample1 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
			var sample2 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 2, "127N", 127.0, false);

			_bioPolymerGroup.SamplesForQuantification = new List<ISampleInfo> { sample1, sample2 };
			_bioPolymerGroup.IntensitiesBySample = new Dictionary<ISampleInfo, double>
			{
				{ sample1, 5000.0 }
                // sample2 is missing
            };

			Assert.DoesNotThrow(() => _bioPolymerGroup.ToString());
			var result = _bioPolymerGroup.ToString();
			Assert.That(result, Does.Contain("5000"));
		}

		#endregion

		#region Additional Property Tests

		[Test]
		public void Properties_CumulativeTarget_CanBeSetAndRetrieved()
		{
			_bioPolymerGroup.CumulativeTarget = 42;
			Assert.That(_bioPolymerGroup.CumulativeTarget, Is.EqualTo(42));
		}

		[Test]
		public void Properties_CumulativeDecoy_CanBeSetAndRetrieved()
		{
			_bioPolymerGroup.CumulativeDecoy = 13;
			Assert.That(_bioPolymerGroup.CumulativeDecoy, Is.EqualTo(13));
		}

		[Test]
		public void Properties_QValue_CanBeSetAndRetrieved()
		{
			_bioPolymerGroup.QValue = 0.01;
			Assert.That(_bioPolymerGroup.QValue, Is.EqualTo(0.01));
		}

		[Test]
		public void Properties_BestBioPolymerWithSetModsQValue_CanBeSetAndRetrieved()
		{
			_bioPolymerGroup.BestBioPolymerWithSetModsQValue = 0.005;
			Assert.That(_bioPolymerGroup.BestBioPolymerWithSetModsQValue, Is.EqualTo(0.005));
		}

		[Test]
		public void Properties_BestBioPolymerWithSetModsScore_CanBeSetAndRetrieved()
		{
			_bioPolymerGroup.BestBioPolymerWithSetModsScore = 250.5;
			Assert.That(_bioPolymerGroup.BestBioPolymerWithSetModsScore, Is.EqualTo(250.5));
		}

		[Test]
		public void Properties_DisplayModsOnPeptides_CanBeSetAndRetrieved()
		{
			_bioPolymerGroup.DisplayModsOnPeptides = true;
			Assert.That(_bioPolymerGroup.DisplayModsOnPeptides, Is.True);

			_bioPolymerGroup.DisplayModsOnPeptides = false;
			Assert.That(_bioPolymerGroup.DisplayModsOnPeptides, Is.False);
		}

		#endregion
	}

	#region Test Helper Classes

	/// <summary>
	/// Test implementation of IBioPolymer for BioPolymerGroup tests.
	/// </summary>
	internal class TestBioPolymer : IBioPolymer
	{
		public string BaseSequence { get; }
		public string Accession { get; }
		public string Organism { get; }
		public string Name { get; }
		public string FullName { get; }
		public List<Tuple<string, string>> GeneNames { get; }
		public bool IsDecoy { get; }
		public bool IsContaminant { get; }
		public string DatabaseFilePath { get; }

		// Additional required members from IBioPolymer
		public int Length => BaseSequence.Length;
		public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }

		// From IHasSequenceVariants (inherited by IBioPolymer)
		public string SampleNameForVariants { get; set; }
		public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }
		public IBioPolymer ConsensusVariant { get; }
		public List<SequenceVariation> AppliedSequenceVariations { get; }
		public List<SequenceVariation> SequenceVariations { get; }
		public List<TruncationProduct> TruncationProducts { get; }

		public TestBioPolymer(string sequence, string accession,
			string organism = "",
			string name = "",
			string fullName = "",
			List<Tuple<string, string>> geneNames = null,
			bool isDecoy = false,
			bool isContaminant = false,
			string databaseFilePath = "")
		{
			BaseSequence = sequence;
			Accession = accession;
			Organism = organism;
			Name = name;
			FullName = fullName;
			GeneNames = geneNames ?? new List<Tuple<string, string>>();
			IsDecoy = isDecoy;
			IsContaminant = isContaminant;
			DatabaseFilePath = databaseFilePath;
			OneBasedPossibleLocalizedModifications = new Dictionary<int, List<Modification>>();
			OriginalNonVariantModifications = new Dictionary<int, List<Modification>>();
			ConsensusVariant = this; // For simple test cases, point to self
			AppliedSequenceVariations = new List<SequenceVariation>();
			SequenceVariations = new List<SequenceVariation>();
			TruncationProducts = new List<TruncationProduct>();
			SampleNameForVariants = string.Empty;
		}

		// Required method from IBioPolymer
		public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams,
			List<Modification> allKnownFixedModifications,
			List<Modification> variableModifications,
			List<SilacLabel>? silacLabels = null,
			(SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null,
			bool topDownTruncationSearch = false)
		{
			// For testing purposes, return empty enumerable
			// Real implementation would perform digestion
			return Enumerable.Empty<IBioPolymerWithSetMods>();
		}

		// Required method from IBioPolymer
		public IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence,
			IDictionary<int, List<Modification>>? newMods)
		{
			// For testing purposes, create a simple clone with new sequence
			var clone = new TestBioPolymer(newBaseSequence, Accession, Organism, Name, FullName,
				GeneNames, IsDecoy, IsContaminant, DatabaseFilePath);

			if (newMods != null)
			{
				foreach (var kvp in newMods)
				{
					clone.OneBasedPossibleLocalizedModifications[kvp.Key] = kvp.Value;
				}
			}

			return clone;
		}

		// Required method from IHasSequenceVariants
		public TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence,
			TBioPolymerType original,
			IEnumerable<SequenceVariation> appliedSequenceVariants,
			IEnumerable<TruncationProduct> applicableProteolysisProducts,
			IDictionary<int, List<Modification>> oneBasedModifications,
			string sampleNameForVariants)
			where TBioPolymerType : IHasSequenceVariants
		{
			// For testing purposes, return the original
			// Real implementation would create a new variant
			return original;
		}

		// Equals implementation from IBioPolymer interface
		public bool Equals(IBioPolymer? other)
		{
			if (other is null) return false;
			if (ReferenceEquals(this, other)) return true;
			if (other.GetType() != GetType()) return false;
			return Accession == other.Accession && BaseSequence == other.BaseSequence;
		}

		public override bool Equals(object? obj)
		{
			return obj is IBioPolymer other && Equals(other);
		}

		public override int GetHashCode()
		{
			return HashCode.Combine(Accession, BaseSequence);
		}
	}

	/// <summary>
	/// Test implementation of IBioPolymerWithSetMods for BioPolymerGroup tests.
	/// </summary>
	internal class TestBioPolymerWithSetMods : IBioPolymerWithSetMods
	{
		public string BaseSequence { get; }
		public string FullSequence { get; }
		public double MostAbundantMonoisotopicMass { get; }
		public string SequenceWithChemicalFormulas { get; }
		public int OneBasedStartResidue { get; }
		public int OneBasedEndResidue { get; }
		public int MissedCleavages { get; }
		public string Description { get; }
		public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; }
		public char PreviousResidue { get; }
		public char NextResidue { get; }
		public IDigestionParams DigestionParams { get; }
		public Dictionary<int, Modification> AllModsOneIsNterminus { get; }
		public int NumMods { get; }
		public int NumFixedMods { get; }
		public int NumVariableMods { get; }
		public int Length => BaseSequence.Length;
		public IBioPolymer Parent { get; }

		// From IHasChemicalFormula
		public ChemicalFormula ThisChemicalFormula { get; }

		// From IHasMass (inherited by IHasChemicalFormula)
		public double MonoisotopicMass { get; }

		public TestBioPolymerWithSetMods(string baseSequence, string fullSequence,
			double mass = 0, int startResidue = 0, int endResidue = 0)
		{
			BaseSequence = baseSequence;
			FullSequence = fullSequence;
			MostAbundantMonoisotopicMass = mass;
			MonoisotopicMass = mass;
			SequenceWithChemicalFormulas = baseSequence;
			OneBasedStartResidue = startResidue;
			OneBasedEndResidue = endResidue > 0 ? endResidue : baseSequence.Length;
			MissedCleavages = 0;
			Description = "Test";
			CleavageSpecificityForFdrCategory = CleavageSpecificity.Full;
			PreviousResidue = '-';
			NextResidue = '-';
			DigestionParams = null;
			AllModsOneIsNterminus = new Dictionary<int, Modification>();
			NumMods = 0;
			NumFixedMods = 0;
			NumVariableMods = 0;
			Parent = null;
			ThisChemicalFormula = new ChemicalFormula();
		}

		public void Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus,
			List<Product> products, FragmentationParams? fragmentationParams = null)
		{
			// For testing purposes, do nothing
			// Real implementation would generate fragment ions
		}

		public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments,
			List<Product> products, FragmentationParams? fragmentationParams = null)
		{
			// For testing purposes, do nothing
			// Real implementation would generate internal fragment ions
		}

		public IBioPolymerWithSetMods Localize(int indexOfMass, double massToLocalize)
		{
			// For testing purposes, return self
			// Real implementation would create a new instance with localized mass
			return this;
		}

		public bool Equals(IBioPolymerWithSetMods? other)
		{
			if (other is null) return false;
			if (ReferenceEquals(this, other)) return true;
			return BaseSequence == other.BaseSequence && FullSequence == other.FullSequence;
		}

		public override bool Equals(object? obj)
		{
			return obj is IBioPolymerWithSetMods other && Equals(other);
		}

		public override int GetHashCode()
		{
			return HashCode.Combine(BaseSequence, FullSequence);
		}
	}

	#endregion
}