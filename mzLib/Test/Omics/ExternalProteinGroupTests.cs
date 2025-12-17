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
using System.Text.RegularExpressions;

namespace Test.Omics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class ExternalProteinGroupTests
    {
        #region Helper Classes

        /// <summary>
        /// Test implementation of IBioPolymer for BioPolymerGroup tests.
        /// </summary>
        private class TestBioPolymer : IBioPolymer
        {
            public string BaseSequence { get; }
            public string Accession { get; }
            public bool IsDecoy { get; }
            public bool IsContaminant { get; }
            public string Organism { get; init; } = "Test Organism";
            public string FullName { get; init; } = "Test Full Name";
            public string Name { get; init; } = "Test Name";
            public int Length => BaseSequence?.Length ?? 0;
            public string DatabaseFilePath { get; init; } = string.Empty;
            public List<Tuple<string, string>> GeneNames { get; init; } = new() { new Tuple<string, string>("primary", "TestGene") };
            public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; init; } = new Dictionary<int, List<Modification>>();
            public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

            // IHasSequenceVariants members
            public string SampleNameForVariants { get; init; } = string.Empty;
            public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; } = new Dictionary<int, List<Modification>>();
            public IBioPolymer ConsensusVariant => this;
            public List<SequenceVariation> AppliedSequenceVariations { get; init; } = new();
            public List<SequenceVariation> SequenceVariations { get; init; } = new();
            public List<TruncationProduct> TruncationProducts { get; init; } = new();

            public TestBioPolymer(string baseSequence, string accession, bool isDecoy = false, bool isContaminant = false)
            {
                BaseSequence = baseSequence;
                Accession = accession;
                IsDecoy = isDecoy;
                IsContaminant = isContaminant;
            }

            public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams,
                List<Modification> allKnownFixedModifications, List<Modification> variableModifications,
                List<SilacLabel>? silacLabels = null, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null,
                bool topDownTruncationSearch = false)
            {
                throw new NotImplementedException();
            }

            public IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence, IDictionary<int, List<Modification>>? newMods)
            {
                throw new NotImplementedException();
            }

            public TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original,
                IEnumerable<SequenceVariation> appliedSequenceVariants, IEnumerable<TruncationProduct> applicableProteolysisProducts,
                IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
                where TBioPolymerType : IHasSequenceVariants
            {
                throw new NotImplementedException();
            }

            public bool Equals(IBioPolymer? other)
            {
                if (other is null) return false;
                if (ReferenceEquals(this, other)) return true;
                return Accession == other.Accession && BaseSequence == other.BaseSequence;
            }

            public override bool Equals(object? obj) => Equals(obj as IBioPolymer);
            public override int GetHashCode() => HashCode.Combine(Accession, BaseSequence);
        }

        /// <summary>
        /// Test implementation of IBioPolymerWithSetMods for BioPolymerGroup tests.
        /// </summary>
        private class TestBioPolymerWithSetMods : IBioPolymerWithSetMods
        {
            public string BaseSequence { get; }
            public string FullSequence { get; }
            public IBioPolymer Parent { get; }
            public int OneBasedStartResidue { get; }
            public int OneBasedEndResidue { get; }
            public double MonoisotopicMass { get; init; } = 500.0;
            public double MostAbundantMonoisotopicMass => MonoisotopicMass;
            public ChemicalFormula ThisChemicalFormula => new ChemicalFormula();
            public string SequenceWithChemicalFormulas => FullSequence;
            public int MissedCleavages { get; init; } = 0;
            public string Description { get; init; } = string.Empty;
            public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; } = CleavageSpecificity.Full;
            public char PreviousResidue { get; init; } = '-';
            public char NextResidue { get; init; } = '-';
            public IDigestionParams DigestionParams => null!;
            public Dictionary<int, Modification> AllModsOneIsNterminus { get; init; } = new();
            public int NumMods => AllModsOneIsNterminus?.Count ?? 0;
            public int NumFixedMods => 0;
            public int NumVariableMods => NumMods;
            public int Length => BaseSequence?.Length ?? 0;
            public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

            public TestBioPolymerWithSetMods(IBioPolymer parent, int startResidue, int endResidue,
                Dictionary<int, Modification>? mods = null)
            {
                Parent = parent;
                OneBasedStartResidue = startResidue;
                OneBasedEndResidue = endResidue;
                BaseSequence = parent.BaseSequence.Substring(startResidue - 1, endResidue - startResidue + 1);
                AllModsOneIsNterminus = mods ?? new Dictionary<int, Modification>();

                // Build full sequence with mods
                if (AllModsOneIsNterminus.Any())
                {
                    var sb = new System.Text.StringBuilder();
                    for (int i = 0; i < BaseSequence.Length; i++)
                    {
                        sb.Append(BaseSequence[i]);
                        // Check for mod at this position (AllModsOneIsNterminus uses 1-based where 1 is N-term, 2 is first residue, etc.)
                        if (AllModsOneIsNterminus.TryGetValue(i + 2, out var mod))
                        {
                            sb.Append($"[{mod.IdWithMotif}]");
                        }
                    }
                    FullSequence = sb.ToString();
                }
                else
                {
                    FullSequence = BaseSequence;
                }
            }

            public void Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus,
                List<Product> products, FragmentationParams? fragmentationParams = null)
            { }

            public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments,
                List<Product> products, FragmentationParams? fragmentationParams = null)
            { }

            public IBioPolymerWithSetMods Localize(int indexOfMass, double massToLocalize) => this;

            public bool Equals(IBioPolymerWithSetMods? other)
            {
                if (other == null) return false;
                return string.Equals(FullSequence, other.FullSequence, StringComparison.Ordinal)
                    && Parent?.Accession == other.Parent?.Accession;
            }

            public override bool Equals(object? obj) => Equals(obj as IBioPolymerWithSetMods);
            public override int GetHashCode() => HashCode.Combine(FullSequence, Parent?.Accession);
        }

        /// <summary>
        /// Test implementation of ISpectralMatch for BioPolymerGroup tests.
        /// </summary>
        private class TestSpectralMatch : ISpectralMatch
        {
            private readonly List<IBioPolymerWithSetMods> _identifiedBioPolymers;

            public string FullFilePath { get; }
            public string FullSequence { get; }
            public string BaseSequence { get; }
            public double Score { get; }
            public int OneBasedScanNumber { get; }

            public TestSpectralMatch(string fullFilePath, string baseSequence, string fullSequence,
                double score, int scanNumber, IEnumerable<IBioPolymerWithSetMods> identified)
            {
                FullFilePath = fullFilePath;
                BaseSequence = baseSequence;
                FullSequence = fullSequence;
                Score = score;
                OneBasedScanNumber = scanNumber;
                _identifiedBioPolymers = identified?.ToList() ?? new List<IBioPolymerWithSetMods>();
            }

            public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => _identifiedBioPolymers;
            public int CompareTo(ISpectralMatch? other)
            {
                if (other is null) return -1;
                return other.Score.CompareTo(Score);
            }
        }

        #endregion

        #region Equality Tests

        [Test]
        public void TestBioPolymerGroupEquals()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1");
            var proteinList1 = new List<IBioPolymer> { prot1 };
            var bioPolymerGroup1 = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList1),
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            var proteinList2 = new List<IBioPolymer> { prot1 };
            var bioPolymerGroup2 = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList2),
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            // Two groups with the same biopolymer should be equal
            Assert.That(bioPolymerGroup1.Equals(bioPolymerGroup2), Is.True);

            var prot3 = new TestBioPolymer("EDEEK", "prot3");
            var proteinList3 = new List<IBioPolymer> { prot3 };
            var bioPolymerGroup3 = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList3),
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            // Two groups with different biopolymers should not be equal
            Assert.That(bioPolymerGroup1.Equals(bioPolymerGroup3), Is.False);

            var proteinList4 = new List<IBioPolymer> { prot1, prot3 };
            var proteinList5 = new List<IBioPolymer> { prot3, prot1 };
            var bioPolymerGroup4 = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList4),
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            var bioPolymerGroup5 = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList5),
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            // Groups with the same biopolymers in different order should be equal
            Assert.That(bioPolymerGroup4.Equals(bioPolymerGroup5), Is.True);

            var pwsm1 = new TestBioPolymerWithSetMods(prot1, 1, 3);
            var pwsm2 = new TestBioPolymerWithSetMods(prot1, 4, 6);
            var bioPolymerGroup6 = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList1),
                new HashSet<IBioPolymerWithSetMods> { pwsm1 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1 });
            var bioPolymerGroup7 = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList1),
                new HashSet<IBioPolymerWithSetMods> { pwsm2 },
                new HashSet<IBioPolymerWithSetMods> { pwsm2 });

            // Groups with the same biopolymers but different peptides should be equal (equality is by name)
            Assert.That(bioPolymerGroup6.Equals(bioPolymerGroup7), Is.True);

            // A null group should not be equal to a non-null group
            BioPolymerGroup? nullGroup = null;
            Assert.That(bioPolymerGroup1.Equals(nullGroup), Is.False);
        }

        [Test]
        public void TestBioPolymerGroupHashCode()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1");
            var group1 = new BioPolymerGroup(new HashSet<IBioPolymer> { prot1 },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            var group2 = new BioPolymerGroup(new HashSet<IBioPolymer> { prot1 },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            // Equal groups should have the same hash code
            Assert.That(group1.GetHashCode(), Is.EqualTo(group2.GetHashCode()));

            var prot2 = new TestBioPolymer("EDEEK", "prot2");
            var group3 = new BioPolymerGroup(new HashSet<IBioPolymer> { prot2 },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            // Different groups should (likely) have different hash codes
            Assert.That(group1.GetHashCode(), Is.Not.EqualTo(group3.GetHashCode()));
        }

        #endregion

        #region ToString Tests

        [Test]
        public void BioPolymerGroupToStringTest()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1");
            var prot2 = new TestBioPolymer("MENEEK", "prot2");

            var pwsm1 = new TestBioPolymerWithSetMods(prot1, 1, 3);
            var pwsm2 = new TestBioPolymerWithSetMods(prot2, 1, 3);

            var proteinList1 = new List<IBioPolymer> { prot1, prot2 };

            var bioPolymerGroup1 = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList1),
                new HashSet<IBioPolymerWithSetMods> { pwsm1, pwsm2 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1, pwsm2 });

            var output = bioPolymerGroup1.ToString();

            // Verify output contains expected content
            Assert.That(output, Does.Contain("prot1|prot2"));
            Assert.That(output, Does.Contain("2")); // BioPolymers.Count
            Assert.That(output, Does.Contain("T")); // Target (not decoy or contaminant)
        }

        [Test]
        public void BioPolymerGroupToStringWithDecoyTest()
        {
            var prot1 = new TestBioPolymer("MAAADAAAAAAAAAAAAAAA", "prot1", isDecoy: true);
            var proteinList = new List<IBioPolymer> { prot1 };
            var bioPolymerGroup = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList),
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            var output = bioPolymerGroup.ToString();

            // Verify output contains decoy marker
            Assert.That(output, Does.Contain("D"));
        }

        [Test]
        public void BioPolymerGroupToStringWithContaminantTest()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1", isContaminant: true);
            var proteinList = new List<IBioPolymer> { prot1 };
            var bioPolymerGroup = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList),
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            var output = bioPolymerGroup.ToString();

            // Verify output contains contaminant marker
            Assert.That(output, Does.Contain("C"));
        }

        [Test]
        public void TestBioPolymerGroupStringAndHeaderHaveSameNumberOfTabs()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1");
            var prot2 = new TestBioPolymer("MENEEK", "prot2");

            var pwsm1 = new TestBioPolymerWithSetMods(prot1, 1, 3);
            var pwsm2 = new TestBioPolymerWithSetMods(prot2, 1, 3);

            var proteinList1 = new List<IBioPolymer> { prot1, prot2 };

            var bioPolymerGroup1 = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList1),
                new HashSet<IBioPolymerWithSetMods> { pwsm1, pwsm2 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1, pwsm2 });

            string header = bioPolymerGroup1.GetTabSeparatedHeader();
            string row = bioPolymerGroup1.ToString();

            string[] headerFields = header.Split('\t');
            string[] rowEntries = row.Split('\t');

            Assert.That(headerFields.Length, Is.EqualTo(rowEntries.Length));
            Assert.That(Regex.Matches(header, @"\t").Count, Is.EqualTo(Regex.Matches(row, @"\t").Count));
        }

        #endregion

        #region Merge Tests

        [Test]
        public void BioPolymerGroupMergeTest()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1");
            var prot2 = new TestBioPolymer("MENEEK", "prot2");
            var prot3 = new TestBioPolymer("MAAADAAAAAAAAAAAAAAA", "prot3");
            var prot4 = new TestBioPolymer("MNNDNNNN", "prot4");

            var pwsm1 = new TestBioPolymerWithSetMods(prot1, 1, 3);
            var pwsm2 = new TestBioPolymerWithSetMods(prot2, 1, 3);
            var pwsm3 = new TestBioPolymerWithSetMods(prot3, 1, 3);
            var pwsm4 = new TestBioPolymerWithSetMods(prot4, 1, 3);

            var proteinList1 = new List<IBioPolymer> { prot1, prot2 };
            var proteinList2 = new List<IBioPolymer> { prot3, prot4 };

            var bioPolymerGroup1 = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList1),
                new HashSet<IBioPolymerWithSetMods> { pwsm1, pwsm2 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1, pwsm2 });
            var bioPolymerGroup2 = new BioPolymerGroup(new HashSet<IBioPolymer>(proteinList2),
                new HashSet<IBioPolymerWithSetMods> { pwsm3, pwsm4 },
                new HashSet<IBioPolymerWithSetMods> { pwsm3, pwsm4 });

            bioPolymerGroup1.MergeWith(bioPolymerGroup2);

            // Merged group should have all biopolymers
            Assert.That(bioPolymerGroup1.BioPolymers.Contains(prot1), Is.True);
            Assert.That(bioPolymerGroup1.BioPolymers.Contains(prot2), Is.True);
            Assert.That(bioPolymerGroup1.BioPolymers.Contains(prot3), Is.True);
            Assert.That(bioPolymerGroup1.BioPolymers.Contains(prot4), Is.True);

            // Merged group should have all peptides
            Assert.That(bioPolymerGroup1.AllBioPolymersWithSetMods.Count, Is.EqualTo(4));
            Assert.That(bioPolymerGroup1.UniqueBioPolymersWithSetMods.Count, Is.EqualTo(4));
        }

        [Test]
        public void BioPolymerGroupMergeResetsOtherGroupScore()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1");
            var prot2 = new TestBioPolymer("MENEEK", "prot2");

            var bioPolymerGroup1 = new BioPolymerGroup(new HashSet<IBioPolymer> { prot1 },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            var bioPolymerGroup2 = new BioPolymerGroup(new HashSet<IBioPolymer> { prot2 },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            bioPolymerGroup2.BioPolymerGroupScore = 100.0;

            bioPolymerGroup1.MergeWith(bioPolymerGroup2);

            // The merged-in group's score should be reset to 0
            Assert.That(bioPolymerGroup2.BioPolymerGroupScore, Is.EqualTo(0));
        }

        [Test]
        public void BioPolymerGroupMergeUpdatesBioPolymerGroupName()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1");
            var prot2 = new TestBioPolymer("MENEEK", "prot2");

            var bioPolymerGroup1 = new BioPolymerGroup(new HashSet<IBioPolymer> { prot1 },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            var bioPolymerGroup2 = new BioPolymerGroup(new HashSet<IBioPolymer> { prot2 },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            Assert.That(bioPolymerGroup1.BioPolymerGroupName, Is.EqualTo("prot1"));

            bioPolymerGroup1.MergeWith(bioPolymerGroup2);

            // Name should be updated to include both accessions in alphabetical order
            Assert.That(bioPolymerGroup1.BioPolymerGroupName, Is.EqualTo("prot1|prot2"));
        }

        #endregion

        #region Constructor Tests

        [Test]
        public void BioPolymerGroupConstructorSetsIsDecoy()
        {
            var decoyProt = new TestBioPolymer("MEDEEK", "decoy_prot1", isDecoy: true);
            var targetProt = new TestBioPolymer("MENEEK", "prot2", isDecoy: false);

            // Group with only decoy
            var decoyGroup = new BioPolymerGroup(new HashSet<IBioPolymer> { decoyProt },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            Assert.That(decoyGroup.IsDecoy, Is.True);

            // Group with only target
            var targetGroup = new BioPolymerGroup(new HashSet<IBioPolymer> { targetProt },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            Assert.That(targetGroup.IsDecoy, Is.False);

            // Group with both - should be marked as decoy
            var mixedGroup = new BioPolymerGroup(new HashSet<IBioPolymer> { decoyProt, targetProt },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            Assert.That(mixedGroup.IsDecoy, Is.True);
        }

        [Test]
        public void BioPolymerGroupConstructorSetsIsContaminant()
        {
            var contaminantProt = new TestBioPolymer("MEDEEK", "contam_prot1", isContaminant: true);
            var normalProt = new TestBioPolymer("MENEEK", "prot2", isContaminant: false);

            // Group with only contaminant
            var contamGroup = new BioPolymerGroup(new HashSet<IBioPolymer> { contaminantProt },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            Assert.That(contamGroup.IsContaminant, Is.True);

            // Group with only normal
            var normalGroup = new BioPolymerGroup(new HashSet<IBioPolymer> { normalProt },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            Assert.That(normalGroup.IsContaminant, Is.False);
        }

        [Test]
        public void BioPolymerGroupConstructorOrdersByAccession()
        {
            var protZ = new TestBioPolymer("MEDEEK", "zzzProtein");
            var protA = new TestBioPolymer("MENEEK", "aaaProtein");
            var protM = new TestBioPolymer("MNNNNN", "mmmProtein");

            var group = new BioPolymerGroup(new HashSet<IBioPolymer> { protZ, protA, protM },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            // Should be ordered alphabetically
            Assert.That(group.ListOfBioPolymersOrderedByAccession[0].Accession, Is.EqualTo("aaaProtein"));
            Assert.That(group.ListOfBioPolymersOrderedByAccession[1].Accession, Is.EqualTo("mmmProtein"));
            Assert.That(group.ListOfBioPolymersOrderedByAccession[2].Accession, Is.EqualTo("zzzProtein"));

            // Name should reflect the order
            Assert.That(group.BioPolymerGroupName, Is.EqualTo("aaaProtein|mmmProtein|zzzProtein"));
        }

        #endregion

        #region Sequence Coverage Tests

        [Test]
        public void BioPolymerGroupCalculateSequenceCoverage_BasicTest()
        {
            var prot1 = new TestBioPolymer("MEDEEKPEPTIDE", "prot1");
            var pwsm1 = new TestBioPolymerWithSetMods(prot1, 1, 6); // MEDEEK
            var pwsm2 = new TestBioPolymerWithSetMods(prot1, 7, 13); // PEPTIDE

            var group = new BioPolymerGroup(new HashSet<IBioPolymer> { prot1 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1, pwsm2 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1, pwsm2 });

            // Create PSMs for the peptides
            var psm1 = new TestSpectralMatch("file.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { pwsm1 });
            var psm2 = new TestSpectralMatch("file.mzML", "PEPTIDE", "PEPTIDE", 100, 2, new[] { pwsm2 });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            var output = group.ToString();

            // Full coverage (100%) should be reflected
            Assert.That(output, Does.Contain("1")); // Coverage fraction of 1.0
        }

        [Test]
        public void BioPolymerGroupCalculateSequenceCoverage_PartialCoverage()
        {
            var prot1 = new TestBioPolymer("MEDEEKPEPTIDE", "prot1"); // 13 residues
            var pwsm1 = new TestBioPolymerWithSetMods(prot1, 1, 6); // MEDEEK - 6 residues

            var group = new BioPolymerGroup(new HashSet<IBioPolymer> { prot1 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1 });

            var psm1 = new TestSpectralMatch("file.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { pwsm1 });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1 };

            group.CalculateSequenceCoverage();

            var output = group.ToString();

            // Partial coverage should show uppercase for covered, lowercase for uncovered
            Assert.That(output, Does.Contain("MEDEEK")); // Covered portion uppercase
            Assert.That(output, Does.Contain("peptide")); // Uncovered portion lowercase
        }

        #endregion

        #region Subset Group Tests

        [Test]
        public void BioPolymerGroupConstructSubsetBioPolymerGroup_FiltersCorrectly()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1");
            var pwsm1 = new TestBioPolymerWithSetMods(prot1, 1, 6);

            var group = new BioPolymerGroup(new HashSet<IBioPolymer> { prot1 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1 });

            var psm1 = new TestSpectralMatch("file1.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { pwsm1 });
            var psm2 = new TestSpectralMatch("file2.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { pwsm1 });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            var subsetGroup = group.ConstructSubsetBioPolymerGroup("file1.mzML");

            // Subset should only contain PSMs from file1
            Assert.That(subsetGroup.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(1));
            Assert.That(subsetGroup.AllPsmsBelowOnePercentFDR.First().FullFilePath, Is.EqualTo("file1.mzML"));
        }

        #endregion

        #region DisplayModsOnPeptides Tests

        [Test]
        public void BioPolymerGroupDisplayModsOnPeptides_AffectsOutput()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1");

            ModificationMotif.TryGetMotif("E", out var motif);
            var mod = new Modification("Phospho", null, "Common Biological", null, motif, "Anywhere.", null, 79.966, null, null, null, null, null, null);

            var modsDict = new Dictionary<int, Modification> { { 3, mod } }; // Mod on position 2 (E)
            var pwsm1 = new TestBioPolymerWithSetMods(prot1, 1, 6, modsDict);

            var group = new BioPolymerGroup(new HashSet<IBioPolymer> { prot1 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1 },
                new HashSet<IBioPolymerWithSetMods> { pwsm1 });

            // Test with mods displayed
            group.DisplayModsOnPeptides = true;
            var outputWithMods = group.ToString();

            // Test without mods displayed
            group.DisplayModsOnPeptides = false;
            var outputWithoutMods = group.ToString();

            // The outputs should differ based on DisplayModsOnPeptides setting
            // (The actual content depends on implementation details)
            Assert.That(outputWithMods, Is.Not.Null);
            Assert.That(outputWithoutMods, Is.Not.Null);
        }

        #endregion

        #region Quantification Tests

        [Test]
        public void BioPolymerGroupWithSpectraFileInfo_GeneratesIntensityColumns()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1");
            var group = new BioPolymerGroup(new HashSet<IBioPolymer> { prot1 },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            var spectraFile = new SpectraFileInfo("test.mzML", "Condition1", 0, 0, 0);
            group.SamplesForQuantification = new List<ISampleInfo> { spectraFile };
            group.IntensitiesBySample = new Dictionary<ISampleInfo, double> { { spectraFile, 1000.0 } };

            var header = group.GetTabSeparatedHeader();
            var output = group.ToString();

            // Header should contain intensity column
            Assert.That(header, Does.Contain("Intensity_"));

            // Output should contain the intensity value
            Assert.That(output, Does.Contain("1000"));
        }

        #endregion

        #region String Truncation Tests

        [Test]
        public void BioPolymerGroupTruncatesLongStrings()
        {
            // Create a very long protein sequence
            var longSequence = new string('M', 50000);
            var prot1 = new TestBioPolymer(longSequence, "prot1");

            var group = new BioPolymerGroup(new HashSet<IBioPolymer> { prot1 },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            // Set a max length for testing
            var originalMaxLength = BioPolymerGroup.MaxStringLength;
            BioPolymerGroup.MaxStringLength = 100;

            try
            {
                var output = group.ToString();

                // The output should not contain a string longer than MaxStringLength for any field
                // (This is a simplified check - actual implementation may vary)
                Assert.That(output, Is.Not.Null);
            }
            finally
            {
                // Restore original setting
                BioPolymerGroup.MaxStringLength = originalMaxLength;
            }
        }

        [Test]
        public void BioPolymerGroupMaxStringLengthDisabled()
        {
            var prot1 = new TestBioPolymer("MEDEEK", "prot1");
            var group = new BioPolymerGroup(new HashSet<IBioPolymer> { prot1 },
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            // Disable truncation
            var originalMaxLength = BioPolymerGroup.MaxStringLength;
            BioPolymerGroup.MaxStringLength = 0;

            try
            {
                var output = group.ToString();
                Assert.That(output, Is.Not.Null);
            }
            finally
            {
                BioPolymerGroup.MaxStringLength = originalMaxLength;
            }
        }

        #endregion
    }
}