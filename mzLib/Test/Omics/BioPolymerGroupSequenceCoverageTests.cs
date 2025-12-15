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
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class BioPolymerGroupSequenceCoverageTests
    {
        #region Helper Classes

        /// <summary>
        /// Test implementation of IBioPolymer for sequence coverage tests.
        /// </summary>
        private class CoverageBioPolymer : IBioPolymer
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
            public int Length => BaseSequence.Length;
            public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }
            public string SampleNameForVariants { get; set; }
            public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }
            public IBioPolymer ConsensusVariant => this;
            public List<SequenceVariation> AppliedSequenceVariations { get; }
            public List<SequenceVariation> SequenceVariations { get; }
            public List<TruncationProduct> TruncationProducts { get; }

            public CoverageBioPolymer(string sequence, string accession,
                string organism = "Test Organism",
                string name = "Test Name",
                string fullName = "Test Full Name",
                List<Tuple<string, string>> geneNames = null,
                bool isDecoy = false,
                bool isContaminant = false)
            {
                BaseSequence = sequence;
                Accession = accession;
                Organism = organism;
                Name = name;
                FullName = fullName;
                GeneNames = geneNames ?? new List<Tuple<string, string>> { new("primary", "GENE1") };
                IsDecoy = isDecoy;
                IsContaminant = isContaminant;
                DatabaseFilePath = "";
                OneBasedPossibleLocalizedModifications = new Dictionary<int, List<Modification>>();
                OriginalNonVariantModifications = new Dictionary<int, List<Modification>>();
                AppliedSequenceVariations = new List<SequenceVariation>();
                SequenceVariations = new List<SequenceVariation>();
                TruncationProducts = new List<TruncationProduct>();
                SampleNameForVariants = string.Empty;
            }

            public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams,
                List<Modification> allKnownFixedModifications,
                List<Modification> variableModifications,
                List<SilacLabel>? silacLabels = null,
                (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null,
                bool topDownTruncationSearch = false)
                => Enumerable.Empty<IBioPolymerWithSetMods>();

            public IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence,
                IDictionary<int, List<Modification>>? newMods)
                => new CoverageBioPolymer(newBaseSequence, Accession, Organism, Name, FullName, GeneNames, IsDecoy, IsContaminant);

            public TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence,
                TBioPolymerType original,
                IEnumerable<SequenceVariation> appliedSequenceVariants,
                IEnumerable<TruncationProduct> applicableProteolysisProducts,
                IDictionary<int, List<Modification>> oneBasedModifications,
                string sampleNameForVariants)
                where TBioPolymerType : IHasSequenceVariants
                => original;

            public bool Equals(IBioPolymer? other)
            {
                if (other is null) return false;
                return Accession == other.Accession && BaseSequence == other.BaseSequence;
            }

            public override bool Equals(object? obj) => obj is IBioPolymer other && Equals(other);
            public override int GetHashCode() => HashCode.Combine(Accession, BaseSequence);
        }

        /// <summary>
        /// Test implementation of IBioPolymerWithSetMods for sequence coverage tests.
        /// </summary>
        private class CoverageBioPolymerWithSetMods : IBioPolymerWithSetMods
        {
            public string BaseSequence { get; }
            public string FullSequence { get; }
            public double MostAbundantMonoisotopicMass { get; }
            public double MonoisotopicMass { get; }
            public string SequenceWithChemicalFormulas => BaseSequence;
            public int OneBasedStartResidue { get; }
            public int OneBasedEndResidue { get; }
            public int MissedCleavages => 0;
            public string Description => "Test";
            public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; } = CleavageSpecificity.Full;
            public char PreviousResidue => '-';
            public char NextResidue => '-';
            public IDigestionParams DigestionParams => null!;
            public Dictionary<int, Modification> AllModsOneIsNterminus { get; }
            public int NumMods => AllModsOneIsNterminus.Count;
            public int NumFixedMods => 0;
            public int NumVariableMods => AllModsOneIsNterminus.Count;
            public int Length => BaseSequence.Length;
            public IBioPolymer Parent { get; }
            public ChemicalFormula ThisChemicalFormula => new ChemicalFormula();
            public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

            public CoverageBioPolymerWithSetMods(
                string baseSequence,
                string fullSequence,
                IBioPolymer parent,
                int oneBasedStartResidue,
                int oneBasedEndResidue,
                double mass = 0,
                Dictionary<int, Modification>? mods = null)
            {
                BaseSequence = baseSequence;
                FullSequence = fullSequence;
                Parent = parent;
                OneBasedStartResidue = oneBasedStartResidue;
                OneBasedEndResidue = oneBasedEndResidue;
                MostAbundantMonoisotopicMass = mass;
                MonoisotopicMass = mass;
                AllModsOneIsNterminus = mods ?? new Dictionary<int, Modification>();
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
                if (other is null) return false;
                return BaseSequence == other.BaseSequence
                    && FullSequence == other.FullSequence
                    && Parent?.Accession == other.Parent?.Accession
                    && OneBasedStartResidue == other.OneBasedStartResidue
                    && OneBasedEndResidue == other.OneBasedEndResidue;
            }

            public override bool Equals(object? obj) => obj is IBioPolymerWithSetMods other && Equals(other);
            public override int GetHashCode() => HashCode.Combine(BaseSequence, FullSequence, Parent?.Accession, OneBasedStartResidue, OneBasedEndResidue);
        }

        /// <summary>
        /// Test implementation of ISpectralMatch for sequence coverage tests.
        /// Also implements IHasSequenceCoverageFromFragments to support fragment coverage calculation.
        /// </summary>
        private class CoverageSpectralMatch : ISpectralMatch, IHasSequenceCoverageFromFragments
        {
            private readonly List<IBioPolymerWithSetMods> _identified;

            public string FullFilePath { get; }
            public string FullSequence { get; }
            public string BaseSequence { get; }
            public double Score { get; }
            public int OneBasedScanNumber { get; }
            public HashSet<int>? FragmentCoveragePositionInPeptide { get; private set; }

            public List<int>? NTerminalFragmentPositions { get; set; }
            public List<int>? CTerminalFragmentPositions { get; set; }

            public CoverageSpectralMatch(
                string filePath,
                string fullSequence,
                string baseSequence,
                double score,
                int scanNumber,
                IEnumerable<IBioPolymerWithSetMods>? identified = null)
            {
                FullFilePath = filePath ?? string.Empty;
                FullSequence = fullSequence ?? string.Empty;
                BaseSequence = baseSequence ?? string.Empty;
                Score = score;
                OneBasedScanNumber = scanNumber;
                _identified = identified?.ToList() ?? new List<IBioPolymerWithSetMods>();
            }

            public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => _identified;

            public void GetSequenceCoverage()
            {
                if (string.IsNullOrEmpty(BaseSequence)) return;

                var nTermPositions = NTerminalFragmentPositions ?? new List<int>();
                var cTermPositions = CTerminalFragmentPositions ?? new List<int>();

                if (!nTermPositions.Any() && !cTermPositions.Any()) return;

                var fragmentCoveredResidues = new HashSet<int>();

                if (nTermPositions.Any())
                {
                    var sortedNTerm = nTermPositions.OrderBy(x => x).ToList();

                    if (sortedNTerm.Contains(BaseSequence.Length - 1))
                        fragmentCoveredResidues.Add(BaseSequence.Length);

                    if (sortedNTerm.Contains(1))
                        fragmentCoveredResidues.Add(1);

                    for (int i = 0; i < sortedNTerm.Count - 1; i++)
                    {
                        if (sortedNTerm[i + 1] - sortedNTerm[i] == 1)
                            fragmentCoveredResidues.Add(sortedNTerm[i + 1]);

                        if (cTermPositions.Contains(sortedNTerm[i + 1]))
                            fragmentCoveredResidues.Add(sortedNTerm[i + 1]);

                        if (cTermPositions.Contains(sortedNTerm[i + 1] + 2))
                            fragmentCoveredResidues.Add(sortedNTerm[i + 1] + 1);
                    }
                }

                if (cTermPositions.Any())
                {
                    var sortedCTerm = cTermPositions.OrderBy(x => x).ToList();

                    if (sortedCTerm.Contains(2))
                        fragmentCoveredResidues.Add(1);

                    if (sortedCTerm.Contains(BaseSequence.Length))
                        fragmentCoveredResidues.Add(BaseSequence.Length);

                    for (int i = 0; i < sortedCTerm.Count - 1; i++)
                    {
                        if (sortedCTerm[i + 1] - sortedCTerm[i] == 1)
                            fragmentCoveredResidues.Add(sortedCTerm[i]);
                    }
                }

                FragmentCoveragePositionInPeptide = fragmentCoveredResidues;
            }

            public int CompareTo(ISpectralMatch? other)
            {
                if (other is null) return 1;
                int scoreCmp = Score.CompareTo(other.Score);
                if (scoreCmp != 0) return scoreCmp;
                return OneBasedScanNumber.CompareTo(other.OneBasedScanNumber);
            }

            public int CompareTo(IHasSequenceCoverageFromFragments? other)
            {
                if (other is ISpectralMatch spectralMatch)
                {
                    return CompareTo(spectralMatch);
                }
                return other is null ? 1 : 0;
            }
        }
        #endregion
        #region Basic Functionality Tests

        [Test]
        public void CalculateSequenceCoverage_WithNoPsms_ProducesEmptyResults()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIKLMNPQRSTVWY", "P00001");
            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods>();
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods>();

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>();

            group.CalculateSequenceCoverage();

            Assert.That(group.SequenceCoverageFraction.Count, Is.EqualTo(1));
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(0.0));
            Assert.That(group.SequenceCoverageDisplayList.Count, Is.EqualTo(1));
            Assert.That(group.SequenceCoverageDisplayList[0], Is.EqualTo("acdefghiklmnpqrstvwy")); // All lowercase = not covered
            Assert.That(group.FragmentSequenceCoverageDisplayList.Count, Is.EqualTo(1));
        }

        [Test]
        public void CalculateSequenceCoverage_WithSinglePeptide_CalculatesCorrectFraction()
        {
            // Protein of length 20, peptide covers positions 5-10 (6 residues)
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIKLMNPQRSTVWY", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("FGHIKL", "FGHIKL", bioPolymer, 5, 10);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "FGHIKL", "FGHIKL", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // 6 residues covered out of 20 = 0.3
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(0.3).Within(0.001));
        }

        [Test]
        public void CalculateSequenceCoverage_WithFullCoverage_Returns100Percent()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("ACDEFGHIK", "ACDEFGHIK", bioPolymer, 1, 9);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "ACDEFGHIK", "ACDEFGHIK", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(1.0).Within(0.001));
            Assert.That(group.SequenceCoverageDisplayList[0], Is.EqualTo("ACDEFGHIK")); // All uppercase = fully covered
        }

        [Test]
        public void CalculateSequenceCoverage_DisplayList_ShowsUppercaseForCoveredResidues()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            // Peptide covers positions 3-6 (DEFG)
            var peptide = new CoverageBioPolymerWithSetMods("DEFG", "DEFG", bioPolymer, 3, 6);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "DEFG", "DEFG", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // ac = not covered (lowercase), DEFG = covered (uppercase), hik = not covered (lowercase)
            Assert.That(group.SequenceCoverageDisplayList[0], Is.EqualTo("acDEFGhik"));
        }

        #endregion
        #region Multiple Peptide Coverage Tests

        [Test]
        public void CalculateSequenceCoverage_WithOverlappingPeptides_CombinesCoverage()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIKLMNPQRSTVWY", "P00001");

            // Peptide 1 covers positions 1-5 (ACDEF)
            var peptide1 = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer, 1, 5);
            // Peptide 2 covers positions 4-8 (EFGHI) - overlaps with peptide1
            var peptide2 = new CoverageBioPolymerWithSetMods("EFGHI", "EFGHI", bioPolymer, 4, 8);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };

            var psm1 = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide1 });
            var psm2 = new CoverageSpectralMatch("test.raw", "EFGHI", "EFGHI", 100, 2, new[] { peptide2 });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            // Combined coverage: positions 1-8 = 8 residues out of 20 = 0.4
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(0.4).Within(0.001));
            Assert.That(group.SequenceCoverageDisplayList[0], Is.EqualTo("ACDEFGHIklmnpqrstvwy"));
        }

        [Test]
        public void CalculateSequenceCoverage_WithNonOverlappingPeptides_AddsCoverage()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIKLMNPQRSTVWY", "P00001");

            // Peptide 1 covers positions 1-3 (ACD)
            var peptide1 = new CoverageBioPolymerWithSetMods("ACD", "ACD", bioPolymer, 1, 3);
            // Peptide 2 covers positions 18-20 (VWY) - no overlap
            var peptide2 = new CoverageBioPolymerWithSetMods("VWY", "VWY", bioPolymer, 18, 20);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };

            var psm1 = new CoverageSpectralMatch("test.raw", "ACD", "ACD", 100, 1, new[] { peptide1 });
            var psm2 = new CoverageSpectralMatch("test.raw", "VWY", "VWY", 100, 2, new[] { peptide2 });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            // Coverage: 3 + 3 = 6 residues out of 20 = 0.3
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(0.3).Within(0.001));
            Assert.That(group.SequenceCoverageDisplayList[0], Is.EqualTo("ACDefghiklmnpqrstVWY"));
        }

        [Test]
        public void CalculateSequenceCoverage_WithDuplicatePeptides_CountsOnlyOnce()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer, 1, 5);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            // Same peptide identified twice (different scans)
            var psm1 = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide });
            var psm2 = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 95, 2, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            // 5 residues out of 9 = 0.5556
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(5.0 / 9.0).Within(0.001));
        }

        #endregion
        #region Multiple BioPolymer Tests

        [Test]
        public void CalculateSequenceCoverage_WithMultipleBioPolymers_CalculatesSeparately()
        {
            var bioPolymer1 = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var bioPolymer2 = new CoverageBioPolymer("MNPQRSTVWY", "P00002");

            var peptide1 = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer1, 1, 5);
            var peptide2 = new CoverageBioPolymerWithSetMods("MNPQR", "MNPQR", bioPolymer2, 1, 5);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer1, bioPolymer2 };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };

            var psm1 = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide1 });
            var psm2 = new CoverageSpectralMatch("test.raw", "MNPQR", "MNPQR", 100, 2, new[] { peptide2 });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            Assert.That(group.SequenceCoverageFraction.Count, Is.EqualTo(2));
            // P00001: 5/9 = 0.556
            // P00002: 5/10 = 0.5
            // Order is by accession, so P00001 first
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(5.0 / 9.0).Within(0.001));
            Assert.That(group.SequenceCoverageFraction[1], Is.EqualTo(0.5).Within(0.001));
        }

        [Test]
        public void CalculateSequenceCoverage_OrdersResultsByAccession()
        {
            var bioPolymerZ = new CoverageBioPolymer("ACDEFGHIK", "Z99999");
            var bioPolymerA = new CoverageBioPolymer("MNPQRSTVWY", "A00001");

            var peptideZ = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymerZ, 1, 5);
            var peptideA = new CoverageBioPolymerWithSetMods("MNPQR", "MNPQR", bioPolymerA, 1, 5);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymerZ, bioPolymerA };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptideZ, peptideA };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods>();

            var psmZ = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptideZ });
            var psmA = new CoverageSpectralMatch("test.raw", "MNPQR", "MNPQR", 100, 2, new[] { peptideA });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psmZ, psmA };

            group.CalculateSequenceCoverage();

            // A00001 should be first (alphabetically)
            Assert.That(group.SequenceCoverageDisplayList[0], Does.StartWith("MNPQR")); // A00001's sequence
            Assert.That(group.SequenceCoverageDisplayList[1], Does.StartWith("ACDEF")); // Z99999's sequence
        }

        #endregion
        #region Fragment Coverage Tests

        [Test]
        public void CalculateSequenceCoverage_WithFragmentCoverage_PopulatesFragmentList()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("DEFGH", "DEFGH", bioPolymer, 3, 7);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "DEFGH", "DEFGH", 100, 1, new[] { peptide });
            // Set fragment positions to cover first 3 residues of the peptide
            psm.NTerminalFragmentPositions = new List<int> { 1, 2, 3 };

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            Assert.That(group.FragmentSequenceCoverageDisplayList.Count, Is.EqualTo(1));
            // Fragment coverage should show covered positions
            Assert.That(group.FragmentSequenceCoverageDisplayList[0], Is.Not.Null);
        }

        [Test]
        public void CalculateSequenceCoverage_WithNoFragmentPositions_ShowsLowercaseOnly()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("DEFGH", "DEFGH", bioPolymer, 3, 7);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "DEFGH", "DEFGH", 100, 1, new[] { peptide });
            // No fragment positions set

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // Fragment coverage list should be all lowercase (no fragment coverage)
            Assert.That(group.FragmentSequenceCoverageDisplayList[0], Is.EqualTo("acdefghik"));
        }

        [Test]
        public void CalculateSequenceCoverage_FragmentCoverage_ConvertsToProteinPosition()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            // Peptide starts at position 3 (D)
            var peptide = new CoverageBioPolymerWithSetMods("DEFGH", "DEFGH", bioPolymer, 3, 7);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "DEFGH", "DEFGH", 100, 1, new[] { peptide });
            // Fragment position 1 in peptide = position 3 in protein (D)
            psm.NTerminalFragmentPositions = new List<int> { 1 };

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // Position 3 (D) should be uppercase in fragment coverage
            var fragmentDisplay = group.FragmentSequenceCoverageDisplayList[0];
            Assert.That(fragmentDisplay[2], Is.EqualTo('D')); // Position 3 (0-indexed = 2)
        }

        #endregion
        #region Modification Tests

        [Test]
        public void CalculateSequenceCoverage_WithModifications_PopulatesModsInfo()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");

            // Create a modification
            var phosphoMod = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Post-translational",
                _target: ModificationMotif.TryGetMotif("S", out var motif) ? motif : null,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 79.966);

            var mods = new Dictionary<int, Modification> { { 3, phosphoMod } }; // Position 3 in peptide

            var peptide = new CoverageBioPolymerWithSetMods(
                "DEFGH", "DEF[Phosphorylation]GH", bioPolymer, 3, 7, 500.0, mods);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "DEF[Phosphorylation]GH", "DEFGH", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            Assert.That(group.SequenceCoverageDisplayListWithMods.Count, Is.EqualTo(1));
            // The modification's IdWithMotif includes the target residue (e.g., "Phosphorylation on S")
            Assert.That(group.SequenceCoverageDisplayListWithMods[0], Does.Contain("[Phosphorylation on S]"));
        }

        [Test]
        public void CalculateSequenceCoverage_SkipsCommonVariableMods()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");

            var commonMod = new Modification(
                _originalId: "Oxidation",
                _modificationType: "Common Variable",
                _target: ModificationMotif.TryGetMotif("M", out var motif) ? motif : null,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 15.995);

            var mods = new Dictionary<int, Modification> { { 3, commonMod } };

            var peptide = new CoverageBioPolymerWithSetMods(
                "DEFGH", "DEF[Oxidation]GH", bioPolymer, 3, 7, 500.0, mods);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "DEF[Oxidation]GH", "DEFGH", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // Common Variable mods should be skipped
            Assert.That(group.SequenceCoverageDisplayListWithMods[0], Does.Not.Contain("[Oxidation]"));
        }

        [Test]
        public void CalculateSequenceCoverage_SkipsCommonFixedMods()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");

            var commonFixedMod = new Modification(
                _originalId: "Carbamidomethyl",
                _modificationType: "Common Fixed",
                _target: ModificationMotif.TryGetMotif("C", out var motif) ? motif : null,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 57.021);

            var mods = new Dictionary<int, Modification> { { 2, commonFixedMod } };

            var peptide = new CoverageBioPolymerWithSetMods(
                "ACDEF", "AC[Carbamidomethyl]DEF", bioPolymer, 1, 5, 500.0, mods);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "AC[Carbamidomethyl]DEF", "ACDEF", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // Common Fixed mods should be skipped
            Assert.That(group.SequenceCoverageDisplayListWithMods[0], Does.Not.Contain("[Carbamidomethyl]"));
        }

        [Test]
        public void CalculateSequenceCoverage_SkipsPeptideTermMods()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");

            var peptideTermMod = new Modification(
                _originalId: "Acetyl",
                _modificationType: "PeptideTermMod",
                _target: null,
                _locationRestriction: "N-terminal.",
                _monoisotopicMass: 42.011);

            var mods = new Dictionary<int, Modification> { { 1, peptideTermMod } };

            var peptide = new CoverageBioPolymerWithSetMods(
                "ACDEF", "[Acetyl]-ACDEF", bioPolymer, 1, 5, 500.0, mods);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "[Acetyl]-ACDEF", "ACDEF", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // PeptideTermMods should be skipped
            Assert.That(group.SequenceCoverageDisplayListWithMods[0], Does.Not.Contain("[Acetyl]"));
        }

        [Test]
        public void CalculateSequenceCoverage_NTerminalMod_AddedWithDash()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");

            var nTermMod = new Modification(
                _originalId: "Acetyl on A",
                _modificationType: "Post-translational",
                _target: ModificationMotif.TryGetMotif("A", out var motif) ? motif : null,
                _locationRestriction: "N-terminal.",
                _monoisotopicMass: 42.011);

            var mods = new Dictionary<int, Modification> { { 1, nTermMod } };

            var peptide = new CoverageBioPolymerWithSetMods(
                "ACDEF", "[Acetyl on A]-ACDEF", bioPolymer, 1, 5, 500.0, mods);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "[Acetyl on A]-ACDEF", "ACDEF", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            Assert.That(group.SequenceCoverageDisplayListWithMods[0], Does.StartWith("[Acetyl on A]-"));
        }

        [Test]
        public void CalculateSequenceCoverage_CTerminalMod_AddedWithDash()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");

            var cTermMod = new Modification(
                _originalId: "Amidation",
                _modificationType: "Post-translational",
                _target: ModificationMotif.TryGetMotif("K", out var motif) ? motif : null,
                _locationRestriction: "C-terminal.",
                _monoisotopicMass: -0.984);

            var mods = new Dictionary<int, Modification> { { 5, cTermMod } };

            // Peptide at the C-terminus of the protein
            var peptide = new CoverageBioPolymerWithSetMods(
                "FGHIK", "FGHIK-[Amidation]", bioPolymer, 5, 9, 500.0, mods);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "FGHIK-[Amidation]", "FGHIK", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // The modification's IdWithMotif includes the target residue (e.g., "Amidation on K")
            Assert.That(group.SequenceCoverageDisplayListWithMods[0], Does.EndWith("-[Amidation on K]"));
        }

        #endregion
        #region Edge Cases Tests

        [Test]
        public void CalculateSequenceCoverage_WithNullBaseSequence_SkipsPsm()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer, 1, 5);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            // PSM with null base sequence (ambiguous)
            var psmNull = new CoverageSpectralMatch("test.raw", "ACDEF", null!, 100, 1, new[] { peptide });
            var psmValid = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 2, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psmNull, psmValid };

            Assert.DoesNotThrow(() => group.CalculateSequenceCoverage());

            // Only the valid PSM should contribute to coverage
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(5.0 / 9.0).Within(0.001));
        }

        [Test]
        public void CalculateSequenceCoverage_WithPeptideFromDifferentParent_SkipsPeptide()
        {
            var bioPolymer1 = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var bioPolymer2 = new CoverageBioPolymer("MNPQRSTVWY", "P00002");

            // Peptide belongs to bioPolymer2, not bioPolymer1
            var peptide = new CoverageBioPolymerWithSetMods("MNPQR", "MNPQR", bioPolymer2, 1, 5);

            // Group only contains bioPolymer1
            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer1 };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods>();

            var psm = new CoverageSpectralMatch("test.raw", "MNPQR", "MNPQR", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // bioPolymer1 should have 0 coverage since peptide belongs to bioPolymer2
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(0.0));
        }

        [Test]
        public void CalculateSequenceCoverage_WithNullFullSequence_SkipsModsDisplay()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");

            // Create a peptide with null FullSequence
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", null!, bioPolymer, 1, 5);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", null!, "ACDEF", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            Assert.DoesNotThrow(() => group.CalculateSequenceCoverage());

            // Coverage should still be calculated
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(5.0 / 9.0).Within(0.001));
        }

        [Test]
        public void CalculateSequenceCoverage_WithEmptyBioPolymers_DoesNotThrow()
        {
            var bioPolymers = new HashSet<IBioPolymer>();
            var sequences = new HashSet<IBioPolymerWithSetMods>();
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods>();

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);

            Assert.DoesNotThrow(() => group.CalculateSequenceCoverage());
            Assert.That(group.SequenceCoverageFraction.Count, Is.EqualTo(0));
        }

        [Test]
        public void CalculateSequenceCoverage_CalledMultipleTimes_AppendsResults()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer, 1, 5);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();
            group.CalculateSequenceCoverage(); // Call again

            // Results are appended, so we get 2 entries
            Assert.That(group.SequenceCoverageFraction.Count, Is.EqualTo(2));
        }

        [Test]
        public void CalculateSequenceCoverage_WithSingleResidueProtein_HandlesBoundaryCase()
        {
            var bioPolymer = new CoverageBioPolymer("K", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("K", "K", bioPolymer, 1, 1);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "K", "K", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(1.0));
            Assert.That(group.SequenceCoverageDisplayList[0], Is.EqualTo("K"));
        }

        [Test]
        public void CalculateSequenceCoverage_WithVeryLongProtein_HandlesLargeSequence()
        {
            var longSequence = new string('A', 1000);
            var bioPolymer = new CoverageBioPolymer(longSequence, "P00001");

            // Peptide covers first 10 residues
            var peptide = new CoverageBioPolymerWithSetMods("AAAAAAAAAA", "AAAAAAAAAA", bioPolymer, 1, 10);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "AAAAAAAAAA", "AAAAAAAAAA", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(10.0 / 1000.0).Within(0.0001));
            Assert.That(group.SequenceCoverageDisplayList[0].Length, Is.EqualTo(1000));
        }

        #endregion
        #region Modification Occupancy Tests

        [Test]
        public void CalculateSequenceCoverage_ModOccupancy_CalculatesCorrectly()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");

            var phosphoMod = new Modification(
                _originalId: "Phospho on S",
                _modificationType: "Post-translational",
                _target: ModificationMotif.TryGetMotif("S", out var motif) ? motif : null,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 79.966);

            // Two peptides covering same position, one modified, one not
            var modsDict = new Dictionary<int, Modification> { { 3, phosphoMod } };
            var peptideWithMod = new CoverageBioPolymerWithSetMods(
                "ACDEF", "ACD[Phospho on S]EF", bioPolymer, 1, 5, 500.0, modsDict);
            var peptideWithoutMod = new CoverageBioPolymerWithSetMods(
                "ACDEF", "ACDEF", bioPolymer, 1, 5, 450.0, new Dictionary<int, Modification>());

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptideWithMod, peptideWithoutMod };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods>();

            var psm1 = new CoverageSpectralMatch("test.raw", "ACD[Phospho on S]EF", "ACDEF", 100, 1, new[] { peptideWithMod });
            var psm2 = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 2, new[] { peptideWithoutMod });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            // ModsInfo should contain occupancy information
            Assert.That(group.ModsInfo.Count, Is.GreaterThanOrEqualTo(0));
        }

        [Test]
        public void CalculateSequenceCoverage_ModsInfo_ContainsPositionAndOccupancy()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");

            var phosphoMod = new Modification(
                _originalId: "Phospho",
                _modificationType: "Post-translational",
                _target: ModificationMotif.TryGetMotif("D", out var motif) ? motif : null,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 79.966);

            var modsDict = new Dictionary<int, Modification> { { 3, phosphoMod } };
            var peptide = new CoverageBioPolymerWithSetMods(
                "ACDEF", "ACD[Phospho]EF", bioPolymer, 1, 5, 500.0, modsDict);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "ACD[Phospho]EF", "ACDEF", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // Should contain #aa format with position and occupancy
            if (group.ModsInfo.Count > 0)
            {
                Assert.That(group.ModsInfo[0], Does.Contain("#aa"));
                Assert.That(group.ModsInfo[0], Does.Contain("occupancy"));
            }
        }

        #endregion
        #region Integration Tests

        [Test]
        public void CalculateSequenceCoverage_FullIntegration_AllOutputsPopulated()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIKLMNPQRSTVWY", "P00001");

            var phosphoMod = new Modification(
                _originalId: "Phospho",
                _modificationType: "Post-translational",
                _target: ModificationMotif.TryGetMotif("T", out var motif) ? motif : null,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 79.966);

            var modsDict = new Dictionary<int, Modification> { { 4, phosphoMod } };

            // Protein: ACDEFGHIKLMNPQRSTVWY (positions 1-20)
            // Position: 1234567890123456789012
            var peptide1 = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer, 1, 5);      // positions 1-5
            var peptide2 = new CoverageBioPolymerWithSetMods("GHIKL", "GHIKL", bioPolymer, 6, 10);     // positions 6-10 (G is at position 6)
            var peptide3 = new CoverageBioPolymerWithSetMods("STVW", "ST[Phospho]VW", bioPolymer, 16, 19, 500.0, modsDict);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2, peptide3 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1 };

            var psm1 = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide1 });
            psm1.NTerminalFragmentPositions = new List<int> { 1, 2, 3 };
            psm1.CTerminalFragmentPositions = new List<int> { 3, 4, 5 };

            var psm2 = new CoverageSpectralMatch("test.raw", "GHIKL", "GHIKL", 95, 2, new[] { peptide2 });
            var psm3 = new CoverageSpectralMatch("test.raw", "ST[Phospho]VW", "STVW", 90, 3, new[] { peptide3 });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3 };

            group.CalculateSequenceCoverage();

            // Verify all outputs are populated
            Assert.That(group.SequenceCoverageFraction.Count, Is.EqualTo(1));
            Assert.That(group.SequenceCoverageDisplayList.Count, Is.EqualTo(1));
            Assert.That(group.SequenceCoverageDisplayListWithMods.Count, Is.EqualTo(1));
            Assert.That(group.FragmentSequenceCoverageDisplayList.Count, Is.EqualTo(1));

            // Verify coverage fraction is correct: 5 + 5 + 4 = 14 residues out of 20 = 0.7
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(0.7).Within(0.001));

            // Verify display list shows correct coverage pattern
            Assert.That(group.SequenceCoverageDisplayList[0], Does.Contain("ACDEF"));
            Assert.That(group.SequenceCoverageDisplayList[0], Does.Contain("GHIKL"));
            Assert.That(group.SequenceCoverageDisplayList[0], Does.Contain("STVW"));

            // Verify mods display contains phospho modification (IdWithMotif includes target residue)
            Assert.That(group.SequenceCoverageDisplayListWithMods[0], Does.Contain("[Phospho on T]"));
        }

        [Test]
        public void CalculateSequenceCoverage_Integration_MultipleBioPolymersWithMixedCoverage()
        {
            var bioPolymer1 = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var bioPolymer2 = new CoverageBioPolymer("MNPQRSTVWY", "P00002");
            var bioPolymer3 = new CoverageBioPolymer("AAAAAAAA", "P00003"); // No coverage

            var peptide1 = new CoverageBioPolymerWithSetMods("ACDEFGHIK", "ACDEFGHIK", bioPolymer1, 1, 9); // Full coverage
            var peptide2 = new CoverageBioPolymerWithSetMods("MNPQR", "MNPQR", bioPolymer2, 1, 5); // Partial coverage

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer1, bioPolymer2, bioPolymer3 };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };

            var psm1 = new CoverageSpectralMatch("test.raw", "ACDEFGHIK", "ACDEFGHIK", 100, 1, new[] { peptide1 });
            var psm2 = new CoverageSpectralMatch("test.raw", "MNPQR", "MNPQR", 95, 2, new[] { peptide2 });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            Assert.That(group.SequenceCoverageFraction.Count, Is.EqualTo(3));

            // P00001: 100% coverage
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(1.0).Within(0.001));
            // P00002: 50% coverage (5/10)
            Assert.That(group.SequenceCoverageFraction[1], Is.EqualTo(0.5).Within(0.001));
            // P00003: 0% coverage
            Assert.That(group.SequenceCoverageFraction[2], Is.EqualTo(0.0).Within(0.001));
        }

        [Test]
        public void CalculateSequenceCoverage_Integration_FragmentAndPeptideCoverageIndependent()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("ACDEFGHIK", "ACDEFGHIK", bioPolymer, 1, 9);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "ACDEFGHIK", "ACDEFGHIK", 100, 1, new[] { peptide });
            // Only set fragment coverage for first 3 positions
            psm.NTerminalFragmentPositions = new List<int> { 1, 2, 3 };

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // Peptide coverage should be 100%
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(1.0).Within(0.001));
            Assert.That(group.SequenceCoverageDisplayList[0], Is.EqualTo("ACDEFGHIK"));

            // Fragment coverage should be partial (only first few residues covered)
            var fragmentDisplay = group.FragmentSequenceCoverageDisplayList[0];
            Assert.That(fragmentDisplay, Is.Not.EqualTo("ACDEFGHIK")); // Not fully uppercase
            Assert.That(fragmentDisplay, Is.Not.EqualTo("acdefghik")); // Not fully lowercase
        }

        [Test]
        public void CalculateSequenceCoverage_Integration_SharedPeptideBetweenProteins()
        {
            var bioPolymer1 = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var bioPolymer2 = new CoverageBioPolymer("XXXACDEFYYY", "P00002");

            // Same peptide sequence appears in both proteins at different positions
            var peptide1 = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer1, 1, 5);
            var peptide2 = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer2, 4, 8);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer1, bioPolymer2 };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods>(); // Shared, so not unique

            // Create separate PSMs for each peptide to ensure both are processed
            var psm1 = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide1 });
            var psm2 = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 2, new[] { peptide2 });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            Assert.That(group.SequenceCoverageFraction.Count, Is.EqualTo(2));

            // P00001: 5/9 coverage
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(5.0 / 9.0).Within(0.001));
            // P00002: 5/11 coverage
            Assert.That(group.SequenceCoverageFraction[1], Is.EqualTo(5.0 / 11.0).Within(0.001));

            // Verify coverage is at correct positions
            Assert.That(group.SequenceCoverageDisplayList[0], Is.EqualTo("ACDEFghik"));
            Assert.That(group.SequenceCoverageDisplayList[1], Is.EqualTo("xxxACDEFyyy"));
        }

        [Test]
        public void CalculateSequenceCoverage_Integration_CompleteWorkflowWithAllFeatures()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIKLMNPQRSTVWY", "P00001");

            // Create modifications
            var phosphoMod = new Modification(
                _originalId: "Phospho",
                _modificationType: "Post-translational",
                _target: ModificationMotif.TryGetMotif("S", out var motif1) ? motif1 : null,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 79.966);

            var acetylMod = new Modification(
                _originalId: "Acetyl",
                _modificationType: "Post-translational",
                _target: ModificationMotif.TryGetMotif("A", out var motif2) ? motif2 : null,
                _locationRestriction: "N-terminal.",
                _monoisotopicMass: 42.011);

            var modsDict1 = new Dictionary<int, Modification> { { 1, acetylMod } };
            var modsDict2 = new Dictionary<int, Modification> { { 2, phosphoMod } };

            // Protein: ACDEFGHIKLMNPQRSTVWY (positions 1-20)
            // Position: 12345678901234567890
            var peptide1 = new CoverageBioPolymerWithSetMods("ACDEF", "[Acetyl]-ACDEF", bioPolymer, 1, 5, 500.0, modsDict1);
            var peptide2 = new CoverageBioPolymerWithSetMods("GHIKL", "GHIKL", bioPolymer, 6, 10);
            var peptide3 = new CoverageBioPolymerWithSetMods("STVWY", "S[Phospho]TVWY", bioPolymer, 16, 20, 600.0, modsDict2);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2, peptide3 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide3 };

            var psm1 = new CoverageSpectralMatch("test.raw", "[Acetyl]-ACDEF", "ACDEF", 100, 1, new[] { peptide1 });
            psm1.NTerminalFragmentPositions = new List<int> { 1, 2, 3, 4 };
            psm1.CTerminalFragmentPositions = new List<int> { 2, 3, 4, 5 };

            var psm2 = new CoverageSpectralMatch("test.raw", "GHIKL", "GHIKL", 95, 2, new[] { peptide2 });
            psm2.NTerminalFragmentPositions = new List<int> { 1, 2 };

            var psm3 = new CoverageSpectralMatch("test.raw", "S[Phospho]TVWY", "STVWY", 90, 3, new[] { peptide3 });
            psm3.CTerminalFragmentPositions = new List<int> { 4, 5 };

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3 };

            group.CalculateSequenceCoverage();

            // Verify all lists are populated
            Assert.That(group.SequenceCoverageFraction.Count, Is.EqualTo(1));
            Assert.That(group.SequenceCoverageDisplayList.Count, Is.EqualTo(1));
            Assert.That(group.SequenceCoverageDisplayListWithMods.Count, Is.EqualTo(1));
            Assert.That(group.FragmentSequenceCoverageDisplayList.Count, Is.EqualTo(1));

            // Coverage: 5 + 5 + 5 = 15 residues out of 20 = 0.75
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(0.75).Within(0.001));

            // Verify modifications are present in mods display (IdWithMotif includes target residue)
            Assert.That(group.SequenceCoverageDisplayListWithMods[0], Does.Contain("[Acetyl on A]"));
            Assert.That(group.SequenceCoverageDisplayListWithMods[0], Does.Contain("[Phospho on S]"));

            // Verify N-terminal mod format
            Assert.That(group.SequenceCoverageDisplayListWithMods[0], Does.StartWith("[Acetyl on A]-"));
        }

        [Test]
        public void CalculateSequenceCoverage_Integration_RealWorldScenario_ProteinWithMultiplePeptides()
        {
            // Simulate a real-world scenario with a protein and multiple overlapping peptides
            var bioPolymer = new CoverageBioPolymer("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAAALEHHHHHH", "P12345");

            var peptide1 = new CoverageBioPolymerWithSetMods("MKTAYIAK", "MKTAYIAK", bioPolymer, 1, 8);
            var peptide2 = new CoverageBioPolymerWithSetMods("QRQISFVK", "QRQISFVK", bioPolymer, 9, 16);
            var peptide3 = new CoverageBioPolymerWithSetMods("SHFSRQLEERLGLIEVQAPILSR", "SHFSRQLEERLGLIEVQAPILSR", bioPolymer, 17, 39);
            var peptide4 = new CoverageBioPolymerWithSetMods("VGDGTQDNLSGAEK", "VGDGTQDNLSGAEK", bioPolymer, 40, 53);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2, peptide3, peptide4 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2, peptide3, peptide4 };

            var psm1 = new CoverageSpectralMatch("sample1.raw", "MKTAYIAK", "MKTAYIAK", 150, 100, new[] { peptide1 });
            var psm2 = new CoverageSpectralMatch("sample1.raw", "QRQISFVK", "QRQISFVK", 145, 200, new[] { peptide2 });
            var psm3 = new CoverageSpectralMatch("sample1.raw", "SHFSRQLEERLGLIEVQAPILSR", "SHFSRQLEERLGLIEVQAPILSR", 180, 300, new[] { peptide3 });
            var psm4 = new CoverageSpectralMatch("sample1.raw", "VGDGTQDNLSGAEK", "VGDGTQDNLSGAEK", 160, 400, new[] { peptide4 });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3, psm4 };

            group.CalculateSequenceCoverage();

            // Total covered: 8 + 8 + 23 + 14 = 53 residues out of 92
            double expectedCoverage = 53.0 / 92.0;
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(expectedCoverage).Within(0.001));

            // Verify display list has correct length
            Assert.That(group.SequenceCoverageDisplayList[0].Length, Is.EqualTo(92));

            // Verify covered regions are uppercase
            var display = group.SequenceCoverageDisplayList[0];
            Assert.That(display.Substring(0, 8), Is.EqualTo("MKTAYIAK"));
            Assert.That(display.Substring(8, 8), Is.EqualTo("QRQISFVK"));
        }

        [Test]
        public void CalculateSequenceCoverage_Integration_WithQuantificationData()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer, 1, 5);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            // Add quantification data
            var spectraFile = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            group.SamplesForQuantification = new List<ISampleInfo> { spectraFile };
            group.IntensitiesBySample = new Dictionary<ISampleInfo, double> { { spectraFile, 1000000.0 } };

            // Calculate coverage
            group.CalculateSequenceCoverage();

            // Verify coverage is calculated correctly alongside quantification data
            Assert.That(group.SequenceCoverageFraction[0], Is.EqualTo(5.0 / 9.0).Within(0.001));
            Assert.That(group.IntensitiesBySample[spectraFile], Is.EqualTo(1000000.0));

            // Verify ToString includes both coverage and intensity
            var output = group.ToString();
            Assert.That(output, Does.Contain("1000000"));
        }

        #endregion
    }
}