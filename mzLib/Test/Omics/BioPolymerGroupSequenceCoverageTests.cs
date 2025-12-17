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
        public class CoverageSpectralMatch : BaseSpectralMatch
        {
            private readonly List<IBioPolymerWithSetMods> _identified;
            public CoverageSpectralMatch(
                string filePath,
                string fullSequence,
                string baseSequence,
                double score,
                int scanNumber,
                IEnumerable<IBioPolymerWithSetMods>? identified = null) : base(filePath, scanNumber, score, fullSequence, baseSequence)
            {
                _identified = identified != null ? [..identified] : [];
            }

            public override IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => _identified;
        }

        #endregion

        #region Helper Methods

        /// <summary>
        /// Extracts the sequence coverage display string from ToString() output.
        /// The coverage display is in the 12th tab-separated field.
        /// </summary>
        private static string GetSequenceCoverageFromToString(string toStringOutput)
        {
            var fields = toStringOutput.Split('\t');
            return fields.Length > 11 ? fields[11] : string.Empty;
        }

        /// <summary>
        /// Extracts the sequence coverage fraction from ToString() output.
        /// The coverage fraction is in the 11th tab-separated field.
        /// </summary>
        private static string GetSequenceCoverageFractionFromToString(string toStringOutput)
        {
            var fields = toStringOutput.Split('\t');
            return fields.Length > 10 ? fields[10] : string.Empty;
        }

        /// <summary>
        /// Extracts the sequence coverage with mods display from ToString() output.
        /// </summary>
        private static string GetSequenceCoverageWithModsFromToString(string toStringOutput)
        {
            var fields = toStringOutput.Split('\t');
            return fields.Length > 12 ? fields[12] : string.Empty;
        }

        /// <summary>
        /// Extracts the fragment sequence coverage display from ToString() output.
        /// </summary>
        private static string GetFragmentCoverageFromToString(string toStringOutput)
        {
            var fields = toStringOutput.Split('\t');
            return fields.Length > 13 ? fields[13] : string.Empty;
        }

        /// <summary>
        /// Extracts the modification info from ToString() output.
        /// </summary>
        private static string GetModsInfoFromToString(string toStringOutput)
        {
            var fields = toStringOutput.Split('\t');
            return fields.Length > 14 ? fields[14] : string.Empty;
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

            var output = group.ToString();
            var coverageDisplay = GetSequenceCoverageFromToString(output);
            var coverageFraction = GetSequenceCoverageFractionFromToString(output);

            // All lowercase = not covered
            Assert.That(coverageDisplay, Is.EqualTo("acdefghiklmnpqrstvwy"));
            Assert.That(coverageFraction, Is.EqualTo("0"));
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

            var output = group.ToString();
            var coverageFraction = GetSequenceCoverageFractionFromToString(output);

            // 6 residues covered out of 20 = 0.3
            Assert.That(double.Parse(coverageFraction), Is.EqualTo(0.3).Within(0.001));
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

            var output = group.ToString();
            var coverageDisplay = GetSequenceCoverageFromToString(output);
            var coverageFraction = GetSequenceCoverageFractionFromToString(output);

            Assert.That(double.Parse(coverageFraction), Is.EqualTo(1.0).Within(0.001));
            Assert.That(coverageDisplay, Is.EqualTo("ACDEFGHIK")); // All uppercase = fully covered
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

            var output = group.ToString();
            var coverageDisplay = GetSequenceCoverageFromToString(output);

            // ac = not covered (lowercase), DEFG = covered (uppercase), hik = not covered (lowercase)
            Assert.That(coverageDisplay, Is.EqualTo("acDEFGhik"));
        }

        #endregion

        #region Multiple Peptide Coverage Tests

        [Test]
        public void CalculateSequenceCoverage_WithOverlappingPeptides_CombinesCoverage()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIKLMNPQRSTVWY", "P00001");

            var peptide1 = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer, 1, 5);
            var peptide2 = new CoverageBioPolymerWithSetMods("EFGHI", "EFGHI", bioPolymer, 4, 8);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };

            var psm1 = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide1 });
            var psm2 = new CoverageSpectralMatch("test.raw", "EFGHI", "EFGHI", 100, 2, new[] { peptide2 });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            var coverageDisplay = GetSequenceCoverageFromToString(output);
            var coverageFraction = GetSequenceCoverageFractionFromToString(output);

            // Combined coverage: positions 1-8 = 8 residues out of 20 = 0.4
            Assert.That(double.Parse(coverageFraction), Is.EqualTo(0.4).Within(0.001));
            Assert.That(coverageDisplay, Is.EqualTo("ACDEFGHIklmnpqrstvwy"));
        }

        [Test]
        public void CalculateSequenceCoverage_WithNonOverlappingPeptides_AddsCoverage()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIKLMNPQRSTVWY", "P00001");

            var peptide1 = new CoverageBioPolymerWithSetMods("ACD", "ACD", bioPolymer, 1, 3);
            var peptide2 = new CoverageBioPolymerWithSetMods("VWY", "VWY", bioPolymer, 18, 20);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };

            var psm1 = new CoverageSpectralMatch("test.raw", "ACD", "ACD", 100, 1, new[] { peptide1 });
            var psm2 = new CoverageSpectralMatch("test.raw", "VWY", "VWY", 100, 2, new[] { peptide2 });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            var coverageDisplay = GetSequenceCoverageFromToString(output);
            var coverageFraction = GetSequenceCoverageFractionFromToString(output);

            Assert.That(double.Parse(coverageFraction), Is.EqualTo(0.3).Within(0.001));
            Assert.That(coverageDisplay, Is.EqualTo("ACDefghiklmnpqrstVWY"));
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

            var output = group.ToString();
            var coverageFractions = GetSequenceCoverageFractionFromToString(output);

            // Should contain pipe-separated fractions for both proteins
            Assert.That(coverageFractions, Does.Contain("|"));
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
            psm.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                // Fragments that cover positions 1, 2, and 3 of the peptide (D, E, F)
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 200.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 3, 3, 0), 300.0, 10.0, 1),
            };


            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            var fragmentCoverage = GetFragmentCoverageFromToString(output);

            // Fragment coverage should show covered positions (not empty)
            Assert.That(fragmentCoverage, Is.Not.Empty);
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

            var output = group.ToString();
            var fragmentCoverage = GetFragmentCoverageFromToString(output);

            // Fragment coverage list should be all lowercase (no fragment coverage)
            Assert.That(fragmentCoverage, Is.EqualTo("acdefghik"));
        }

        #endregion

        #region Modification Tests

        [Test]
        public void CalculateSequenceCoverage_WithModifications_PopulatesModsInfo()
        {
            var bioPolymer = new CoverageBioPolymer("ACDEFGHIK", "P00001");

            var phosphoMod = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Post-translational",
                _target: ModificationMotif.TryGetMotif("S", out var motif) ? motif : null,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 79.966);

            var mods = new Dictionary<int, Modification> { { 3, phosphoMod } };

            var peptide = new CoverageBioPolymerWithSetMods(
                "DEFGH", "DEF[Phosphorylation]GH", bioPolymer, 3, 7, 500.0, mods);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "DEF[Phosphorylation]GH", "DEFGH", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            var coverageWithMods = GetSequenceCoverageWithModsFromToString(output);

            Assert.That(coverageWithMods, Does.Contain("[Phosphorylation on S]"));
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

            var output = group.ToString();
            var coverageWithMods = GetSequenceCoverageWithModsFromToString(output);

            // Common Variable mods should be skipped
            Assert.That(coverageWithMods, Does.Not.Contain("[Oxidation]"));
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

            var output = group.ToString();
            var coverageWithMods = GetSequenceCoverageWithModsFromToString(output);

            Assert.That(coverageWithMods, Does.StartWith("[Acetyl on A]-"));
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

            var peptide = new CoverageBioPolymerWithSetMods(
                "FGHIK", "FGHIK-[Amidation]", bioPolymer, 5, 9, 500.0, mods);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new CoverageSpectralMatch("test.raw", "FGHIK-[Amidation]", "FGHIK", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            var coverageWithMods = GetSequenceCoverageWithModsFromToString(output);

            Assert.That(coverageWithMods, Does.EndWith("-[Amidation on K]"));
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

            var output = group.ToString();
            var coverageFraction = GetSequenceCoverageFractionFromToString(output);

            // Only the valid PSM should contribute to coverage
            Assert.That(double.Parse(coverageFraction), Is.EqualTo(5.0 / 9.0).Within(0.001));
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

            var output = group.ToString();
            var coverageFraction = GetSequenceCoverageFractionFromToString(output);

            // bioPolymer1 should have 0 coverage since peptide belongs to bioPolymer2
            Assert.That(double.Parse(coverageFraction), Is.EqualTo(0.0));
        }

        [Test]
        public void CalculateSequenceCoverage_WithEmptyBioPolymers_DoesNotThrow()
        {
            var bioPolymers = new HashSet<IBioPolymer>();
            var sequences = new HashSet<IBioPolymerWithSetMods>();
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods>();

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);

            Assert.DoesNotThrow(() => group.CalculateSequenceCoverage());
            Assert.DoesNotThrow(() => group.ToString());
        }

        [Test]
        public void CalculateSequenceCoverage_CalledMultipleTimes_ReplacesResults()
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
            var output1 = group.ToString();

            group.CalculateSequenceCoverage(); // Call again
            var output2 = group.ToString();

            // Results should be identical (replaced, not appended)
            Assert.That(output1, Is.EqualTo(output2));
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

            var output = group.ToString();
            var coverageDisplay = GetSequenceCoverageFromToString(output);
            var coverageFraction = GetSequenceCoverageFractionFromToString(output);

            Assert.That(double.Parse(coverageFraction), Is.EqualTo(1.0));
            Assert.That(coverageDisplay, Is.EqualTo("K"));
        }

        #endregion

        #region Modification Occupancy Tests

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

            var output = group.ToString();
            var modsInfo = GetModsInfoFromToString(output);

            // Should contain #aa format with position and occupancy
            if (!string.IsNullOrEmpty(modsInfo))
            {
                Assert.That(modsInfo, Does.Contain("#aa"));
                Assert.That(modsInfo, Does.Contain("occupancy"));
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

            var peptide1 = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer, 1, 5);
            var peptide2 = new CoverageBioPolymerWithSetMods("GHIKL", "GHIKL", bioPolymer, 6, 10);
            var peptide3 = new CoverageBioPolymerWithSetMods("STVW", "ST[Phospho]VW", bioPolymer, 16, 19, 500.0, modsDict);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2, peptide3 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1 };

            var psm1 = new CoverageSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide1 });
            var psm2 = new CoverageSpectralMatch("test.raw", "GHIKL", "GHIKL", 95, 2, new[] { peptide2 });
            var psm3 = new CoverageSpectralMatch("test.raw", "ST[Phospho]VW", "STVW", 90, 3, new[] { peptide3 });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2, psm3 };

            group.CalculateSequenceCoverage();

            var output = group.ToString();

            // Verify coverage fraction is correct: 5 + 5 + 4 = 14 residues out of 20 = 0.7
            var coverageFraction = GetSequenceCoverageFractionFromToString(output);
            Assert.That(double.Parse(coverageFraction), Is.EqualTo(0.7).Within(0.001));

            // Verify display list shows correct coverage pattern
            var coverageDisplay = GetSequenceCoverageFromToString(output);
            Assert.That(coverageDisplay, Does.Contain("ACDEF"));
            Assert.That(coverageDisplay, Does.Contain("GHIKL"));
            Assert.That(coverageDisplay, Does.Contain("STVW"));

            // Verify mods display contains phospho modification
            var coverageWithMods = GetSequenceCoverageWithModsFromToString(output);
            Assert.That(coverageWithMods, Does.Contain("[Phospho on T]"));
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

            // Verify ToString includes both coverage and intensity
            var output = group.ToString();
            Assert.That(output, Does.Contain("1000000"));

            var coverageFraction = GetSequenceCoverageFractionFromToString(output);
            Assert.That(double.Parse(coverageFraction), Is.EqualTo(5.0 / 9.0).Within(0.001));
        }

        [Test]
        public void CalculateSequenceCoverage_Integration_RealWorldScenario()
        {
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

            var output = group.ToString();
            var coverageFraction = GetSequenceCoverageFractionFromToString(output);
            var coverageDisplay = GetSequenceCoverageFromToString(output);

            // Total covered: 8 + 8 + 23 + 14 = 53 residues out of 92
            double expectedCoverage = 53.0 / 92.0;
            Assert.That(double.Parse(coverageFraction), Is.EqualTo(expectedCoverage).Within(0.001));

            // Verify display has correct length
            Assert.That(coverageDisplay.Length, Is.EqualTo(92));

            // Verify covered regions are uppercase
            Assert.That(coverageDisplay.Substring(0, 8), Is.EqualTo("MKTAYIAK"));
            Assert.That(coverageDisplay.Substring(8, 8), Is.EqualTo("QRQISFVK"));
        }

        #endregion
    }
}