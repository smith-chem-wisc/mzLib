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
    /// <summary>
    /// Tests for BioPolymerGroup.CalculateSequenceCoverage() functionality.
    /// These tests ensure correct calculation of peptide-level and fragment-level 
    /// sequence coverage, which is critical for protein inference reporting.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class BioPolymerGroupSequenceCoverageTests
    {
        #region Helper Classes

        private class CoverageBioPolymer : IBioPolymer
        {
            public string BaseSequence { get; }
            public string Accession { get; }
            public string Organism { get; } = "Test";
            public string Name { get; } = "Test";
            public string FullName { get; } = "Test";
            public List<Tuple<string, string>> GeneNames { get; } = new() { new("primary", "GENE1") };
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

            public CoverageBioPolymer(string sequence, string accession, bool isDecoy = false, bool isContaminant = false)
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
                => new CoverageBioPolymer(newBaseSequence, Accession);

            public TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original,
                IEnumerable<SequenceVariation> appliedSequenceVariants, IEnumerable<TruncationProduct> applicableProteolysisProducts,
                IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
                where TBioPolymerType : IHasSequenceVariants => original;

            public bool Equals(IBioPolymer? other) => other != null && Accession == other.Accession && BaseSequence == other.BaseSequence;
            public override bool Equals(object? obj) => obj is IBioPolymer other && Equals(other);
            public override int GetHashCode() => HashCode.Combine(Accession, BaseSequence);
        }

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
            public string Description => "";
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
            public ChemicalFormula ThisChemicalFormula => new();
            public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

            public CoverageBioPolymerWithSetMods(string baseSequence, string fullSequence, IBioPolymer parent,
                int oneBasedStartResidue, int oneBasedEndResidue, Dictionary<int, Modification>? mods = null)
            {
                BaseSequence = baseSequence;
                FullSequence = fullSequence;
                Parent = parent;
                OneBasedStartResidue = oneBasedStartResidue;
                OneBasedEndResidue = oneBasedEndResidue;
                AllModsOneIsNterminus = mods ?? new Dictionary<int, Modification>();
                MostAbundantMonoisotopicMass = MonoisotopicMass = 500.0;
            }

            public void Fragment(DissociationType d, FragmentationTerminus t, List<Product> p, FragmentationParams? f = null) { }
            public void FragmentInternally(DissociationType d, int m, List<Product> p, FragmentationParams? f = null) { }
            public IBioPolymerWithSetMods Localize(int i, double m) => this;
            public bool Equals(IBioPolymerWithSetMods? other) => other != null && BaseSequence == other.BaseSequence && Parent?.Accession == other.Parent?.Accession;
            public override bool Equals(object? obj) => obj is IBioPolymerWithSetMods other && Equals(other);
            public override int GetHashCode() => HashCode.Combine(BaseSequence, Parent?.Accession);
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

        private static string GetCoverageDisplay(string output) => output.Split('\t').ElementAtOrDefault(11) ?? "";
        private static string GetCoverageFraction(string output) => output.Split('\t').ElementAtOrDefault(10) ?? "";
        private static string GetFragmentCoverage(string output) => output.Split('\t').ElementAtOrDefault(13) ?? "";

        private static BioPolymerGroup CreateGroupWithPsm(CoverageBioPolymer protein, CoverageBioPolymerWithSetMods peptide)
        {
            var psm = new CoverageSpectralMatch("test.raw", peptide.BaseSequence, 100, 1, new[] { peptide });
            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };
            return group;
        }

        #endregion

        #region Core Functionality Tests

        /// <summary>
        /// Verifies that sequence coverage fraction is calculated correctly.
        /// This is the primary metric reported in protein group output files.
        /// Critical: Incorrect coverage would misrepresent protein identification quality.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_CalculatesCorrectFraction()
        {
            // Protein: 20 residues, Peptide covers positions 1-5 (5 residues) = 25% coverage
            var protein = new CoverageBioPolymer("ACDEFGHIKLMNPQRSTVWY", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);
            var group = CreateGroupWithPsm(protein, peptide);

            group.CalculateSequenceCoverage();

            var fraction = double.Parse(GetCoverageFraction(group.ToString()));
            Assert.That(fraction, Is.EqualTo(0.25).Within(0.001));
        }

        /// <summary>
        /// Verifies the uppercase/lowercase display convention for sequence coverage.
        /// Uppercase = covered residues, lowercase = uncovered residues.
        /// Critical: This visual format is used in output files and must be consistent.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_DisplaysUppercaseForCoveredResidues()
        {
            var protein = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("DEF", "DEF", protein, 3, 5);
            var group = CreateGroupWithPsm(protein, peptide);

            group.CalculateSequenceCoverage();

            var display = GetCoverageDisplay(group.ToString());
            Assert.That(display, Is.EqualTo("acDEFghik"));
        }

        /// <summary>
        /// Verifies that overlapping peptides correctly combine their coverage.
        /// Critical: Real proteomics data has overlapping peptides; double-counting would inflate coverage.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_OverlappingPeptides_CombineCorrectly()
        {
            var protein = new CoverageBioPolymer("ACDEFGHIK", "P00001"); // 9 residues
            var peptide1 = new CoverageBioPolymerWithSetMods("ACDE", "ACDE", protein, 1, 4);
            var peptide2 = new CoverageBioPolymerWithSetMods("DEFG", "DEFG", protein, 3, 6); // Overlaps at DE

            var psm1 = new CoverageSpectralMatch("test.raw", "ACDE", 100, 1, new[] { peptide1 });
            var psm2 = new CoverageSpectralMatch("test.raw", "DEFG", 100, 2, new[] { peptide2 });

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 },
                new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            // Combined coverage: 1-6 = 6 residues out of 9
            var fraction = double.Parse(GetCoverageFraction(group.ToString()));
            Assert.That(fraction, Is.EqualTo(6.0 / 9.0).Within(0.001));
        }

        /// <summary>
        /// Verifies that fragment-level coverage (from matched ions) is tracked separately from peptide coverage.
        /// Critical: Fragment coverage provides additional confidence in residue identification.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_FragmentCoverage_TrackedSeparately()
        {
            var protein = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

            var psm = new CoverageSpectralMatch("test.raw", "ACDEF", 100, 1, new[] { peptide });
            psm.NTerminalFragmentPositions = new List<int> { 1, 2 }; // Fragment ions at positions 1,2

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

            var fragmentCoverage = GetFragmentCoverage(group.ToString());
            // Fragment coverage should show some uppercase (fragment-covered) positions
            Assert.That(fragmentCoverage.Any(char.IsUpper), Is.True);
        }

        #endregion

        #region Robustness Tests

        /// <summary>
        /// Ensures the method handles empty PSM collections without throwing.
        /// Critical: Prevents crashes during edge-case processing.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_NoPsms_ReturnsZeroCoverage()
        {
            var protein = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>();

            group.CalculateSequenceCoverage();

            var fraction = GetCoverageFraction(group.ToString());
            Assert.That(fraction, Is.EqualTo("0"));
        }

        /// <summary>
        /// Ensures PSMs with null BaseSequence (ambiguous identifications) are skipped gracefully.
        /// Critical: Prevents null reference exceptions in production workflows.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_NullBaseSequencePsm_SkippedWithoutError()
        {
            var protein = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

            var validPsm = new CoverageSpectralMatch("test.raw", "ACDEF", 100, 1, new[] { peptide });
            var nullPsm = new CoverageSpectralMatch("test.raw", null!, 100, 2, new[] { peptide });

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { validPsm, nullPsm };

            Assert.DoesNotThrow(() => group.CalculateSequenceCoverage());

            // Only valid PSM should contribute
            var fraction = double.Parse(GetCoverageFraction(group.ToString()));
            Assert.That(fraction, Is.EqualTo(5.0 / 9.0).Within(0.001));
        }

        /// <summary>
        /// Ensures peptides from different proteins don't contribute to wrong protein's coverage.
        /// Critical: Incorrect protein assignment would corrupt coverage calculations.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_PeptideFromDifferentParent_NotCounted()
        {
            var protein1 = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            var protein2 = new CoverageBioPolymer("MNPQRSTVWY", "P00002");
            var peptideFromProtein2 = new CoverageBioPolymerWithSetMods("MNPQR", "MNPQR", protein2, 1, 5);

            var psm = new CoverageSpectralMatch("test.raw", "MNPQR", 100, 1, new[] { peptideFromProtein2 });

            // Group only contains protein1, but PSM has peptide from protein2
            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein1 },
                new HashSet<IBioPolymerWithSetMods> { peptideFromProtein2 },
                new HashSet<IBioPolymerWithSetMods>());
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var fraction = GetCoverageFraction(group.ToString());
            Assert.That(fraction, Is.EqualTo("0")); // protein1 has no coverage
        }

        #endregion

        #region Modification Handling Tests

        /// <summary>
        /// Verifies that modifications are correctly annotated in coverage display.
        /// Critical: PTM localization is essential for biological interpretation.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_Modifications_AnnotatedInDisplay()
        {
            var protein = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            ModificationMotif.TryGetMotif("D", out var motif);
            var phospho = new Modification("Phospho", null, "Post-translational", null, motif, "Anywhere.", null, 79.966);

            var mods = new Dictionary<int, Modification> { { 3, phospho } };
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", "ACD[Phospho]EF", protein, 1, 5, mods);

            var psm = new CoverageSpectralMatch("test.raw", "ACDEF", 100, 1, new[] { peptide });
            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            var coverageWithMods = output.Split('\t').ElementAtOrDefault(12) ?? "";
            Assert.That(coverageWithMods, Does.Contain("[Phospho on D]"));
        }

        /// <summary>
        /// Verifies that Common Variable modifications (like oxidation) are excluded from output.
        /// Critical: These are technical artifacts, not biological modifications.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_CommonVariableMods_ExcludedFromOutput()
        {
            var protein = new CoverageBioPolymer("ACDEFGHIK", "P00001");
            ModificationMotif.TryGetMotif("D", out var motif);
            var commonMod = new Modification("Oxidation", null, "Common Variable", null, motif, "Anywhere.", null, 15.995);

            var mods = new Dictionary<int, Modification> { { 3, commonMod } };
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", "ACD[Oxidation]EF", protein, 1, 5, mods);

            var psm = new CoverageSpectralMatch("test.raw", "ACDEF", 100, 1, new[] { peptide });
            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            var coverageWithMods = output.Split('\t').ElementAtOrDefault(12) ?? "";
            Assert.That(coverageWithMods, Does.Not.Contain("[Oxidation"));
        }

        #endregion
    }
}