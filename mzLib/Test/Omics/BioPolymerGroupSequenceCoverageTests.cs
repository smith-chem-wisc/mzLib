using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectralMatch;
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
        #region Helper Methods

        private static string GetCoverageDisplay(string output) => output.Split('\t').ElementAtOrDefault(11) ?? "";
        private static string GetCoverageFraction(string output) => output.Split('\t').ElementAtOrDefault(10) ?? "";
        private static string GetFragmentCoverage(string output) => output.Split('\t').ElementAtOrDefault(13) ?? "";

        private static BioPolymerGroup CreateGroupWithPsm(MockBioPolymer protein, MockBioPolymerWithSetMods peptide)
        {
            var psm = new MockSpectralMatch("test.raw", peptide.BaseSequence, peptide.FullSequence, 100, 1, new[] { peptide });
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
            var protein = new MockBioPolymer("ACDEFGHIKLMNPQRSTVWY", "P00001", false);
            var peptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);
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
            var protein = new MockBioPolymer("ACDEFGHIK", "P00001", false);
            var peptide = new MockBioPolymerWithSetMods("DEF", "DEF", protein, 3, 5);
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
            var protein = new MockBioPolymer("ACDEFGHIK", "P00001"); // 9 residues
            var peptide1 = new MockBioPolymerWithSetMods("ACDE", "ACDE", protein, 1, 4);
            var peptide2 = new MockBioPolymerWithSetMods("DEFG", "DEFG", protein, 3, 6); // Overlaps at DE

            var psm1 = new MockSpectralMatch("test.raw", "ACDE","ACDE", 100, 1, new[] { peptide1 });
            var psm2 = new MockSpectralMatch("test.raw", "DEFG","DEFG", 100, 2, new[] { peptide2 });

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

        #endregion

        #region Robustness Tests

        /// <summary>
        /// Ensures the method handles empty PSM collections without throwing.
        /// Critical: Prevents crashes during edge-case processing.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_NoPsms_ReturnsZeroCoverage()
        {
            var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
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
            var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

            var validPsm = new MockSpectralMatch("test.raw", "ACDEF","ACDEF", 100, 1, new[] { peptide });
            var nullPsm = new MockSpectralMatch("test.raw", null!, null!, 100, 2, new[] { peptide });

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
            var protein1 = new MockBioPolymer("ACDEFGHIK", "P00001");
            var protein2 = new MockBioPolymer("MNPQRSTVWY", "P00002");
            var peptideFromProtein2 = new MockBioPolymerWithSetMods("MNPQR", "MNPQR", protein2, 1, 5);

            var psm = new MockSpectralMatch("test.raw", "MNPQR", "MNPQR", 100, 1, new[] { peptideFromProtein2 });

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


        [Test]
        public void CalculateSequenceCoverage_WithNonOverlappingPeptides_AddsCoverage()
        {
            var bioPolymer = new MockBioPolymer("ACDEFGHIKLMNPQRSTVWY", "P00001");

            var peptide1 = new MockBioPolymerWithSetMods("ACD", "ACD", bioPolymer, 1, 3);
            var peptide2 = new MockBioPolymerWithSetMods("VWY", "VWY", bioPolymer, 18, 20);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };

            var psm1 = new MockSpectralMatch("test.raw", "ACD", "ACD", 100, 1, new[] { peptide1 });
            psm1.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                // Fragments that cover positions 1, 2, and 3 of the peptide (A,B,C)
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 1, 1, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 2, 2, 0), 200.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 10, 3, 3, 0), 300.0, 10.0, 1),
            };
            var psm2 = new MockSpectralMatch("test.raw", "VWY", "VWY", 100, 2, new[] { peptide2 });
            psm2.MatchedFragmentIons = new List<MatchedFragmentIon>
            {
                // Fragments that cover positions 1, 2, and 3 of the peptide (V,W,Y)
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 1, 3, 0), 100.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 2, 2, 0), 200.0, 10.0, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 10, 3, 1, 0), 300.0, 10.0, 1),
            };

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            var coverageDisplay = group.CoverageResult.FragmentSequenceCoverageDisplayList.FirstOrDefault();
            var coverageFraction = group.CoverageResult.SequenceCoverageFraction.FirstOrDefault();
            Assert.That(coverageFraction, Is.EqualTo(0.3).Within(0.001));
            Assert.That(coverageDisplay, Is.EqualTo("ACDefghiklmnpqrstVWY"));
        }

        #endregion

        #region Multiple BioPolymer Tests

        [Test]
        public void CalculateSequenceCoverage_WithMultipleBioPolymers_CalculatesSeparately()
        {
            var bioPolymer1 = new MockBioPolymer("ACDEFGHIK", "P00001");
            var bioPolymer2 = new MockBioPolymer("MNPQRSTVWY", "P00002");

            var peptide1 = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", bioPolymer1, 1, 5);
            var peptide2 = new MockBioPolymerWithSetMods("MNPQR", "MNPQR", bioPolymer2, 1, 5);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer1, bioPolymer2 };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 };

            var psm1 = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide1 });
            var psm2 = new MockSpectralMatch("test.raw", "MNPQR", "MNPQR", 100, 2, new[] { peptide2 });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            // Should contain coverage fractions for both proteins
            Assert.That(group.CoverageResult.SequenceCoverageFraction.Count, Is.GreaterThan(1));
            Assert.That(group.CoverageResult.SequenceCoverageDisplayList.Count, Is.GreaterThan(1));
            Assert.That(group.CoverageResult.SequenceCoverageDisplayListWithMods.Count, Is.GreaterThan(1));
            Assert.That(group.CoverageResult.FragmentSequenceCoverageDisplayList.Count, Is.GreaterThan(1));
            Assert.That(group.CoverageResult.ModsInfo.Count, Is.EqualTo(0));
        }

        #endregion

        #region Fragment Coverage Tests

        /// <summary>
        /// Verifies that Common Variable modifications (like oxidation) are excluded from output.
        /// Critical: These are technical artifacts, not biological modifications.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_CommonVariableMods_ExcludedFromOutput()
        {
            var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
            ModificationMotif.TryGetMotif("D", out var motif);
            var commonMod = new Modification("Oxidation", null, "Common Variable", null, motif, "Anywhere.", null, 15.995);

            var mods = new Dictionary<int, Modification> { { 3, commonMod } };
            var peptide = new MockBioPolymerWithSetMods("ACDEF", "ACD[Oxidation]EF", protein, 1, 5, mods);

            var psm = new MockSpectralMatch("test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide });
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