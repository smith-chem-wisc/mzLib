using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class PeptideWithSetModsDecoyTests
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        /// <summary>
        /// CRITICAL: Tests reverse decoy peptide generation from target peptides.
        /// Validates correct sequence reversal, modification position remapping,
        /// cleavage site preservation, and handling of palindromic sequences.
        /// Essential for target-decoy FDR estimation in all database searches.
        /// </summary>
        [Test]
        public static void TestReverseDecoyFromTarget()
        {
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif_p);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif_p, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

            ModificationMotif.TryGetMotif("K", out ModificationMotif motif_k);
            Modification acetylation = new Modification(_originalId: "acetyl", _modificationType: "CommonBiological", _target: motif_k, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

            Modification nTermAcet = new Modification(_originalId: "n-acetyl", _modificationType: "CommonBiological", _target: null, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"));

            Dictionary<int, Modification> allmodsoneisnterminus = new Dictionary<int, Modification> { { 1, nTermAcet }, { 6, phosphorylation }, { 9, acetylation } };

            PeptideWithSetModifications p = new PeptideWithSetModifications(new Protein("PEPTIDEK", "ACCESSIION"), new DigestionParams(), 1, 8, CleavageSpecificity.Full, null, 0, allmodsoneisnterminus, 0, null);

            int[] newAminoAcidPositions = new int["PEPTIDEK".Length];
            PeptideWithSetModifications reverse = p.GetReverseDecoyFromTarget(newAminoAcidPositions);

            string targetSequence = p.FullSequence;
            string decoySequence = reverse.FullSequence; 
            Assert.AreEqual(reverse.PairedTargetDecoySequence, targetSequence);
            Assert.AreEqual(p.PairedTargetDecoySequence, decoySequence);
            Assert.AreEqual("EDITPEPK", reverse.BaseSequence);
            Assert.AreEqual(new int[] { 6, 5, 4, 3, 2, 1, 0, 7 }, newAminoAcidPositions);
            Assert.IsTrue(reverse.Protein.IsDecoy);
            Assert.AreEqual(p.Protein.BaseSequence.Length, reverse.Protein.BaseSequence.Length);//we replaced the original with the new so the protein should have the same length
            Assert.AreEqual(reverse.PeptideDescription, p.FullSequence);

            List<Product> decoyProducts = new List<Product>();
            reverse.Fragment(MassSpectrometry.DissociationType.HCD, FragmentationTerminus.Both, decoyProducts);

            Assert.AreEqual(14, decoyProducts.Count);

            //  Arg-C -- Cleave after R
            newAminoAcidPositions = new int["RPEPTIREK".Length];
            PeptideWithSetModifications p_argC = new PeptideWithSetModifications(new Protein("RPEPTIREK", "DECOY_ARGC"), new DigestionParams(protease: "Arg-C"), 1, 9, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications p_argC_reverse = p_argC.GetReverseDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual("RKEITPREP", p_argC_reverse.BaseSequence);
            Assert.AreEqual(new int[] { 0, 8, 7, 5, 4, 3, 6, 2, 1 }, newAminoAcidPositions);
            Assert.AreEqual(p_argC_reverse.PeptideDescription, p_argC.FullSequence);

            //  Asp-N -- Cleave before D
            newAminoAcidPositions = new int["DPEPTIDEK".Length];
            PeptideWithSetModifications p_aspN = new PeptideWithSetModifications(new Protein("DPEPTIDEK", "DECOY_ASPN"), new DigestionParams(protease: "Asp-N"), 1, 9, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications p_aspN_reverse = p_aspN.GetReverseDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual("DKEITPDEP", p_aspN_reverse.BaseSequence);
            Assert.AreEqual(p_aspN.FullSequence, p_aspN_reverse.PeptideDescription);

            //  chymotrypsin|P
            newAminoAcidPositions = new int["FKFPRWAWPSYGYPG".Length];
            PeptideWithSetModifications p_chymoP = new PeptideWithSetModifications(new Protein("FKFPRWAWPSYGYPG", "DECOY_CHYMOP"), new DigestionParams(protease: "chymotrypsin|P", maxMissedCleavages: 10), 1, 15, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications p_chymoP_reverse = p_chymoP.GetReverseDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual("FGPYGWSPWAYRPFK", p_chymoP_reverse.BaseSequence);
            Assert.AreEqual(p_chymoP.FullSequence, p_chymoP_reverse.PeptideDescription);

            //  CNBr cleave after M
            newAminoAcidPositions = new int["MPEPTIMEK".Length];
            PeptideWithSetModifications p_cnbr = new PeptideWithSetModifications(new Protein("MPEPTIMEK", "DECOY_CNBR"), new DigestionParams(protease: "CNBr"), 1, 9, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications p_cnbr_reverse = p_cnbr.GetReverseDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual("MKEITPMEP", p_cnbr_reverse.BaseSequence);
            Assert.AreEqual(p_cnbr.FullSequence, p_cnbr_reverse.PeptideDescription);

            //  elastase cleave after A, V, S, G, L, I,
            newAminoAcidPositions = new int["KAYVPSRGHLDIN".Length];
            PeptideWithSetModifications p_elastase = new PeptideWithSetModifications(new Protein("KAYVPSRGHLDIN", "DECOY_ELASTASE"), new DigestionParams(protease: "elastase|P"), 1, 13, CleavageSpecificity.Semi, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications p_elastase_reverse = p_elastase.GetReverseDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual("NADHRSPGVLYIK", p_elastase_reverse.BaseSequence);

            //  top-down
            newAminoAcidPositions = new int["RPEPTIREK".Length];
            PeptideWithSetModifications p_topdown = new PeptideWithSetModifications(new Protein("RPEPTIREK", "DECOY_TD"), new DigestionParams(protease: "top-down"), 1, 9, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications p_topdown_reverse = p_topdown.GetReverseDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual("KERITPEPR", p_topdown_reverse.BaseSequence);

            //  Arg-C -- Cleave after R
            newAminoAcidPositions = new int["PEPTGPYGPYIDE".Length];
            PeptideWithSetModifications p_coll = new PeptideWithSetModifications(new Protein("PEPTGPYGPYIDE", "DECOY_COL"), new DigestionParams(protease: "collagenase"), 1, 13, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications p_coll_reverse = p_coll.GetReverseDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual("EDITGPYGPYPEP", p_coll_reverse.BaseSequence);

            // Tyrpsin -- Reverse Decoy is identical so get mirror
            Dictionary<int, Modification> VTIRTVR_modsDictionary = new Dictionary<int, Modification> { { 1, nTermAcet }, { 3, phosphorylation } };
            newAminoAcidPositions = new int["VTIRTVR".Length];
            PeptideWithSetModifications p_tryp = new PeptideWithSetModifications(new Protein("VTIRTVR", "DECOY_TRYP"), new DigestionParams(protease: "trypsin"), 1, 7, CleavageSpecificity.Full, null, 0, VTIRTVR_modsDictionary, 0, null);
            PeptideWithSetModifications p_tryp_reverse = p_tryp.GetReverseDecoyFromTarget(newAminoAcidPositions);
            // Hash code corresponding to the target sequence, should be PairedTargetDecoyHash for reverse
            string mirrorTarget = p_tryp.FullSequence;
            // Hash code corresponding to the decoy sequence, should be PairedTargetDecoyHash for target
            string mirrorDecoy = p_tryp_reverse.FullSequence;
            Assert.AreEqual(mirrorTarget, p_tryp_reverse.PairedTargetDecoySequence);
            Assert.AreEqual(mirrorDecoy, p_tryp.PairedTargetDecoySequence);
            Assert.AreEqual("RVTRITV", p_tryp_reverse.BaseSequence);
            Assert.AreEqual(new int[] { 6, 5, 4, 3, 2, 1, 0 }, newAminoAcidPositions);
            Assert.IsTrue(p_tryp_reverse.AllModsOneIsNterminus.ContainsKey(1));//n-term acetyl
            Assert.IsTrue(p_tryp_reverse.AllModsOneIsNterminus.ContainsKey(7));//moved T-phospho from 3 to 7
            Assert.IsTrue(p_tryp_reverse.Protein.IsDecoy);
            Assert.AreEqual(p_tryp.FullSequence, p_tryp_reverse.PeptideDescription);

        }
        /// <summary>
        /// CRITICAL: Tests scrambled decoy peptide generation from target peptides.
        /// Scrambled decoys are an alternative to reversed decoys that may provide
        /// better FDR estimation for certain peptide populations. Tests modification
        /// position mapping and cleavage site preservation during scrambling.
        /// </summary>
        [Test]
        public static void TestScrambledDecoyFromTarget()
        {
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif_p);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif_p, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

            ModificationMotif.TryGetMotif("K", out ModificationMotif motif_k);
            Modification acetylation = new Modification(_originalId: "acetyl", _modificationType: "CommonBiological", _target: motif_k, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

            Modification nTermAcet = new Modification(_originalId: "n-acetyl", _modificationType: "CommonBiological", _target: null, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"));

            Dictionary<int, Modification> allmodsoneisnterminus = new Dictionary<int, Modification> { { 1, nTermAcet }, { 5, phosphorylation }, { 9, acetylation } };

            PeptideWithSetModifications p = new PeptideWithSetModifications(new Protein("PEPTIDEK", "ACCESSIION"), new DigestionParams(), 1, 8, CleavageSpecificity.Full, null, 0, allmodsoneisnterminus, 0, null);
            int[] newAminoAcidPositions = new int["PEPTIDEK".Length];
            PeptideWithSetModifications testScrambled = p.GetScrambledDecoyFromTarget(newAminoAcidPositions);

            string targetSequence = p.FullSequence;
            string decoySequence = testScrambled.FullSequence;
            Assert.AreEqual(testScrambled.PairedTargetDecoySequence, targetSequence);
            Assert.AreEqual(p.PairedTargetDecoySequence, decoySequence);
            Assert.AreEqual("IDEETPPK", testScrambled.BaseSequence);
            Assert.AreEqual(new int[] { 4, 5, 6, 1, 3, 0, 2, 7 }, newAminoAcidPositions);
            // Check n-term acetyl
            Assert.AreEqual(p.AllModsOneIsNterminus[1], testScrambled.AllModsOneIsNterminus[1]); 
            // Check phosphorylated Thr originally at position 5
            Assert.AreEqual(p.AllModsOneIsNterminus[5], testScrambled.AllModsOneIsNterminus[6]);
            Assert.IsTrue(testScrambled.Protein.IsDecoy);
            Assert.AreEqual(p.Protein.BaseSequence.Length, testScrambled.Protein.BaseSequence.Length);
            // Check that the scrambled PeptideDescription is equivalent to the original peptide's full sequence
            Assert.AreEqual(testScrambled.PeptideDescription, p.FullSequence);

            // Test for complex cleavage motif patterns
            Dictionary<int, Modification> complexTrypticMods = new Dictionary<int, Modification> { { 1, nTermAcet }, { 89, acetylation } };
            string insulinFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.fasta");
            string insulinSeq = ProteinDbLoader.LoadProteinFasta(insulinFile, true, DecoyType.None, false, out var dbError)[0].BaseSequence;
            PeptideWithSetModifications insulinTryptic = new PeptideWithSetModifications(new Protein(insulinSeq, "COMPLEX_TRYP"), new DigestionParams(), 1, 110, CleavageSpecificity.Full, null, 0, complexTrypticMods, 0, null);
            int[] insulinScrambledOrder = new int[insulinSeq.Length];
            PeptideWithSetModifications insulinScrambled = insulinTryptic.GetScrambledDecoyFromTarget(insulinScrambledOrder);
            // Check that cleavage motif positions are preserved in the more complex case
            Assert.AreEqual(insulinTryptic.BaseSequence[5], insulinScrambled[5]);
            Assert.AreEqual(insulinScrambled.BaseSequence[55], insulinTryptic.BaseSequence[55]);          
            // Check that modified cleavage motif residues are preserved
            Assert.AreEqual(insulinScrambled.AllModsOneIsNterminus[89], insulinTryptic.AllModsOneIsNterminus[89]);

            // Test with Arg-C as the protease
            newAminoAcidPositions = new int["RPEPTIREAVLKK".Length];
            PeptideWithSetModifications testArgC = new PeptideWithSetModifications(new Protein("RPEPTIREAVLKK", "DECOY_ARGC"), new DigestionParams(protease: "Arg-C"), 1, 13, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications scrambledArgC = testArgC.GetScrambledDecoyFromTarget(newAminoAcidPositions);
            // Check for preserved cleavage motif positions
            Assert.AreEqual(testArgC.BaseSequence[6], scrambledArgC.BaseSequence[6]);

            // Tests that peptide is mirrored when the maximum number of scramble attempts is reached
            newAminoAcidPositions = new int["AVLRRRKKRDEF".Length];
            PeptideWithSetModifications forceMirror = new PeptideWithSetModifications(new Protein("AVLRRRKKRDEF", "ACCESSIION"), new DigestionParams(), 1, 12, CleavageSpecificity.Full, null, 0, allmodsoneisnterminus, 0, null);
            PeptideWithSetModifications mirroredTarget = forceMirror.GetScrambledDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual(new int[] { 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 }, newAminoAcidPositions);
        }
        /// <summary>
        /// CRITICAL: Tests decoy generation from peptides digested from real XML databases.
        /// Ensures decoy generation works correctly with database-derived modifications,
        /// sequence variants, and other annotations. No unchanged peptides should remain
        /// after decoy generation (except true palindromes).
        /// </summary>
        [Test]
        public static void TestReverseDecoyFromPeptideFromProteinXML()
        {
            //Just making sure there are no snafus when creating decoy peptides from an xml,which will have mods in various places, etc.
            //sequence variants, modifications
            Dictionary<string, Modification> un = new Dictionary<string, Modification>();
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "cRAP_databaseGPTMD.xml"), true, DecoyType.None, UniProtPtms, false, new string[] { "exclude_me" }, out un);

            List<Modification> fixedMods = new List<Modification>();
            List<Modification> variableMods = new List<Modification>();
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif_C);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif_M);

            fixedMods.Add(new Modification(_originalId: "resMod_C", _target: motif_C, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: PeriodicTable.GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "resMod_M", _target: motif_C, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("O"), _monoisotopicMass: PeriodicTable.GetElement(8).PrincipalIsotope.AtomicMass));

            int unchangedPeptides = 0;
            int totalPeptides = 0;

            foreach (Protein p in proteins)
            {
                List<PeptideWithSetModifications> targetPeptides = p.Digest(new DigestionParams(), fixedMods, variableMods, null, null).ToList();
                foreach (PeptideWithSetModifications targetPeptide in targetPeptides)
                {
                    totalPeptides++;
                    int[] newAminoAcidPositions = new int[targetPeptide.BaseSequence.Length];
                    PeptideWithSetModifications decoyPeptide = targetPeptide.GetReverseDecoyFromTarget(newAminoAcidPositions);

                    if (decoyPeptide.BaseSequence == targetPeptide.BaseSequence)
                    {
                        unchangedPeptides++;
                    }
                }
            }

            Assert.AreEqual(0, unchangedPeptides);
        }

        /// <summary>
        /// INFORMATIONAL: Counts how many target peptides have matching decoy sequences.
        /// This metric helps assess decoy database quality - too many matches indicate
        /// potential FDR estimation issues. Used for database quality validation.
        /// </summary>
        [Test]
        public static void CountTargetsWithMatchingDecoys()
        {
            Dictionary<string, Modification> un = new Dictionary<string, Modification>();
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "cRAP_databaseGPTMD.xml"), true, DecoyType.None, UniProtPtms, false, new string[] { "exclude_me" }, out un);

            List<Modification> fixedMods = new List<Modification>();
            List<Modification> variableMods = new List<Modification>();
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif_C);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif_M);

            fixedMods.Add(new Modification(_originalId: "resMod_C", _target: motif_C, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: PeriodicTable.GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "resMod_M", _target: motif_C, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("O"), _monoisotopicMass: PeriodicTable.GetElement(8).PrincipalIsotope.AtomicMass));

            Dictionary<string, int> targets = new Dictionary<string, int>();

            foreach (Protein p in proteins)
            {
                List<PeptideWithSetModifications> targetPeptides = p.Digest(new DigestionParams(), fixedMods, variableMods, null, null).ToList();

                foreach (PeptideWithSetModifications targetPeptide in targetPeptides)
                {
                    if (targets.ContainsKey(targetPeptide.BaseSequence))
                    {
                        targets[targetPeptide.BaseSequence]++;
                    }
                    else
                    {
                        targets.Add(targetPeptide.BaseSequence, 1);
                    }
                }
            }

            int matchingDecoys = 0;
            foreach (Protein p in proteins)
            {
                List<PeptideWithSetModifications> targetPeptides = p.Digest(new DigestionParams(), fixedMods, variableMods, null, null).ToList();

                foreach (PeptideWithSetModifications target in targetPeptides)
                {
                    int[] newAminoAcidPositions = new int[target.BaseSequence.Length];
                    string decoySequence = target.GetReverseDecoyFromTarget(newAminoAcidPositions).BaseSequence;

                    if (targets.ContainsKey(decoySequence))
                    {
                        matchingDecoys++;
                    }
                }
            }
        }
    }
}
