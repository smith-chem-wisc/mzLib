using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class PeptideWithSetModsDigestionTests
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
        /// CRITICAL: Tests semi-specific digestion doesn't crash with various protein sequences.
        /// Semi-specific search is essential for identifying endogenous peptides and degradation
        /// products. This test guards against edge cases that could crash the digestion algorithm.
        /// </summary>
        [Test]
        public static void TestSemiFewCleavages()
        {
            Protein protein = new("MQLLRCFSIFSVIASVLAQELTTICEQIPSPTLESTPYSLSTTTILANGKAMQGVFEYYKSVTFVSNCGSHPSTTSKGSPINTQYVF", "P32781");
            DigestionParams nParams = new("trypsin", 2, 7, 50, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.N);
            DigestionParams cParams = new("trypsin", 2, 7, 50, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.C);

            //Unit test is not crashing
            protein.Digest(nParams, null, null).ToList();
            protein.Digest(cParams, null, null).ToList();

            Protein protein2 = new Protein("MQFSTVASVAFVALANFVAAESAAAISQITDGQIQATTTATTEATTTAAPSSTVETVSPSSTETISQQTENGAAKAAVGMGAGALAAAAMLL", "P43497");
            protein2.Digest(nParams, null, null).ToList();
            protein2.Digest(cParams, null, null).ToList();

            List<TruncationProduct> proteolysisProducts = new List<TruncationProduct>
            {
                new TruncationProduct(5, 25, "asdf")
            };

            //speedy
            Protein protein3 = new Protein("MQFSTVASVAFVALANFVAAESAAAISQITDGQIQATTTATTEATTTAAPSSTVETVSPSSTETISQQTENGAAKAAVGMGAGALAAAAMLL", "P43497", proteolysisProducts: proteolysisProducts);
            protein3.Digest(nParams, null, null).ToList();
            protein3.Digest(cParams, null, null).ToList();
            cParams = new DigestionParams("trypsin", 0, 7, 9, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.C, initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave);
            protein3.Digest(cParams, null, null).ToList();
        }

        /// <summary>
        /// CRITICAL: Tests that MaxLength constraints are respected in semi/non-specific digestion.
        /// Without proper length enforcement, search space would explode and include peptides
        /// that cannot physically be detected, wasting computational resources and increasing FDR.
        /// </summary>
        [Test]
        public static void TestSpeedyNonAndSemiSpecificMaxLength()
        {
            Protein Q07065 = new Protein("MPSAKQRGSKGGHGAASPSEKGAHPSGGADDV" +
                "AKKPPPAPQQPPPPPAPHPQQHPQQHPQNQAHGKGGHRGGGGGGGKSSSSSSASAAAAAA" +
                "AASSSASCSRRLGRALNFLFYLALVAAAAFSGWCVHHVLEEVQQVRRSHQDFSRQREELGQ" +
                "GLQGVEQKVQSLQATFGTFESILRSSQHKQDLTEKAVKQGESEVSRISEVLQKLQNEILKDL" +
                "SDGIHVVKDARERDFTSLENTVEERLTELTKSINDNIAIFTEVQKRSQKEINDMKAKVASLEE" +
                "SEGNKQDLKALKEAVKEIQTSAKSREWDMEALRSTLQTMESDIYTEVRELVSLKQEQQAFKEA" +
                "ADTERLALQALTEKLLRSEESVSRLPEEIRRLEEELRQLKSDSHGPKEDGGFRHSEAFEALQQK" +
                "SQGLDSRLQHVEDGVLSMQVASARQTESLESLLSKSQEHEQRLAALQGRLEGLGSSEADQDGLAST" +
                "VRSLGETQLVLYGDVEELKRSVGELPSTVESLQKVQEQVHTLLSQDQAQAARLPPQDFLDRLSSLD" +
                "NLKASVSQVEADLKMLRTAVDSLVAYSVKIETNENNLESAKGLLDDLRNDLDRLFVKVEKIHEKV", "Q07065");

            //Semi
            DigestionParams semiNParams = new DigestionParams("Asp-N", 3, 7, 50, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.N);
            DigestionParams semiCParams = new DigestionParams("Asp-N", 3, 7, 50, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.C);
            List<PeptideWithSetModifications> nPwsms = Q07065.Digest(semiNParams, null, null).ToList();
            List<PeptideWithSetModifications> cPwsms = Q07065.Digest(semiCParams, null, null).ToList();
            Assert.IsFalse(nPwsms.Any(x => x.Length > semiNParams.MaxLength));
            Assert.IsFalse(cPwsms.Any(x => x.Length > semiCParams.MaxLength));
            Assert.IsTrue(nPwsms.Any(x => x.Length == semiNParams.MaxPeptideLength));
            Assert.IsTrue(cPwsms.Any(x => x.Length == semiCParams.MaxPeptideLength));

            //Non
            DigestionParams nonNParams = new DigestionParams("Asp-N", 20, 7, 50, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N); //more missed cleavages here so we can test the end
            DigestionParams nonCParams = new DigestionParams("Asp-N", 3, 7, 50, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            nPwsms = Q07065.Digest(nonNParams, null, null).ToList();
            cPwsms = Q07065.Digest(nonCParams, null, null).ToList();
            Assert.IsFalse(nPwsms.Any(x => x.Length > nonNParams.MaxLength));
            Assert.IsFalse(cPwsms.Any(x => x.Length > nonCParams.MaxLength));
            Assert.IsTrue(nPwsms.Any(x => x.Length == nonNParams.MaxLength));
            Assert.IsTrue(cPwsms.Any(x => x.Length == nonCParams.MaxLength));
            Assert.IsTrue(nPwsms.Any(x => x.Length == nonNParams.MinPeptideLength));
            Assert.IsTrue(cPwsms.Any(x => x.Length == nonCParams.MinPeptideLength));
        }

        /// <summary>
        /// CRITICAL: Comprehensive test for non-specific and semi-specific digestion correctness.
        /// Tests peptide generation, cleavage specificity labeling, duplicate prevention, and
        /// methionine behavior. Essential for all non-tryptic and semi-tryptic searches.
        /// </summary>
        [Test]
        public static void TestNonAndSemiSpecificDigests()
        {
            Protein fiveCleavages = new Protein("MAAKCCKDDKEEKFFKGG", "fiveCleavages"); //protein with 5 K's
            List<Tuple<string, FragmentationTerminus>> trypticSequencesInducingClevage = new List<Tuple<string, FragmentationTerminus>>
            {
                new Tuple<string, FragmentationTerminus>("K",FragmentationTerminus.C ) //cleave at C terminus of K
            };
            List<Tuple<string, FragmentationTerminus>> trypticSequencesPreventingClevage = new List<Tuple<string, FragmentationTerminus>>();

            //make two identical proteases, but one is fully specific and one is semi specific

            var motifList = DigestionMotif.ParseDigestionMotifsFromString("K|");
            Protease trypsinForTestNonAndSemiSpecificDigests = new Protease("trypsinForTestNonAndSemiSpecificDigests", CleavageSpecificity.Full, "asdf", "asdf", motifList);
            Protease semiTrypsinForTestNonAndSemiSpecificDigests = new Protease("semitrypsinForTestNonAndSemiSpecificDigests", CleavageSpecificity.Semi, "asdf", "asdf", motifList);

            //add these made up proteases to the dictionary
            ProteaseDictionary.Dictionary.Add(trypsinForTestNonAndSemiSpecificDigests.Name, trypsinForTestNonAndSemiSpecificDigests);
            ProteaseDictionary.Dictionary.Add(semiTrypsinForTestNonAndSemiSpecificDigests.Name, semiTrypsinForTestNonAndSemiSpecificDigests);

            //Digest with the full
            DigestionParams fullyDigestParams = new DigestionParams(trypsinForTestNonAndSemiSpecificDigests.Name, 3, 2);
            var fiveCleavageProductsTrypsin = fiveCleavages.Digest(fullyDigestParams, null, null).ToList();
            Assert.AreEqual(22, fiveCleavageProductsTrypsin.Count);

            //digests with the semi (variable methionine)
            DigestionParams semiDigestionParams = new DigestionParams(semiTrypsinForTestNonAndSemiSpecificDigests.Name, 3, 2);
            var fiveCleavageProductsSemiTrypsin = fiveCleavages.Digest(semiDigestionParams, null, null).ToList();

            //This is a partial list of the full peptides. From this, we can GENERATE every possible semi that we would expect to see
            List<string> expectedProductsSemiFiveCleavages = new List<string>
            {
                "MAAKCCKDDKEEK",
                "AAKCCKDDKEEK",
                "CCKDDKEEKFFK",
                "DDKEEKFFKGG",
                "EEKFFKGG",
                "FFKGG",
                "AAK",
                "AAKCCK",
                "AAKCCKDDK",
            };

            //Check that, when we digested with semi, we made all possible semi sequences, labeled full and semi correctly, and have no duplicates
            foreach (string s in expectedProductsSemiFiveCleavages) //foreach precursor peptide
            {
                for (int i = 0; i < s.Length - semiDigestionParams.MinLength; i++) //cleave it to be semi
                {
                    string sToFind = s.Substring(i); //get a peptide from this precursor (fixed C)
                    var peps = fiveCleavageProductsSemiTrypsin.Where(x => x.BaseSequence.Equals(sToFind)).ToArray(); //find the peptide in the digested list
                    Assert.IsTrue(peps.Length == 1); //There should be exactly one! More than that means there are duplicates, fewer means we didn't generate it!
                    var pep = peps[0]; //get that single peptide
                    //if it's a full sequence (cleaved at both indexes (including termini))
                    if ((pep.BaseSequence[0] == pep.BaseSequence[1] || pep.BaseSequence[0] == 'M') && (pep.BaseSequence[pep.BaseSequence.Length - 1] == 'K'
                        || (pep.BaseSequence[pep.BaseSequence.Length - 1] == 'G' && pep.BaseSequence[pep.BaseSequence.Length - 2] == 'G')))
                    {
                        Assert.IsTrue(pep.CleavageSpecificityForFdrCategory == CleavageSpecificity.Full);
                    }
                    else
                    {
                        Assert.IsTrue(pep.CleavageSpecificityForFdrCategory == CleavageSpecificity.Semi);
                    }
                    //try to remake the pwsm with unknown specificity... was it assigned the correct specificity?
                    PeptideWithSetModifications pwsmRemake = new PeptideWithSetModifications(fiveCleavages, semiDigestionParams, pep.OneBasedStartResidueInProtein, pep.OneBasedEndResidueInProtein, CleavageSpecificity.Unknown, "", 3, pep.AllModsOneIsNterminus, 0);
                    Assert.IsTrue(pwsmRemake.CleavageSpecificityForFdrCategory == pep.CleavageSpecificityForFdrCategory);

                    //Repeat the above going from the other direction (fixed N)
                    sToFind = s.Substring(0, semiDigestionParams.MinLength + i); //get a peptide from this precursor (fixed N)
                    peps = fiveCleavageProductsSemiTrypsin.Where(x => x.BaseSequence.Equals(sToFind)).ToArray();//find the peptide in the digested list
                    Assert.IsTrue(peps.Length == 1);//There should be exactly one! More than that means there are duplicates, fewer means we didn't generate it!
                    pep = peps[0];//get that single peptide
                    //if it's a full sequence (cleaved at both indexes (including termini))
                    if ((pep.BaseSequence[0] == pep.BaseSequence[1] || pep.BaseSequence[0] == 'M') && pep.BaseSequence.Last() == 'K')
                    {
                        Assert.IsTrue(pep.CleavageSpecificityForFdrCategory == CleavageSpecificity.Full);
                    }
                    else
                    {
                        Assert.IsTrue(pep.CleavageSpecificityForFdrCategory == CleavageSpecificity.Semi);
                    }
                    //try to remake the pwsm with unknown specificity... was it assigned the correct specificity?
                    pwsmRemake = new PeptideWithSetModifications(fiveCleavages, semiDigestionParams, pep.OneBasedStartResidueInProtein, pep.OneBasedEndResidueInProtein, CleavageSpecificity.Unknown, "", 3, pep.AllModsOneIsNterminus, 0);
                    Assert.IsTrue(pwsmRemake.CleavageSpecificityForFdrCategory == pep.CleavageSpecificityForFdrCategory);
                }
            }
            //confirm there were 85 peptides generated by the semi
            Assert.AreEqual(85, fiveCleavageProductsSemiTrypsin.Count);

            //The rest of the tests are less intense
            //check semi when methionine is cleaved
            DigestionParams semiCleaveDigestionParams = new DigestionParams(semiTrypsinForTestNonAndSemiSpecificDigests.Name, 3, 2, initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave);
            var fiveCleavageProductsSemiTrypsinCleave = fiveCleavages.Digest(semiCleaveDigestionParams, null, null).ToList();
            int numVariableWithMet = fiveCleavageProductsSemiTrypsin.Where(x => x.BaseSequence[0] == 'M').Count(); //how many had methionine in the variable digestion?
            Assert.AreEqual(fiveCleavageProductsSemiTrypsin.Count, fiveCleavageProductsSemiTrypsinCleave.Count + numVariableWithMet); //there should be the same number of sequences as before, minus the amount of methionine peptides

            //check semi when methionine is retained
            DigestionParams semiRetainDigestionParams = new DigestionParams(semiTrypsinForTestNonAndSemiSpecificDigests.Name, 3, 2, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var fiveCleavageProductsSemiTrypsinRetain = fiveCleavages.Digest(semiRetainDigestionParams, null, null).ToList();
            int numNotRetained = fiveCleavageProductsSemiTrypsin.Where(x => x.BaseSequence[0] == 'A' && x.BaseSequence[1] == 'A' //how many had Met cleaved in the variable digestion?
            && (x.BaseSequence[x.BaseSequence.Length - 1] != 'K' && !(x.BaseSequence[x.BaseSequence.Length - 1] == 'G' && x.BaseSequence[x.BaseSequence.Length - 2] == 'G'))).Count();
            Assert.AreEqual(fiveCleavageProductsSemiTrypsinRetain.Count + numNotRetained, fiveCleavageProductsSemiTrypsin.Count); //there should be the same number of sequences as before, minus the amount of cleaved peptides

            //Check the speedy semi-specific search (the previous ones were the slow classic)
            //Fixed N
            DigestionParams modernSemiDigestionParamsN = new DigestionParams(trypsinForTestNonAndSemiSpecificDigests.Name, 3, 2, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.N);
            var fiveCleavageProductsModernSemiTrypsinN = fiveCleavages.Digest(modernSemiDigestionParamsN, null, null).ToList();
            Assert.AreEqual(7, fiveCleavageProductsModernSemiTrypsinN.Count);

            //Fixed C
            DigestionParams modernSemiDigestionParamsC = new DigestionParams(trypsinForTestNonAndSemiSpecificDigests.Name, 3, 2, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.C);
            var fiveCleavageProductsModernSemiTrypsinC = fiveCleavages.Digest(modernSemiDigestionParamsC, null, null).ToList();
            Assert.AreEqual(6, fiveCleavageProductsModernSemiTrypsinC.Count);

            //test the maxPeptideLength for both singleN and SingleC (variable methionine)
            //Single N max peptide length
            var modernNonSpecificN = new DigestionParams("singleN", 4, 2, 4, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N);
            var fiveCleavageProductsModernNonSpecificN = fiveCleavages.Digest(modernNonSpecificN, null, null).ToList();
            Assert.AreEqual(17, fiveCleavageProductsModernNonSpecificN.Count);
            foreach (var pep in fiveCleavageProductsModernNonSpecificN)
            {
                Assert.IsTrue(pep.BaseSequence.Length <= 4 && pep.BaseSequence.Length >= 2);
            }

            //Single C max peptide length
            var modernNonSpecificC = new DigestionParams("singleC", 4, 2, 4, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            var fiveCleavageProductsModernNonSpecificC = fiveCleavages.Digest(modernNonSpecificC, null, null).ToList();
            Assert.AreEqual(17, fiveCleavageProductsModernNonSpecificC.Count);
            foreach (var pep in fiveCleavageProductsModernNonSpecificC)
            {
                Assert.IsTrue(pep.BaseSequence.Length <= 4 && pep.BaseSequence.Length >= 2);
            }

            //test speedy nonspecific with variable methionine
            TestSingleProteases(fiveCleavages, InitiatorMethionineBehavior.Variable, FragmentationTerminus.N, 17);
            TestSingleProteases(fiveCleavages, InitiatorMethionineBehavior.Variable, FragmentationTerminus.C, 17);

            //test speedy nonspecific with cleaved methionine
            TestSingleProteases(fiveCleavages, InitiatorMethionineBehavior.Cleave, FragmentationTerminus.N, 16);
            TestSingleProteases(fiveCleavages, InitiatorMethionineBehavior.Cleave, FragmentationTerminus.C, 16);

            //test speedy nonspecific with retained methionine
            TestSingleProteases(fiveCleavages, InitiatorMethionineBehavior.Retain, FragmentationTerminus.N, 17);
            TestSingleProteases(fiveCleavages, InitiatorMethionineBehavior.Retain, FragmentationTerminus.C, 17);

            //test classic nonspecific
            DigestionParams classicNonspecificDigest = new DigestionParams("non-specific", 50);
            List<PeptideWithSetModifications> classicNonspecificPeptides = fiveCleavages.Digest(classicNonspecificDigest, null, null).ToList();
            Assert.IsTrue(classicNonspecificPeptides.Count == 78);
        }

        private static void TestSingleProteases(Protein protein, InitiatorMethionineBehavior initiatorMethionineBehavior, FragmentationTerminus fragmentationTerminus, int numSequencesExpected)
        {
            string protease = FragmentationTerminus.N == fragmentationTerminus ? "singleN" : "singleC";
            DigestionParams digestionParams = new DigestionParams(protease, 50, 2, searchModeType: CleavageSpecificity.None, initiatorMethionineBehavior: initiatorMethionineBehavior, fragmentationTerminus: fragmentationTerminus);
            var products = protein.Digest(digestionParams, null, null).ToList();
            Assert.AreEqual(numSequencesExpected, products.Count);
        }

        /// <summary>
        /// CRITICAL: Tests that terminal modifications are correctly handled with single proteases.
        /// For SingleN/SingleC searches, terminal mods that don't affect fragment ions should be
        /// excluded to prevent redundant peptide generation. Essential for spectral library searches.
        /// </summary>
        [Test]
        public static void TestSingleProteasesWithTerminalMods()
        {
            //we actually don't want C-terminal mods on SingleN or N-terminal mods on SingleC, because they don't influence the fragment ion series and just create a redundant peptide
            //the modified peptides are found using precursor mass matching after scoring, but that happens downstream in MetaMorpheus
            ModificationMotif.TryGetMotif("A", out ModificationMotif motif);
            Modification nTermMod = new Modification(_originalId: "acetylation", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "N-terminal.");
            Modification cTermMod = new Modification(_originalId: "amide", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "C-terminal.");

            Protein proteinWithMods = new Protein("MAGIAAKLAKDREAAEGLGSHA", "testProtein",
                oneBasedModifications: new Dictionary<int, List<Modification>>
                {
                    { 2, new List<Modification>{nTermMod } },
                    { 22, new List<Modification>{cTermMod } }
                });

            DigestionParams singleN = new DigestionParams(protease: "singleN", searchModeType: CleavageSpecificity.SingleN, fragmentationTerminus: FragmentationTerminus.N);
            DigestionParams singleC = new DigestionParams(protease: "singleC", searchModeType: CleavageSpecificity.SingleC, fragmentationTerminus: FragmentationTerminus.C);

            List<Modification> empty = new List<Modification>();
            List<Modification> allMods = new List<Modification> { nTermMod, cTermMod };
            List<PeptideWithSetModifications> nPeps = proteinWithMods.Digest(singleN, empty, empty).ToList();
            List<PeptideWithSetModifications> cPeps = proteinWithMods.Digest(singleC, empty, empty).ToList();
            Assert.IsTrue(nPeps.Count == cPeps.Count);
            Assert.IsTrue(cPeps.Count == 17);

            Protein proteinWithoutMods = new Protein("MAGIAAKLAKDREAAEGLGSHA", "testProtein");

            //Test that variable mods are removed
            nPeps = proteinWithoutMods.Digest(singleN, empty, allMods).ToList();
            cPeps = proteinWithoutMods.Digest(singleC, empty, allMods).ToList();
            Assert.IsTrue(nPeps.Count == cPeps.Count);
            Assert.IsTrue(cPeps.Count == 17);

            //Test that fixed mods are NOT removed
            nPeps = proteinWithoutMods.Digest(singleN, allMods, empty).ToList();
            cPeps = proteinWithoutMods.Digest(singleC, allMods, empty).ToList();
            Assert.IsTrue(nPeps.Count == cPeps.Count);
            Assert.IsTrue(nPeps.All(x => x.FullSequence.Contains("testModType:amide on A")));
            Assert.IsTrue(nPeps.Last().FullSequence.Contains("testModType:amide on A"));
            Assert.IsTrue(cPeps.Count == 16);

            //Test single proteases with specific protease
            DigestionParams specificNonN = new DigestionParams(protease: "Asp-N", searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N);
            DigestionParams specificNonC = new DigestionParams(protease: "Asp-N", searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            List<PeptideWithSetModifications> nSpecificPeps = proteinWithMods.Digest(specificNonN, empty, empty).ToList();
            List<PeptideWithSetModifications> cSpecificPeps = proteinWithMods.Digest(specificNonC, empty, empty).ToList();
            Assert.IsTrue(nSpecificPeps.Count == cSpecificPeps.Count);
            Assert.IsTrue(cSpecificPeps.Count == 17);

            //try again with no missed cleavages
            specificNonN = new DigestionParams(protease: "Asp-N", maxMissedCleavages: 0, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N);
            specificNonC = new DigestionParams(protease: "Asp-N", maxMissedCleavages: 0, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            nSpecificPeps = proteinWithMods.Digest(specificNonN, empty, empty).ToList();
            cSpecificPeps = proteinWithMods.Digest(specificNonC, empty, empty).ToList();
            Assert.IsTrue(nSpecificPeps.Count == 11);
            Assert.IsTrue(cSpecificPeps.Count == 11);
        }

        /// <summary>
        /// CRITICAL: Tests single proteases (N/C terminus fixed) with specific cleavage enzymes.
        /// Validates peptide counts and methionine cleavage behavior for non-specific searches
        /// that still respect enzyme cleavage sites. Essential for hybrid search strategies.
        /// </summary>
        [Test]
        public static void TestSingleProteasesWithSpecificProteases()
        {
            Protein tinyProteinWithCleavages = new Protein("ACDREFGHIKLMNPQRST", "tiny");
            Protein tinyProteinWithoutCleavages = new Protein("ACDEFGHILMNPQST", "tinier");
            Protein bigProteinWithStretchOfNoCleavages = new Protein("ACDREFGHIKLMNPQRSTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGACDREFGHIKLMNPQRST", "big");

            DigestionParams dpN = new DigestionParams("trypsin", 1, 5, 20, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N);
            DigestionParams dpC = new DigestionParams("trypsin", 1, 5, 20, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            List<Modification> empty = new List<Modification>();

            //SingleN tests

            List<PeptideWithSetModifications> peptides = tinyProteinWithCleavages.Digest(dpN, empty, empty).ToList();
            Assert.IsTrue(peptides.Count == 14);
            peptides = tinyProteinWithoutCleavages.Digest(dpN, empty, empty).ToList();
            Assert.IsTrue(peptides.Count == 11);
            peptides = bigProteinWithStretchOfNoCleavages.Digest(dpN, empty, empty).ToList();
            Assert.IsTrue(peptides.Count == 63);

            //SingleC tests
            peptides = tinyProteinWithCleavages.Digest(dpC, empty, empty).ToList();
            Assert.IsTrue(peptides.Count == 14);
            peptides = tinyProteinWithoutCleavages.Digest(dpC, empty, empty).ToList();
            Assert.IsTrue(peptides.Count == 11);
            peptides = bigProteinWithStretchOfNoCleavages.Digest(dpC, empty, empty).ToList();
            Assert.IsTrue(peptides.Count == 63);

            //Methionine cleavage fringe test
            Protein methionineProtein = new Protein("MDBCEFGDHIKLMNODPQRST", "tiny");
            dpN = new DigestionParams("Asp-N", 2, 5, 20, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N, initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave);
            dpC = new DigestionParams("Asp-N", 2, 5, 20, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C, initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave);
            peptides = methionineProtein.Digest(dpN, empty, empty).ToList();
            Assert.IsTrue(peptides.Count == 16);
            peptides = methionineProtein.Digest(dpC, empty, empty).ToList();
            Assert.IsTrue(peptides.Count == 16);
        }

        /// <summary>
        /// EDGE CASE: Tests single proteases on small proteins where peptide count
        /// should equal protein length minus minimum peptide length plus one.
        /// Guards against off-by-one errors in peptide enumeration.
        /// </summary>
        [Test]
        public static void TestSingleProteasesTinyProtein()
        {
            Protein P56381 = new Protein("MVAYWRQAGLSYIRYSQICAKAVRDALKTEFKANAEKTSGSNVKIVKVKKE", "P56381");
            DigestionParams singleN = new DigestionParams(protease: "Asp-N", maxMissedCleavages: 3, minPeptideLength: 7, maxPeptideLength: 50, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N);
            DigestionParams singleC = new DigestionParams(protease: "Asp-N", maxMissedCleavages: 3, minPeptideLength: 7, maxPeptideLength: 50, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            List<PeptideWithSetModifications> nPwsms = P56381.Digest(singleN, null, null).ToList();
            List<PeptideWithSetModifications> cPwsms = P56381.Digest(singleC, null, null).ToList();
            Assert.IsTrue(nPwsms.Count == cPwsms.Count);
            Assert.IsTrue(nPwsms.Count == P56381.Length - singleN.MinLength + 1);
        }

        /// <summary>
        /// CRITICAL: Tests cleavage specificity assignment for FDR categorization.
        /// Peptides must be correctly labeled as Full/Semi/None specificity for proper
        /// FDR estimation. Incorrect labeling leads to inflated or deflated FDR.
        /// </summary>
        [Test]
        public static void TestUpdateCleavageSpecificity()
        {
            Protein protein = new Protein("MACDEFGHIKLMNPQRST", "test");
            DigestionParams dpVariable = new DigestionParams();
            DigestionParams dpRetain = new DigestionParams(initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            Dictionary<int, Modification> empty = new Dictionary<int, Modification>();

            //Test with varying Methionine
            PeptideWithSetModifications fullCleavageVariableMet = new PeptideWithSetModifications(protein, dpVariable, 2, 10, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(fullCleavageVariableMet.CleavageSpecificityForFdrCategory == CleavageSpecificity.Full);
            PeptideWithSetModifications fullCleavageRetainMet = new PeptideWithSetModifications(protein, dpRetain, 2, 10, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(fullCleavageRetainMet.CleavageSpecificityForFdrCategory == CleavageSpecificity.Semi);
            PeptideWithSetModifications semiCleavageVariableMet = new PeptideWithSetModifications(protein, dpVariable, 2, 9, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(semiCleavageVariableMet.CleavageSpecificityForFdrCategory == CleavageSpecificity.Semi);
            PeptideWithSetModifications semiCleavageRetainMet = new PeptideWithSetModifications(protein, dpRetain, 2, 9, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(semiCleavageRetainMet.CleavageSpecificityForFdrCategory == CleavageSpecificity.None);
            PeptideWithSetModifications noneCleavageVariableMet = new PeptideWithSetModifications(protein, dpVariable, 3, 9, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(noneCleavageVariableMet.CleavageSpecificityForFdrCategory == CleavageSpecificity.None);
            PeptideWithSetModifications noneCleavageRetainMet = new PeptideWithSetModifications(protein, dpRetain, 3, 9, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(noneCleavageRetainMet.CleavageSpecificityForFdrCategory == CleavageSpecificity.None);

            //Test with proteolytic cleavages
            protein = new Protein("MACDEFGHIKLMNPQRST", "test", proteolysisProducts: new List<TruncationProduct> { new TruncationProduct(3, 9, "chain") });
            PeptideWithSetModifications fullProteolytic = new PeptideWithSetModifications(protein, dpVariable, 3, 9, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(fullProteolytic.CleavageSpecificityForFdrCategory == CleavageSpecificity.Full);
            fullProteolytic = new PeptideWithSetModifications(protein, dpVariable, 3, 10, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(fullProteolytic.CleavageSpecificityForFdrCategory == CleavageSpecificity.Full);
            PeptideWithSetModifications semiProteolytic = new PeptideWithSetModifications(protein, dpVariable, 3, 6, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(semiProteolytic.CleavageSpecificityForFdrCategory == CleavageSpecificity.Semi);
            semiProteolytic = new PeptideWithSetModifications(protein, dpVariable, 5, 9, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(semiProteolytic.CleavageSpecificityForFdrCategory == CleavageSpecificity.Semi);
        }

        /// <summary>
        /// CRITICAL: Tests splice site intersection logic for variant peptide analysis.
        /// Correctly identifying which peptides span splice sites is essential for
        /// alternative splicing detection and isoform-specific peptide annotation.
        /// </summary>
        [Test]
        public static void TestIncludeSpliceSiteRanges()
        {
            Protein protein = new Protein("MACDEFGHIKLMNPQRST", "test");
            PeptideWithSetModifications pepe = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 10, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);
            SpliceSite ss1Before = new SpliceSite(1, 1, "");
            SpliceSite ss2BeginningBefore = new SpliceSite(1, 2, "");
            SpliceSite ss3 = new SpliceSite(2, 3, "");
            SpliceSite ss4 = new SpliceSite(3, 4, "");
            SpliceSite ss5 = new SpliceSite(9, 10, "");
            SpliceSite ss6EndAfter = new SpliceSite(10, 11, "");
            SpliceSite ss7After = new SpliceSite(11, 12, "");
            Assert.IsFalse(pepe.IncludesSpliceSite(ss1Before));
            Assert.IsFalse(pepe.IncludesSpliceSite(ss2BeginningBefore));
            Assert.IsTrue(pepe.IncludesSpliceSite(ss3));
            Assert.IsTrue(pepe.IncludesSpliceSite(ss4));
            Assert.IsTrue(pepe.IncludesSpliceSite(ss5));
            Assert.IsFalse(pepe.IncludesSpliceSite(ss6EndAfter));
            Assert.IsFalse(pepe.IncludesSpliceSite(ss7After));
        }

        /// <summary>
        /// CRITICAL: Tests sequence variation intersection and identification logic.
        /// Validates both 'intersects' (peptide overlaps variant) and 'identifies'
        /// (peptide uniquely confirms variant) for SAV, MNV, indels, and truncations.
        /// Essential for variant peptide identification in proteogenomics.
        /// </summary>
        [Test]
        public static void TestIntersectsSequenceVariations()
        {
            Protein protein = new Protein("MACDEFGHIK", "test");
            PeptideWithSetModifications pepe = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 10, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);

            // The weird thing here is that IntersectsWithVariation takes in applied variations,
            // so these are constructed as if already applied
            SequenceVariation sv1Before = new SequenceVariation(1, 1, "A", "M", ""); // before peptide (not identified)
            SequenceVariation sv2Synonymous = new SequenceVariation(2, 2, "A", "A", ""); // no change (intersects because peptide crosses entire variant but is not truly "identified")
            SequenceVariation sv4MissenseBeginning = new SequenceVariation(2, 2, "V", "A", ""); // missense at beginning
            SequenceVariation sv5InsertionAtEnd = new SequenceVariation(7, 9, "GHI", "GHIK", ""); // insertion or stop loss
            SequenceVariation sv6Deletion = new SequenceVariation(2, 3, "AC", "A", ""); // deletion
            SequenceVariation sv66Truncation = new SequenceVariation(10, 20, "KAAAAAAAAAA", "K", ""); // truncation or stop gain (identified because peptide crosses entire variant)
            SequenceVariation sv7MNP = new SequenceVariation(2, 3, "AA", "AC", ""); // mnp
            SequenceVariation sv77MNP = new SequenceVariation(2, 3, "AC", "AC", ""); // synonymous mnp (identified because peptide crosses entire variant)
            SequenceVariation sv9MissenseInRange = new SequenceVariation(3, 3, "C", "V", ""); // missense in range
            SequenceVariation sv10MissenseRangeEdge = new SequenceVariation(10, 10, "K", "R", ""); // missense at end
            SequenceVariation sv11After = new SequenceVariation(11, 11, "L", "V", ""); // after peptide (not identified)

            Assert.IsFalse(pepe.IntersectsAndIdentifiesVariation(sv1Before).intersects);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv2Synonymous).intersects);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv4MissenseBeginning).intersects);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv5InsertionAtEnd).intersects);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv6Deletion).intersects);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv66Truncation).intersects);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv7MNP).intersects);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv77MNP).intersects);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv9MissenseInRange).intersects);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv10MissenseRangeEdge).intersects);
            Assert.IsFalse(pepe.IntersectsAndIdentifiesVariation(sv11After).intersects);

            Assert.IsFalse(pepe.IntersectsAndIdentifiesVariation(sv1Before).identifies);
            Assert.IsFalse(pepe.IntersectsAndIdentifiesVariation(sv2Synonymous).identifies);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv4MissenseBeginning).identifies);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv5InsertionAtEnd).identifies);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv6Deletion).identifies);
            Assert.IsFalse(pepe.IntersectsAndIdentifiesVariation(sv66Truncation).identifies);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv7MNP).identifies);
            Assert.IsFalse(pepe.IntersectsAndIdentifiesVariation(sv77MNP).identifies);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv9MissenseInRange).identifies);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv10MissenseRangeEdge).identifies);
            Assert.IsFalse(pepe.IntersectsAndIdentifiesVariation(sv11After).identifies);

            PeptideWithSetModifications pepe2 = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);
            Assert.IsTrue(pepe2.IntersectsAndIdentifiesVariation(sv5InsertionAtEnd).intersects); // this only intersects GHI, which is the same in GHI -> GHIK
            Assert.IsFalse(pepe2.IntersectsAndIdentifiesVariation(sv5InsertionAtEnd).identifies);
        }

        /// <summary>
        /// CRITICAL: Tests detection of variant-containing peptides.
        /// The IsVariantPeptide() method must correctly distinguish peptides that
        /// span applied sequence variations from those that don't. Essential for
        /// filtering and reporting variant peptides.
        /// </summary>
        [Test]
        public static void TestIsVariantPeptide()
        {
            Protein protein = new Protein("MPEPTIDENEWPEPTIDE", "protein0", appliedSequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "V", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) });

            PeptideWithSetModifications pepe = new PeptideWithSetModifications(protein, new DigestionParams(), 1, 8, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications notPepe = new PeptideWithSetModifications(protein, new DigestionParams(), 9, 18, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);

            Assert.IsTrue(pepe.IsVariantPeptide());
            Assert.IsFalse(notPepe.IsVariantPeptide());
        }

        /// <summary>
        /// EDGE CASE: Tests peptide behavior when parent protein is null or when
        /// accessing flanking residues. PreviousResidue/NextResidue should return '-'
        /// at protein termini or when parent is null. Important for de novo peptides.
        /// </summary>
        [Test]
        public static void TestPeptideWithSetModsNoParentProtein()
        {
            // null parent
            DigestionParams dParams = new DigestionParams();
            var pwsm = new PeptideWithSetModifications("P", null,
                digestionParams: dParams, p: null);
            Assert.AreEqual('-', pwsm.PreviousAminoAcid);
            Assert.AreEqual('-', pwsm.PreviousResidue);
            Assert.AreEqual('-', pwsm.NextAminoAcid);
            Assert.AreEqual('-', pwsm.NextResidue);

            // non-null parent
            Protein protein = new("MQLLRCFSIFSVIASVLAQELTTICEQIPSPTLESTPYSLSTTTILANGKAMQGVFEYYKSVTFVSNCGSHPSTTSKGSPINTQYVF", "P32781");
            var pwsMods = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).ToList();

            var first = pwsMods.First(p => p.BaseSequence == "MQLLRCFSIFSVIASVLAQELTTICEQIPSPTLESTPYSLSTTTILANGK");
            Assert.AreEqual('-', first.PreviousAminoAcid);
            Assert.AreEqual('-', first.PreviousResidue);
            Assert.AreEqual('A', first.NextAminoAcid);
            Assert.AreEqual('A', first.NextResidue);

            var middle = pwsMods.First(p => p.BaseSequence == "SVTFVSNCGSHPSTTSK");
            Assert.AreEqual('K', middle.PreviousAminoAcid);
            Assert.AreEqual('K',middle.PreviousResidue);
            Assert.AreEqual('G',middle.NextAminoAcid);
            Assert.AreEqual('G',middle.NextResidue);

            var last = pwsMods.First(p => p.BaseSequence == "GSPINTQYVF");
            Assert.AreEqual('K', last.PreviousAminoAcid);
            Assert.AreEqual('K', last.PreviousResidue);
            Assert.AreEqual('-', last.NextAminoAcid);
            Assert.AreEqual('-', last.NextResidue);
        }
    }
}
