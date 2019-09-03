using Chemistry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    public static class TestPeptideWithSetMods
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
        /// The purpose of this test is to ensure that two peptides digested from two different proteases are not equal even if their sequences are equal
        /// This is important for multiprotease parsimony in MetaMorpheus
        /// </summary>
        [Test]
        public static void TestDifferentProteaseEquals()
        {
            Protein myProtein = new Protein("SEQUENCEK", "accession");

            DigestionParams digest1 = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            DigestionParams digest2 = new DigestionParams(protease: "Lys-C (cleave before proline)", maxMissedCleavages: 0, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            PeptideWithSetModifications pep1 = myProtein.Digest(digest1, new List<Modification>(), new List<Modification>()).First();
            PeptideWithSetModifications pep2 = myProtein.Digest(digest2, new List<Modification>(), new List<Modification>()).First();

            Assert.That(pep1.FullSequence.Equals(pep2.FullSequence));
            Assert.That(pep1.Protein.Equals(pep2.Protein));
            Assert.That(!pep1.DigestionParams.Protease.Equals(pep2.DigestionParams.Protease));
            Assert.That(!pep1.Equals(pep2));
            Assert.That(!pep1.GetHashCode().Equals(pep2.GetHashCode()));
        }

        [Test]
        public static void TestSemiFewCleavages()
        {
            Protein protein = new Protein("MQLLRCFSIFSVIASVLAQELTTICEQIPSPTLESTPYSLSTTTILANGKAMQGVFEYYKSVTFVSNCGSHPSTTSKGSPINTQYVF", "P32781");
            DigestionParams nParams = new DigestionParams("trypsin", 2, 7, 50, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.N);
            DigestionParams cParams = new DigestionParams("trypsin", 2, 7, 50, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.C);

            //Unit test is not crashing
            protein.Digest(nParams, null, null).ToList();
            protein.Digest(cParams, null, null).ToList();

            Protein protein2 = new Protein("MQFSTVASVAFVALANFVAAESAAAISQITDGQIQATTTATTEATTTAAPSSTVETVSPSSTETISQQTENGAAKAAVGMGAGALAAAAMLL", "P43497");
            protein2.Digest(nParams, null, null).ToList();
            protein2.Digest(cParams, null, null).ToList();

            List<ProteolysisProduct> proteolysisProducts = new List<ProteolysisProduct>
            {
                new ProteolysisProduct(5, 25, "asdf")
            };

            //speedy
            Protein protein3 = new Protein("MQFSTVASVAFVALANFVAAESAAAISQITDGQIQATTTATTEATTTAAPSSTVETVSPSSTETISQQTENGAAKAAVGMGAGALAAAAMLL", "P43497", proteolysisProducts: proteolysisProducts);
            protein3.Digest(nParams, null, null).ToList();
            protein3.Digest(cParams, null, null).ToList();
            cParams = new DigestionParams("trypsin", 0, 7, 9, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.C, initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave);
            protein3.Digest(cParams, null, null).ToList();

            //classic
            DigestionParams classicSemi = new DigestionParams("semi-trypsin", 2, 7, 50);
            protein3.Digest(classicSemi, null, null).ToList();
        }

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
            Assert.IsFalse(nPwsms.Any(x => x.Length > semiNParams.MaxPeptideLength));
            Assert.IsFalse(cPwsms.Any(x => x.Length > semiCParams.MaxPeptideLength));
            Assert.IsTrue(nPwsms.Any(x => x.Length == semiNParams.MaxPeptideLength));
            Assert.IsTrue(cPwsms.Any(x => x.Length == semiCParams.MaxPeptideLength));

            //Non
            DigestionParams nonNParams = new DigestionParams("Asp-N", 20, 7, 50, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N); //more missed cleavages here so we can test the end
            DigestionParams nonCParams = new DigestionParams("Asp-N", 3, 7, 50, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            nPwsms = Q07065.Digest(nonNParams, null, null).ToList();
            cPwsms = Q07065.Digest(nonCParams, null, null).ToList();
            Assert.IsFalse(nPwsms.Any(x => x.Length > nonNParams.MaxPeptideLength));
            Assert.IsFalse(cPwsms.Any(x => x.Length > nonCParams.MaxPeptideLength));
            Assert.IsTrue(nPwsms.Any(x => x.Length == nonNParams.MaxPeptideLength));
            Assert.IsTrue(cPwsms.Any(x => x.Length == nonCParams.MaxPeptideLength));
            Assert.IsTrue(nPwsms.Any(x => x.Length == nonNParams.MinPeptideLength));
            Assert.IsTrue(cPwsms.Any(x => x.Length == nonCParams.MinPeptideLength));
        }

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
                for (int i = 0; i < s.Length - semiDigestionParams.MinPeptideLength; i++) //cleave it to be semi
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
                    sToFind = s.Substring(0, semiDigestionParams.MinPeptideLength + i); //get a peptide from this precursor (fixed N)
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

        [Test]
        public static void TestHardToParseModifiedSequence()
        {
            string fullSequence = "PE[Metal:Cation:Fe[III] on X]PTIDE";

            ModificationMotif.TryGetMotif("X", out var motif);

            Modification mod = new Modification(_originalId: "Cation:Fe[III]", _modificationType: "Metal",
                _monoisotopicMass: 1, _locationRestriction: "Anywhere.", _target: motif);

            Dictionary<string, Modification> mods = new Dictionary<string, Modification> { { "Cation:Fe[III] on X", mod } };

            PeptideWithSetModifications pep = new PeptideWithSetModifications(fullSequence, mods);

            Assert.That(pep.AllModsOneIsNterminus.Count == 1);
            var annotatedMod = pep.AllModsOneIsNterminus.First();
            Assert.That(annotatedMod.Key == 3);
            Assert.That(annotatedMod.Value.IdWithMotif == "Cation:Fe[III] on X");
            Assert.That(annotatedMod.Value.OriginalId == "Cation:Fe[III]");
            Assert.That(annotatedMod.Value.ModificationType == "Metal");

            fullSequence = "[Metal:Cation:Fe[III] on X]PE[Metal:Cation:Fe[III] on X]PTIDE[Metal:Cation:Fe[III] on X]";
            pep = new PeptideWithSetModifications(fullSequence, mods);
            Assert.That(pep.AllModsOneIsNterminus.Count == 3);
            Assert.That(pep.AllModsOneIsNterminus.Keys.ToList().SequenceEqual(new int[] { 1, 3, 8 }));
        }

        [Test]
        public static void TestCTermAndLastSideChainModParsing()
        {
            string fullSequence = "PEPTIDE[Mod:MyMod on E][PeptideCTermMod:MyCTermMod on E]";

            ModificationMotif.TryGetMotif("E", out var motif);

            Modification mod = new Modification(_originalId: "MyMod", _modificationType: "Mod",
                _monoisotopicMass: 1, _locationRestriction: "Anywhere.", _target: motif);

            Modification cTermMod = new Modification(_originalId: "MyCTermMod", _modificationType: "PeptideCTermMod",
                _monoisotopicMass: 1, _locationRestriction: "Peptide C-terminal.", _target: motif);

            Dictionary<string, Modification> mods = new Dictionary<string, Modification>
            {
                { "MyMod on E", mod },
                { "MyCTermMod on E", cTermMod }
            };

            PeptideWithSetModifications pep = new PeptideWithSetModifications(fullSequence, mods);

            Assert.That(pep.AllModsOneIsNterminus.Count == 2);
            Assert.That(pep.AllModsOneIsNterminus.Keys.SequenceEqual(new int[] { 8, 9 }));
        }

        [Test]
        public static void TestPeptideWithSetMod_GetHashCode()
        {
            PeptideWithSetModifications pep1 = new PeptideWithSetModifications("SEQUENCEK", new Dictionary<string, Modification>());
            int oneHashCode = pep1.GetHashCode();

            //if digestion params is not defined, the peptidewithsetmods should still return a hashcode.
            Assert.IsNotNull(oneHashCode);

            Protein myProtein = new Protein("SEQUENCEK", "accession");

            DigestionParams digest1 = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            PeptideWithSetModifications pep2 = myProtein.Digest(digest1, new List<Modification>(), new List<Modification>()).First();

            int twoHashCode = pep2.GetHashCode();

            //if digestion params IS defined, the peptidewithsetmods should  return a hashcode.
            Assert.IsNotNull(twoHashCode);
        }

        [Test]
        public static void TestTopDownDigestion()
        {
            List<ProteolysisProduct> proteolysisProducts = new List<ProteolysisProduct>
            {
                new ProteolysisProduct(5, 20, "asdf")
            };
            Protein protein = new Protein("MACDEFGHIKLMNOPQRSTVWYMACDEFGHIKLMNOPQRSTVWYMACDEFGHIKLMNOPQRSTVWY", "testProtein", "Mus", proteolysisProducts: proteolysisProducts);
            DigestionParams topdownParams = new DigestionParams("top-down");
            List<PeptideWithSetModifications> peptides = protein.Digest(topdownParams, null, null).ToList();
            Assert.IsTrue(peptides.Count == 3);
        }

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
            protein = new Protein("MACDEFGHIKLMNPQRST", "test", proteolysisProducts: new List<ProteolysisProduct> { new ProteolysisProduct(3, 9, "chain") });
            PeptideWithSetModifications fullProteolytic = new PeptideWithSetModifications(protein, dpVariable, 3, 9, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(fullProteolytic.CleavageSpecificityForFdrCategory == CleavageSpecificity.Full);
            fullProteolytic = new PeptideWithSetModifications(protein, dpVariable, 3, 10, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(fullProteolytic.CleavageSpecificityForFdrCategory == CleavageSpecificity.Full);
            PeptideWithSetModifications semiProteolytic = new PeptideWithSetModifications(protein, dpVariable, 3, 6, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(semiProteolytic.CleavageSpecificityForFdrCategory == CleavageSpecificity.Semi);
            semiProteolytic = new PeptideWithSetModifications(protein, dpVariable, 5, 9, CleavageSpecificity.Unknown, "", 0, empty, 0);
            Assert.IsTrue(semiProteolytic.CleavageSpecificityForFdrCategory == CleavageSpecificity.Semi);
        }

        [Test]
        public static void TestSingleProteasesTinyProtein()
        {
            Protein P56381 = new Protein("MVAYWRQAGLSYIRYSQICAKAVRDALKTEFKANAEKTSGSNVKIVKVKKE", "P56381");
            DigestionParams singleN = new DigestionParams(protease: "Asp-N", maxMissedCleavages: 3, minPeptideLength: 7, maxPeptideLength: 50, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N);
            DigestionParams singleC = new DigestionParams(protease: "Asp-N", maxMissedCleavages: 3, minPeptideLength: 7, maxPeptideLength: 50, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            List<PeptideWithSetModifications> nPwsms = P56381.Digest(singleN, null, null).ToList();
            List<PeptideWithSetModifications> cPwsms = P56381.Digest(singleC, null, null).ToList();
            Assert.IsTrue(nPwsms.Count == cPwsms.Count);
            Assert.IsTrue(nPwsms.Count == P56381.Length - singleN.MinPeptideLength + 1);
        }

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

        [Test]
        public static void TestIntersectsSequenceVariations()
        {
            Protein protein = new Protein("MACDEFGHIK", "test");
            PeptideWithSetModifications pepe = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 10, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);

            // The weird thing here is that IntersectsWithVariation takes in applied variations,
            // so these are constructed as if already applied
            SequenceVariation sv1Before = new SequenceVariation(1, 1, "A", "M", ""); // before peptide (not identified)
            SequenceVariation sv2Synonymous = new SequenceVariation(2, 2, "A", "A", ""); // no change (identified because peptide crosses entire variant)
            SequenceVariation sv4MissenseBeginning = new SequenceVariation(2, 2, "V", "A", ""); // missense at beginning
            SequenceVariation sv5InsertionAtEnd = new SequenceVariation(7, 9, "GHI", "GHIK", ""); // insertion or stop loss
            SequenceVariation sv6Deletion = new SequenceVariation(2, 3, "AC", "A", ""); // deletion
            SequenceVariation sv66Truncation = new SequenceVariation(10, 20, "KAAAAAAAAAA", "K", ""); // truncation or stop gain (identified because peptide crosses entire variant)
            SequenceVariation sv7MNP = new SequenceVariation(2, 3, "AA", "AC", ""); // mnp
            SequenceVariation sv77MNP = new SequenceVariation(2, 3, "AC", "AC", ""); // synonymous mnp (identified because peptide crosses entire variant)
            SequenceVariation sv9MissenseInRange = new SequenceVariation(3, 3, "C", "V", ""); // missense in range
            SequenceVariation sv10MissenseRangeEdge = new SequenceVariation(10, 10, "K", "R", ""); // missense at end
            SequenceVariation sv11After = new SequenceVariation(11, 11, "L", "V", ""); // after peptide (not identified)

            Assert.IsFalse(pepe.IntersectsAndIdentifiesVariation(sv1Before).Item2);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv2Synonymous).Item2);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv4MissenseBeginning).Item2);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv5InsertionAtEnd).Item2);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv6Deletion).Item2);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv66Truncation).Item2);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv7MNP).Item2);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv77MNP).Item2);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv9MissenseInRange).Item2);
            Assert.IsTrue(pepe.IntersectsAndIdentifiesVariation(sv10MissenseRangeEdge).Item2);
            Assert.IsFalse(pepe.IntersectsAndIdentifiesVariation(sv11After).Item2);

            PeptideWithSetModifications pepe2 = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);
            Assert.IsFalse(pepe2.IntersectsAndIdentifiesVariation(sv5InsertionAtEnd).Item2); // this only intersects GHI, which is the same in GHI -> GHIK
        }

        [Test]
        public static void TestSeqVarString()
        {
            Protein protein = new Protein("MACDEFGHIK", "test");

            // mod on N-terminus
            PeptideWithSetModifications pepe = new PeptideWithSetModifications(protein, new DigestionParams(), 1, 10, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification> { { 1, new Modification("mod on M", "mod", "mod", "mod") } }, 0);
            SequenceVariation sv1Before = new SequenceVariation(1, 1, "A", "M", ""); // n-terminal mod goes before the sequence
            Assert.AreEqual("A1[mod:mod on M]M", pepe.SequenceVariantString(sv1Before,true));

            // mod in middle
            PeptideWithSetModifications pepe2 = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 10, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification> { { 2, new Modification("mod on A", "mod", "mod", "mod") } }, 0);
            SequenceVariation sv4MissenseBeginning = new SequenceVariation(2, 2, "V", "A", ""); // missense at beginning
            Assert.AreEqual("V2A[mod:mod on A]", pepe2.SequenceVariantString(sv4MissenseBeginning,true));

            // truncated seqvar doesn't truncate in string report (using applied variation correctly)
            PeptideWithSetModifications pepe3 = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);
            SequenceVariation svvvv = new SequenceVariation(7, 10, "GHM", "GHIK", ""); // insertion
            Assert.AreEqual("GHM7GHIK", pepe3.SequenceVariantString(svvvv,true));
        }

        [Test]
        public static void TestIdentifyandStringMethods()
        {
            ModificationMotif.TryGetMotif("V", out ModificationMotif motifV);
            Modification mv = new Modification("mod", null, "type", null, motifV, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            Modification mp = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            Dictionary<int, Modification> modV = new Dictionary<int, Modification>();
            modV.Add(4, mv);
            Dictionary<int, Modification> modP = new Dictionary<int, Modification>();
            modP.Add(5, mp);
           
            List<Protein> proteins = new List<Protein>
            {
                new Protein("MPEPTIDE", "protein0", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "V", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein1", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 5, "PT", "KT", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein2", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "PPP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPPPTIDE", "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "PPP", "P", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPKPKTIDE", "protein4", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 7, "PKPK", "PK", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTAIDE", "protein5",sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "PTA", "KT", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEKKAIDE", "protein6", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "KKA", "K", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein7", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "V", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", new Dictionary<int, List<Modification>> {{ 4, new[] { mv }.ToList() } }) }),
                new Protein("MPEPTIDE", "protein8",sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "PPP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", new Dictionary<int, List<Modification>> {{ 5, new[] { mp }.ToList() } }) }),
                new Protein("MPEPTIDEPEPTIDE", "protein9", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 15, "PTIDEPEPTIDE", "PPP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
            };

            DigestionParams dp = new DigestionParams(minPeptideLength: 2);            
            var protein0_variant = proteins.ElementAt(0).GetVariantProteins().ElementAt(0);
            var protein1_variant = proteins.ElementAt(1).GetVariantProteins().ElementAt(0);
            var protein2_variant = proteins.ElementAt(2).GetVariantProteins().ElementAt(0);
            var protein3_variant = proteins.ElementAt(3).GetVariantProteins().ElementAt(0);
            var protein4_variant = proteins.ElementAt(4).GetVariantProteins().ElementAt(0);
            var protein5_variant = proteins.ElementAt(5).GetVariantProteins().ElementAt(0);
            var protein6_variant = proteins.ElementAt(6).GetVariantProteins().ElementAt(0);
            var protein7_variant = proteins.ElementAt(7).GetVariantProteins().ElementAt(0);
            var protein8_variant = proteins.ElementAt(8).GetVariantProteins().ElementAt(0);
            var protein9_variant = proteins.ElementAt(9).GetVariantProteins().ElementAt(0);

            List<Modification> digestMods = new List<Modification>();

            var protein0_peptide = protein0_variant.Digest(dp, digestMods, digestMods).ElementAt(0);
            var protein1_peptide = protein1_variant.Digest(dp, digestMods, digestMods).ElementAt(2);
            var protein2_peptide = protein2_variant.Digest(dp, digestMods, digestMods).ElementAt(0);
            var protein3_peptide = protein3_variant.Digest(dp, digestMods, digestMods).ElementAt(0);
            var protein4_peptide = protein4_variant.Digest(dp, digestMods, digestMods).ElementAt(2);
            var protein5_peptide = protein5_variant.Digest(dp, digestMods, digestMods).ElementAt(2);
            var protein6_peptide = protein6_variant.Digest(dp, digestMods, digestMods).ElementAt(2);
            var protein7_peptide = protein7_variant.Digest(dp, digestMods, digestMods).ElementAt(1);
            var protein8_peptide = protein8_variant.Digest(dp, digestMods, digestMods).ElementAt(1);
            var protein9_peptide = protein9_variant.Digest(dp, digestMods, digestMods).ElementAt(0);            

            Assert.AreEqual((true,true), protein0_peptide.IntersectsAndIdentifiesVariation(protein0_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein1_peptide.IntersectsAndIdentifiesVariation(protein1_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein2_peptide.IntersectsAndIdentifiesVariation(protein2_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein3_peptide.IntersectsAndIdentifiesVariation(protein3_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((false, false), protein4_peptide.IntersectsAndIdentifiesVariation(protein4_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein5_peptide.IntersectsAndIdentifiesVariation(protein5_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((false, true), protein6_peptide.IntersectsAndIdentifiesVariation(protein6_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein7_peptide.IntersectsAndIdentifiesVariation(protein7_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein8_peptide.IntersectsAndIdentifiesVariation(protein8_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein9_peptide.IntersectsAndIdentifiesVariation(protein9_variant.AppliedSequenceVariations.ElementAt(0)));

            Assert.AreEqual("P4V", protein0_peptide.SequenceVariantString(protein0_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("PT4KT", protein1_peptide.SequenceVariantString(protein1_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("P4PPP", protein2_peptide.SequenceVariantString(protein2_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("PPP4P", protein3_peptide.SequenceVariantString(protein3_variant.AppliedSequenceVariations.ElementAt(0), true));            
            Assert.AreEqual("PTA4KT", protein5_peptide.SequenceVariantString(protein5_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("KKA4K", protein6_peptide.SequenceVariantString(protein6_variant.AppliedSequenceVariations.ElementAt(0), false));
            Assert.AreEqual("P4V[type:mod on V]", protein7_peptide.SequenceVariantString(protein7_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("P4PP[type:mod on P]P", protein8_peptide.SequenceVariantString(protein8_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("PTIDEPEPTIDE4PPP", protein9_peptide.SequenceVariantString(protein9_variant.AppliedSequenceVariations.ElementAt(0), true));
        }        

        [Test]
        public static void BreakDeserializationMethod()
        {
            Assert.Throws<MzLibUtil.MzLibException>(() => new PeptideWithSetModifications("|", new Dictionary<string, Modification>())); // ambiguous
            Assert.Throws<MzLibUtil.MzLibException>(() => new PeptideWithSetModifications("[]", new Dictionary<string, Modification>())); // bad mod
            Assert.Throws<MzLibUtil.MzLibException>(() => new PeptideWithSetModifications("A[:mod]", new Dictionary<string, Modification>())); // nonexistent mod
        }
    }
}