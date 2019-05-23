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
            specificNonN = new DigestionParams(protease: "Asp-N", 0, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N);
            specificNonC = new DigestionParams(protease: "Asp-N", 0, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            nSpecificPeps = proteinWithMods.Digest(specificNonN, empty, empty).ToList();
            cSpecificPeps = proteinWithMods.Digest(specificNonC, empty, empty).ToList();
            Assert.IsTrue(nSpecificPeps.Count == 11);
            Assert.IsTrue(cSpecificPeps.Count == 10);
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
            DigestionParams singleN = new DigestionParams("Asp-N", 3, 7, 50, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N);
            DigestionParams singleC = new DigestionParams("Asp-N", 3, 7, 50, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            List<PeptideWithSetModifications> nPwsms = P56381.Digest(singleN, null, null).ToList();
            List<PeptideWithSetModifications> cPwsms = P56381.Digest(singleC, null, null).ToList();
            Assert.IsTrue(nPwsms.Count == cPwsms.Count);
            Assert.IsTrue(nPwsms.Count == P56381.Length - singleN.MinPeptideLength + 1);
        }
    }
}