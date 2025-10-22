using Chemistry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MzLibUtil;
using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;
using Transcriptomics;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
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
            Assert.That(pep1.Parent.Equals(pep2.Parent));
            Assert.That(!pep1.DigestionParams.DigestionAgent.Equals(pep2.DigestionParams.DigestionAgent));
            Assert.That(!pep1.Equals(pep2));
            Assert.That(!pep1.Equals((object)pep2));
            Assert.That(!pep1.GetHashCode().Equals(pep2.GetHashCode()));
        }

        [Test]
        public static void TestPeptideOligoEquality()
        {
            var oligo = new OligoWithSetMods("GUACUG", []);
            var peptide = new PeptideWithSetModifications("PEPTIDE", []);

            Assert.That(!oligo.Equals(peptide));
            Assert.That(!peptide.Equals(oligo));
            Assert.That(!((IBioPolymerWithSetMods)oligo).Equals(peptide));
            Assert.That(!((IBioPolymerWithSetMods)peptide).Equals(oligo));
            Assert.That(!((object)oligo).Equals(peptide));
            Assert.That(!((object)peptide).Equals(oligo));
        }

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
            string fullSequenceBothMods = "PEPTIDE[Mod:MyMod on E]-[PeptideCTermMod:MyCTermMod on E]";
            string fullSequenceCTermOnly = "PEPTIDE-[PeptideCTermMod:MyCTermMod on E]";
            string fullSequenceSideChainOnly = "PEPTIDE[Mod:MyMod on E]";

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

            PeptideWithSetModifications pepBothMods = new PeptideWithSetModifications(fullSequenceBothMods, mods);
            PeptideWithSetModifications pepCterm = new PeptideWithSetModifications(fullSequenceCTermOnly, mods);
            PeptideWithSetModifications pepSideChain = new PeptideWithSetModifications(fullSequenceSideChainOnly, mods);

            Assert.That(pepBothMods.AllModsOneIsNterminus.Count == 2);
            Assert.That(pepBothMods.AllModsOneIsNterminus.Keys.SequenceEqual(new int[] { 8, 9 }));
            Assert.That(pepBothMods.AllModsOneIsNterminus[8].IdWithMotif == "MyMod on E");
            Assert.That(pepBothMods.AllModsOneIsNterminus[9].IdWithMotif == "MyCTermMod on E");
            Assert.That(pepCterm.AllModsOneIsNterminus.Count == 1);
            Assert.That(pepCterm.AllModsOneIsNterminus.Keys.SequenceEqual(new int[] { 9 }));
            Assert.That(pepCterm.AllModsOneIsNterminus[9].IdWithMotif == "MyCTermMod on E");
            Assert.That(pepSideChain.AllModsOneIsNterminus.Count == 1);
            Assert.That(pepSideChain.AllModsOneIsNterminus.Keys.SequenceEqual(new int[] { 8 }));
            Assert.That(pepSideChain.AllModsOneIsNterminus[8].IdWithMotif == "MyMod on E");
        }

        [Test]
        public static void TestPeptideWithSetMods_GetHashCode()
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
            List<TruncationProduct> proteolysisProducts = new List<TruncationProduct>
            {
                new TruncationProduct(5, 20, "asdf")
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
            // Protein:  M  A  C  D  E  F  G  H  I  K
            // Position: 1  2  3  4  5  6  7  8  9  10
            var protein = new Protein("MACDEFGHIK", "test");

            // Peptide covering residues 2–10 (A..K)
            var pepFull = new PeptideWithSetModifications(
                protein, new DigestionParams(), 2, 10,
                CleavageSpecificity.Unknown, "", 0,
                new Dictionary<int, Modification>(), 0);

            // Shorter peptide (2–9) to exercise non-intersect terminal logic with a downstream stop gain
            var pepShort = new PeptideWithSetModifications(
                protein, new DigestionParams(), 2, 9,
                CleavageSpecificity.Unknown, "", 0,
                new Dictionary<int, Modification>(), 0);

            // 1. Missense BEFORE peptide start (pos 1: M -> A)
            var vBefore = new SequenceVariation(1, 1, "M", "A", "missense_before");
            // 2. Missense AT peptide start (pos 2: A -> V)
            var vBegin = new SequenceVariation(2, 2, "A", "V", "missense_begin");
            // 3. Internal insertion / expansion (pos 5: E -> EQK; expansion length +2)
            var vInsertion = new SequenceVariation(5, 5, "E", "EQK", "insertion_expansion");
            // 4. Internal deletion / contraction (pos 7–8: GH -> G; net -1)
            var vDeletion = new SequenceVariation(7, 8, "GH", "G", "internal_deletion");
            // 5. Stop gain at last residue (pos 10: K -> * )
            var vStopEnd = new SequenceVariation(10, 10, "K", "*", "stop_gain_terminal");
            // 6. Same stop gain evaluated against shorter peptide (should not intersect, but can identify via terminal logic)
            var vStopBeyondShort = vStopEnd; // reuse object

            // Assertions for pepFull (2–10)
            Assert.AreEqual((false, false), pepFull.IntersectsAndIdentifiesVariation(vBefore), "Missense before peptide should neither intersect nor identify.");
            Assert.AreEqual((true, true), pepFull.IntersectsAndIdentifiesVariation(vBegin), "Missense at peptide start should intersect & identify.");
            Assert.AreEqual((true, true), pepFull.IntersectsAndIdentifiesVariation(vInsertion), "Insertion expansion should intersect & identify.");
            Assert.AreEqual((true, true), pepFull.IntersectsAndIdentifiesVariation(vDeletion), "Internal deletion should intersect & identify (length contraction).");
            Assert.AreEqual((true, true), pepFull.IntersectsAndIdentifiesVariation(vStopEnd), "Terminal stop gain inside span should intersect & identify.");

            // Assertions for pepShort (2–9)
            // Stop gain at position 10 is exactly one residue beyond pepShort end (9);
            // Intersects = false, but identification can occur if a new protease site / termination is introduced.
            var shortResult = pepShort.IntersectsAndIdentifiesVariation(vStopBeyondShort);
            Assert.IsFalse(shortResult.intersects, "Stop gain beyond shorter peptide should not intersect.");
            Assert.IsTrue(shortResult.identifies, "Stop gain just beyond peptide end should identify (terminal change).");
        }
        [Test]
        public static void TestIsVariantPeptide()
        {
            Protein protein = new Protein("MPEPTIDENEWPEPTIDE", "protein0", appliedSequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "V", "substitution", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) });

            PeptideWithSetModifications pepe = new PeptideWithSetModifications(protein, new DigestionParams(), 1, 8, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications notPepe = new PeptideWithSetModifications(protein, new DigestionParams(), 9, 18, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);

            Assert.IsTrue(pepe.IsVariantPeptide());
            Assert.IsFalse(notPepe.IsVariantPeptide());
        }
        [Test]
        public static void TestSeqVarString()
        {
            // Protein baseline
            Protein protein = new Protein("MACDEFGHIK", "test");

            // 1. Substitution at N-terminus with variant-specific modification (M -> A + mod on A)
            var subMod = new Modification("mod on A", "mod", "mod", "mod");
            var vSubNterm = new SequenceVariation(
                oneBasedBeginPosition: 1,
                oneBasedEndPosition: 1,
                originalSequence: "M",
                variantSequence: "A",
                description: "nterm_substitution_with_variant_mod",
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { subMod } } });

            var pepFull = new PeptideWithSetModifications(
                protein, new DigestionParams(), 1, 10,
                CleavageSpecificity.Unknown, "", 0,
                new Dictionary<int, Modification>(), 0);

            Assert.AreEqual("M1A[mod:mod on A]", pepFull.SequenceVariantString(vSubNterm, true));

            // 2. Missense at peptide position 2 with variant-specific modification (A -> V)
            var pos2Mod = new Modification("mod on V", "mod", "mod", "mod");
            var vMissense = new SequenceVariation(
                oneBasedBeginPosition: 2,
                oneBasedEndPosition: 2,
                originalSequence: "A",
                variantSequence: "V",
                description: "missense_with_variant_mod",
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 2, new List<Modification> { pos2Mod } } });

            var pep2toEnd = new PeptideWithSetModifications(
                protein, new DigestionParams(), 2, 10,
                CleavageSpecificity.Unknown, "", 0,
                new Dictionary<int, Modification>(), 0);

            Assert.AreEqual("A2V[mod:mod on V]", pep2toEnd.SequenceVariantString(vMissense, true));

            // 3. Insertion / expansion: positions 7–9 (GHI -> GHIK)
            // Original segment (7–9) == GHI; variant adds K
            var vInsertion = new SequenceVariation(
                oneBasedBeginPosition: 7,
                oneBasedEndPosition: 9,
                originalSequence: "GHI",
                variantSequence: "GHIK",
                description: "insertion_extension");
            var pepMid = new PeptideWithSetModifications(
                protein, new DigestionParams(), 2, 10,
                CleavageSpecificity.Unknown, "", 0,
                new Dictionary<int, Modification>(), 0);
            Assert.AreEqual("GHI7GHIK", pepMid.SequenceVariantString(vInsertion, true));

            // 4. Frameshift/large replacement: full span (1–10) replaced by longer sequence
            var vFrameshift = new SequenceVariation(
                oneBasedBeginPosition: 1,
                oneBasedEndPosition: 10,
                originalSequence: "MACDEFGHIK",
                variantSequence: "MABCDEFGHIJKLMNOP",
                description: "frameshift_extension");
            Assert.AreEqual("MACDEFGHIK1MABCDEFGHIJKLMNOP", pepFull.SequenceVariantString(vFrameshift, true));

            // 5. Synonymous with variant-specific mod (no sequence change but mod should appear)
            var synMod = new Modification("mod on C", "mod", "mod", "mod");
            var vSynonymous = new SequenceVariation(
                oneBasedBeginPosition: 3,
                oneBasedEndPosition: 3,
                originalSequence: "C",
                variantSequence: "C",
                description: "synonymous_with_variant_mod",
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 3, new List<Modification> { synMod } } });
            Assert.AreEqual("C3C[mod:mod on C]", pepFull.SequenceVariantString(vSynonymous, true));
        }
        [Test]
        public static void BreakDeserializationMethod()
        {
            Assert.Throws<MzLibException>(() => new PeptideWithSetModifications("|", new Dictionary<string, Modification>())); // ambiguous
            Assert.Throws<MzLibException>(() => new PeptideWithSetModifications("[]", new Dictionary<string, Modification>())); // bad mod
            Assert.Throws<MzLibException>(() => new PeptideWithSetModifications("A[:mod]", new Dictionary<string, Modification>())); // nonexistent mod
        }

        [Test]
        public static void TestIdentifyandStringMethodsRevised()
        {
            // Picks a peptide that fully covers the variant span in the applied proteoform (if possible).
            // Notes:
            // - For insertions/deletions, the "effective variant end" includes the length delta.
            // - If the effective end would move before the begin (e.g., contraction past the begin), we clamp to begin.
            //   That specific clamp is exercising the branch:
            //     if (effectiveVariantEnd < appliedVariation.OneBasedBeginPosition)
            //         effectiveVariantEnd = appliedVariation.OneBasedBeginPosition;
            // - This helper is meant to deterministically hit the "peptide fully covers variant" path.
            static PeptideWithSetModifications PickCoveringPeptide(
                Protein variantProteoform,
                DigestionParams dp,
                SequenceVariation v)
            {
                var peps = variantProteoform
                    .Digest(dp, new List<Modification>(), new List<Modification>())
                    .OfType<PeptideWithSetModifications>()
                    .OrderBy(p => p.Length)
                    .ThenBy(p => p.OneBasedStartResidueInProtein)
                    .ToList();

                if (!peps.Any())
                    Assert.Fail($"No peptides produced for {variantProteoform.Accession}.");

                // Compute effective end of the variant after accounting for length difference
                // (e.g., insertions/expansions push the end right; deletions/contractions pull it left).
                int lengthDiff = v.VariantSequence.Length - v.OriginalSequence.Length;
                int effectiveVariantEnd = v.OneBasedEndPosition + lengthDiff;

                // Clamp if the effective end "overshot" left of the begin due to contraction.
                // This is the branch under test:
                //   if (effectiveVariantEnd < appliedVariation.OneBasedBeginPosition)
                //       effectiveVariantEnd = appliedVariation.OneBasedBeginPosition;
                if (effectiveVariantEnd < v.OneBasedBeginPosition)
                    effectiveVariantEnd = v.OneBasedBeginPosition;

                // Prefer the shortest peptide that fully contains [variantBegin..effectiveVariantEnd]
                var covering = peps
                    .Where(p => p.OneBasedStartResidueInProtein <= v.OneBasedBeginPosition
                             && p.OneBasedEndResidueInProtein >= effectiveVariantEnd)
                    .OrderBy(p => p.Length)
                    .ThenBy(p => p.OneBasedStartResidueInProtein)
                    .FirstOrDefault();

                // Fallback (no full-cover peptide): return the shortest overall
                return covering ?? peps.First();
            }

            // Picks a peptide around the variant begin anchor (or a requested index) to exercise
            // intersect/non-intersect and terminal-cleavage identification logic deterministically.
            static PeptideWithSetModifications PickPeptide(
                Protein variantProteoform,
                DigestionParams dp,
                SequenceVariation v,
                int? requestedIndex)
            {
                var peps = variantProteoform
                    .Digest(dp, new List<Modification>(), new List<Modification>())
                    .OfType<PeptideWithSetModifications>()
                    .OrderBy(p => p.OneBasedStartResidueInProtein)
                    .ThenBy(p => p.Length)
                    .ToList();

                if (!peps.Any())
                    Assert.Fail($"No peptides produced for {variantProteoform.Accession}.");

                if (requestedIndex.HasValue && requestedIndex.Value < peps.Count)
                    return peps[requestedIndex.Value];

                // Anchor near the variant begin in the applied proteoform coordinate space
                int anchor = Math.Min(v.OneBasedBeginPosition, variantProteoform.BaseSequence.Length);

                // Choose a peptide that spans the anchor residue if possible
                var covering = peps.FirstOrDefault(p =>
                    p.OneBasedStartResidueInProtein <= anchor &&
                    p.OneBasedEndResidueInProtein >= Math.Min(anchor, variantProteoform.BaseSequence.Length));

                return covering ?? peps.First();
            }

            // Build two simple mods used by some variant cases (variant-specific PTMs)
            ModificationMotif.TryGetMotif("V", out var motifV);
            ModificationMotif.TryGetMotif("P", out var motifP);
            var mv = new Modification("mod", null, "type", null, motifV, "Anywhere.", null, 42.01,
                new Dictionary<string, IList<string>>(), null, null, null, null, null);
            var mp = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01,
                new Dictionary<string, IList<string>>(), null, null, null, null, null);

            // Protein-level PTM on P(4) for testing combined variant/PTM handling
            var proteinPMods = new Dictionary<int, List<Modification>> { { 4, new List<Modification> { mp } } };

            // Each protein has a single variant. Some are substitutions (equal-length),
            // some are insertions (expansion), deletions (contraction), or stops.
            // These cover the major branches inside IntersectsAndIdentifiesVariation:
            //  - Intersection determination (original and effective windows)
            //  - Equal-length substitution identification (per-residue differences)
            //  - Insertion/deletion identification rules
            //  - Terminal changes (stop gains/losses) affecting cleavage identification when non-intersecting
            var proteins = new List<(string Label, Protein Protein)>
            {
                ("protein0",  new Protein("MPEPTIDE","protein0",
                    sequenceVariations: new(){ new SequenceVariation(4,4,"P","V","substitution") })),
                // protein1 is a multi-residue equal-length substitution (PT->KT, span 4..5).
                // For reliable identification, we construct a peptide that spans exactly 4..5.
                ("protein1",  new Protein("MPEPTIDE","protein1",
                    sequenceVariations: new(){ new SequenceVariation(4,5,"PT","KT","mnp") })),
                ("protein2",  new Protein("MPEPTIDE","protein2",
                    sequenceVariations: new(){ new SequenceVariation(4,4,"P","PPP","insertion") })),
                ("protein3",  new Protein("MPEPPPTIDE","protein3",
                    sequenceVariations: new(){ new SequenceVariation(4,6,"PPP","P","deletion") })),
                ("protein4",  new Protein("MPEPKPKTIDE","protein4",
                    sequenceVariations: new(){ new SequenceVariation(4,7,"PKPK","PK","internal_deletion") })),
                ("protein5",  new Protein("MPEPTAIDE","protein5",
                    sequenceVariations: new(){ new SequenceVariation(4,6,"PTA","KT","mnp") })),
                ("protein6",  new Protein("MPEKKAIDE","protein6",
                    sequenceVariations: new(){ new SequenceVariation(4,6,"KKA","K","deletion") })),
                // Variant-specific mod added at pos 4 (post-variation coordinates)
                ("protein7",  new Protein("MPEPTIDE","protein7",
                    sequenceVariations: new(){ new SequenceVariation(4,4,"P","V","substitution_with_variant_mod",
                        oneBasedModifications: new Dictionary<int,List<Modification>>{{4,new(){mv}}}) })),
                // Insertion with a variant mod located within the inserted region (post-variation pos 5)
                ("protein8",  new Protein("MPEPTIDE","protein8",
                    sequenceVariations: new(){ new SequenceVariation(4,4,"P","PPP","insertion_with_variant_mod",
                        oneBasedModifications: new Dictionary<int,List<Modification>>{{5,new(){mp}}}) })),
                ("protein9",  new Protein("MPEPTIDEPEPTIDE","protein9",
                    sequenceVariations: new(){ new SequenceVariation(4,15,"PTIDEPEPTIDE","PPP","replacement_contraction") })),
                // Protein-level PTM co-exists with a substitution at position 4
                ("protein10", new Protein("MPEPTIDE","protein10",
                    oneBasedModifications: proteinPMods,
                    sequenceVariations: new(){ new SequenceVariation(4,4,"P","V","substitution_with_protein_mod") })),
                // Stop gain inside peptide span (intersect=false can still identify via terminal logic for flanks; here intersecting case)
                ("protein11", new Protein("MPEPTIDE","protein11",
                    sequenceVariations: new(){ new SequenceVariation(5,5,"T","*","stop_gain_identifying") })),
                // Same stop but in a context that should not identify for chosen peptide
                ("protein12", new Protein("MPEKTIDE","protein12",
                    sequenceVariations: new(){ new SequenceVariation(5,5,"T","*","stop_gain_non_identifying") })),
                ("protein13", new Protein("MPEPTIPEPEPTIPE","protein13",
                    sequenceVariations: new(){ new SequenceVariation(7,7,"P","D","missense") })),
                // Extension at position 8 (E->EK) tests insertion-like behavior at a single position
                ("protein14", new Protein("MPEPTIDE","protein14",
                    sequenceVariations: new(){ new SequenceVariation(8,8,"E","EK","extension") })),
                // Stop loss extension beyond peptide end; used to assert non-identifying for certain flanks
                ("protein15", new Protein("MPEPTIDE","protein15",
                    sequenceVariations: new(){ new SequenceVariation(9,9,"*","KMPEP","stop_loss_extension") }))
            };

            // Expected variant-string encodings for a subset of the above (label -> expected).
            // These strings summarize: OriginalSubstr + BeginIndex + VariantSubstr (+ [mod annotations] if present).
            var expectedVariantStrings = new Dictionary<string, string>
            {
                {"protein0","P4V"},
                {"protein1","PT4KT"},
                {"protein2","P4PPP"}, // insertion keeps full variant (no compression)
                {"protein3","PPP4P"},
                {"protein5","PTA4KT"},
                {"protein6","KKA4K"},
                {"protein7","P4V[type:mod on V]"},
                {"protein8","P4PP[type:mod on P]P"},
                {"protein9","PTIDEPEPTIDE4PPP"},
                {"protein10","P4V"},
                {"protein11","T5*"},
                {"protein13","P7D"}
            };

            var dpTrypsin = new DigestionParams(minPeptideLength: 2);
            var dpAspN = new DigestionParams(protease: "Asp-N", minPeptideLength: 2);
            var dpLysN = new DigestionParams(protease: "Lys-N", minPeptideLength: 2);

            // Build a map of label -> (applied proteoform, applied variant)
            // If a proteoform with AppliedSequenceVariations exists, we prefer it for testing the "applied space".
            int autoApplied = 0;
            var appliedMap = new Dictionary<string, (Protein Prot, SequenceVariation V)>();

            foreach (var (label, prot) in proteins)
            {
                var variant = prot.SequenceVariations.Single();
                var applied = prot
                    .GetVariantBioPolymers(maxSequenceVariantIsoforms: 50)
                    .OfType<Protein>()
                    .FirstOrDefault(p => p.AppliedSequenceVariations.Any());

                if (applied != null)
                {
                    autoApplied++;
                    appliedMap[label] = (applied, applied.AppliedSequenceVariations.First());
                }
                else
                {
                    appliedMap[label] = (prot, variant);
                }
            }

            TestContext.WriteLine($"[INFO] Variant application summary: autoApplied={autoApplied}, total={appliedMap.Count}");

            // protein0: simple point substitution P->V at pos 4, covered under both proteases
            (Protein p0v, var v0) = appliedMap["protein0"];
            var p0_pep = PickCoveringPeptide(p0v, dpTrypsin, v0);
            Assert.AreEqual((true, true), p0_pep.IntersectsAndIdentifiesVariation(v0));
            var p0_pep2 = PickCoveringPeptide(p0v, dpAspN, v0);
            Assert.AreEqual((true, true), p0_pep2.IntersectsAndIdentifiesVariation(v0));

            // protein1: multi-residue equal-length substitution PT(4..5)->KT.
            // To avoid ambiguity from digestion (e.g., tryptic cleavage near K), construct a peptide from the
            // non-applied proteoform that spans exactly 4..5 and test against the "raw" variant. This ensures a
            // full-window overlap and deterministic identification via per-residue difference (P!=K).
            (Protein p1v, var v1) = appliedMap["protein1"];
            var v1Raw = proteins.First(p => p.Label == "protein1").Protein.SequenceVariations.Single();
            var p1_origin = proteins.First(p => p.Label == "protein1").Protein; // non-applied proteoform

            var p1_pep = new PeptideWithSetModifications(
                p1_origin,
                dpTrypsin,
                oneBasedStartResidueInProtein: 4,
                oneBasedEndResidueInProtein: 5, // exactly the variant window
                CleavageSpecificity.Full,
                peptideDescription: "",
                missedCleavages: 0,
                allModsOneIsNterminus: new Dictionary<int, Modification>(),
                numFixedMods: 0);

            // Expected: (intersects=true, identifies=true) because equal-length substitution differs inside overlap.
            // Also exercises downstream string building for multi-residue substitutions.
            Assert.AreEqual((true, true), p1_pep.IntersectsAndIdentifiesVariation(v1Raw));

            // protein7: substitution with a variant-specific PTM (annotation should still identify)
            (Protein p7v, var v7) = appliedMap["protein7"];
            var p7_pep = PickCoveringPeptide(p7v, dpTrypsin, v7);
            Assert.AreEqual((true, true), p7_pep.IntersectsAndIdentifiesVariation(v7));

            // protein10: substitution co-existing with a protein-level PTM at the same site
            (Protein p10v, var v10) = appliedMap["protein10"];
            var p10_pep = PickCoveringPeptide(p10v, dpTrypsin, v10);
            Assert.AreEqual((true, true), p10_pep.IntersectsAndIdentifiesVariation(v10));

            // protein2/protein3/protein4/protein5: insertion and deletion flavors
            // Insertions/expansions or deletions/contractions that overlap are identifying.
            (Protein p2v, var v2) = appliedMap["protein2"];
            var p2_pep = PickCoveringPeptide(p2v, dpTrypsin, v2);
            Assert.AreEqual((true, true), p2_pep.IntersectsAndIdentifiesVariation(v2));

            (Protein p3v, var v3) = appliedMap["protein3"];
            var p3_pep = PickCoveringPeptide(p3v, dpTrypsin, v3);
            Assert.AreEqual((true, true), p3_pep.IntersectsAndIdentifiesVariation(v3));

            (Protein p4v, var v4) = appliedMap["protein4"];
            var p4_pep = PickCoveringPeptide(p4v, dpTrypsin, v4);
            Assert.AreEqual((true, true), p4_pep.IntersectsAndIdentifiesVariation(v4));

            (Protein p5v, var v5) = appliedMap["protein5"];
            var p5_pep = PickCoveringPeptide(p5v, dpTrypsin, v5);
            Assert.AreEqual((true, true), p5_pep.IntersectsAndIdentifiesVariation(v5));

            // protein6: deletion; even partial overlapping deletions are considered identifying once intersecting.
            (Protein p6v, var v6) = appliedMap["protein6"];
            var p6_pep = PickPeptide(p6v, dpTrypsin, v6, 2);
            Assert.AreEqual((true, true), p6_pep.IntersectsAndIdentifiesVariation(v6));

            // protein8/protein9: insertion-with-mod and replacement-contraction cases
            (Protein p8v, var v8) = appliedMap["protein8"];
            var p8_pep = PickCoveringPeptide(p8v, dpTrypsin, v8);
            Assert.AreEqual((true, true), p8_pep.IntersectsAndIdentifiesVariation(v8));

            (Protein p9v, var v9) = appliedMap["protein9"];
            var p9_pep = PickCoveringPeptide(p9v, dpTrypsin, v9);
            Assert.AreEqual((true, true), p9_pep.IntersectsAndIdentifiesVariation(v9));

            // protein11: stop gain that can be identified even when the chosen peptide doesn’t overlap,
            // via terminal-cleavage logic (new terminal introduced). We assert (false, true) using two proteases.
            (Protein p11v, var v11) = appliedMap["protein11"];
            var p11_pep_AspN = PickPeptide(p11v, dpAspN, v11, 0);
            Assert.AreEqual((false, true), p11_pep_AspN.IntersectsAndIdentifiesVariation(v11));
            var p11_pep_Tryp = PickPeptide(p11v, dpTrypsin, v11, 0);
            Assert.AreEqual((false, true), p11_pep_Tryp.IntersectsAndIdentifiesVariation(v11));

            // protein12: stop gain in a context that should not identify for the peptide chosen
            (Protein p12v, var v12) = appliedMap["protein12"];
            var p12_pep = PickPeptide(p12v, dpTrypsin, v12, 0);
            Assert.AreEqual((false, false), p12_pep.IntersectsAndIdentifiesVariation(v12));

            // protein13: missense away from anchor, demonstrate non-intersecting but identifying due to rules
            (Protein p13v, var v13) = appliedMap["protein13"];
            var p13_pep = PickPeptide(p13v, dpAspN, v13, 0);
            Assert.AreEqual((false, true), p13_pep.IntersectsAndIdentifiesVariation(v13));

            // protein14: single-position extension (E->EK) treated like insertion at that coordinate
            (Protein p14v, var v14) = appliedMap["protein14"];
            var p14_pep = PickPeptide(p14v, dpLysN, v14, 0);
            Assert.AreEqual((true, true), p14_pep.IntersectsAndIdentifiesVariation(v14));
            AssertVariantStringIfExpected("protein14", p14_pep, v14, true);

            // protein15: stop loss extension beyond peptide end in a context that should not identify
            (Protein p15v, var v15) = appliedMap["protein15"];
            var p15_pep = PickPeptide(p15v, dpLysN, v15, 0);
            Assert.AreEqual((false, false), p15_pep.IntersectsAndIdentifiesVariation(v15));

            // Helper for asserting variant-string outputs only when expected is provided.
            // The boolean intersectsFlag is the legacy "intersects" parameter used by SequenceVariantString overloads.
            void AssertVariantStringIfExpected(string label, PeptideWithSetModifications pep, SequenceVariation v, bool intersectsFlag)
            {
                if (!expectedVariantStrings.TryGetValue(label, out var expected))
                    return;
                var actual = pep.SequenceVariantString(v, intersectsFlag);
                Assert.AreEqual(expected, actual, $"Variant string mismatch for {label} (intersectsFlag={intersectsFlag})");
            }

            // Validate the human-readable variant strings for selected cases
            AssertVariantStringIfExpected("protein0", p0_pep, v0, true);
            AssertVariantStringIfExpected("protein0", p0_pep2, v0, true);
            AssertVariantStringIfExpected("protein1", p1_pep, v1, true);
            AssertVariantStringIfExpected("protein2", p2_pep, v2, true);
            AssertVariantStringIfExpected("protein3", p3_pep, v3, true);
            AssertVariantStringIfExpected("protein5", p5_pep, v5, true);
            AssertVariantStringIfExpected("protein6", p6_pep, v6, true);
            AssertVariantStringIfExpected("protein7", p7_pep, v7, true);
            AssertVariantStringIfExpected("protein8", p8_pep, v8, true);
            AssertVariantStringIfExpected("protein9", p9_pep, v9, true);
            AssertVariantStringIfExpected("protein10", p10_pep, v10, true);
            AssertVariantStringIfExpected("protein11", p11_pep_AspN, v11, false);
            AssertVariantStringIfExpected("protein11", p11_pep_Tryp, v11, false);
            AssertVariantStringIfExpected("protein13", p13_pep, v13, false);

            TestContext.WriteLine("[INFO] TestIdentifyandStringMethods completed (deletion overlaps now intersect & identify).");
        }


        [Test]
        public static void TestIdentifyandStringMethods()
        {
            static PeptideWithSetModifications PickCoveringPeptide(
                Protein variantProteoform,
                DigestionParams dp,
                SequenceVariation v)
            {
                var peps = variantProteoform
                    .Digest(dp, new List<Modification>(), new List<Modification>())
                    .OfType<PeptideWithSetModifications>()
                    .OrderBy(p => p.Length)
                    .ThenBy(p => p.OneBasedStartResidueInProtein)
                    .ToList();

                if (!peps.Any())
                    Assert.Fail($"No peptides produced for {variantProteoform.Accession}.");

                int lengthDiff = v.VariantSequence.Length - v.OriginalSequence.Length;
                int effectiveVariantEnd = v.OneBasedEndPosition + lengthDiff;
                if (effectiveVariantEnd < v.OneBasedBeginPosition)
                    effectiveVariantEnd = v.OneBasedBeginPosition;

                var covering = peps
                    .Where(p => p.OneBasedStartResidueInProtein <= v.OneBasedBeginPosition
                             && p.OneBasedEndResidueInProtein >= effectiveVariantEnd)
                    .OrderBy(p => p.Length)
                    .ThenBy(p => p.OneBasedStartResidueInProtein)
                    .FirstOrDefault();

                return covering ?? peps.First();
            }

            static PeptideWithSetModifications PickPeptide(
                Protein variantProteoform,
                DigestionParams dp,
                SequenceVariation v,
                int? requestedIndex)
            {
                var peps = variantProteoform
                    .Digest(dp, new List<Modification>(), new List<Modification>())
                    .OfType<PeptideWithSetModifications>()
                    .OrderBy(p => p.OneBasedStartResidueInProtein)
                    .ThenBy(p => p.Length)
                    .ToList();

                if (!peps.Any())
                    Assert.Fail($"No peptides produced for {variantProteoform.Accession}.");

                if (requestedIndex.HasValue && requestedIndex.Value < peps.Count)
                    return peps[requestedIndex.Value];

                int anchor = Math.Min(v.OneBasedBeginPosition, variantProteoform.BaseSequence.Length);

                var covering = peps.FirstOrDefault(p =>
                    p.OneBasedStartResidueInProtein <= anchor &&
                    p.OneBasedEndResidueInProtein >= Math.Min(anchor, variantProteoform.BaseSequence.Length));

                return covering ?? peps.First();
            }

            ModificationMotif.TryGetMotif("V", out var motifV);
            ModificationMotif.TryGetMotif("P", out var motifP);
            var mv = new Modification("mod", null, "type", null, motifV, "Anywhere.", null, 42.01,
                new Dictionary<string, IList<string>>(), null, null, null, null, null);
            var mp = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01,
                new Dictionary<string, IList<string>>(), null, null, null, null, null);

            var proteinPMods = new Dictionary<int, List<Modification>> { { 4, new List<Modification> { mp } } };

            var proteins = new List<(string Label, Protein Protein)>
            {
                ("protein0",  new Protein("MPEPTIDE","protein0",
                    sequenceVariations: new(){ new SequenceVariation(4,4,"P","V","substitution") })),
                ("protein1",  new Protein("MPEPTIDE","protein1",
                    sequenceVariations: new(){ new SequenceVariation(4,5,"PT","KT","mnp") })),
                ("protein2",  new Protein("MPEPTIDE","protein2",
                    sequenceVariations: new(){ new SequenceVariation(4,4,"P","PPP","insertion") })),
                ("protein3",  new Protein("MPEPPPTIDE","protein3",
                    sequenceVariations: new(){ new SequenceVariation(4,6,"PPP","P","deletion") })),
                ("protein4",  new Protein("MPEPKPKTIDE","protein4",
                    sequenceVariations: new(){ new SequenceVariation(4,7,"PKPK","PK","internal_deletion") })),
                ("protein5",  new Protein("MPEPTAIDE","protein5",
                    sequenceVariations: new(){ new SequenceVariation(4,6,"PTA","KT","mnp") })),
                ("protein6",  new Protein("MPEKKAIDE","protein6",
                    sequenceVariations: new(){ new SequenceVariation(4,6,"KKA","K","deletion") })),
                ("protein7",  new Protein("MPEPTIDE","protein7",
                    sequenceVariations: new(){ new SequenceVariation(4,4,"P","V","substitution_with_variant_mod",
                        oneBasedModifications: new Dictionary<int,List<Modification>>{{4,new(){mv}}}) })),
                ("protein8",  new Protein("MPEPTIDE","protein8",
                    sequenceVariations: new(){ new SequenceVariation(4,4,"P","PPP","insertion_with_variant_mod",
                        oneBasedModifications: new Dictionary<int,List<Modification>>{{5,new(){mp}}}) })),
                ("protein9",  new Protein("MPEPTIDEPEPTIDE","protein9",
                    sequenceVariations: new(){ new SequenceVariation(4,15,"PTIDEPEPTIDE","PPP","replacement_contraction") })),
                ("protein10", new Protein("MPEPTIDE","protein10",
                    oneBasedModifications: proteinPMods,
                    sequenceVariations: new(){ new SequenceVariation(4,4,"P","V","substitution_with_protein_mod") })),
                ("protein11", new Protein("MPEPTIDE","protein11",
                    sequenceVariations: new(){ new SequenceVariation(5,5,"T","*","stop_gain_identifying") })),
                ("protein12", new Protein("MPEKTIDE","protein12",
                    sequenceVariations: new(){ new SequenceVariation(5,5,"T","*","stop_gain_non_identifying") })),
                ("protein13", new Protein("MPEPTIPEPEPTIPE","protein13",
                    sequenceVariations: new(){ new SequenceVariation(7,7,"P","D","missense") })),
                ("protein14", new Protein("MPEPTIDE","protein14",
                    sequenceVariations: new(){ new SequenceVariation(8,8,"E","EK","extension") })),
                ("protein15", new Protein("MPEPTIDE","protein15",
                    sequenceVariations: new(){ new SequenceVariation(9,9,"*","KMPEP","stop_loss_extension") }))
            };

            var expectedVariantStrings = new Dictionary<string, string>
            {
                {"protein0","P4V"},
                {"protein1","PT4KT"},
                {"protein2","P4PPP"}, // restored full insertion (no compression)
                {"protein3","PPP4P"},
                {"protein5","PTA4KT"},
                {"protein6","KKA4K"},
                {"protein7","P4V[type:mod on V]"},
                {"protein8","P4PP[type:mod on P]P"},
                {"protein9","PTIDEPEPTIDE4PPP"},
                {"protein10","P4V"},
                {"protein11","T5*"},
                {"protein13","P7D"}
            };

            var dpTrypsin = new DigestionParams(minPeptideLength: 2);
            var dpAspN = new DigestionParams(protease: "Asp-N", minPeptideLength: 2);
            var dpLysN = new DigestionParams(protease: "Lys-N", minPeptideLength: 2);

            int autoApplied = 0;
            var appliedMap = new Dictionary<string, (Protein Prot, SequenceVariation V)>();

            foreach (var (label, prot) in proteins)
            {
                var variant = prot.SequenceVariations.Single();
                var applied = prot
                    .GetVariantBioPolymers(maxSequenceVariantIsoforms: 50)
                    .OfType<Protein>()
                    .FirstOrDefault(p => p.AppliedSequenceVariations.Any());

                if (applied != null)
                {
                    autoApplied++;
                    appliedMap[label] = (applied, applied.AppliedSequenceVariations.First());
                }
                else
                {
                    appliedMap[label] = (prot, variant);
                }
            }

            TestContext.WriteLine($"[INFO] Variant application summary: autoApplied={autoApplied}, total={appliedMap.Count}");

            (Protein p0v, var v0) = appliedMap["protein0"];
            var p0_pep = PickCoveringPeptide(p0v, dpTrypsin, v0);
            Assert.AreEqual((true, true), p0_pep.IntersectsAndIdentifiesVariation(v0));
            var p0_pep2 = PickCoveringPeptide(p0v, dpAspN, v0);
            Assert.AreEqual((true, true), p0_pep2.IntersectsAndIdentifiesVariation(v0));

            (Protein p1v, var v1) = appliedMap["protein1"];
            // Use the raw (non-applied) variant and construct a peptide that exactly spans the variant window (4..5).
            var v1Raw = proteins.First(p => p.Label == "protein1").Protein.SequenceVariations.Single();
            var p1_origin = proteins.First(p => p.Label == "protein1").Protein; // non-applied proteoform
            var p1_pep = new PeptideWithSetModifications(
                p1_origin,
                dpTrypsin,
                oneBasedStartResidueInProtein: 4,
                oneBasedEndResidueInProtein: 5,
            CleavageSpecificity.Full,
            peptideDescription: "",
            missedCleavages: 0,
            allModsOneIsNterminus: new Dictionary<int, Modification>(),
            numFixedMods: 0);
            Assert.AreEqual((true, true), p1_pep.IntersectsAndIdentifiesVariation(v1Raw));

            (Protein p7v, var v7) = appliedMap["protein7"];
            var p7_pep = PickCoveringPeptide(p7v, dpTrypsin, v7);
            Assert.AreEqual((true, true), p7_pep.IntersectsAndIdentifiesVariation(v7));

            (Protein p10v, var v10) = appliedMap["protein10"];
            var p10_pep = PickCoveringPeptide(p10v, dpTrypsin, v10);
            Assert.AreEqual((true, true), p10_pep.IntersectsAndIdentifiesVariation(v10));

            (Protein p2v, var v2) = appliedMap["protein2"];
            var p2_pep = PickCoveringPeptide(p2v, dpTrypsin, v2);
            Assert.AreEqual((true, true), p2_pep.IntersectsAndIdentifiesVariation(v2));

            (Protein p3v, var v3) = appliedMap["protein3"];
            var p3_pep = PickCoveringPeptide(p3v, dpTrypsin, v3);
            Assert.AreEqual((true, true), p3_pep.IntersectsAndIdentifiesVariation(v3));

            (Protein p4v, var v4) = appliedMap["protein4"];
            var p4_pep = PickCoveringPeptide(p4v, dpTrypsin, v4);
            Assert.AreEqual((true, true), p4_pep.IntersectsAndIdentifiesVariation(v4));

            (Protein p5v, var v5) = appliedMap["protein5"];
            var p5_pep = PickCoveringPeptide(p5v, dpTrypsin, v5);
            Assert.AreEqual((true, true), p5_pep.IntersectsAndIdentifiesVariation(v5));

            (Protein p6v, var v6) = appliedMap["protein6"];
            // Updated expectation: deletion overlap ⇒ (true,true)
            var p6_pep = PickPeptide(p6v, dpTrypsin, v6, 2);
            Assert.AreEqual((true, true), p6_pep.IntersectsAndIdentifiesVariation(v6));

            (Protein p8v, var v8) = appliedMap["protein8"];
            var p8_pep = PickCoveringPeptide(p8v, dpTrypsin, v8);
            Assert.AreEqual((true, true), p8_pep.IntersectsAndIdentifiesVariation(v8));

            (Protein p9v, var v9) = appliedMap["protein9"];
            var p9_pep = PickCoveringPeptide(p9v, dpTrypsin, v9);
            Assert.AreEqual((true, true), p9_pep.IntersectsAndIdentifiesVariation(v9));

            (Protein p11v, var v11) = appliedMap["protein11"];
            var p11_pep_AspN = PickPeptide(p11v, dpAspN, v11, 0);
            Assert.AreEqual((false, true), p11_pep_AspN.IntersectsAndIdentifiesVariation(v11));
            var p11_pep_Tryp = PickPeptide(p11v, dpTrypsin, v11, 0);
            Assert.AreEqual((false, true), p11_pep_Tryp.IntersectsAndIdentifiesVariation(v11));

            (Protein p12v, var v12) = appliedMap["protein12"];
            var p12_pep = PickPeptide(p12v, dpTrypsin, v12, 0);
            Assert.AreEqual((false, false), p12_pep.IntersectsAndIdentifiesVariation(v12));

            (Protein p13v, var v13) = appliedMap["protein13"];
            var p13_pep = PickPeptide(p13v, dpAspN, v13, 0);
            Assert.AreEqual((false, true), p13_pep.IntersectsAndIdentifiesVariation(v13));

            (Protein p14v, var v14) = appliedMap["protein14"];
            var p14_pep = PickPeptide(p14v, dpLysN, v14, 0);
            Assert.AreEqual((true, true), p14_pep.IntersectsAndIdentifiesVariation(v14));
            AssertVariantStringIfExpected("protein14", p14_pep, v14, true); // if you decide to include it in expected strings

            (Protein p15v, var v15) = appliedMap["protein15"];
            var p15_pep = PickPeptide(p15v, dpLysN, v15, 0);
            Assert.AreEqual((false, false), p15_pep.IntersectsAndIdentifiesVariation(v15));

            void AssertVariantStringIfExpected(string label, PeptideWithSetModifications pep, SequenceVariation v, bool intersectsFlag)
            {
                if (!expectedVariantStrings.TryGetValue(label, out var expected))
                    return;
                var actual = pep.SequenceVariantString(v, intersectsFlag);
                Assert.AreEqual(expected, actual, $"Variant string mismatch for {label} (intersectsFlag={intersectsFlag})");
            }

            AssertVariantStringIfExpected("protein0", p0_pep, v0, true);
            AssertVariantStringIfExpected("protein0", p0_pep2, v0, true);
            AssertVariantStringIfExpected("protein1", p1_pep, v1, true);
            AssertVariantStringIfExpected("protein2", p2_pep, v2, true);
            AssertVariantStringIfExpected("protein3", p3_pep, v3, true);
            AssertVariantStringIfExpected("protein5", p5_pep, v5, true);
            AssertVariantStringIfExpected("protein6", p6_pep, v6, true); // intersects now true
            AssertVariantStringIfExpected("protein7", p7_pep, v7, true);
            AssertVariantStringIfExpected("protein8", p8_pep, v8, true);
            AssertVariantStringIfExpected("protein9", p9_pep, v9, true);
            AssertVariantStringIfExpected("protein10", p10_pep, v10, true);
            AssertVariantStringIfExpected("protein11", p11_pep_AspN, v11, false);
            AssertVariantStringIfExpected("protein11", p11_pep_Tryp, v11, false);
            AssertVariantStringIfExpected("protein13", p13_pep, v13, false);

            TestContext.WriteLine("[INFO] TestIdentifyandStringMethods completed (deletion overlaps now intersect & identify).");
        }
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

            //  chymotrypsin (don't cleave before proline)
            newAminoAcidPositions = new int["FKFPRWAWPSYGYPG".Length];
            PeptideWithSetModifications p_chymoP = new PeptideWithSetModifications(new Protein("FKFPRWAWPSYGYPG", "DECOY_CHYMOP"), new DigestionParams(protease: "chymotrypsin (don't cleave before proline)", maxMissedCleavages: 10), 1, 15, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications p_chymoP_reverse = p_chymoP.GetReverseDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual("FGPYGWSPWAYRPFK", p_chymoP_reverse.BaseSequence);
            Assert.AreEqual(p_chymoP.FullSequence, p_chymoP_reverse.PeptideDescription);

            //  chymotrypsin (cleave before proline)
            newAminoAcidPositions = new int["FKFPRWAWPSYGYPG".Length];
            PeptideWithSetModifications p_chymo = new PeptideWithSetModifications(new Protein("FKFPRWAWPSYGYPG", "DECOY_CHYMO"), new DigestionParams(protease: "chymotrypsin (cleave before proline)", maxMissedCleavages: 10), 1, 15, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications p_chymo_reverse = p_chymo.GetReverseDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual("FGFPGWSWPAYRYPK", p_chymo_reverse.BaseSequence);
            Assert.AreEqual(p_chymo.FullSequence, p_chymo_reverse.PeptideDescription);

            //  CNBr cleave after M
            newAminoAcidPositions = new int["MPEPTIMEK".Length];
            PeptideWithSetModifications p_cnbr = new PeptideWithSetModifications(new Protein("MPEPTIMEK", "DECOY_CNBR"), new DigestionParams(protease: "CNBr"), 1, 9, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications p_cnbr_reverse = p_cnbr.GetReverseDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual("MKEITPMEP", p_cnbr_reverse.BaseSequence);
            Assert.AreEqual(p_cnbr.FullSequence, p_cnbr_reverse.PeptideDescription);

            //  elastase cleave after A, V, S, G, L, I,
            newAminoAcidPositions = new int["KAYVPSRGHLDIN".Length];
            PeptideWithSetModifications p_elastase = new PeptideWithSetModifications(new Protein("KAYVPSRGHLDIN", "DECOY_ELASTASE"), new DigestionParams(protease: "elastase"), 1, 13, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            PeptideWithSetModifications p_elastase_reverse = p_elastase.GetReverseDecoyFromTarget(newAminoAcidPositions);
            Assert.AreEqual("NADVHSRGPLYIK", p_elastase_reverse.BaseSequence);

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
        // Helper: make a minimal peptide from a protein interval
        private static PeptideWithSetModifications MakePep(Protein prot, int begin, int end)
        {
            var dp = new DigestionParams(); // default protease/settings are fine; not used in these branches
            return new PeptideWithSetModifications(
                protein: prot,
                digestionParams: dp,
                oneBasedStartResidueInProtein: begin,
                oneBasedEndResidueInProtein: end,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "unit-test",
                missedCleavages: 0,
                allModsOneIsNterminus: new Dictionary<int, Modification>(),
                numFixedMods: 0);
        }

        [Test]
        public static void IntersectsAndIdentifiesVariation_EffectiveVariantEndClamped_And_EffectiveDegenerate_EarlyReturn()
        {
            // Protein indices:        1 2 3 4 5 6 7 8 9 10 11 12 ...
            // Sequence (20 AAs):     A C D E F G H I K  L  M  N  P  Q  R  S  T  V  W  Y
            // Variant: deletion of 5..10 ("FGHIKL") → VariantSequence == "" (lengthDiff negative)
            // Peptide under test: 8..12 (overlaps the original region but, after effective clamp, becomes degenerate)
            // Why this triggers the clamp:
            // - effectiveVariantEnd = end + (len(variant) - len(original)) = 10 + (0 - 6) = 4 < begin(=5) → clamped to 5
            // - intersectStartEff = max(pepStart=8, varBegin=5) = 8; intersectEndEff = min(pepEnd=12, effEnd=5) = 5
            //   -> intersectEndEff (5) < intersectStartEff (8) → effectiveDegenerate == true → early return
            var prot = new Protein("ACDEFGHIKLMNPQRSTVWY", "P1");
            var pep = MakePep(prot, begin: 8, end: 12);

            var deletion = new SequenceVariation(
                oneBasedBeginPosition: 5,
                oneBasedEndPosition: 10,
                originalSequence: "FGHIKL",
                variantSequence: string.Empty, // deletion
                description: "del 5..10");

            var (intersects, identifies) = pep.IntersectsAndIdentifiesVariation(deletion);

            // For deletions, code sets identifiesFlag = true; and early return occurs due to effectiveDegenerate
            Assert.Multiple(() =>
            {
                Assert.That(intersects, Is.True, "Expected 'intersects' == true (original region overlaps the peptide).");
                Assert.That(identifies, Is.True, "Deletion should set identifiesFlag = true.");
            });

            TestContext.WriteLine("Early-return path hit: effectiveVariantEnd clamped below begin and effectiveDegenerate == true");
        }

        [Test]
        public static void IntersectsAndIdentifiesVariation_NoClamp_NonDegenerate_ContinuesAndIdentifies()
        {
            // Same protein as above. Use a same-length substitution 5..7 where sequences differ.
            // Variant: 5..7 original "FGH" replaced with "YYY" (lengthDiff = 0) → no clamp.
            // Peptide under test: 5..7 (fully covers the variant span).
            // Since we cross the entire effective variant and Original != Variant over that window, identifiesFlag becomes true.
            var prot = new Protein("ACDEFGHIKLMNPQRSTVWY", "P2");
            var pep = MakePep(prot, begin: 5, end: 7);

            var substitution = new SequenceVariation(
                oneBasedBeginPosition: 5,
                oneBasedEndPosition: 7,
                originalSequence: "FGH",
                variantSequence: "YYY",
                description: "sub 5..7 FGH->YYY");

            var (intersects, identifies) = pep.IntersectsAndIdentifiesVariation(substitution);

            // No clamp (lengthDiff == 0), effectiveDegenerate == false (non-empty overlap), and sequences differ across full window -> identifies == true
            Assert.Multiple(() =>
            {
                Assert.That(intersects, Is.True, "Expected 'intersects' == true.");
                Assert.That(identifies, Is.True, "Expected identify due to full-span substitution with differing sequence.");
            });

            TestContext.WriteLine("Non-degenerate path hit: no clamp, full-span substitution identified correctly");
        }
    }
}