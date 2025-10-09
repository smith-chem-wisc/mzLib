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
            Assert.Throws<MzLibUtil.MzLibException>(() => new PeptideWithSetModifications("|", new Dictionary<string, Modification>())); // ambiguous
            Assert.Throws<MzLibUtil.MzLibException>(() => new PeptideWithSetModifications("[]", new Dictionary<string, Modification>())); // bad mod
            Assert.Throws<MzLibUtil.MzLibException>(() => new PeptideWithSetModifications("A[:mod]", new Dictionary<string, Modification>())); // nonexistent mod
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

            (Protein p7v, var v7) = appliedMap["protein7"];
            var p7_pep = PickCoveringPeptide(p7v, dpTrypsin, v7);
            Assert.AreEqual((true, true), p7_pep.IntersectsAndIdentifiesVariation(v7));

            (Protein p10v, var v10) = appliedMap["protein10"];
            var p10_pep = PickCoveringPeptide(p10v, dpTrypsin, v10);
            Assert.AreEqual((true, true), p10_pep.IntersectsAndIdentifiesVariation(v10));

            (Protein p1v, var v1) = appliedMap["protein1"];
            var p1_pep = PickCoveringPeptide(p1v, dpTrypsin, v1);
            Assert.AreEqual((true, true), p1_pep.IntersectsAndIdentifiesVariation(v1));

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

            //  chymotrypsin (don't cleave before proline)
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
        [Test]
        public static void TestReverseDecoyFromPeptideFromProteinXML()
        {
            //Just making sure there are no snafus when creating decoy peptides from an xml,which will have mods in various places, etc.
            //sequence variants, modifications
            Dictionary<string, Modification> un = new Dictionary<string, Modification>();
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "cRAP_databaseGPTMD.xml"),
                true, DecoyType.None, UniProtPtms, false, new string[] { "exclude_me" }, out un,
                maxSequenceVariantsPerIsoform: 4, minAlleleDepth: 1, maxSequenceVariantIsoforms: 1);

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
        [Test]
        public static void CountTargetsWithMatchingDecoys()
        {
            Dictionary<string, Modification> un = new Dictionary<string, Modification>();
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            // Pin legacy LoadProteinXML defaults to avoid new default behavior
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "cRAP_databaseGPTMD.xml"),
                true, DecoyType.None, UniProtPtms, false, new[] { "exclude_me" }, out un,
                maxSequenceVariantsPerIsoform: 4, minAlleleDepth: 1, maxSequenceVariantIsoforms: 1);

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
        [Test]
        public static void TestPeptideWithSetModsReturnsTruncationsInTopDown()
        {
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.xml");

            // Pin legacy variant-expansion defaults to restore previous behavior
            Protein insulin = ProteinDbLoader.LoadProteinXML(
                xmlDatabase,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: null,
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out var unknownModifications,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1,
                addTruncations: true)[0];

            Protease protease = new Protease("top-down", CleavageSpecificity.None, "", "", new List<DigestionMotif>(), null);
            List<PeptideWithSetModifications> insulinTruncations = insulin.Digest(new DigestionParams(protease: protease.Name), new List<Modification>(), new List<Modification>(), topDownTruncationSearch: true).ToList();
            Assert.AreEqual(68, insulinTruncations.Count);
        }
        [Test]
        public static void TestPeptideWithSetModsReturnsDecoyTruncationsInTopDown()
        {
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.xml");
            List<Protein> insulinProteins = ProteinDbLoader.LoadProteinXML(
                xmlDatabase,
                generateTargets: true,
                decoyType: DecoyType.Reverse,
                allKnownModifications: null,
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out var unknownModifications,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1,
                addTruncations: true);

            Protease protease = new Protease("top-down", CleavageSpecificity.None, "", "", new List<DigestionMotif>(), null);
            List<PeptideWithSetModifications> insulintTargetTruncations = insulinProteins.Where(p => !p.IsDecoy).First()
                .Digest(new DigestionParams(protease: protease.Name), new List<Modification>(), new List<Modification>(), topDownTruncationSearch: true).ToList();
            Assert.AreEqual(68, insulintTargetTruncations.Count);
            List<PeptideWithSetModifications> insulintDecoyTruncations = insulinProteins.Where(p => p.IsDecoy).First()
                .Digest(new DigestionParams(protease: protease.Name), new List<Modification>(), new List<Modification>(), topDownTruncationSearch: true).ToList();
            Assert.AreEqual(68, insulintDecoyTruncations.Count);
        }
        [Test]
        public static void CheckFullChemicalFormula()
        {
            PeptideWithSetModifications small_pep = new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION"), new DigestionParams(protease: "trypsin"), 1, 7, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            ChemicalFormula small_pep_cf = ChemicalFormula.ParseFormula("C34H53N7O15");
            Assert.AreEqual(small_pep.FullChemicalFormula, small_pep_cf);

            PeptideWithSetModifications large_pep = new PeptideWithSetModifications(new Protein("PEPTIDEKRNSPEPTIDEKECUEIRQUV", "ACCESSION"), new DigestionParams(protease: "trypsin"), 1, 28, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            ChemicalFormula large_pep_cf = ChemicalFormula.ParseFormula("C134H220N38O50S1Se2");
            Assert.AreEqual(large_pep.FullChemicalFormula, large_pep_cf);

            ModificationMotif.TryGetMotif("S", out ModificationMotif motif_s);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif_s, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            Dictionary<int, Modification> modDict_small = new Dictionary<int, Modification>();
            modDict_small.Add(4, phosphorylation);

            PeptideWithSetModifications small_pep_mod = new PeptideWithSetModifications(new Protein("PEPSIDE", "ACCESSION"), new DigestionParams(protease: "trypsin"), 1, 7, CleavageSpecificity.Full, null, 0, modDict_small, 0, null);
            ChemicalFormula small_pep_mod_cf = ChemicalFormula.ParseFormula("C33H52N7O18P1");
            Assert.AreEqual(small_pep_mod.FullChemicalFormula, small_pep_mod_cf);

            ModificationMotif.TryGetMotif("K", out ModificationMotif motif_k);
            Modification acetylation = new Modification(_originalId: "acetyl", _modificationType: "CommonBiological", _target: motif_k, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("C2H3O"));
            Dictionary<int, Modification> modDict_large = new Dictionary<int, Modification>();
            modDict_large.Add(4, phosphorylation);
            modDict_large.Add(11, phosphorylation);
            modDict_large.Add(8, acetylation);

            PeptideWithSetModifications large_pep_mod = new PeptideWithSetModifications(new Protein("PEPSIDEKRNSPEPTIDEKECUEIRQUV", "ACCESSION"), new DigestionParams(protease: "trypsin"), 1, 28, CleavageSpecificity.Full, null, 0, modDict_large, 0, null);
            ChemicalFormula large_pep_mod_cf = ChemicalFormula.ParseFormula("C135H223N38O57P2S1Se2");
            Assert.AreEqual(large_pep_mod.FullChemicalFormula, large_pep_mod_cf);

            ModificationMotif.TryGetMotif("C", out var motif_c);
            ModificationMotif.TryGetMotif("G", out var motif_g);
            Dictionary<string, Modification> modDict =
                new()
                {
                    { "Carbamidomethyl on C", new Modification(_originalId: "Carbamidomethyl", _modificationType: "Common Fixed",
                        _target: motif_c, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("C2H3ON")) },
                    { "BS on G" , new Modification(_originalId: "BS on G", _modificationType: "BS", _target: motif_g, _monoisotopicMass: 96.0875)}
                };
            PeptideWithSetModifications pwsmWithMissingCfMods = new PeptideWithSetModifications(
                "ENQGDETQG[Speculative:BS on G]C[Common Fixed:Carbamidomethyl on C]PPQR", modDict, p: new Protein("ENQGDETQGCPPQR", "FakeProtein"), digestionParams: new DigestionParams(),
                oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 14);
            Assert.Null(pwsmWithMissingCfMods.FullChemicalFormula);
        }

        [Test]
        public static void CheckMostAbundantMonoisotopicMass()
        {
            PeptideWithSetModifications small_pep = new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION"), new DigestionParams(protease: "trypsin"), 1, 7, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            double small_pep_most_abundant_mass_prospector = 800.36724 - 1.0079;
            Assert.That(small_pep.MostAbundantMonoisotopicMass, Is.EqualTo(small_pep_most_abundant_mass_prospector).Within(0.01));

            PeptideWithSetModifications large_pep = new PeptideWithSetModifications(new Protein("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", "ACCESSION"), new DigestionParams(protease: "trypsin"), 1, 42, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            double large_pep_most_abundant_mass_prospector = 4709.12020 - 1.0079;
            Assert.That(large_pep.MostAbundantMonoisotopicMass, Is.EqualTo(large_pep_most_abundant_mass_prospector).Within(0.01));
        }

        [Test]
        public static void TestPeptideWithSetModsEssentialSequence()
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            Dictionary<string, int> modsToWrite = new Dictionary<string, int>();
            modsToWrite.Add("UniProt", 0);

            var proteinXml = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "humanGAPDH.xml"),
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: UniProtPtms,
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out var unknownMod,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);
            var gapdh = proteinXml[0];

            var gapdhPeptides = gapdh.Digest(
                new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Variable),
                UniProtPtms,
                new List<Modification>());

            List<string> allSequences = new List<string>();
            foreach (var peptide in gapdhPeptides)
            {
                allSequences.Add(peptide.EssentialSequence(modsToWrite));
            }

            var expectedFullStrings = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "essentialSequences.txt"));

            CollectionAssert.AreEquivalent(expectedFullStrings, allSequences.ToArray());
        }
        [Test]
        public static void TestPeptideWithSetModsFullSequence()
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();
            var proteinXml = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "humanGAPDH.xml"),
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: UniProtPtms,
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out var unknownMod,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);
            var gapdh = proteinXml[0];

            var gapdhPeptides = gapdh.Digest(
                new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Variable),
                UniProtPtms,
                new List<Modification>());

            List<string> allSequences = new List<string>();
            foreach (var peptide in gapdhPeptides)
            {
                allSequences.Add(peptide.FullSequence);
            }

            var expectedFullStrings = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "fullSequences.txt"));
            CollectionAssert.AreEquivalent(expectedFullStrings, allSequences.ToArray());

            allSequences.Clear();
            foreach (var peptide in gapdhPeptides)
            {
                allSequences.Add(peptide.FullSequenceWithMassShift());
            }

            var expectedFullStringsWithMassShifts = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "fullSequencesWithMassShift.txt"));
            CollectionAssert.AreEquivalent(expectedFullStringsWithMassShifts, allSequences.ToArray());
        }
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

        [Test]
        public static void TestPeptideWithSetModsEquals()
        {
            // Create two proteins
            Protein protein1 = new Protein("SEQUENCEK", "accession1");
            Protein protein2 = new Protein("SEQUENCEK", "accession2");

            // Create digestion parameters
            DigestionParams digestionParams = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            // Digest the proteins to get peptides
            PeptideWithSetModifications peptide1 = protein1.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
            PeptideWithSetModifications peptide2 = protein2.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();

            // Test equality - same peptide
            Assert.IsTrue(peptide1.Equals(peptide1));

            // different peptide
            Assert.IsTrue(!peptide1.Equals(peptide2));
            Assert.IsTrue(!peptide1.Equals((object)peptide2));
            Assert.IsTrue(!peptide1.Equals((IBioPolymerWithSetMods)peptide2));
            Assert.AreNotEqual(peptide1.GetHashCode(), peptide2.GetHashCode());

            // Test inequality with different start residue
            PeptideWithSetModifications peptide3 = new PeptideWithSetModifications(protein1, digestionParams, 2, 9, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            Assert.IsFalse(peptide1.Equals(peptide3));

            // Test inequality with different parent accession
            PeptideWithSetModifications peptide4 = new PeptideWithSetModifications(protein2, digestionParams, 1, 9, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            Assert.IsFalse(peptide1.Equals(peptide4));

            // all fail on null
            Assert.That(!peptide1.Equals(null));
            Assert.That(!peptide1.Equals((object)null));
            Assert.That(!peptide1.Equals((PeptideWithSetModifications)null));
        }

        [Test]
        public static void TestIBioPolymerWithSetModsModificationFromFullSequence()
        {
            Dictionary<string, Modification> un = new Dictionary<string, Modification>();
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"),
                    formalChargesDictionary).ToList();

            // Pin legacy LoadProteinXML defaults to restore previous behavior
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "cRAP_databaseGPTMD.xml"),
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: UniProtPtms,
                isContaminant: false,
                modTypesToExclude: new string[] { "exclude_me" },
                unknownModifications: out un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);

            var allKnownModDict = UniProtPtms.ToDictionary(p => p.IdWithMotif, p => p);
            var digestionParameters = new DigestionParams(maxModsForPeptides: 3);

            foreach (Protein p in proteins)
            {
                List<PeptideWithSetModifications> digestedPeptides =
                    p.Digest(digestionParameters, [], [], null, null).ToList();
                // take the most modified peptide by base sequence and ensure all methods function properly
                foreach (var targetPeptide in digestedPeptides
                             .Where(pep => pep.FullSequence.Contains('['))
                             .GroupBy(pep => pep.BaseSequence)
                             .Select(pepGroup => pepGroup.MaxBy(pep => pep.AllModsOneIsNterminus.Count)))
                {
                    var startResidue = targetPeptide.OneBasedStartResidue;
                    var endResidue = targetPeptide.OneBasedEndResidue;

                    int expectedModCount = 0;
                    foreach (var modDictEntry in p.OneBasedPossibleLocalizedModifications
                                 .Where(mod => mod.Key >= startResidue && mod.Key <= endResidue))
                    {
                        if (modDictEntry.Value.Count > 1)
                        {
                            var locRestrictions = modDictEntry.Value.Select(mod => mod.LocationRestriction).ToList();

                            if (locRestrictions.AllSame())
                            {
                                if (locRestrictions.First() == "Anywhere.")
                                    expectedModCount++;
                                else if (locRestrictions.First() == "N-terminal." && modDictEntry.Key == startResidue)
                                    expectedModCount++;
                            }
                            else if (modDictEntry.Value.Select(mod => mod.LocationRestriction).Contains("Anywhere.")
                                     && modDictEntry.Value.Select(mod => mod.LocationRestriction)
                                         .Contains("N-terminal."))
                            {
                                expectedModCount++;
                                if (modDictEntry.Key == startResidue)
                                    expectedModCount++;
                            }
                        }
                        else
                        {
                            switch (modDictEntry.Value.First().LocationRestriction)
                            {
                                case "Anywhere.":
                                case "N-terminal." when modDictEntry.Key == startResidue:
                                    expectedModCount++;
                                    break;
                            }
                        }
                    }

                    expectedModCount = Math.Min(expectedModCount, digestionParameters.MaxMods);

                    var expectedModifications = p.OneBasedPossibleLocalizedModifications.Where(mod =>
                        mod.Key >= startResidue &&
                        mod.Key <= endResidue).SelectMany(mod => mod.Value).ToList();

                    // Parse modifications from PWSM and two IBioPolymerWithSetMods methods
                    var pwsmModDict = targetPeptide.AllModsOneIsNterminus;
                    var bpwsmModDict = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(targetPeptide.FullSequence, allKnownModDict);
                    var bpwsmModList = IBioPolymerWithSetMods.GetModificationsFromFullSequence(targetPeptide.FullSequence, allKnownModDict);

                    // Ensure all methods are in agreement by modification count
                    Assert.AreEqual(pwsmModDict.Count, expectedModCount);
                    Assert.AreEqual(bpwsmModDict.Count, expectedModCount);
                    Assert.AreEqual(bpwsmModList.Count, expectedModCount);

                    // Ensure all methods are in agreement by modification identify
                    foreach (var pwsmModification in pwsmModDict.Values)
                        Assert.Contains(pwsmModification, expectedModifications);
                    foreach (var pwsmModification in bpwsmModDict.Values)
                        Assert.Contains(pwsmModification, expectedModifications);
                    foreach (var pwsmModification in bpwsmModList)
                        Assert.Contains(pwsmModification, expectedModifications);
                }
            }
        }
        [Test]
        public static void TestGetSubstitutedFullSequence()
        {
            //It should take care of multiple substitutions
            string test1 = "F[1 nucleotide substitution:F->Y on F]SIMGGGLA[1 nucleotide substitution:A->S on A]DR";
            string expected1 = "YSIMGGGLSDR";
            var actual1 = IBioPolymerWithSetMods.ParseSubstitutedFullSequence(test1);
            Assert.That(actual1, Is.EqualTo(expected1));

            //It should not change other modifications
            string test2 = "SANH[1 nucleotide substitution:H->L on H]M[Common Variable:Oxidation on M]AGHWVAISGAAGGLGSLAVQYAK";
            string expected2 = "SANLM[Common Variable:Oxidation on M]AGHWVAISGAAGGLGSLAVQYAK";
            var actual2 = IBioPolymerWithSetMods.ParseSubstitutedFullSequence(test2);
            Assert.That(actual2, Is.EqualTo(expected2));

            //It should work on 2 nucleotide substitutions
            string test3 = "S[2+ nucleotide substitution:S->E on S]AAADRLNLTSGHLNAGR";
            string expected3 = "EAAADRLNLTSGHLNAGR";
            var actual3 = IBioPolymerWithSetMods.ParseSubstitutedFullSequence(test3);
            Assert.That(actual3, Is.EqualTo(expected3));
        }
        private static SequenceVariation MakePointVariant(int pos, char original, char variant)
            => new SequenceVariation(
                oneBasedBeginPosition: pos,
                oneBasedEndPosition: pos,
                originalSequence: original.ToString(),
                variantSequence: variant.ToString(),
                description: $"{original}{pos}{variant}");

        private static Protein MakeOriginalProtein(string seq, string accession = "P1")
            => new Protein(sequence: seq, accession: accession);

        private static Protein MakeVariantProtein(Protein original, string variantSequence, SequenceVariation variation)
            => new Protein(variantSequence, original, new[] { variation }, applicableProteolysisProducts: new List<TruncationProduct>(),
                           oneBasedModifications: new Dictionary<int, List<Modification>>(), sampleNameForVariants: null);

        [Test]
        public static void IntersectsAndIdentifiesVariation_NewCTermCleavageSite_SetsIdentifiesTrue()
        {
            // Original sequence (position 5 = A, not a trypsin cleavage residue)
            // Index: 1 2 3 4 5 6 7 8 9
            //        P E P T A I D E K
            string originalSeq = "PEPTAIDEK";
            var originalProtein = MakeOriginalProtein(originalSeq);

            // Variant changes A5 -> K5 creating a new potential C-terminal cleavage site before peptide start
            var variation = MakePointVariant(5, 'A', 'K');
            string variantSeq = "PEPTKIDEK";
            var variantProtein = MakeVariantProtein(originalProtein, variantSeq, variation);

            // Peptide starts immediately after the variant (residues 6-8: IDE)
            var dp = new DigestionParams(protease: "trypsin");
            var peptide = new PeptideWithSetModifications(
                variantProtein,
                dp,
                oneBasedStartResidueInProtein: 6,
                oneBasedEndResidueInProtein: 8,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "test",
                missedCleavages: 0,
                allModsOneIsNterminus: new Dictionary<int, Modification>(),
                numFixedMods: 0,
                baseSequence: "IDE");

            // Act
            var (intersects, identifies) = peptide.IntersectsAndIdentifiesVariation(variation);

            // Assert: variant is immediately upstream (no intersection) but creates a new cleavage site => identifies == true
            Assert.Multiple(() =>
            {
                Assert.That(intersects, Is.False, "Expected no positional overlap with the variant");
                Assert.That(identifies, Is.True, "Expected identification of new upstream cleavage site (A->K)");
            });
        }

        [Test]
        public static void IntersectsAndIdentifiesVariation_NoNewCleavageSite_IdentifiesFalse()
        {
            // Original sequence (position 5 = A)
            string originalSeq = "PEPTAIDEK";
            var originalProtein = MakeOriginalProtein(originalSeq);

            // Variant changes A5 -> V5 (neither A nor V is a trypsin cleavage residue => no new site)
            var variation = MakePointVariant(5, 'A', 'V');
            string variantSeq = "PEPTVIDEK";
            var variantProtein = MakeVariantProtein(originalProtein, variantSeq, variation);

            var dp = new DigestionParams(protease: "trypsin");
            var peptide = new PeptideWithSetModifications(
                variantProtein,
                dp,
                oneBasedStartResidueInProtein: 6,
                oneBasedEndResidueInProtein: 8,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "test-noneg",
                missedCleavages: 0,
                allModsOneIsNterminus: new Dictionary<int, Modification>(),
                numFixedMods: 0,
                baseSequence: "IDE");

            var (intersects, identifies) = peptide.IntersectsAndIdentifiesVariation(variation);

            Assert.Multiple(() =>
            {
                Assert.That(intersects, Is.False, "Expected no intersection");
                Assert.That(identifies, Is.False, "No new cleavage site introduced (A->V) so identifies should be false");
            });
        }
        // Helper: build original protein
        private static Protein MakeProtein(string seq, string acc = "PVAR") => new Protein(seq, acc);

        // Helper: apply variation to produce variant base sequence
        private static (SequenceVariation variation, string variantBase) MakeDeletionVariation(
            string originalSeq, int begin, int end, string variantInserted)
        {
            string originalSegment = originalSeq.Substring(begin - 1, end - begin + 1);
            string prefix = originalSeq.Substring(0, begin - 1);
            string suffix = originalSeq.Substring(end); // after end
            string variantBase = prefix + variantInserted + suffix;

            var sv = new SequenceVariation(
                oneBasedBeginPosition: begin,
                oneBasedEndPosition: end,
                originalSequence: originalSegment,
                variantSequence: variantInserted,
                description: $"del_{begin}_{end}_len{variantInserted.Length}");

            return (sv, variantBase);
        }

        private static PeptideWithSetModifications MakePeptide(
            Protein variantProtein,
            int start,
            int end,
            string baseSeq,
            DigestionParams dp)
        {
            return new PeptideWithSetModifications(
                variantProtein,
                dp,
                oneBasedStartResidueInProtein: start,
                oneBasedEndResidueInProtein: end,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "test-pep",
                missedCleavages: 0,
                allModsOneIsNterminus: new Dictionary<int, Modification>(),
                numFixedMods: 0,
                baseSequence: baseSeq);
        }

        private const string OriginalProteinSeq = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"; // length 40

        // Matrix of scenarios:
        // EVC (effectiveVariantEnd correction) & effectiveDegenerate combinations

        [Test]
        public static void IntersectsAndIdentifiesVariation_FullDeletion_EVCTrue_DegenerateTrue()
        {
            // Deletion remove 10-20 entirely (variant sequence empty)
            int begin = 10;
            int end = 20;

            var originalProtein = MakeProtein(OriginalProteinSeq);
            var (variation, variantBase) = MakeDeletionVariation(OriginalProteinSeq, begin, end, variantInserted: "");
            // Variant protein (shorter by 11 aa)
            var variantProtein = new Protein(originalProtein, variantBase);

            // Peptide starts AFTER the corrected effectiveVariantEnd (= begin) so degenerate
            // In variant coordinates: positions after deletion are compressed.
            // Choose start 15 end 18 (no actual overlap in effective span → degenerate).
            var dp = new DigestionParams(protease: "trypsin");

            // Derive base sequence from variant
            string pepBase = variantBase.Substring(15 - 1, 18 - 15 + 1);
            var peptide = MakePeptide(variantProtein, 15, 18, pepBase, dp);

            var (intersects, identifies) = peptide.IntersectsAndIdentifiesVariation(variation);

            Assert.Multiple(() =>
            {
                Assert.That(intersects, Is.True, "Deletion path still reports intersects tuple true.");
                Assert.That(identifies, Is.True, "Full deletion sets identifiesFlag true.");
            });
        }

        [Test]
        public static void IntersectsAndIdentifiesVariation_FullDeletion_EVCTrue_DegenerateFalse()
        {
            int begin = 10;
            int end = 20;
            var originalProtein = MakeProtein(OriginalProteinSeq);
            var (variation, variantBase) = MakeDeletionVariation(OriginalProteinSeq, begin, end, variantInserted: "");
            var variantProtein = new Protein(originalProtein, variantBase);
            var dp = new DigestionParams(protease: "trypsin");

            // Peptide spans original prefix (variant coords 9..11)
            // start 9 -> before deletion; end 11 -> after junction (compressed) ensures intersectEndEff == startEff (not degenerate)
            string pepBase = variantBase.Substring(9 - 1, 11 - 9 + 1);
            var peptide = MakePeptide(variantProtein, 9, 11, pepBase, dp);

            var (intersects, identifies) = peptide.IntersectsAndIdentifiesVariation(variation);

            Assert.Multiple(() =>
            {
                Assert.That(intersects, Is.True);
                Assert.That(identifies, Is.True, "Deletion still marks identifiesFlag.");
            });
        }

        [Test]
        public static void IntersectsAndIdentifiesVariation_PartialDeletion_EVCFalse_DegenerateTrue()
        {
            int begin = 10;
            int end = 20;
            // Partial deletion: replace 11-length region with 5 aa
            string inserted = "KLMNP";
            var originalProtein = MakeProtein(OriginalProteinSeq);
            var (variation, variantBase) = MakeDeletionVariation(OriginalProteinSeq, begin, end, inserted);
            var variantProtein = new Protein(originalProtein, variantBase);
            var dp = new DigestionParams(protease: "trypsin");

            // Choose peptide start AFTER effectiveVariantEnd (which will be end + (lenDiff) = 20 -6 =14)
            // Variant coordinate 15..17 -> degenerate (intersectEndEff < intersectStartEff)
            string pepBase = variantBase.Substring(15 - 1, 17 - 15 + 1);
            var peptide = MakePeptide(variantProtein, 15, 17, pepBase, dp);

            var (intersects, identifies) = peptide.IntersectsAndIdentifiesVariation(variation);

            Assert.Multiple(() =>
            {
                Assert.That(intersects, Is.True);
                Assert.That(identifies, Is.True, "Deletion (partial) sets identifiesFlag.");
            });
        }

        [Test]
        public static void IntersectsAndIdentifiesVariation_PartialDeletion_EVCFalse_DegenerateFalse()
        {
            int begin = 10;
            int end = 20;
            string inserted = "KLMNP";
            var originalProtein = MakeProtein(OriginalProteinSeq);
            var (variation, variantBase) = MakeDeletionVariation(OriginalProteinSeq, begin, end, inserted);
            var variantProtein = new Protein(originalProtein, variantBase);
            var dp = new DigestionParams(protease: "trypsin");

            // Peptide 9..12 (variant coords) => intersects effective variant span (effectiveVariantEnd=14) producing non-degenerate overlap
            string pepBase = variantBase.Substring(9 - 1, 12 - 9 + 1);
            var peptide = MakePeptide(variantProtein, 9, 12, pepBase, dp);

            var (intersects, identifies) = peptide.IntersectsAndIdentifiesVariation(variation);

            Assert.Multiple(() =>
            {
                Assert.That(intersects, Is.True);
                Assert.That(identifies, Is.True);
            });
        }
    }
}