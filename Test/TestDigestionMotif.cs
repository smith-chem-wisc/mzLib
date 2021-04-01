using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestDigestionMotif
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setup()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void TestParseProtease()
        {
            var argn = DigestionMotif.ParseDigestionMotifsFromString("|D");
            Assert.AreEqual(argn.Count, 1);

            var c = argn[0];
            Assert.AreEqual(c.InducingCleavage, "D");
            Assert.AreEqual(c.PreventingCleavage, null);
            Assert.AreEqual(c.CutIndex, 0);

            var chymotrypsin = DigestionMotif.ParseDigestionMotifsFromString("F[P]|,W[P]|,Y[P]|");
            Assert.AreEqual(chymotrypsin.Count, 3);
        }

        [Test]
        public static void TestBasicProtease1()
        {
            List<Modification> empty = new();
            DigestionParams myDigestionParams = new(minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new("PROTEIN", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual(first, "PR");
            Assert.AreEqual(last, "OTEIN");
        }

        [Test]
        public static void TestBasicProtease2()
        {
            List<Modification> empty = new();
            DigestionParams myDigestionParams = new("Lys-C (don't cleave before proline)", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new("MKPKPKPMKA", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual("MKPKPKPMK", first);
            Assert.AreEqual("A", last);
        }

        [Test]
        public static void TestWildCardExclusion()
        {
            var empty = new List<Modification>();
            var digestionmotifs = DigestionMotif.ParseDigestionMotifsFromString("RX{P}|");
            Protease multiletter = new("multiletter", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(multiletter.Name, multiletter);

            DigestionParams myDigestionParams = new("multiletter", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new("PROPRPPM", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual(first, "PRO");
            Assert.AreEqual(last, "PRPPM");
        }

        [Test]
        public static void TestMultiLetterNTerm()
        {
            var empty = new List<Modification>();
            var digestionmotifs = DigestionMotif.ParseDigestionMotifsFromString("|AAA");
            Protease multiletter = new("multi-custom", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(multiletter.Name, multiletter);

            DigestionParams myDigestionParams = new("multi-custom", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new("FAAAMAAM", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual(first, "F");
            Assert.AreEqual(last, "AAAMAAM");
        }

        [Test]
        public static void TestMultiLetterProtease()
        {
            List<Modification> empty = new();
            DigestionParams myDigestionParams = new("collagenase", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new("ABCGPXGPMFKCGPMKK", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual(first, "ABCGPX");
            Assert.AreEqual(last, "GPMFKCGPMKK");
        }

        [Test]
        public static void TestNTerminusProtease()
        {
            List<Modification> empty = new();
            DigestionParams myDigestionParams = new("Asp-N", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new("PADDMSKDPDMMAASMDJSSM", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual(first, "PA");
            Assert.AreEqual(last, "DJSSM");
        }

        [Test]
        public static void TestWrongSyntax()
        {
            Assert.Throws<MzLibException>(() =>
            {
                var protease = DigestionMotif.ParseDigestionMotifsFromString("X[Y,P]");
                Assert.Fail("Exception shold be thrown for incorrect syntax.");
            });
        }

        [Test]
        public static void TestCutIndexDifferentSyntax()
        {
            var empty = new List<Modification>();
            var digestionmotifs = DigestionMotif.ParseDigestionMotifsFromString("K|[P]"); // same as K[P]|
            Protease protease = new("lys-c", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);

            DigestionParams myDigestionParams = new("lys-c", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new("PROKPKMKP", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual(first, "PROKPK");
            Assert.AreEqual(last, "MKP");
        }

        [Test]
        public static void TestEndSequenceCTerm()
        {
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new("chymotrypsin (don't cleave before proline)", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new("AASFPWDJSSMF", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual(first, "AASFPW");
            Assert.AreEqual(last, "DJSSMF");
        }

        [Test]
        public static void TestNonSpecificProtease()
        {
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new("non-specific", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new("PRO", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            Assert.AreEqual(myPeptides.Count, 3);
        }

        [Test]
        public static void TestSpecificProteaseWriting()
        {
            //check that the specific protease is the one written for indexing
            //This is needed for speedy non-specific searches to have the stratified full/semi/none peptide cleavages
            //If the protease is written instead of the specific protease, then the specific protease information is lost upon deserialization.

            //check for nonspecific
            DigestionParams dp = new(protease: "Arg-C", searchModeType: CleavageSpecificity.None);
            string proteaseString = dp.ToString().Split(',')[6];
            Assert.IsTrue(proteaseString.Equals("Arg-C"));

            //Check for semi
            dp = new DigestionParams(protease: "Arg-C", searchModeType: CleavageSpecificity.Semi);
            proteaseString = dp.ToString().Split(',')[6];
            Assert.IsTrue(proteaseString.Equals("Arg-C"));

            //check for normal
            dp = new DigestionParams(protease: "Arg-C"); //default searchModeType is Full
            proteaseString = dp.ToString().Split(',')[6];
            Assert.IsTrue(proteaseString.Equals("Arg-C"));
        }

        [Test]
        public static void TestEndSequenceNTerm()
        {
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new("Lys-N", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new("KKPROTEIN", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual(first, "K");
            Assert.AreEqual(myPeptides.Count, 2);
            Assert.IsNotNull(last);
        }

        [Test]
        public static void TestOneMotifMultiplePreventing()
        {
            List<Modification> empty = new();
            var digestionmotifs = DigestionMotif.ParseDigestionMotifsFromString("N[M]|,N[C]|,N[A]|");
            Protease customProtease = new("custom", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(customProtease.Name, customProtease);

            DigestionParams myDigestionParams = new("custom", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new("PRONFNMMHFHAA", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual(myPeptides.Count, 2);
            Assert.AreEqual(first, "PRON");
            Assert.IsNotNull(last);
        }

        [Test]
        public static void TestProteolyticDigestion()
        {
            Protein humanInsulin = new("MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN", "P01308",
                proteolysisProducts: new List<ProteolysisProduct>
                {
                    new ProteolysisProduct(1, 24, ""),
                    new ProteolysisProduct(25, 54, ""),
                    new ProteolysisProduct(57, 87, ""),
                    new ProteolysisProduct(90, 110, "")
                });
            DigestionParams dp = new(maxMissedCleavages: 10, minPeptideLength: 1, maxPeptideLength: 120); //this should allow for all peptides to be generated
            List<PeptideWithSetModifications> pwsms = humanInsulin.Digest(dp, null, null).ToList();
            HashSet<PeptideWithSetModifications> hashset = new(pwsms);

            //check that all the full length proteolysis products were made
            Assert.IsTrue(pwsms.Any(x => x.OneBasedStartResidueInProtein == 1 && x.OneBasedEndResidueInProtein == 24));
            Assert.IsTrue(pwsms.Any(x => x.OneBasedStartResidueInProtein == 25 && x.OneBasedEndResidueInProtein == 54));
            Assert.IsTrue(pwsms.Any(x => x.OneBasedStartResidueInProtein == 57 && x.OneBasedEndResidueInProtein == 87));
            Assert.IsTrue(pwsms.Any(x => x.OneBasedStartResidueInProtein == 90 && x.OneBasedEndResidueInProtein == 110));

            //check that all the correct peptides were made
            Assert.IsTrue(hashset.Count == 52);

            //check that there are no duplicates
            Assert.IsTrue(pwsms.Count == hashset.Count);

            //Speedy semi specific test
            DigestionParams speedySemiN = new("trypsin", 10, 29, 30, 1024, InitiatorMethionineBehavior.Retain, 2, CleavageSpecificity.Semi, Proteomics.Fragmentation.FragmentationTerminus.N);
            DigestionParams speedySemiC = new("trypsin", 10, 29, 30, 1024, InitiatorMethionineBehavior.Retain, 2, CleavageSpecificity.Semi, Proteomics.Fragmentation.FragmentationTerminus.C);
            List<PeptideWithSetModifications> pwsmsN = humanInsulin.Digest(speedySemiN, null, null).ToList();
            List<PeptideWithSetModifications> pwsmsC = humanInsulin.Digest(speedySemiC, null, null).ToList();
            Assert.IsTrue(pwsmsN.Count == 7);
            Assert.IsTrue(pwsmsC.Count == 9);
            Assert.IsFalse(pwsmsN.Any(x => x.Length > speedySemiN.MaxPeptideLength));
            Assert.IsFalse(pwsmsC.Any(x => x.Length > speedySemiC.MaxPeptideLength));
            Assert.IsFalse(pwsmsN.Any(x => x.Length < speedySemiN.MinPeptideLength));
            Assert.IsFalse(pwsmsC.Any(x => x.Length < speedySemiC.MinPeptideLength));
        }
    }
}