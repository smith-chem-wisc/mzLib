﻿using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
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
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams(minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new Protein("PROTEIN", "myAccession");

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
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams("Lys-C (don't cleave before proline)", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new Protein("MKPKPKPMKA", "myAccession");

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
            Protease multiletter = new Protease("multiletter", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(multiletter.Name, multiletter);

            DigestionParams myDigestionParams = new DigestionParams("multiletter", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new Protein("PROPRPPM", "myAccession");

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
            Protease multiletter = new Protease("multi-custom", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(multiletter.Name, multiletter);

            DigestionParams myDigestionParams = new DigestionParams("multi-custom", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new Protein("FAAAMAAM", "myAccession");

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
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams("collagenase", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new Protein("ABCGPXGPMFKCGPMKK", "myAccession");

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
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams("Asp-N", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new Protein("PADDMSKDPDMMAASMDJSSM", "myAccession");

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
            Protease protease = new Protease("lys-c", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);

            DigestionParams myDigestionParams = new DigestionParams("lys-c", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new Protein("PROKPKMKP", "myAccession");

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
            DigestionParams myDigestionParams = new DigestionParams("chymotrypsin (don't cleave before proline)", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new Protein("AASFPWDJSSMF", "myAccession");

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
            DigestionParams myDigestionParams = new DigestionParams("non-specific", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new Protein("PRO", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            Assert.AreEqual(myPeptides.Count(), 3);
        }

        [Test]
        public static void TestSpecificProteaseWriting()
        {
            //check that the specific protease is the one written for indexing
            //This is needed for speedy non-specific searches to have the stratified full/semi/none peptide cleavages
            //If the protease is written instead of the specific protease, then the specific protease information is lost upon deserialization.

            //check for nonspecific
            DigestionParams dp = new DigestionParams(protease: "Arg-C", searchModeType: CleavageSpecificity.None);
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
        public static void TestDigestionParamsSerializeDeserialize()
        {
            DigestionParams weirdDigestionParams = new DigestionParams("Asp-N", 77, 88, 99, 69, InitiatorMethionineBehavior.Cleave, 420, CleavageSpecificity.Unknown, Proteomics.Fragmentation.FragmentationTerminus.None, false);
            string serializedString = weirdDigestionParams.ToString();
            DigestionParams deserializedDigestionParams = DigestionParams.FromString(serializedString);
            Assert.AreEqual(weirdDigestionParams, deserializedDigestionParams);
        }

        [Test]
        public static void TestEndSequenceNTerm()
        {
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams("Lys-N", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new Protein("KKPROTEIN", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual(first, "K");
            Assert.AreEqual(myPeptides.Count(), 2);
        }

        [Test]
        public static void TestOneMotifMultiplePreventing()
        {
            var empty = new List<Modification>();
            var digestionmotifs = DigestionMotif.ParseDigestionMotifsFromString("N[M]|,N[C]|,N[A]|");
            Protease customProtease = new Protease("custom", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(customProtease.Name, customProtease);

            DigestionParams myDigestionParams = new DigestionParams("custom", minPeptideLength: 1, maxMissedCleavages: 0);

            // create a protein
            Protein myProtein = new Protein("PRONFNMMHFHAA", "myAccession");

            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();

            Assert.AreEqual(myPeptides.Count, 2);
            Assert.AreEqual(first, "PRON");
        }

        [Test]
        public static void TestProteolyticDigestion()
        {
            Protein humanInsulin = new Protein("MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN", "P01308",
                proteolysisProducts: new List<ProteolysisProduct>
                {
                    new ProteolysisProduct(1, 24, ""),
                    new ProteolysisProduct(25, 54, ""),
                    new ProteolysisProduct(57, 87, ""),
                    new ProteolysisProduct(90, 110, "")
                });
            DigestionParams dp = new DigestionParams(maxMissedCleavages: 10, minPeptideLength: 1, maxPeptideLength: 120); //this should allow for all peptides to be generated
            List<PeptideWithSetModifications> pwsms = humanInsulin.Digest(dp, null, null).ToList();
            HashSet<PeptideWithSetModifications> hashset = new HashSet<PeptideWithSetModifications>(pwsms);

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
            DigestionParams speedySemiN = new DigestionParams("trypsin", 10, 29, 30, 1024, InitiatorMethionineBehavior.Retain, 2, CleavageSpecificity.Semi, Proteomics.Fragmentation.FragmentationTerminus.N);
            DigestionParams speedySemiC = new DigestionParams("trypsin", 10, 29, 30, 1024, InitiatorMethionineBehavior.Retain, 2, CleavageSpecificity.Semi, Proteomics.Fragmentation.FragmentationTerminus.C);
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