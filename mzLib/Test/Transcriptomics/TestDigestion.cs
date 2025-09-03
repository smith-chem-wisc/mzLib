using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace Test.Transcriptomics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestDigestion
    {
        public record RnaDigestionTestCase(string BaseSequence, string Enzyme, int MissedCleavages, int MinLength,
            int MaxLength, int DigestionProductCount,
            double[] MonoMasses, string[] Sequences);

        public static IEnumerable<RnaDigestionTestCase> GetTestCases()
        {
            // 6bp Top Down
            yield return new RnaDigestionTestCase("GUACUG", "top-down",
                0, 1, 6, 1,
                new[] { 1874.28 },
                new[] { "GUACUG" });
            // 6bp Rnase T1, normal
            yield return new RnaDigestionTestCase("GUACUG", "RNase T1",
                0, 1, 6, 2,
                new[] { 363.057, 1529.234 },
                new[] { "G", "UACUG" });
            // 6bp Cusativin, normal
            yield return new RnaDigestionTestCase("GUACUG", "Cusativin",
                0, 1, 6, 2,
                new[] { 1303.175, 589.116 },
                new[] { "GUAC", "UG" });
            // 6bp Rnase T1, one product too short
            yield return new RnaDigestionTestCase("GUACUG", "RNase T1",
                0, 3, 6, 1,
                new[] { 1529.234 },
                new[] { "UACUG" });
            // 6bp Rnase T1, one product too long
            yield return new RnaDigestionTestCase("GUACUG", "RNase T1",
                0, 1, 2, 1,
                new[] { 363.057 },
                new[] { "G" });
            // 6bp Rnase T1, 1 missed cleavage
            yield return new RnaDigestionTestCase("GUACUG", "RNase T1",
                1, 1, 6, 3,
                new[] { 363.057, 1529.234, 1874.28 },
                new[] { "G", "UACUG", "GUACUG" });
            // 6bp Rnase A
            yield return new RnaDigestionTestCase("GUACUG", "RNase A",
                0, 1, 6, 4,
                new[] { 669.082, 652.103, 324.035, 283.091 },
                new[] { "GU", "AC", "U", "G" });
            // 6bp Rnase A, 1 missed cleavage
            yield return new RnaDigestionTestCase("GUACUG", "RNase A",
                1, 1, 6, 7,
                new[] { 669.082, 652.103, 324.035, 283.091, 1303.175, 958.128, 589.116 },
                new[] { "GU", "AC", "U", "G", "GUAC", "ACU", "UG" });
            // 6bp Rnase A, 2 missed cleavages
            yield return new RnaDigestionTestCase("GUACUG", "RNase A",
                2, 1, 6, 9,
                new[] { 669.082, 652.103, 324.035, 283.091, 1303.175, 958.128, 589.116, 1609.200, 1223.209 },
                new[] { "GU", "AC", "U", "G", "GUAC", "ACU", "UG", "GUACU", "ACUG" });
            // 20bp top-down
            yield return new RnaDigestionTestCase("GUACUGCCUCUAGUGAAGCA", "top-down",
                0, 1, int.MaxValue, 1,
                new[] { 6363.871 },
                new[] { "GUACUGCCUCUAGUGAAGCA" });
            // 20bp Rnase T1, normal
            yield return new RnaDigestionTestCase("GUACUGCCUCUAGUGAAGCA", "RNase T1",
                0, 1, int.MaxValue, 6,
                new[] { 363.057, 1609.200, 2219.282, 669.082, 1021.161, 572.137 },
                new[] { "G", "UACUG", "CCUCUAG", "UG", "AAG", "CA" });
        }

        public static string rnaseTsvpath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"Digestion\rnases.tsv");

        [OneTimeSetUp]
        public void OneTimeSetup()
        {
            RnaseDictionary.Dictionary = RnaseDictionary.LoadRnaseDictionary(rnaseTsvpath);
        }

        #region Rnase

        [Test]
        public void TestRnaseDictionaryLoading()
        {
            var rnaseCountFromTsv = File.ReadAllLines(rnaseTsvpath).Length - 1;
            Assert.That(RnaseDictionary.Dictionary.Count, Is.EqualTo(rnaseCountFromTsv));
        }

        [Test]
        [TestCaseSource(nameof(GetTestCases))]
        public void TestRnase_GetUnmodifiedOligos_Counts(RnaDigestionTestCase testCase)
        {
            RNA rna = new RNA(testCase.BaseSequence);
            Rnase rnase = RnaseDictionary.Dictionary[testCase.Enzyme];
            var digestionProducts =
                rnase.GetUnmodifiedOligos(rna, testCase.MissedCleavages, testCase.MinLength, testCase.MaxLength);

            Assert.That(digestionProducts.Count(), Is.EqualTo(testCase.DigestionProductCount));
        }

        [Test]
        [TestCaseSource(nameof(GetTestCases))]
        public void TestRnase_GetUnmodifiedOligo_Sequence(RnaDigestionTestCase testCase)
        {
            RNA rna = new RNA(testCase.BaseSequence);
            Rnase rnase = RnaseDictionary.Dictionary[testCase.Enzyme];
            var digestionProducts =
                rnase.GetUnmodifiedOligos(rna, testCase.MissedCleavages, testCase.MinLength, testCase.MaxLength).ToList();

            Assert.That(digestionProducts.Count, Is.EqualTo(testCase.Sequences.Length));
            for (var i = 0; i < digestionProducts.Count; i++)
            {
                var product = digestionProducts[i];
                var testCaseCaseSequence = testCase.Sequences[i];
                Assert.That(product.BaseSequence == testCaseCaseSequence);
            }
        }

        [Test]
        public void TestRnaseEqualityProperties()
        {
            Rnase t1 = RnaseDictionary.Dictionary["RNase T1"];
            Rnase t1Duplicate = RnaseDictionary.Dictionary["RNase T1"];
            Rnase t2 = RnaseDictionary.Dictionary["RNase T2"];

            Assert.That(t1.Equals(t1Duplicate));
            Assert.That(t1.Equals(t1));
            Assert.That(!t1.Equals(t2));
            Assert.That(!t1.Equals(null));
            Assert.That(t1.GetHashCode(), Is.EqualTo(t1Duplicate.GetHashCode()));
            Assert.That(t1.GetHashCode(), Is.Not.EqualTo(t2.GetHashCode()));
            Assert.That(t1.Equals((object)t1Duplicate));
            Assert.That(t1.Equals((object)t1));
            Assert.That(!t1.Equals((object)t2));
            Assert.That(!t1.Equals((object)null));
            // ReSharper disable once SuspiciousTypeConversion.Global
            Assert.That(!t1.Equals((object)new RNA("GUA")));
        }

        [Test]
        public void TestRnase_UnmodifiedOligos_Exception()
        {
            Rnase rnase = new Rnase("Bad", CleavageSpecificity.SingleC, new List<DigestionMotif>());
            Assert.Throws<ArgumentException>(() => { rnase.GetUnmodifiedOligos(new RNA("GUACUG"), 0, 1, 6); });
        }

        #endregion

        #region NucleolyticOligo

        [Test]
        public void TestNucleolyticOligoProperties_FivePrimeDigestionProduct()
        {
            RNA rna = new("GUACUG");
            Rnase rnase = RnaseDictionary.Dictionary["RNase U2"];
            var digestionProducts = rnase.GetUnmodifiedOligos(rna, 0, 1, 6).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(3));

            var oligo = digestionProducts[0];
            Assert.That(oligo.BaseSequence, Is.EqualTo("G"));
            Assert.That(oligo.OneBasedStartResidue, Is.EqualTo(1));
            Assert.That(oligo.OneBasedEndResidue, Is.EqualTo(1));
            Assert.That(oligo.MissedCleavages, Is.EqualTo(0));
            Assert.That(oligo.CleavageSpecificityForFdrCategory, Is.EqualTo(CleavageSpecificity.Full));
            Assert.That(oligo.NextResidue, Is.EqualTo('U'));
            Assert.That(oligo.PreviousResidue, Is.EqualTo('-'));
            Assert.That(oligo.ToString(), Is.EqualTo(oligo.BaseSequence));
        }

        [Test]
        public void TestNucleolyticOligoProperties_ThreePrimeDigestionProduct()
        {
            RNA rna = new("GUACUG");
            Rnase rnase = RnaseDictionary.Dictionary["RNase U2"];
            var digestionProducts = rnase.GetUnmodifiedOligos(rna, 0, 1, 6).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(3));

            NucleolyticOligo oligo = digestionProducts[2];
            Assert.That(oligo.BaseSequence, Is.EqualTo("CUG"));
            Assert.That(oligo.OneBasedStartResidue, Is.EqualTo(4));
            Assert.That(oligo.OneBasedEndResidue, Is.EqualTo(6));
            Assert.That(oligo.MissedCleavages, Is.EqualTo(0));
            Assert.That(oligo.CleavageSpecificityForFdrCategory, Is.EqualTo(CleavageSpecificity.Full));
            Assert.That(oligo.NextResidue, Is.EqualTo('-'));
            Assert.That(oligo.PreviousResidue, Is.EqualTo('A'));
            Assert.That(oligo.ToString(), Is.EqualTo(oligo.BaseSequence));
        }

        [Test]
        public void TestNucleolyticOligoProperties_InternalDigestionProduct()
        {
            RNA rna = new("GUACUG");
            Rnase rnase = RnaseDictionary.Dictionary["RNase U2"];
            var digestionProducts = rnase.GetUnmodifiedOligos(rna, 0, 1, 6).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(3));

            NucleolyticOligo oligo = digestionProducts[1];
            Assert.That(oligo.BaseSequence, Is.EqualTo("UA"));
            Assert.That(oligo.OneBasedStartResidue, Is.EqualTo(2));
            Assert.That(oligo.OneBasedEndResidue, Is.EqualTo(3));
            Assert.That(oligo.MissedCleavages, Is.EqualTo(0));
            Assert.That(oligo.CleavageSpecificityForFdrCategory, Is.EqualTo(CleavageSpecificity.Full));
            Assert.That(oligo.NextResidue, Is.EqualTo('C'));
            Assert.That(oligo.PreviousResidue, Is.EqualTo('G'));
            Assert.That(oligo.ToString(), Is.EqualTo(oligo.BaseSequence));
        }

        [Test]
        public void TestNucleolyticOligoProperties_TopDownDigestionProduct()
        {
            RNA rna = new("GUACUG");
            Rnase rnase = RnaseDictionary.Dictionary["top-down"];
            var digestionProducts = rnase.GetUnmodifiedOligos(rna, 0, 1, 6).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(1));

            NucleolyticOligo oligo = digestionProducts[0];
            Assert.That(oligo.BaseSequence, Is.EqualTo("GUACUG"));
            Assert.That(oligo.OneBasedStartResidue, Is.EqualTo(1));
            Assert.That(oligo.OneBasedEndResidue, Is.EqualTo(6));
            Assert.That(oligo.MissedCleavages, Is.EqualTo(0));
            Assert.That(oligo.CleavageSpecificityForFdrCategory, Is.EqualTo(CleavageSpecificity.Full));
            Assert.That(oligo.NextResidue, Is.EqualTo('-'));
            Assert.That(oligo.PreviousResidue, Is.EqualTo('-'));
            Assert.That(oligo.ToString(), Is.EqualTo(oligo.BaseSequence));
        }

        #endregion

        #region OligoWithSetMods

        private static (string Sequence, int FragmentNumber, ProductType Type, double Mass)[] DigestFragmentTestCases =>
            new (string Sequence, int FragmentNumber, ProductType Type, double Mass)[]
            {
                ("UAG", 0, ProductType.M, 998.134),
                ("UAG", 1, ProductType.aBaseLoss, 114.031), ("UAG", 2, ProductType.aBaseLoss, 420.056),
                ("UAG", 1, ProductType.c, 308.031), ("UAG", 2, ProductType.c, 637.093),
                ("UAG", 1, ProductType.dWaterLoss, 306.025), ("UAG", 2, ProductType.dWaterLoss, 635.077),
                ("UAG", 1, ProductType.w, 443.023), ("UAG", 2, ProductType.w, 772.075),
                ("UAG", 1, ProductType.y,  363.057), ("UAG", 2, ProductType.y, 692.109),
                ("UAG", 1, ProductType.yWaterLoss,  345.047), ("UAG", 2, ProductType.yWaterLoss, 674.100),

                ("UCG", 0, ProductType.M, 974.123),
                ("UCG", 1, ProductType.aBaseLoss, 114.031), ("UCG", 2, ProductType.aBaseLoss, 420.056),
                ("UCG", 1, ProductType.c, 308.040), ("UCG", 2, ProductType.c, 613.082),
                ("UCG", 1, ProductType.dWaterLoss, 306.025), ("UCG", 2, ProductType.dWaterLoss, 611.066),
                ("UCG", 1, ProductType.w, 443.023), ("UCG", 2, ProductType.w, 748.064),
                ("UCG", 1, ProductType.y,  363.057), ("UCG", 2, ProductType.y, 668.098),
                ("UCG", 1, ProductType.yWaterLoss,  345.047), ("UCG", 2, ProductType.yWaterLoss, 650.089),

                ("UUG", 0, ProductType.M, 975.107),
                ("UUG", 1, ProductType.aBaseLoss, 114.031), ("UUG", 2, ProductType.aBaseLoss, 420.056),
                ("UUG", 1, ProductType.c, 308.041), ("UUG", 2, ProductType.c, 614.066),
                ("UUG", 1, ProductType.dWaterLoss, 306.025), ("UUG", 2, ProductType.dWaterLoss, 612.050),
                ("UUG", 1, ProductType.w, 443.023), ("UUG", 2, ProductType.w, 749.048),
                ("UUG", 1, ProductType.y,  363.057), ("UUG", 2, ProductType.y, 669.082),
                ("UUG", 1, ProductType.yWaterLoss,  345.047), ("UUG", 2, ProductType.yWaterLoss, 651.073),

                ("AUAG", 0, ProductType.M, 1247.220),
                ("AUAG", 1, ProductType.aBaseLoss, 114.031), ("AUAG", 2, ProductType.aBaseLoss, 443.083), ("AUAG", 3, ProductType.aBaseLoss, 749.108),
                ("AUAG", 1, ProductType.c, 331.068), ("AUAG", 2, ProductType.c, 637.093), ("AUAG", 3, ProductType.c, 966.146),
                ("AUAG", 1, ProductType.dWaterLoss, 329.052), ("AUAG", 2, ProductType.dWaterLoss, 635.077), ("AUAG", 3, ProductType.dWaterLoss, 964.129),
                ("AUAG", 1, ProductType.w, 363.057), ("AUAG", 2, ProductType.w, 692.109), ("AUAG", 3, ProductType.w, 998.134),
                ("AUAG", 1, ProductType.y,  283.091), ("AUAG", 2, ProductType.y, 612.143), ("AUAG", 3, ProductType.y, 918.168),
                ("AUAG", 1, ProductType.yWaterLoss,  265.081), ("AUAG", 2, ProductType.yWaterLoss, 594.134), ("AUAG", 3, ProductType.yWaterLoss, 900.159),
            };

        [Test] // test values calculated with http://rna.rega.kuleuven.be/masspec/mongo.htm
        [TestCase("UAGUCGUUGAUAG", 4140.555, new[] { "UAG", "UCG", "UUG", "AUAG" },
            new[] { 998.134, 974.123, 975.107, 1247.220 })]
        public static void TestDigestionAndFragmentation(string sequence, double monoMass,
            string[] digestionProductSequences, double[] digestionProductMasses)
        {
            RNA rna = new(sequence);
            Assert.That(rna.MonoisotopicMass, Is.EqualTo(monoMass).Within(0.01));

            // digest RNA
            var digestionParams = new RnaDigestionParams("RNase T1");
            var products = rna.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .Select(p => (OligoWithSetMods)p).ToList();
            Assert.That(products.Count, Is.EqualTo(digestionProductSequences.Length));

            // ensure digestion sequence and masses are correct
            for (var index = 0; index < products.Count; index++)
            {
                var digestionProduct = products[index];
                Assert.That(digestionProduct.BaseSequence, Is.EqualTo(digestionProductSequences[index]));
                Assert.That(digestionProduct.MonoisotopicMass, Is.EqualTo(digestionProductMasses[index]).Within(0.01));

                List<Product> fragments = new();
                digestionProduct.Fragment(DissociationType.CID, FragmentationTerminus.Both, fragments);

                // test that fragments are correct
                var fragmentsToCompare = DigestFragmentTestCases
                    .Where(p => p.Sequence.Equals(digestionProduct.BaseSequence)).ToList();
                for (var i = 0; i < fragments.Count; i++)
                {
                    var fragment = fragments[i];
                    var theoreticalFragment = fragmentsToCompare.FirstOrDefault(p =>
                        p.FragmentNumber == fragment.FragmentNumber && p.Type == fragment.ProductType);
                    if (theoreticalFragment.Mass is 0.0 ) continue;
                    Assert.That(fragment.MonoisotopicMass, Is.EqualTo(theoreticalFragment.Mass).Within(0.01));
                    Assert.That(fragment.FragmentNumber, Is.EqualTo(theoreticalFragment.FragmentNumber));
                    Assert.That(fragment.ProductType, Is.EqualTo(theoreticalFragment.Type));
                    Assert.That(fragment.FragmentNumber, Is.EqualTo(theoreticalFragment.FragmentNumber));
                    if (fragment.Terminus == FragmentationTerminus.FivePrime)
                        Assert.That(fragment.AminoAcidPosition, Is.EqualTo(theoreticalFragment.FragmentNumber));
                    else if (fragment.Terminus == FragmentationTerminus.None)
                        Assert.That(fragment.FragmentNumber, Is.EqualTo(0));
                    else
                        Assert.That(fragment.AminoAcidPosition, Is.EqualTo(digestionProductSequences[index].Length - theoreticalFragment.FragmentNumber));
                }
            }
        }

        [Test]
        [TestCase("UAGUCGUUGAUAG", new[] { "UAG", "UCG", "UUG", "AUAG" },
            new[] { 1, 4, 7, 10 }, new[] { 3, 6, 9, 13 }, new[] { '-', 'G', 'G', 'G' },
            new[] { 'U', 'U', 'A', '-' })]
        public static void TestOligoWithSetMods_AAPositions(string sequence, string[] digestionProductSequences,
        int[] startResidue, int[] endResidue, char[] preciousResidue, char[] nextResidue)
        {
            RNA rna = new RNA(sequence);
            var digestionProducts = rna.Digest(new RnaDigestionParams("RNase T1"), new List<Modification>(),
                new List<Modification>()).Select(p => (OligoWithSetMods)p).ToList();

            Assert.That(digestionProducts.All(p => p.DigestionParams.DigestionAgent.Name == "RNase T1"));
            for (var index = 0; index < digestionProducts.Count; index++)
            {
                var digestionProduct = digestionProducts[index];
                Assert.That(digestionProduct.BaseSequence, Is.EqualTo(digestionProductSequences[index]));
                Assert.That(digestionProduct.OneBasedStartResidue, Is.EqualTo(startResidue[index]));
                Assert.That(digestionProduct.OneBasedEndResidue, Is.EqualTo(endResidue[index]));
                Assert.That(digestionProduct.PreviousResidue, Is.EqualTo(preciousResidue[index]));
                Assert.That(digestionProduct.NextResidue, Is.EqualTo(nextResidue[index]));
            }
        }

        [Test]
        public static void TestTermini_ThreePrimeCyclicPhosphate()
        {
            string sequence = "UAGUCGUUGAUAG";
            RNA rna = new RNA(sequence);
            var oligoCyclicPhosphate = PtmListLoader.ReadModsFromString(
                "ID   Cyclic Phosphate\r\nTG   X\r\nPP   Oligo 3'-terminal.\r\nMT   Digestion Termini\r\nCF   H-2 O-1\r\nDR   Unimod; 280.\r\n//",
                out List<(Modification, string)> errors).First();
            var nucleicAcidCyclicPhosphate = PtmListLoader.ReadModsFromString(
                "ID   Cyclic Phosphate\r\nTG   X\r\nPP   3'-terminal.\r\nMT   Digestion Termini\r\nCF   H-2 O-1\r\nDR   Unimod; 280.\r\n//",
                out errors).First();
            Assert.That(!errors.Any());

            // top-down digestion, 3' terminal modification
            var variableMods = new List<Modification> { nucleicAcidCyclicPhosphate };
            var digestionParams = new RnaDigestionParams("top-down");
            var digestionProducts = rna.Digest(digestionParams, new List<Modification>(), variableMods)
                .Select(p => (OligoWithSetMods)p).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(2));
            Assert.That(digestionProducts[0].FullSequence, Is.EqualTo("UAGUCGUUGAUAG"));
            Assert.That(digestionProducts[0].SequenceWithChemicalFormulas, Is.EqualTo("UAGUCGUUGAUAG"));
            Assert.That(digestionProducts[0].FullSequenceWithMassShift(), Is.EqualTo("UAGUCGUUGAUAG"));
            
            Assert.That(digestionProducts[1].FullSequence, Is.EqualTo("UAGUCGUUGAUAG-[Digestion Termini:Cyclic Phosphate on X]"));
            Assert.That(digestionProducts[1].SequenceWithChemicalFormulas, Is.EqualTo("UAGUCGUUGAUAG-[H-2O-1]"));
            Assert.That(digestionProducts[1].FullSequenceWithMassShift(), Is.EqualTo("UAGUCGUUGAUAG-[-18.010565]"));

            // top-down digestion, 3' oligo terminal modification
            variableMods = new List<Modification> { oligoCyclicPhosphate };
            digestionProducts = rna.Digest(digestionParams, new List<Modification>(), variableMods)
                .Select(p => (OligoWithSetMods)p).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(2));
            Assert.That(digestionProducts[0].FullSequence, Is.EqualTo("UAGUCGUUGAUAG"));
            Assert.That(digestionProducts[1].FullSequence, Is.EqualTo("UAGUCGUUGAUAG-[Digestion Termini:Cyclic Phosphate on X]"));

            // RNase T1 digestion, 3' terminal modification
            digestionParams = new RnaDigestionParams("RNase T1");
            variableMods = new List<Modification> { nucleicAcidCyclicPhosphate };
            digestionProducts = rna.Digest(digestionParams, new List<Modification>(), variableMods)
                .Select(p => (OligoWithSetMods)p).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(5));
            var expected = new List<string>()
            {
                "UAG", "UCG", "UUG", "AUAG", "AUAG-[Digestion Termini:Cyclic Phosphate on X]"
            };
            for (int i = 0; i < expected.Count; i++)
            {
                Assert.That(digestionProducts[i].FullSequence, Is.EqualTo(expected[i]));
            }

            // RNase T1 digestion, 3' oligo terminal modification 
            variableMods = new List<Modification> { oligoCyclicPhosphate };
            digestionProducts = rna.Digest(digestionParams, new List<Modification>(), variableMods)
                .Select(p => (OligoWithSetMods)p).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(8));
            expected = new List<string>()
            {
                "UAG", "UAG-[Digestion Termini:Cyclic Phosphate on X]",
                "UCG", "UCG-[Digestion Termini:Cyclic Phosphate on X]",
                "UUG", "UUG-[Digestion Termini:Cyclic Phosphate on X]",
                "AUAG","AUAG-[Digestion Termini:Cyclic Phosphate on X]"
            };

            for (int i = 0; i < expected.Count; i++)
            {
                Assert.That(digestionProducts[i].FullSequence, Is.EqualTo(expected[i]));
            }
        }

        [Test]
        public static void TestTermini_FivePrimeLargeMod()
        {
            string sequence = "UAGUCGUUGAUAG";
            RNA rna = new RNA(sequence);
            var oligoLargeMod = PtmListLoader.ReadModsFromString(
                "ID   Pfizer 5'-Cap\r\nTG   X\r\nPP   Oligo 5'-terminal.\r\nMT   Standard\r\nCF   C13H22N5O14P3\r\nDR   Unimod; 280.\r\n//",
                out List<(Modification, string)> errors).First();
            var nucleicAcidLargeMod = PtmListLoader.ReadModsFromString(
                "ID   Pfizer 5'-Cap\r\nTG   X\r\nPP   5'-terminal.\r\nMT   Standard\r\nCF   C13H22N5O14P3\r\nDR   Unimod; 280.\r\n//",
                out errors).First();
            Assert.That(!errors.Any());

            // top-down digestion, 5' terminal modification, expect two products
            var variableMods = new List<Modification> { nucleicAcidLargeMod };
            var digestionParams = new RnaDigestionParams("top-down");
            var digestionProducts = rna.Digest(digestionParams, new List<Modification>(), variableMods)
                .Select(p => (OligoWithSetMods)p).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(2));
            Assert.That(digestionProducts[0].FullSequence, Is.EqualTo("UAGUCGUUGAUAG"));
            Assert.That(digestionProducts[1].FullSequence, Is.EqualTo("[Standard:Pfizer 5'-Cap on X]UAGUCGUUGAUAG"));

            // top-down digestion, 5' oligo terminal modification, expect two products
            variableMods = new List<Modification> { oligoLargeMod };
            digestionProducts = rna.Digest(digestionParams, new List<Modification>(), variableMods)
                .Select(p => (OligoWithSetMods)p).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(2));
            Assert.That(digestionProducts[0].FullSequence, Is.EqualTo("UAGUCGUUGAUAG"));
            Assert.That(digestionProducts[1].FullSequence, Is.EqualTo("[Standard:Pfizer 5'-Cap on X]UAGUCGUUGAUAG"));

            // RNase T1 digestion, 5' terminal modification
            digestionParams = new RnaDigestionParams("RNase T1");
            variableMods = new List<Modification> { nucleicAcidLargeMod };
            digestionProducts = rna.Digest(digestionParams, new List<Modification>(), variableMods)
                .Select(p => (OligoWithSetMods)p).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(5));
            var expected = new List<string>()
            {
                "UAG", "[Standard:Pfizer 5'-Cap on X]UAG", "UCG", "UUG", "AUAG"
            };
            for (int i = 0; i < expected.Count; i++)
            {
                Assert.That(digestionProducts[i].FullSequence, Is.EqualTo(expected[i]));
            }

            // RNase T1 digestion, 5' oligo terminal modification 
            variableMods = new List<Modification> { oligoLargeMod };
            digestionProducts = rna.Digest(digestionParams, new List<Modification>(), variableMods)
                .Select(p => (OligoWithSetMods)p).ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(8));
            expected = new List<string>()
            {
                "UAG", "[Standard:Pfizer 5'-Cap on X]UAG",
                "UCG", "[Standard:Pfizer 5'-Cap on X]UCG",
                "UUG", "[Standard:Pfizer 5'-Cap on X]UUG",
                "AUAG", "[Standard:Pfizer 5'-Cap on X]AUAG"
            };

            for (int i = 0; i < expected.Count; i++)
            {
                Assert.That(digestionProducts[i].FullSequence, Is.EqualTo(expected[i]));
            }
        }

        [Test]
        [TestCase("UAGUCGUUGAUAG")]
        public static void TestOligoWithSetMods_PropertiesWithTopDownDigestion(string sequence)
        {
            var rna = new RNA(sequence);
            var oligoWithSetMods =
                rna.Digest(new RnaDigestionParams(), new List<Modification>(), new List<Modification>())
                        .First() as OligoWithSetMods ?? throw new NullReferenceException();

            Assert.That(rna.BaseSequence, Is.EqualTo(oligoWithSetMods.BaseSequence));
            Assert.That(rna.ThreePrimeTerminus, Is.EqualTo(oligoWithSetMods.ThreePrimeTerminus));
            Assert.That(rna.FivePrimeTerminus, Is.EqualTo(oligoWithSetMods.FivePrimeTerminus));
            Assert.That(rna.ThisChemicalFormula, Is.EqualTo(oligoWithSetMods.ThisChemicalFormula));
            Assert.That(rna.Length, Is.EqualTo(oligoWithSetMods.Length));
        }

        [Test]
        public static void OligoWithSetMods_CalculatedValues()
        {
            var rna = new RNA("GUACUG");
            var rnaFormula = rna.ThisChemicalFormula;

            string modText = "ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A\r\nCF   Na1H-1\r\n" + @"//";
            var sodiumAdduct = PtmListLoader.ReadModsFromString(modText, out List<(Modification, string)> mods).First();
            var oligoWithSetMods =
                rna.Digest(new RnaDigestionParams(), new List<Modification>() { sodiumAdduct }, new List<Modification>())
                    .First() as OligoWithSetMods ?? throw new NullReferenceException();

            Assert.That(oligoWithSetMods.NumMods, Is.EqualTo(1));
            Assert.That(oligoWithSetMods.NumFixedMods, Is.EqualTo(1));
            Assert.That(oligoWithSetMods.NumVariableMods, Is.EqualTo(0));
            Assert.That(oligoWithSetMods.CleavageSpecificityForFdrCategory, Is.EqualTo(CleavageSpecificity.Full));

            var formula = oligoWithSetMods.ThisChemicalFormula;
            Assert.That(formula, Is.EqualTo(rnaFormula + sodiumAdduct.ChemicalFormula));

            var formulaToAdd = ChemicalFormula.ParseFormula("H");
            var deltaMass = formulaToAdd.MonoisotopicMass;
            var oldMonoMass = oligoWithSetMods.MonoisotopicMass;
            var oldMostAbundantMass = oligoWithSetMods.MostAbundantMonoisotopicMass;

            oligoWithSetMods.FivePrimeTerminus = formulaToAdd + oligoWithSetMods.FivePrimeTerminus;

            Assert.That(oligoWithSetMods.MonoisotopicMass, Is.EqualTo(oldMonoMass + deltaMass).Within(0.01));
            Assert.That(oligoWithSetMods.MostAbundantMonoisotopicMass, Is.EqualTo(oldMostAbundantMass + deltaMass).Within(0.01));
            Assert.That(oligoWithSetMods.ThisChemicalFormula, Is.EqualTo(formula + formulaToAdd));

            oldMonoMass = oligoWithSetMods.MonoisotopicMass;
            oldMostAbundantMass = oligoWithSetMods.MostAbundantMonoisotopicMass;
            oligoWithSetMods.ThreePrimeTerminus = formulaToAdd + oligoWithSetMods.ThreePrimeTerminus;

            Assert.That(oligoWithSetMods.MonoisotopicMass, Is.EqualTo(oldMonoMass + deltaMass).Within(0.01));
            Assert.That(oligoWithSetMods.MostAbundantMonoisotopicMass, Is.EqualTo(oldMostAbundantMass + deltaMass).Within(0.01));
            Assert.That(oligoWithSetMods.ThisChemicalFormula, Is.EqualTo(formula + formulaToAdd + formulaToAdd));

            Assert.Throws<MzLibException>(() =>
            {
                var oligo = new OligoWithSetMods("GUA|GAUGUC", new Dictionary<string, Modification>());
            });
        }

        #endregion

        #region DigestionParams

        [Test]
        [TestCaseSource(nameof(GetTestCases))]
        public void TestDigestionParams_Properties(RnaDigestionTestCase testCase)
        {
            var rna = new RNA(testCase.BaseSequence);
            var digestionParams = new RnaDigestionParams(testCase.Enzyme, testCase.MissedCleavages, testCase.MinLength,
                testCase.MaxLength);

            Assert.That(digestionParams.DigestionAgent, Is.EqualTo(RnaseDictionary.Dictionary[testCase.Enzyme]));
            Assert.That(digestionParams.MaxMissedCleavages, Is.EqualTo(testCase.MissedCleavages));
            Assert.That(digestionParams.MinLength, Is.EqualTo(testCase.MinLength));
            Assert.That(digestionParams.MaxLength, Is.EqualTo(testCase.MaxLength));

            digestionParams.MaxModificationIsoforms = 2048;
            digestionParams.MaxMods = 3;
            Assert.That(digestionParams.MaxModificationIsoforms, Is.EqualTo(2048));
            Assert.That(digestionParams.MaxMods, Is.EqualTo(3));

            var digestionProducts = rna.Digest(digestionParams, new List<Modification>(), new List<Modification>());
            Assert.That(digestionProducts.Count(), Is.EqualTo(testCase.DigestionProductCount));
        }

        [Test]
        public void TestDigestionParamsClone()
        {
            var digestionParams = new RnaDigestionParams("top-down", 0, 3, 20000);
            var cloned = digestionParams.Clone(FragmentationTerminus.C);

            // set new terminus, all values except terminus are retained
            Assert.That(digestionParams.DigestionAgent, Is.EqualTo(cloned.DigestionAgent));
            Assert.That(digestionParams.MaxMissedCleavages, Is.EqualTo(cloned.MaxMissedCleavages));
            Assert.That(digestionParams.MinLength, Is.EqualTo(cloned.MinLength));
            Assert.That(digestionParams.MaxLength, Is.EqualTo(cloned.MaxLength));
            Assert.That(digestionParams.MaxMods, Is.EqualTo(cloned.MaxMods));
            Assert.That(digestionParams.FragmentationTerminus, Is.Not.EqualTo(cloned.FragmentationTerminus));
            Assert.That(digestionParams.SearchModeType, Is.EqualTo(CleavageSpecificity.Full));
            Assert.That(cloned.FragmentationTerminus, Is.EqualTo(FragmentationTerminus.C));

            // do not set new terminus, all values are retained
            cloned = digestionParams.Clone();
            Assert.That(digestionParams.DigestionAgent, Is.EqualTo(cloned.DigestionAgent));
            Assert.That(digestionParams.MaxMissedCleavages, Is.EqualTo(cloned.MaxMissedCleavages));
            Assert.That(digestionParams.MinLength, Is.EqualTo(cloned.MinLength));
            Assert.That(digestionParams.MaxLength, Is.EqualTo(cloned.MaxLength));
            Assert.That(digestionParams.MaxMods, Is.EqualTo(cloned.MaxMods));
            Assert.That(digestionParams.FragmentationTerminus, Is.EqualTo(cloned.FragmentationTerminus));
            Assert.That(digestionParams.SearchModeType, Is.EqualTo(CleavageSpecificity.Full));
            Assert.That(cloned.FragmentationTerminus, Is.EqualTo(FragmentationTerminus.Both));
        }

        #endregion

        #region NucleicAcid


        [Test]
        [TestCaseSource(nameof(GetTestCases))]
        public void TestNucleicAcid_Digestion_WithoutMods_Counts(RnaDigestionTestCase testCase)
        {
            var rna = new RNA(testCase.BaseSequence);
            var digestionParams = new RnaDigestionParams(testCase.Enzyme, testCase.MissedCleavages, testCase.MinLength,
                testCase.MaxLength);

            var digestionProducts = rna.Digest(digestionParams, new List<Modification>(), new List<Modification>());
            Assert.That(digestionProducts.Count(), Is.EqualTo(testCase.DigestionProductCount));
        }

        [Test]
        [TestCaseSource(nameof(GetTestCases))]
        public void TestNucleicAcid_Digestion_WithoutMods_Sequences(RnaDigestionTestCase testCase)
        {
            var rna = new RNA(testCase.BaseSequence);
            var digestionParams = new RnaDigestionParams(testCase.Enzyme, testCase.MissedCleavages, testCase.MinLength,
                testCase.MaxLength);

            var digestionProducts = rna.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .ToList();

            Assert.That(digestionProducts.Count, Is.EqualTo(testCase.Sequences.Length));
            for (var i = 0; i < digestionProducts.Count; i++)
            {
                var product = digestionProducts[i];
                var testCaseCaseSequence = testCase.Sequences[i];
                Assert.That(product.BaseSequence, Is.EqualTo(testCaseCaseSequence));
                Assert.That(product.FullSequence, Is.EqualTo(testCaseCaseSequence));
            }
        }

        [Test]
        [TestCaseSource(nameof(GetTestCases))]
        public void TestNucleicAcid_Digestion_WithoutMods_MonoMasses(RnaDigestionTestCase testCase)
        {
            var rna = new RNA(testCase.BaseSequence);
            var digestionParams = new RnaDigestionParams(testCase.Enzyme, testCase.MissedCleavages, testCase.MinLength,
                testCase.MaxLength);

            var digestionProducts = rna.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .ToList();

            Assert.That(digestionProducts.Count, Is.EqualTo(testCase.Sequences.Length));
            for (var i = 0; i < digestionProducts.Count; i++)
            {
                var productMass = digestionProducts[i].MonoisotopicMass;
                var testCaseCaseMass = testCase.MonoMasses[i];
                Assert.That(productMass, Is.EqualTo(testCaseCaseMass).Within(0.01));
            }
        }

        [Test]
        public static void TestNucleicAcid_Digestion_Exception()
        {
            IDigestionParams digestionParams = new Proteomics.ProteolyticDigestion.DigestionParams();
            var rna = new RNA("GUACUGGUACUG");

            try
            {
                var result = rna.Digest(digestionParams, new List<Modification>(), new List<Modification>());
            }
            catch (Exception e)
            {
                Assert.That(e, Is.TypeOf<MzLibException>());
                Assert.That(e.InnerException, Is.TypeOf<ArgumentException>());
            }
        }

        #endregion

        #region Digestion with Modifications

        public static List<Modification> SodiumAdducts =>
            PtmListLoader.ReadModsFromString("ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A or C or G or U\r\nCF   Na1H-1\r\n" + @"//",
                    out List<(Modification, string)> mods).ToList();

        public static List<Modification> PotassiumAdducts =>
            PtmListLoader.ReadModsFromString("ID   Potassium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A or C or G or U\r\nCF   K1H-1\r\n" + @"//",
                    out List<(Modification, string)> mods).ToList();

        public static List<Modification> TerminalSodiumAdducts =>
            PtmListLoader.ReadModsFromString("ID   Sodium\r\nMT   Metal\r\nPP   3'-terminal.\r\nTG   A or C or G or U\r\nCF   Na1H-1\r\n" + @"//",
            out List<(Modification, string)> mods).ToList();

        public static List<Modification> TerminalPotassiumAdducts =>
            PtmListLoader.ReadModsFromString("ID   Potassium\r\nMT   Metal\r\nPP   5'-terminal.\r\nTG   A or C or G or U\r\nCF   K1H-1\r\n" + @"//",
            out List<(Modification, string)> mods).ToList();

        [Test]
        public static void TestVariableModsCountCorrect()
        {
            var rna = new RNA("GUACUG");
            var rnaDigestionParams = new RnaDigestionParams()
            {
                MaxMods = 1,
            };

            var precursors = rna.Digest(rnaDigestionParams, new List<Modification>(), SodiumAdducts)
                .ToList();
            Assert.That(precursors.Count, Is.EqualTo(7));
            var fullSequences = precursors.Select(p => p.FullSequence).ToList();
            Assert.That(fullSequences.Contains("GUACUG"));
            Assert.That(fullSequences.Contains("G[Metal:Sodium on G]UACUG"));
            Assert.That(fullSequences.Contains("GU[Metal:Sodium on U]ACUG"));
            Assert.That(fullSequences.Contains("GUA[Metal:Sodium on A]CUG"));
            Assert.That(fullSequences.Contains("GUAC[Metal:Sodium on C]UG"));
            Assert.That(fullSequences.Contains("GUACU[Metal:Sodium on U]G"));
            Assert.That(fullSequences.Contains("GUACUG[Metal:Sodium on G]"));

            rnaDigestionParams.MaxMods = 2;
            precursors = rna.Digest(rnaDigestionParams, new List<Modification>(), SodiumAdducts)
                .ToList();
            Assert.That(precursors.Count, Is.EqualTo(22));
            fullSequences = precursors.Select(p => p.FullSequence).ToList();
            Assert.That(fullSequences.Contains("GUACUG"));
            Assert.That(fullSequences.Contains("G[Metal:Sodium on G]UACUG"));
            Assert.That(fullSequences.Contains("G[Metal:Sodium on G]U[Metal:Sodium on U]ACUG"));
            Assert.That(fullSequences.Contains("G[Metal:Sodium on G]UA[Metal:Sodium on A]CUG"));
            Assert.That(fullSequences.Contains("G[Metal:Sodium on G]UAC[Metal:Sodium on C]UG"));
            Assert.That(fullSequences.Contains("G[Metal:Sodium on G]UACU[Metal:Sodium on U]G"));
            Assert.That(fullSequences.Contains("G[Metal:Sodium on G]UACUG[Metal:Sodium on G]"));
            Assert.That(fullSequences.Contains("GU[Metal:Sodium on U]ACUG"));
            Assert.That(fullSequences.Contains("GU[Metal:Sodium on U]A[Metal:Sodium on A]CUG"));
            Assert.That(fullSequences.Contains("GU[Metal:Sodium on U]AC[Metal:Sodium on C]UG"));
            Assert.That(fullSequences.Contains("GU[Metal:Sodium on U]ACU[Metal:Sodium on U]G"));
            Assert.That(fullSequences.Contains("GU[Metal:Sodium on U]ACUG[Metal:Sodium on G]"));
            Assert.That(fullSequences.Contains("GUA[Metal:Sodium on A]CUG"));
            Assert.That(fullSequences.Contains("GUA[Metal:Sodium on A]C[Metal:Sodium on C]UG"));
            Assert.That(fullSequences.Contains("GUA[Metal:Sodium on A]CU[Metal:Sodium on U]G"));
            Assert.That(fullSequences.Contains("GUA[Metal:Sodium on A]CUG[Metal:Sodium on G]"));
            Assert.That(fullSequences.Contains("GUAC[Metal:Sodium on C]UG"));
            Assert.That(fullSequences.Contains("GUAC[Metal:Sodium on C]U[Metal:Sodium on U]G"));
            Assert.That(fullSequences.Contains("GUAC[Metal:Sodium on C]UG[Metal:Sodium on G]"));
            Assert.That(fullSequences.Contains("GUACU[Metal:Sodium on U]G"));
            Assert.That(fullSequences.Contains("GUACU[Metal:Sodium on U]G[Metal:Sodium on G]"));
            Assert.That(fullSequences.Contains("GUACUG[Metal:Sodium on G]"));
        }

        [Test]
        public static void TestFixedModsCountCorrect()
        {
            var sodiumAdduct = new List<Modification>() { SodiumAdducts[0] };

            var rna = new RNA("GUACUG");
            var rnaDigestionParams = new RnaDigestionParams()
            {
                MaxMods = 1,
            };
            var precursors = rna.Digest(rnaDigestionParams, sodiumAdduct, new List<Modification>())
                .ToList();
            Assert.That(precursors.Count, Is.EqualTo(1));
            Assert.That(precursors.First().NumFixedMods, Is.EqualTo(1));
            Assert.That(precursors.First().FullSequence, Is.EqualTo("GUA[Metal:Sodium on A]CUG"));
            Assert.That(precursors.First().MonoisotopicMass, Is.EqualTo(1896.26).Within(0.01));

            sodiumAdduct = new List<Modification>() { SodiumAdducts[2] };

            precursors = rna.Digest(rnaDigestionParams, sodiumAdduct, new List<Modification>())
                .ToList();
            Assert.That(precursors.Count, Is.EqualTo(1));
            Assert.That(precursors.First().NumFixedMods, Is.EqualTo(2));
            Assert.That(precursors.First().FullSequence, Is.EqualTo("G[Metal:Sodium on G]UACUG[Metal:Sodium on G]"));
            Assert.That(precursors.First().MonoisotopicMass, Is.EqualTo(1918.25).Within(0.01));
        }

        [Test]
        public static void TestFixedAndVariableMods()
        {
            var rna = new RNA("GUACUG");
            var rnaDigestionParams = new RnaDigestionParams();

            rnaDigestionParams.MaxMods = 1;
            var fixedMods = new List<Modification> { PotassiumAdducts[0] }; // A
            var variableMods = new List<Modification> { SodiumAdducts[1] }; // C
            var precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                .ToList();

            var fullSequences = precursors.Select(p => p.FullSequence).ToList();
            Assert.That(precursors.Count, Is.EqualTo(2));
            Assert.That(precursors.All(p => p.NumFixedMods == 1));
            Assert.That(fullSequences.Contains("GUA[Metal:Potassium on A]CUG"));
            Assert.That(fullSequences.Contains("GUA[Metal:Potassium on A]C[Metal:Sodium on C]UG"));

            var oneOfEach = precursors.First(p => p.FullSequence.Equals("GUA[Metal:Potassium on A]C[Metal:Sodium on C]UG"));
            Assert.That(oneOfEach.NumFixedMods, Is.EqualTo(1));
            Assert.That(oneOfEach.NumVariableMods, Is.EqualTo(1));
            Assert.That(oneOfEach.NumMods, Is.EqualTo(2));

            fixedMods = new List<Modification> { PotassiumAdducts[2] }; // G
            variableMods = new List<Modification> { SodiumAdducts[1] }; // C
            precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                .ToList();
            fullSequences = precursors.Select(p => p.FullSequence).ToList();
            Assert.That(precursors.Count, Is.EqualTo(2));
            Assert.That(precursors.All(p => p.NumFixedMods == 2));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UACUG[Metal:Potassium on G]"));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UAC[Metal:Sodium on C]UG[Metal:Potassium on G]"));

            fixedMods = new List<Modification> { PotassiumAdducts[2] }; // G
            variableMods = new List<Modification> { SodiumAdducts[1], SodiumAdducts[3] }; // C, U
            precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                .ToList();
            fullSequences = precursors.Select(p => p.FullSequence).ToList();
            Assert.That(precursors.Count, Is.EqualTo(4));
            Assert.That(precursors.All(p => p.NumFixedMods == 2));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UACUG[Metal:Potassium on G]"));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UAC[Metal:Sodium on C]UG[Metal:Potassium on G]"));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UACU[Metal:Sodium on U]G[Metal:Potassium on G]"));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]U[Metal:Sodium on U]ACUG[Metal:Potassium on G]"));

            rnaDigestionParams.MaxMods = 2;
            precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                .ToList();
            fullSequences = precursors.Select(p => p.FullSequence).ToList();
            Assert.That(precursors.Count, Is.EqualTo(7));
            Assert.That(precursors.All(p => p.NumFixedMods == 2));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UACUG[Metal:Potassium on G]"));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UAC[Metal:Sodium on C]UG[Metal:Potassium on G]"));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]U[Metal:Sodium on U]ACUG[Metal:Potassium on G]"));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UACU[Metal:Sodium on U]G[Metal:Potassium on G]"));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]U[Metal:Sodium on U]ACU[Metal:Sodium on U]G[Metal:Potassium on G]"));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UAC[Metal:Sodium on C]U[Metal:Sodium on U]G[Metal:Potassium on G]"));
            Assert.That(fullSequences.Contains("G[Metal:Potassium on G]U[Metal:Sodium on U]AC[Metal:Sodium on C]UG[Metal:Potassium on G]"));
        }

        /// <summary>
        /// Test when one fixed and one variable mod are used and share a localization
        /// expect two results, one with the fixed, and one with the variable
        /// </summary>
        [Test]
        public static void TestFixedAndVariableMods_LocalizationOverlap()
        {
            var rna = new RNA("GUACUG");
            var rnaDigestionParams = new RnaDigestionParams();

            for (int i = 1; i < 3; i++)
            {
                rnaDigestionParams.MaxMods = i;
                var fixedMods = new List<Modification> { PotassiumAdducts[1] }; // C
                var variableMods = new List<Modification> { SodiumAdducts[1] }; // C
                var precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                    .ToList();

                var fullSequences = precursors.Select(p => p.FullSequence).ToList();
                Assert.That(precursors.Count, Is.EqualTo(2));
                Assert.That(precursors.Any(p => p.NumFixedMods == 1));
                Assert.That(precursors.Any(p => p.NumVariableMods == 1));
                Assert.That(precursors.Any(p => p.NumFixedMods == 0));
                Assert.That(precursors.Any(p => p.NumVariableMods == 0));
                Assert.That(precursors.All(p => p.NumMods == 1));
                Assert.That(fullSequences.Contains("GUAC[Metal:Potassium on C]UG"));
                Assert.That(fullSequences.Contains("GUAC[Metal:Sodium on C]UG"));
            }
        }

        /// <summary>
        /// Test when two variable mods are used and share a localization
        /// expect three results, one unmodified, and two singly modified
        /// </summary>
        [Test]
        public static void TestVariableMods_LocalizationOverlap()
        {
            var rna = new RNA("GUACUG");
            var rnaDigestionParams = new RnaDigestionParams();

            for (int i = 1; i < 3; i++)
            {
                rnaDigestionParams.MaxMods = i;
                var fixedMods = new List<Modification> { }; // C
                var variableMods = new List<Modification> { PotassiumAdducts[1], SodiumAdducts[1] };
                var precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                    .ToList();

                // expect three results, one unmodified, and two singly modified
                var fullSequences = precursors.Select(p => p.FullSequence).ToList();
                Assert.That(precursors.Count, Is.EqualTo(3));
                Assert.That(precursors.Any(p => p.NumFixedMods == 0));
                Assert.That(precursors.Any(p => p.NumVariableMods == 1));
                Assert.That(fullSequences.Contains("GUACUG"));
                Assert.That(fullSequences.Contains("GUAC[Metal:Potassium on C]UG"));
                Assert.That(fullSequences.Contains("GUAC[Metal:Sodium on C]UG"));
            }
        }

        /// <summary>
        /// Test when one modification is annotated in the database, out of bounds
        /// expect two results, one with the fixed, and one with the variable
        /// </summary>
        [Test]
        public static void TestDatabaseAnnotatedMods_OutOfBounds()
        {
            var rnaDigestionParams = new RnaDigestionParams();
            var oneBasedModifications = new Dictionary<int, List<Modification>>()
            {
                { 23, new List<Modification>() { PotassiumAdducts[1] } }
            };
            var rna = new RNA("GUACUG", oneBasedModifications);
            
            for (int i = 1; i < 3; i++)
            {
                rnaDigestionParams.MaxMods = i;
                var fixedMods = new List<Modification> {  }; // C
                var variableMods = new List<Modification> { }; // C
                var precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                    .ToList();

                var fullSequences = precursors.Select(p => p.FullSequence).ToList();
                Assert.That(precursors.Count, Is.EqualTo(1));
                Assert.That(precursors.All(p => p.NumFixedMods == 0));
                Assert.That(precursors.All(p => p.NumVariableMods == 0));
                Assert.That(precursors.All(p => p.NumMods == 0));
                Assert.That(fullSequences.Contains("GUACUG"));
            }
        }

        /// <summary>
        /// Test when one modification is annotated in the database
        /// expect two results, one unmodified, and one singly modified
        /// </summary>
        [Test]
        public static void TestDatabaseAnnotatedMods_SingleModification()
        {
            var rnaDigestionParams = new RnaDigestionParams();
            var oneBasedModifications = new Dictionary<int, List<Modification>>()
            {
                { 4, new List<Modification>() { PotassiumAdducts[1] } }
            };
            var rna = new RNA("GUACUG", oneBasedPossibleLocalizedModifications: oneBasedModifications);

            
            for (int i = 1; i < 3; i++)
            {
                rnaDigestionParams.MaxMods = i;
                var fixedMods = new List<Modification> { }; // C
                var variableMods = new List<Modification> { };
                var precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                    .ToList();

                
                var fullSequences = precursors.Select(p => p.FullSequence).ToList();
                Assert.That(precursors.Count, Is.EqualTo(2));
                Assert.That(precursors[0].NumMods, Is.EqualTo(0));
                Assert.That(precursors[1].NumMods, Is.EqualTo(1));
                Assert.That(precursors[1].NumVariableMods, Is.EqualTo(1));
                Assert.That(fullSequences.Contains("GUACUG"));
                Assert.That(fullSequences.Contains("GUAC[Metal:Potassium on C]UG"));
            }
        }

        /// <summary>
        /// Test when two modifications are annotated in the database at the same location
        /// expect three results, one unmodified, and two singly modified
        /// </summary>
        [Test]
        public static void TestDatabaseAnnotatedMods_LocalizationOverlap()
        {
            var rnaDigestionParams = new RnaDigestionParams();
            var oneBasedModifications = new Dictionary<int, List<Modification>>()
            {
                { 4, new List<Modification>() { PotassiumAdducts[1], SodiumAdducts[1] } }
            };
            var rna = new RNA("GUACUG", oneBasedPossibleLocalizedModifications: oneBasedModifications);

            for (int i = 1; i < 3; i++)
            {
                rnaDigestionParams.MaxMods = i;
                var fixedMods = new List<Modification> { };
                var variableMods = new List<Modification> { };
                var precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                    .ToList();

                var fullSequences = precursors.Select(p => p.FullSequence).ToList();
                Assert.That(precursors.Count, Is.EqualTo(3));
                Assert.That(precursors.Any(p => p.NumFixedMods == 0));
                Assert.That(precursors.Any(p => p.NumVariableMods == 1));
                Assert.That(fullSequences.Contains("GUACUG"));
                Assert.That(fullSequences.Contains("GUAC[Metal:Potassium on C]UG"));
                Assert.That(fullSequences.Contains("GUAC[Metal:Sodium on C]UG"));
            }
        }

        /// <summary>
        /// Test when two terminal modifications are annotated in the database
        /// MaxMods 1: expect three results, one unmodified, and two singly modified
        /// MaxMods 2: expect four results, one unmodified, and two singly modified, and one double modified
        /// </summary>
        [Test]
        public static void TestDatabaseAnnotatedMods_TerminalMods()
        {
            var rnaDigestionParams = new RnaDigestionParams();
            var oneBasedModifications = new Dictionary<int, List<Modification>>()
            {
                { 1, new List<Modification>() { TerminalPotassiumAdducts[2]} },
                { 6, new List<Modification>() { TerminalSodiumAdducts[2]} }
            };
            var rna = new RNA("GUACUG", oneBasedPossibleLocalizedModifications: oneBasedModifications);

            // Test when two terminal modifications are annotated in the database
            for (int i = 1; i < 3; i++)
            {
                rnaDigestionParams.MaxMods = i;
                var fixedMods = new List<Modification> { };
                var variableMods = new List<Modification> { };
                var precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                    .ToList();

                // expect three results, one unmodified, and two singly modified
                var fullSequences = precursors.Select(p => p.FullSequence).ToList();
                Assert.That(precursors.Count, Is.EqualTo(2 + i));
                Assert.That(precursors.Any(p => p.NumFixedMods == 0));
                Assert.That(precursors.Any(p => p.NumVariableMods == 1));
                Assert.That(fullSequences.Contains("GUACUG"));
                Assert.That(fullSequences.Contains("[Metal:Potassium on G]GUACUG"));
                Assert.That(fullSequences.Contains("GUACUG-[Metal:Sodium on G]"));

                if (rnaDigestionParams.MaxMods != 2) continue;
                Assert.That(precursors.Any(p => p.NumVariableMods == 2));
                Assert.That(fullSequences.Contains("[Metal:Potassium on G]GUACUG-[Metal:Sodium on G]"));
            }
        }

        /// <summary>
        /// Test when two terminal modifications are annotated in the database and one database mod on first residue
        /// MaxMods 1: expect four results, one unmodified, and three singly modified
        /// MaxMods 2: expect seven results, one unmodified, and three singly modified, and three double modified
        /// MaxMods 3: expect eight results, one unmodified, and three singly modified, and three double modified, and one triply modified
        /// </summary>
        [Test]
        public static void TestDatabaseAnnotatedMods_TerminalMods_WithFirstResidueDatabaseMod()
        {
            var rnaDigestionParams = new RnaDigestionParams();
            var oneBasedModifications = new Dictionary<int, List<Modification>>()
            {
                { 1, new List<Modification>() { TerminalPotassiumAdducts[2], PotassiumAdducts[2] } },
                { 6, new List<Modification>() { TerminalSodiumAdducts[2]} }
            };
            var rna = new RNA("GUACUG", oneBasedPossibleLocalizedModifications: oneBasedModifications);

            // Test when two terminal modifications are annotated in the database and one database mod on first residue
            for (int i = 1; i < 4; i++)
            {
                rnaDigestionParams.MaxMods = i;
                var fixedMods = new List<Modification> { };
                var variableMods = new List<Modification> { };

                var precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                    .ToList();

                var fullSequences = precursors.Select(p => p.FullSequence).ToList();
                Assert.That(precursors.All(p => p.NumFixedMods == 0));

                switch (rnaDigestionParams.MaxMods)
                {
                    case 1:
                        Assert.That(precursors.Count(), Is.EqualTo(4));
                        Assert.That(precursors.Skip(1).All(p => p.NumVariableMods == 1));
                        break;
                    case 2:
                        Assert.That(precursors.Count(), Is.EqualTo(7));
                        Assert.That(precursors.Skip(1).All(p => p.NumVariableMods >= 1));
                        break;

                    case 3:
                        Assert.That(precursors.Count(), Is.EqualTo(8));
                        Assert.That(precursors.Skip(1).All(p => p.NumVariableMods >= 1));
                        break;
                }

                if (rnaDigestionParams.MaxMods >= 1)
                {
                    
                    Assert.That(fullSequences.Contains("GUACUG"));
                    Assert.That(fullSequences.Contains("[Metal:Potassium on G]GUACUG"));
                    Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UACUG"));
                    Assert.That(fullSequences.Contains("GUACUG-[Metal:Sodium on G]"));
                }
                else if (rnaDigestionParams.MaxMods >= 2)
                {
                    Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UACUG-[Metal:Sodium on G]"));
                    Assert.That(fullSequences.Contains("[Metal:Potassium on G]GUACUG-[Metal:Sodium on G]"));
                    Assert.That(fullSequences.Contains("[Metal:Potassium on G]G[Metal:Potassium on G]UACUG"));
                }
                else if (rnaDigestionParams.MaxMods >= 3)
                {
                    Assert.That(fullSequences.Contains("[Metal:Potassium on G]G[Metal:Potassium on G]UACUG[Metal:Sodium on G]"));
                }
            }
        }

        /// <summary>
        /// Test when two terminal modifications are annotated in the database and one database mod on first residue
        /// MaxMods 1: expect five results, one unmodified, and four singly modified
        /// MaxMods 2: expect eleven results, one unmodified, and four singly modified, and six double modified
        /// MaxMods 3: expect fifteen results, one unmodified, and four singly modified, and six double modified, and four triply modified
        /// </summary>
        [Test]
        public static void TestDatabaseAnnotatedMods_TerminalMods_WithFirstResidueVariableMod()
        {
            var rnaDigestionParams = new RnaDigestionParams();
            var oneBasedModifications = new Dictionary<int, List<Modification>>()
            {
                { 1, new List<Modification>() { TerminalPotassiumAdducts[2] } },
                { 6, new List<Modification>() { TerminalSodiumAdducts[2]} }
            };
            var rna = new RNA("GUACUG", oneBasedPossibleLocalizedModifications: oneBasedModifications);

            // Test when two terminal modifications are annotated in the database and one database mod on first residue
            for (int i = 1; i < 4; i++)
            {
                rnaDigestionParams.MaxMods = i;
                var fixedMods = new List<Modification> { };
                var variableMods = new List<Modification> { PotassiumAdducts[2] };

                var precursors = rna.Digest(rnaDigestionParams, fixedMods, variableMods)
                    .ToList();

                var fullSequences = precursors.Select(p => p.FullSequence).ToList();
                Assert.That(precursors.All(p => p.NumFixedMods == 0));

                switch (rnaDigestionParams.MaxMods)
                {
                    case 1:
                        Assert.That(precursors.Count(), Is.EqualTo(5));
                        Assert.That(precursors.Skip(1).All(p => p.NumVariableMods == 1));
                        break;
                    case 2:
                        Assert.That(precursors.Count(), Is.EqualTo(11));
                        Assert.That(precursors.Skip(1).All(p => p.NumVariableMods >= 1));
                        break;

                    case 3:
                        Assert.That(precursors.Count(), Is.EqualTo(15));
                        Assert.That(precursors.Skip(1).All(p => p.NumVariableMods >= 1));
                        break;
                }

                if (rnaDigestionParams.MaxMods >= 1)
                {

                    Assert.That(fullSequences.Contains("GUACUG"));
                    Assert.That(fullSequences.Contains("GUACUG[Metal:Potassium on G]"));
                    Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UACUG"));
                    Assert.That(fullSequences.Contains("GUACUG-[Metal:Sodium on G]"));
                    Assert.That(fullSequences.Contains("[Metal:Potassium on G]GUACUG"));
                }
                else if (rnaDigestionParams.MaxMods >= 2)
                {
                    Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UACUG[Metal:Potassium on G]"));
                    Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UACUG-[Metal:Sodium on G]"));
                    Assert.That(fullSequences.Contains("GUACUG[Metal:Potassium on G]-[Metal:Sodium on G]"));
                    Assert.That(fullSequences.Contains("[Metal:Potassium on G]GUACUG-[Metal:Sodium on G]"));
                    Assert.That(fullSequences.Contains("[Metal:Potassium on G]G[Metal:Potassium on G]UACUG"));
                    Assert.That(fullSequences.Contains("[Metal:Potassium on G]GUACUG[Metal:Potassium on G]"));
                }
                else if (rnaDigestionParams.MaxMods >= 3)
                {
                    Assert.That(fullSequences.Contains("[Metal:Potassium on G]G[Metal:Potassium on G]UACUG-[Metal:Sodium on G]"));
                    Assert.That(fullSequences.Contains("[Metal:Potassium on G]G[Metal:Potassium on G]UACUG[Metal:Potassium on G]"));

                    Assert.That(fullSequences.Contains("G[Metal:Potassium on G]UACUG[Metal:Potassium on G]-[Metal:Sodium on G]"));
                    Assert.That(fullSequences.Contains("[Metal:Potassium on G]GUACUG[Metal:Potassium on G]-[Metal:Sodium on G]"));
                }
            }
        }

        [Test]
        public static void TestDigestionMaxIsoforms()
        {
            var rna = new RNA("GUACUAGACUACAUGGUACAUCA");
            var rnaDigestionParams = new RnaDigestionParams();
            var variableMods = SodiumAdducts.Concat(PotassiumAdducts)
                .Concat(TerminalPotassiumAdducts).Concat(TerminalSodiumAdducts).ToList();

            var digestionProducts = rna.Digest(rnaDigestionParams, new List<Modification>(), variableMods)
                .ToList();
            Assert.That(digestionProducts.Count, Is.EqualTo(rnaDigestionParams.MaxModificationIsoforms));
        }

        #endregion
    }
}
