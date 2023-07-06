using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using NUnit.Framework;
using Transcriptomics;
using Transcriptomics.Fragmentation;
using UsefulProteomicsDatabases;

namespace Test.Transcriptomics
{
    /// <summary>
    /// Test Data generated with  http://rna.rega.kuleuven.be/masspec/mongo.htm
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestNucleicAcid
    {

        internal record SixmerTestCase(string Sequence, FragmentType Type, double[] NeutralMasses, string[] ChemicalFormulas);

        internal static IEnumerable<SixmerTestCase> GetSixmerTestCases()
        {
            Loaders.LoadElements();

            yield return new SixmerTestCase("GUACUG", FragmentType.a,
                new[] { 267.089, 573.114, 902.167, 1207.208, 1513.233 }, 
                new[] { "C10H13N5O4", "C19H24N7O12P", "C29H36N12O18P2", "C38H48N15O25P3", "C47H59N17O33P4" });
            yield return new SixmerTestCase("GUACUG", FragmentType.b, 
            new[] { 283.084, 589.109, 918.162, 1223.203, 1529.228 },
            new[] { "C10H13N5O5", "C19H24N7O13P", "C29H36N12O19P2", "C38H48N15O26P3", "C47H59N17O34P4" });
            yield return new SixmerTestCase("GUACUG", FragmentType.c,
                new[] { 347.055, 653.081, 982.133, 1287.174, 1593.2 },
            new[] { "C10H14N5O7P", "C19H25N7O15P2", "C29H37N12O21P3", "C38H49N15O28P4", "C47H60N17O36P5", });
            yield return new SixmerTestCase("GUACUG", FragmentType.d,
                new[] { 363.05, 669.075, 998.128, 1303.169, 1609.195 },
            new[] { "C10H14N5O8P", "C19H25N7O16P2", "C29H37N12O22P3", "C38H49N15O29P4", "C47H60N17O37P5", });
            yield return new SixmerTestCase("GUACUG", FragmentType.dH2O,
                new[] { 345.039, 651.064, 980.116, 1285.157, 1591.184 },
            new[] { "C10H12N5O7P", "C19H23N7O15P2", "C29H35N12O21P3", "C38H47N15O28P4", "C47H58N17O36P5", });
            yield return new SixmerTestCase("GUACUG", FragmentType.w, 
                new[] { 363.049, 669.074, 974.115, 1303.169, 1609.195 },
            new[] { "C10H14N5O8P", "C19H25N7O16P2", "C28H37N10O23P3", "C38H49N15O29P4", "C47H60N17O37P5", });
            yield return new SixmerTestCase("GUACUG", FragmentType.x,
                new[] { 347.055, 653.081, 958.122, 1287.174, 1593.2 },
            new[] { "C10H14N5O7P", "C19H25N7O15P2", "C28H37N10O22P3", "C38H49N15O28P4", "C47H60N17O36P5" });
            yield return new SixmerTestCase("GUACUG", FragmentType.y,
                new[] { 283.084, 589.109, 894.15, 1223.203, 1529.228 },
            new[] { "C10H13N5O5", "C19H24N7O13P", "C28H36N10O20P2", "C38H48N15O26P3", "C47H59N17O34P4", });
            yield return new SixmerTestCase("GUACUG", FragmentType.z,
                new[] { 267.089, 573.124, 878.156, 1207.208, 1513.233 },
            new[] { "C10H13N5O4", "C19H24N7O12P", "C28H36N10O19P2", "C38H48N15O25P3", "C47H59N17O33P4", });


            yield return new SixmerTestCase("GUACUG", FragmentType.aBase,
                new[] { 114.03, 459.07, 765.095, 1094.147, 1399.198 },
                new[] { "C5H6O3", "C15H18N5O10P", "C24H29N7O18P2", "C34H41N12O24P3", "C43H53N15O31P4" });
            yield return new SixmerTestCase("GUACUG", FragmentType.bBase,
                new[] { 130.027, 475.074, 781.099, 1110.152, 1415.193 }, 
                new[] { "C5H6O4", "C15H18N5O11P", "C24H29N7O19P2", "C34H41N12O25P3", "C43H53N15O32P4" });
            yield return new SixmerTestCase("GUACUG", FragmentType.cBase, 
                new[] { 193.998, 539.045, 845.071, 1174.123, 1479.164 },
                new[] { "C5H7O6P", "C15H19N5O13P2", "C24H30N7O21P3", "C34H42N12O27P4", "C43H54N15O34P5" });
            yield return new SixmerTestCase("GUACUG", FragmentType.dBase,
                new[] { 209.993, 555.04, 861.066, 1190.118, 1495.16 },
                new[] { "C5H7O7P", "C15H19N5O14P2", "C24H30N7O22P3", "C34H42N12O28P4", "C43H54N15O35P5" });



            //yield return new SixmerTestCase("GUACUG", FragmentType.wBase, new[] { 459.078, 765.103, 1094.155, 1399.196 },
            //new[] { "C10H13N5O4", "C19H24N7O12P", "C28H36N10O19P2", "C38H48N15O25P3", "C47H59N17O33P4", });
            //yield return new SixmerTestCase("GUACUG", FragmentType.xBase, new[] { 459.078, 765.103, 1094.155, 1399.196 },
            //new[] { "C10H13N5O4", "C19H24N7O12P", "C28H36N10O19P2", "C38H48N15O25P3", "C47H59N17O33P4", });
            //yield return new SixmerTestCase("GUACUG", FragmentType.yBase, new[] { 459.078, 765.103, 1094.155, 1399.196 },
            //new[] { "C10H13N5O4", "C19H24N7O12P", "C28H36N10O19P2", "C38H48N15O25P3", "C47H59N17O33P4", });
            //yield return new SixmerTestCase("GUACUG", FragmentType.zBase, new[] { 459.078, 765.103, 1094.155, 1399.196 },
            //new[] { "C10H13N5O4", "C19H24N7O12P", "C28H36N10O19P2", "C38H48N15O25P3", "C47H59N17O33P4", });

            // TODO: Add dot
        }


        [Test]
        [TestCase("GUACUG", 1874.281)]
        [TestCase("A", 267.096)]
        [TestCase("C", 243.085)]
        [TestCase("U", 244.069)]
        [TestCase("G", 283.091)]
        [TestCase("GU", 589.116)]
        [TestCase("AAA", 925.200)]
        [TestCase("CCC", 853.166)]
        [TestCase("UUU", 856.119)]
        [TestCase("GGG", 973.185)]
        public void TestConstructorsAndEquality(string sequence, double monoMass)
        {
            // test constructors and equality
            RNA rna = new RNA(sequence);

            Assert.IsNotNull(rna);
            Assert.That(rna.Length, Is.EqualTo(sequence.Length));
            Assert.That(rna.MonoisotopicMass, Is.EqualTo(monoMass).Within(0.01));
            Assert.That(rna.FivePrimeTerminus.Equals(NucleicAcid.DefaultFivePrimeTerminus));
            Assert.That(rna.ThreePrimeTerminus.Equals(NucleicAcid.DefaultThreePrimeTerminus));
            List<Nucleotide> nucList = new();
            foreach (var nucleotide in sequence)
            {
                nucList.Add(Nucleotide.GetResidue(nucleotide));
            }
            Assert.That(rna.NucleicAcidArray.SequenceEqual(nucList.ToArray()));

            var rna2 = new RNA(sequence, NucleicAcid.DefaultFivePrimeTerminus, NucleicAcid.DefaultThreePrimeTerminus);

            Assert.IsNotNull(rna2);
            Assert.That(rna2.Length, Is.EqualTo(sequence.Length));
            Assert.That(rna2.MonoisotopicMass, Is.EqualTo(monoMass).Within(0.01));
            Assert.That(rna.FivePrimeTerminus.Equals(NucleicAcid.DefaultFivePrimeTerminus));
            Assert.That(rna.ThreePrimeTerminus.Equals(NucleicAcid.DefaultThreePrimeTerminus));
            nucList.Clear();
            foreach (var nucleotide in sequence)
            {
                nucList.Add(Nucleotide.GetResidue(nucleotide));
            }
            Assert.That(rna.NucleicAcidArray.SequenceEqual(nucList.ToArray()));

            Assert.That(rna.Equals(rna2));
            Assert.That(rna.Equals((object)rna2));
        }


        [Test]
        [TestCase("GUACUG", new[] {-1, -2, -3, -4, -5}, new[] {1873.273, 936.133, 623.752, 467.562, 373.848})]
        public void TestElectroSpraySeries(string sequence, int[] charges, double[] mzs)
        {
            RNA rna = new(sequence);

            int i = 0;
            foreach (var ion in rna.GetElectrospraySeries(charges.First(), charges.Last()))
            {
                Assert.That(ion, Is.EqualTo(mzs[i]).Within(0.001));
                i++;
            }
        }

        [Test]
        [TestCaseSource(nameof(GetSixmerTestCases))]
        public void TestGetNeutralFragments(SixmerTestCase testCase)
        {
            RNA rna = new(testCase.Sequence);

            var neutralFragments = rna.GetNeutralFragments(testCase.Type).ToList();
            //Assert.That(neutralFragments.Count, Is.EqualTo(neutralMasses.Length));
            for (int i = 1; i < neutralFragments.Count; i++)
            {
                Assert.That(neutralFragments[i].MonoisotopicMass, Is.EqualTo(testCase.NeutralMasses[i]).Within(0.01));
            }
        }


        [Test]
        [TestCaseSource(nameof(GetSixmerTestCases))]
        public void TestGetNeutralFragmentEquations(SixmerTestCase testCase)
        {
            RNA rna = new(testCase.Sequence);
            var neutralFragments = rna.GetNeutralFragments(testCase.Type, true).Select(p => (ChemicalFormulaFragment)p).ToList();
            for (int i = 0; i < neutralFragments.Count; i++)
            {
                var theoreticalFormula = ChemicalFormula.ParseFormula(testCase.ChemicalFormulas[i]);
                var calculatedFormula = neutralFragments[i].ThisChemicalFormula;

                Assert.That(calculatedFormula, Is.EqualTo(theoreticalFormula));
                Assert.That(neutralFragments[i].MonoisotopicMass,
                    Is.EqualTo(neutralFragments[i].ThisChemicalFormula.MonoisotopicMass));

            }
        }


        [Test]
        public static void TESTNAME()
        {
            var rna = new RNA("GUACUG");
            var a = rna.GetNeutralFragments(FragmentType.a).ToList();
            var adot = rna.GetNeutralFragments(FragmentType.adot).ToList();
            var abase = rna.GetNeutralFragments(FragmentType.aBase).ToList();
        }
    }
}
