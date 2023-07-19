using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using MassSpectrometry;
using Microsoft.VisualBasic;
using NUnit.Framework;
using Proteomics.Fragmentation;
using Readers;
using Transcriptomics;
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
        internal record SixmerTestCase(string Sequence, ProductType Type, double[] NeutralMasses, string[] ChemicalFormulas);

        internal static IEnumerable<SixmerTestCase> GetSixmerIndividualFragmentTypeTestCases()
        {
            Loaders.LoadElements();

            yield return new SixmerTestCase("GUACUG", ProductType.a,
                new[] { 267.089, 573.114, 902.167, 1207.208, 1513.233 }, 
                new[] { "C10H13N5O4", "C19H24N7O12P", "C29H36N12O18P2", "C38H48N15O25P3", "C47H59N17O33P4" });
            yield return new SixmerTestCase("GUACUG", ProductType.b, 
            new[] { 283.084, 589.109, 918.162, 1223.203, 1529.228 },
            new[] { "C10H13N5O5", "C19H24N7O13P", "C29H36N12O19P2", "C38H48N15O26P3", "C47H59N17O34P4" });
            yield return new SixmerTestCase("GUACUG", ProductType.c,
                new[] { 347.055, 653.081, 982.133, 1287.174, 1593.2 },
            new[] { "C10H14N5O7P", "C19H25N7O15P2", "C29H37N12O21P3", "C38H49N15O28P4", "C47H60N17O36P5", });
            yield return new SixmerTestCase("GUACUG", ProductType.d,
                new[] { 363.05, 669.075, 998.128, 1303.169, 1609.195 },
            new[] { "C10H14N5O8P", "C19H25N7O16P2", "C29H37N12O22P3", "C38H49N15O29P4", "C47H60N17O37P5", });
            yield return new SixmerTestCase("GUACUG", ProductType.dH2O,
                new[] { 345.039, 651.064, 980.116, 1285.157, 1591.184 },
            new[] { "C10H12N5O7P", "C19H23N7O15P2", "C29H35N12O21P3", "C38H47N15O28P4", "C47H58N17O36P5", });
            yield return new SixmerTestCase("GUACUG", ProductType.w, 
                new[] { 363.049, 669.074, 974.115, 1303.169, 1609.195 },
            new[] { "C10H14N5O8P", "C19H25N7O16P2", "C28H37N10O23P3", "C38H49N15O29P4", "C47H60N17O37P5", });
            yield return new SixmerTestCase("GUACUG", ProductType.x,
                new[] { 347.055, 653.081, 958.122, 1287.174, 1593.2 },
            new[] { "C10H14N5O7P", "C19H25N7O15P2", "C28H37N10O22P3", "C38H49N15O28P4", "C47H60N17O36P5" });
            yield return new SixmerTestCase("GUACUG", ProductType.y,
                new[] { 283.084, 589.109, 894.15, 1223.203, 1529.228 },
            new[] { "C10H13N5O5", "C19H24N7O13P", "C28H36N10O20P2", "C38H48N15O26P3", "C47H59N17O34P4", });
            yield return new SixmerTestCase("GUACUG", ProductType.z,
                new[] { 267.089, 573.124, 878.156, 1207.208, 1513.233 },
            new[] { "C10H13N5O4", "C19H24N7O12P", "C28H36N10O19P2", "C38H48N15O25P3", "C47H59N17O33P4", });


            yield return new SixmerTestCase("GUACUG", ProductType.aBase,
                new[] { 114.03, 459.07, 765.095, 1094.147, 1399.198 },
                new[] { "C5H6O3", "C15H18N5O10P", "C24H29N7O18P2", "C34H41N12O24P3", "C43H53N15O31P4" });
            yield return new SixmerTestCase("GUACUG", ProductType.bBase,
                new[] { 130.027, 475.074, 781.099, 1110.152, 1415.193 }, 
                new[] { "C5H6O4", "C15H18N5O11P", "C24H29N7O19P2", "C34H41N12O25P3", "C43H53N15O32P4" });
            yield return new SixmerTestCase("GUACUG", ProductType.cBase, 
                new[] { 193.998, 539.045, 845.071, 1174.123, 1479.164 },
                new[] { "C5H7O6P", "C15H19N5O13P2", "C24H30N7O21P3", "C34H42N12O27P4", "C43H54N15O34P5" });
            yield return new SixmerTestCase("GUACUG", ProductType.dBase,
                new[] { 209.993, 555.04, 861.066, 1190.118, 1495.16 },
                new[] { "C5H7O7P", "C15H19N5O14P2", "C24H30N7O22P3", "C34H42N12O28P4", "C43H54N15O35P5" });

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

        #region Product Type Retreival

        [Test]
        [TestCase(DissociationType.HCD, new[] {ProductType.y, ProductType.w, ProductType.aBase, ProductType.dH2O} )]
        [TestCase(DissociationType.CID, new[] {ProductType.y, ProductType.w, ProductType.aBase, ProductType.dH2O} )]
        public void TestProductTypes_Dissociation(DissociationType dissociation, ProductType[] products)
        {
            CollectionAssert.AreEquivalent(products, dissociation.GetRnaProductTypesFromDissociationType());
        }

        [Test]
        [TestCase(FragmentationTerminus.FivePrime, new[]
        {
            ProductType.a, ProductType.adot, ProductType.aBase,
            ProductType.b, ProductType.bdot, ProductType.bBase,
            ProductType.c, ProductType.cdot, ProductType.cBase,
            ProductType.d, ProductType.ddot, ProductType.dBase, ProductType.dH2O
        })]
        [TestCase(FragmentationTerminus.ThreePrime, new[]
        {
            ProductType.w, ProductType.wdot, ProductType.wBase,
            ProductType.x, ProductType.xdot, ProductType.xBase,
            ProductType.y, ProductType.ydot, ProductType.yBase,
            ProductType.z, ProductType.zDot, ProductType.zBase,
        })]
        public void TestProductTypes_Terminus(FragmentationTerminus terminus, ProductType[] products)
        {
            CollectionAssert.AreEquivalent(products, terminus.GetRnaTerminusSpecificProductTypes());
        }

        [Test]
        [TestCase(DissociationType.HCD, FragmentationTerminus.FivePrime, new[] {ProductType.aBase, ProductType.dH2O})]
        [TestCase(DissociationType.HCD, FragmentationTerminus.ThreePrime, new[] {ProductType.w, ProductType.y})]
        [TestCase(DissociationType.HCD, FragmentationTerminus.Both, new[] {ProductType.w, ProductType.y, ProductType.aBase, ProductType.dH2O })]
        [TestCase(DissociationType.CID, FragmentationTerminus.FivePrime, new[] {ProductType.aBase, ProductType.dH2O})]
        [TestCase(DissociationType.CID, FragmentationTerminus.ThreePrime, new[] {ProductType.w, ProductType.y})]
        [TestCase(DissociationType.CID, FragmentationTerminus.Both, new[] {ProductType.w, ProductType.y, ProductType.aBase, ProductType.dH2O })]
        public void TestProductTypes_TerminusAndDissociation(DissociationType dissociation, FragmentationTerminus terminus, ProductType[] products)
        {
            CollectionAssert.AreEquivalent(products, dissociation.GetRnaTerminusSpecificProductTypesFromDissociation(terminus));
        }

        [Test]
        public static void Test_NeutralMassShiftFromProductType()
        {
            foreach (ProductType p in Enum.GetValues(typeof(ProductType)))
            {
                double mass = 0;
                switch (p)
                {
                    case ProductType.a:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.b:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("OH").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.c:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O3H2P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.x:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-1H").MonoisotopicMass , 2).Value, mass);
                        break;

                    case ProductType.y:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-3P-1").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.zDot:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-4HP-1").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.adot:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H2").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.aBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H-2").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.bdot:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("OH2").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.bBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O1H-2").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.cdot:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O3H3P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.cBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O3H-1P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.d:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O4H2P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.ddot:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O4H3P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.dBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O4H-1P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.dH2O:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O3P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.w:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.wdot:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H2").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.xdot:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-1H2").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.ydot:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-3HP-1").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.z:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-4P-1").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.wBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H-2").MonoisotopicMass, 2).Value, mass);
                        break;
                    case ProductType.xBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-1H-2").MonoisotopicMass, 2).Value, mass);
                        break;
                    case ProductType.yBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-3H-2P-1").MonoisotopicMass, 2).Value, mass);
                        break;
                    case ProductType.zBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-4H-3P-1").MonoisotopicMass, 2).Value, mass);
                        break;
                }
            }
        }

        [Test]
        public void TestProductTypes_GetRnaTerminusType()
        {
            foreach (var type in Enum.GetValues<ProductType>())
            {
                switch (type)
                {
                    case ProductType.a:
                    case ProductType.adot:
                    case ProductType.aBase:
                    case ProductType.b:
                    case ProductType.bdot:
                    case ProductType.bBase:
                    case ProductType.c:
                    case ProductType.cdot:
                    case ProductType.cBase:
                    case ProductType.d:
                    case ProductType.ddot:
                    case ProductType.dBase:
                    case ProductType.dH2O:
                        Assert.That(type.GetRnaTerminusType(), Is.EqualTo(FragmentationTerminus.FivePrime));
                        break;

                    case ProductType.w:
                    case ProductType.wdot:
                    case ProductType.wBase:
                    case ProductType.x:
                    case ProductType.xdot:
                    case ProductType.xBase:
                    case ProductType.y:
                    case ProductType.ydot:
                    case ProductType.yBase:
                    case ProductType.z:
                    case ProductType.zDot:
                    case ProductType.zBase:
                        Assert.That(type.GetRnaTerminusType(), Is.EqualTo(FragmentationTerminus.ThreePrime));
                        break;

                    case ProductType.M:
                    case ProductType.aStar:
                    case ProductType.bAmmoniaLoss:
                    case ProductType.bWaterLoss:
                    case ProductType.D:
                    case ProductType.Ycore:
                    case ProductType.Y:
                    case ProductType.aDegree:
                    case ProductType.yAmmoniaLoss:
                    case ProductType.zPlusOne:
                    case ProductType.yWaterLoss:
                        Assert.Throws<ArgumentOutOfRangeException>(() => type.GetRnaTerminusType());
                        break;
                }
            }
        }

        #endregion

        #region Fragmentation

        [Test]
        [TestCaseSource(nameof(GetSixmerIndividualFragmentTypeTestCases))]
        public void TestGetNeutralFragments(SixmerTestCase testCase)
        {
            RNA rna = new(testCase.Sequence);

            var neutralFragments = rna.GetNeutralFragments(testCase.Type).ToList();
            for (int i = 1; i < neutralFragments.Count; i++)
            {
                Assert.That(neutralFragments[i].NeutralMass, Is.EqualTo(testCase.NeutralMasses[i]).Within(0.01));
            }
        }

        [Test]
        public void TestFragmentation_Unmodified()
        {
            List<IProduct> productsToTest = new List<IProduct>();
            var cidProductTypes = DissociationType.CID.GetRnaProductTypesFromDissociationType();
            foreach (var testCase in GetSixmerIndividualFragmentTypeTestCases().Where(p => cidProductTypes.Contains(p.Type)))
            {
                for (int i = 0; i < testCase.NeutralMasses.Length; i++)
                {
                    var terminus = testCase.Type.GetRnaTerminusType();
                    double neutralMass = testCase.NeutralMasses[i];
                    int fragmentNum = i;
                    double neutralLoss = 0;
                    if (testCase.Type.ToString().Contains("Base"))
                    {
                       // neutralLoss = previousNucleotide.BaseChemicalFormula.MonoisotopicMass;
                    }

                    var t = new RnaProduct(testCase.Type, terminus, neutralMass, i + 1, i + 1, neutralLoss);
                    productsToTest.Add(t);
                }
            }


            List<IProduct> products = new();
            RNA rna = new RNA("GUACUG");
            rna.Fragment(DissociationType.CID, FragmentationTerminus.Both, products);
            CollectionAssert.AreEquivalent(products, productsToTest);

        }

        #endregion
    }
}
