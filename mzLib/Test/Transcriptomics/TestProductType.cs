using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Transcriptomics;

namespace Test.Transcriptomics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestProductType
    {
        [Test]
        [TestCase(DissociationType.HCD, new[] { ProductType.y, ProductType.w, ProductType.aBase, ProductType.dWaterLoss })]
        [TestCase(DissociationType.CID, new[] { ProductType.y, ProductType.w, ProductType.aBase, ProductType.dWaterLoss })]
        public void TestProductTypes_Dissociation(DissociationType dissociation, ProductType[] products)
        {
            CollectionAssert.AreEquivalent(products, dissociation.GetRnaProductTypesFromDissociationType());
        }

        [Test]
        [TestCase(FragmentationTerminus.FivePrime, new[]
        {
            ProductType.a, ProductType.aWaterLoss, ProductType.aBase,
            ProductType.b, ProductType.bWaterLoss, ProductType.bBase,
            ProductType.c, ProductType.cWaterLoss, ProductType.cBase,
            ProductType.d, ProductType.dWaterLoss, ProductType.dBase,
        })]
        [TestCase(FragmentationTerminus.ThreePrime, new[]
        {
            ProductType.w, ProductType.wWaterLoss, ProductType.wBase,
            ProductType.x, ProductType.xWaterLoss, ProductType.xBase,
            ProductType.y, ProductType.yWaterLoss, ProductType.yBase,
            ProductType.z, ProductType.zWaterLoss, ProductType.zBase,
        })]
        public void TestProductTypes_Terminus(FragmentationTerminus terminus, ProductType[] products)
        {
            CollectionAssert.AreEquivalent(products, terminus.GetRnaTerminusSpecificProductTypes());
        }

        [Test]
        [TestCase(DissociationType.HCD, FragmentationTerminus.FivePrime, new[] { ProductType.aBase, ProductType.dWaterLoss })]
        [TestCase(DissociationType.HCD, FragmentationTerminus.ThreePrime, new[] { ProductType.w, ProductType.y })]
        [TestCase(DissociationType.HCD, FragmentationTerminus.Both, new[] { ProductType.w, ProductType.y, ProductType.aBase, ProductType.dWaterLoss })]
        [TestCase(DissociationType.CID, FragmentationTerminus.FivePrime, new[] { ProductType.aBase, ProductType.dWaterLoss })]
        [TestCase(DissociationType.CID, FragmentationTerminus.ThreePrime, new[] { ProductType.w, ProductType.y })]
        [TestCase(DissociationType.CID, FragmentationTerminus.Both, new[] { ProductType.w, ProductType.y, ProductType.aBase, ProductType.dWaterLoss })]
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
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-1H").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.y:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-3P-1").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.zWaterLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-5H-2P-1").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.aWaterLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H-1O-1").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.aBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H-2").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.bBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O1H-2").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.cWaterLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O2P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.cBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O3H-1P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.d:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O4H2P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.dWaterLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O3P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.dBase:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O4H-1P").MonoisotopicMass, 2).Value, mass);
                        break;


                    case ProductType.w:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.wWaterLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H-1O-1").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.xWaterLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-2H-1").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.yWaterLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-4H-2P-1").MonoisotopicMass, 2).Value, mass);
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
                    case ProductType.aWaterLoss:
                    case ProductType.aBase:
                    case ProductType.b:
                    case ProductType.bWaterLoss:
                    case ProductType.bBase:
                    case ProductType.c:
                    case ProductType.cWaterLoss:
                    case ProductType.cBase:
                    case ProductType.d:
                    case ProductType.dWaterLoss:
                    case ProductType.dBase:
                        Assert.That(type.GetRnaTerminusType(), Is.EqualTo(FragmentationTerminus.FivePrime));
                        break;

                    case ProductType.w:
                    case ProductType.wWaterLoss:
                    case ProductType.wBase:
                    case ProductType.x:
                    case ProductType.xWaterLoss:
                    case ProductType.xBase:
                    case ProductType.y:
                    case ProductType.yWaterLoss:
                    case ProductType.yBase:
                    case ProductType.z:
                    case ProductType.zWaterLoss:
                    case ProductType.zBase:
                        Assert.That(type.GetRnaTerminusType(), Is.EqualTo(FragmentationTerminus.ThreePrime));
                        break;

                    case ProductType.M:
                    case ProductType.aStar:
                    case ProductType.bAmmoniaLoss:
                    case ProductType.D:
                    case ProductType.Ycore:
                    case ProductType.Y:
                    case ProductType.aDegree:
                    case ProductType.yAmmoniaLoss:
                    case ProductType.zPlusOne:
                    case ProductType.zDot:
                        Assert.Throws<ArgumentOutOfRangeException>(() => type.GetRnaTerminusType());
                        break;
                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }
        }

        [Test]
        [TestCase(ProductType.a, ProductType.aWaterLoss)]
        [TestCase(ProductType.b, ProductType.bWaterLoss)]
        [TestCase(ProductType.c, ProductType.cWaterLoss)]
        [TestCase(ProductType.d, ProductType.dWaterLoss)]
        [TestCase(ProductType.w, ProductType.wWaterLoss)]
        [TestCase(ProductType.x, ProductType.xWaterLoss)]
        [TestCase(ProductType.y, ProductType.yWaterLoss)]
        [TestCase(ProductType.z, ProductType.zWaterLoss)]
        public void EnsureWaterLossMassesAreCorrect(ProductType normal, ProductType waterLoss)
        {
            RNA rna = new("GUACUG");

            List<IProduct> normalFragments = rna.GetNeutralFragments(normal).ToList();
            List<IProduct> waterLossFragments = rna.GetNeutralFragments(waterLoss).ToList();
            for (var index = 0; index < waterLossFragments.Count; index++)
            {
                var waterLossFragment = waterLossFragments[index];
                var normalFragment = normalFragments[index];
                var watermass = 2 * Constants.ProtonMass + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

                Assert.That(normalFragment.MonoisotopicMass, Is.EqualTo(waterLossFragment.MonoisotopicMass + watermass).Within(0.01));
            }
        }
    }
}
