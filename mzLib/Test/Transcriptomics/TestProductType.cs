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
using Omics.Fragmentation;
using Omics.Fragmentation.Oligo;
using Omics.Modifications;
using Transcriptomics.Digestion;

namespace Test.Transcriptomics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestProductType
    {
        [Test]
        [TestCase(DissociationType.HCD, new[] { ProductType.M, ProductType.y, ProductType.w, ProductType.aBaseLoss, ProductType.dWaterLoss })]
        [TestCase(DissociationType.CID, new[] { ProductType.M, ProductType.y, ProductType.yWaterLoss, ProductType.c, ProductType.w, ProductType.aBaseLoss, ProductType.dWaterLoss })]
        public void TestProductTypes_Dissociation(DissociationType dissociation, ProductType[] products)
        {
            CollectionAssert.AreEquivalent(products, dissociation.GetRnaProductTypesFromDissociationType());
        }

        [Test]
        [TestCase(FragmentationTerminus.FivePrime, new[]
        {
            ProductType.a, ProductType.aWaterLoss, ProductType.aBaseLoss,
            ProductType.b, ProductType.bWaterLoss, ProductType.bBaseLoss,
            ProductType.c, ProductType.cWaterLoss, ProductType.cBaseLoss,
            ProductType.d, ProductType.dWaterLoss, ProductType.dBaseLoss,
        })]
        [TestCase(FragmentationTerminus.ThreePrime, new[]
        {
            ProductType.w, ProductType.wWaterLoss, ProductType.wBaseLoss,
            ProductType.x, ProductType.xWaterLoss, ProductType.xBaseLoss,
            ProductType.y, ProductType.yWaterLoss, ProductType.yBaseLoss,
            ProductType.z, ProductType.zWaterLoss, ProductType.zBaseLoss,
        })]
        public void TestProductTypes_Terminus(FragmentationTerminus terminus, ProductType[] products)
        {
            CollectionAssert.AreEquivalent(products, terminus.GetRnaTerminusSpecificProductTypes());
        }

        [Test]
        [TestCase(DissociationType.HCD, FragmentationTerminus.FivePrime, new[] { ProductType.aBaseLoss, ProductType.dWaterLoss })]
        [TestCase(DissociationType.HCD, FragmentationTerminus.ThreePrime, new[] { ProductType.w, ProductType.y })]
        [TestCase(DissociationType.HCD, FragmentationTerminus.Both, new[] { ProductType.M, ProductType.y, ProductType.w, ProductType.aBaseLoss, ProductType.dWaterLoss })]
        [TestCase(DissociationType.CID, FragmentationTerminus.FivePrime, new[] { ProductType.aBaseLoss, ProductType.dWaterLoss, ProductType.c })]
        [TestCase(DissociationType.CID, FragmentationTerminus.ThreePrime, new[] { ProductType.w, ProductType.y, ProductType.yWaterLoss })]
        [TestCase(DissociationType.CID, FragmentationTerminus.Both, new[] { ProductType.M, ProductType.y, ProductType.yWaterLoss, ProductType.c, ProductType.w, ProductType.aBaseLoss, ProductType.dWaterLoss })]
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

                    case ProductType.aBaseLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H-2").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.bBaseLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O1H-2").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.cWaterLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O2P").MonoisotopicMass, 2).Value, mass);
                        break;

                    case ProductType.cBaseLoss:
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

                    case ProductType.dBaseLoss:
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

                    case ProductType.wBaseLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("H-2").MonoisotopicMass, 2).Value, mass);
                        break;
                    case ProductType.xBaseLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-1H-2").MonoisotopicMass, 2).Value, mass);
                        break;
                    case ProductType.yBaseLoss:
                        mass = ClassExtensions.RoundedDouble(p.GetRnaMassShiftFromProductType(), 2).Value;
                        Assert.AreEqual(ClassExtensions.RoundedDouble(ChemicalFormula.ParseFormula("O-3H-2P-1").MonoisotopicMass, 2).Value, mass);
                        break;
                    case ProductType.zBaseLoss:
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
                    case ProductType.aBaseLoss:
                    case ProductType.b:
                    case ProductType.bWaterLoss:
                    case ProductType.bBaseLoss:
                    case ProductType.c:
                    case ProductType.cWaterLoss:
                    case ProductType.cBaseLoss:
                    case ProductType.d:
                    case ProductType.dWaterLoss:
                    case ProductType.dBaseLoss:
                        Assert.That(type.GetRnaTerminusType(), Is.EqualTo(FragmentationTerminus.FivePrime));
                        break;

                    case ProductType.w:
                    case ProductType.wWaterLoss:
                    case ProductType.wBaseLoss:
                    case ProductType.x:
                    case ProductType.xWaterLoss:
                    case ProductType.xBaseLoss:
                    case ProductType.y:
                    case ProductType.yWaterLoss:
                    case ProductType.yBaseLoss:
                    case ProductType.z:
                    case ProductType.zWaterLoss:
                    case ProductType.zBaseLoss:
                        Assert.That(type.GetRnaTerminusType(), Is.EqualTo(FragmentationTerminus.ThreePrime));
                        break;

                    case ProductType.M:
                        Assert.That(type.GetRnaTerminusType(), Is.EqualTo(FragmentationTerminus.None));
                        break;

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
            var rna = new RNA("GUACUG")
                .Digest(new RnaDigestionParams(), new List<Modification>(), new List<Modification>())
                .First() as OligoWithSetMods ?? throw new NullReferenceException();

            List<Product> normalFragments = rna.GetNeutralFragments(normal).ToList();
            List<Product> waterLossFragments = rna.GetNeutralFragments(waterLoss).ToList();
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
