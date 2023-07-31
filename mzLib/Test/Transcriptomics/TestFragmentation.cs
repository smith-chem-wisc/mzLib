using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Test.Transcriptomics.TestNucleicAcid;
using Transcriptomics;
using MassSpectrometry;
using Proteomics.Fragmentation;
using UsefulProteomicsDatabases.Generated;

namespace Test.Transcriptomics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestFragmentation
    {

        internal static IEnumerable<SixmerTestCase> GetSixMerIndividualFragmentTypeTestCases() =>
            TestNucleicAcid.GetSixmerIndividualFragmentTypeTestCases();

        [Test]
        [TestCaseSource(nameof(GetSixMerIndividualFragmentTypeTestCases))]
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
            RNA rna = new RNA("GUACUG");
            List<IProduct> products = new();

            // both termini
            rna.Fragment(DissociationType.CID, FragmentationTerminus.Both, products);
            Assert.That(products.Count, Is.EqualTo(20));
            Assert.That(products.All(p => p.ProductType is ProductType.w or ProductType.dWaterLoss or ProductType.y or ProductType.aBaseLoss));

            // only 5'
            rna.Fragment(DissociationType.CID, FragmentationTerminus.FivePrime, products);
            Assert.That(products.Count, Is.EqualTo(10));
            Assert.That(products.All(p => p.ProductType is ProductType.dWaterLoss or  ProductType.aBaseLoss));

            // only 3'
            rna.Fragment(DissociationType.CID, FragmentationTerminus.ThreePrime, products);
            Assert.That(products.Count, Is.EqualTo(10));
            Assert.That(products.All(p => p.ProductType is ProductType.w or ProductType.y ));
        }

        [Test]
        [TestCaseSource(nameof(GetSixMerIndividualFragmentTypeTestCases))]
        public void TestRnaFragments(SixmerTestCase testCase)
        {
            RNA rna = new("GUACUG");
            List<RnaProduct> products = rna.GetNeutralFragments(testCase.Type).Select(p => (RnaProduct)p).ToList();

            for (int i = 0; i < products.Count; i++)
            {
                var product = products[i];
                Assert.That(testCase.Type, Is.EqualTo(product.ProductType));
                Assert.That(testCase.Type.GetRnaTerminusType(), Is.EqualTo(product.Terminus));
                Assert.That(testCase.NeutralMasses[i], Is.EqualTo(product.NeutralMass).Within(0.01));
                Assert.That(testCase.NeutralMasses[i], Is.EqualTo(product.MonoisotopicMass).Within(0.01));
                Assert.That(0, Is.EqualTo(product.NeutralLoss));
                Assert.That(null, Is.EqualTo(product.SecondaryProductType));
                Assert.That(null, Is.EqualTo(product.Parent));
                Assert.That(0, Is.EqualTo(product.SecondaryFragmentNumber));

                string annotation = $"{product.ProductType}{product.FragmentNumber}";
                Assert.That(annotation, Is.EqualTo(product.Annotation));
                string toString =
                    $"{product.ProductType}{product.FragmentNumber};{product.NeutralMass:F5}-{product.NeutralLoss:0.##}";
                Assert.That(toString, Is.EqualTo(product.ToString()));
            }
        }

        [Test]
        [TestCaseSource(nameof(GetSixMerIndividualFragmentTypeTestCases))]
        public void TestRnaFragmentNumbers(SixmerTestCase testCase)
        {
            RNA rna = new("GUACUG");
            List<RnaProduct> products = rna.GetNeutralFragments(testCase.Type).Select(p => (RnaProduct)p).ToList();

            for (int i = 0; i < products.Count; i++)
            {
                var product = products[i];
                bool isThreePrime = product.ProductType.GetRnaTerminusType() == FragmentationTerminus.ThreePrime;

                int fragmentNumber = i + 1;
                int residuePosition = isThreePrime ? rna.Length - fragmentNumber : fragmentNumber;

                Assert.That(product.FragmentNumber, Is.EqualTo(fragmentNumber));
                Assert.That(product.ResiduePosition, Is.EqualTo(residuePosition));
            }

        }

        [Test]
        public void TestConstructorAndEquality()
        {
            RnaProduct product1 = new RnaProduct(ProductType.d, FragmentationTerminus.FivePrime, 200, 4, 4, 0.0);
            RnaProduct product2 = new RnaProduct(ProductType.d, FragmentationTerminus.FivePrime, 200, 4, 4, 0.0);
            RnaProduct uniqueProduct = new RnaProduct(ProductType.a, FragmentationTerminus.FivePrime, 200, 4, 4, 0.0);

            Assert.That(product1.Equals(product1));
            Assert.That(product1.Equals(product2));
            Assert.That(product1.GetHashCode(), Is.EqualTo(product2.GetHashCode()));
            Assert.That(!product1.Equals(uniqueProduct));
            Assert.That(!product1.Equals(null));
            Assert.That(product1.GetHashCode(), Is.Not.EqualTo(uniqueProduct.GetHashCode()));

            Assert.That(product1.Equals((object)product1));
            Assert.That(product1.Equals((object)product2));
            Assert.That(!product1.Equals((object)uniqueProduct));
            Assert.That(!product1.Equals((object)new Product(ProductType.d, FragmentationTerminus.N, 200, 4, 4, 0.0)));
            Assert.That(!product1.Equals((object)null));
        }
    }
}
