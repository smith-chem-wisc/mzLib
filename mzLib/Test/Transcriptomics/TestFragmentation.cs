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
using Omics.Fragmentation;
using Omics.Fragmentation.Oligo;
using Omics.Modifications;
using UsefulProteomicsDatabases;

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
            var rna = new RNA("GUACUG")
                .Digest(new RnaDigestionParams(), new List<Modification>(), new List<Modification>())
                .First() as OligoWithSetMods ?? throw new NullReferenceException();

            var neutralFragments = rna.GetNeutralFragments(testCase.Type).ToList();
            for (int i = 1; i < neutralFragments.Count; i++)
            {
                Assert.That(neutralFragments[i].NeutralMass, Is.EqualTo(testCase.NeutralMasses[i]).Within(0.01));
            }
        }


        [Test]
        public void TestFragmentation_Unmodified()
        {
            var rna = new RNA("GUACUG")
                .Digest(new RnaDigestionParams(), new List<Modification>(), new List<Modification>())
                .First();
            List<Product> products = new();

            // both termini
            rna.Fragment(DissociationType.CID, FragmentationTerminus.Both, products);
            Assert.That(products.Count, Is.EqualTo(31));
            Assert.That(products.All(p => p.ProductType is ProductType.M or ProductType.w or ProductType.dWaterLoss or ProductType.y or ProductType.yWaterLoss or ProductType.aBaseLoss or ProductType.c));

            // only 5'
            rna.Fragment(DissociationType.CID, FragmentationTerminus.FivePrime, products);
            Assert.That(products.Count, Is.EqualTo(15));
            Assert.That(products.All(p => p.ProductType is ProductType.dWaterLoss or  ProductType.aBaseLoss or ProductType.c));

            // only 3'
            rna.Fragment(DissociationType.CID, FragmentationTerminus.ThreePrime, products);
            Assert.That(products.Count, Is.EqualTo(15));
            Assert.That(products.All(p => p.ProductType is ProductType.w or ProductType.y or ProductType.yWaterLoss ));

        }

        [Test]
        [TestCase("GUACUG", "ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A\r\nCF   Na1H-1\r\n" + @"//",
            "GUA[Metal:Sodium on A]CUG", 1874.28, 1896.26, ProductType.a,
            new[] { 267.089, 573.114, 902.167, 1207.208, 1513.233 },
            new[] { 267.089, 573.114, 902.167 + 21.982, 1207.208 + 21.982, 1513.233 + 21.982 })]
        [TestCase("GUACUG", "ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A\r\nCF   Na1H-1\r\n" + @"//",
            "GUA[Metal:Sodium on A]CUG", 1874.28, 1896.26, ProductType.b,
            new[] { 283.084, 589.109, 918.162, 1223.203, 1529.228 },
            new[] { 283.084, 589.109, 918.162 + 21.982, 1223.203 + 21.982, 1529.228 + 21.982 })]
        [TestCase("GUACUG", "ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A\r\nCF   Na1H-1\r\n" + @"//",
            "GUA[Metal:Sodium on A]CUG", 1874.28, 1896.26, ProductType.c,
            new[] { 347.055, 653.081, 982.133, 1287.174, 1593.2 },
            new[] { 347.055, 653.081, 982.133 + 21.982, 1287.174 + 21.982, 1593.2 + 21.982 })]
        [TestCase("GUACUG", "ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A\r\nCF   Na1H-1\r\n" + @"//",
            "GUA[Metal:Sodium on A]CUG", 1874.28, 1896.26, ProductType.d,
            new[] { 363.05, 669.075, 998.128, 1303.169, 1609.195 },
            new[] { 363.05, 669.075, 998.128 + 21.982, 1303.169 + 21.982, 1609.195 + 21.982 })]
        [TestCase("GUACUG", "ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A\r\nCF   Na1H-1\r\n" + @"//",
            "GUA[Metal:Sodium on A]CUG", 1874.28, 1896.26, ProductType.dWaterLoss,
            new[] { 345.039, 651.064, 980.116, 1285.157, 1591.184 },
            new[] { 345.039, 651.064, 980.116 + 21.982, 1285.157 + 21.982, 1591.184 + 21.982 })]
        [TestCase("GUACUG", "ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A\r\nCF   Na1H-1\r\n" + @"//",
            "GUA[Metal:Sodium on A]CUG", 1874.28, 1896.26, ProductType.w,
            new[] { 363.049, 669.074, 974.115, 1303.169, 1609.195 },
            new[] { 363.049, 669.074, 974.115, 1303.169 + 21.982, 1609.195 + 21.982 })]
        [TestCase("GUACUG", "ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A\r\nCF   Na1H-1\r\n" + @"//",
            "GUA[Metal:Sodium on A]CUG", 1874.28, 1896.26, ProductType.x,
            new[] { 347.055, 653.081, 958.122, 1287.174, 1593.2 },
            new[] { 347.055, 653.081, 958.122, 1287.174 + 21.982, 1593.2 + 21.982 })]
        [TestCase("GUACUG", "ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A\r\nCF   Na1H-1\r\n" + @"//",
            "GUA[Metal:Sodium on A]CUG", 1874.28, 1896.26, ProductType.y,
            new[] { 283.084, 589.109, 894.15, 1223.203, 1529.228 },
            new[] { 283.084, 589.109, 894.15, 1223.203 + 21.982, 1529.228 + 21.982 })]
        [TestCase("GUACUG", "ID   Sodium\r\nMT   Metal\r\nPP   Anywhere.\r\nTG   A\r\nCF   Na1H-1\r\n" + @"//",
            "GUA[Metal:Sodium on A]CUG", 1874.28, 1896.26, ProductType.z,
            new[] { 267.089, 573.124, 878.156, 1207.208, 1513.233 },
            new[] { 267.089, 573.124, 878.156, 1207.208 + 21.982, 1513.233 + 21.982 })]
        public void TestFragmentation_Modified(string sequence, string modString, string fullSequence, double unmodifiedMass, double modifiedMass,
            ProductType productType, double[] unmodifiedFragmentMass, double[] modifiedFragmentMasses)
        {
            var mods = PtmListLoader.ReadModsFromString(modString, out List<(Modification, string)> modsOut).ToList();
            var rna = new RNA(sequence);

            var unmodifiedOligo = rna.Digest(new RnaDigestionParams(), new List<Modification>(), new List<Modification>())
                .First() as OligoWithSetMods ?? throw new NullReferenceException();
            Assert.That(unmodifiedOligo.AllModsOneIsNterminus.Count, Is.EqualTo(0));
            Assert.That(unmodifiedOligo.FullSequence, Is.EqualTo(sequence));
            Assert.That(unmodifiedOligo.MonoisotopicMass, Is.EqualTo(unmodifiedMass).Within(0.01));

            var modifiedOligo = rna.Digest(new RnaDigestionParams(), mods, new List<Modification>())
                .First() as OligoWithSetMods ?? throw new NullReferenceException();
            Assert.That(modifiedOligo.AllModsOneIsNterminus.Count, Is.EqualTo(mods.Count));
            Assert.That(modifiedOligo.FullSequence, Is.EqualTo(fullSequence));
            Assert.That(modifiedOligo.MonoisotopicMass, Is.EqualTo(modifiedMass).Within(0.01));

            var unmodifiedProducts = unmodifiedOligo.GetNeutralFragments(productType).ToList();
            Assert.That(unmodifiedProducts.Count, Is.EqualTo(5));
            var modifiedProducts = modifiedOligo.GetNeutralFragments(productType).ToList();
            Assert.That(modifiedProducts.Count, Is.EqualTo(5));


            for (int i = 0; i < unmodifiedProducts.Count; i++)
            {
                var unModifedProduct = unmodifiedProducts[i];
                var modifiedProduct = modifiedProducts[i];

                Assert.That(unModifedProduct.NeutralMass, Is.EqualTo(unmodifiedFragmentMass[i]).Within(0.01));
                Assert.That(modifiedProduct.NeutralMass, Is.EqualTo(modifiedFragmentMasses[i]).Within(0.01));
            }
        }


        [Test]
        [TestCaseSource(nameof(GetSixMerIndividualFragmentTypeTestCases))]
        public void TestRnaFragments(SixmerTestCase testCase)
        {
            var rna = new RNA("GUACUG")
                .Digest(new RnaDigestionParams(), new List<Modification>(), new List<Modification>())
                .First() as OligoWithSetMods ?? throw new NullReferenceException();
            List<Product> products = rna.GetNeutralFragments(testCase.Type).Select(p => (Product)p).ToList();

            for (int i = 0; i < products.Count; i++)
            {
                var product = products[i];
                Assert.That(testCase.Type, Is.EqualTo(product.ProductType));
                Assert.That(testCase.Type.GetRnaTerminusType(), Is.EqualTo(product.Terminus));
                Assert.That(testCase.NeutralMasses[i], Is.EqualTo(product.NeutralMass).Within(0.01));
                Assert.That(testCase.NeutralMasses[i], Is.EqualTo(product.MonoisotopicMass).Within(0.01));
                Assert.That(0, Is.EqualTo(product.NeutralLoss));
                Assert.That(null, Is.EqualTo(product.SecondaryProductType));
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
            var rna = new RNA("GUACUG")
                .Digest(new RnaDigestionParams(), new List<Modification>(), new List<Modification>())
                .First() as OligoWithSetMods ?? throw new NullReferenceException();
            List<Product> products = rna.GetNeutralFragments(testCase.Type).Select(p => (Product)p).ToList();

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
            Product product1 = new Product(ProductType.d, FragmentationTerminus.FivePrime, 200, 4, 4, 0.0);
            Product product2 = new Product(ProductType.d, FragmentationTerminus.FivePrime, 200, 4, 4, 0.0);
            Product uniqueProduct = new Product(ProductType.a, FragmentationTerminus.FivePrime, 200, 4, 4, 0.0);

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
