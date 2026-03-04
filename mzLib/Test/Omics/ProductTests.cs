using NUnit.Framework;
using Omics.Fragmentation;
using System.Diagnostics.CodeAnalysis;

namespace Test.Omics
{
    [ExcludeFromCodeCoverage]
    public class ProductTests
    {
        [Test]
        public void TestProductConstructor()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            Assert.That(product.NeutralMass, Is.EqualTo(100.0));
            Assert.That(product.ProductType, Is.EqualTo(ProductType.b));
            Assert.That(product.NeutralLoss, Is.EqualTo(0.0));
            Assert.That(product.Terminus, Is.EqualTo(FragmentationTerminus.N));
            Assert.That(product.FragmentNumber, Is.EqualTo(1));
            Assert.That(product.ResiduePosition, Is.EqualTo(1));
            Assert.That(product.SecondaryProductType, Is.Null);
            Assert.That(product.SecondaryFragmentNumber, Is.EqualTo(0));
            Assert.That(product.IsInternalFragment, Is.False);
        }

        [Test]
        public void TestProductAnnotation()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            Assert.That(product.Annotation, Is.EqualTo("b1"));
            Assert.That(product.IsInternalFragment, Is.False);

            var productWithLoss = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 18.0);
            Assert.That(productWithLoss.Annotation, Is.EqualTo("b1-18.00"));

            var internalProduct = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0, ProductType.y, 2);
            Assert.That(internalProduct.Annotation, Is.EqualTo("bIy[1-2]"));
            Assert.That(internalProduct.IsInternalFragment, Is.True);
        }

        [Test]
        public void TestProductToString()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            Assert.That(product.ToString(), Is.EqualTo("b1;100.00000-0"));
            Assert.That(product.IsInternalFragment, Is.False);

            var productWithLoss = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 18.0);
            Assert.That(productWithLoss.ToString(), Is.EqualTo("b1;100.00000-18"));

            var internalProduct = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0, ProductType.y, 2);
            Assert.That(internalProduct.ToString(), Is.EqualTo("bIy[1-2];100.00000-0"));
            Assert.That(internalProduct.IsInternalFragment, Is.True);
        }

        [Test]
        public void TestProductWithCacheAnnotation()
        {
            var productWithCache = new ProductWithCache(ProductType.b, FragmentationTerminus.N, 100.0, 3, 5, 0.0);
            Assert.That(productWithCache.Annotation, Is.EqualTo("b3"));

            var productWithCacheWithLoss = new ProductWithCache(ProductType.b, FragmentationTerminus.N, 100.0, 3, 5, 18.0);
            Assert.That(productWithCacheWithLoss.Annotation, Is.EqualTo("b3-18.00"));

            var internalProductWithCache = new ProductWithCache(ProductType.y, FragmentationTerminus.C, 200.0, 18, 5, 0.0, ProductType.b, 36);
            Assert.That(internalProductWithCache.Annotation, Is.EqualTo("yIb[18-36]"));
        }

        [Test]
        public void TestProductEquality()
        {
            var product1 = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var product2 = new ProductWithCache(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var product3 = new Product(ProductType.y, FragmentationTerminus.C, 200.0, 2, 2, 0.0);

            Assert.That(product1.Equals(product2), Is.True);
            Assert.That(product1.Equals(product3), Is.False);
            Assert.That(product1.Equals(null), Is.False);
        }

        [Test]
        public void TestProductHashCode()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            Assert.That(product.GetHashCode(), Is.EqualTo(product.NeutralMass.GetHashCode()));


            Product P = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Assert.That(1072693248, Is.EqualTo(P.GetHashCode()));
        }

        [Test]
        public static void Test_ProductMonoisotopicMass()
        {
            Product P = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Assert.That(P.MonoisotopicMass.Equals(P.NeutralMass));
        }
    }
}
