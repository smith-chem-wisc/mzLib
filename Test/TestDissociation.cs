using MassSpectrometry;
using NUnit.Framework;
using System.Collections.Generic;
using Proteomics.Fragmentation;

namespace Test
{
    [TestFixture]
    public sealed class TestDissociation
    {
        [Test]
        public void TestHCD_ProductTypes()
        {

            List<ProductType> d = DissociationTypeCollection.ProductsFromDissociationType[DissociationType.HCD];
            Assert.IsTrue(d.Contains(ProductType.b));
            Assert.IsTrue(d.Contains(ProductType.y));
        }

        [Test]
        public void TestECD_ProductTypes()
        {
            List<ProductType> d = DissociationTypeCollection.ProductsFromDissociationType[DissociationType.ECD];
            Assert.IsTrue(d.Contains(ProductType.c));
            Assert.IsTrue(d.Contains(ProductType.y));
            Assert.IsTrue(d.Contains(ProductType.zPlusOne));
        }

        [Test]
        public void TestETD_ProductTypes()
        {
            List<ProductType> d = DissociationTypeCollection.ProductsFromDissociationType[DissociationType.ETD];
            Assert.IsTrue(d.Contains(ProductType.c));
            Assert.IsTrue(d.Contains(ProductType.y));
            Assert.IsTrue(d.Contains(ProductType.zPlusOne));
        }

        [Test]
        public void TestEThCD_ProductTypes()
        {
            List<ProductType> d = DissociationTypeCollection.ProductsFromDissociationType[DissociationType.EThcD];
            Assert.IsTrue(d.Contains(ProductType.b));
            Assert.IsTrue(d.Contains(ProductType.y));
        }

        [Test]
        public void TestAny_ProductTypes()
        {
            List<ProductType> d = DissociationTypeCollection.ProductsFromDissociationType[DissociationType.AnyActivationType];
            Assert.IsTrue(d.Contains(ProductType.b));
            Assert.IsTrue(d.Contains(ProductType.y));
        }

        [Test]
        public void TestCustom_ProductTypes()
        {

        }

        [Test]
        public void TestHCD_ProductTypes_IgnorePassedList()
        {

        }
    }
}