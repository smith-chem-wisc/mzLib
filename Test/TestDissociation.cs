using MassSpectrometry;
using NUnit.Framework;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public sealed class TestDissociation
    {
        [Test]
        public void TestHCD_ProductTypes()
        {

            List<ProductType> d = DissociationTypeCollection.ProductsFromDissociationType[DissociationType.HCD];
            Assert.IsTrue(d.Contains(ProductType.BnoB1ions));
            Assert.IsTrue(d.Contains(ProductType.Y));
        }

        [Test]
        public void TestECD_ProductTypes()
        {
            List<ProductType> d = DissociationTypeCollection.ProductsFromDissociationType[DissociationType.ECD];
            Assert.IsTrue(d.Contains(ProductType.C));
            Assert.IsTrue(d.Contains(ProductType.Y));
            Assert.IsTrue(d.Contains(ProductType.Zdot));
        }

        [Test]
        public void TestETD_ProductTypes()
        {
            List<ProductType> d = DissociationTypeCollection.ProductsFromDissociationType[DissociationType.ETD];
            Assert.IsTrue(d.Contains(ProductType.C));
            Assert.IsTrue(d.Contains(ProductType.Y));
            Assert.IsTrue(d.Contains(ProductType.Zdot));
        }

        [Test]
        public void TestEThCD_ProductTypes()
        {
            List<ProductType> d = DissociationTypeCollection.ProductsFromDissociationType[DissociationType.EThCD];
            Assert.IsTrue(d.Contains(ProductType.BnoB1ions));
            Assert.IsTrue(d.Contains(ProductType.Y));
        }

        [Test]
        public void TestAny_ProductTypes()
        {
            List<ProductType> d = DissociationTypeCollection.ProductsFromDissociationType[DissociationType.AnyActivationType];
            Assert.IsTrue(d.Contains(ProductType.B));
            Assert.IsTrue(d.Contains(ProductType.Y));
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