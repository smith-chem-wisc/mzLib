using MassSpectrometry;
using NUnit.Framework;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public sealed class TestDissassociation
    {
        [Test]
        public void TestHCD_ProductTypes()
        {
            DissassociationTypeCollection d = new DissassociationTypeCollection(DissassociationType.HCD);
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.BnoB1ions));
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.Y));
        }

        [Test]
        public void TestECD_ProductTypes()
        {
            DissassociationTypeCollection d = new DissassociationTypeCollection(DissassociationType.ECD);
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.C));
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.Y));
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.Z));
        }

        [Test]
        public void TestETD_ProductTypes()
        {
            DissassociationTypeCollection d = new DissassociationTypeCollection(DissassociationType.ETD);
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.C));
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.Y));
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.Z));
        }

        [Test]
        public void TestEThCD_ProductTypes()
        {
            DissassociationTypeCollection d = new DissassociationTypeCollection(DissassociationType.EThCD);
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.BnoB1ions));
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.Y));
        }

        [Test]
        public void TestAny_ProductTypes()
        {
            DissassociationTypeCollection d = new DissassociationTypeCollection(DissassociationType.Any);
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.B));
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.Y));
        }

        [Test]
        public void TestCustom_ProductTypes()
        {
            DissassociationTypeCollection d = new DissassociationTypeCollection(DissassociationType.Custom, new List<ProductType> { ProductType.BnoB1ions, ProductType.Y });
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.BnoB1ions));
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.Y));
        }

        [Test]
        public void TestHCD_ProductTypes_IgnorePassedList()
        {
            DissassociationTypeCollection d = new DissassociationTypeCollection(DissassociationType.HCD, new List<ProductType> { ProductType.C, ProductType.Z });
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.BnoB1ions));
            Assert.IsTrue(d.ProductTypes.Contains(ProductType.Y));
            Assert.IsFalse(d.ProductTypes.Contains(ProductType.C));
            Assert.IsFalse(d.ProductTypes.Contains(ProductType.Z));
        }
    }
}