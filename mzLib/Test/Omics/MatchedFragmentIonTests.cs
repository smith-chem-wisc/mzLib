using System;
using Omics.Fragmentation;
using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;

namespace Test.Omics
{
    [ExcludeFromCodeCoverage]
    public class MatchedFragmentIonTests
    {
        [Test]
        public void TestMatchedFragmentIonConstructor() 
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var ion = new MatchedFragmentIon(product, 101.0, 200.0, 1);

            Assert.That(ion.NeutralTheoreticalProduct, Is.EqualTo(product));
            Assert.That(ion.Mz, Is.EqualTo(101.0));
            Assert.That(ion.Intensity, Is.EqualTo(200.0));
            Assert.That(ion.Charge, Is.EqualTo(1));
        }

        [Test]
        public void TestIsInternalFragment()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0, ProductType.aStar, 2);
            var ion = new MatchedFragmentIon(product, 101.0, 200.0, 1);

            Assert.That(ion.IsInternalFragment, Is.True);

            product = new ProductWithCache(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0, ProductType.aStar, 2);
            ion = new MatchedFragmentIonWithCache(product, 101.0, 200.0, 1);

            Assert.That(ion.IsInternalFragment, Is.True);
        }

        [Test]
        public void TestMassErrorDa()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var ion = new MatchedFragmentIon(product, 101.0, 200.0, 1);
            Assert.That(ion.MassErrorDa, Is.EqualTo(1.0).Within(10));

            product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            ion = new MatchedFragmentIonWithCache(product, 101.0, 200.0, 1);
            Assert.That(ion.MassErrorDa, Is.EqualTo(1.0).Within(10));
            Assert.That(ion.MassErrorDa, Is.EqualTo(1.0).Within(10));

            product = new ProductWithCache(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            ion = new MatchedFragmentIon(product, 101.0, 200.0, 1);
            Assert.That(ion.MassErrorDa, Is.EqualTo(1.0).Within(10));
            Assert.That(ion.MassErrorDa, Is.EqualTo(1.0).Within(10));

            product = new ProductWithCache(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            ion = new MatchedFragmentIonWithCache(product, 101.0, 200.0, 1);
            Assert.That(ion.MassErrorDa, Is.EqualTo(1.0).Within(10));
            Assert.That(ion.MassErrorDa, Is.EqualTo(1.0).Within(10));
        }

        [Test]
        public void TestMassErrorPpm()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var ion = new MatchedFragmentIon(product, 101.0, 200.0, 1);
            Assert.That(ion.MassErrorPpm, Is.EqualTo(-72.764).Within(0.001));

            product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            ion = new MatchedFragmentIonWithCache(product, 101, 200.0, 1);
            Assert.That(ion.MassErrorPpm, Is.EqualTo(-72.764).Within(0.001));
            Assert.That(ion.MassErrorPpm, Is.EqualTo(-72.764).Within(0.001));

            product = new ProductWithCache(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            ion = new MatchedFragmentIon(product, 101, 200.0, 1);
            ion = new MatchedFragmentIon(product, 101, 200.0, 1);
            Assert.That(ion.MassErrorPpm, Is.EqualTo(-72.764).Within(0.001));

            product = new ProductWithCache(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            ion = new MatchedFragmentIonWithCache(product, 101, 200.0, 1);
            Assert.That(ion.MassErrorPpm, Is.EqualTo(-72.764).Within(0.001));
            Assert.That(ion.MassErrorPpm, Is.EqualTo(-72.764).Within(0.001));
        }

        [Test]
        public void TestAnnotation()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var ion = new MatchedFragmentIon(product, 101.0, 200.0, 1);

            Assert.That(ion.Annotation, Is.EqualTo("b1+1"));
        }

        [Test]
        public void TestToString()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var ion = new MatchedFragmentIon(product, 101.0, 200.0, 1);

            Assert.That(ion.ToString(), Is.EqualTo("b1+1\t;100"));

            Product P = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            MatchedFragmentIon m = new MatchedFragmentIon(P, 1, 1, 1);
            Assert.That(m.ToString(), Is.EqualTo("b1+1\t;1"));
        }

        [Test]
        public void TestEquals()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var ion1 = new MatchedFragmentIon(product, 101.0, 200.0, 1);
            var ion2 = new MatchedFragmentIon(product, 101.0, 200.0, 1);

            Assert.That(ion1.Equals(ion2), Is.True);
        }

        [Test]
        public static void TestMatchedFragmentAndCachedIonEquals()
        {
            Product P = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            MatchedFragmentIon ion1 = new MatchedFragmentIon(P, experMz: 150, experIntensity: 99.99999999999, charge: 2);
            MatchedFragmentIon ion2 = new MatchedFragmentIonWithCache(P, experMz: 149.99999999999, experIntensity: 100, charge: 2);
            Assert.That(ion1, Is.EqualTo(ion2));
        }

        [Test]
        public void TestGetHashCode()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var ion = new MatchedFragmentIon(product, 101.0, 200.0, 1);

            Assert.That(ion.GetHashCode(), Is.EqualTo(HashCode.Combine(product.GetHashCode(), 1.GetHashCode(), Math.Round(101.0, 10).GetHashCode(), Math.Round(200.0, 6).GetHashCode())));
        }

        [Test]
        public static void Test_MatchedFragmentGetHashCode()
        {
            Product P = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product pPrime = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            MatchedFragmentIon m = new MatchedFragmentIon(P, 1, 1, 1);
            MatchedFragmentIon mPrime = new MatchedFragmentIonWithCache(pPrime, 1, 1, 1);
            Assert.That(P.GetHashCode(), Is.EqualTo(pPrime.GetHashCode()));
            Assert.That(mPrime.GetHashCode(), Is.EqualTo(m.GetHashCode()));
        }
    }
}
