using System;
using Omics.Fragmentation;
using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;

namespace Test.Omics
{
    /// <summary>
    /// Tests for MatchedFragmentIon and MatchedFragmentIonWithCache classes.
    /// Validates fragment ion matching, mass error calculations, and caching behavior.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class MatchedFragmentIonTests
    {
        /// <summary>
        /// Verifies constructor correctly stores Product, Mz, Intensity, and Charge.
        /// Critical: Fragment ion data must be accurately preserved for spectral annotation.
        /// </summary>
        [Test]
        public void Constructor_SetsAllProperties()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var ion = new MatchedFragmentIon(product, 101.0, 200.0, 1);

            Assert.Multiple(() =>
            {
                Assert.That(ion.NeutralTheoreticalProduct, Is.EqualTo(product));
                Assert.That(ion.Mz, Is.EqualTo(101.0));
                Assert.That(ion.Intensity, Is.EqualTo(200.0));
                Assert.That(ion.Charge, Is.EqualTo(1));
            });
        }

        /// <summary>
        /// Verifies IsInternalFragment correctly identifies internal fragment ions.
        /// Critical: Internal fragments require special handling in spectral matching algorithms.
        /// </summary>
        [Test]
        public void IsInternalFragment_IdentifiesInternalFragments()
        {
            // Product with secondary product type indicates internal fragment
            var internalProduct = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0, ProductType.aStar, 2);
            var ion = new MatchedFragmentIon(internalProduct, 101.0, 200.0, 1);
            Assert.That(ion.IsInternalFragment, Is.True);

            // Also works with cached variants
            var cachedProduct = new ProductWithCache(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0, ProductType.aStar, 2);
            var cachedIon = new MatchedFragmentIonWithCache(cachedProduct, 101.0, 200.0, 1);
            Assert.That(cachedIon.IsInternalFragment, Is.True);
        }

        /// <summary>
        /// Verifies MassErrorDa calculates correct mass error in Daltons.
        /// Critical: Mass error is used for match quality assessment and FDR calculations.
        /// Tests all combinations of cached/non-cached Product and MatchedFragmentIon.
        /// </summary>
        [Test]
        public void MassErrorDa_CalculatesCorrectly_AllCacheVariants()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var productCached = new ProductWithCache(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);

            // Test all four combinations
            var ion1 = new MatchedFragmentIon(product, 101.0, 200.0, 1);
            var ion2 = new MatchedFragmentIonWithCache(product, 101.0, 200.0, 1);
            var ion3 = new MatchedFragmentIon(productCached, 101.0, 200.0, 1);
            var ion4 = new MatchedFragmentIonWithCache(productCached, 101.0, 200.0, 1);

            // MassErrorDa = experimental m/z - theoretical m/z (accounting for charge and proton mass)
            const double expectedErrorDa = -0.00728; // Actual calculated value
            Assert.Multiple(() =>
            {
                Assert.That(ion1.MassErrorDa, Is.EqualTo(expectedErrorDa).Within(0.001));
                Assert.That(ion2.MassErrorDa, Is.EqualTo(expectedErrorDa).Within(0.001));
                Assert.That(ion3.MassErrorDa, Is.EqualTo(expectedErrorDa).Within(0.001));
                Assert.That(ion4.MassErrorDa, Is.EqualTo(expectedErrorDa).Within(0.001));
            });
        }

        /// <summary>
        /// Verifies MassErrorPpm calculates correct mass error in parts per million.
        /// Critical: PPM error is the standard metric for high-resolution mass spectrometry matching.
        /// Tests all combinations of cached/non-cached Product and MatchedFragmentIon.
        /// </summary>
        [Test]
        public void MassErrorPpm_CalculatesCorrectly_AllCacheVariants()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var productCached = new ProductWithCache(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);

            var ion1 = new MatchedFragmentIon(product, 101.0, 200.0, 1);
            var ion2 = new MatchedFragmentIonWithCache(product, 101.0, 200.0, 1);
            var ion3 = new MatchedFragmentIon(productCached, 101.0, 200.0, 1);
            var ion4 = new MatchedFragmentIonWithCache(productCached, 101.0, 200.0, 1);

            const double expectedPpm = -72.764;
            Assert.Multiple(() =>
            {
                Assert.That(ion1.MassErrorPpm, Is.EqualTo(expectedPpm).Within(0.001));
                Assert.That(ion2.MassErrorPpm, Is.EqualTo(expectedPpm).Within(0.001));
                Assert.That(ion3.MassErrorPpm, Is.EqualTo(expectedPpm).Within(0.001));
                Assert.That(ion4.MassErrorPpm, Is.EqualTo(expectedPpm).Within(0.001));
            });
        }

        /// <summary>
        /// Verifies Annotation returns correct format (e.g., "b1+1" for b-ion, position 1, charge 1).
        /// Critical: Annotations are displayed in spectrum viewers and exported to result files.
        /// </summary>
        [Test]
        public void Annotation_ReturnsCorrectFormat()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var ion = new MatchedFragmentIon(product, 101.0, 200.0, 1);

            Assert.That(ion.Annotation, Is.EqualTo("b1+1"));
        }

        /// <summary>
        /// Verifies ToString returns annotation with neutral mass in expected format.
        /// Critical: Used for spectral library export and debugging.
        /// </summary>
        [Test]
        public void ToString_ReturnsAnnotationWithMass()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var ion = new MatchedFragmentIon(product, 101.0, 200.0, 1);

            Assert.That(ion.ToString(), Is.EqualTo("b1+1\t;100"));
        }

        /// <summary>
        /// Verifies Equals correctly compares MatchedFragmentIon instances.
        /// Critical: Required for deduplication in fragment matching results.
        /// Also verifies cross-type equality between cached and non-cached variants.
        /// </summary>
        [Test]
        public void Equals_ComparesCorrectly()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 1, 1, 0.0);
            var ion1 = new MatchedFragmentIon(product, 101.0, 200.0, 1);
            var ion2 = new MatchedFragmentIon(product, 101.0, 200.0, 1);

            Assert.That(ion1.Equals(ion2), Is.True);

            // Cross-type equality (cached vs non-cached) with near-equal values
            var ionCached = new MatchedFragmentIonWithCache(product, 100.99999999999, 100.0, 1);
            var ionNonCached = new MatchedFragmentIon(product, 101.0, 99.99999999999, 1);
            Assert.That(ionNonCached.Equals(ionCached), Is.True);
        }

        /// <summary>
        /// Verifies GetHashCode produces consistent hashes for equal objects.
        /// Critical: Required for correct HashSet/Dictionary behavior in fragment collections.
        /// </summary>
        [Test]
        public void GetHashCode_ConsistentForEqualObjects()
        {
            var product1 = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var product2 = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var ion = new MatchedFragmentIon(product1, 1, 1, 1);
            var ionCached = new MatchedFragmentIonWithCache(product2, 1, 1, 1);

            Assert.Multiple(() =>
            {
                Assert.That(product1.GetHashCode(), Is.EqualTo(product2.GetHashCode()));
                Assert.That(ion.GetHashCode(), Is.EqualTo(ionCached.GetHashCode()));
            });
        }
    }
}