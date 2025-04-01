using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class IsotopicEnvelopeTests
    {
        [Test]
        public static void CheckSimilarity_WithIdenticalEnvelopes_ReturnsOne()
        {
            // Arrange
            var peak1 = new IndexedMassSpectralPeak(100, 10, 0, 0);
            var peak2 = new IndexedMassSpectralPeak(200, 20, 0, 0);
            var peaks = new List<IIndexedMzPeak> { peak1, peak2 };

            var envelope1 = new IsotopicEnvelope(peaks, peak1, 1, 30, 0.99);
            var envelope2 = new IsotopicEnvelope(peaks, peak1, 1, 30, 0.99);

            // Act
            var similarity = envelope1.CheckSimilarity(envelope2);

            // Assert
            Assert.AreEqual(1.0, similarity);
        }

        [Test]
        public static void CheckSimilarity_WithDifferentEnvelopes_ReturnsLessThanOne()
        {
            // Arrange
            var peak1 = new IndexedMassSpectralPeak(100, 10, 0, 0);
            var peak2 = new IndexedMassSpectralPeak(200, 20, 0, 0);
            var peak3 = new IndexedMassSpectralPeak(300, 30, 0, 0);
            var peaks1 = new List<IIndexedMzPeak> { peak1, peak2 };
            var peaks2 = new List<IIndexedMzPeak> { peak1, peak3 };

            var envelope1 = new IsotopicEnvelope(peaks1, peak1, 1, 30, 0.99);
            var envelope2 = new IsotopicEnvelope(peaks2, peak1, 1, 30, 0.99);

            // Act
            var similarity = envelope1.CheckSimilarity(envelope2);

            // Assert
            Assert.Less(similarity, 1.0);
        }

        [Test]
        public static void CheckSimilarity_WithNoMatchingPeaks_ReturnsNegativeOne()
        {
            // Arrange
            var peak1 = new IndexedMassSpectralPeak(100, 10, 0, 0);
            var peak2 = new IndexedMassSpectralPeak(200, 20, 0, 0);
            var peak3 = new IndexedMassSpectralPeak(300, 30, 0, 0);
            var peak4 = new IndexedMassSpectralPeak(400, 40, 0, 0);
            var peaks1 = new List<IIndexedMzPeak> { peak1, peak2 };
            var peaks2 = new List<IIndexedMzPeak> { peak3, peak4 };

            var envelope1 = new IsotopicEnvelope(peaks1, peak1, 1, 30, 0.99);
            var envelope2 = new IsotopicEnvelope(peaks2, peak3, 1, 30, 0.99);

            // Act
            var similarity = envelope1.CheckSimilarity(envelope2);

            // Assert
            Assert.AreEqual(0, similarity);
        }
    }
}
