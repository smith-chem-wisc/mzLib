using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestSpectralSimilarity
    {
        [Test]
        public void TestAllSpectrumSimilarities()
        {
            double ppmTolerance = 10;
            MzSpectrum primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            MzSpectrum secondary = new MzSpectrum(new double[] { 1 }, new double[] { 2 }, false);

            SpectralSimilarity s = new SpectralSimilarity(primary, secondary, ppmTolerance);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(1).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.29).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.17).Within(0.01));
        }
    }
}