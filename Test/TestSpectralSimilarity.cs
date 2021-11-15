using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
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

            SpectralSimilarity s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(1).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.5).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.33).Within(0.01));

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(1).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.29).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.17).Within(0.01));

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(1).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.24).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.39).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.24).Within(0.01));

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.unnormalized, ppmTolerance);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(1).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.24).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.39).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.24).Within(0.01));

            //Make sure there is no crash w/ empty array
            primary = new MzSpectrum(new double[] {  }, new double[] {  }, false);
            secondary = new MzSpectrum(new double[] {  }, new double[] {  }, false);

            Assert.Throws<MzLibException>(() =>
            {
                s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance);
            }, "Empty YArray in spectrum.");

            //We should have any zero intensity YArrays but just to be sure
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1 }, new double[] { 0 }, false);

            Assert.Throws<MzLibException>(() =>
            {
                s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance);
            }, "Spectrum has no intensity.");

            //What happens when there are no intensity pairs
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 4 }, new double[] { 2 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.41).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));
        }
    }
}