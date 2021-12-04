using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using NUnit.Framework;
using System;

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

            //Test all different similarity calculations
            MzSpectrum primary = new MzSpectrum(new double[] { 1, 2, 3, 4, 5 }, new double[] { 2, 4, 6, 8, 10 }, false);
            MzSpectrum secondary = new MzSpectrum(new double[] { 3, 4, 5, 6, 7 }, new double[] { 9, 7, 5, 3, 1 }, false);
            SpectralSimilarity s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, true);
            Assert.AreEqual(7, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.8).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.59).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.70).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.66).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(0.42).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.17).Within(0.01));

            //Test all normalization schemes
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, true);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.37).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.22).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.33).Within(0.01));

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-.03).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.17).Within(0.01));

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, true);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.41).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.07).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.24).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.90).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.24).Within(0.01));

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.unnormalized, ppmTolerance, true);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.41).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.07).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.24).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.90).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.24).Within(0.01));

            //Make sure there is no crash w/ empty array
            primary = new MzSpectrum(new double[] { }, new double[] { }, false);
            secondary = new MzSpectrum(new double[] { }, new double[] { }, false);

            Assert.Throws<MzLibException>(() =>
            {
                s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            }, "Empty YArray in spectrum.");

            //We should have any zero intensity YArrays but just to be sure
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1 }, new double[] { 0 }, false);

            Assert.Throws<MzLibException>(() =>
            {
                s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            }, "Spectrum has no intensity.");

            //What happens when all intensity pairs include a zero
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 4 }, new double[] { 2 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.18).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.77).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));

            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 0, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 2, 0, 2, 2 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.48).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.32).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.33).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.20).Within(0.01));

            //Test what happens when all intensity pairs include 1 zero
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 0, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 4, 5 }, new double[] { 2, 0 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            Assert.AreEqual(5, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.23).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.40).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));

            //explore bounds of binary search
            primary = new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 1, 2, 3, 4 }, false);
            secondary = new MzSpectrum(new double[] { 1.000009, 1.99999, 3.00004, 3.99995 }, new double[] { 1, 2, 3, 4 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            Assert.AreEqual(6, s.intensityPairs.Count);

            //Test alternate constructor
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, true);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.37).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.22).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.33).Within(0.01));

            //Test alternate constructor only library peaks. Since library has one peak, and primary has three peaks, we get only one intensity pair
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false);
            Assert.AreEqual(1, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.50).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1.0).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.33).Within(0.01));

            //Test cosine similarity when there are no peaks from spectrum one matching spectrum 2
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 4,6,8 }, new double[] { 2,4,6 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false);
            Assert.AreEqual(3, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
        }
    }
}