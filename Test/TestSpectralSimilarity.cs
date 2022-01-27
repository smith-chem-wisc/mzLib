using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestSpectralSimilarity
    {
        [Test]
        public void TestAllSpectrumSimilaritiesWithoutMzFilter()
        {
            double ppmTolerance = 10;

            //Test all different similarity calculations
            MzSpectrum primary = new MzSpectrum(new double[] { 1, 2, 3, 4, 5 }, new double[] { 2, 4, 6, 8, 10 }, false);
            MzSpectrum secondary = new MzSpectrum(new double[] { 3, 4, 5, 6, 7 }, new double[] { 9, 7, 5, 3, 1 }, false);
            SpectralSimilarity s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, true,0);
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
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, true,0);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.37).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.22).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.33).Within(0.01));

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true,0);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-.03).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.17).Within(0.01));

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, true,0);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.41).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.07).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.24).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.90).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.24).Within(0.01));

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.unnormalized, ppmTolerance, true,0);
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
                s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true, 0);
            }, "Empty YArray in spectrum.");

            //We should have any zero intensity YArrays but just to be sure
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1 }, new double[] { 0 }, false);

            Assert.Throws<MzLibException>(() =>
            {
                s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true, 0);
            }, "Spectrum has no intensity.");

            //What happens when all intensity pairs include a zero
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 4 }, new double[] { 2 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true,0);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.18).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.77).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));

            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 0, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 2, 0, 2, 2 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true,0);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.48).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.32).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.33).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.20).Within(0.01));

            //Test what happens when all intensity pairs include 1 zero
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 0, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 4, 5 }, new double[] { 2, 0 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true,0);
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

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true,0);
            Assert.AreEqual(6, s.intensityPairs.Count);

            //Test alternate constructor
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, true,0);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.37).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.22).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.33).Within(0.01));

            //Test alternate constructor only library peaks. Since library has one peak, and primary has three peaks, we get only one intensity pair
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false,0);
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
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false,0);
            Assert.AreEqual(3, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
        }

        [Test]
        public void TestAllSpectrumSimilaritiesWithDefaultedMzFilter()
        {
            double ppmTolerance = 10;
            MzSpectrum primary = new MzSpectrum(new double[] { 100, 200, 300, 400, 500 ,600,700}, new double[] { 2, 4, 6, 8, 10,12,14 }, false);
            MzSpectrum secondary = new MzSpectrum(new double[] { 200, 300, 500, 600 ,800}, new double[] { 9, 7, 5, 3, 1 }, false);
            //Test when using all peaks of primary(experimental) and secondary(theoretical)  spectra (bool allpeaks is true) and mz cut off is 0 (no cut off)
            SpectralSimilarity s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, true,0);
            Assert.AreEqual(8, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.68).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.48).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.65).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.56).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(0.02).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.13).Within(0.01));

            //Test when using all peaks of primary(experimental) and secondary(theoretical) spectra (bool allpeaks is true) and mz cut off is 300 (default cut off)
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, true,300);
            Assert.AreEqual(6, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.70).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.49).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.72).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.60).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.13).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.09).Within(0.01));

            //Test when not using all peaks of primary(experimental) spectra (bool allpeaks is false) and mz cut off is is 0 (no cut off)
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, false, 0);
            Assert.AreEqual(5, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.92).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.75).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.79).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.73).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(0.80).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.12).Within(0.01));

            //Test when not using all peaks of primary(experimental) spectra (bool allpeaks is false) and mz cut off is is 300 (default cut off)
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, false,300);
            Assert.AreEqual(4, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.89).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.69).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.79).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.71).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(0.89).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.09).Within(0.01));

            //Test all different similarity calculations
            primary = new MzSpectrum(new double[] { 100, 200, 300, 400, 500 }, new double[] { 2, 4, 6, 8, 10 }, false);
            secondary = new MzSpectrum(new double[] { 300, 400, 500, 600, 700 }, new double[] { 9, 7, 5, 3, 1 }, false);
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, false);
            Assert.AreEqual(5, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.89).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.70).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.79).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.77).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(0.81).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.17).Within(0.01));


            //Test all normalization schemes
            primary = new MzSpectrum(new double[] { 1000, 2000, 3000 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1000 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, true);
            Assert.AreEqual(3, s.intensityPairs.Count);
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


            //What happens when all intensity pairs include a zero
            primary = new MzSpectrum(new double[] { 100, 200, 300 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 400, 500 }, new double[] { 2, 4 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            Assert.AreEqual(3, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.1).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));

            primary = new MzSpectrum(new double[] { 1000, 2000, 3000 }, new double[] { 0, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1000, 2000, 3000, 4000 }, new double[] { 2, 0, 2, 2 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            Assert.AreEqual(4, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.48).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.32).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.33).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.20).Within(0.01));

            //Test what happens when all intensity pairs include 1 zero
            primary = new MzSpectrum(new double[] { 1000, 2000, 3000 }, new double[] { 0, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 4000, 5000 }, new double[] { 2, 0 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            Assert.AreEqual(5, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.23).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.4).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, false);
            Assert.AreEqual(2, s.intensityPairs.Count);



            //Test alternate constructor
            primary = new MzSpectrum(new double[] { 100, 350, 400 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 350 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, true);
            Assert.AreEqual(2, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.55).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.37).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.05).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.50).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.67).Within(0.01));

            //Test alternate constructor only library peaks. Since library has one peak, and primary has three peaks, we get only one intensity pair
            primary = new MzSpectrum(new double[] { 350, 400, 500 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 400 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false);
            Assert.AreEqual(1, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.67).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.80).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1.0).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.67).Within(0.01));

            //Test cosine similarity when there are no peaks from spectrum one matching spectrum 2
            primary = new MzSpectrum(new double[] { 450, 650, 850 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 400, 600, 800 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false);
            Assert.AreEqual(3, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));

            //Test when intensityPairs is null 
            //Test when all mz are less than 300 in experimental peaks
            primary = new MzSpectrum(new double[] { 100, 150, 200 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 400, 600, 800 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false);
            Assert.AreEqual(1, s.intensityPairs.Count);
            Assert.AreEqual(new List<(double, double)> { (-1, -1) }, s.intensityPairs);
            Assert.That(s.CosineSimilarity(), Is.Null);
            Assert.That(s.SpectralContrastAngle(), Is.Null);
            Assert.That(s.EuclideanDistance(), Is.Null);
            Assert.That(s.BrayCurtis(), Is.Null);
            Assert.That(s.PearsonsCorrelation(), Is.Null);
            Assert.That(s.DotProduct(), Is.Null);

            //Test when all mz are less than 300 in theoretical peaks
            primary = new MzSpectrum(new double[] { 150, 250, 350 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 100, 150, 200 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false);
            Assert.AreEqual(1, s.intensityPairs.Count);
            Assert.AreEqual(new List<(double, double)> { (-1, -1) }, s.intensityPairs);
            Assert.That(s.CosineSimilarity(), Is.Null);
            Assert.That(s.SpectralContrastAngle(), Is.Null);
            Assert.That(s.EuclideanDistance(), Is.Null);
            Assert.That(s.BrayCurtis(), Is.Null);
            Assert.That(s.PearsonsCorrelation(), Is.Null);
            Assert.That(s.DotProduct(), Is.Null);
        }
    }
}