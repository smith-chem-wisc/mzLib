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
            MzSpectrum primary = new(new double[] { 1, 2, 3, 4, 5 }, new double[] { 2, 4, 6, 8, 10 }, false);
            MzSpectrum secondary = new(new double[] { 3, 4, 5, 6, 7 }, new double[] { 9, 7, 5, 3, 1 }, false);
            SpectralSimilarity s = new(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, true,0);

            //mz pairs in tolerance are (3,3), (4,4), (5,5). Since we are using all peaks, we get 7 intensity pairs with 1,2,6 and 7 intensities being paired zero
            Assert.AreEqual(7, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.8).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.59).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.70).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.66).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(0.42).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.SearleSimilarity(), Is.EqualTo(2.4391).Within(0.01));
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true, 0);
            Assert.That(s.SpectralEntropy(), Is.EqualTo(0.79).Within(0.01));

            //Test all normalization schemes
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, true,0);
            //mz pairs in tolerance are (1,1). Since we are using all peaks, we get 3 intensity pairs with 2 and 3 intensities being paired zero
            Assert.AreEqual(3, s.intensityPairs.Count);
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
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-6.21).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.29).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(4).Within(0.01));

            //Make sure there is no crash w/ empty array
            primary = new MzSpectrum(Array.Empty<double>(), new double[] { }, false);
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
            //there are no mz pairs in tolerance, but since we have all peaks, we get 4 intensity pairs with 1,2,3 and 4 intensities being paired zero
            Assert.AreEqual(4, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.18).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.77).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));

            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 0, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 2, 0, 2, 2 }, false);
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true,0);
            //mz pairs in tolerance are (1,1), (2,2), (3,3) and (4,).
            Assert.AreEqual(4, s.intensityPairs.Count);
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
            //there are no mz pairs in tolerance, but since we have all peaks, we get 5 intensity pairs with 1,2,3,4 and 5 intensities being paired zero, The pair with mz 1 should be (0,0) as should the pair with mz 5
            Assert.AreEqual(5, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.23).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.40).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));

            //explore bounds of binary search
            primary = new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 1, 2, 3, 4 }, false);
            secondary = new MzSpectrum(new double[] { 1.000011, 1.99997, 3.000031, 3.99995 }, new double[] { 1, 2, 3, 4 }, false);
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true,0);
            //The ppm difference between the 4 closest pairs are 11, 15,10.3 and 12.5. These are all beyond the 10ppm tolerance that is allowed. Therefore we get 8 intensity pairs with all intensities being paired zero
            Assert.AreEqual(8, s.intensityPairs.Count);

            //Test alternate constructor
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, true,0);
            //mz pairs in tolerance are (1,1). Since we are using all peaks, we get 3 intensity pairs with 2 and 3 intensities being paired zero
            Assert.AreEqual(3, s.intensityPairs.Count);
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
            //mz pairs in tolerance are (1,1). Since we are NOT using all peaks, we get only 1 intensity pair
            Assert.AreEqual(1, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.50).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1.0).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.33).Within(0.01));

            //Test cosine similarity when there are no peaks from spectrum one matching spectrum 2
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 4, 6, 8 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false,0);
            //There are no mz pairs in tolerance. Since we are NOT using all peaks, we get 0 intensity pairs.
            //however, since there are zero pairs, we return as a pair (-1,-1) so we know that and there is no crash.
            Assert.AreEqual(1, s.intensityPairs.Count);
            Assert.AreEqual((-1, -1), s.intensityPairs[0]);
            Assert.IsNull(s.CosineSimilarity());
            Assert.IsNull(s.SpectralContrastAngle());

            //Test SearleSimilarity with both spectra are identical
            primary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false, 0);

            Assert.AreEqual(s.SearleSimilarity(), double.MaxValue);
        }

        [Test]
        public void TestAllSpectrumSimilaritiesWithDefaultedMzFilter()
        {
            double ppmTolerance = 10;
            MzSpectrum primary = new MzSpectrum(new double[] { 100, 200, 300, 400, 500 , 600, 700}, new double[] { 2, 4, 6, 8, 10, 12, 14 }, false);
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
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, true);
            Assert.AreEqual(6, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.70).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.49).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.61).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.58).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.04).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.16).Within(0.01));

            //Test when not using all peaks of primary(experimental) spectra (bool allpeaks is false) and mz cut off is is 0 (no cut off)
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, false, 0);
            //mz pairs in tolerance are (200,200), (300,300), (500,500), (600,600). Since we are NOT using all peaks, we get the four corresponding intensity pairs.
            Assert.AreEqual(4, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.922).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.747).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.780).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.757).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.98).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.126).Within(0.01));

            //Test when not using all peaks of primary(experimental) spectra (bool allpeaks is false) and mz cut off is is 300 (default cut off)
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, false);
            //mz pairs in tolerance are (200,200), (300,300), (500,500), (600,600). Since we are NOT using all peaks but we ARE using cut off we get only 3 intensity pairs.
            Assert.AreEqual(3, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.954).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.806).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.788).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.73).Within(0.801));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-.958).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.164).Within(0.01));

            //Test all different similarity calculations
            primary = new MzSpectrum(new double[] { 100, 200, 300, 400, 500 }, new double[] { 2, 4, 6, 8, 10 }, false);
            secondary = new MzSpectrum(new double[] { 300, 400, 500, 600, 700 }, new double[] { 9, 7, 5, 3, 1 }, false);
            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, ppmTolerance, false);
            //mz pairs in tolerance are (300,300), (400,400), (500,500). Since we are NOT using all peaks, we get 3 intensity pairs. We are using the default 300mz cut off but that has no effect here
            Assert.AreEqual(3, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.976).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.859).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.815).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.852).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-.997).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.244).Within(0.01));


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
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-6.21).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.29).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(4).Within(0.01));


            //What happens when all intensity pairs include a zero
            primary = new MzSpectrum(new double[] { 100, 200, 300 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 400, 500 }, new double[] { 2, 4 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            //there are no mz pairs in tolerance, but since we have all peaks and 300 default cutoff, we get 3 intensity pairs with 300, 400, 500 intensities being paired zero
            Assert.AreEqual(3, s.intensityPairs.Count);

            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.247).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.867).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));


            //Test what happens when some initial peak intensities are zero.
            primary = new MzSpectrum(new double[] { 1000, 2000, 3000 }, new double[] { 0, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 1000, 2000, 3000, 4000 }, new double[] { 2, 0, 2, 2 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);

            //First, we throw out zero intensity peaks.
            //primary spectrum becomes 2000, 3000 and secondary spectrum becomes 1000, 3000, 4000
            //mz pairs are (1000,1000), (2000,2000) and (3000,3000). Since we are using all peaks, we get an additional pairs with 4000 intensity being paired zero
            Assert.AreEqual(4, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.48).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.319).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.327).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.333).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.333).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.20).Within(0.01));

            //Test what happens when all intensity pairs include 1 zero
            primary = new MzSpectrum(new double[] { 1000, 2000, 3000 }, new double[] { 0, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 4000, 5000 }, new double[] { 2, 0 }, false);

            s = new SpectralSimilarity(primary, secondary, SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true);
            //There are no mz pairs within tolerance. Since we are using all peaks, we get 5 intensity pairs with 1000, 2000, 3000, 4000 and 5000 intensities being paired zero
            Assert.AreEqual(5, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.232).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.395).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));

            //Test alternate constructor
            primary = new MzSpectrum(new double[] { 100, 350, 400 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 350 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, true);
            //mz pairs in tolerance are (350,350). Since we are using all peaks and default mz cut off of 300 we get 1 additional intensity pairs with 400 intensity being paired zero
            Assert.AreEqual(2, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.555).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.374).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.054).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.50).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.667).Within(0.01));

            //Test alternate constructor only library peaks. Since library has one peak, and primary has three peaks, we get only one intensity pair
            primary = new MzSpectrum(new double[] { 350, 400, 500 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 400 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false);
            //mz pairs in tolerance are (400,400). Since we are NOT using all peaks, we get only 1 intensity pair
            Assert.AreEqual(1, s.intensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.667).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.80).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1.0).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.667).Within(0.01));

            //Test cosine similarity when there are no peaks from spectrum one matching spectrum 2
            primary = new MzSpectrum(new double[] { 450, 650, 850 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 400, 600, 800 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false);
            //There are no mz pairs in tolerance. Since we are NOT using all peaks, we get 0 intensity pairs.
            //However, since there are zero pairs, we return as a pair (-1,-1) so we know that and there is no crash.
            Assert.AreEqual(1, s.intensityPairs.Count);
            Assert.AreEqual((-1, -1), s.intensityPairs[0]);
            Assert.IsNull(s.CosineSimilarity());
            Assert.IsNull(s.SpectralContrastAngle());
            Assert.IsNull(s.EuclideanDistance());
            Assert.IsNull(s.BrayCurtis());
            Assert.IsNull(s.PearsonsCorrelation());
            Assert.IsNull(s.DotProduct());

            //Test when intensityPairs is null 
            //Test when all mz are less than 300 in experimental peaks
            primary = new MzSpectrum(new double[] { 100, 150, 200 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 400, 600, 800 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false);
            //There are no mz pairs in tolerance. It doesn't matter but 3 peaks are also below mz300 cut off.
            //Since we are NOT using all peaks, we get 0 intensity pairs.
            //However, since there are zero pairs, we return as a pair (-1,-1) so we know that and there is no crash.
            Assert.AreEqual(1, s.intensityPairs.Count);
            Assert.AreEqual((-1, -1), s.intensityPairs[0]);
            Assert.IsNull(s.CosineSimilarity());
            Assert.IsNull(s.SpectralContrastAngle());
            Assert.IsNull(s.EuclideanDistance());
            Assert.IsNull(s.BrayCurtis());
            Assert.IsNull(s.PearsonsCorrelation());
            Assert.IsNull(s.DotProduct());

            //Test when all mz are less than 300 in theoretical peaks
            primary = new MzSpectrum(new double[] { 150, 250, 350 }, new double[] { 2, 4, 6 }, false);
            secondary = new MzSpectrum(new double[] { 100, 150, 200 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(primary, secondary.XArray, secondary.YArray, SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, ppmTolerance, false);
            //There are no mz pairs in tolerance. It doesn't matter but all 6 peaks are also below mz300 cut off.
            //Since we are NOT using all peaks, we get 0 intensity pairs.
            //However, since there are zero pairs, we return as a pair (-1,-1) so we know that and there is no crash.
            Assert.AreEqual(1, s.intensityPairs.Count);
            Assert.AreEqual((-1, -1), s.intensityPairs[0]);
            Assert.IsNull(s.CosineSimilarity());
            Assert.IsNull(s.SpectralContrastAngle());
            Assert.IsNull(s.EuclideanDistance());
            Assert.IsNull(s.BrayCurtis());
            Assert.IsNull(s.PearsonsCorrelation());
            Assert.IsNull(s.DotProduct());
        }

        [Test]
        public void TestKullbackLeiblerDivergence()
        {
            double ppmTolerance = 10;
            double[] p_XArray = new double[] { 1, 2, 3 };
            double[] p_YArray = new double[] { 9.0/25.0, 12.0/25.0, 4.0/25.0 };
            double[] q_XArray = new double[] { 1, 2, 3 };
            double[] q_YArray = new double[] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

            //Test when using all peaks of primary(experimental) and secondary(theoretical)  spectra(bool allpeaks is true) and mz cut off is 0(no cut off)
            SpectralSimilarity s = new(p_XArray, p_YArray, q_XArray, q_YArray, SpectralSimilarity.SpectrumNormalizationScheme.unnormalized, ppmTolerance, true, 0);
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.EqualTo(0.0853).Within(0.001));

            // ignore negative intensity
            p_XArray = new double[] { 1, 2, 3, 4 };
            p_YArray = new double[] { 9.0 / 25.0, 12.0 / 25.0, 4.0 / 25.0, -1.0 / 25.0 };
            q_XArray = new double[] { 1, 2, 3 };
            q_YArray = new double[] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

            //Test when using all peaks of primary(experimental) and secondary(theoretical)  spectra (bool allpeaks is true) and mz cut off is 0 (no cut off)
            s = new(p_XArray, p_YArray, q_XArray, q_YArray, SpectralSimilarity.SpectrumNormalizationScheme.unnormalized, ppmTolerance, true, 0);
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.EqualTo(0.0853).Within(0.001));
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.EqualTo(s.KullbackLeiblerDivergence_P_Q(correctionConstant: 0)).Within(0.001));

            // ignore negative mz
            p_XArray = new double[] { 1, 2, 3, -4.0 };
            p_YArray = new double[] { 9.0 / 25.0, 12.0 / 25.0, 4.0 / 25.0, 1.0 / 25.0 };
            q_XArray = new double[] { 1, 2, 3 };
            q_YArray = new double[] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

            //Test when using all peaks of primary(experimental) and secondary(theoretical)  spectra (bool allpeaks is true) and mz cut off is 0 (no cut off)
            s = new(p_XArray, p_YArray, q_XArray, q_YArray, SpectralSimilarity.SpectrumNormalizationScheme.unnormalized, ppmTolerance, true, 0);
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.EqualTo(0.0853).Within(0.001));

            // correct for 0 intensity values
            p_XArray = new double[] { 1, 2, 3, 4, 5 };
            p_YArray = new double[] { 0.0, 0.0, 9.0 / 25.0, 12.0 / 25.0, 4.0 / 25.0 };
            q_XArray = new double[] { 1, 2, 3, 4, 5 };
            q_YArray = new double[] { 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 0.0 };

            //Test when using all peaks of primary(experimental) and secondary(theoretical)  spectra (bool allpeaks is true) and mz cut off is 0 (no cut off)
            s = new(p_XArray, p_YArray, q_XArray, q_YArray,
                SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true, 0);
            // With correction, this should increase divergence for missing peaks
            //double? corrected = s.KullbackLeiblerDivergence_P_Q(correctionConstant: 1e-8);
            //double? uncorrected = s.KullbackLeiblerDivergence_P_Q(correctionConstant: 0);
            //double? testCorrect = s.KullbackLeiblerDivergence_P_Q(correctionConstant: 1.0/50.0);
            //double? defaultValue = s.KullbackLeiblerDivergence_P_Q(correctionConstant: 1e-9);
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.GreaterThan(3));
            Assert.That(s.KullbackLeiblerDivergence_P_Q() > s.KullbackLeiblerDivergence_P_Q(correctionConstant: 0));

            // correct for 0 intensity values
            p_XArray = new double[] { 1, 2, 3, 4, 5 };
            p_YArray = new double[] { 0.0, 4.0/25.0, 9.0 / 25.0, 8.0 / 25.0, 4.0 / 25.0 };
            q_XArray = new double[] { 1, 2, 3, 4, 5 };
            q_YArray = new double[] { 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 0.0 };

            //Test when using all peaks of primary(experimental) and secondary(theoretical)  spectra (bool allpeaks is true) and mz cut off is 0 (no cut off)
            s = new(p_XArray, p_YArray, q_XArray, q_YArray,
                SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, ppmTolerance, true, 0);
            // With correction, this should increase divergence for missing peaks
            //corrected = s.KullbackLeiblerDivergence_P_Q(correctionConstant: 1e-8);
            //uncorrected = s.KullbackLeiblerDivergence_P_Q(correctionConstant: 0);
            //testCorrect = s.KullbackLeiblerDivergence_P_Q(correctionConstant: 1.0 / 50.0);
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.GreaterThan(3));
            Assert.That(s.KullbackLeiblerDivergence_P_Q() > s.KullbackLeiblerDivergence_P_Q(correctionConstant: 0));

            // Test for no overlapping peaks
            p_XArray = new double[] { 1, 2, 3, 4 };
            p_YArray = new double[] { 9.0 / 25.0, 12.0 / 25.0, 0.0 / 25.0, 0.0 };
            q_XArray = new double[] { 1, 2, 3, 4 };
            q_YArray = new double[] { 0.0 / 3.0, 0.0 / 25.0, 1.0 / 3.0, 8.0 / 25.0 };

            s = new(p_XArray, p_YArray, q_XArray, q_YArray, SpectralSimilarity.SpectrumNormalizationScheme.unnormalized, ppmTolerance, true, 0);
            // With correction, this should increase divergence for missing peaks
            Assert.That(s.KullbackLeiblerDivergence_P_Q() == null);

        }
    }
}