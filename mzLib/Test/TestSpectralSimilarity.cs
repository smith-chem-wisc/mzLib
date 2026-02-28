using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using System;
using System.Collections.Generic;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

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
            MzSpectrum experimentalSpectrum = new(new double[] { 1, 2, 3, 4, 5 }, new double[] { 2, 4, 6, 8, 10 }, false);
            MzSpectrum theoreticalSpectrum = new(new double[] { 3, 4, 5, 6, 7 }, new double[] { 9, 7, 5, 3, 1 }, false);
            //
            SpectralSimilarity s = new(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, ppmTolerance, true, 0);
            //mz pairs in tolerance are (3,3), (4,4), (5,5). Since we are using all peaks, we get 7 intensity pairs with 1,2,6 and 7 intensities being paired zero
            Assert.AreEqual(7, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.8).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.59).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.70).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.66).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(0.42).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.SearleSimilarity(), Is.EqualTo(2.4391).Within(0.01));
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true, 0);
            Assert.That(s.SpectralEntropy(), Is.EqualTo(0.79).Within(0.01));

            //Test all normalization schemes
            experimentalSpectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 1 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, ppmTolerance, true,0);
            //mz pairs in tolerance are (1,1). Since we are creating pairs for all experimental peaks we get additional pairs for 2 and 3 with zero intensities
            Assert.AreEqual(3, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.37).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.22).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.33).Within(0.01));

            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true, 0);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-.03).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.17).Within(0.01));

            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, ppmTolerance, true, 0);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.41).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.07).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.24).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.90).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.24).Within(0.01));

            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.Unnormalized, ppmTolerance, true, 0);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-6.21).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.29).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(4).Within(0.01));

            //Make sure there is no crash w/ empty array
            experimentalSpectrum = new MzSpectrum(Array.Empty<double>(), new double[] { }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { }, new double[] { }, false);

            Assert.Throws<MzLibException>(() =>
            {
                s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true, 0);
            }, "Empty YArray in spectrum.");

            //We should have any zero intensity YArrays but just to be sure
            experimentalSpectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 1 }, new double[] { 0 }, false);

            Assert.Throws<MzLibException>(() =>
            {
                s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true, 0);
            }, "Spectrum has no intensity.");

            //What happens when all intensity pairs include a zero
            experimentalSpectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 4 }, new double[] { 2 }, false);

            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true, 0);
            //There are no mz pairs in tolerance. But we automatically get a pair for 4 because it is a theoretical mz. And we get pairs for 1,2,3 with zero intensities because we are using all experimental peaks.
            Assert.AreEqual(4, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.18).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.77).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));

            experimentalSpectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 0, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 2, 0, 2, 2 }, false);

            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true, 0);
            //mz pairs in tolerance are (1,1), (2,2), (3,3). We also automatically get a pair for 4 because it is a theoretical mz.
            Assert.AreEqual(4, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.48).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.32).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.33).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.20).Within(0.01));

            //Test what happens when all intensity pairs include 1 zero
            experimentalSpectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 0, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 4, 5 }, new double[] { 2, 0 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true, 0);
            //There are no mz pairs in tolerance. We eliminate the peak at 3 because it has zero intensity. We also remove the theoretical pair at 5 because it has zero intensity.
            //We get 3 pairs for 2, 3 and 4.
            Assert.AreEqual(3, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.23).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.945).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));

            //explore bounds of binary search
            experimentalSpectrum = new MzSpectrum(new double[] { 1, 2, 3, 4 }, new double[] { 1, 2, 3, 4 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 1.000011, 1.99997, 3.000031, 3.99995 }, new double[] { 1, 2, 3, 4 }, false);

            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true, 0);
            //The ppm difference between the 4 closest pairs are 11, 15,10.3 and 12.5. These are all beyond the 10ppm tolerance that is allowed. Therefore we get 8 intensity pairs with all intensities being paired zero
            Assert.AreEqual(8, s.IntensityPairs.Count);

            //Test alternate constructor
            experimentalSpectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 1 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, ppmTolerance, true,0);
            //mz pairs in tolerance are (1,1). Since we are using all peaks, we get 3 intensity pairs with 2 and 3 intensities being paired z
            Assert.AreEqual(3, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.37).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.22).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.33).Within(0.01));

            //Test alternate constructor only library peaks. Since library has one peak, and primary has three peaks, we get only one intensity pair
            experimentalSpectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 1 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, ppmTolerance, false,0);
            //mz pairs in tolerance are (1,1). Since we are NOT using all peaks, we get only 1 intensity pair
            Assert.AreEqual(1, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.50).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1.0).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.33).Within(0.01));

            //Test cosine similarity when there are no peaks from spectrum one matching spectrum 2
            experimentalSpectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 4,6,8 }, new double[] { 2,4,6 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, ppmTolerance, false, 0);
            //There are no mz pairs in tolerance. But we keep all three theoretical peaks so we get three intesity pairs with all intensities being paired zero
            Assert.AreEqual(3, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));

            //Test SearleSimilarity with both spectra are identical
            experimentalSpectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, ppmTolerance, false, 0);
            //there are 3 mz pairs in tolerance. 
            Assert.AreEqual(s.SearleSimilarity(), double.MaxValue);
        }

        [Test]
        public void TestAllSpectrumSimilaritiesWithDefaultedMzFilter()
        {
            double ppmTolerance = 10;
            MzSpectrum experimentalSpectrum = new MzSpectrum(new double[] { 100, 200, 300, 400, 500 ,600,700}, new double[] { 2, 4, 6, 8, 10,12,14 }, false);
            MzSpectrum theoreticalSpectrum = new MzSpectrum(new double[] { 200, 300, 500, 600 ,800}, new double[] { 9, 7, 5, 3, 1 }, false);
            //Test when using all peaks of primary(experimental) and secondary(theoretical)  spectra (bool all experimental peaks is true) and mz cut off is 0 (no cut off)
            //we are keeping everything and there are 8 mz peaks with more than zero intnsity so we get 8 intensity pairs
            SpectralSimilarity s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, ppmTolerance, true, 0);
            Assert.AreEqual(8, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.68).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.48).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.65).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.56).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(0.02).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.13).Within(0.01));

            //Test when using all peaks of primary(experimental) and secondary(theoretical) spectra (bool all experimental peaks is true) and mz cut off is 300 (default cut off)
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, ppmTolerance, true);
            //similary to above but we remove peaks below 300. That leaves 6 intensity pairs.
            Assert.AreEqual(6, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.70).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.49).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.61).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.58).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.04).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.16).Within(0.01));

            //Test when not using all peaks of primary(experimental) spectra (bool all experimental peaks is false) and mz cut off is is 0 (no cut off)
            //experimental xArray 100, 200, 300, 400, 500, 600, 700 and theoretical xArray 200, 300, 500, 600, 800. So, 200, 300, 500, 600 are common. We keep 800 from the theoretical spectrum
            //because the default is to keep all theoretical peaks.
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, ppmTolerance, false, 0);
            Assert.AreEqual(5, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.903).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.718).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.76).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.712).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(0.475).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.126).Within(0.01));

            //Test when not using all peaks of primary(experimental) spectra (bool all experimental peaks is false) and mz cut off is is 300 (default cut off)
            //primary xArray 100, 200, 300, 400, 500, 600, 700 and secondary xArray 200, 300, 500, 600, 800. So, 200, 300, 500, 600 are common. But with 300 cut off, only 300, 500, 600 are common
            //we keep 800 because it is a theoretical peak.
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, ppmTolerance, false);
            Assert.AreEqual(4, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.924).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.75).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.751).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.734).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(0.681).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.164).Within(0.01));

            //Test all different similarity calculations
            experimentalSpectrum = new MzSpectrum(new double[] { 100, 200, 300, 400, 500 }, new double[] { 2, 4, 6, 8, 10 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 300, 400, 500, 600, 700 }, new double[] { 9, 7, 5, 3, 1 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, ppmTolerance, false);
            //all experimental peaks is false and there is no mz cut off.
            //Therefore we get a intensity pair for each theoretical mz for a total of 5.
            Assert.AreEqual(5, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.894).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.704).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.736).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.743).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(0.812).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.244).Within(0.01));

            //Test all normalization schemes
            experimentalSpectrum = new MzSpectrum(new double[] { 1000, 2000, 3000 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 1000 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, ppmTolerance, true);
            //there is one mz pair in tolerance (1000,1000). Since we are using all experimental peaks, we get additional pairs for 2000 and 3000 with zero intensities
            Assert.AreEqual(3, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.267).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.172).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.374).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.222).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.866).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.333).Within(0.01));

            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-.03).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.17).Within(0.01));

            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, ppmTolerance, true);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.41).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.07).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.24).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.90).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.24).Within(0.01));

            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.Unnormalized, ppmTolerance, true);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.27).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.17).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-6.21).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.29).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.87).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(4).Within(0.01));


            //What happens when all intensity pairs include a zero
            experimentalSpectrum = new MzSpectrum(new double[] { 100, 200, 300 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 400, 500 }, new double[] { 2, 4 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true);
            //There are no mz pairs in tolerence. But, we are keeping all experimental peaks. With the theoretical peaks with mz >= 300. That gives us 3 intensity pairs
            Assert.AreEqual(3, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.247).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.866).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));

            experimentalSpectrum = new MzSpectrum(new double[] { 1000, 2000, 3000 }, new double[] { 0, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 1000, 2000, 3000, 4000 }, new double[] { 2, 0, 2, 2 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true);
            //There are 3 mz pairs in tolerance plus an extra from theoretical totaling 4 intensity pairs
            Assert.AreEqual(4, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.48).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.32).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.33).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.33).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.20).Within(0.01));

            //Test what happens when all intensity pairs include 1 zero
            experimentalSpectrum = new MzSpectrum(new double[] { 1000, 2000, 3000 }, new double[] { 0, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 4000, 5000 }, new double[] { 2, 0 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true);
            //We toss the experimental peak at 1000 and the theoretical peak at 5000 because they have zero intensity.
            //there are no mz pairs in tolerance. But we are keeping all experimental and theoretical peaks. We get 3 intensity pairs for 2000, 3000 and 4000
            Assert.AreEqual(3, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.23).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-0.949).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0).Within(0.01));

            //Test alternate constructor
            experimentalSpectrum = new MzSpectrum(new double[] { 100, 350, 400 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 350 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, ppmTolerance, true);
            //there is one mz pair with mz > 300 (350,350). We keep all experimental peaks > 300. That give one more intensity pair for 400.
            Assert.AreEqual(2, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0.55).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0.37).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(-0.05).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.50).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.67).Within(0.01));

            //Test alternate constructor only library peaks. Since library has one peak, and primary has three peaks, we get only one intensity pair
            experimentalSpectrum = new MzSpectrum(new double[] { 350, 400, 500 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 400 }, new double[] { 2 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, ppmTolerance, false);
            Assert.AreEqual(1, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(1.0).Within(0.01));
            Assert.That(s.EuclideanDistance(), Is.EqualTo(0.67).Within(0.01));
            Assert.That(s.BrayCurtis(), Is.EqualTo(0.80).Within(0.01));
            Assert.That(s.PearsonsCorrelation(), Is.EqualTo(-1.0).Within(0.01));
            Assert.That(s.DotProduct(), Is.EqualTo(0.67).Within(0.01));

            //Test cosine similarity when there are no peaks from spectrum one matching spectrum 2
            experimentalSpectrum = new MzSpectrum(new double[] { 450, 650, 850 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 400, 600, 800 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, ppmTolerance, false);
            Assert.AreEqual(3, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.EqualTo(0).Within(0.01));
            Assert.That(s.SpectralContrastAngle(), Is.EqualTo(0).Within(0.01));

            //Test when intensityPairs is null 
            //Test when all mz are less than 300 in experimental peaks
            experimentalSpectrum = new MzSpectrum(new double[] { 100, 150, 200 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 400, 600, 800 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, ppmTolerance, false);
            Assert.AreEqual(1, s.IntensityPairs.Count);
            Assert.AreEqual(new List<(double, double)> { (-1, -1) }, s.IntensityPairs);
            Assert.That(s.CosineSimilarity(), Is.Null);
            Assert.That(s.SpectralContrastAngle(), Is.Null);
            Assert.That(s.EuclideanDistance(), Is.Null);
            Assert.That(s.BrayCurtis(), Is.Null);
            Assert.That(s.PearsonsCorrelation(), Is.Null);
            Assert.That(s.DotProduct(), Is.Null);

            //Test when all mz are less than 300 in theoretical peaks
            experimentalSpectrum = new MzSpectrum(new double[] { 150, 250, 350 }, new double[] { 2, 4, 6 }, false);
            theoreticalSpectrum = new MzSpectrum(new double[] { 100, 150, 200 }, new double[] { 2, 4, 6 }, false);
            s = new SpectralSimilarity(experimentalSpectrum, theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, ppmTolerance, false);
            Assert.AreEqual(1, s.IntensityPairs.Count);
            Assert.AreEqual(new List<(double, double)> { (-1, -1) }, s.IntensityPairs);
            Assert.That(s.CosineSimilarity(), Is.Null);
            Assert.That(s.SpectralContrastAngle(), Is.Null);
            Assert.That(s.EuclideanDistance(), Is.Null);
            Assert.That(s.BrayCurtis(), Is.Null);
            Assert.That(s.PearsonsCorrelation(), Is.Null);
            Assert.That(s.DotProduct(), Is.Null);
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
            SpectralSimilarity s = new(p_XArray, p_YArray, q_XArray, q_YArray, SpectralSimilarity.SpectrumNormalizationScheme.Unnormalized, ppmTolerance, true, 0);
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.EqualTo(0.0853).Within(0.001));

            // ignore negative intensity
            p_XArray = new double[] { 1, 2, 3, 4 };
            p_YArray = new double[] { 9.0 / 25.0, 12.0 / 25.0, 4.0 / 25.0, -1.0 / 25.0 };
            q_XArray = new double[] { 1, 2, 3 };
            q_YArray = new double[] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

            //Test when using all peaks of primary(experimental) and secondary(theoretical)  spectra (bool allpeaks is true) and mz cut off is 0 (no cut off)
            s = new(p_XArray, p_YArray, q_XArray, q_YArray, SpectralSimilarity.SpectrumNormalizationScheme.Unnormalized, ppmTolerance, true, 0);
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.EqualTo(0.0853).Within(0.001));
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.EqualTo(s.KullbackLeiblerDivergence_P_Q(correctionConstant: 0)).Within(0.001));

            // ignore negative mz
            p_XArray = new double[] { 1, 2, 3, -4.0 };
            p_YArray = new double[] { 9.0 / 25.0, 12.0 / 25.0, 4.0 / 25.0, 1.0 / 25.0 };
            q_XArray = new double[] { 1, 2, 3 };
            q_YArray = new double[] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

            //Test when using all peaks of primary(experimental) and secondary(theoretical)  spectra (bool allpeaks is true) and mz cut off is 0 (no cut off)
            s = new(p_XArray, p_YArray, q_XArray, q_YArray, SpectralSimilarity.SpectrumNormalizationScheme.Unnormalized, ppmTolerance, true, 0);
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.EqualTo(0.0853).Within(0.001));

            // correct for 0 intensity values
            p_XArray = new double[] { 3, 4, 5 };
            p_YArray = new double[] {9.0 / 25.0, 12.0 / 25.0, 4.0 / 25.0 };
            q_XArray = new double[] { 1, 2, 3, 4, 5 };
            q_YArray = new double[] { 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 0.0 };

            //Test when using all peaks of primary(experimental) and secondary(theoretical)  spectra (bool all experimental peaks is true) and mz cut off is 0 (no cut off)
            s = new(p_XArray, p_YArray, q_XArray, q_YArray,
                SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true, 0);
            // With correction, this should increase divergence for missing peaks
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.EqualTo(3.467).Within(0.01));
            Assert.That(s.KullbackLeiblerDivergence_P_Q() > s.KullbackLeiblerDivergence_P_Q(correctionConstant: 0));

            // correct for 0 intensity values
            p_XArray = new double[] { 1, 2, 3, 4, 5 };
            p_YArray = new double[] { 0.0, 4.0/25.0, 9.0 / 25.0, 8.0 / 25.0, 4.0 / 25.0 };
            q_XArray = new double[] { 1, 2, 3, 4, 5 };
            q_YArray = new double[] { 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 0.0 };

            //Test when using all peaks of primary(experimental) and secondary(theoretical)  spectra (bool allpeaks is true) and mz cut off is 0 (no cut off)
            s = new(p_XArray, p_YArray, q_XArray, q_YArray,
                SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum, ppmTolerance, true, 0);
            // With correction, this should increase divergence for missing peaks
            Assert.That(s.KullbackLeiblerDivergence_P_Q(), Is.GreaterThan(3));
            Assert.That(s.KullbackLeiblerDivergence_P_Q() > s.KullbackLeiblerDivergence_P_Q(correctionConstant: 0));

            // Test for no overlapping peaks
            p_XArray = new double[] { 1, 2, 3, 4 };
            p_YArray = new double[] { 9.0 / 25.0, 12.0 / 25.0, 0.0 / 25.0, 0.0 };
            q_XArray = new double[] { 1, 2, 3, 4 };
            q_YArray = new double[] { 0.0 / 3.0, 0.0 / 25.0, 1.0 / 3.0, 8.0 / 25.0 };

            s = new(p_XArray, p_YArray, q_XArray, q_YArray, SpectralSimilarity.SpectrumNormalizationScheme.Unnormalized, ppmTolerance, true, 0);
            // With correction, this should increase divergence for missing peaks
            Assert.That(s.KullbackLeiblerDivergence_P_Q() == null);

        }
        [Test]
        public void TestSpectralSimilarityWithInternalFragmentIons()
        {
            double ppmTolerance = 10;

            // Create a spectrum with both standard and internal fragment ions
            // Simulating what would come from a library spectrum with internal ions
            // Standard ions: b1, b2, y1, y2
            // Internal ions: bIb[2-4], bIb[3-5]
            double[] spectrum1Mz = new double[] { 100.0, 200.0, 300.0, 400.0, 250.0, 350.0 };
            double[] spectrum1Intensity = new double[] { 1000, 2000, 1500, 2500, 800, 1200 };

            // Second spectrum with similar peaks (simulating another acquisition of the same peptide)
            double[] spectrum2Mz = new double[] { 100.0, 200.0, 300.0, 400.0, 250.0, 350.0 };
            double[] spectrum2Intensity = new double[] { 1100, 1900, 1600, 2400, 750, 1300 };

            // Test that similarity calculations work with internal ion m/z values
            SpectralSimilarity s = new SpectralSimilarity(
                spectrum1Mz, spectrum1Intensity,
                spectrum2Mz, spectrum2Intensity,
                SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
                ppmTolerance, true, 0);

            // All 6 peaks should match
            Assert.AreEqual(6, s.IntensityPairs.Count);

            // High similarity expected since spectra are nearly identical
            Assert.That(s.CosineSimilarity(), Is.GreaterThan(0.95));
            Assert.That(s.SpectralContrastAngle(), Is.GreaterThan(0.85));
            Assert.That(s.PearsonsCorrelation(), Is.GreaterThan(0.95));

            // Test with partial overlap - some internal ions missing
            double[] spectrum3Mz = new double[] { 100.0, 200.0, 300.0, 400.0 }; // Missing internal ions
            double[] spectrum3Intensity = new double[] { 1000, 2000, 1500, 2500 };

            s = new SpectralSimilarity(
                spectrum1Mz, spectrum1Intensity,
                spectrum3Mz, spectrum3Intensity,
                SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
                ppmTolerance, true, 0);

            // 6 peaks total (4 matching + 2 unmatched internal ions paired with zero)
            Assert.AreEqual(6, s.IntensityPairs.Count);
            // Similarity should be lower due to missing internal ions
            Assert.That(s.CosineSimilarity(), Is.LessThan(0.95));
            Assert.That(s.CosineSimilarity(), Is.GreaterThan(0.7));

            // Test using only library peaks (allPeaks = false)
            s = new SpectralSimilarity(
                spectrum1Mz, spectrum1Intensity,
                spectrum3Mz, spectrum3Intensity,
                SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
                ppmTolerance, false, 0);

            // Only 4 theoretical peaks considered
            Assert.AreEqual(4, s.IntensityPairs.Count);
            // Higher similarity since we're only comparing matched peaks
            Assert.That(s.CosineSimilarity(), Is.EqualTo(1.0).Within(0.01));
        }

        [Test]
        public void TestSpectralSimilarityWithInternalIonsFromLibrarySpectrum()
        {
            // Create MatchedFragmentIons including internal fragments
            // Standard b-ion
            Product bIon = new Product(ProductType.b, FragmentationTerminus.N, 100.0, 2, 2, 0);
            // Standard y-ion
            Product yIon = new Product(ProductType.y, FragmentationTerminus.C, 200.0, 3, 3, 0);
            // Internal fragment ion: bIb[2-5]
            Product internalIon1 = new Product(
                ProductType.b, FragmentationTerminus.None, 300.0,
                fragmentNumber: 2, residuePosition: 2, neutralLoss: 0,
                secondaryProductType: ProductType.b, secondaryFragmentNumber: 5);
            // Internal fragment ion: bIy[3-6]
            Product internalIon2 = new Product(
                ProductType.b, FragmentationTerminus.None, 400.0,
                fragmentNumber: 3, residuePosition: 3, neutralLoss: 0,
                secondaryProductType: ProductType.y, secondaryFragmentNumber: 6);

            // Verify internal ion properties
            Assert.That(bIon.IsInternalFragment, Is.False);
            Assert.That(yIon.IsInternalFragment, Is.False);
            Assert.That(internalIon1.IsInternalFragment, Is.True);
            Assert.That(internalIon2.IsInternalFragment, Is.True);

            // Create matched fragment ions
            var matchedIons1 = new List<MatchedFragmentIon>
    {
        new MatchedFragmentIon(bIon, 100.0, 1000, 1),
        new MatchedFragmentIon(yIon, 200.0, 2000, 1),
        new MatchedFragmentIon(internalIon1, 300.0, 1500, 1),
        new MatchedFragmentIon(internalIon2, 400.0, 2500, 1)
    };

            var matchedIons2 = new List<MatchedFragmentIon>
    {
        new MatchedFragmentIon(bIon, 100.0, 1100, 1),
        new MatchedFragmentIon(yIon, 200.0, 1900, 1),
        new MatchedFragmentIon(internalIon1, 300.0, 1600, 1),
        new MatchedFragmentIon(internalIon2, 400.0, 2400, 1)
    };

            // Create library spectra
            var librarySpectrum1 = new LibrarySpectrum("PEPTIDE", 500.0, 2, matchedIons1, 10.0);
            var librarySpectrum2 = new LibrarySpectrum("PEPTIDE", 500.0, 2, matchedIons2, 10.0);

            // Verify internal ions are in the spectrum
            Assert.That(librarySpectrum1.MatchedFragmentIons.Count(m => m.NeutralTheoreticalProduct.IsInternalFragment), Is.EqualTo(2));
            Assert.That(librarySpectrum2.MatchedFragmentIons.Count(m => m.NeutralTheoreticalProduct.IsInternalFragment), Is.EqualTo(2));

            // Test spectral similarity using the library spectrum's m/z and intensity arrays
            SpectralSimilarity s = new SpectralSimilarity(
                librarySpectrum1.XArray, librarySpectrum1.YArray,
                librarySpectrum2.XArray, librarySpectrum2.YArray,
                SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
                toleranceInPpm: 10, allPeaks: true, filterOutBelowThisMz: 0);

            Assert.AreEqual(4, s.IntensityPairs.Count);
            Assert.That(s.CosineSimilarity(), Is.GreaterThan(0.95));
            Assert.That(s.SpectralContrastAngle(), Is.GreaterThan(0.85));

            // Test using the LibrarySpectrum's built-in method
            string spectralAngle = librarySpectrum1.CalculateSpectralAngleOnTheFly(matchedIons2);
            Assert.That(spectralAngle, Is.Not.EqualTo("N/A"));
            double parsedAngle = double.Parse(spectralAngle);
            Assert.That(parsedAngle, Is.GreaterThan(0.85));
        }

        [Test]
        public void TestSpectralEntropyWithInternalFragmentIons()
        {
            double ppmTolerance = 10;

            // Create spectra with internal ions for spectral entropy calculation
            double[] spectrum1Mz = new double[] { 100.0, 200.0, 300.0, 400.0, 250.0, 350.0 };
            double[] spectrum1Intensity = new double[] { 0.2, 0.3, 0.15, 0.2, 0.08, 0.07 }; // Normalized

            double[] spectrum2Mz = new double[] { 100.0, 200.0, 300.0, 400.0, 250.0, 350.0 };
            double[] spectrum2Intensity = new double[] { 0.2, 0.3, 0.15, 0.2, 0.08, 0.07 };

            SpectralSimilarity s = new SpectralSimilarity(
                spectrum1Mz, spectrum1Intensity,
                spectrum2Mz, spectrum2Intensity,
                SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum,
                ppmTolerance, true, 0);

            // Identical spectra should have spectral entropy = 1.0
            Assert.That(s.SpectralEntropy(), Is.EqualTo(1.0).Within(0.01));

            // Test with different intensity distributions (internal ions have different intensities)
            double[] spectrum3Intensity = new double[] { 0.1, 0.4, 0.1, 0.25, 0.1, 0.05 };

            s = new SpectralSimilarity(
                spectrum1Mz, spectrum1Intensity,
                spectrum2Mz, spectrum3Intensity,
                SpectralSimilarity.SpectrumNormalizationScheme.SpectrumSum,
                ppmTolerance, true, 0);

            // Spectral entropy should be less than 1.0 for different spectra
            Assert.That(s.SpectralEntropy(), Is.LessThan(1.0));
            Assert.That(s.SpectralEntropy(), Is.GreaterThan(0.8));
        }
    }
}