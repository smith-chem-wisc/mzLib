using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MathNet.Numerics.Distributions;
using Microsoft.VisualBasic.FileIO;
using MzLibSpectralAveraging;
using NUnit.Framework;

namespace Test.AveragingTests
{
    [ExcludeFromCodeCoverage]
    public class NoiseEstimatorTests
    {
        [Test]
        [TestCase(WaveletType.Haar, 
            new double[]{0.5d, 0.5d}, 
            new double[] {0.5d,-0.5d})]
        [TestCase(WaveletType.Db4, new[] {
                0.1629,    
                0.5055, 
                0.4461,
                -0.0198,
                -0.1323,
                0.0218,
                0.0233,
                -0.0075
            }, 
            new[] 
            {
                -0.0075,
                -0.0233,
                0.0218,
                0.1323,
                -0.0198, 
                -0.4461,
                0.5055,  
                -0.1629

            }
            )]
        // This test makes sure that the wavelets coefficients and scales are correctly 
        // calculated. 
        public void TestCreateWavelet(WaveletType type, 
            double[] expectedScaling, 
            double[] expectedWavelet)
        {
            var wflt = new WaveletFilter(); 
            wflt.CreateFiltersFromCoeffs(type);
            Assert.That(wflt.WaveletCoefficients, Is.EqualTo(expectedWavelet).Within(0.01));
            Assert.That(wflt.ScalingCoefficients, Is.EqualTo(expectedScaling).Within(0.01));
        }
        
        [Test]
        [TestCase(0d, 0.5)]
        [TestCase(1d, 0.75d)]
        [TestCase(2d, 0.5d)]
        public void TestSimulatedNoiseSineWave(double mean, double stdev)
        {
            Normal normalDist = new Normal(mean, stdev);
            
            // generate sine wave
            double[] sineWave = Enumerable.Range(0, 1000)
                .Select(i => Math.Sin((double)i/(2*Math.PI)))
                .Select(i => 
                    5d*i + Normal.Sample(normalDist.RandomSource, normalDist.Mean, normalDist.StdDev))
                .ToArray();
            double noiseEstimate = BasicStatistics.CalculateStandardDeviation(sineWave.ToList());
            bool success = MRSNoiseEstimator.MRSNoiseEstimation(sineWave, epsilon: 0.1, 
                out double mrsEstimate, 
                waveletType:WaveletType.Db4);
            Console.WriteLine("mrsEstimate: {0}; stddev estimate: {1}", mrsEstimate, noiseEstimate);

            // Explanation of what is tested: 
            /*
             * A sine wave oscillates around zero. The standard deviation is going to be a
             * good enough approximation of the mean of the noise because the amplitude of the sine wave is going
             * to be the primary driver of the standard deviation.
             *
             * To check the results against the standard deviation, we are going to check that the
             * MRS noise estimate 1) succeeds, 2) is less than the standard deviation, and 3)
             * the values given are no greater than 0.1 apart. 
             */
            Assert.True(success);
            Assert.That(mrsEstimate, Is.EqualTo(noiseEstimate).Within(0.1));
        }

        [Test]
        [TestCase(500d, 10, WaveletType.Haar)]
        public void TestSimulatedNoiseGaussianPeak(double mean, double stddev, WaveletType type)
        {
            double frontTerm = 1 / (stddev * Math.Sqrt(2 * Math.PI));
            
            Normal normalDist = new Normal(20, 1);
            
            double[] gaussian = Enumerable.Range(0, 1000)
                .Select(i => 50*frontTerm * Math.Exp(-0.5 * (i - mean) * (i - mean) / (stddev * stddev)))
                .Select(i =>
                    100*i + Normal.Sample(normalDist.RandomSource, normalDist.Mean, normalDist.StdDev))
                .ToArray();
            double noiseEstimateStddev = BasicStatistics.CalculateStandardDeviation(gaussian);
            bool estimateSuccess = MRSNoiseEstimator.MRSNoiseEstimation(gaussian, epsilon: 0.001,
                out double mrsNoiseEstimate,
                waveletType:type);
            Console.WriteLine("Stddev: {0}\nMrs: {1}", noiseEstimateStddev, mrsNoiseEstimate);
            /*
             * Explanation of test:
             * Unlike the sine wave, this single has two components: the first is a noise baseline generated
             * from random Gaussian noise. The second is the gaussian peak at x = 500 that has a standard deviation of 10.
             * This means that the peak is very sharp relative to the signal. This is precisely the situation we are
             * looking at in mass spectrometry data: a sharp peak and a relatively flat noise baseline. I have added
             * a noise to the signal. This noise is Gaussian and has a mean of 20, with a standard deviation of 20.
             *
             * When we run the standard deviation, we get a value of 26, which is a 30% error. However, we get a value
             * 20.7 for the MRS noise estimate, which is only a 3.5% error, and we must also remember that the Gaussian noise
             * values have a standard deviation of 1, so likely, this is as good of an error estimate as we could possibly get.
             *
             * We assert that the MRS noise estimate is lower than the standard deviation noise estimate.
             *
             * DB4 wavelet gives a higher standard deviation, probably because the because width of the wavelet is as longer
             * or longer than the width of peak. So in this case, we should probably not use the DB4 wavelet. 
             */
            Assert.True(estimateSuccess);
            Assert.That(mrsNoiseEstimate, Is.LessThan(noiseEstimateStddev));
            Assert.That(mrsNoiseEstimate, Is.EqualTo(20.7).Within(0.25));
            Assert.That(noiseEstimateStddev, Is.EqualTo(26.5).Within(0.5));
        }

        [Test]
        [TestCase(WaveletType.Haar)]
        [TestCase(WaveletType.Db4)]
        public void TestRealData(WaveletType waveType)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "AveragingTestData",
                "RealScanCarbonicAnhydrase.csv");  
            List<double> mzVals = new();
            List<double> intensityVals = new();
            using (TextFieldParser csvParser = new TextFieldParser(path))
            {
                csvParser.CommentTokens = new[] { "#" };
                csvParser.SetDelimiters(new string[] { "," });
                csvParser.HasFieldsEnclosedInQuotes = false;
                while (!csvParser.EndOfData)
                {
                    string[] fields = csvParser.ReadFields();
                    mzVals.Add(Convert.ToDouble(fields[0]));
                    intensityVals.Add(Convert.ToDouble(fields[1]));
                }
            }
            WaveletFilter wflt = new WaveletFilter();
            wflt.CreateFiltersFromCoeffs(WaveletType.Db4);
            double[] signal = intensityVals.ToArray();
            double noiseStddev = BasicStatistics.CalculateStandardDeviation(signal); 
            bool noiseEstimationSuccess = MRSNoiseEstimator.MRSNoiseEstimation(signal,
                0.1, out double noiseEstimate, waveletType:waveType);
            Console.WriteLine("Standard Deviation: {0}; MRS: {1}", noiseStddev, noiseEstimate);
            /*
             * For the real data test, we assert that the MRS noise estimate is lower than the
             * standard deviation noise estimate.
             *
             * The haar wavelet gives a higher estimate than the db4 wavelet. But the haar wavelet
             * is significantly faster than the haar wavelet. Interesting. 
             */
            Assert.True(noiseEstimationSuccess);
            Assert.That(noiseEstimate, Is.LessThan(noiseStddev));
        }
    }
}
