using System;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MathNet.Numerics.Distributions;
using MzLibSpectralAveraging;
using NUnit.Framework;

namespace Test.AveragingTests
{
    [ExcludeFromCodeCoverage]
    public class TestScanCombination
    {
        private double[][] xArrays; 
        private double[][] yArrays;
        private double[] tics;
        private double[] averagedArrayFin; 

        [OneTimeSetUp]
        public void OneTimeSetup()
        {
            #region Old, Easy test
            //// make ten arrays of the following: 
            //// make x array from 500 to 2000 m/z, spaced at 0.001 m/z apart
            //// make y array random 

            //double[][] xarrays = new double[10][];
            //double[][] yarrays = new double[10][];
            //double numberSteps = (2000 - 500) / 0.01;
            //double[] xarray = new double[(int)numberSteps];
            //for (int i = 0; i < numberSteps; i++)
            //{
            //    xarray[i] = 500 + i * 0.01;
            //}

            //// create random values with seed starting at 1551. 
            //int initialSeed = 1551; 

            //for (int i = 0; i < xarrays.GetLength(0); i++)
            //{
            //    xarrays[i] = xarray; 
            //    yarrays[i] = new double[(int)numberSteps];

            //    Random rnd = new Random(initialSeed + i);
            //    for (int j = 0; j < yarrays[i].Length; j++)
            //    {
            //        yarrays[i][j] = rnd.NextDouble() * 10000; 
            //    }
            //}
            //_xarrays = xarrays;
            //_yarrays = yarrays;
            //_tics = new double[xarrays.GetLength(0)];
            //for (int i = 0; i < _tics.Length; i++)
            //{
            //    _tics[i] = yarrays[i].Sum(); 
            //}
            #endregion

            // this tests more mass spectrometry-like data. 
            double stddev = 50;
            double mean = 500; 
            double frontTerm = 1 / (stddev * Math.Sqrt(2 * Math.PI));

            Normal normalDist = new Normal(20, 10);
            // normal distribution to shift the m/z values by
            Normal mzShifts = new Normal(0, 0.001);

            xArrays = new double[25][];
            yArrays = new double[25][];
            tics = new double[25];
            for (int m = 0; m < xArrays.GetLength(0); m++)
            {
                double[] gaussian = Enumerable.Range(1, 1000)
                    .Select(i => 100 * frontTerm * Math.Exp(-0.5 * (i - mean) * (i - mean) / (stddev * stddev)))
                    .Select(i =>
                        100 * i + Normal.Sample(normalDist.RandomSource, normalDist.Mean, normalDist.StdDev))
                    .ToArray();
                double[] xAxis = Enumerable.Range(1, 1000)
                    .Select(i => i + Normal.Sample(mzShifts.RandomSource, mzShifts.Mean, mzShifts.StdDev))
                    .ToArray();

                xArrays[m] = new double[gaussian.Length];
                yArrays[m] = new double[gaussian.Length];

                yArrays[m] = gaussian;
                xArrays[m] = xAxis;
                tics[m] = gaussian.Sum();
            }

            double[] averagedArray = new double[yArrays[0].Length];

            for (int i = 0; i < yArrays.Length; i++)
            {
                double tic = yArrays[i].Sum();
                for (int j = 0; j < yArrays[0].Length; j++)
                {
                    averagedArray[j] += yArrays[i][j] / tic;
                }
            }
            averagedArrayFin = averagedArray.Select(i => i / yArrays.GetLength(0))
                .ToArray();
        }

        [Test, 
        TestCase(RejectionType.SigmaClipping), 
        TestCase(RejectionType.WinsorizedSigmaClipping), 
        TestCase(RejectionType.AveragedSigmaClipping), 
        TestCase(RejectionType.MinMaxClipping), 
        TestCase(RejectionType.PercentileClipping)]
        public void TestCombination(RejectionType rejection)
        {
            using (StreamWriter sr = new("noisyOutput.txt"))
            {
                for (int i = 0; i < xArrays[0].Length; i++)
                {
                    sr.WriteLine("{0},{1}", xArrays[0][i], yArrays[0][i]);
                }
                sr.Flush();
            }

            Stopwatch sw = Stopwatch.StartNew();
            SpectralAveragingOptions options = new SpectralAveragingOptions();
            options.SetDefaultValues();
            options.BinSize = 1.0; 
            options.SpectrumMergingType = SpectrumMergingType.MrsNoiseEstimate;
            options.RejectionType = rejection; 
            // TODO: Generate multiple binned spectra from a set of xArrays given the 
            // number of spectra. 
            double[][] results = SpectralMerging.CombineSpectra(xArrays, yArrays,
                tics, 25, options);
            sw.Stop();
            Console.WriteLine(sw.ElapsedMilliseconds);

            using (StreamWriter sr = new("averagedOutputs.txt"))
            {
                for (int i = 0; i < results[0].Length; i++)
                {
                    sr.WriteLine("{0},{1}", results[0][i], results[1][i]);
                }
                sr.Flush();
            }

            MRSNoiseEstimator.MRSNoiseEstimation(results[1][..300], 0.01, out double noiseEstimateAveraged);
            // normalize reference spectra to tics, then estimate noise 
            double ticReference = yArrays[0].Sum();
            double[] referenceSpectrum = yArrays[0].Select(i => i / ticReference).ToArray();
            MRSNoiseEstimator.MRSNoiseEstimation(referenceSpectrum[..300], 0.01, out double noiseEstimateReference); 
            // get the effective noise reduction function (ENR)
            // calculate it for parts that you know are dominated by noise
            // ENR = sigma_final / (scale_ref / scale_final * sigma_ref)

            // TODO: WEIGHTS ARE NOT BEING CALCULATED OR USED CORRECTLY

            MRSNoiseEstimator.MRSNoiseEstimation(averagedArrayFin, 0.01, out double estimateSimpleAverage); 
            Console.WriteLine("Simple Average: {0}; Integration {1}; Reference {2}", 
                estimateSimpleAverage, noiseEstimateAveraged, noiseEstimateReference);
        }
        private double MedianAbsoluteDeviationFromMedian(double[] array)
        {
            double arrayMedian = BasicStatistics.CalculateMedian(array);
            double[] results = new double[array.Length];
            for (int j = 0; j < array.Length; j++)
            {
                results[j] = Math.Abs(array[j] - arrayMedian);
            }

            return BasicStatistics.CalculateMedian(results);
        }

        private double BiweightMidvariance(double[] array)
        {
            double[] y_i = new double[array.Length];
            double[] a_i = new double[array.Length];
            double MAD_X = MedianAbsoluteDeviationFromMedian(array);
            double median = BasicStatistics.CalculateMedian(array);
            for (int i = 0; i < y_i.Length; i++)
            {
                y_i[i] = (array[i] - median) / (9d * MAD_X);
                if (y_i[i] < 1d)
                {
                    a_i[i] = 1d;
                }
                else
                {
                    a_i[i] = 0;
                }
            }

            // biweight midvariance calculation

            double denomSum = 0;
            double numeratorSum = 0;
            for (int i = 0; i < y_i.Length; i++)
            {
                numeratorSum += a_i[i] * Math.Pow(array[i] - median, 2) * Math.Pow(1 - y_i[i] * y_i[i], 4);
                denomSum += a_i[i] * (1 - 5 * y_i[i] * y_i[i]) * (1 - y_i[i] * y_i[i]);
            }

            return (double)y_i.Length * numeratorSum / Math.Pow(Math.Abs(denomSum), 2);
        }
    }
}
