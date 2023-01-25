using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;
using Chemistry;
using Easy.Common.Extensions;
using FlashLFQ;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Providers.LinearAlgebra;
using NUnit;
using NUnit.Framework;
using OxyPlot.Wpf;
using Plotly.NET.CSharp;
using Proteomics.AminoAcidPolymer;
using SimulatedData;
using SpectralAveraging;
using GenericChartExtensions = Plotly.NET.CSharp.GenericChartExtensions;
using Peptide = Proteomics.AminoAcidPolymer.Peptide;

namespace Test.SimulatedData
{
    public class SimualtedDataTests
    {
        // demonstrates simple Gaussian peak and ensures peak is at the correct index (500). 
        [Test]
        [TestCase(500, 20)]
        public void TestGaussianPeak(double mean, double stddev)
        {
            var gaussianPeak = new GaussianPeakSpectra(mean, stddev, 10000,
                1000, 0, 1.0);
            GenericChartExtensions.Show(gaussianPeak.Plot());
            Assert.AreEqual(mean, gaussianPeak.Yarray.IndexOf(gaussianPeak.Yarray.Max()));
        }

        [Test]
        [TestCase(500, 20)]
        public void TestAddHighFrequencyNoiseToGaussian(double mean, double stddev)
        {
            var gaussianPeak = new GaussianPeakSpectra(mean, stddev, 10000,
                1000, 0, 1.0);
            // max value is ~200 for the intensity multiple of 10000. 
            // therefore, add a mean of 75, stddev of 2.5. 
            Normal hfNoiseDistribution = new(75, 2.5);
            gaussianPeak.AddHighFrequencyNoise(hfNoiseDistribution);
            GenericChartExtensions.Show(gaussianPeak.Plot());

            // because we are adding a random noise distribution, we need to check that the index 
            // of the max value is within a range of possible max values. 
            Assert.That(mean, Is.EqualTo(gaussianPeak.Yarray.IndexOf(gaussianPeak.Yarray.Max())).Within(10));

        }
        // demonstrate adding low frequency noise
        [Test]
        [TestCase(500, 20)]
        public void TestAddLowFrequencyNoiseToGaussian(double mean, double stddev)
        {
            var gaussianPeak = new GaussianPeakSpectra(mean, stddev, 10000,
                1000, 0, 1.0);
            LowFrequencyNoiseParameters lfNoiseParam = new(10, 100, 
                200, 800, 1, 2, 
                5, 100, (450, 550));

            gaussianPeak.AddLowFrequencyNoise(lfNoiseParam);
            GenericChartExtensions.Show(gaussianPeak.Plot());
        }
        // demonstrate combining low frequency noise and high frequency noise 
        [Test]
        [TestCase(500, 20)]
        public void TestAddLowAndHighFreqNoiseToGauss(double mean, double stddev)
        {
            var gaussianPeak = new GaussianPeakSpectra(mean, stddev, 10000,
                1000, 0, 1.0);
            
            Normal hfNoiseDistribution = new(75, 2.5);
            gaussianPeak.AddHighFrequencyNoise(hfNoiseDistribution);

            LowFrequencyNoiseParameters lfNoiseParam = new(10, 100,
                200, 800, 1, 2,
                5, 100, (450, 550));
            gaussianPeak.AddLowFrequencyNoise(lfNoiseParam);
            GenericChartExtensions.Show(gaussianPeak.Plot());
        }
        // demonstrate simulated isotopic envelope
        [Test]
        [TestCase()]
        public void TestSimulatedProtein()
        {
            Peptide peptide = new Peptide("AAAAAAAAA");
            var distribution = IsotopicDistribution.GetDistribution(peptide.GetChemicalFormula());
            double mzLow = 500;
            double mzHigh = 2000;
            double stepSize = 0.1;
            
            int steps = (int)((mzHigh - mzLow) / stepSize); 

            var simProtein = new SimulatedProtein(distribution, 500, 2000,
                steps, stepSize); 
            Console.WriteLine(distribution.Masses.First());
            GenericChartExtensions.Show(simProtein.Plot());
            double yArrayMax = simProtein.Yarray.Max();
            int indexOfMaxYarray = simProtein.Yarray.IndexOf(yArrayMax);

            int expectedIndex = (int)((657.3 - mzLow) / stepSize); // cause zero-based indexing 

            Assert.That(expectedIndex, Is.EqualTo(indexOfMaxYarray - 1));
        }

        [Test]
        [TestCase()]
        public void TestNormalizeByTic()
        {
            Peptide peptide = new Peptide("AAAAAAAAA");
            var distribution = IsotopicDistribution.GetDistribution(peptide.GetChemicalFormula());
            double mzLow = 500;
            double mzHigh = 2000;
            double stepSize = 0.1;

            int steps = (int)((mzHigh - mzLow) / stepSize);

            var simProtein = new SimulatedProtein(distribution, mzLow, mzHigh,
                steps, stepSize);

            simProtein.NormalizeByTic();
            // asserts that the sum of the y array is now equal to 1. 
            Assert.That(1.00, Is.EqualTo(simProtein.Yarray.Sum()).Within(0.01));
        }
        // demonstrates making a charge state envelope
        [Test]
        [TestCase(-0.75, 0.50)]
        [TestCase(-0.75, 0.18)]
        [TestCase(-0.70, 0.18)]
        [TestCase(-0.65, 0.18)]
        public void TestSimulatedChargeStateEnvelope(double mu, double sigma)
        {
            Peptide peptide = new Peptide("MASPDWGYDDKNGPEQWSKLYPIANGNNQSPVD" +
                                          "IKTSETKHDTSLKPISVSYNPATAKEIINVGHSFHVNFEDNDNR" +
                                          "SVLKGGPFSDSYRLFQFHFHWGSTNEHGSEHTVDGVKYSAELHVAHWN" +
                                          "SAKYSSLAEAASKADGLAVIGVLMKVGEANPKLQKVLDALQAIKTKGKRAP" +
                                          "FTNFDPSTLLPSSLDFWTYPGSLTHPPLYESVTWIICKESISVSSEQLAQF" +
                                          "RSLLSNVEGDNAVPMQHNNRPTQPLKGRTVRASF");
            var distribution = IsotopicDistribution.GetDistribution(peptide.GetChemicalFormula());
            double mzLow = 500;
            double mzHigh = 2000;
            double stepSize = 0.01;

            int steps = (int)((mzHigh - mzLow) / stepSize);
            //SimulatedChargeStateEnvelope envelope = new(mzLow, mzHigh, stepSize, 
            //    steps, 15, 60, distribution);
            SimulatedChargeStateEnvelope envelope = new(mzLow, mzHigh, stepSize, 
                steps, 15, 60, distribution, (mu, sigma));
            envelope.NormalizeByTic();
            GenericChartExtensions.Show(envelope.Plot());
        }

        [Test]
        [TestCase(25)]
        public void GaussianDistributionHfNoiseOnly(int numberPeaks)
        {
            List<GaussianPeakSpectra> peaksList = new List<GaussianPeakSpectra>();

            while (numberPeaks > 0)
            {
                GaussianPeakSpectra gaussian = new GaussianPeakSpectra(500, 25, 10000, 1000, 0, 1.0);
                Normal hfNoise = new Normal(100, 50);
                gaussian.AddHighFrequencyNoise(hfNoise);
                gaussian.NormalizeByTic();
                peaksList.Add(gaussian);
                
                numberPeaks--; 
            }

            // perform averaging without rejection
            double[] naiveAverage = peaksList.NaiveAverage();
            double[] xArray = peaksList[0].Xarray;
            
            var averagedPlot = naiveAverage.Plot(xArray);
            averagedPlot.WithTraceInfo("25 spectra, naive average"); 
            
            // get a representative plot
            var originalPlot = peaksList[0].Yarray.Plot(peaksList[0].Xarray);
            originalPlot.WithTraceInfo("Representative plot, unaveraged");

            // perform averaging with rejection
            SpectralAveragingParameters averagingParms = new SpectralAveragingParameters()
            {
                MinSigmaValue = 10.0,
                MaxSigmaValue = 10.0,
                BinSize = 1.0,
                NormalizationType = NormalizationType.AbsoluteToTic,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping, 
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };
            double[][] averagingWithRejectionOutput = peaksList.AverageWithRejection(averagingParms);
            var rejectionPlot = averagingWithRejectionOutput[1].Plot(averagingWithRejectionOutput[0]);

            List<Plotly.NET.GenericChart.GenericChart> plotList = new()
            {
                //originalPlot,
                averagedPlot,
                rejectionPlot
            };
            Plotly.NET.GenericChart.GenericChart charge = Chart.Combine(plotList);
            // create the combined plot: 
            charge.Show();
        }

        [Test]
        [TestCase(25)]
        public void TestEffectOfMinMaxOnSigmaClipping(int numberPeaks = 25)
        {
            SpectralAveragingParameters averagingParmsRelative = new SpectralAveragingParameters()
            {
                MinSigmaValue = 3.0,
                MaxSigmaValue = 3.0,
                BinSize = 1.0,
                NormalizationType = NormalizationType.RelativeToTics,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.SigmaClipping,
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };
            SpectralAveragingParameters averagingParmsAbsolute = new SpectralAveragingParameters()
            {
                MinSigmaValue = 3.0,
                MaxSigmaValue = 3.0,
                BinSize = 1.0,
                NormalizationType = NormalizationType.AbsoluteToTic,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.SigmaClipping,
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };
            SpectralAveragingParameters averagingParmsNone = new SpectralAveragingParameters()
            {
                MinSigmaValue = 3.0,
                MaxSigmaValue = 3.0,
                BinSize = 1.0,
                NormalizationType = NormalizationType.NoNormalization,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.SigmaClipping,
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };
            List<GaussianPeakSpectra> peaksList = new List<GaussianPeakSpectra>();

            while (numberPeaks > 0)
            {
                GaussianPeakSpectra gaussian = new GaussianPeakSpectra(500, 25, 10000, 1000, 0, 1.0);
                Normal hfNoise = new Normal(100, 50);
                gaussian.AddHighFrequencyNoise(hfNoise);
                gaussian.NormalizeByTic();
                peaksList.Add(gaussian);

                numberPeaks--;
            }

            double[][] relative = peaksList.AverageWithRejection(averagingParmsRelative);
            double[][] absolute = peaksList.AverageWithRejection(averagingParmsAbsolute);
            double[][] none = peaksList.AverageWithRejection(averagingParmsNone);
            double[] naive = peaksList.NaiveAverage(); 
            List<Plotly.NET.GenericChart.GenericChart> plotList = new()
            {
                { relative[1].Plot(relative[0]).WithTraceInfo("Relative") },
                { absolute[1].Plot(absolute[0]).WithTraceInfo("Absolute") },
                { none[1].Plot(none[0]).WithTraceInfo("None") },
                naive.Plot(peaksList[0].Xarray).WithTraceInfo("Naive")

            };

            Plotly.NET.GenericChart.GenericChart stackedPlot = Chart.Combine(plotList);
            stackedPlot.Show();

        }
        [Test]
        [Repeat(10)]
        [TestCase()]
        public void TestSigmaClipplingVsNaive(int numberPeaks = 25)
        {
            SpectralAveragingParameters sigma = new SpectralAveragingParameters()
            {
                MinSigmaValue = 3.0,
                MaxSigmaValue = 3.0,
                BinSize = 1.0,
                NormalizationType = NormalizationType.RelativeToTics,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.SigmaClipping,
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };
            SpectralAveragingParameters winsorized = new SpectralAveragingParameters()
            {
                MinSigmaValue = 3.0,
                MaxSigmaValue = 3.0,
                BinSize = 1.0,
                NormalizationType = NormalizationType.RelativeToTics,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.WinsorizedSigmaClipping,
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };
            SpectralAveragingParameters averageSigma = new SpectralAveragingParameters()
            {
                MinSigmaValue = 3.0,
                MaxSigmaValue = 3.0,
                BinSize = 1.0,
                NormalizationType = NormalizationType.RelativeToTics,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping,
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };
            SpectralAveragingParameters noRejection = new SpectralAveragingParameters()
            {
                MinSigmaValue = 3.0,
                MaxSigmaValue = 3.0,
                BinSize = 1.0,
                NormalizationType = NormalizationType.RelativeToTics,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.NoRejection,
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };
            List<GaussianPeakSpectra> peaksList = new List<GaussianPeakSpectra>();

            while (numberPeaks > 0)
            {
                GaussianPeakSpectra gaussian = new GaussianPeakSpectra(500, 10, 10000, 1000, 0, 1.0);
                Normal hfNoise = new Normal(100, 10);
                gaussian.AddHighFrequencyNoise(hfNoise);
                gaussian.NormalizeByTic();
                peaksList.Add(gaussian);

                numberPeaks--;
            }

            double[][] sigmaClipping = peaksList.AverageWithRejection(sigma);
            double[][] winsorizedClipping = peaksList.AverageWithRejection(winsorized);
            double[][] averageSigmaClipping = peaksList.AverageWithRejection(averageSigma);
            double[][] noRejectionClipping = peaksList.AverageWithRejection(noRejection);
            double[] naive = peaksList.NaiveAverage();
            List<Plotly.NET.GenericChart.GenericChart> plotList = new()
            {
                { sigmaClipping[1].Plot(sigmaClipping[0]).WithTraceInfo("Sigma") },
                { winsorizedClipping[1].Plot(winsorizedClipping[0]).WithTraceInfo("Winsorized") },
                { averageSigmaClipping[1].Plot(averageSigmaClipping[0]).WithTraceInfo("Average Sigma") },
                noRejectionClipping[1].Plot(noRejectionClipping[0]).WithTraceInfo("No rejection"),
                naive.Plot(peaksList[0].Xarray).WithTraceInfo("Naive"),
                peaksList[0].Plot().WithTraceInfo("Original")
            };

            Plotly.NET.GenericChart.GenericChart stackedPlot = Chart.Combine(plotList);
            var layoutPlotly = Plotly.NET.Layout.init<double>(Width: 2000);
            Plotly.NET.GenericChartExtensions.WithLayout(stackedPlot, layoutPlotly).Show();

            List<double[]> yarrays = new()
            {
                sigmaClipping[1], 
                winsorizedClipping[1], 
                averageSigmaClipping[1], 
                noRejectionClipping[1], 
                naive,
                peaksList[0].Yarray
            };
            List<double> noiseEstimates = new();
            List<double> stddevNoiseEstimate = new();
            
            foreach (double[] yarray in yarrays)
            {
                MRSNoiseEstimator.MRSNoiseEstimation(yarray[..400], 0.1, out double noiseEstimate); 
                noiseEstimates.Add(noiseEstimate);
                stddevNoiseEstimate.Add(BasicStatistics.CalculateStandardDeviation(yarray[..400]));
            }
            // TODO: Calculate the ENR instead of improvement ratio. Need to account for scale. 
            List<double> improvementRatio = noiseEstimates.Select(i => noiseEstimates.Last() / i).ToList();
            List<string> method = new()
            {
                "Sigma Clipping", 
                "Winsorized Sigma Clipping", 
                "AveragedSigma Clipping", 
                "No Rejection", 
                "Averaging Without Rejection", 
                "Raw Reference Spectrum"
            }; 
            Console.WriteLine("Method\tMRS Noise Estimate\tStddev Noise Estimate\tImprovement Ratio (Raw Noise Estimate / Averaging Method Noise Estimate)");
            for (int i = 0; i < noiseEstimates.Count; i++)
            {
                Console.WriteLine("{0}\t{1}\t{2}\t{3}",
                    method[i], noiseEstimates[i], stddevNoiseEstimate[i], improvementRatio[i]);
            }
        }
    }

    public static class ConvenienceExtensions
    {
        public static Plotly.NET.GenericChart.GenericChart Plot(this double[] yArray, double[] xArray)
        {
            return Chart.Line<double, double, string>(x: xArray, y: yArray)
                .WithXAxisStyle<double, double, string>("m/z")
                .WithYAxisStyle<double, double, string>("Intensity"); 
        }
    }
}
