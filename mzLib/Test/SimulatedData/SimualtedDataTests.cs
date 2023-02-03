using System;
using System.Collections.Generic;
using System.DirectoryServices.ActiveDirectory;
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
using MathNet.Numerics.Statistics;
using Newtonsoft.Json.Linq;
using NUnit;
using NUnit.Framework;
using OxyPlot.Wpf;
using Plotly.NET;
using Plotly.NET.CSharp;
using Plotly.NET.TraceObjects;
using Proteomics.AminoAcidPolymer;
using SimulatedData;
using SpectralAveraging;
using Test.SimulatedData;
using Chart = Plotly.NET.CSharp.Chart;
using GenericChart = Plotly.NET.GenericChart;
using GenericChartExtensions = Plotly.NET.CSharp.GenericChartExtensions;
using Peptide = Proteomics.AminoAcidPolymer.Peptide;

namespace Test
{
    public class SimualtedDataTests
    {
        const string carbonicAnhydrase = "MASPDWGYDDKNGPEQWSKLYPIANGNNQS" +
                                         "PVDIKTSETKHDTSLKPISVSYNPATAKEIINVGHSFH" +
                                         "VNFEDNDNRSVLKGGPFSDSYRLFQFHFHWGSTNEHGSEHTV" +
                                         "DGVKYSAELHVAHWNSAKYSSLAEAASKADGLAVIGVLMKVG" +
                                         "EANPKLQKVLDALQAIKTKGKRAPFTNFDPSTLLPSSLDFWTY" +
                                         "PGSLTHPPLYESVTWIICKESISVSSEQLAQFRSLLSNVEG" +
                                         "DNAVPMQHNNRPTQPLKGRTVRASF";
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
            GenericChartExtensions.WithTraceInfo(averagedPlot, "25 spectra, naive average"); 
            
            // get a representative plot
            var originalPlot = peaksList[0].Yarray.Plot(peaksList[0].Xarray);
            GenericChartExtensions.WithTraceInfo(originalPlot, "Representative plot, unaveraged");

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
            GenericChartExtensions.Show(charge);
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
                { GenericChartExtensions.WithTraceInfo(relative[1].Plot(relative[0]), "Relative") },
                { GenericChartExtensions.WithTraceInfo(absolute[1].Plot(absolute[0]), "Absolute") },
                { GenericChartExtensions.WithTraceInfo(none[1].Plot(none[0]), "None") },
                GenericChartExtensions.WithTraceInfo(naive.Plot(peaksList[0].Xarray), "Naive")

            };

            Plotly.NET.GenericChart.GenericChart stackedPlot = Chart.Combine(plotList);
            GenericChartExtensions.Show(stackedPlot);

        }
        [Test]
        [Repeat(1)]
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
                Normal hfNoise = new Normal(1000, 500);
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
                { GenericChartExtensions.WithTraceInfo(sigmaClipping[1].Plot(sigmaClipping[0]), "Sigma") },
                { GenericChartExtensions.WithTraceInfo(winsorizedClipping[1].Plot(winsorizedClipping[0]), "Winsorized") },
                { GenericChartExtensions.WithTraceInfo(averageSigmaClipping[1].Plot(averageSigmaClipping[0]), "Average Sigma") },
                GenericChartExtensions.WithTraceInfo(noRejectionClipping[1].Plot(noRejectionClipping[0]), "No rejection"),
                GenericChartExtensions.WithTraceInfo(naive.Plot(peaksList[0].Xarray), "Naive"),
                GenericChartExtensions.WithTraceInfo(peaksList[0].Plot(), "Original")
            };

            Plotly.NET.GenericChart.GenericChart stackedPlot = Chart.Combine(plotList);
            var layoutPlotly = Plotly.NET.Layout.init<double>(Width: 1000);
            GenericChartExtensions.Show(Plotly.NET.GenericChartExtensions.WithLayout(stackedPlot, layoutPlotly));

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
                MRSNoiseEstimator.MRSNoiseEstimation(yarray[..400], 0.01, out double noiseEstimate); 
                noiseEstimates.Add(noiseEstimate);
                stddevNoiseEstimate.Add(BasicStatistics.CalculateStandardDeviation(yarray[..400]));
            }

            double refNoiseEstimate = noiseEstimates.Last();
            double refScaleEstimate = Math.Sqrt(BasicStatistics.BiweightMidvariance(yarrays.Last()));
            List<double> enrList = yarrays
                .Select((i, index) =>
                    i.CalculateEnr(refScaleEstimate, noiseEstimates[index], refNoiseEstimate))
                .ToList();
            List<string> method = new()
            {
                "Sigma Clipping", 
                "Winsorized Sigma Clipping", 
                "AveragedSigma Clipping", 
                "No Rejection", 
                "Averaging Without Rejection", 
                "Raw Reference Spectrum"
            }; 
            Console.WriteLine("Method\tMRS Noise Estimate\tStddev Noise Estimate\tENR");
            for (int i = 0; i < noiseEstimates.Count; i++)
            {
                Console.WriteLine("{0}\t{1}\t{2}\t{3}",
                    method[i], noiseEstimates[i], stddevNoiseEstimate[i], enrList[i]);
            }

            // create a histogram of the values
            double[] testHistogramYarray = yarrays[0];
            Plotly.NET.GenericChart.GenericChart histogramSigmaClipping =
                Plotly.NET.CSharp.Chart.Histogram<double, double, string>(testHistogramYarray[..400]);

            var histogramRaw = Chart.Histogram<double, double, string>(yarrays.Last()[..400]);

            var stackedHistograms = Chart.Combine(new Plotly.NET.GenericChart.GenericChart[]
                { histogramSigmaClipping, histogramRaw });
            GenericChartExtensions.Show(stackedHistograms);

            double sigmaClippingStdev = BasicStatistics.CalculateStandardDeviation(yarrays[0][..400]); 
            double rawStdev = BasicStatistics.CalculateStandardDeviation(yarrays.Last()[..400]); 
            double ratioStddev = rawStdev / sigmaClippingStdev;
            Console.WriteLine("Sigma Clipping: {0}\tRaw: {1}\tImprovement: {2}", 
                sigmaClippingStdev, rawStdev, ratioStddev);
            List<double> noiseStddevsAll =
                peaksList.Select(i => i.Yarray[..400].StandardDeviation()).ToList(); 
            Console.WriteLine(noiseStddevsAll.Mean() +"\t" + noiseStddevsAll.StandardDeviation());
            // need to flesh out more that the 
        }

        [Test]
        [TestCase()]
        public void TestLfNoise(int numberPeaks = 25)
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

            LowFrequencyNoiseParameters lfNoiseParams = new LowFrequencyNoiseParameters(1, 15, 200, 800,
                0.25, 0.75, 200, 1000); 
            while (numberPeaks > 0)
            {
                GaussianPeakSpectra gaussian = new GaussianPeakSpectra(500, 10, 10000, 1000, 0, 1.0);
                Normal hfNoise = new Normal(1000, 10);
                gaussian.AddHighFrequencyNoise(hfNoise);
                gaussian.AddLowFrequencyNoise(lfNoiseParams);
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
                { GenericChartExtensions.WithTraceInfo(sigmaClipping[1].Plot(sigmaClipping[0]), "Sigma") },
                { GenericChartExtensions.WithTraceInfo(winsorizedClipping[1].Plot(winsorizedClipping[0]), "Winsorized") },
                { GenericChartExtensions.WithTraceInfo(averageSigmaClipping[1].Plot(averageSigmaClipping[0]), "Average Sigma") },
                GenericChartExtensions.WithTraceInfo(noRejectionClipping[1].Plot(noRejectionClipping[0]), "No rejection"),
                GenericChartExtensions.WithTraceInfo(naive.Plot(peaksList[0].Xarray), "Naive"),
                GenericChartExtensions.WithTraceInfo(peaksList[0].Plot(), "Original")
            };

            Plotly.NET.GenericChart.GenericChart stackedPlot = Chart.Combine(plotList);
            var layoutPlotly = Plotly.NET.Layout.init<double>(Width: 1000);
            GenericChartExtensions.Show(Plotly.NET.GenericChartExtensions.WithLayout(stackedPlot, layoutPlotly));

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
                MRSNoiseEstimator.MRSNoiseEstimation(yarray[200..400], 0.01, out double noiseEstimate);
                noiseEstimates.Add(noiseEstimate);
                stddevNoiseEstimate.Add(BasicStatistics.CalculateStandardDeviation(yarray[200..400]));
            }

            double refNoiseEstimate = noiseEstimates.Last();
            double refScaleEstimate = Math.Sqrt(BasicStatistics.BiweightMidvariance(yarrays.Last()));
            List<double> enrList = yarrays
                .Select((i, index) =>
                    i.CalculateEnr(refScaleEstimate, noiseEstimates[index], refNoiseEstimate))
                .ToList();
            List<string> method = new()
            {
                "Sigma Clipping",
                "Winsorized Sigma Clipping",
                "AveragedSigma Clipping",
                "No Rejection",
                "Averaging Without Rejection",
                "Raw Reference Spectrum"
            };
            Console.WriteLine("Method\tMRS Noise Estimate\tStddev Noise Estimate\tENR");
            for (int i = 0; i < noiseEstimates.Count; i++)
            {
                Console.WriteLine("{0}\t{1}\t{2}\t{3}",
                    method[i], noiseEstimates[i], stddevNoiseEstimate[i], enrList[i]);
            }
            // based on the experimentation in these tests, it looks like the best use case of averaging without rejection is when 
            // there are few low-frequency noise peaks relative to the total data points. 
        }

        [Test]
        [Repeat(5)]
        public void LfNoiseRejection()
        {

            SpectralAveragingParameters sigma = new SpectralAveragingParameters()
            {
                MinSigmaValue = 1.45,
                MaxSigmaValue = 1.45,
                BinSize = 1.0,
                NormalizationType = NormalizationType.RelativeToTics,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.SigmaClipping,
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };

            Dictionary<int, LowFrequencyNoiseParameters> lfParamsDict = new()
            {
                {1, new LowFrequencyNoiseParameters(1,1,10,990,0.25,0.75,50, 200)}, 
                {2, new LowFrequencyNoiseParameters(1,2,10,990,0.25,0.75,50, 200)}, 
                {3, new LowFrequencyNoiseParameters(1,3,10,990,0.25,0.75,50, 200)}, 
                {4, new LowFrequencyNoiseParameters(1,4,10,990,0.25,0.75,50, 200)}, 
                {5, new LowFrequencyNoiseParameters(1,5,10,990,0.25,0.75,50, 200)}, 
                {6, new LowFrequencyNoiseParameters(1,6,10,990,0.25,0.75,50, 200)}, 
                {7, new LowFrequencyNoiseParameters(1,7,10,990,0.25,0.75,50, 200)}, 
                {8, new LowFrequencyNoiseParameters(1,8,10,990,0.25,0.75,50, 200)}, 
                {9, new LowFrequencyNoiseParameters(1,9,10,990,0.25,0.75,50, 200)}, 
                {10, new LowFrequencyNoiseParameters(1,10,10,990,0.25,0.75,50, 200)}, 
                {11, new LowFrequencyNoiseParameters(1,11,10,990,0.25,0.75,50, 200)}, 
                {12, new LowFrequencyNoiseParameters(1,12,10,990,0.25,0.75,50, 200)}, 
                {13, new LowFrequencyNoiseParameters(1,13,10,990,0.25,0.75,50, 200)}, 
                {14, new LowFrequencyNoiseParameters(1,14,10,990,0.25,0.75,50, 200)}, 
                {15, new LowFrequencyNoiseParameters(1,15,10,990,0.25,0.75,50, 200)},
                {20, new LowFrequencyNoiseParameters(1,20,10,990,0.25,0.75,50, 200)},
                {25, new LowFrequencyNoiseParameters(1,25,10,990,0.25,0.75,50, 200)},
                {30, new LowFrequencyNoiseParameters(1,30,10,990,0.25,0.75,50, 200)},
                {35, new LowFrequencyNoiseParameters(1,35,10,990,0.25,0.75,50, 200)},
                {40, new LowFrequencyNoiseParameters(1,40,10,990,0.25,0.75,50, 200)},
                {45, new LowFrequencyNoiseParameters(1,45,10,990,0.25,0.75,50, 200)},
                {50, new LowFrequencyNoiseParameters(1,50,10,990,0.25,0.75,50, 200)},
                {55, new LowFrequencyNoiseParameters(1,55,10,990,0.25,0.75,50, 200)},
                {60, new LowFrequencyNoiseParameters(1,60,10,990,0.25,0.75,50, 200)},
                {65, new LowFrequencyNoiseParameters(1,65,10,990,0.25,0.75,50, 200)},
                {70, new LowFrequencyNoiseParameters(1,70,10,990,0.25,0.75,50, 200)},
                {75, new LowFrequencyNoiseParameters(1,75,10,990,0.25,0.75,50, 200)},
                {80, new LowFrequencyNoiseParameters(1,80,10,990,0.25,0.75,50, 200)},
                {85, new LowFrequencyNoiseParameters(1,85,10,990,0.25,0.75,50, 200)},
                {90, new LowFrequencyNoiseParameters(1,90,10,990,0.25,0.75,50, 200)},
            };

            List<double[][]> rejectionAveragedOutputs = new();
            List<double[]> naiveAveragedOutputs = new();
            List<List<GaussianPeakSpectra>> peaksList = new();
            List<(double refScaleEstimate, double refNoiseEstimate)> scaleNoiseEstimates =
                new List<(double refScaleEstimate, double refNoiseEstimate)>(); 
            foreach (var kvp in lfParamsDict)
            {
                int numberPeaks = 25;
                List<GaussianPeakSpectra> innerTempList = new();
                while (numberPeaks > 0)
                {
                    GaussianPeakSpectra gaussian = new GaussianPeakSpectra(500, 10, 10000, 1000, 0, 1.0);
                    Normal hfNoise = new Normal(1000, 10);
                    gaussian.AddHighFrequencyNoise(hfNoise);
                    gaussian.AddLowFrequencyNoise(kvp.Value);
                    gaussian.NormalizeByTic();
                    innerTempList.Add(gaussian);

                    numberPeaks--;
                }
                peaksList.Add(innerTempList);
                rejectionAveragedOutputs.Add(innerTempList.AverageWithRejection(sigma));
                naiveAveragedOutputs.Add(innerTempList.NaiveAverage());
                MRSNoiseEstimator.MRSNoiseEstimation(innerTempList[0].Yarray, 0.01, out double noise); 
                scaleNoiseEstimates.Add((Math.Sqrt(BasicStatistics.BiweightMidvariance(innerTempList[0].Yarray)), noise)); 
            }
            
            // calculate the ENR for each set of lf peaks added
            List< (double rejection, double noRejection)> enrs = new();
            int j = 0;
            while (j < rejectionAveragedOutputs.Count)
            {
                //var rejectionPlot = rejectionAveragedOutputs[j][1].Plot(rejectionAveragedOutputs[j][0]);
                //var naivePlot = naiveAveragedOutputs[j].Plot(rejectionAveragedOutputs[j][0]);
                //var rawPlot = peaksList[j][0].Yarray.Plot(peaksList[j][0].Xarray);
                //List<Plotly.NET.GenericChart.GenericChart> plotList = new()
                //{
                //    { GenericChartExtensions.WithTraceInfo(rejectionPlot, "with Rejection") },
                //    { GenericChartExtensions.WithTraceInfo(naivePlot, "without Rejection") },
                //    { GenericChartExtensions.WithTraceInfo(rawPlot, "raw Spectra") }
                //};
                //var stackedPlot = Chart.Combine(plotList);
                //stackedPlot.Show();

                double enr = rejectionAveragedOutputs[j][1].CalculateEnr(scaleNoiseEstimates[j].refScaleEstimate, null,
                        scaleNoiseEstimates[j].refNoiseEstimate);
                double naiveEnr = naiveAveragedOutputs[j].CalculateEnr(scaleNoiseEstimates[j].refScaleEstimate,
                    null, scaleNoiseEstimates[j].refNoiseEstimate);
                enrs.Add((enr, naiveEnr));
                j++;
            }

            int maxPeaksIndexer = 0;
            int[] maxPeaks = Enumerable.Range(1,30).ToArray(); 
            foreach (var valueTuple in enrs)
            {
                Console.WriteLine("{0}\t{1}\t{2}", maxPeaks[maxPeaksIndexer], valueTuple.rejection, valueTuple.noRejection);
                if (++maxPeaksIndexer == maxPeaks.Length)
                {
                    maxPeaksIndexer = 0; 
                }
            }
        }

        [Test]
        [TestCase(-0.65, 0.18, true)]
        [Repeat(5)]
        public void ProteinSpectraNoiseTest(double mu, double sigma, bool plotting)
        {
            const string carbonicAnhydrase = "MASPDWGYDDKNGPEQWSKLYPIANGNNQS" +
                                             "PVDIKTSETKHDTSLKPISVSYNPATAKEIINVGHSFH" +
                                             "VNFEDNDNRSVLKGGPFSDSYRLFQFHFHWGSTNEHGSEHTV" +
                                             "DGVKYSAELHVAHWNSAKYSSLAEAASKADGLAVIGVLMKVG" +
                                             "EANPKLQKVLDALQAIKTKGKRAPFTNFDPSTLLPSSLDFWTY" +
                                             "PGSLTHPPLYESVTWIICKESISVSSEQLAQFRSLLSNVEG" +
                                             "DNAVPMQHNNRPTQPLKGRTVRASF"; 
            Peptide peptide = new Peptide(carbonicAnhydrase); // from uniprot 
            var distribution = IsotopicDistribution.GetDistribution(peptide.GetChemicalFormula(), 0.005, 0.001);
            double mzLow = distribution.Masses.Min() - 5d;
            double mzHigh = distribution.Masses.Max() + 5d; 
            double stepSize = 0.005;
            SpectralAveragingParameters sigmaParams = new SpectralAveragingParameters()
            {
                MinSigmaValue = 1.0,
                MaxSigmaValue = 3.0,
                BinSize = 0.0049,
                NormalizationType = NormalizationType.RelativeToTics,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.SigmaClipping,
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };

            int steps = (int)((mzHigh - mzLow) / stepSize);

            Dictionary<int, LowFrequencyNoiseParameters> lfParamsDict = new()
            {
                {1, new LowFrequencyNoiseParameters(1, (int)(steps * 0.001),mzLow,mzHigh,
                    0.005,0.001,0.000001, 0.0005)},
                {2, new LowFrequencyNoiseParameters(1,(int)(steps * 0.01),mzLow,mzHigh,
                    0.005,0.001,0.000001, 0.0005)},
                {3, new LowFrequencyNoiseParameters(1,(int)(steps * 0.02),mzLow,mzHigh,
                    0.005,0.001,0.000001, 0.0005)},
                {4, new LowFrequencyNoiseParameters(1,(int)(steps * 0.05),mzLow,mzHigh,
                    0.005,0.001,0.000001, 0.0005)},
                {5, new LowFrequencyNoiseParameters(1,(int)(steps * 0.1),mzLow,mzHigh,
                    0.005,0.001,0.000001, 0.0005)},
                {6, new LowFrequencyNoiseParameters(1,(int)(steps * 0.25),mzLow,mzHigh,
                    0.005,0.001,0.000001, 0.0005)},
                {7, new LowFrequencyNoiseParameters(1,(int)(steps * 0.50),mzLow,mzHigh,
                    0.005,0.001,0.000001, 0.0005)},
                {8, new LowFrequencyNoiseParameters(1,(int)(steps * 0.75),mzLow,mzHigh,
                    0.005,0.001,0.000001, 0.0005)},
                {9, new LowFrequencyNoiseParameters(1,(int)(steps * 0.9),mzLow,mzHigh,
                    0.005,0.001,0.000001, 0.0005)}
            };
            List<double[][]> rejectionAveragedOutputs = new();
            List<double[]> naiveAveragedOutputs = new();
            List<List<SimulatedChargeStateEnvelope>> peaksList = new();
            List<(double refScaleEstimate, double refNoiseEstimate)> scaleNoiseEstimates =
                new List<(double refScaleEstimate, double refNoiseEstimate)>();
            foreach (var kvp in lfParamsDict)
            {
                int numberPeaks = 10;
                List<SimulatedChargeStateEnvelope> innerTempList = new();
                while (numberPeaks > 0)
                {
                    SimulatedChargeStateEnvelope cse = new(mzLow, mzHigh, stepSize, steps, 
                        1, 2, distribution);
                    Normal hfNoise = new Normal(0.001, 0.001);
                    cse.AddHighFrequencyNoise(hfNoise);
                    cse.AddLowFrequencyNoise(kvp.Value);
                    cse.NormalizeByTic();
                    innerTempList.Add(cse);

                    numberPeaks--;
                }
                peaksList.Add(innerTempList);
                rejectionAveragedOutputs.Add(innerTempList.AverageWithRejection(sigmaParams));
                naiveAveragedOutputs.Add(innerTempList.NaiveAverage());
                // normalize y array 
                double[] normalizedYarray = NormalizeArray(innerTempList[0].Yarray, 0, 1);

                MRSNoiseEstimator.MRSNoiseEstimation(normalizedYarray, 0.01, out double noise);
                scaleNoiseEstimates.Add((Math.Sqrt(BasicStatistics.BiweightMidvariance(normalizedYarray)), noise));
            }

            // calculate the ENR for each set of lf peaks added
            List<(double rejection, double noRejection)> enrs = new();
            int j = 0;
            while (j < rejectionAveragedOutputs.Count)
            {
                #region Plotting

                if (plotting)
                {
                    var rejectionPlot = NormalizeArray(rejectionAveragedOutputs[j][1], 0, 1)
                        .Plot(rejectionAveragedOutputs[j][0]);
                    var naivePlot = NormalizeArray(naiveAveragedOutputs[j], 0,1)
                        .Plot(peaksList[j][0].Xarray);
                    var rawPlot = NormalizeArray(peaksList[j][0].Yarray, 0,1)
                        .Plot(peaksList[j][0].Xarray);

                    List<Plotly.NET.GenericChart.GenericChart> plotList = new()
                    {
                        { GenericChartExtensions.WithTraceInfo(rejectionPlot, "with Rejection") },
                        { GenericChartExtensions.WithTraceInfo(naivePlot, "without Rejection") },
                        { GenericChartExtensions.WithTraceInfo(rawPlot, "raw Spectra") },
                    };

                    // combine theoretical with each plot
                    List<GenericChart.GenericChart> plotListCombined = new();
                    int maximumListIndexer = 0; 
                    foreach (var plot in plotList)
                    {
                        Peptide ca = new Peptide(carbonicAnhydrase);
                        var distributionInner = IsotopicDistribution.GetDistribution(ca.GetChemicalFormula(), 0.05, 0.001);
                        var distributionMax = distributionInner.Intensities.Max(); 

                        Plotly.NET.Color markerColor = Color.fromKeyword(ColorKeyword.Red);
                        var markerSymbol = StyleParam.MarkerSymbol.NewModified(StyleParam.MarkerSymbol.Circle, StyleParam.SymbolStyle.Open);
                        Marker marker = Marker.init(Symbol: markerSymbol, Color: markerColor);

                        var theoreticalPlot = Chart.Point<double, double, string>(distributionInner.Masses,
                            distributionInner.Intensities.Select(z => z / distributionMax), Marker: marker);
                        GenericChartExtensions.WithTraceInfo(theoreticalPlot, "Theoretical Isotopic Distribution"); 
                        plotListCombined.Add(Chart.Combine(new []{ theoreticalPlot, plot })); 
                    }

                    GenericChartExtensions.Show(Plotly.NET.Chart.Grid<IEnumerable<Plotly.NET.GenericChart.GenericChart>>(3, 1).Invoke(plotListCombined)
                            .WithSize(1000));

                    //var stackedPlot = Chart.Combine(plotList);
                    //stackedPlot.Show();
                }

                #endregion

                double enr = NormalizeArray(rejectionAveragedOutputs[j][1], 0, 1)
                    .CalculateEnr(scaleNoiseEstimates[j].refScaleEstimate, null,
                        scaleNoiseEstimates[j].refNoiseEstimate);
                double naiveEnr = NormalizeArray(naiveAveragedOutputs[j],0,1)
                    .CalculateEnr(scaleNoiseEstimates[j].refScaleEstimate,
                    null, scaleNoiseEstimates[j].refNoiseEstimate);
                enrs.Add((enr, naiveEnr));
                j++;
            }

            int maxPeaksIndexer = 0;
            int[] maxPeaks = Enumerable.Range(1, lfParamsDict.Count).ToArray();
            foreach (var valueTuple in enrs)
            {
                Console.WriteLine("{0}\t{1}\t{2}", maxPeaks[maxPeaksIndexer], valueTuple.rejection, valueTuple.noRejection);
                if (++maxPeaksIndexer == maxPeaks.Length)
                {
                    maxPeaksIndexer = 0;
                }
            }
        }

        [Test]
        [TestCase(true, 1)]
        [TestCase(true, 5)]
        [TestCase(true, 7)]
        public void LowResolutionProtein(bool plotting, int iterations)
        {
            // from uniprot 

            Peptide peptide = new Peptide(carbonicAnhydrase); 
            var distribution = IsotopicDistribution.GetDistribution(peptide.GetChemicalFormula(), 0.001, 0.001);
            double mzLow = distribution.Masses.Min() - 5d;
            double mzHigh = distribution.Masses.Max() + 5d;
            double stepSize = 0.01;
            SpectralAveragingParameters sigmaParams = new SpectralAveragingParameters()
            {
                MinSigmaValue = 1.0,
                MaxSigmaValue = 3.0,
                BinSize = 0.0049,
                NormalizationType = NormalizationType.RelativeToTics,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.SigmaClipping,
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };

            int steps = (int)((mzHigh - mzLow) / stepSize);
            SimulatedChargeStateEnvelope cse = new(mzLow, mzHigh, stepSize, steps,
                1, 2, distribution);

            cse.Blur(75, 0.75, iterations);
            cse.NormalizeToMaxIntensity();
            Normal hfNoise = new Normal(0.01, 0.005);
            cse.AddHighFrequencyNoise(hfNoise);
            var lfNoiseParams = new LowFrequencyNoiseParameters(
                1, (int)(steps * 0.5),
                mzLow, mzHigh,
                0.01, 0.02,
                0.00001, 0.01);
            cse.AddLowFrequencyNoise(lfNoiseParams); 
            if (plotting)
                GenericChartExtensions.Show(cse.Plot().WithTitle(string.Join(" = ", "Blur iterations", iterations)));

        }
        [Test]
        [TestCase(true, 5)]
        [Repeat(5)]
        public void LowResolutionProteinRejection(bool plotting, int iterations)
        {
            // from uniprot 

            Peptide peptide = new Peptide(carbonicAnhydrase);
            var distribution = IsotopicDistribution.GetDistribution(peptide.GetChemicalFormula(), 0.001, 0.001);
            double mzLow = distribution.Masses.Min() - 5d;
            double mzHigh = distribution.Masses.Max() + 5d;
            double stepSize = 0.001;
            SpectralAveragingParameters sigmaParams = new SpectralAveragingParameters()
            {
                MinSigmaValue = 1.5,
                MaxSigmaValue = 1.5,
                BinSize = 0.0049,
                NormalizationType = NormalizationType.RelativeToTics,
                NumberOfScansToAverage = 25,
                OutlierRejectionType = OutlierRejectionType.SigmaClipping,
                SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation
            };

            int steps = (int)((mzHigh - mzLow) / stepSize);
            SimulatedChargeStateEnvelope cse = new(mzLow, mzHigh, stepSize, steps,
                1, 2, distribution);
           
            Dictionary<int, LowFrequencyNoiseParameters> lfNoiseDict = new(); 
            var lfNoiseParams = new LowFrequencyNoiseParameters(
                1, (int)(steps * 0.01),
                mzLow, mzHigh,
                0.01, 0.02,
                0.00001, 0.01);
            lfNoiseDict.Add(1, lfNoiseParams);
            var cseParams = new SimulatedChargeStateEnvelopeParams(SimulatedDataParamsType.SimulatedChargeStateEnvelope,
                (mzLow, mzHigh), stepSize, steps, (1, 2), distribution, null); 

            SimulatePlotWrite(cseParams, lfNoiseDict, sigmaParams, plotting, 
                new Normal(0.1, 0.001), lfNoiseParams, (503,0.75,iterations), true, 
                CreateSimulatedChargeState);
        }
        public global::SimulatedData.SimulatedData CreateSimulatedChargeState(SimulatedDataParams cseParams, 
            Normal? hfNoise, LowFrequencyNoiseParameters? lfNoiseParams, (int width, double sigma, int iterations)? blurParams, bool normalize)
        {
            var cse = new SimulatedChargeStateEnvelope(cseParams.MzRange.Item1,
                cseParams.MzRange.Item2, cseParams.StepSize, cseParams.Length, cseParams.ChargeStateRange.Item1,
                cseParams.ChargeStateRange.Item2, cseParams.ParentDistribution);
            if (blurParams != null)
            {
                cse.Blur(blurParams.Value.width, blurParams.Value.sigma, blurParams.Value.iterations);
            }

            if (normalize)
            {
                cse.NormalizeToMaxIntensity();
            }
            
            if (hfNoise != null)
            {
                cse.AddHighFrequencyNoise(hfNoise);
            }

            if (lfNoiseParams != null)
            {
                cse.AddLowFrequencyNoise(lfNoiseParams);
            }

            return cse; 
        }
        
        public void SimulatePlotWrite(SimulatedDataParams cseParams, 
            Dictionary<int,LowFrequencyNoiseParameters> lfParamsDict, 
            SpectralAveragingParameters averagingParams, 
            bool plotting, 
            Normal? hfNoise, LowFrequencyNoiseParameters? lfNoise, (int, double, int)? blurParams,
            bool normalize, 
            Func<SimulatedDataParams, Normal?, LowFrequencyNoiseParameters?, 
                (int, double, int)?, bool, global::SimulatedData.SimulatedData> dataCreator)
        {
            List<double[][]> rejectionAveragedOutputs = new();
            List<double[]> naiveAveragedOutputs = new();
            List<List<global::SimulatedData.SimulatedData>> peaksList = new();
            List<(double refScaleEstimate, double refNoiseEstimate)> scaleNoiseEstimates =
                new List<(double refScaleEstimate, double refNoiseEstimate)>();
            foreach (var kvp in lfParamsDict)
            {
                int numberPeaks = 10;
                List<global::SimulatedData.SimulatedData> innerTempList = new();
                while (numberPeaks > 0)
                {
                    var output = dataCreator.Invoke(cseParams, new Normal(hfNoise.Mean, hfNoise.StdDev), 
                        lfNoise, blurParams, normalize);
                    innerTempList.Add(output);
                    numberPeaks--;
                }
                peaksList.Add(innerTempList);
                rejectionAveragedOutputs.Add(innerTempList.AverageWithRejection(averagingParams));
                naiveAveragedOutputs.Add(innerTempList.NaiveAverage());
                
                double[] normalizedYarray = NormalizeArray(innerTempList[0].Yarray, 0.1, 1);
                MRSNoiseEstimator.MRSNoiseEstimation(normalizedYarray, 0.01, out double noise);
                scaleNoiseEstimates.Add((Math.Sqrt(BasicStatistics.BiweightMidvariance(normalizedYarray)), noise));
            }

            // calculate the ENR for each set of lf peaks added
            List<(double rejection, double noRejection)> enrs = new();
            int j = 0;
            while (j < rejectionAveragedOutputs.Count)
            {
                #region Plotting

                if (plotting)
                {
                    var rejectionPlot = NormalizeArray(rejectionAveragedOutputs[j][1], 0, 1)
                        .Plot(rejectionAveragedOutputs[j][0]);
                    var naivePlot = NormalizeArray(naiveAveragedOutputs[j], 0, 1)
                        .Plot(peaksList[j][0].Xarray);
                    var rawPlot = NormalizeArray(peaksList[j][0].Yarray, 0, 1)
                        .Plot(peaksList[j][0].Xarray);

                    List<Plotly.NET.GenericChart.GenericChart> plotList = new()
                    {
                        { GenericChartExtensions.WithTraceInfo(rejectionPlot, "with Rejection") },
                        { GenericChartExtensions.WithTraceInfo(naivePlot, "without Rejection") },
                        { GenericChartExtensions.WithTraceInfo(rawPlot, "raw Spectra") },
                    };

                    // combine theoretical with each plot
                    List<GenericChart.GenericChart> plotListCombined = new();
                    int maximumListIndexer = 0;
                    foreach (var plot in plotList)
                    {
                        Peptide ca = new Peptide(carbonicAnhydrase);
                        var distributionInner = IsotopicDistribution.GetDistribution(ca.GetChemicalFormula(), 0.05, 0.001);
                        var distributionMax = distributionInner.Intensities.Max();

                        Plotly.NET.Color markerColor = Color.fromKeyword(ColorKeyword.Red);
                        var markerSymbol = StyleParam.MarkerSymbol.NewModified(StyleParam.MarkerSymbol.Circle, StyleParam.SymbolStyle.Open);
                        Marker marker = Marker.init(Symbol: markerSymbol, Color: markerColor);

                        var theoreticalPlot = Chart.Point<double, double, string>(distributionInner.Masses,
                            distributionInner.Intensities.Select(z => z / distributionMax), Marker: marker);
                        GenericChartExtensions.WithTraceInfo(theoreticalPlot, "Theoretical Isotopic Distribution");
                        plotListCombined.Add(Chart.Combine(new[] { theoreticalPlot, plot }));
                    }

                    GenericChartExtensions.Show(Plotly.NET.Chart.Grid<IEnumerable<Plotly.NET.GenericChart.GenericChart>>(3, 1).Invoke(plotListCombined)
                            .WithSize(1000));

                    //var stackedPlot = Chart.Combine(plotList);
                    //stackedPlot.Show();
                }

                #endregion

                double enr = NormalizeArray(rejectionAveragedOutputs[j][1], 0.1, 1)
                    .CalculateEnr(scaleNoiseEstimates[j].refScaleEstimate, null,
                        scaleNoiseEstimates[j].refNoiseEstimate);
                double naiveEnr = NormalizeArray(naiveAveragedOutputs[j], 0.1, 1)
                    .CalculateEnr(scaleNoiseEstimates[j].refScaleEstimate,
                    null, scaleNoiseEstimates[j].refNoiseEstimate);

                enrs.Add((enr, naiveEnr));
                j++;
            }

            int maxPeaksIndexer = 0;
            int[] maxPeaks = Enumerable.Range(1, lfParamsDict.Count).ToArray();
            foreach (var valueTuple in enrs)
            {
                Console.WriteLine("{0}\t{1}\t{2}", maxPeaks[maxPeaksIndexer], valueTuple.rejection, valueTuple.noRejection);
                if (++maxPeaksIndexer == maxPeaks.Length)
                {
                    maxPeaksIndexer = 0;
                }
            }
        }


        [Test]

        private double[] NormalizeArray(double[] yarray, double normMin, double normMax)
        {
            double max = yarray.Max();
            double min = yarray.Min();
            double range = max - min;

            return yarray.Select(d => (d - min) / range)
                .Select(n => (1d - n) * normMin + n * normMax)
                .ToArray(); 
        }
        #region Demonstrate Noise Types
        [Test]
        public void DemonstrateHighFrequencyNoise()
        {
            

            double[] xarray = Enumerable.Range(0, 1000).Select(i => (double)i).ToArray();
            Random rnd = new Random(1551);
            double[] yarray = new double[xarray.Length];
            yarray = yarray.Select(_ => rnd.NextDouble()).ToArray();
            // yeah, this is an incredibly annoying way to generate the ylimits for plotly. Wish there was a better way, 
            // but I wasn't able to figure anything out. 
            var ylimits = new Optional<Tuple<double, double>>(new Tuple<double, double>(0.0, 2.0), true);
            GenericChartExtensions.Show(Chart.Line<double,double,string>(x: xarray, y: yarray)
                .WithYAxisStyle<double, double, string>(MinMax:ylimits)
                .WithTitle("Random, high frequency noise"));

        }

        [Test]
        public void DemonstrateLowFrequencyNoise()
        {
            Random rnd = new Random(1551);
            double[] xarray = Enumerable.Range(0, 1000).Select(i => (double)i).ToArray();
            double[] yarray = new double[xarray.Length];
            // add the high-frequency noise. 
            yarray = yarray.Select(_ => rnd.NextDouble()).ToArray();
            
            int numberOfRandomPeaks = 10;
            for (int i = 0; i < numberOfRandomPeaks; i++)
            {
                int indexOfRandomPeak = rnd.Next(1, 999); 
                yarray[indexOfRandomPeak] *= 10; 
            }
            var ylimits = new Optional<Tuple<double, double>>(new Tuple<double, double>(0, 10), true);
            GenericChartExtensions.Show(Chart.Line<double, double, string>(x: xarray, y: yarray)
                .WithYAxisStyle<double, double, string>(MinMax: ylimits)
                .WithTitle("Random, low frequency noise"));

        }
        #endregion
    }

    public enum SimulatedDataParamsType
    {
        SimulatedChargeStateEnvelope, 
        SimulatedProtein
    }
    public class SimulatedDataParams
    {
        public (double, double) MzRange { get; set; }
        public double StepSize { get; set; }
        public int Length { get; set; }
        public (int, int) ChargeStateRange { get; set; }
        public IsotopicDistribution ParentDistribution { get; set; }
        public (double mu, double sigma)? EnvelopeDistrParams { get; set; }
        public double AdductMass { get; set; }
        private SimulatedDataParamsType ParamsType { get; set; }
        protected SimulatedDataParams(SimulatedDataParamsType paramsType, (double, double) mzRange, double stepSize, int length,
            (int, int) chargeStateRange, IsotopicDistribution parentIsotopicDistribution,
            (double mu, double sigma)? envelopeDistribution = null, double adductMass = 1.007276)
        {
            MzRange = mzRange;
            StepSize = stepSize;
            Length = length;
            ChargeStateRange = chargeStateRange;
            ParentDistribution = parentIsotopicDistribution;
            EnvelopeDistrParams = envelopeDistribution;
            AdductMass = adductMass;
            ParamsType = paramsType; 
        }
    }
    public class SimulatedChargeStateEnvelopeParams : SimulatedDataParams
    {
        public SimulatedChargeStateEnvelopeParams(SimulatedDataParamsType paramsType, (double, double) mzRange, double stepSize, int length,
            (int, int) chargeStateRange, IsotopicDistribution parentIsotopicDistribution,
            (double mu, double sigma)? envelopeDistribution = null, double adductMass = 1.007276)
            : base(paramsType, mzRange, stepSize, length, chargeStateRange, parentIsotopicDistribution, envelopeDistribution,
                adductMass)
        {

        }

    }

}
