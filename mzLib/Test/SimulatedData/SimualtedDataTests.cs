using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;
using FlashLFQ;
using MathNet.Numerics.Distributions;
using NUnit;
using NUnit.Framework;
using Plotly.NET;
using Proteomics.AminoAcidPolymer;
using SimulatedData;
using SpectralAveraging;
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
            gaussianPeak.Plot().Show();
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
            gaussianPeak.Plot().Show();

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
            gaussianPeak.Plot().Show();
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
            gaussianPeak.Plot().Show();
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
            simProtein.Plot().Show();
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
            simProtein.Plot().Show(); 
        }

        [Test]
        [TestCase()]
        public void TestGetNoiseEstimation()
        {

        }
    }
}
