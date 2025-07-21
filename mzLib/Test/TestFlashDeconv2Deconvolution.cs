using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Algorithms;
using MassSpectrometry.Deconvolution.Parameters;
using MzLibUtil;
using NUnit.Framework;
using Test.FileReadingTests;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class TestFlashDeconv2Deconvolution
    {
        private FlashDeconv2 _flashDeconv2;

        [SetUp]
        public void Setup()
        {
            var deconParams = new FlashDeconvDeconvolutionParameters(1,60);
            _flashDeconv2 = new FlashDeconv2(deconParams);
        }
        [Test]
        public void FlashDeconvWithArticialMs1Spectrum()
        {
            MsDataScan[] Scans = new MsDataScan[1];
            double selectedIonMz = 850.24;
            int selectedIonChargeStateGuess = 15;
            double selectedIonIntensity = 13.4;
            double isolationMz = 850.24; // This is the isolation m/z for the selected ion, which is the most intense proteoform in this test case.

            string Ms1SpectrumPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\artificialProteoform.txt");
            string[] spectrumLines = File.ReadAllLines(Ms1SpectrumPath);

            int mzIntensityPairsCount = spectrumLines.Length;
            double[] ms1mzs = new double[mzIntensityPairsCount];
            double[] ms1intensities = new double[mzIntensityPairsCount];

            for (int i = 0; i < mzIntensityPairsCount; i++)
            {
                string[] pair = spectrumLines[i].Split('\t');
                ms1mzs[i] = Convert.ToDouble(pair[0], CultureInfo.InvariantCulture);
                ms1intensities[i] = Convert.ToDouble(pair[1], CultureInfo.InvariantCulture);
            }

            MzSpectrum spectrum = new MzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MsDataScan(spectrum, 1, 1, false, Polarity.Positive, 1.0, new MzRange(740, 990), "first spectrum", MZAnalyzerType.Unknown, spectrum.SumOfAllY, null, null, null, selectedIonMz, selectedIonChargeStateGuess, selectedIonIntensity, isolationMz, 4);

            var myMsDataFile = new FakeMsDataFile(Scans);

            MsDataScan scan = myMsDataFile.GetAllScansList()[0];

            // The ones marked 2 are for checking an overload method

            DeconvolutionParameters deconParameters = new FlashDeconvDeconvolutionParameters(13, 17);

            FlashDeconv2 alg = new FlashDeconv2(deconParameters);
            List<IsotopicEnvelope> allMasses = alg.Deconvolute(scan.MassSpectrum, new MzRange((double)scan.MassSpectrum.FirstX, (double)scan.MassSpectrum.LastX)).ToList();

            Assert.That(allMasses.Count, Is.EqualTo(1));
            Assert.That(allMasses[0].Charge, Is.EqualTo(1));
            Assert.That(allMasses[0].MonoisotopicMass, Is.EqualTo(12730.5057551409).Within(0.01));
            Assert.That(allMasses[0].MostAbundantObservedIsotopicMass, Is.EqualTo(12738.53003).Within(0.01));
            Assert.That(allMasses[0].Score, Is.EqualTo(5111.138226).Within(0.01));
            Assert.That(allMasses[0].TotalIntensity, Is.EqualTo(499.45).Within(0.01));
        }
        [Test]
        public void FlashDeconvWithRealMs1Spectrum()
        {
            MsDataScan[] Scans = new MsDataScan[1];
            double selectedIonMz = 850.24;
            int selectedIonChargeStateGuess = 15;
            double selectedIonIntensity = 13.4;
            double isolationMz = 850.24; // This is the isolation m/z for the selected ion, which is the most intense proteoform in this test case.

            string Ms1SpectrumPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\realProteoform.txt");
            string[] spectrumLines = File.ReadAllLines(Ms1SpectrumPath);

            int mzIntensityPairsCount = spectrumLines.Length;
            double[] ms1mzs = new double[mzIntensityPairsCount];
            double[] ms1intensities = new double[mzIntensityPairsCount];

            for (int i = 0; i < mzIntensityPairsCount; i++)
            {
                string[] pair = spectrumLines[i].Split('\t');
                ms1mzs[i] = Convert.ToDouble(pair[0], CultureInfo.InvariantCulture);
                ms1intensities[i] = Convert.ToDouble(pair[1], CultureInfo.InvariantCulture);
            }

            MzSpectrum spectrum = new MzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MsDataScan(spectrum, 1, 1, false, Polarity.Positive, 1.0, new MzRange(740, 990), "first spectrum", MZAnalyzerType.Unknown, spectrum.SumOfAllY, null, null, null, selectedIonMz, selectedIonChargeStateGuess, selectedIonIntensity, isolationMz, 4);

            var myMsDataFile = new FakeMsDataFile(Scans);

            MsDataScan scan = myMsDataFile.GetAllScansList()[0];

            // The ones marked 2 are for checking an overload method

            DeconvolutionParameters deconParameters = new FlashDeconvDeconvolutionParameters(13, 17);

            FlashDeconv2 alg = new FlashDeconv2(deconParameters);
            List<IsotopicEnvelope> allMasses = alg.Deconvolute(scan.MassSpectrum, new MzRange((double)scan.MassSpectrum.FirstX, (double)scan.MassSpectrum.LastX)).ToList();

            Assert.That(allMasses.Count, Is.EqualTo(1));
            Assert.That(allMasses[0].Charge, Is.EqualTo(1));
            Assert.That(allMasses[0].MonoisotopicMass, Is.EqualTo(12730.499718360948).Within(0.01));
            Assert.That(allMasses[0].MostAbundantObservedIsotopicMass, Is.EqualTo(12738.532541223694).Within(0.01));
            Assert.That(allMasses[0].Score, Is.EqualTo(91500544022.111084).Within(0.01));
            Assert.That(allMasses[0].TotalIntensity, Is.EqualTo(1602985825.25).Within(0.01));
        }
        [Test]
        // Tests that LogTransformSpectrum filters out low-intensity peaks and applies log transform to X values.
        public void LogTransformSpectrum_FiltersAndTransformsCorrectly()
        {
            var x = new[] { 100.0, 200.0, 300.0 };
            var y = new[] { 0.005, 0.02, 0.5 };
            var spectrum = new MzSpectrum(x, y, false);

            var result = _flashDeconv2.GetType()
                .GetMethod("LogTransformSpectrum", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(_flashDeconv2, new object[] { spectrum, 0.01 }) as MzSpectrum;

            Assert.That(result.XArray.Length, Is.EqualTo(2));
            Assert.That(result.YArray.Length, Is.EqualTo(2));
            Assert.That(result.XArray, Is.All.GreaterThan(0));
        }



        [Test]
        // Tests that FilterMassIntensityGroupsByPpmTolerance separates groups into likely correct and incorrect based on ppm tolerance.
        public void FilterMassIntensityGroupsByPpmTolerance_FiltersCorrectly()
        {
            var groups = new List<(double[] neutralMass, double[] intensity)>
            {
                (new double[] { 1000, 1000.01 }, new double[] { 10, 20 }),
                (new double[] { 1000, 2000 }, new double[] { 10, 20 })
            };

            FlashDeconv2.FilterMassIntensityGroupsByPpmTolerance(
                groups,
                out var likelyCorrect,
                out var likelyIncorrect,
                correctPpmTolerance: 20_000, // Large tolerance to force first group as correct
                incorrectPpmTolerance: 10);

            Assert.That(likelyCorrect.Count, Is.EqualTo(1));
            Assert.That(likelyIncorrect.Count, Is.EqualTo(1));
        }

        [Test]
        // Tests that GetMostCommonNeutralMassAndSummedIntensity finds the mode cluster and sums the correct intensities.
        public void GetMostCommonNeutralMassAndSummedIntensity_ReturnsExpected()
        {
            var groups = new List<(double[] neutralMass, double[] intensity)>
            {
                (new double[] { 1000, 1000.01, 2000 }, new double[] { 10, 20, 30 })
            };

            var result = FlashDeconv2.GetMostCommonNeutralMassAndSummedIntensity(groups, ppmTolerance: 20);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].summedIntensity, Is.EqualTo(30)); // Only 2000 is outside 20 ppm, so 10+20=30
        }

        [Test]
        // Tests that CreateNeutralMassIntensityGroups creates a group with the expected number of neutral masses.
        public void CreateNeutralMassIntensityGroups_CreatesExpectedGroups()
        {
            var groups = new List<(double[] X, double[] Y, int[] ChargeState)>
            {
                (new double[] { Math.Log(100) }, new double[] { 10 }, new int[] { 1 })
            };

            var result = typeof(FlashDeconv2)
                .GetMethod("CreateNeutralMassIntensityGroups", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static)
                .Invoke(null, new object[] { groups }) as List<(double[] neutralMass, double[] intensity)>;

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].neutralMass.Length, Is.EqualTo(1));
        }

        [Test]
        // Tests that LogMzDependentTolerance returns a positive tolerance value for a given log(m/z).
        public void LogMzDependentTolerance_ReturnsExpectedTolerance()
        {
            double logMz = Math.Log(1000);
            var result = typeof(FlashDeconv2)
                .GetMethod("LogMzDependentTolerance", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static)
                .Invoke(null, new object[] { logMz, 250.0 });

            Assert.That(result, Is.TypeOf<double>());
            Assert.That((double)result, Is.GreaterThan(0));
        }

        [Test]
        // Tests that Deconvolute returns an empty result when given an empty spectrum.
        public void Deconvolute_ReturnsEmptyListForEmptySpectrum()
        {
            var spectrum = new MzSpectrum(new double[0], new double[0], false);
            var range = new MzRange(0, 1000);

            var result = _flashDeconv2.Deconvolute(spectrum, range);

            Assert.That(result, Is.Empty);
        }
    }
}
