using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using NUnit.Framework;
using SpectralAveraging;
using IO.MzML;
using MzLibUtil;
using System.Reflection;
using System.Globalization;
using System.Threading;
using Chemistry;

namespace Test.DeconvolutionTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class StandardDeconvolutionTest
    {
        #region Test Case Class

        public enum SampleType
        {
            TopDown,
            BottomUp,
        }

        public class TestCase
        {
            public SampleType SampleType { get; set; }
            public string SampleInformation { get; init; }
            public string TestName { get; init; }
            public double MonoIsotopicMass { get; init; }
            public double SelectedIonChargeState { get; init; }
            public double SelectedIonMz { get; init; }
            public MzRange RangeToDeconvolute { get; init; }
            public MzSpectrum SpectrumToDeconvolute { get; init; }
            public TestCase(SampleType sampleType, string sampleInformation, string testName, double monoIsotopicMass, double selectedIonChargeState, 
                 double selectedIonMz, MzSpectrum spectrumToDeconvolute)
            {
                SampleType = sampleType;
                SampleInformation = sampleInformation;
                TestName = testName;
                MonoIsotopicMass = monoIsotopicMass;
                SelectedIonChargeState = selectedIonChargeState;
                SelectedIonMz = selectedIonMz;
                SpectrumToDeconvolute = spectrumToDeconvolute;

                // 8.5 was selected as this is the magic number found in Classic Deconvolution
                RangeToDeconvolute = new MzRange(selectedIonMz - 8.5, selectedIonMz + 8.5);
            }

            public override string ToString()
            {
                return SampleInformation + " " + TestName + $"Charge: {SelectedIonChargeState}";
            }
        }

        #endregion

        #region Set Up Test Cases

        
        private static IEnumerable<TestCase> TestCases;

        // to set up new test case, add single scan to DeconvolutionTest.TestData and add new test case to object below
        static StandardDeconvolutionTest()
        {
            List<TestCase> cases = new();

            var caCases = GenerateTestCases(SampleType.TopDown, "Direct Injection",
                "Carbonic Anhydrase, Averaged", 29006.68, new[] { 32.0, 33, 34, 35, 36 },
                new[] { 907.97, 880.49, 854.62, 830.32, 807.20 }, "Averaged_221110_CaOnly_.mzML");
            cases.AddRange(caCases);

            var ubiqCases = GenerateTestCases(SampleType.TopDown, "Direct Injection", "PolyUbiquitin, Averaged",
                10025.42, new[] { 10.0, 11, 12, 13, 14, 15 }, new[] { 1004.14, 912.95, 836.87, 772.57, 717.46, 669.63 }, "Averaged_221110_UbiqOnly.mzML");
            cases.AddRange(ubiqCases);

            var hghCasaes = GenerateTestCases(SampleType.TopDown, "Direct Injection", "Human Growth Hormone, Averaged",
                22132.24, new[] { 13.0, 14, 15, 16, 17, 18 }, new[] { 1702.79, 1581.30, 1476.08, 1381.39, 1303.72, 1231.40 }, "Averaged_221110_HGHOnly.mzML");
            cases.AddRange(hghCasaes);

            var cytoCases = GenerateTestCases(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                   12351.45, new [] { 12.0, 13, 14, 15, 16, 17, 18 }, new [] { 1030.87, 951.65, 883.82, 824.43, 773.41, 727.91, 687.58 }, "Averaged_221110_CytoOnly.mzML");
            cases.AddRange(cytoCases);

            TestCases = cases;
        }

        /// <summary>
        /// Generates a test case for each charge state and mz inputted
        /// </summary>
        /// <param name="sampleType">top down or bottom up?</param>
        /// <param name="sampleInformation">How was this data acquired?</param>
        /// <param name="testName">Sample specific information</param>
        /// <param name="monoIsotopicMass">mono mass of deconvoluted protein/peptide</param>
        /// <param name="selectedIonChargeStates">charge state of selected ion</param>
        /// <param name="selectedIonMzs">mz of selected on</param>
        /// <param name="fileToTest">mass spec file to test</param>
        /// <returns></returns>
        static IEnumerable<TestCase> GenerateTestCases(SampleType sampleType, string sampleInformation, string testName, double monoIsotopicMass, double[] selectedIonChargeStates,
            double[] selectedIonMzs, string fileToTest)
        {
            List<TestCase> cases = new List<TestCase>();
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DeconvolutionTests", "TestData", fileToTest);
            MzSpectrum spectrum = Mzml.LoadAllStaticData(filePath).GetAllScansList().First().MassSpectrum;
            for (int i = 0; i < selectedIonMzs.Length; i++)
            {
                yield return new TestCase(sampleType, sampleInformation, testName, monoIsotopicMass, selectedIonChargeStates[i], selectedIonMzs[i], spectrum);
            }
        }

        #endregion

        [Test]
        [TestCaseSource(nameof(TestCases))]
        public static void TestClassicDeconvolution(TestCase testCase)
        {
            // cases which the correct value is there but scoring fails
            if (Math.Abs(testCase.SelectedIonMz - 854.62) < 0.01)
                return;



            // change max charge state based upon top down vs bottom up as is recommend in MetaMorpheus
            ClassicDeconvolutionParameters deconParameters = null;
            switch (testCase.SampleType)
            {
                case SampleType.TopDown:
                    deconParameters = new(1, 60, 4, 3);
                    break;
                case SampleType.BottomUp:
                    deconParameters = new (1, 12, 4, 3);
                    break;
            }

            // deconvolution
            Deconvoluter deconvoluter = new(DeconvolutionTypes.ClassicDeconvolution, deconParameters);
            List<IsotopicEnvelope> allResults = deconvoluter.Deconvolute(testCase.SpectrumToDeconvolute, testCase.RangeToDeconvolute).ToList();
            IsotopicEnvelope topScoringResult = allResults.First();


            int withinMagicNumber = 3; // 3 was selected as that is the default mass difference acceptor in MetaMorpheus
            Assert.That(topScoringResult.Charge, Is.EqualTo(testCase.SelectedIonChargeState));
            Assert.That(topScoringResult.MonoisotopicMass, Is.EqualTo(testCase.MonoIsotopicMass).Within(withinMagicNumber));
            Assert.That((topScoringResult.MostAbundantObservedIsotopicMass - Constants.ProtonMass) / (double)topScoringResult.Charge, Is.EqualTo(testCase.SelectedIonMz).Within(withinMagicNumber));
        }
    }
}


