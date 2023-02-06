using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using NUnit.Framework;
using MzLibUtil;
using System.Reflection;
using System.Globalization;
using System.Threading;
using Chemistry;

namespace Development.Deconvolution
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class StandardDeconvolutionTest
    {
        

        #region Set Up Test Cases


        private static IEnumerable<DeconvolutionTestCase> TestCases;

        // to set up new test case, add single scan to DeconvolutionTest.TestData and add new test case to object below
        // Can either add a single Deconvolution Test Case, or generate a list of them with numerous charges and m/zs using DeconvolutionTestCase.GenerateTestCases
        static StandardDeconvolutionTest()
        {
            List<DeconvolutionTestCase> cases = new();

            var ubiqCases = DeconvolutionTestCase.GenerateTestCases(SampleType.TopDown, "Direct Injection", "PolyUbiquitin, Averaged",
                10030.54, new[] { 8.0, 9, 10, 11, 12, 13, 14, 15, 16 }, new[] { 1254.8, 1115.49, 1004.14, 912.86, 836.95, 772.57, 717.46, 669.70, 627.84 }, "Averaged_221110_UbiqOnly.mzML");
            cases.AddRange(ubiqCases);

            // These test cases fail as the top scoring result is off by 24 Da, but the correct answer is still within the results of isotopic envelopes
            // This indicates that our scoring function needs some work
            var hghCasaes = DeconvolutionTestCase.GenerateTestCases(SampleType.TopDown, "Direct Injection", "Human Growth Hormone, Averaged",
                22124.41, new[] { 11.0, 12, 13, 14, 15, 16, 17, 18, 19 }, new[] { 2012.29, 1844.69, 1702.87, 1581.38, 1475.95, 1383.77, 1302.43, 1230.13, 1165.44 }, "Averaged_221110_HGHOnly.mzML");
            cases.AddRange(hghCasaes);

            var cytoCases = DeconvolutionTestCase.GenerateTestCases(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                   12358.47, new[] { 9.0, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20 }, new[] { 1374.16, 1236.74, 1124.40, 1030.87, 951.65, 883.82, 824.90, 773.41, 727.91, 687.58, 619.03 }, "Averaged_221110_CytoOnly.mzML");
            cases.AddRange(cytoCases);

            TestCases = cases;
        }

        

        #endregion

        [Test]
        [TestCaseSource(nameof(TestCases))]
        public static void TestClassicDeconvolution(DeconvolutionTestCase testCase)
        {
            // change max charge state based upon top down vs bottom up as is recommend in MetaMorpheus
            ClassicDeconvolutionParameters deconParameters = null;
            switch (testCase.SampleType)
            {
                case SampleType.TopDown:
                    deconParameters = new(1, 60, 4, 3);
                    break;
                case SampleType.BottomUp:
                    deconParameters = new(1, 12, 4, 3);
                    break;
            }

            // deconvolution
            Deconvoluter deconvoluter = new(DeconvolutionTypes.ClassicDeconvolution, deconParameters);
            List<IsotopicEnvelope> allResults = deconvoluter.Deconvolute(testCase.SpectrumToDeconvolute, testCase.RangeToDeconvolute).ToList();
            IsotopicEnvelope topScoringResult = allResults.First();

            // 3 was selected as it is the default mass difference acceptor in MetaMorpheus
            int withinMagicNumber = 3;
            Assert.That(topScoringResult.Charge, Is.EqualTo(testCase.SelectedIonChargeState));
            Assert.That(topScoringResult.MostAbundantObservedIsotopicMass, Is.EqualTo(testCase.MostAbundantMass).Within(withinMagicNumber));
            Assert.That((topScoringResult.MostAbundantObservedIsotopicMass - Constants.ProtonMass) / topScoringResult.Charge, Is.EqualTo(testCase.SelectedIonMz).Within(withinMagicNumber));
        }
    }
}


