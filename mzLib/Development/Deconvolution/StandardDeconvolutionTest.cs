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

            string ubiquitinPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Deconvolution", "TestData", "Averaged_221110_UbiqOnly.mzML");
            string hghPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Deconvolution", "TestData", "Averaged_221110_HGHOnly.mzML");
            string cytoPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Deconvolution", "TestData", "Averaged_221110_CytoOnly.mzML");
            List<DeconvolutionTestCase> cases = new()
            {
                // uniquitin, direct injection
                new(SampleType.TopDown, "Direct Injection", "PolyUbiquitin, Averaged",
                    ubiquitinPath, 10038.4, 8, 1254.8, 20),
                new(SampleType.TopDown, "Direct Injection", "PolyUbiquitin, Averaged",
                    ubiquitinPath, 10039.41, 9, 1115.49, 20),
                new(SampleType.TopDown, "Direct Injection", "PolyUbiquitin, Averaged",
                    ubiquitinPath, 10041.4, 10, 1004.14, 20),
                new(SampleType.TopDown, "Direct Injection", "PolyUbiquitin, Averaged",
                    ubiquitinPath, 10041.46, 11, 912.86, 20),
                new(SampleType.TopDown, "Direct Injection", "PolyUbiquitin, Averaged",
                    ubiquitinPath, 10043.4, 12, 836.95, 20),
                new(SampleType.TopDown, "Direct Injection", "PolyUbiquitin, Averaged",
                    ubiquitinPath, 10043.41, 13, 772.57, 20),
                new(SampleType.TopDown, "Direct Injection", "PolyUbiquitin, Averaged",
                    ubiquitinPath, 10044.44, 14, 717.46, 20),
                new(SampleType.TopDown, "Direct Injection", "PolyUbiquitin, Averaged",
                    ubiquitinPath, 10045.5, 15, 669.70, 20),
                new(SampleType.TopDown, "Direct Injection", "PolyUbiquitin, Averaged",
                    ubiquitinPath, 10045.44, 16, 627.84, 20),

                // hgh, direct injection
                // These test cases fail as the top scoring result is off by 24 Da, but the correct answer is still within the results of isotopic envelopes
                // This indicates that our scoring function needs some work
                new(SampleType.TopDown, "Direct Injection", "Human Growth Hormone, Averaged",
                    hghPath, 22139.41, 11, 2012.29, 20),
                new(SampleType.TopDown, "Direct Injection", "Human Growth Hormone, Averaged",
                    hghPath, 22136.28, 12, 1844.69, 20),
                new(SampleType.TopDown, "Direct Injection", "Human Growth Hormone, Averaged",
                    hghPath, 22137.31, 13, 1702.87, 20),
                new(SampleType.TopDown, "Direct Injection", "Human Growth Hormone, Averaged",
                    hghPath, 22139.32, 14, 1581.38, 20),
                new(SampleType.TopDown, "Direct Injection", "Human Growth Hormone, Averaged",
                    hghPath, 22139.25, 15, 1475.95, 20),
                new(SampleType.TopDown, "Direct Injection", "Human Growth Hormone, Averaged",
                    hghPath, 22140.32, 16, 1383.77, 20),
                new(SampleType.TopDown, "Direct Injection", "Human Growth Hormone, Averaged",
                    hghPath, 22141.31, 17, 1302.43, 20),
                new(SampleType.TopDown, "Direct Injection", "Human Growth Hormone, Averaged",
                    hghPath, 22142.34, 18, 1230.13, 20),
                new(SampleType.TopDown, "Direct Injection", "Human Growth Hormone, Averaged",
                    hghPath, 22143.36, 19, 1165.44, 20),

                // cytochrome c, direct injection 
                new(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                    cytoPath, 12367.44, 9, 1374.16, 20),
                new(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                    cytoPath, 12367.4, 10, 1236.74, 20),
                new(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                    cytoPath, 12368.4, 11, 1124.40, 20),
                new(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                    cytoPath, 12370.44, 12, 1030.87, 20),
                new(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                    cytoPath, 12371.45, 13, 951.65, 20),
                new(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                    cytoPath, 12373.48, 14, 883.82, 20),
                new(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                    cytoPath, 12373.5, 15, 824.90, 20),
                new(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                    cytoPath, 12374.56, 16, 773.41, 20),
                new(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                    cytoPath, 12374.47, 17, 727.91, 20),
                new(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                    cytoPath, 12376.44, 18, 687.58, 20),
                new(SampleType.TopDown, "Direct Injection", "Cytochrome C, Averaged",
                    cytoPath, 12360.6, 20, 619.03, 20)
            };

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

            Assert.That(topScoringResult.Charge, Is.EqualTo(testCase.SelectedIonChargeState));

            var acceptableDistanceFromTheoreticalWithinTestCaseTolerance =
                testCase.DeconvolutionPPmTolerance.GetMaximumValue(testCase.MostAbundantMass) -
                testCase.MostAbundantMass;
            Assert.That(topScoringResult.MostAbundantObservedIsotopicMass, Is.EqualTo(testCase.MostAbundantMass).Within(acceptableDistanceFromTheoreticalWithinTestCaseTolerance));

            acceptableDistanceFromTheoreticalWithinTestCaseTolerance =
                testCase.DeconvolutionPPmTolerance.GetMaximumValue(testCase.SelectedIonMz) -
                testCase.SelectedIonMz;
            Assert.That(topScoringResult.MostAbundantObservedIsotopicMass / topScoringResult.Charge, Is.EqualTo(testCase.SelectedIonMz).Within(acceptableDistanceFromTheoreticalWithinTestCaseTolerance));
        }
    }
}


