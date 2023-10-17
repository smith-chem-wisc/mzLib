using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;
using NUnit.Framework;
using UsefulProteomicsDatabases;

namespace Development.Deconvolution
{
    /// <summary>
    /// This class is designed as a playground for deconvolution developers.
    /// <remarks>
    /// Test Cases are instantiated in the static constructor found within the Set Up Test Cases Region
    /// New test cases can be added to the static constructor and will automatically be ran in all previously implemented development tests
    /// New testing methods can call the private test case collections as the TestCaseSource to iterate through all test cases
    /// </remarks>
    /// </summary>
    [TestFixture]
    [Ignore("Only needed when developing deconvolution methods")]
    [ExcludeFromCodeCoverage]
    public class StandardDeconvolutionTest
    {
        /// <summary>
        /// All Test Cases where a single peak is deconvoluted
        /// </summary>
        private static IEnumerable<SinglePeakDeconvolutionTestCase> _singlePeakTestCases;

        /// <summary>
        /// All Test Cases where an entire spectrum is deconvoluted
        /// </summary>
        private static IEnumerable<WholeSpectrumDeconvolutionTestCase> _wholeSpectrumDeconvolutionTestCases;

        #region Set Up Test Cases

        /// <summary>
        /// Instantiates all test cases to be ran
        /// </summary>
        static StandardDeconvolutionTest()
        {
            Loaders.LoadElements();

            // define paths to spectra
            var ubiquitinPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Deconvolution", "TestData",
                "Averaged_221110_UbiqOnly.mzML");
            var hghPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Deconvolution", "TestData",
                "Averaged_221110_HGHOnly.mzML");
            var cytoPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Deconvolution", "TestData",
                "Averaged_221110_CytoOnly.mzML");

            // set up deconvoluters to be utilized by test cases
            Deconvoluter classicTopDownDeconvoluter = new Deconvoluter(DeconvolutionType.ClassicDeconvolution,
                new ClassicDeconvolutionParameters(1, 60, 4, 3));
            Deconvoluter classicBottomUpDeconvoluter = new Deconvoluter(DeconvolutionType.ClassicDeconvolution,
                new ClassicDeconvolutionParameters(1, 12, 4, 3));

            // Add Individual peak test cases
            List<SinglePeakDeconvolutionTestCase> singlePeakDeconvolutionTestCases = new()
            {
                // uniquitin, direct injection
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10038.4, 8, 1254.8, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10039.41, 9, 1115.49, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10041.4, 10, 1004.14, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10041.46, 11, 912.86, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10043.4, 12, 836.95, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10043.41, 13, 772.57, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10044.44, 14, 717.46, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10045.5, 15, 669.70, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10045.44, 16, 627.84, 20),

                // hgh, direct injection
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22139.41, 11, 2012.29, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22136.28, 12, 1844.69, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22137.31, 13, 1702.87, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22139.32, 14, 1581.38, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22139.25, 15, 1475.95, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter,
                    "Direct Injectio Human Growth Hormone, Averaged",
                    hghPath, 1, 22140.32, 16, 1383.77, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22141.31, 17, 1302.43, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22142.34, 18, 1230.13, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22143.36, 19, 1165.44, 20),

                // cytochrome c, direct injection 
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12367.44, 9, 1374.16, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12367.4, 10, 1236.74, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12368.4, 11, 1124.40, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12370.44, 12, 1030.87, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12371.45, 13, 951.65, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12373.48, 14, 883.82, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12373.5, 15, 824.90, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12374.56, 16, 773.41, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12374.47, 17, 727.91, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12376.44, 18, 687.58, 20),
                new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12360.6, 20, 619.03, 20)
            };
            _singlePeakTestCases = singlePeakDeconvolutionTestCases;

            // Add whole spectrum test cases
            List<WholeSpectrumDeconvolutionTestCase> wholeSpectrumDeconvolutionTestCases = new()
            {
                new WholeSpectrumDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 20,
                    new[] { 10038.4, 10039.41, 10041.4, 10041.46, 10043.4, 10043.41, 10044.44, 10045.5, 10045.44, },
                    new[] { 8, 9, 10, 11, 12, 13, 14, 15, 16 },
                    new[] { 1254.8, 1115.49, 1004.14, 912.86, 836.95, 772.57, 717.46, 669.70, 627.84 }),

                new WholeSpectrumDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 20,
                    new []{22139.41, 22136.28, 22137.31, 22139.32, 22139.25, 22140.32, 22141.31, 22142.34, 22143.36},
                    new []{11, 12, 13, 14, 15, 16, 17, 18, 19},
                    new []{2012.29, 1844.69, 1702.87, 1581.38, 1475.95, 1383.77, 1302.43, 1230.13, 1165.44}),

                new WholeSpectrumDeconvolutionTestCase(classicTopDownDeconvoluter, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 20,
                    new []{12367.44, 12367.4, 12368.4, 12370.44, 12371.45, 12373.48, 12373.5, 12374.56, 12374.47, 12376.44, 12360.6},
                    new []{9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20},
                    new []{1374.16, 1236.74, 1124.40, 1030.87, 951.65, 883.82, 824.90, 773.41, 727.91, 687.58, 619.03}),
            };
            _wholeSpectrumDeconvolutionTestCases = wholeSpectrumDeconvolutionTestCases;
        }

        #endregion

        /// <summary>
        /// Tests fails if deconvolution top scoring result does not equal the expected
        /// </summary>
        /// <param name="testCase"></param>
        [Test]
        [TestCaseSource(nameof(_singlePeakTestCases))]
        public static void SinglePeak_TopScoringResult_IsCorrectWithinTolerance(SinglePeakDeconvolutionTestCase testCase)
        {
            // deconvolution
            List<IsotopicEnvelope> allResults = testCase.Deconvoluter
                .Deconvolute(testCase.SpectrumToDeconvolute, testCase.RangeToDeconvolute).ToList();

            // get top scoring result
            var topScoringResult = allResults.First();

            // test deconvolution results
            Assert.That(topScoringResult.Charge, Is.EqualTo(testCase.ExpectedIonChargeState));

            var acceptableDistanceFromTheoreticalWithinTestCaseTolerance = 
                testCase.DeconvolutionPPmTolerance.GetMaximumValue(testCase.ExpectedMostAbundantObservedIsotopicMass) -
                testCase.ExpectedMostAbundantObservedIsotopicMass;
            Assert.That(topScoringResult.MostAbundantObservedIsotopicMass,
                Is.EqualTo(testCase.ExpectedMostAbundantObservedIsotopicMass)
                    .Within(acceptableDistanceFromTheoreticalWithinTestCaseTolerance));

            acceptableDistanceFromTheoreticalWithinTestCaseTolerance =
                testCase.DeconvolutionPPmTolerance.GetMaximumValue(testCase.SelectedIonMz) -
                testCase.SelectedIonMz;
            Assert.That(topScoringResult.MostAbundantObservedIsotopicMass / topScoringResult.Charge,
                Is.EqualTo(testCase.SelectedIonMz).Within(acceptableDistanceFromTheoreticalWithinTestCaseTolerance));
        }

        /// <summary>
        /// Test fails if deconvolution does not have any result that matches the expected
        /// </summary>
        /// <param name="testCase"></param>
        [Test]
        [TestCaseSource(nameof(_singlePeakTestCases))]
        public static void SinglePeak_AnyResult_IsCorrectWithinTolerance(SinglePeakDeconvolutionTestCase testCase)
        {
            // deconvolution
            List<IsotopicEnvelope> allResults = testCase.Deconvoluter
                .Deconvolute(testCase.SpectrumToDeconvolute, testCase.RangeToDeconvolute).ToList();

            // extract tested properties from IsotopicEnvelopeResults
            var resultChargeStates = allResults.Select(p => p.Charge);
            var resultMostAbundantObservedIsotopicMass = allResults.Select(p => p.MostAbundantObservedIsotopicMass);
            var resultSelectedIonMz = allResults.Select(p => p.MostAbundantObservedIsotopicMass / p.Charge);

            // test deconvolution results
            Assert.That(resultChargeStates, Has.Some.EqualTo(testCase.ExpectedIonChargeState));

            var acceptableDistanceFromTheoreticalWithinTestCaseTolerance =
                testCase.DeconvolutionPPmTolerance.GetMaximumValue(testCase.ExpectedMostAbundantObservedIsotopicMass) -
                testCase.ExpectedMostAbundantObservedIsotopicMass;
            Assert.That(resultMostAbundantObservedIsotopicMass,
                Has.Some.EqualTo(testCase.ExpectedMostAbundantObservedIsotopicMass)
                    .Within(acceptableDistanceFromTheoreticalWithinTestCaseTolerance));

            acceptableDistanceFromTheoreticalWithinTestCaseTolerance =
                testCase.DeconvolutionPPmTolerance.GetMaximumValue(testCase.SelectedIonMz) -
                testCase.SelectedIonMz;
            Assert.That(resultSelectedIonMz,
                Has.Some.EqualTo(testCase.SelectedIonMz)
                    .Within(acceptableDistanceFromTheoreticalWithinTestCaseTolerance));
        }

        /// <summary>
        /// Test fails if deconvolution results do not contain the entire set of expected results
        /// </summary>
        /// <param name="testCase"></param>
        [Test]
        [TestCaseSource(nameof(_wholeSpectrumDeconvolutionTestCases))]
        public static void WholeSpectrum_ResultContainsAllExpected(WholeSpectrumDeconvolutionTestCase testCase)
        {
            // deconvolution
            List<IsotopicEnvelope> allResults = testCase.Deconvoluter
                .Deconvolute(testCase.SpectrumToDeconvolute).ToList();

            // extract tested properties from IsotopicEnvelopeResults
            var resultChargeStates = allResults.Select(p => p.Charge);
            var resultMostAbundantObservedIsotopicMass = allResults.Select(p => p.MostAbundantObservedIsotopicMass);
            var resultSelectedIonMz = allResults.Select(p => p.MostAbundantObservedIsotopicMass / p.Charge);

            // test deconvolution results
            for (int i = 0; i < testCase.Count; i++)
            {
                Assert.That(resultChargeStates, Has.Some.EqualTo(testCase.ExpectedIonChargeStates[i]));

                var acceptableDistanceFromTheoreticalWithinTestCaseTolerance =
                    testCase.DeconvolutionPPmTolerance.GetMaximumValue(testCase.ExpectedMostAbundantObservedIsotopicMasses[i]) -
                    testCase.ExpectedMostAbundantObservedIsotopicMasses[i];
                Assert.That(resultMostAbundantObservedIsotopicMass,
                    Has.Some.EqualTo(testCase.ExpectedMostAbundantObservedIsotopicMasses[i])
                        .Within(acceptableDistanceFromTheoreticalWithinTestCaseTolerance));

                acceptableDistanceFromTheoreticalWithinTestCaseTolerance =
                    testCase.DeconvolutionPPmTolerance.GetMaximumValue(testCase.SelectedIonMzs[i]) -
                    testCase.SelectedIonMzs[i];
                Assert.That(resultSelectedIonMz,
                    Has.Some.EqualTo(testCase.SelectedIonMzs[i])
                        .Within(acceptableDistanceFromTheoreticalWithinTestCaseTolerance));
            }
        }
    }
}


