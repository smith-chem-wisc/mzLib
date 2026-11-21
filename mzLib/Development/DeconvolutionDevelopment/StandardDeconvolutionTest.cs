using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
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
            // define paths to spectra
            var ubiquitinPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DeconvolutionDevelopment", "TestData",
                "Averaged_221110_UbiqOnly.mzML");
            var hghPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DeconvolutionDevelopment", "TestData",
                "Averaged_221110_HGHOnly.mzML");
            var cytoPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DeconvolutionDevelopment", "TestData",
                "Averaged_221110_CytoOnly.mzML");

            // Resolve the FLASHDeconv executable once. The wrapper algorithm needs
            // an explicit path on its parameters; if FLASHDeconv isn't installed
            // on the dev machine, skip its test cases instead of crashing every
            // test in the fixture with "FLASHDeconvExePath is not set".
            string flashDeconvExePath = null;
            try { flashDeconvExePath = FlashDeconvExePathRegistry.Resolve(); }
            catch { /* FLASHDeconv not installed -- skip those cases below */ }

            // set up deconvoluters to be utilized by test cases
            List<DeconvolutionParameters> topDownDeconvolutionParametersToTest =
            [
                new ClassicDeconvolutionParameters(1, 60, 4, 3),
                new IsoDecDeconvolutionParameters(),
            ];

            List<DeconvolutionParameters> bottomUpDeconvolutionParametersToTest =
            [
                new ClassicDeconvolutionParameters(1, 12, 4, 3),
                new IsoDecDeconvolutionParameters(),
            ];

            if (flashDeconvExePath != null)
            {
                topDownDeconvolutionParametersToTest.Add(
                    new RealFLASHDeconvolutionParameters(1, 60, tolerancePpm: 4, flashDeconvExePath: flashDeconvExePath));
                bottomUpDeconvolutionParametersToTest.Add(
                    new RealFLASHDeconvolutionParameters(1, 12, tolerancePpm: 4, flashDeconvExePath: flashDeconvExePath));
            }



            // Add Individual peak test cases for top down
            List<SinglePeakDeconvolutionTestCase> singlePeakDeconvolutionTestCases = new();
            foreach (var deconParams in topDownDeconvolutionParametersToTest)
            {
                // uniquitin, direct injection
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10038.4, 8, 1254.8, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10039.41, 9, 1115.49, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10041.4, 10, 1004.14, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10041.46, 11, 912.86, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10043.4, 12, 836.95, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10043.41, 13, 772.57, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10044.44, 14, 717.46, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10045.5, 15, 669.70, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 10045.44, 16, 627.84, 20));

                // hgh, direct injection
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22139.41, 11, 2012.29, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22136.28, 12, 1844.69, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22137.31, 13, 1702.87, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22139.32, 14, 1581.38, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22139.25, 15, 1475.95, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injectio Human Growth Hormone, Averaged",
                    hghPath, 1, 22140.32, 16, 1383.77, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22141.31, 17, 1302.43, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22142.34, 18, 1230.13, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 22143.36, 19, 1165.44, 20));

                // cytochrome c, direct injection 
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12367.44, 9, 1374.16, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12367.4, 10, 1236.74, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12368.4, 11, 1124.40, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12370.44, 12, 1030.87, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12371.45, 13, 951.65, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12373.48, 14, 883.82, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12373.5, 15, 824.90, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12374.56, 16, 773.41, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12374.47, 17, 727.91, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12376.44, 18, 687.58, 20));
                singlePeakDeconvolutionTestCases.Add(new SinglePeakDeconvolutionTestCase(deconParams,
                    "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 12360.6, 20, 619.03, 20));
            }
            _singlePeakTestCases = singlePeakDeconvolutionTestCases;

            // Add whole spectrum test cases for top down
            List<WholeSpectrumDeconvolutionTestCase> wholeSpectrumDeconvolutionTestCases = new();
            foreach (var deconParams in topDownDeconvolutionParametersToTest)
            {
                wholeSpectrumDeconvolutionTestCases.Add(new WholeSpectrumDeconvolutionTestCase(deconParams, "Direct Injection PolyUbiquitin, Averaged",
                    ubiquitinPath, 1, 20,
                    new[] { 10038.4, 10039.41, 10041.4, 10041.46, 10043.4, 10043.41, 10044.44, 10045.5, 10045.44, },
                    new[] { 8, 9, 10, 11, 12, 13, 14, 15, 16 },
                    new[] { 1254.8, 1115.49, 1004.14, 912.86, 836.95, 772.57, 717.46, 669.70, 627.84 }));

                wholeSpectrumDeconvolutionTestCases.Add(new WholeSpectrumDeconvolutionTestCase(deconParams, "Direct Injection Human Growth Hormone, Averaged",
                    hghPath, 1, 20,
                    new []{22139.41, 22136.28, 22137.31, 22139.32, 22139.25, 22140.32, 22141.31, 22142.34, 22143.36},
                    new []{11, 12, 13, 14, 15, 16, 17, 18, 19},
                    new []{2012.29, 1844.69, 1702.87, 1581.38, 1475.95, 1383.77, 1302.43, 1230.13, 1165.44}));

                wholeSpectrumDeconvolutionTestCases.Add(new WholeSpectrumDeconvolutionTestCase(deconParams, "Direct Injection Cytochrome C, Averaged",
                    cytoPath, 1, 20,
                    new []{12367.44, 12367.4, 12368.4, 12370.44, 12371.45, 12373.48, 12373.5, 12374.56, 12374.47, 12376.44, 12360.6},
                    new []{9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20},
                    new []{1374.16, 1236.74, 1124.40, 1030.87, 951.65, 883.82, 824.90, 773.41, 727.91, 687.58, 619.03}));
            };
            _wholeSpectrumDeconvolutionTestCases = wholeSpectrumDeconvolutionTestCases;

            // TODO: Add cases for bottom up deconvolution
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
            List<IsotopicEnvelope> allResults = Deconvoluter
                .Deconvolute(testCase.SpectrumToDeconvolute, testCase.DeconvolutionParameters, testCase.RangeToDeconvolute).ToList();

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
            List<IsotopicEnvelope> allResults = Deconvoluter
                .Deconvolute(testCase.SpectrumToDeconvolute, testCase.DeconvolutionParameters, testCase.RangeToDeconvolute).ToList();

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
            List<IsotopicEnvelope> allResults = Deconvoluter
                .Deconvolute(testCase.SpectrumToDeconvolute, testCase.DeconvolutionParameters).ToList();

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

        /// <summary>
        /// Diagnostic test. Runs each available deconvolution algorithm against a
        /// multi-scan top-down yeast LC-MS snippet (lvsYeastTopDownSnip.mzML) and
        /// reports the total number of MS1 envelopes each produces, written to the
        /// test log.
        ///
        /// FLASHDeconv is exercised in two modes:
        ///   - Per-scan (what Deconvoluter.Deconvolute uses): one in-memory MzSpectrum
        ///     at a time. Documented to be sparse because no multi-scan feature
        ///     tracing happens.
        ///   - Whole-file (RealFLASHDeconvolutionAlgorithm.DeconvoluteFile): FLASHDeconv's
        ///     intended best mode. The numbers from this row are the meaningful ones
        ///     for comparing FLASHDeconv against Classic/IsoDec on the same input.
        ///
        /// FLASHDeconv rows are skipped silently if the exe isn't installed on the
        /// dev machine. Each printed row gets a > 0 sanity check; the headline output
        /// is the count itself.
        /// </summary>
        [Test]
        public static void TopDownSnip_ReportEnvelopeCounts()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "DeconvolutionDevelopment", "TestData", "lvsYeastTopDownSnip.mzML");
            Assume.That(File.Exists(path), $"Top-down snippet not found at {path}");

            var dataFile = MsDataFileReader.GetDataFile(path).LoadAllStaticData();
            var ms1Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToList();
            TestContext.Out.WriteLine($"Input: {Path.GetFileName(path)} -- {ms1Scans.Count} MS1 scan(s)");

            string flashExe = null;
            try { flashExe = FlashDeconvExePathRegistry.Resolve(); }
            catch { /* not installed -- skip those rows */ }

            // Algorithms invoked per-scan via the Deconvoluter factory.
            var perScan = new List<(string Name, DeconvolutionParameters Params)>
            {
                ("Classic top-down (1-60, 4ppm)", new ClassicDeconvolutionParameters(1, 60, 4, 3)),
                ("IsoDec (defaults)",             new IsoDecDeconvolutionParameters()),
            };
            if (flashExe != null)
                perScan.Add(("FLASHDeconv per-scan (1-60, 4ppm)",
                             new RealFLASHDeconvolutionParameters(1, 60, tolerancePpm: 4, flashDeconvExePath: flashExe)));

            foreach (var (name, prms) in perScan)
            {
                int total = ms1Scans.Sum(scan =>
                    Deconvoluter.Deconvolute(scan.MassSpectrum, prms).Count());
                TestContext.Out.WriteLine($"  {name,-40} : {total} envelope(s) across {ms1Scans.Count} MS1 scan(s)");
                Assert.That(total, Is.GreaterThan(0), $"{name} produced zero envelopes; check params.");
            }

            // FLASHDeconv's intended whole-file mode: bypasses the per-scan factory
            // and the WriteSingleScanMzml round-trip, lets FLASHDeconv do multi-scan
            // feature tracing on the original mzML. This is the apples-to-apples
            // row for "how does FLASHDeconv really do on this data."
            if (flashExe != null)
            {
                var p = new RealFLASHDeconvolutionParameters(1, 60, tolerancePpm: 4, flashDeconvExePath: flashExe);
                var byScan = RealFLASHDeconvolutionAlgorithm.DeconvoluteFile(path, p);
                int total = byScan.Values.Sum(list => list.Count);
                TestContext.Out.WriteLine($"  {"FLASHDeconv whole-file (1-60, 4ppm)",-40} : {total} envelope(s) across {byScan.Count} scan(s)");
                Assert.That(total, Is.GreaterThan(0), "FLASHDeconv whole-file mode produced zero envelopes; check params.");
            }
        }

        /// <summary>
        /// Diagnostic test. Asks "do the three algorithms agree on which monoisotopic
        /// masses are present in the top-down snippet?" by:
        ///   1. Running each algorithm on every MS1 scan.
        ///   2. Collecting the union of envelope monoisotopic masses, tagged by algorithm.
        ///   3. Clustering them within 10 ppm so that "same mass" (within mass-accuracy
        ///      tolerance) maps to one bucket regardless of which algorithm reported it.
        ///   4. Tallying: how many clusters were found by all 3 algorithms, by 2, by 1.
        ///
        /// Pure reporting -- the only assertion is that there is at least one cluster
        /// (sanity that something deconvolved). Output goes to the test log.
        ///
        /// Skipped quietly if FLASHDeconv isn't installed; otherwise the same-scan
        /// envelopes from per-scan and whole-file FLASHDeconv modes are merged into
        /// one set since they target the same algorithm's idea of "what's there."
        /// </summary>
        [Test]
        public static void TopDownSnip_ConsensusAcrossAlgorithms()
        {
            const double ClusterPpm = 10.0;

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "DeconvolutionDevelopment", "TestData", "lvsYeastTopDownSnip.mzML");
            Assume.That(File.Exists(path), $"Top-down snippet not found at {path}");

            var dataFile = MsDataFileReader.GetDataFile(path).LoadAllStaticData();
            var ms1Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToList();

            string flashExe = null;
            try { flashExe = FlashDeconvExePathRegistry.Resolve(); } catch { }

            // (algorithm-label, masses-observed-anywhere-in-the-file)
            var observed = new List<(string Algo, List<double> Masses)>();

            var classic = new ClassicDeconvolutionParameters(1, 60, 4, 3);
            observed.Add(("Classic", ms1Scans
                .SelectMany(s => Deconvoluter.Deconvolute(s.MassSpectrum, classic))
                .Select(e => e.MonoisotopicMass).ToList()));

            var isoDec = new IsoDecDeconvolutionParameters();
            observed.Add(("IsoDec", ms1Scans
                .SelectMany(s => Deconvoluter.Deconvolute(s.MassSpectrum, isoDec))
                .Select(e => e.MonoisotopicMass).ToList()));

            if (flashExe != null)
            {
                var fd = new RealFLASHDeconvolutionParameters(1, 60, tolerancePpm: 4, flashDeconvExePath: flashExe);
                var byScan = RealFLASHDeconvolutionAlgorithm.DeconvoluteFile(path, fd);
                observed.Add(("FLASHDeconv", byScan.Values
                    .SelectMany(list => list)
                    .Select(e => e.MonoisotopicMass).ToList()));
            }

            // Flatten to one big sorted list tagged with algorithm.
            var tagged = observed
                .SelectMany(o => o.Masses.Select(m => (Mass: m, Algo: o.Algo)))
                .OrderBy(x => x.Mass)
                .ToList();

            // Greedy single-pass clustering: a point joins the current cluster if it
            // is within ClusterPpm of the cluster's running max mass.
            var clusters = new List<List<(double Mass, string Algo)>>();
            foreach (var p in tagged)
            {
                if (clusters.Count == 0)
                {
                    clusters.Add(new List<(double, string)> { p });
                    continue;
                }
                var last = clusters[^1];
                double anchor = last[^1].Mass;
                double delta = anchor * ClusterPpm * 1e-6;
                if (p.Mass - anchor <= delta)
                    last.Add(p);
                else
                    clusters.Add(new List<(double, string)> { p });
            }

            // Tally consensus: how many distinct algorithms contributed to each cluster.
            int totalAlgos = observed.Count;
            var perLevel = Enumerable.Range(1, totalAlgos).ToDictionary(k => k, _ => 0);
            var privateTo = observed.ToDictionary(o => o.Algo, _ => 0);
            foreach (var c in clusters)
            {
                var algos = c.Select(x => x.Algo).Distinct().ToList();
                perLevel[algos.Count] += 1;
                if (algos.Count == 1)
                    privateTo[algos[0]] += 1;
            }

            TestContext.Out.WriteLine($"Input: {Path.GetFileName(path)} -- {ms1Scans.Count} MS1 scan(s), {totalAlgos} algorithm(s), {ClusterPpm}ppm clustering");
            TestContext.Out.WriteLine($"Total distinct mass clusters union: {clusters.Count}");
            for (int k = totalAlgos; k >= 1; k--)
                TestContext.Out.WriteLine($"  Found by {k} algorithm(s) : {perLevel[k],6}");
            TestContext.Out.WriteLine($"Private (found by only one algorithm):");
            foreach (var (algo, count) in privateTo.OrderByDescending(kv => kv.Value))
                TestContext.Out.WriteLine($"  {algo,-12} : {count,6}");

            // Per-algorithm reach: how many of the all-clusters did this algorithm hit?
            TestContext.Out.WriteLine($"Per-algorithm cluster reach:");
            foreach (var (algo, _) in observed)
            {
                int hit = clusters.Count(c => c.Any(x => x.Algo == algo));
                double pct = 100.0 * hit / clusters.Count;
                TestContext.Out.WriteLine($"  {algo,-12} : {hit,6} / {clusters.Count} ({pct:F1}%)");
            }

            Assert.That(clusters.Count, Is.GreaterThan(0), "No clusters formed -- nothing deconvolved.");
        }

        // ── Consensus-recall infrastructure ──────────────────────────────────
        //
        // The idea: define a fixed reference set of established deconvolution
        // algorithms whose intersection (within a ppm tolerance) defines the
        // "true" mass list for lvsYeastTopDownSnip.mzML. Then run every algorithm
        // -- reference and new -- and measure how many of those consensus masses
        // it recovers. Reference algorithms must recover 100% (they're the ones
        // that defined the consensus). New algorithms are measured against the
        // fixed bar.
        //
        // To add a new algorithm to the comparison: add one entry to
        // _algorithmsToTest below, with IsReference = false and a sensible
        // MinRecallPct threshold. No other code changes required.

        private const double ConsensusClusterPpm = 10.0;

        private sealed record AlgorithmRecallEntry(
            string Label,
            Func<string, DeconvolutionParameters> ParamFactory,
            double MinRecallPct,
            bool IsReference);

        private static readonly List<AlgorithmRecallEntry> _algorithmsToTest =
        [
            new(Label: "Classic top-down (1-60, 4ppm)",
                ParamFactory: _ => new ClassicDeconvolutionParameters(1, 60, 4, 3),
                MinRecallPct: 100.0,
                IsReference: true),
            new(Label: "IsoDec (defaults)",
                ParamFactory: _ => new IsoDecDeconvolutionParameters(),
                MinRecallPct: 100.0,
                IsReference: true),
            new(Label: "FLASHDeconv whole-file (1-60, 4ppm)",
                ParamFactory: exe => new RealFLASHDeconvolutionParameters(1, 60, tolerancePpm: 4, flashDeconvExePath: exe),
                MinRecallPct: 100.0,
                IsReference: true),
            // Add new algorithms here. Example shape:
            // new(Label: "MyNewAlgo top-down",
            //     ParamFactory: _ => new MyNewAlgoParameters(...),
            //     MinRecallPct: 80.0,
            //     IsReference: false),
        ];

        /// <summary>
        /// Establishes a consensus mass list on lvsYeastTopDownSnip.mzML from the
        /// reference algorithms and asserts each registered algorithm recovers at
        /// least its declared MinRecallPct of that consensus. Skipped quietly if
        /// FLASHDeconv isn't installed (the consensus requires all references).
        /// </summary>
        [Test]
        public static void TopDownSnip_ConsensusRecall()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "DeconvolutionDevelopment", "TestData", "lvsYeastTopDownSnip.mzML");
            Assume.That(File.Exists(path), $"Top-down snippet not found at {path}");

            string flashExe = null;
            try { flashExe = FlashDeconvExePathRegistry.Resolve(); } catch { }
            Assume.That(flashExe, Is.Not.Null,
                "FLASHDeconv not installed on this machine; consensus requires all reference algorithms.");

            var dataFile = MsDataFileReader.GetDataFile(path).LoadAllStaticData();
            var ms1Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1).ToList();
            Assume.That(ms1Scans, Is.Not.Empty, "No MS1 scans in input file.");

            // Pre-run every registered algorithm once, collect the masses they
            // produce anywhere in the file. Cached so we don't repeat the work
            // when computing both the consensus AND each algorithm's recall.
            var massesByAlgo = new Dictionary<string, List<double>>();
            foreach (var entry in _algorithmsToTest)
            {
                var prms = entry.ParamFactory(flashExe);
                List<double> masses = prms is RealFLASHDeconvolutionParameters fdParams
                    ? RealFLASHDeconvolutionAlgorithm.DeconvoluteFile(path, fdParams)
                        .Values.SelectMany(list => list).Select(e => e.MonoisotopicMass).ToList()
                    : ms1Scans.SelectMany(s => Deconvoluter.Deconvolute(s.MassSpectrum, prms))
                        .Select(e => e.MonoisotopicMass).ToList();
                massesByAlgo[entry.Label] = masses;
            }

            // Build consensus from reference algorithms only.
            var referenceLabels = _algorithmsToTest.Where(a => a.IsReference).Select(a => a.Label).ToList();
            var consensus = BuildConsensusMasses(massesByAlgo, referenceLabels, ConsensusClusterPpm);

            TestContext.Out.WriteLine($"Input: {Path.GetFileName(path)} -- {ms1Scans.Count} MS1 scan(s)");
            TestContext.Out.WriteLine($"Reference algorithms ({referenceLabels.Count}): {string.Join(", ", referenceLabels)}");
            TestContext.Out.WriteLine($"Consensus masses (all references agree within {ConsensusClusterPpm}ppm): {consensus.Count}");
            TestContext.Out.WriteLine($"Recall per algorithm:");

            Assert.That(consensus, Is.Not.Empty,
                "Empty consensus -- reference algorithms found no common masses; cannot evaluate.");

            Assert.Multiple(() =>
            {
                foreach (var entry in _algorithmsToTest)
                {
                    int matched = CountMassesWithinPpm(consensus, massesByAlgo[entry.Label], ConsensusClusterPpm);
                    double recallPct = 100.0 * matched / consensus.Count;
                    string flag = entry.IsReference ? "[ref]" : "[new]";
                    TestContext.Out.WriteLine(
                        $"  {flag} {entry.Label,-40} : {matched,5} / {consensus.Count} = {recallPct,6:F2}%   (min {entry.MinRecallPct:F1}%)");
                    Assert.That(recallPct, Is.GreaterThanOrEqualTo(entry.MinRecallPct),
                        $"{entry.Label} recall {recallPct:F2}% is below its declared minimum {entry.MinRecallPct:F1}%.");
                }
            });
        }

        /// <summary>
        /// Returns the set of monoisotopic masses where every algorithm in
        /// <paramref name="referenceLabels"/> contributed at least one mass within
        /// <paramref name="ppmTolerance"/> of the cluster. Cluster representative
        /// is the mean of the contributing masses.
        /// </summary>
        private static List<double> BuildConsensusMasses(
            Dictionary<string, List<double>> massesByAlgo,
            List<string> referenceLabels,
            double ppmTolerance)
        {
            var tagged = referenceLabels
                .SelectMany(label => massesByAlgo[label].Select(m => (Mass: m, Algo: label)))
                .OrderBy(x => x.Mass)
                .ToList();
            if (tagged.Count == 0) return new List<double>();

            var clusters = new List<List<(double Mass, string Algo)>>();
            foreach (var p in tagged)
            {
                if (clusters.Count == 0)
                {
                    clusters.Add(new List<(double, string)> { p });
                    continue;
                }
                var last = clusters[^1];
                double anchor = last[^1].Mass;
                double delta = anchor * ppmTolerance * 1e-6;
                if (p.Mass - anchor <= delta)
                    last.Add(p);
                else
                    clusters.Add(new List<(double, string)> { p });
            }

            var consensus = new List<double>();
            int refCount = referenceLabels.Count;
            foreach (var c in clusters)
            {
                if (c.Select(x => x.Algo).Distinct().Count() == refCount)
                    consensus.Add(c.Average(x => x.Mass));
            }
            return consensus;
        }

        /// <summary>
        /// For each consensus mass, returns the count of consensus entries that
        /// have at least one match in <paramref name="candidates"/> within
        /// <paramref name="ppmTolerance"/>. Uses two-pointer scan after sorting
        /// candidates so the lookup is O((n+m) log m).
        /// </summary>
        private static int CountMassesWithinPpm(
            IReadOnlyList<double> consensus,
            IReadOnlyList<double> candidates,
            double ppmTolerance)
        {
            if (consensus.Count == 0 || candidates.Count == 0) return 0;
            var sorted = candidates.OrderBy(m => m).ToList();
            int matched = 0;
            foreach (double target in consensus)
            {
                double delta = target * ppmTolerance * 1e-6;
                int lo = LowerBound(sorted, target - delta);
                if (lo < sorted.Count && sorted[lo] <= target + delta)
                    matched++;
            }
            return matched;
        }

        private static int LowerBound(List<double> sorted, double key)
        {
            int lo = 0, hi = sorted.Count;
            while (lo < hi)
            {
                int mid = lo + (hi - lo) / 2;
                if (sorted[mid] < key) lo = mid + 1;
                else hi = mid;
            }
            return lo;
        }
    }
}


