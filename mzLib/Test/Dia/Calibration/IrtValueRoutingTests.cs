// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using NUnit.Framework;

namespace Test.DiaTests
{
    /// <summary>
    /// Tests that IrtValue is correctly used wherever LibraryPrecursorInput is consumed
    /// in production mzLib code.
    ///
    /// Three production sites are covered:
    ///
    ///   1. DiaLibraryQueryGenerator.Generate()
    ///      Bug: only checked p.RetentionTime — if IrtValue was set and RetentionTime was null,
    ///      every precursor fell back to the full-run RT range.
    ///      Fix: resolve IrtValue ?? RetentionTime.
    ///
    ///   2. DiaLibraryQueryGenerator.AssembleResults()
    ///      Bug: passed input.RetentionTime as libraryRetentionTime on DiaSearchResult — null
    ///      for iRT libraries, so RecalibrateRtDeviations skipped every result.
    ///      Fix: pass input.IrtValue ?? input.RetentionTime.
    ///
    ///   3. DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring()
    ///      Same bug and same fix as (2).
    ///
    /// Placement: mzLib/Test/Dia/IrtValueRoutingTests.cs
    /// </summary>
    [TestFixture]
    public class IrtValueRoutingTests
    {
        // ─────────────────────────────────────────────────────────────────────
        //  Helpers
        // ─────────────────────────────────────────────────────────────────────

        private static LibraryPrecursorInput MakePrecursor(
            double precursorMz,
            double? retentionTime = null,
            double? irtValue = null,
            bool isDecoy = false)
        {
            return new LibraryPrecursorInput(
                sequence: "PEPTIDEK",
                precursorMz: precursorMz,
                chargeState: 2,
                retentionTime: retentionTime,
                isDecoy: isDecoy,
                fragmentMzs: new float[] { 100f, 200f, 300f },
                fragmentIntensities: new float[] { 0.9f, 0.5f, 0.1f },
                irtValue: irtValue);
        }

        /// <summary>
        /// Builds a minimal DiaScanIndex with one DIA window covering [windowMzLow, windowMzHigh].
        /// One MS2 scan is included mid-run so GetGlobalRtMin/Max return rtMin/rtMax.
        /// Uses DiaScanIndexBuilder.Build(), which is the only public construction path.
        /// </summary>
        private static DiaScanIndex MakeMinimalScanIndex(
            double rtMin, double rtMax, double windowMzLow, double windowMzHigh)
        {
            double windowCenter = (windowMzLow + windowMzHigh) / 2.0;
            double windowWidth = windowMzHigh - windowMzLow;
            double midRt = (rtMin + rtMax) / 2.0;

            // Two scans at rtMin and rtMax so GetGlobalRtMin/Max return the expected values
            var scans = new[]
            {
                new MsDataScan(
                    massSpectrum:    new MzSpectrum(new double[] { 100.0 }, new double[] { 1000.0 }, false),
                    oneBasedScanNumber: 1,
                    msnOrder:        2,
                    isCentroid:      true,
                    polarity:        Polarity.Positive,
                    retentionTime:   rtMin,
                    scanWindowRange: new MzRange(100, 2000),
                    scanFilter:      "FTMS",
                    mzAnalyzer:      MZAnalyzerType.Orbitrap,
                    totalIonCurrent: 1000.0,
                    injectionTime:   20.0,
                    noiseData:       null,
                    nativeId:        "scan=1",
                    isolationMZ:     windowCenter,
                    isolationWidth:  windowWidth,
                    dissociationType: DissociationType.HCD),

                new MsDataScan(
                    massSpectrum:    new MzSpectrum(new double[] { 100.0 }, new double[] { 1000.0 }, false),
                    oneBasedScanNumber: 2,
                    msnOrder:        2,
                    isCentroid:      true,
                    polarity:        Polarity.Positive,
                    retentionTime:   rtMax,
                    scanWindowRange: new MzRange(100, 2000),
                    scanFilter:      "FTMS",
                    mzAnalyzer:      MZAnalyzerType.Orbitrap,
                    totalIonCurrent: 1000.0,
                    injectionTime:   20.0,
                    noiseData:       null,
                    nativeId:        "scan=2",
                    isolationMZ:     windowCenter,
                    isolationWidth:  windowWidth,
                    dissociationType: DissociationType.HCD),
            };

            return DiaScanIndexBuilder.Build(scans);
        }

        /// <summary>
        /// Calls AssembleResults with a single precursor and one data point per fragment
        /// so MeetsMinFragments(1) passes. Only LibraryRetentionTime is asserted.
        /// FragmentResult signature: (queryId, dataPointCount, rtBufferOffset,
        ///                            intensityBufferOffset, totalIntensity)
        /// </summary>
        private static DiaSearchResult BuildSingleResult(LibraryPrecursorInput precursor)
        {
            var precursors = new List<LibraryPrecursorInput> { precursor };

            var group = new DiaLibraryQueryGenerator.PrecursorQueryGroup(
                inputIndex: 0,
                queryOffset: 0,
                queryCount: precursor.FragmentCount,
                windowId: 0,
                rtMin: 40f,
                rtMax: 50f);

            var genResult = new DiaLibraryQueryGenerator.GenerationResult(
                queries: new FragmentQuery[precursor.FragmentCount],
                precursorGroups: new[] { group },
                skippedNoWindow: 0,
                skippedNoFragments: 0);

            var extractionResults = new FragmentResult[precursor.FragmentCount];
            for (int i = 0; i < extractionResults.Length; i++)
                extractionResults[i] = new FragmentResult(
                    queryId: i,
                    dataPointCount: 1,
                    rtBufferOffset: i,
                    intensityBufferOffset: i,
                    totalIntensity: 100f);

            var parameters = new DiaSearchParameters
            {
                MinFragmentsRequired = 1,
                MinScoreThreshold = 0f
            };

            var results = DiaLibraryQueryGenerator.AssembleResults(
                precursors, genResult, extractionResults, parameters);

            Assert.That(results.Count, Is.EqualTo(1),
                "Expected exactly one result from AssembleResults");
            return results[0];
        }

        /// <summary>
        /// Calls AssembleResultsWithTemporalScoring with zero data points per fragment.
        /// The "maxPts == 0" early-out branch adds the result before any scoring, so
        /// LibraryRetentionTime is set correctly regardless.
        /// MinFragmentsRequired = 0 bypasses the fragment-count pre-filter.
        /// ExtractionResult signature: (FragmentResult[] results, float[] rtBuffer,
        ///                              float[] intensityBuffer, int totalDataPoints)
        /// </summary>
        private static DiaSearchResult BuildSingleResultTemporal(LibraryPrecursorInput precursor)
        {
            var precursors = new List<LibraryPrecursorInput> { precursor };

            var group = new DiaLibraryQueryGenerator.PrecursorQueryGroup(
                inputIndex: 0,
                queryOffset: 0,
                queryCount: precursor.FragmentCount,
                windowId: 0,
                rtMin: 40f,
                rtMax: 50f);

            var genResult = new DiaLibraryQueryGenerator.GenerationResult(
                queries: new FragmentQuery[precursor.FragmentCount],
                precursorGroups: new[] { group },
                skippedNoWindow: 0,
                skippedNoFragments: 0);

            int fragCount = precursor.FragmentCount;
            var fragResults = new FragmentResult[fragCount];
            for (int i = 0; i < fragCount; i++)
                fragResults[i] = new FragmentResult(
                    queryId: i,
                    dataPointCount: 0,
                    rtBufferOffset: 0,
                    intensityBufferOffset: 0,
                    totalIntensity: 0f);

            var extractionResult = new ExtractionResult(
                results: fragResults,
                rtBuffer: System.Array.Empty<float>(),
                intensityBuffer: System.Array.Empty<float>(),
                totalDataPoints: 0);

            var parameters = new DiaSearchParameters
            {
                MinFragmentsRequired = 0,   // bypass fragment-count filter
                MinScoreThreshold = 0f
            };

            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursors, genResult, extractionResult, parameters);

            Assert.That(results.Count, Is.EqualTo(1),
                "Expected one result (MinFragmentsRequired=0 bypasses fragment count filter)");
            return results[0];
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Group 1: LibraryPrecursorInput field population
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void LibraryPrecursorInput_NativeRt_RetentionTimeSetIrtNull()
        {
            var p = MakePrecursor(550.3, retentionTime: 45.7, irtValue: null);
            Assert.That(p.RetentionTime, Is.EqualTo(45.7).Within(0.001));
            Assert.That(p.IrtValue, Is.Null);
        }

        [Test]
        public void LibraryPrecursorInput_IrtLibrary_IrtValueSetRetentionTimeNull()
        {
            var p = MakePrecursor(550.3, retentionTime: null, irtValue: -12.5);
            Assert.That(p.IrtValue, Is.EqualTo(-12.5).Within(0.001));
            Assert.That(p.RetentionTime, Is.Null);
        }

        [Test]
        public void LibraryPrecursorInput_BothNull_BothNull()
        {
            var p = MakePrecursor(550.3);
            Assert.That(p.RetentionTime, Is.Null);
            Assert.That(p.IrtValue, Is.Null);
        }

        [Test]
        public void LibraryPrecursorInput_BothSet_BothPreserved()
        {
            var p = MakePrecursor(550.3, retentionTime: 45.7, irtValue: -12.5);
            Assert.That(p.RetentionTime, Is.EqualTo(45.7).Within(0.001));
            Assert.That(p.IrtValue, Is.EqualTo(-12.5).Within(0.001));
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Group 2: AssembleResults — LibraryRetentionTime on DiaSearchResult
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void AssembleResults_NativeRt_LibraryRetentionTimeIsRetentionTime()
        {
            var result = BuildSingleResult(
                MakePrecursor(550.3, retentionTime: 45.7, irtValue: null));

            Assert.That(result.LibraryRetentionTime, Is.EqualTo(45.7).Within(0.001),
                "Native-RT library: LibraryRetentionTime must come from RetentionTime");
        }

        [Test]
        public void AssembleResults_IrtLibrary_LibraryRetentionTimeIsIrtValue()
        {
            // Before fix: RetentionTime was null → LibraryRetentionTime = null
            //             → RecalibrateRtDeviations silently skipped every result
            // After fix:  IrtValue ?? RetentionTime = -12.5 → correctly propagated
            var result = BuildSingleResult(
                MakePrecursor(550.3, retentionTime: null, irtValue: -12.5));

            Assert.That(result.LibraryRetentionTime, Is.EqualTo(-12.5).Within(0.001),
                "iRT library: LibraryRetentionTime must come from IrtValue, not be null");
        }

        [Test]
        public void AssembleResults_BothNull_LibraryRetentionTimeIsNull()
        {
            var result = BuildSingleResult(
                MakePrecursor(550.3));

            Assert.That(result.LibraryRetentionTime, Is.Null);
        }

        [Test]
        public void AssembleResults_IrtValuePreferredOverRetentionTime()
        {
            // IrtValue wins when both are set — it is what calibration was fitted on
            var result = BuildSingleResult(
                MakePrecursor(550.3, retentionTime: 45.7, irtValue: -12.5));

            Assert.That(result.LibraryRetentionTime, Is.EqualTo(-12.5).Within(0.001),
                "When both fields are set, IrtValue takes priority");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Group 3: AssembleResultsWithTemporalScoring — same fix
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void AssembleResultsWithTemporalScoring_IrtLibrary_LibraryRetentionTimeIsIrtValue()
        {
            var result = BuildSingleResultTemporal(
                MakePrecursor(550.3, retentionTime: null, irtValue: 88.4));

            Assert.That(result.LibraryRetentionTime, Is.EqualTo(88.4).Within(0.001),
                "Temporal assembly: LibraryRetentionTime must come from IrtValue for iRT libraries");
        }

        [Test]
        public void AssembleResultsWithTemporalScoring_NativeRt_LibraryRetentionTimeIsRetentionTime()
        {
            var result = BuildSingleResultTemporal(
                MakePrecursor(550.3, retentionTime: 45.7, irtValue: null));

            Assert.That(result.LibraryRetentionTime, Is.EqualTo(45.7).Within(0.001));
        }

        [Test]
        public void AssembleResultsWithTemporalScoring_BothNull_LibraryRetentionTimeIsNull()
        {
            var result = BuildSingleResultTemporal(
                MakePrecursor(550.3));

            Assert.That(result.LibraryRetentionTime, Is.Null);
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Group 4: Generate() — RT window uses IrtValue ?? RetentionTime
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void Generate_NativeRt_WindowCenteredOnRetentionTime()
        {
            using var scanIndex = MakeMinimalScanIndex(0.0, 120.0, 540.0, 560.0);
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor(550.3, retentionTime: 45.7, irtValue: null)
            };
            var parameters = new DiaSearchParameters
            {
                RtToleranceMinutes = 5f,
                PpmTolerance = 10f,
                MinFragmentsRequired = 1
            };

            var result = DiaLibraryQueryGenerator.Generate(precursors, scanIndex, parameters);

            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(1));
            var group = result.PrecursorGroups[0];
            Assert.That(group.RtMin, Is.EqualTo(40.7f).Within(0.01f),
                "RT window min must be RetentionTime - tolerance");
            Assert.That(group.RtMax, Is.EqualTo(50.7f).Within(0.01f),
                "RT window max must be RetentionTime + tolerance");
        }

        [Test]
        public void Generate_IrtLibrary_WindowCenteredOnIrtValue_NotFullRun()
        {
            // Before fix: RetentionTime = null → window = [0, 120] (full run, width 120)
            // After fix:  IrtValue = -12.5  → window = [-17.5, -7.5] (width 10)
            using var scanIndex = MakeMinimalScanIndex(0.0, 120.0, 540.0, 560.0);
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor(550.3, retentionTime: null, irtValue: -12.5)
            };
            var parameters = new DiaSearchParameters
            {
                RtToleranceMinutes = 5f,
                PpmTolerance = 10f,
                MinFragmentsRequired = 1
            };

            var result = DiaLibraryQueryGenerator.Generate(precursors, scanIndex, parameters);

            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(1));
            var group = result.PrecursorGroups[0];

            // Regression guard: window must not span the full run
            Assert.That(group.RtMax - group.RtMin, Is.LessThan(100f),
                "Bug regression: window must not span the full run when IrtValue is set");

            // Correct behaviour: centred on IrtValue ± tolerance
            Assert.That(group.RtMin, Is.EqualTo(-17.5f).Within(0.01f),
                "Window min must be IrtValue - tolerance");
            Assert.That(group.RtMax, Is.EqualTo(-7.5f).Within(0.01f),
                "Window max must be IrtValue + tolerance");
        }

        [Test]
        public void Generate_BothNull_WindowIsFullRun()
        {
            using var scanIndex = MakeMinimalScanIndex(0.0, 120.0, 540.0, 560.0);
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor(550.3)
            };
            var parameters = new DiaSearchParameters
            {
                RtToleranceMinutes = 5f,
                PpmTolerance = 10f,
                MinFragmentsRequired = 1
            };

            var result = DiaLibraryQueryGenerator.Generate(precursors, scanIndex, parameters);

            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(1));
            var group = result.PrecursorGroups[0];
            Assert.That(group.RtMin, Is.EqualTo(0f).Within(0.001f),
                "No RT: window must start at global RT min");
            Assert.That(group.RtMax, Is.EqualTo(120f).Within(0.001f),
                "No RT: window must end at global RT max");
        }

        [Test]
        public void Generate_IrtValuePreferredOverRetentionTime()
        {
            // When both are set, IrtValue wins (consistent with AssembleResults)
            using var scanIndex = MakeMinimalScanIndex(0.0, 120.0, 540.0, 560.0);
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor(550.3, retentionTime: 45.7, irtValue: -12.5)
            };
            var parameters = new DiaSearchParameters
            {
                RtToleranceMinutes = 5f,
                PpmTolerance = 10f,
                MinFragmentsRequired = 1
            };

            var result = DiaLibraryQueryGenerator.Generate(precursors, scanIndex, parameters);

            var group = result.PrecursorGroups[0];
            Assert.That(group.RtMin, Is.EqualTo(-17.5f).Within(0.01f),
                "IrtValue must take priority over RetentionTime in Generate()");
            Assert.That(group.RtMax, Is.EqualTo(-7.5f).Within(0.01f));
        }
    }
}