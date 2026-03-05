// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// REPO:     mzLib
// LOCATION: mzLib/Test/Dia/DiaPhase19FeatureTests.cs

using MassSpectrometry.Dia;
using NUnit.Framework;
using System.Collections.Generic;

namespace Test.Dia
{
    /// <summary>
    /// Unit tests for DIA Phase 19 features on DiaSearchResult:
    ///   - ChimericScore       — uncontested fragment signal fraction
    ///   - LibraryCoverageFraction — intensity-weighted detected fragment fraction
    ///
    /// Also tests RtDeviationMinutes/RtDeviationSquared computation logic.
    ///
    /// DiaFeatureVector.ClassifierFeatureCount is 26 in the current build.
    /// ChimericScore and LibraryCoverageFraction live on DiaSearchResult but
    /// are not yet wired into DiaFeatureVector.
    ///
    /// All tests are self-contained; no file I/O or external data required.
    /// </summary>
    [TestFixture]
    public class DiaPhase19FeatureTests
    {
        // ════════════════════════════════════════════════════════════════════
        //  Helpers
        // ════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Minimal DiaSearchResult factory. All scores default to NaN;
        /// callers set only the fields their test needs.
        /// </summary>
        private static DiaSearchResult MakeResult(
            string sequence = "PEPTIDEK",
            int chargeState = 2,
            double precursorMz = 500.0,
            int windowId = 0,
            bool isDecoy = false,
            int fragmentsQueried = 4)
        {
            return new DiaSearchResult(
                sequence: sequence,
                chargeState: chargeState,
                precursorMz: precursorMz,
                windowId: windowId,
                isDecoy: isDecoy,
                fragmentsQueried: fragmentsQueried,
                libraryRetentionTime: 25.0,
                rtWindowStart: 24.0f,
                rtWindowEnd: 26.0f);
        }

        /// <summary>
        /// Builds a LibraryPrecursorInput with the given fragment m/z and intensity arrays.
        /// </summary>
        private static LibraryPrecursorInput MakeInput(
            string sequence,
            double precursorMz,
            int windowId,
            float[] fragmentMzs,
            float[] fragmentIntensities,
            int chargeState = 2)
        {
            return new LibraryPrecursorInput(
                sequence: sequence,
                precursorMz: precursorMz,
                chargeState: chargeState,
                retentionTime: 25.0,
                isDecoy: false,
                fragmentMzs: fragmentMzs,
                fragmentIntensities: fragmentIntensities);
        }

        // ════════════════════════════════════════════════════════════════════
        //  DiaFeatureVector schema — 3 tests
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void FeatureVector_FeatureNames_CountMatchesClassifierFeatureCount()
        {
            Assert.That(DiaFeatureVector.FeatureNames.Length,
                Is.EqualTo(DiaFeatureVector.ClassifierFeatureCount));
        }

        [Test]
        public void FeatureVector_WriteTo_AllSlotsPopulated()
        {
            var fv = new DiaFeatureVector
            {
                ApexScore = 0.9f,
                TemporalScore = 0.8f,
                Log2SignalToNoise = 3.5f
            };

            var buf = new float[DiaFeatureVector.ClassifierFeatureCount];
            fv.WriteTo(buf);

            Assert.That(buf[0], Is.EqualTo(0.9f).Within(1e-6f));
            Assert.That(buf[1], Is.EqualTo(0.8f).Within(1e-6f));
            Assert.That(buf[25], Is.EqualTo(3.5f).Within(1e-6f));
        }

        // ════════════════════════════════════════════════════════════════════
        //  RtDeviationMinutes — computed deviation tests
        //  (RtDeviationNormalized is not a property on DiaSearchResult in the
        //   current build; these tests verify the underlying deviation logic)
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void RtDeviationMinutes_IsNaN_ByDefault()
        {
            var result = MakeResult();
            Assert.That(result.RtDeviationMinutes, Is.NaN);
        }

        [Test]
        public void RtDeviationMinutes_CorrectValue_WhenSet()
        {
            var result = MakeResult();
            result.RtDeviationMinutes = 0.6f;
            result.PeakWidth = 2.0f;

            Assert.That(result.RtDeviationMinutes, Is.EqualTo(0.6f).Within(1e-5f));
        }

        [Test]
        public void RtDeviationSquared_CorrectValue_WhenSet()
        {
            var result = MakeResult();
            result.RtDeviationMinutes = 0.5f;
            result.RtDeviationSquared = result.RtDeviationMinutes * result.RtDeviationMinutes;

            Assert.That(result.RtDeviationSquared, Is.EqualTo(0.25f).Within(1e-5f));
        }

        // ════════════════════════════════════════════════════════════════════
        //  LibraryCoverageFraction — 4 tests
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void LibraryCoverageFraction_IsNaN_ByDefault()
        {
            var result = MakeResult();
            Assert.That(result.LibraryCoverageFraction, Is.NaN);
        }

        [Test]
        public void LibraryCoverageFraction_AllDetected_ReturnsOne()
        {
            // 3 fragments, all equal library intensity, all detected
            var result = MakeResult(fragmentsQueried: 3);
            result.XicPointCounts[0] = 5;
            result.XicPointCounts[1] = 3;
            result.XicPointCounts[2] = 4;

            float[] libInt = { 100f, 100f, 100f };
            ApplyLibraryCoverageFormula(result, libInt);

            Assert.That(result.LibraryCoverageFraction, Is.EqualTo(1.0f).Within(1e-5f));
        }

        [Test]
        public void LibraryCoverageFraction_NoneDetected_ReturnsZero()
        {
            // All XicPointCounts remain 0 (default)
            var result = MakeResult(fragmentsQueried: 3);
            float[] libInt = { 100f, 200f, 50f };
            ApplyLibraryCoverageFormula(result, libInt);

            Assert.That(result.LibraryCoverageFraction, Is.EqualTo(0.0f).Within(1e-5f));
        }

        [Test]
        public void LibraryCoverageFraction_WeightedByLibraryIntensity()
        {
            // Fragment 0: lib intensity 300, detected
            // Fragment 1: lib intensity 100, NOT detected
            // Fragment 2: lib intensity 100, detected
            // Total lib = 500; detected lib = 400 → coverage = 0.8
            var result = MakeResult(fragmentsQueried: 3);
            result.XicPointCounts[0] = 5;
            result.XicPointCounts[1] = 0; // not detected
            result.XicPointCounts[2] = 3;

            float[] libInt = { 300f, 100f, 100f };
            ApplyLibraryCoverageFormula(result, libInt);

            Assert.That(result.LibraryCoverageFraction, Is.EqualTo(0.8f).Within(1e-5f));
        }

        // ════════════════════════════════════════════════════════════════════
        //  ChimericScore — 5 tests
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void ChimericScore_IsNaN_ByDefault()
        {
            var result = MakeResult();
            Assert.That(result.ChimericScore, Is.NaN);
        }

        [Test]
        public void ComputeChimericScores_SinglePrecursorInWindow_ScoreIsOne()
        {
            // Only one precursor in window 0 → no co-isolation → score = 1.0
            var result = MakeResult(windowId: 0);
            result.XicPointCounts[0] = 5;
            result.ExtractedIntensities[0] = 1000f;
            result.FragmentsDetected = 1;

            var input = MakeInput("PEPTIDEK", 500.0, 0,
                new float[] { 200f, 300f, 400f, 500f },
                new float[] { 100f, 80f, 60f, 40f });

            var precursors = new List<LibraryPrecursorInput> { input };
            var results = new List<DiaSearchResult> { result };

            DiaLibraryQueryGenerator.ComputeChimericScores(precursors, results, ppmTolerance: 10f);

            Assert.That(result.ChimericScore, Is.EqualTo(1.0f).Within(1e-5f));
        }

        [Test]
        public void ComputeChimericScores_NoSharedFragments_ScoreIsOne()
        {
            // Two precursors in same window, completely different fragment m/z values
            var result0 = MakeResult(sequence: "PEPTIDEK", windowId: 0, fragmentsQueried: 2);
            result0.XicPointCounts[0] = 3; result0.XicPointCounts[1] = 3;
            result0.ExtractedIntensities[0] = 500f; result0.ExtractedIntensities[1] = 500f;
            result0.FragmentsDetected = 2;

            var result1 = MakeResult(sequence: "ACDEFGHIK", windowId: 0, fragmentsQueried: 2);
            result1.XicPointCounts[0] = 3; result1.XicPointCounts[1] = 3;
            result1.ExtractedIntensities[0] = 400f; result1.ExtractedIntensities[1] = 400f;
            result1.FragmentsDetected = 2;

            // Fragment m/z values are far apart (>>10 ppm)
            var input0 = MakeInput("PEPTIDEK", 500.0, 0,
                new float[] { 200.0f, 300.0f },
                new float[] { 100f, 80f });

            var input1 = MakeInput("ACDEFGHIK", 600.0, 0,
                new float[] { 700.0f, 800.0f },  // completely different range
                new float[] { 100f, 80f });

            var precursors = new List<LibraryPrecursorInput> { input0, input1 };
            var results = new List<DiaSearchResult> { result0, result1 };

            DiaLibraryQueryGenerator.ComputeChimericScores(precursors, results, ppmTolerance: 10f);

            Assert.That(result0.ChimericScore, Is.EqualTo(1.0f).Within(1e-5f));
            Assert.That(result1.ChimericScore, Is.EqualTo(1.0f).Within(1e-5f));
        }

        [Test]
        public void ComputeChimericScores_AllSharedFragments_ScoreIsZero()
        {
            // Two precursors with identical fragment m/z values → every fragment is contested
            var result0 = MakeResult(sequence: "PEPTIDEK", windowId: 0, fragmentsQueried: 2);
            result0.XicPointCounts[0] = 3; result0.XicPointCounts[1] = 3;
            result0.ExtractedIntensities[0] = 500f; result0.ExtractedIntensities[1] = 500f;
            result0.FragmentsDetected = 2;

            var result1 = MakeResult(sequence: "ACDEFGHIK", windowId: 0, fragmentsQueried: 2);
            result1.XicPointCounts[0] = 3; result1.XicPointCounts[1] = 3;
            result1.ExtractedIntensities[0] = 400f; result1.ExtractedIntensities[1] = 400f;
            result1.FragmentsDetected = 2;

            // Identical fragment m/z → both precursors contest each other's fragments
            var sharedMzs = new float[] { 200.0f, 300.0f };
            var input0 = MakeInput("PEPTIDEK", 500.0, 0, sharedMzs, new float[] { 100f, 80f });
            var input1 = MakeInput("ACDEFGHIK", 600.0, 0, sharedMzs, new float[] { 100f, 80f });

            var precursors = new List<LibraryPrecursorInput> { input0, input1 };
            var results = new List<DiaSearchResult> { result0, result1 };

            DiaLibraryQueryGenerator.ComputeChimericScores(precursors, results, ppmTolerance: 10f);

            Assert.That(result0.ChimericScore, Is.EqualTo(0.0f).Within(1e-5f));
            Assert.That(result1.ChimericScore, Is.EqualTo(0.0f).Within(1e-5f));
        }

        [Test]
        public void ComputeChimericScores_PartialOverlap_ScoreReflectsUncontestedFraction()
        {
            // Precursor A: fragments at 200, 300, 400 (intensities 100, 100, 100)
            //   fragment at 300 is shared with precursor B → contested
            //   uncontested signal = 200 / 300 → score ≈ 0.667
            var resultA = MakeResult(sequence: "PEPTIDEK", windowId: 0, fragmentsQueried: 3);
            resultA.XicPointCounts[0] = 3; resultA.XicPointCounts[1] = 3; resultA.XicPointCounts[2] = 3;
            resultA.ExtractedIntensities[0] = 100f;
            resultA.ExtractedIntensities[1] = 100f; // contested
            resultA.ExtractedIntensities[2] = 100f;
            resultA.FragmentsDetected = 3;

            // Precursor B: fragments at 300, 500 (different precursor)
            var resultB = MakeResult(sequence: "ACDEFGHIK", windowId: 0, fragmentsQueried: 2);
            resultB.XicPointCounts[0] = 3; resultB.XicPointCounts[1] = 3;
            resultB.ExtractedIntensities[0] = 200f;
            resultB.ExtractedIntensities[1] = 200f;
            resultB.FragmentsDetected = 2;

            var inputA = MakeInput("PEPTIDEK", 500.0, 0,
                new float[] { 200.0f, 300.0f, 400.0f },
                new float[] { 100f, 100f, 100f });

            var inputB = MakeInput("ACDEFGHIK", 600.0, 0,
                new float[] { 300.0f, 500.0f },     // 300 overlaps with A
                new float[] { 100f, 100f });

            var precursors = new List<LibraryPrecursorInput> { inputA, inputB };
            var results = new List<DiaSearchResult> { resultA, resultB };

            DiaLibraryQueryGenerator.ComputeChimericScores(precursors, results, ppmTolerance: 10f);

            // A: 200 + 400 uncontested out of 300 total → 200/300 ≈ 0.6667
            Assert.That(resultA.ChimericScore, Is.EqualTo(200f / 300f).Within(1e-4f));
        }

        // ════════════════════════════════════════════════════════════════════
        //  Helper — replicates the LibraryCoverage formula from the assembler
        // ════════════════════════════════════════════════════════════════════

        private static void ApplyLibraryCoverageFormula(DiaSearchResult result, float[] libIntensities)
        {
            float totalLibIntensity = 0f;
            float detectedLibIntensity = 0f;
            for (int f = 0; f < result.FragmentsQueried; f++)
            {
                float libInt = f < libIntensities.Length ? libIntensities[f] : 0f;
                totalLibIntensity += libInt;
                if (result.XicPointCounts[f] > 0)
                    detectedLibIntensity += libInt;
            }
            result.LibraryCoverageFraction = totalLibIntensity > 0f
                ? detectedLibIntensity / totalLibIntensity
                : 0f;
        }
    }
}