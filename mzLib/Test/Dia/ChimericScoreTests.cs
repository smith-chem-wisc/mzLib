// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System.Collections.Generic;
using MassSpectrometry.Dia;
using NUnit.Framework;

namespace Test.DiaTests
{
    /// <summary>
    /// Tests for the PrecursorIndex assembly bug and its downstream effect on ChimericScore.
    ///
    /// Root cause (as of this writing):
    ///   AssembleResultsWithTemporalScoring() creates each DiaSearchResult but never sets
    ///   result.PrecursorIndex = group.InputIndex. Every result enters ComputeChimericScores
    ///   with PrecursorIndex == -1, triggering the early-out branch:
    ///
    ///       if (piA < 0) { ra.ChimericScore = 1.0f; continue; }
    ///
    ///   … so ChimericScore is always 1.0f (dead feature, zero separation).
    ///
    /// Fix (one line in AssembleResultsWithTemporalScoring, after new DiaSearchResult(...)):
    ///   result.PrecursorIndex = group.InputIndex;
    ///
    /// Test structure:
    ///   Layer 1 — Assembly sets PrecursorIndex correctly           [RED until fix]
    ///   Layer 2 — ComputeChimericScores logic is correct when      [GREEN already]
    ///             PrecursorIndex is set manually
    ///   Layer 3 — End-to-end: assembly bug causes all scores=1.0f  [documents current behavior;
    ///             (regression guard — must go GREEN after fix)       flips GREEN after fix]
    ///
    /// Placement: mzLib/Test/Dia/ChimericScoreTests.cs
    /// </summary>
    [TestFixture]
    public class ChimericScoreTests
    {
        // ─────────────────────────────────────────────────────────────────────
        //  Shared helpers
        // ─────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Builds a minimal LibraryPrecursorInput at the given precursor m/z with
        /// the specified fragment m/z values.
        /// </summary>
        private static LibraryPrecursorInput MakePrecursor(
            double precursorMz,
            float[] fragmentMzs,
            bool isDecoy = false,
            string sequence = "PEPTIDEK")
        {
            var intensities = new float[fragmentMzs.Length];
            for (int i = 0; i < intensities.Length; i++)
                intensities[i] = 1000f - i * 100f; // descending, all positive

            return new LibraryPrecursorInput(
                sequence: sequence,
                precursorMz: precursorMz,
                chargeState: 2,
                retentionTime: 30.0,
                isDecoy: isDecoy,
                fragmentMzs: fragmentMzs,
                fragmentIntensities: intensities,
                irtValue: null);
        }

        /// <summary>
        /// Builds a minimal GenerationResult for a single precursor so we can call
        /// AssembleResultsWithTemporalScoring without a real extraction pipeline.
        /// All fragments get zero data points → the maxPts==0 early-out branch fires,
        /// which still creates and adds the result (so PrecursorIndex is observable).
        /// MinFragmentsRequired must be 0 to bypass the fragment-count pre-filter.
        /// </summary>
        private static (DiaLibraryQueryGenerator.GenerationResult genResult,
                        ExtractionResult extractionResult,
                        DiaSearchParameters parameters)
            BuildMinimalInputs(LibraryPrecursorInput precursor, int inputIndex = 0)
        {
            var group = new DiaLibraryQueryGenerator.PrecursorQueryGroup(
                inputIndex: inputIndex,
                queryOffset: 0,
                queryCount: precursor.FragmentCount,
                windowId: 7,
                rtMin: 28f,
                rtMax: 32f);

            var genResult = new DiaLibraryQueryGenerator.GenerationResult(
                queries: new FragmentQuery[precursor.FragmentCount],
                precursorGroups: new[] { group },
                skippedNoWindow: 0,
                skippedNoFragments: 0);

            // Zero data points per fragment → maxPts==0 → early-out adds result immediately
            var fragResults = new FragmentResult[precursor.FragmentCount];
            for (int i = 0; i < fragResults.Length; i++)
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
                MinFragmentsRequired = 0,
                MinScoreThreshold = 0f,
                PpmTolerance = 10f
            };

            return (genResult, extractionResult, parameters);
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 1: Assembly sets PrecursorIndex
        //  These tests are RED until the one-line fix is applied.
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        [Description("RED until fix: AssembleResultsWithTemporalScoring must set PrecursorIndex = group.InputIndex")]
        public void AssembleResultsWithTemporalScoring_SetsPrecursorIndex_ForFirstGroup()
        {
            var precursor = MakePrecursor(550.3, new float[] { 100f, 200f, 300f });
            var precursors = new List<LibraryPrecursorInput> { precursor };
            var (genResult, extractionResult, parameters) = BuildMinimalInputs(precursor, inputIndex: 0);

            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursors, genResult, extractionResult, parameters);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].PrecursorIndex, Is.EqualTo(0),
                "PrecursorIndex must equal group.InputIndex (0). " +
                "Fails today because AssembleResultsWithTemporalScoring never sets it.");
        }

        [Test]
        [Description("RED until fix: PrecursorIndex must reflect the actual InputIndex, not default -1")]
        public void AssembleResultsWithTemporalScoring_SetsPrecursorIndex_NonZeroInputIndex()
        {
            // Simulate a precursor that lives at position 5 in the combined precursor list
            // (e.g. 5th entry in targets+decoys). InputIndex must survive into the result.
            var precursor = MakePrecursor(550.3, new float[] { 100f, 200f, 300f });
            var precursors = new List<LibraryPrecursorInput>
            {
                // Pad with dummy entries so index 5 is valid
                MakePrecursor(400.0, new float[] { 99f }, sequence: "AAAK"),
                MakePrecursor(401.0, new float[] { 99f }, sequence: "BBAK"),
                MakePrecursor(402.0, new float[] { 99f }, sequence: "CCAK"),
                MakePrecursor(403.0, new float[] { 99f }, sequence: "DDAK"),
                MakePrecursor(404.0, new float[] { 99f }, sequence: "EEAK"),
                precursor,   // index 5
            };
            var (genResult, extractionResult, parameters) = BuildMinimalInputs(precursor, inputIndex: 5);

            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursors, genResult, extractionResult, parameters);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].PrecursorIndex, Is.EqualTo(5),
                "PrecursorIndex must equal group.InputIndex (5). " +
                "Fails today because AssembleResultsWithTemporalScoring never sets it.");
        }

        [Test]
        [Description("RED until fix: default PrecursorIndex=-1 causes ChimericScore=1.0 even with real overlap")]
        public void AssembleResultsWithTemporalScoring_PrecursorIndex_NotMinusOne()
        {
            // This is the simplest possible assertion: just that the default -1 was overwritten.
            var precursor = MakePrecursor(550.3, new float[] { 100f, 200f });
            var precursors = new List<LibraryPrecursorInput> { precursor };
            var (genResult, extractionResult, parameters) = BuildMinimalInputs(precursor, inputIndex: 0);

            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursors, genResult, extractionResult, parameters);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].PrecursorIndex, Is.Not.EqualTo(-1),
                "PrecursorIndex must not be -1 after assembly. " +
                "Currently always -1 because the assignment is missing.");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 2: ComputeChimericScores logic is correct when PrecursorIndex
        //  is set manually. These tests are GREEN already — they pin the
        //  chimeric scoring logic independent of the assembly bug.
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        [Description("GREEN: When no fragments overlap, ChimericScore must be 1.0 (all uncontested)")]
        public void ComputeChimericScores_NoOverlap_ScoreIsOne()
        {
            // Two precursors in the same window with completely non-overlapping fragments.
            // 400 ppm apart at ~200 Da → 0.08 Da separation, far exceeds 10 ppm tolerance.
            var precA = MakePrecursor(550.3, new float[] { 100.0f, 200.0f, 300.0f });
            var precB = MakePrecursor(560.3, new float[] { 150.0f, 250.0f, 350.0f });
            var precursors = new List<LibraryPrecursorInput> { precA, precB };

            // Build results manually with PrecursorIndex already set correctly
            var resultA = new DiaSearchResult(
                sequence: precA.Sequence, chargeState: 2, precursorMz: precA.PrecursorMz,
                windowId: 0, isDecoy: false, fragmentsQueried: 3,
                libraryRetentionTime: 30.0, rtWindowStart: 28f, rtWindowEnd: 32f);
            resultA.PrecursorIndex = 0;
            resultA.ExtractedIntensities[0] = 500f;
            resultA.ExtractedIntensities[1] = 300f;
            resultA.ExtractedIntensities[2] = 100f;

            var resultB = new DiaSearchResult(
                sequence: precB.Sequence, chargeState: 2, precursorMz: precB.PrecursorMz,
                windowId: 0, isDecoy: false, fragmentsQueried: 3,
                libraryRetentionTime: 30.0, rtWindowStart: 28f, rtWindowEnd: 32f);
            resultB.PrecursorIndex = 1;
            resultB.ExtractedIntensities[0] = 400f;
            resultB.ExtractedIntensities[1] = 200f;
            resultB.ExtractedIntensities[2] = 100f;

            var results = new List<DiaSearchResult> { resultA, resultB };

            DiaLibraryQueryGenerator.ComputeChimericScores(precursors, results, ppmTolerance: 10f);

            Assert.That(resultA.ChimericScore, Is.EqualTo(1.0f).Within(1e-4f),
                "No fragment overlap → all signal uncontested → score must be 1.0");
            Assert.That(resultB.ChimericScore, Is.EqualTo(1.0f).Within(1e-4f),
                "No fragment overlap → all signal uncontested → score must be 1.0");
        }

        [Test]
        [Description("GREEN: Fully overlapping fragments → ChimericScore must be 0.0")]
        public void ComputeChimericScores_FullOverlap_ScoreIsZero()
        {
            // Both precursors have identical fragment m/z values → every fragment is contested.
            // With co-eluting RT windows, all signal is contested → score = 0.
            float[] sharedFragments = { 100.0f, 200.0f, 300.0f };

            var precA = MakePrecursor(550.3, sharedFragments, sequence: "PEPTIDEK");
            var precB = MakePrecursor(560.3, sharedFragments, sequence: "KEPTIDEM");
            var precursors = new List<LibraryPrecursorInput> { precA, precB };

            var resultA = new DiaSearchResult(
                sequence: precA.Sequence, chargeState: 2, precursorMz: precA.PrecursorMz,
                windowId: 0, isDecoy: false, fragmentsQueried: 3,
                libraryRetentionTime: 30.0, rtWindowStart: 28f, rtWindowEnd: 32f);
            resultA.PrecursorIndex = 0;
            resultA.ExtractedIntensities[0] = 500f;
            resultA.ExtractedIntensities[1] = 300f;
            resultA.ExtractedIntensities[2] = 100f;

            var resultB = new DiaSearchResult(
                sequence: precB.Sequence, chargeState: 2, precursorMz: precB.PrecursorMz,
                windowId: 0, isDecoy: false, fragmentsQueried: 3,
                libraryRetentionTime: 30.0, rtWindowStart: 28f, rtWindowEnd: 32f);
            resultB.PrecursorIndex = 1;
            resultB.ExtractedIntensities[0] = 400f;
            resultB.ExtractedIntensities[1] = 200f;
            resultB.ExtractedIntensities[2] = 100f;

            var results = new List<DiaSearchResult> { resultA, resultB };

            DiaLibraryQueryGenerator.ComputeChimericScores(precursors, results, ppmTolerance: 10f);

            Assert.That(resultA.ChimericScore, Is.EqualTo(0.0f).Within(1e-4f),
                "All fragments contested → score must be 0.0");
            Assert.That(resultB.ChimericScore, Is.EqualTo(0.0f).Within(1e-4f),
                "All fragments contested → score must be 0.0");
        }

        [Test]
        [Description("GREEN: Partial overlap → ChimericScore is the intensity-weighted uncontested fraction")]
        public void ComputeChimericScores_PartialOverlap_ScoreIsUncontestedFraction()
        {
            // precA has 3 fragments: 100, 200, 300.
            // precB has fragment at 200 only (overlaps precA's middle fragment).
            // precA intensities: 100f, 300f, 200f → total = 600f
            //   fragment 100: uncontested → 100f
            //   fragment 200: contested   →   0f
            //   fragment 300: uncontested → 200f
            // expected score for precA = 300f / 600f = 0.5f

            var precA = MakePrecursor(550.3, new float[] { 100.0f, 200.0f, 300.0f });
            var precB = MakePrecursor(560.3, new float[] { 200.0f });
            var precursors = new List<LibraryPrecursorInput> { precA, precB };

            var resultA = new DiaSearchResult(
                sequence: precA.Sequence, chargeState: 2, precursorMz: precA.PrecursorMz,
                windowId: 0, isDecoy: false, fragmentsQueried: 3,
                libraryRetentionTime: 30.0, rtWindowStart: 28f, rtWindowEnd: 32f);
            resultA.PrecursorIndex = 0;
            resultA.ExtractedIntensities[0] = 100f;
            resultA.ExtractedIntensities[1] = 300f;
            resultA.ExtractedIntensities[2] = 200f;

            var resultB = new DiaSearchResult(
                sequence: precB.Sequence, chargeState: 2, precursorMz: precB.PrecursorMz,
                windowId: 0, isDecoy: false, fragmentsQueried: 1,
                libraryRetentionTime: 30.0, rtWindowStart: 28f, rtWindowEnd: 32f);
            resultB.PrecursorIndex = 1;
            resultB.ExtractedIntensities[0] = 400f;

            var results = new List<DiaSearchResult> { resultA, resultB };

            DiaLibraryQueryGenerator.ComputeChimericScores(precursors, results, ppmTolerance: 10f);

            // precA: uncontested = 100f + 200f = 300f, total = 600f → 0.5f
            Assert.That(resultA.ChimericScore, Is.EqualTo(0.5f).Within(1e-4f),
                "precA: fragments 100 and 300 are uncontested (300/600 = 0.5)");

            // precB: its only fragment (200) is contested by precA → score = 0.0f
            Assert.That(resultB.ChimericScore, Is.EqualTo(0.0f).Within(1e-4f),
                "precB: sole fragment is contested by precA → score = 0.0");
        }

        [Test]
        [Description("GREEN: Results in different windows never contest each other")]
        public void ComputeChimericScores_DifferentWindows_NeverContest()
        {
            // Same fragment m/z but different WindowId → they are never grouped together
            // → each result sees no co-eluters → both scores = 1.0
            float[] sharedFragments = { 100.0f, 200.0f };

            var precA = MakePrecursor(550.3, sharedFragments, sequence: "PEPTIDEK");
            var precB = MakePrecursor(560.3, sharedFragments, sequence: "KEPTIDEM");
            var precursors = new List<LibraryPrecursorInput> { precA, precB };

            var resultA = new DiaSearchResult(
                sequence: precA.Sequence, chargeState: 2, precursorMz: precA.PrecursorMz,
                windowId: 0, isDecoy: false, fragmentsQueried: 2,      // window 0
                libraryRetentionTime: 30.0, rtWindowStart: 28f, rtWindowEnd: 32f);
            resultA.PrecursorIndex = 0;
            resultA.ExtractedIntensities[0] = 500f;
            resultA.ExtractedIntensities[1] = 300f;

            var resultB = new DiaSearchResult(
                sequence: precB.Sequence, chargeState: 2, precursorMz: precB.PrecursorMz,
                windowId: 1, isDecoy: false, fragmentsQueried: 2,      // window 1
                libraryRetentionTime: 30.0, rtWindowStart: 28f, rtWindowEnd: 32f);
            resultB.PrecursorIndex = 1;
            resultB.ExtractedIntensities[0] = 400f;
            resultB.ExtractedIntensities[1] = 200f;

            var results = new List<DiaSearchResult> { resultA, resultB };

            DiaLibraryQueryGenerator.ComputeChimericScores(precursors, results, ppmTolerance: 10f);

            Assert.That(resultA.ChimericScore, Is.EqualTo(1.0f).Within(1e-4f),
                "Different windows → no co-eluters → score must be 1.0");
            Assert.That(resultB.ChimericScore, Is.EqualTo(1.0f).Within(1e-4f),
                "Different windows → no co-eluters → score must be 1.0");
        }

        [Test]
        [Description("GREEN: Non-overlapping RT windows prevent contestation even with matching fragment m/z")]
        public void ComputeChimericScores_NonOverlappingRtWindows_NeverContest()
        {
            // Same window ID, same fragment m/z, but RT windows don't overlap.
            // precA elutes at 10–14 min; precB at 20–24 min.
            // The RT guard:  if (rtEnd[b] <= rtStart[a] || rtStart[b] >= rtEnd[a]) continue;
            // ensures they don't contest each other.
            float[] sharedFragments = { 100.0f, 200.0f };

            var precA = MakePrecursor(550.3, sharedFragments, sequence: "PEPTIDEK");
            var precB = MakePrecursor(560.3, sharedFragments, sequence: "KEPTIDEM");
            var precursors = new List<LibraryPrecursorInput> { precA, precB };

            var resultA = new DiaSearchResult(
                sequence: precA.Sequence, chargeState: 2, precursorMz: precA.PrecursorMz,
                windowId: 0, isDecoy: false, fragmentsQueried: 2,
                libraryRetentionTime: 12.0, rtWindowStart: 10f, rtWindowEnd: 14f);
            resultA.PrecursorIndex = 0;
            resultA.ExtractedIntensities[0] = 500f;
            resultA.ExtractedIntensities[1] = 300f;

            var resultB = new DiaSearchResult(
                sequence: precB.Sequence, chargeState: 2, precursorMz: precB.PrecursorMz,
                windowId: 0, isDecoy: false, fragmentsQueried: 2,
                libraryRetentionTime: 22.0, rtWindowStart: 20f, rtWindowEnd: 24f);
            resultB.PrecursorIndex = 1;
            resultB.ExtractedIntensities[0] = 400f;
            resultB.ExtractedIntensities[1] = 200f;

            var results = new List<DiaSearchResult> { resultA, resultB };

            DiaLibraryQueryGenerator.ComputeChimericScores(precursors, results, ppmTolerance: 10f);

            Assert.That(resultA.ChimericScore, Is.EqualTo(1.0f).Within(1e-4f),
                "Non-overlapping RT → no contestation despite shared fragment m/z");
            Assert.That(resultB.ChimericScore, Is.EqualTo(1.0f).Within(1e-4f),
                "Non-overlapping RT → no contestation despite shared fragment m/z");
        }

        [Test]
        [Description("GREEN: PrecursorIndex=-1 → ChimericScore=1.0 (early-out behavior is well-defined)")]
        public void ComputeChimericScores_NegativePrecursorIndex_ScoreIsOne()
        {
            // This is the CURRENT behavior for all results from AssembleResultsWithTemporalScoring.
            // After the fix, no results should have PrecursorIndex=-1. But the fallback
            // must remain defined and return 1.0 (no evidence of contamination = neutral).
            var precA = MakePrecursor(550.3, new float[] { 100.0f, 200.0f });
            var precursors = new List<LibraryPrecursorInput> { precA };

            var result = new DiaSearchResult(
                sequence: precA.Sequence, chargeState: 2, precursorMz: precA.PrecursorMz,
                windowId: 0, isDecoy: false, fragmentsQueried: 2,
                libraryRetentionTime: 30.0, rtWindowStart: 28f, rtWindowEnd: 32f);
            // PrecursorIndex deliberately left at default -1
            result.ExtractedIntensities[0] = 500f;
            result.ExtractedIntensities[1] = 300f;

            var results = new List<DiaSearchResult> { result };

            DiaLibraryQueryGenerator.ComputeChimericScores(precursors, results, ppmTolerance: 10f);

            Assert.That(result.ChimericScore, Is.EqualTo(1.0f).Within(1e-4f),
                "PrecursorIndex=-1 must produce ChimericScore=1.0 (neutral fallback)");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 3: End-to-end regression guard
        //  RED today (documents the bug). GREEN after the one-line fix.
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        [Description("REGRESSION GUARD: After fix, assembled results must have PrecursorIndex>=0 " +
                     "so ChimericScore reflects real interference, not always 1.0")]
        public void EndToEnd_AfterFix_ChimericScoreReflectsRealInterference()
        {
            // Two precursors in the same window with fully shared fragment m/z.
            // After the fix, both should get ChimericScore=0.0.
            // Before the fix, both get ChimericScore=1.0 (all early-out due to PrecursorIndex=-1).
            float[] sharedFragments = { 100.0f, 200.0f, 300.0f };

            var precA = MakePrecursor(550.3, sharedFragments, sequence: "PEPTIDEK");
            var precB = MakePrecursor(560.3, sharedFragments, sequence: "KEPTIDEM");
            var precursors = new List<LibraryPrecursorInput> { precA, precB };

            // Build a GenResult with both precursors assigned to window 0
            var groupA = new DiaLibraryQueryGenerator.PrecursorQueryGroup(
                inputIndex: 0, queryOffset: 0, queryCount: precA.FragmentCount,
                windowId: 0, rtMin: 28f, rtMax: 32f);
            var groupB = new DiaLibraryQueryGenerator.PrecursorQueryGroup(
                inputIndex: 1, queryOffset: precA.FragmentCount, queryCount: precB.FragmentCount,
                windowId: 0, rtMin: 28f, rtMax: 32f);

            var genResult = new DiaLibraryQueryGenerator.GenerationResult(
                queries: new FragmentQuery[precA.FragmentCount + precB.FragmentCount],
                precursorGroups: new[] { groupA, groupB },
                skippedNoWindow: 0,
                skippedNoFragments: 0);

            int totalFrags = precA.FragmentCount + precB.FragmentCount;
            var fragResults = new FragmentResult[totalFrags];
            for (int i = 0; i < totalFrags; i++)
                fragResults[i] = new FragmentResult(
                    queryId: i, dataPointCount: 0,
                    rtBufferOffset: 0, intensityBufferOffset: 0, totalIntensity: 0f);

            var extractionResult = new ExtractionResult(
                results: fragResults,
                rtBuffer: System.Array.Empty<float>(),
                intensityBuffer: System.Array.Empty<float>(),
                totalDataPoints: 0);

            var parameters = new DiaSearchParameters
            {
                MinFragmentsRequired = 0,
                MinScoreThreshold = 0f,
                PpmTolerance = 10f
            };

            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursors, genResult, extractionResult, parameters);

            Assert.That(results.Count, Is.EqualTo(2), "Both precursors must assemble");

            // After fix: PrecursorIndex is set correctly
            Assert.That(results[0].PrecursorIndex, Is.EqualTo(0),
                "After fix: first result must have PrecursorIndex=0");
            Assert.That(results[1].PrecursorIndex, Is.EqualTo(1),
                "After fix: second result must have PrecursorIndex=1");

            // Manually set intensities so ChimericScore has something to compute
            results[0].ExtractedIntensities[0] = 500f;
            results[0].ExtractedIntensities[1] = 300f;
            results[0].ExtractedIntensities[2] = 100f;
            results[1].ExtractedIntensities[0] = 400f;
            results[1].ExtractedIntensities[1] = 200f;
            results[1].ExtractedIntensities[2] = 100f;

            DiaLibraryQueryGenerator.ComputeChimericScores(precursors, results, ppmTolerance: 10f);

            // After fix: both scores reflect real interference (shared m/z → 0.0)
            Assert.That(results[0].ChimericScore, Is.EqualTo(0.0f).Within(1e-4f),
                "After fix: fully shared fragments → ChimericScore must be 0.0, not 1.0");
            Assert.That(results[1].ChimericScore, Is.EqualTo(0.0f).Within(1e-4f),
                "After fix: fully shared fragments → ChimericScore must be 0.0, not 1.0");
        }
    }
}