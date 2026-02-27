// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;
using System.Collections.Generic;
using NUnit.Framework;

namespace MassSpectrometry.Dia.Tests
{
    /// <summary>
    /// Tests for DiaTemporalScorer — the RT-resolved scoring engine.
    /// 
    /// These tests validate that:
    ///   1. Each scoring strategy produces correct results for known inputs
    ///   2. Temporal strategies outscore summed strategy for clean chromatographic peaks
    ///   3. Edge cases (no data, single time point, missing fragments) are handled
    ///   4. The RT alignment logic correctly merges fragments with slightly different RT grids
    ///   5. The nonlinear transform has the expected effect on scores
    /// </summary>
    [TestFixture]
    public class DiaTemporalScorerTests
    {
        // ─────────────────────────────────────────────────────────────────────
        //  Test helper: builds ExtractionResult-style buffers from a
        //  time × fragment intensity matrix plus RT values
        // ─────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Builds FragmentResult[], rtBuffer, and intensityBuffer from a
        /// matrix[time][fragment] layout. This simulates what the extraction
        /// engine would produce.
        /// 
        /// In the real engine, each fragment has its own XIC (RT[], Intensity[]).
        /// All fragments for a precursor share the same scans (same window),
        /// so their RT arrays are identical.
        /// </summary>
        private static (FragmentResult[] results, float[] rtBuffer, float[] intensityBuffer)
            BuildXicBuffers(float[] rts, float[,] intensityMatrix)
        {
            int timePoints = intensityMatrix.GetLength(0);
            int fragmentCount = intensityMatrix.GetLength(1);

            // Layout: fragment 0's XIC first, then fragment 1's, etc.
            // Each fragment has timePoints entries in rtBuffer and intensityBuffer.
            int totalPoints = timePoints * fragmentCount;
            float[] rtBuffer = new float[totalPoints];
            float[] intensityBuffer = new float[totalPoints];
            FragmentResult[] results = new FragmentResult[fragmentCount];

            int offset = 0;
            for (int f = 0; f < fragmentCount; f++)
            {
                float totalIntensity = 0f;
                int nonZeroPoints = 0;

                // Only write data points where intensity > 0 (simulates real extraction
                // where a fragment may not be found in every scan)
                int fragOffset = offset;
                for (int t = 0; t < timePoints; t++)
                {
                    float intensity = intensityMatrix[t, f];
                    if (intensity > 0f)
                    {
                        rtBuffer[offset] = rts[t];
                        intensityBuffer[offset] = intensity;
                        totalIntensity += intensity;
                        nonZeroPoints++;
                        offset++;
                    }
                }

                results[f] = new FragmentResult(
                    queryId: f,
                    dataPointCount: nonZeroPoints,
                    rtBufferOffset: fragOffset,
                    intensityBufferOffset: fragOffset, // aligned 1:1 with RT
                    totalIntensity: totalIntensity);
            }

            // Trim buffers to actual size
            float[] trimmedRt = new float[offset];
            float[] trimmedInt = new float[offset];
            Array.Copy(rtBuffer, trimmedRt, offset);
            Array.Copy(intensityBuffer, trimmedInt, offset);

            return (results, trimmedRt, trimmedInt);
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Test 1: Summed scoring produces same result as original IScorer
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void SummedScoringMatchesOriginalBehavior()
        {
            // 5 time points, 3 fragments
            // Library: [100, 200, 300]
            float[] library = { 100f, 200f, 300f };

            float[] rts = { 10f, 11f, 12f, 13f, 14f };
            float[,] matrix = new float[,]
            {
                { 10f, 20f, 30f },  // t=10: proportional to library
                { 20f, 40f, 60f },  // t=11: same ratio, 2x intensity
                { 30f, 60f, 90f },  // t=12: apex, same ratio, 3x
                { 20f, 40f, 60f },  // t=13: falling
                { 10f, 20f, 30f },  // t=14: back to baseline
            };

            var (results, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            // Summed intensities: [90, 180, 270] — exactly 0.9 × library
            // So cosine between [100,200,300] and [90,180,270] should be 1.0
            // (cosine is scale-invariant)

            var scorer = new DiaTemporalScorer(ScoringStrategy.Summed);
            var score = scorer.ScorePrecursor(library, 3, results, rtBuf, intBuf);

            Assert.That(score.IsValid, Is.True, "Score should be valid");
            Assert.That(score.DotProductScore, Is.EqualTo(1.0f).Within(0.001f),
                "Summed scoring of proportional intensities should give cosine ≈ 1.0");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Test 2: Consensus apex scoring
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void ConsensusApexScoringFindsCorrectApex()
        {
            float[] library = { 100f, 200f, 300f };

            float[] rts = { 10f, 11f, 12f, 13f, 14f };
            float[,] matrix = new float[,]
            {
                { 5f,  10f,  15f },   // t=10: low signal
                { 50f, 100f, 150f },  // t=11: rising
                { 100f, 200f, 300f }, // t=12: APEX — exactly matches library
                { 50f, 100f, 150f },  // t=13: falling
                { 5f,  10f,  15f },   // t=14: low
            };

            var (results, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            var scorer = new DiaTemporalScorer(ScoringStrategy.ConsensusApex);
            var score = scorer.ScorePrecursor(library, 3, results, rtBuf, intBuf);

            Assert.That(score.IsValid, Is.True);
            Assert.That(score.ApexTimeIndex, Is.EqualTo(2),
                "Apex should be at index 2 (highest total intensity)");
            Assert.That(score.DotProductScore, Is.EqualTo(1.0f).Within(0.001f),
                "Apex vector exactly matches library, so cosine should be 1.0");
        }

        [Test]
        public void ConsensusApexWithInterference_ScoresBetterThanSummed()
        {
            // Library: [100, 200, 300]
            // At the apex (t=12), fragments match the library pattern.
            // But at other times, fragment 0 has interference (extra signal)
            // that distorts the summed ratio.
            float[] library = { 100f, 200f, 300f };

            float[] rts = { 10f, 11f, 12f, 13f, 14f };
            float[,] matrix = new float[,]
            {
                { 500f, 10f,  15f },   // t=10: fragment 0 has interference
                { 500f, 100f, 150f },  // t=11: interference continues
                { 100f, 200f, 300f },  // t=12: APEX — clean, matches library
                { 500f, 100f, 150f },  // t=13: interference again
                { 500f, 10f,  15f },   // t=14: interference
            };

            var (results, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            var summedScorer = new DiaTemporalScorer(ScoringStrategy.Summed);
            var apexScorer = new DiaTemporalScorer(ScoringStrategy.ConsensusApex);

            var summedScore = summedScorer.ScorePrecursor(library, 3, results, rtBuf, intBuf);
            var apexScore = apexScorer.ScorePrecursor(library, 3, results, rtBuf, intBuf);

            Assert.That(apexScore.DotProductScore, Is.GreaterThan(summedScore.DotProductScore),
                "Apex scoring should beat summed scoring when interference is present at non-apex times");

            // The apex itself is clean and matches the library perfectly
            Assert.That(apexScore.DotProductScore, Is.EqualTo(1.0f).Within(0.001f));
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Test 3: Temporal cosine scoring
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void TemporalCosine_PerfectMatch_ReturnsOne()
        {
            // Every time point has the same ratio as the library
            float[] library = { 100f, 200f, 300f };

            float[] rts = { 10f, 11f, 12f };
            float[,] matrix = new float[,]
            {
                { 10f, 20f, 30f },
                { 50f, 100f, 150f },
                { 100f, 200f, 300f },
            };

            var (results, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            var scorer = new DiaTemporalScorer(ScoringStrategy.TemporalCosine);
            var score = scorer.ScorePrecursor(library, 3, results, rtBuf, intBuf);

            Assert.That(score.DotProductScore, Is.EqualTo(1.0f).Within(0.001f),
                "All time points match library ratio → average cosine should be 1.0");
            Assert.That(score.TimePointsUsed, Is.EqualTo(3));
        }

        [Test]
        public void TemporalCosine_WithInterference_ScoresBetterThanSummed()
        {
            // The key test: temporal scoring should down-weight time points
            // where interference distorts the fragment ratios.
            float[] library = { 100f, 200f, 300f };

            // 10 time points. The middle 4 are clean (match library).
            // The outer 6 have interference on fragment 0.
            float[] rts = new float[10];
            for (int i = 0; i < 10; i++) rts[i] = 10f + i;

            float[,] matrix = new float[10, 3];
            for (int t = 0; t < 10; t++)
            {
                if (t >= 3 && t <= 6) // Clean peak region
                {
                    float scale = (t == 4 || t == 5) ? 3.0f : 1.5f;
                    matrix[t, 0] = 100f * scale;
                    matrix[t, 1] = 200f * scale;
                    matrix[t, 2] = 300f * scale;
                }
                else // Interference region
                {
                    matrix[t, 0] = 800f; // strong interference on fragment 0
                    matrix[t, 1] = 20f;
                    matrix[t, 2] = 30f;
                }
            }

            var (results, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            var summedScorer = new DiaTemporalScorer(ScoringStrategy.Summed);
            var temporalScorer = new DiaTemporalScorer(ScoringStrategy.TemporalCosine);

            var summedScore = summedScorer.ScorePrecursor(library, 3, results, rtBuf, intBuf);
            var temporalScore = temporalScorer.ScorePrecursor(library, 3, results, rtBuf, intBuf);

            Assert.That(temporalScore.DotProductScore, Is.GreaterThan(summedScore.DotProductScore),
                $"Temporal cosine ({temporalScore.DotProductScore:F4}) should exceed summed ({summedScore.DotProductScore:F4}) when interference affects some but not all time points");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Test 4: Weighted temporal cosine with nonlinear transform
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void WeightedTemporalWithTransform_AmplifiesHighScores()
        {
            float[] library = { 100f, 200f, 300f };

            float[] rts = { 10f, 11f, 12f };
            float[,] matrix = new float[,]
            {
                { 10f, 20f, 30f },    // Perfect ratio, low intensity
                { 50f, 100f, 150f },  // Perfect ratio, medium intensity
                { 100f, 200f, 300f }, // Perfect ratio, high intensity
            };

            var (results, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            var noTransform = new DiaTemporalScorer(ScoringStrategy.TemporalCosine);
            var withTransform = new DiaTemporalScorer(ScoringStrategy.WeightedTemporalCosineWithTransform, 3.0f);

            var noTransformScore = noTransform.ScorePrecursor(library, 3, results, rtBuf, intBuf);
            var withTransformScore = withTransform.ScorePrecursor(library, 3, results, rtBuf, intBuf);

            // When all time points have cosine=1.0, the transform doesn't change the score
            // (1.0^3 = 1.0). But raw cosine values from weighting might differ slightly.
            Assert.That(noTransformScore.DotProductScore, Is.EqualTo(1.0f).Within(0.001f));
            Assert.That(withTransformScore.DotProductScore, Is.EqualTo(1.0f).Within(0.001f));
        }

        [Test]
        public void NonlinearTransform_SuppressesModerateCosine()
        {
            // Create a scenario where temporal cosine is moderate (~0.7)
            // and verify the transform suppresses it
            float[] library = { 100f, 200f, 300f };

            float[] rts = { 10f, 11f, 12f };
            float[,] matrix = new float[,]
            {
                { 100f, 200f, 300f },  // Perfect match
                { 300f, 200f, 100f },  // Reversed — poor match
                { 100f, 200f, 300f },  // Perfect match
            };

            var (results, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            var temporal = new DiaTemporalScorer(ScoringStrategy.TemporalCosine);
            var transformed = new DiaTemporalScorer(ScoringStrategy.WeightedTemporalCosineWithTransform, 3.0f);

            var rawScore = temporal.ScorePrecursor(library, 3, results, rtBuf, intBuf);
            var transScore = transformed.ScorePrecursor(library, 3, results, rtBuf, intBuf);

            // Raw score should be moderate (average of high and low cosines)
            Assert.That(rawScore.DotProductScore, Is.GreaterThan(0.3f).And.LessThan(0.95f),
                "Raw temporal cosine should be moderate due to mixed time points");

            // Transformed score should be lower (cos^3 suppresses moderate values)
            Assert.That(transScore.DotProductScore, Is.LessThan(rawScore.DotProductScore),
                "Nonlinear transform should suppress moderate cosine values");

            // The raw cosine should be stored
            Assert.That(transScore.RawCosine, Is.GreaterThan(transScore.DotProductScore),
                "RawCosine should be higher than the transformed DotProductScore");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Test 5: Edge cases
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void InsufficientData_ReturnsNaN()
        {
            // Only 1 fragment (minimum is 2)
            float[] library = { 100f };
            float[] rts = { 10f };
            float[,] matrix = new float[,] { { 100f } };

            var (results, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            var scorer = new DiaTemporalScorer(ScoringStrategy.TemporalCosine);
            var score = scorer.ScorePrecursor(library, 1, results, rtBuf, intBuf);

            Assert.That(score.IsValid, Is.False, "Should return NaN for < 2 fragments");
        }

        [Test]
        public void NoXicData_ReturnsInsufficientScore()
        {
            float[] library = { 100f, 200f, 300f };

            // All fragments have 0 data points
            var results = new FragmentResult[]
            {
                new FragmentResult(0, 0, 0, 0, 0f),
                new FragmentResult(1, 0, 0, 0, 0f),
                new FragmentResult(2, 0, 0, 0, 0f),
            };

            var scorer = new DiaTemporalScorer(ScoringStrategy.TemporalCosine);
            var score = scorer.ScorePrecursor(library, 3, results, Array.Empty<float>(), Array.Empty<float>());

            Assert.That(score.IsValid, Is.False);
        }

        [Test]
        public void SingleTimePoint_StillScoresCorrectly()
        {
            float[] library = { 100f, 200f, 300f };

            float[] rts = { 12f };
            float[,] matrix = new float[,]
            {
                { 100f, 200f, 300f },
            };

            var (results, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            var scorer = new DiaTemporalScorer(ScoringStrategy.TemporalCosine);
            var score = scorer.ScorePrecursor(library, 3, results, rtBuf, intBuf);

            Assert.That(score.IsValid, Is.True);
            Assert.That(score.DotProductScore, Is.EqualTo(1.0f).Within(0.001f));
            Assert.That(score.TimePointsUsed, Is.EqualTo(1));
        }

        [Test]
        public void MissingFragment_HandledGracefully()
        {
            // Fragment 1 has no XIC data (not detected)
            float[] library = { 100f, 200f, 300f };

            float[] rts = { 10f, 11f, 12f };
            // Fragment 1 (index 1) is zero everywhere → simulates "not detected"
            float[,] matrix = new float[,]
            {
                { 100f, 0f, 300f },
                { 100f, 0f, 300f },
                { 100f, 0f, 300f },
            };

            var (results, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            var scorer = new DiaTemporalScorer(ScoringStrategy.TemporalCosine);
            var score = scorer.ScorePrecursor(library, 3, results, rtBuf, intBuf);

            Assert.That(score.IsValid, Is.True, "Should still score even with one missing fragment");
            // The cosine of [100, 0, 300] vs [100, 200, 300] should be > 0 but < 1
            Assert.That(score.DotProductScore, Is.GreaterThan(0.5f).And.LessThan(1.0f));
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Test 6: Strategy comparison (the key validation)
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void StrategyComparison_TemporalBeatssSummedWithInterference()
        {
            // Simulate a realistic scenario: a peptide peak with interference
            // from a co-eluting peptide on fragment 0. The interference
            // affects scans outside the peak region.
            float[] library = { 100f, 500f, 300f, 200f, 150f, 80f };

            // 20 time points, peak centered at index 10
            int timeCount = 20;
            float[] rts = new float[timeCount];
            for (int i = 0; i < timeCount; i++) rts[i] = 20f + i * 0.1f;

            float[,] matrix = new float[timeCount, 6];
            for (int t = 0; t < timeCount; t++)
            {
                // Gaussian-shaped peak centered at t=10
                float peakShape = MathF.Exp(-0.5f * (t - 10f) * (t - 10f) / (2f * 2f));
                for (int f = 0; f < 6; f++)
                    matrix[t, f] = library[f] * peakShape * 10f; // Scale up

                // Add interference to fragment 0 at non-peak times
                if (t < 5 || t > 15)
                    matrix[t, 0] += 2000f; // Strong interference
            }

            var (results, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            var summed = new DiaTemporalScorer(ScoringStrategy.Summed);
            var apex = new DiaTemporalScorer(ScoringStrategy.ConsensusApex);
            var temporal = new DiaTemporalScorer(ScoringStrategy.TemporalCosine);
            var weighted = new DiaTemporalScorer(ScoringStrategy.WeightedTemporalCosineWithTransform, 3.0f);

            var summedScore = summed.ScorePrecursor(library, 6, results, rtBuf, intBuf);
            var apexScore = apex.ScorePrecursor(library, 6, results, rtBuf, intBuf);
            var temporalScore = temporal.ScorePrecursor(library, 6, results, rtBuf, intBuf);
            var weightedScore = weighted.ScorePrecursor(library, 6, results, rtBuf, intBuf);

            Console.WriteLine($"Strategy comparison with interference:");
            Console.WriteLine($"  Summed:            {summedScore.DotProductScore:F4}");
            Console.WriteLine($"  ConsensusApex:     {apexScore.DotProductScore:F4}");
            Console.WriteLine($"  TemporalCosine:    {temporalScore.DotProductScore:F4}");
            Console.WriteLine($"  WeightedTemporal:  {weightedScore.DotProductScore:F4} (raw={weightedScore.RawCosine:F4})");

            // All temporal strategies should beat summed
            Assert.That(apexScore.DotProductScore, Is.GreaterThan(summedScore.DotProductScore),
                "Apex should beat summed");
            Assert.That(temporalScore.DotProductScore, Is.GreaterThan(summedScore.DotProductScore),
                "Temporal cosine should beat summed");

            // Apex should be very high because the peak center is clean
            Assert.That(apexScore.DotProductScore, Is.GreaterThan(0.95f),
                "Apex should be near-perfect at peak center");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Test 7: RT alignment with slightly mismatched grids
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void RtAlignment_HandlesSlightRtOffsets()
        {
            // Simulate fragments from the same window but with tiny RT offsets
            // (e.g., due to different DIA windows in the same cycle)
            float[] library = { 100f, 200f, 300f };

            // Fragment 0: RTs at exactly [10, 11, 12]
            // Fragment 1: RTs at [10.001, 11.001, 12.001] (0.001 min offset)
            // Fragment 2: RTs at [10.005, 11.005, 12.005] (0.005 min offset)
            // All within the 0.01 min alignment tolerance

            // Build manually to control per-fragment RT values
            float[] rtBuffer = new float[]
            {
                10.0f, 11.0f, 12.0f,       // fragment 0
                10.001f, 11.001f, 12.001f,  // fragment 1
                10.005f, 11.005f, 12.005f,  // fragment 2
            };
            float[] intBuffer = new float[]
            {
                100f, 100f, 100f,  // fragment 0: constant
                200f, 200f, 200f,  // fragment 1: constant
                300f, 300f, 300f,  // fragment 2: constant
            };

            var results = new FragmentResult[]
            {
                new FragmentResult(0, 3, 0, 0, 300f),
                new FragmentResult(1, 3, 3, 3, 600f),
                new FragmentResult(2, 3, 6, 6, 900f),
            };

            var scorer = new DiaTemporalScorer(ScoringStrategy.TemporalCosine);
            var score = scorer.ScorePrecursor(library, 3, results, rtBuffer, intBuffer);

            Assert.That(score.IsValid, Is.True);
            Assert.That(score.DotProductScore, Is.EqualTo(1.0f).Within(0.01f),
                "Slight RT offsets should be handled by alignment tolerance");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Test 8: AssembleResultsWithTemporalScoring integration
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void AssembleResultsWithTemporalScoring_ProducesResults()
        {
            // Build a minimal but complete pipeline test:
            // 1 precursor, 3 fragments, 5 time points
            var precursors = new List<LibraryPrecursorInput>
            {
                new LibraryPrecursorInput(
                    sequence: "PEPTIDE",
                    precursorMz: 500.0,
                    chargeState: 2,
                    retentionTime: 12.0,
                    isDecoy: false,
                    fragmentMzs: new float[] { 200f, 300f, 400f },
                    fragmentIntensities: new float[] { 100f, 200f, 300f })
            };

            // Simulate extraction results with 5 time points per fragment
            float[] rts = { 10f, 11f, 12f, 13f, 14f };
            float[,] matrix = new float[,]
            {
                { 10f, 20f, 30f },
                { 50f, 100f, 150f },
                { 100f, 200f, 300f }, // apex
                { 50f, 100f, 150f },
                { 10f, 20f, 30f },
            };

            var (fragResults, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            // Build the PrecursorQueryGroup
            var groups = new DiaLibraryQueryGenerator.PrecursorQueryGroup[]
            {
                new DiaLibraryQueryGenerator.PrecursorQueryGroup(
                    inputIndex: 0, queryOffset: 0, queryCount: 3,
                    windowId: 0, rtMin: 10f, rtMax: 14f)
            };

            var genResult = new DiaLibraryQueryGenerator.GenerationResult(
                queries: new FragmentQuery[3], // Not used by AssembleResults
                precursorGroups: groups,
                skippedNoWindow: 0,
                skippedNoFragments: 0);

            var extractionResult = new ExtractionResult(fragResults, rtBuf, intBuf, rtBuf.Length);

            var parameters = new DiaSearchParameters
            {
                ScoringStrategy = ScoringStrategy.TemporalCosine,
                MinFragmentsRequired = 2
            };

            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursors, genResult, extractionResult, parameters);

            Assert.That(results.Count, Is.EqualTo(1));
            var r = results[0];
            Assert.That(r.Sequence, Is.EqualTo("PEPTIDE"));
            Assert.That(r.DotProductScore, Is.EqualTo(1.0f).Within(0.01f),
                "All time points have perfect library match");
            Assert.That(r.ScoringStrategyUsed, Is.EqualTo(ScoringStrategy.TemporalCosine));
            Assert.That(r.TimePointsUsed, Is.GreaterThan(0));
            Assert.That(r.FragmentsDetected, Is.EqualTo(3));
            Assert.That(!float.IsNaN(r.SpectralAngleScore), "Spectral angle should be computed");
            Assert.That(!float.IsNaN(r.RawCosine), "RawCosine should be populated");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Test 9: Verify all strategies work through the assembly pipeline
        // ─────────────────────────────────────────────────────────────────────

        [TestCase(ScoringStrategy.Summed)]
        [TestCase(ScoringStrategy.ConsensusApex)]
        [TestCase(ScoringStrategy.TemporalCosine)]
        [TestCase(ScoringStrategy.WeightedTemporalCosineWithTransform)]
        public void AllStrategiesWork_ThroughAssemblyPipeline(ScoringStrategy strategy)
        {
            var precursors = new List<LibraryPrecursorInput>
            {
                new LibraryPrecursorInput("PEPTIDE", 500.0, 2, 12.0, false,
                    new float[] { 200f, 300f, 400f },
                    new float[] { 100f, 200f, 300f })
            };

            float[] rts = { 10f, 11f, 12f };
            float[,] matrix = new float[,]
            {
                { 100f, 200f, 300f },
                { 100f, 200f, 300f },
                { 100f, 200f, 300f },
            };

            var (fragResults, rtBuf, intBuf) = BuildXicBuffers(rts, matrix);

            var groups = new DiaLibraryQueryGenerator.PrecursorQueryGroup[]
            {
                new DiaLibraryQueryGenerator.PrecursorQueryGroup(0, 0, 3, 0, 10f, 12f)
            };
            var genResult = new DiaLibraryQueryGenerator.GenerationResult(
                new FragmentQuery[3], groups, 0, 0);
            var extractionResult = new ExtractionResult(fragResults, rtBuf, intBuf, rtBuf.Length);

            var parameters = new DiaSearchParameters
            {
                ScoringStrategy = strategy,
                MinFragmentsRequired = 2,
                NonlinearPower = 3.0f
            };

            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursors, genResult, extractionResult, parameters);

            Assert.That(results.Count, Is.EqualTo(1), $"Strategy {strategy} should produce 1 result");
            Assert.That(results[0].DotProductScore, Is.EqualTo(1.0f).Within(0.01f),
                $"Strategy {strategy} should give ~1.0 for perfect match");
            Assert.That(results[0].ScoringStrategyUsed, Is.EqualTo(strategy));
        }
    }
}