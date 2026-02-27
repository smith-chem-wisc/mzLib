// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Lightweight input describing one library precursor for DIA query generation.
    /// This struct exists so that DiaLibraryQueryGenerator has NO dependency on Omics
    /// (LibrarySpectrum, MatchedFragmentIon, etc.), avoiding circular project references.
    /// 
    /// The caller (MetaMorpheus DiaEngine facade) converts LibrarySpectrum → LibraryPrecursorInput.
    /// </summary>
    public readonly struct LibraryPrecursorInput
    {
        /// <summary>Peptide sequence with modifications (for result reporting)</summary>
        public readonly string Sequence;

        /// <summary>Precursor m/z (used to find the DIA isolation window)</summary>
        public readonly double PrecursorMz;

        /// <summary>Precursor charge state</summary>
        public readonly int ChargeState;

        /// <summary>
        /// Predicted/library retention time in minutes. Null if unavailable,
        /// in which case the full run RT range is used as fallback.
        /// </summary>
        public readonly double? RetentionTime;

        /// <summary>Whether this is a decoy entry</summary>
        public readonly bool IsDecoy;

        /// <summary>
        /// Indexed retention time (iRT) value from the library, if available.
        /// Used by the iRT calibration system (IrtLibraryIndex) to build
        /// calibration models that map iRT → experimental RT.
        /// Null if the library entry has no iRT annotation.
        /// </summary>
        public readonly double? IrtValue;

        /// <summary>
        /// Fragment ion m/z values from the library spectrum.
        /// Length determines the number of queries generated for this precursor.
        /// </summary>
        public readonly float[] FragmentMzs;

        /// <summary>
        /// Fragment ion intensities from the library spectrum (parallel to FragmentMzs).
        /// Used downstream for scoring (library vs extracted intensity comparison).
        /// </summary>
        public readonly float[] FragmentIntensities;

        public LibraryPrecursorInput(
            string sequence,
            double precursorMz,
            int chargeState,
            double? retentionTime,
            bool isDecoy,
            float[] fragmentMzs,
            float[] fragmentIntensities,
            double? irtValue = null)
        {
            Sequence = sequence ?? throw new ArgumentNullException(nameof(sequence));
            PrecursorMz = precursorMz;
            ChargeState = chargeState;
            RetentionTime = retentionTime;
            IsDecoy = isDecoy;
            IrtValue = irtValue;
            FragmentMzs = fragmentMzs ?? throw new ArgumentNullException(nameof(fragmentMzs));
            FragmentIntensities = fragmentIntensities ?? throw new ArgumentNullException(nameof(fragmentIntensities));

            if (fragmentMzs.Length != fragmentIntensities.Length)
                throw new ArgumentException(
                    $"FragmentMzs length ({fragmentMzs.Length}) must equal FragmentIntensities length ({fragmentIntensities.Length})");
        }

        /// <summary>Number of fragment ions in this library entry</summary>
        public int FragmentCount => FragmentMzs.Length;
    }

    /// <summary>
    /// Converts library precursor inputs into FragmentQuery[] for the DIA extraction engine.
    /// 
    /// For each LibraryPrecursorInput:
    ///   1. Maps precursor m/z → window ID via DiaScanIndex.FindWindowForPrecursorMz()
    ///   2. Computes RT window from RetentionTime ± tolerance
    ///   3. Generates one FragmentQuery per fragment ion
    /// 
    /// Also produces PrecursorQueryGroup[] metadata that records which queries belong to
    /// which precursor, enabling downstream result assembly.
    /// 
    /// This class has NO dependency on Omics (LibrarySpectrum, MatchedFragmentIon, etc.).
    /// </summary>
    public static class DiaLibraryQueryGenerator
    {
        /// <summary>
        /// Metadata linking a precursor to its generated queries.
        /// Used after extraction to assemble DiaSearchResult objects.
        /// </summary>
        public readonly struct PrecursorQueryGroup
        {
            /// <summary>Index into the original input list</summary>
            public readonly int InputIndex;

            /// <summary>Start index into the FragmentQuery[] array</summary>
            public readonly int QueryOffset;

            /// <summary>Number of fragment queries for this precursor</summary>
            public readonly int QueryCount;

            /// <summary>Window ID this precursor maps to</summary>
            public readonly int WindowId;

            /// <summary>RT window lower bound used for queries</summary>
            public readonly float RtMin;

            /// <summary>RT window upper bound used for queries</summary>
            public readonly float RtMax;

            public PrecursorQueryGroup(int inputIndex, int queryOffset, int queryCount,
                int windowId, float rtMin, float rtMax)
            {
                InputIndex = inputIndex;
                QueryOffset = queryOffset;
                QueryCount = queryCount;
                WindowId = windowId;
                RtMin = rtMin;
                RtMax = rtMax;
            }
        }

        /// <summary>
        /// Result of query generation: the queries themselves plus grouping metadata.
        /// </summary>
        public readonly struct GenerationResult
        {
            public readonly FragmentQuery[] Queries;
            public readonly PrecursorQueryGroup[] PrecursorGroups;
            public readonly int SkippedNoWindow;
            public readonly int SkippedNoFragments;

            public GenerationResult(FragmentQuery[] queries, PrecursorQueryGroup[] precursorGroups,
                int skippedNoWindow, int skippedNoFragments)
            {
                Queries = queries;
                PrecursorGroups = precursorGroups;
                SkippedNoWindow = skippedNoWindow;
                SkippedNoFragments = skippedNoFragments;
            }
        }

        /// <summary>
        /// Generates FragmentQuery[] from library precursor inputs.
        /// Two-pass: first counts for a single allocation, then fills.
        /// </summary>
        public static GenerationResult Generate(
            IList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaSearchParameters parameters)
        {
            if (precursors == null) throw new ArgumentNullException(nameof(precursors));
            if (scanIndex == null) throw new ArgumentNullException(nameof(scanIndex));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));

            float rtTolerance = parameters.RtToleranceMinutes;
            float ppmTolerance = parameters.PpmTolerance;

            // Fallback RT bounds for precursors with null RT
            float globalRtMin = scanIndex.GetGlobalRtMin();
            float globalRtMax = scanIndex.GetGlobalRtMax();

            // First pass: count
            int totalQueryCount = 0;
            int skippedNoWindow = 0;
            int skippedNoFragments = 0;

            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0)
                {
                    skippedNoWindow++;
                    continue;
                }
                if (p.FragmentCount == 0)
                {
                    skippedNoFragments++;
                    continue;
                }
                totalQueryCount += p.FragmentCount;
            }

            // Allocate
            var queries = new FragmentQuery[totalQueryCount];
            var groups = new List<PrecursorQueryGroup>(
                precursors.Count - skippedNoWindow - skippedNoFragments);

            // Second pass: fill
            int queryIndex = 0;
            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0 || p.FragmentCount == 0)
                    continue;

                float rtMin, rtMax;
                if (p.RetentionTime.HasValue)
                {
                    float rt = (float)p.RetentionTime.Value;
                    rtMin = rt - rtTolerance;
                    rtMax = rt + rtTolerance;
                }
                else
                {
                    rtMin = globalRtMin;
                    rtMax = globalRtMax;
                }

                int groupOffset = queryIndex;

                for (int f = 0; f < p.FragmentCount; f++)
                {
                    queries[queryIndex] = new FragmentQuery(
                        targetMz: p.FragmentMzs[f],
                        tolerancePpm: ppmTolerance,
                        rtMin: rtMin,
                        rtMax: rtMax,
                        windowId: windowId,
                        queryId: queryIndex
                    );
                    queryIndex++;
                }

                groups.Add(new PrecursorQueryGroup(
                    inputIndex: i,
                    queryOffset: groupOffset,
                    queryCount: p.FragmentCount,
                    windowId: windowId,
                    rtMin: rtMin,
                    rtMax: rtMax
                ));
            }

            return new GenerationResult(queries, groups.ToArray(), skippedNoWindow, skippedNoFragments);
        }

        /// <summary>
        /// Generates FragmentQuery[] using a calibrated RT model instead of a fixed RT tolerance.
        /// 
        /// For each precursor with an iRT value, the calibration model predicts the
        /// experimental RT and computes a per-precursor window based on the model's
        /// residual sigma: predicted_RT ± (sigma × CalibratedWindowSigmaMultiplier).
        /// 
        /// Precursors without iRT values fall back to the fixed RtToleranceMinutes window.
        /// This two-pass approach was validated to narrow windows by ~23× (from ±20 min
        /// to ±0.9 min) in synthetic benchmarks with iRT-based libraries.
        /// </summary>
        /// <param name="precursors">Library precursor inputs (must have IrtValue for calibration benefit).</param>
        /// <param name="scanIndex">The DIA scan index.</param>
        /// <param name="parameters">Search parameters (CalibratedWindowSigmaMultiplier controls window width).</param>
        /// <param name="calibrationModel">
        /// Fitted RT calibration model (iRT → experimental RT). Obtained from RtCalibrationFitter.
        /// </param>
        /// <returns>GenerationResult with per-precursor calibrated RT windows.</returns>
        public static GenerationResult GenerateCalibrated(
            IList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaSearchParameters parameters,
            RtCalibrationModel calibrationModel)
        {
            if (precursors == null) throw new ArgumentNullException(nameof(precursors));
            if (scanIndex == null) throw new ArgumentNullException(nameof(scanIndex));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (calibrationModel == null) throw new ArgumentNullException(nameof(calibrationModel));

            float ppmTolerance = parameters.PpmTolerance;
            double sigmaMultiplier = parameters.CalibratedWindowSigmaMultiplier;
            float fallbackRtTolerance = parameters.RtToleranceMinutes;

            // Calibrated window half-width based on model residual sigma (in minutes)
            float calibratedHalfWindow = (float)calibrationModel.GetMinutesWindowHalfWidth(sigmaMultiplier);

            // Fallback RT bounds for precursors with no RT info at all
            float globalRtMin = scanIndex.GetGlobalRtMin();
            float globalRtMax = scanIndex.GetGlobalRtMax();

            // First pass: count
            int totalQueryCount = 0;
            int skippedNoWindow = 0;
            int skippedNoFragments = 0;

            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0) { skippedNoWindow++; continue; }
                if (p.FragmentCount == 0) { skippedNoFragments++; continue; }
                totalQueryCount += p.FragmentCount;
            }

            // Allocate
            var queries = new FragmentQuery[totalQueryCount];
            var groups = new List<PrecursorQueryGroup>(
                precursors.Count - skippedNoWindow - skippedNoFragments);

            // Second pass: fill with calibrated RT windows
            int queryIndex = 0;
            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0 || p.FragmentCount == 0)
                    continue;

                float rtMin, rtMax;
                if (p.IrtValue.HasValue)
                {
                    // Use calibration model to predict experimental RT from iRT
                    float predictedRt = (float)calibrationModel.ToMinutes(p.IrtValue.Value);
                    rtMin = predictedRt - calibratedHalfWindow;
                    rtMax = predictedRt + calibratedHalfWindow;
                }
                else if (p.RetentionTime.HasValue)
                {
                    // No iRT, but have library RT. Use the calibration model to predict
                    // observed RT from library RT (the model was fitted on library RT vs
                    // observed RT, so library RT serves as the iRT input here).
                    float predictedRt = (float)calibrationModel.ToMinutes(p.RetentionTime.Value);
                    rtMin = predictedRt - calibratedHalfWindow;
                    rtMax = predictedRt + calibratedHalfWindow;
                }
                else
                {
                    rtMin = globalRtMin;
                    rtMax = globalRtMax;
                }

                int groupOffset = queryIndex;

                for (int f = 0; f < p.FragmentCount; f++)
                {
                    queries[queryIndex] = new FragmentQuery(
                        targetMz: p.FragmentMzs[f],
                        tolerancePpm: ppmTolerance,
                        rtMin: rtMin,
                        rtMax: rtMax,
                        windowId: windowId,
                        queryId: queryIndex
                    );
                    queryIndex++;
                }

                groups.Add(new PrecursorQueryGroup(
                    inputIndex: i,
                    queryOffset: groupOffset,
                    queryCount: p.FragmentCount,
                    windowId: windowId,
                    rtMin: rtMin,
                    rtMax: rtMax
                ));
            }

            return new GenerationResult(queries, groups.ToArray(), skippedNoWindow, skippedNoFragments);
        }

        /// <summary>
        /// Original result assembly using summed intensities and IScorer interface.
        /// Retained for backward compatibility and benchmarking.
        /// 
        /// Uses TotalIntensity from FragmentResult (collapses all XIC points to one value per fragment).
        /// </summary>
        public static List<DiaSearchResult> AssembleResults(
            IList<LibraryPrecursorInput> precursors,
            GenerationResult generationResult,
            FragmentResult[] extractionResults,
            DiaSearchParameters parameters,
            IScorer dotProductScorer = null,
            IScorer spectralAngleScorer = null)
        {
            var results = new List<DiaSearchResult>(generationResult.PrecursorGroups.Length);

            for (int g = 0; g < generationResult.PrecursorGroups.Length; g++)
            {
                var group = generationResult.PrecursorGroups[g];
                var input = precursors[group.InputIndex];

                var result = new DiaSearchResult(
                    sequence: input.Sequence,
                    chargeState: input.ChargeState,
                    precursorMz: input.PrecursorMz,
                    windowId: group.WindowId,
                    isDecoy: input.IsDecoy,
                    fragmentsQueried: group.QueryCount,
                    libraryRetentionTime: input.RetentionTime,
                    rtWindowStart: group.RtMin,
                    rtWindowEnd: group.RtMax
                );
                result.ScoringStrategyUsed = ScoringStrategy.Summed;

                int detected = 0;
                for (int f = 0; f < group.QueryCount; f++)
                {
                    int qi = group.QueryOffset + f;
                    var fr = extractionResults[qi];
                    result.ExtractedIntensities[f] = fr.TotalIntensity;
                    result.XicPointCounts[f] = fr.DataPointCount;
                    if (fr.DataPointCount > 0)
                        detected++;
                }
                result.FragmentsDetected = detected;

                if (!result.MeetsMinFragments(parameters.MinFragmentsRequired))
                    continue;

                if (dotProductScorer != null)
                {
                    result.DotProductScore = dotProductScorer.Score(
                        input.FragmentIntensities, result.ExtractedIntensities);
                    result.RawCosine = result.DotProductScore;
                }

                if (spectralAngleScorer != null)
                    result.SpectralAngleScore = spectralAngleScorer.Score(
                        input.FragmentIntensities, result.ExtractedIntensities);

                if (!float.IsNaN(result.DotProductScore) && result.DotProductScore < parameters.MinScoreThreshold)
                    continue;

                results.Add(result);
            }

            return results;
        }

        /// <summary>
        /// Temporal-scoring result assembly using RT-resolved XIC data.
        /// 
        /// This is the new scoring pipeline that replaces the summed-intensity approach.
        /// Instead of collapsing each fragment's XIC to a single summed value, it passes
        /// the full RT × Intensity data to DiaTemporalScorer which computes time-resolved
        /// cosine similarity.
        /// 
        /// The ExtractionResult contains per-fragment XIC buffers (RT[] and Intensity[])
        /// that are accessed via FragmentResult offsets. The temporal scorer aligns these
        /// XICs to a common time grid and applies the configured scoring strategy.
        /// 
        /// For backward compatibility, the summed ExtractedIntensities[] and XicPointCounts[]
        /// are still populated on each DiaSearchResult.
        /// </summary>
        /// <param name="precursors">The library precursor inputs (same as used for query generation).</param>
        /// <param name="generationResult">Output from Generate() — links precursors to query indices.</param>
        /// <param name="extractionResult">Full extraction output including XIC buffers.</param>
        /// <param name="parameters">Search parameters including scoring strategy.</param>
        /// <returns>Scored DiaSearchResult list, filtered by minimum fragments and score threshold.</returns>
        public static List<DiaSearchResult> AssembleResultsWithTemporalScoring(
            IList<LibraryPrecursorInput> precursors,
            GenerationResult generationResult,
            ExtractionResult extractionResult,
            DiaSearchParameters parameters)
        {
            if (precursors == null) throw new ArgumentNullException(nameof(precursors));
            if (extractionResult == null) throw new ArgumentNullException(nameof(extractionResult));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));

            // Create scorers for ALL strategies we want as features
            var apexScorer = new DiaTemporalScorer(ScoringStrategy.ConsensusApex);
            var temporalScorer = new DiaTemporalScorer(ScoringStrategy.TemporalCosine);

            // The "primary" scorer is whatever the user requested
            DiaTemporalScorer primaryScorer = parameters.ScoringStrategy switch
            {
                ScoringStrategy.ConsensusApex => apexScorer,
                ScoringStrategy.TemporalCosine => temporalScorer,
                ScoringStrategy.WeightedTemporalCosineWithTransform =>
                    new DiaTemporalScorer(ScoringStrategy.WeightedTemporalCosineWithTransform, parameters.NonlinearPower),
                _ => new DiaTemporalScorer(ScoringStrategy.Summed),
            };

            var results = new List<DiaSearchResult>(generationResult.PrecursorGroups.Length);

            ReadOnlySpan<float> rtBuffer = extractionResult.RtBuffer.AsSpan();
            ReadOnlySpan<float> intensityBuffer = extractionResult.IntensityBuffer.AsSpan();

            for (int g = 0; g < generationResult.PrecursorGroups.Length; g++)
            {
                var group = generationResult.PrecursorGroups[g];
                var input = precursors[group.InputIndex];

                var result = new DiaSearchResult(
                    sequence: input.Sequence,
                    chargeState: input.ChargeState,
                    precursorMz: input.PrecursorMz,
                    windowId: group.WindowId,
                    isDecoy: input.IsDecoy,
                    fragmentsQueried: group.QueryCount,
                    libraryRetentionTime: input.RetentionTime,
                    rtWindowStart: group.RtMin,
                    rtWindowEnd: group.RtMax
                );
                result.ScoringStrategyUsed = parameters.ScoringStrategy;

                // Populate summed intensities and XIC point counts (backward compat + diagnostics)
                int detected = 0;
                for (int f = 0; f < group.QueryCount; f++)
                {
                    int qi = group.QueryOffset + f;
                    var fr = extractionResult.Results[qi];
                    result.ExtractedIntensities[f] = fr.TotalIntensity;
                    result.XicPointCounts[f] = fr.DataPointCount;
                    if (fr.DataPointCount > 0)
                        detected++;
                }
                result.FragmentsDetected = detected;

                if (!result.MeetsMinFragments(parameters.MinFragmentsRequired))
                    continue;

                // ── Get fragment results span for this precursor ────────────
                ReadOnlySpan<FragmentResult> fragmentResults =
                    extractionResult.Results.AsSpan(group.QueryOffset, group.QueryCount);
                ReadOnlySpan<float> libIntensities = input.FragmentIntensities.AsSpan();

                // ── Always compute apex score ───────────────────────────────
                DiaTemporalScore apexScore = apexScorer.ScorePrecursor(
                    libIntensities, group.QueryCount, fragmentResults, rtBuffer, intensityBuffer);
                result.ApexDotProductScore = apexScore.DotProductScore;
                result.ApexTimeIndex = apexScore.ApexTimeIndex;

                // ── Always compute temporal cosine score ─────────────────────
                DiaTemporalScore temporalScore = temporalScorer.ScorePrecursor(
                    libIntensities, group.QueryCount, fragmentResults, rtBuffer, intensityBuffer);
                result.TemporalCosineScore = temporalScore.DotProductScore;
                result.TimePointsUsed = temporalScore.TimePointsUsed;

                // ── Set primary DotProductScore from the requested strategy ──
                if (parameters.ScoringStrategy == ScoringStrategy.ConsensusApex)
                {
                    result.DotProductScore = apexScore.DotProductScore;
                    result.RawCosine = apexScore.RawCosine;
                }
                else if (parameters.ScoringStrategy == ScoringStrategy.TemporalCosine)
                {
                    result.DotProductScore = temporalScore.DotProductScore;
                    result.RawCosine = temporalScore.RawCosine;
                }
                else if (parameters.ScoringStrategy == ScoringStrategy.WeightedTemporalCosineWithTransform)
                {
                    DiaTemporalScore weightedScore = primaryScorer.ScorePrecursor(
                        libIntensities, group.QueryCount, fragmentResults, rtBuffer, intensityBuffer);
                    result.DotProductScore = weightedScore.DotProductScore;
                    result.RawCosine = weightedScore.RawCosine;
                }
                else // Summed
                {
                    DiaTemporalScore summedScore = primaryScorer.ScorePrecursor(
                        libIntensities, group.QueryCount, fragmentResults, rtBuffer, intensityBuffer);
                    result.DotProductScore = summedScore.DotProductScore;
                    result.RawCosine = summedScore.RawCosine;
                }

                // Compute spectral angle from the raw cosine (pre-transform)
                if (!float.IsNaN(result.RawCosine))
                {
                    float clampedCosine = Math.Clamp(result.RawCosine, 0f, 1f);
                    result.SpectralAngleScore = 1.0f - (2.0f / MathF.PI) * MathF.Acos(clampedCosine);
                }

                // Apply score threshold filter
                if (!float.IsNaN(result.DotProductScore) && result.DotProductScore < parameters.MinScoreThreshold)
                    continue;

                results.Add(result);
            }

            return results;
        }
    }
}