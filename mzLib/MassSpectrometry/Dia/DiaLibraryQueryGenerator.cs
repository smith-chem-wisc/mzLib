// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;
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
        /// Temporal-scoring result assembly with Phase 12 peak group detection.
        /// 
        /// For each precursor:
        ///   1. Builds the time × fragment intensity matrix ONCE (shared across all scoring)
        ///   2. Runs peak group detection to find the chromatographic peak boundaries
        ///   3. Computes temporal cosine and apex scores WITHIN peak boundaries
        ///   4. Computes fragment correlations WITHIN peak boundaries
        ///   5. Also computes full-window scores for backward compatibility
        /// 
        /// This replaces the Phase 9/11 approach where the matrix was built twice
        /// (once for temporal scoring, once for correlations) and no peak boundaries
        /// were used. The refactoring both improves scoring quality (by restricting to
        /// the actual elution peak) and improves performance (single matrix build).
        /// </summary>
        public static List<DiaSearchResult> AssembleResultsWithTemporalScoring(
            IList<LibraryPrecursorInput> precursors,
            GenerationResult generationResult,
            ExtractionResult extractionResult,
            DiaSearchParameters parameters)
        {
            if (precursors == null) throw new ArgumentNullException(nameof(precursors));
            if (extractionResult == null) throw new ArgumentNullException(nameof(extractionResult));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));

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

                // ── Populate summed intensities and XIC point counts ────────
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

                // ── Get fragment results and library intensities ─────────────
                ReadOnlySpan<FragmentResult> fragmentResults =
                    extractionResult.Results.AsSpan(group.QueryOffset, group.QueryCount);
                ReadOnlySpan<float> libIntensities = input.FragmentIntensities.AsSpan();

                // ══════════════════════════════════════════════════════════════
                //  Phase 12: Build matrix ONCE, detect peak, score within peak
                // ══════════════════════════════════════════════════════════════

                // Find reference fragment (most data points) for the RT grid
                int refFragIdx = -1;
                int maxPts = 0;
                for (int f = 0; f < group.QueryCount; f++)
                {
                    int pts = extractionResult.Results[group.QueryOffset + f].DataPointCount;
                    if (pts > maxPts) { maxPts = pts; refFragIdx = f; }
                }

                if (maxPts == 0)
                {
                    // No data at all — set NaN scores and add
                    result.DotProductScore = float.NaN;
                    result.ApexDotProductScore = float.NaN;
                    result.TemporalCosineScore = float.NaN;
                    results.Add(result);
                    continue;
                }

                int refQi = group.QueryOffset + refFragIdx;
                var refResult = extractionResult.Results[refQi];
                ReadOnlySpan<float> refRts = rtBuffer.Slice(refResult.RtBufferOffset, refResult.DataPointCount);
                int timePointCount = refRts.Length;

                int fragmentCount = group.QueryCount;
                int matrixSize = timePointCount * fragmentCount;

                float[] matrixRented = ArrayPool<float>.Shared.Rent(matrixSize);
                Span<float> matrix = matrixRented.AsSpan(0, matrixSize);
                matrix.Clear();

                try
                {
                    // ── Build the time × fragment matrix ────────────────────
                    const float rtTol = 0.01f;
                    for (int f = 0; f < fragmentCount; f++)
                    {
                        var fr = fragmentResults[f];
                        if (fr.DataPointCount == 0) continue;

                        ReadOnlySpan<float> fragRts = rtBuffer.Slice(fr.RtBufferOffset, fr.DataPointCount);
                        ReadOnlySpan<float> fragInts = intensityBuffer.Slice(fr.IntensityBufferOffset, fr.DataPointCount);

                        // Two-pointer alignment
                        int fragPtr = 0;
                        for (int t = 0; t < timePointCount && fragPtr < fragRts.Length; t++)
                        {
                            float refRt = refRts[t];
                            while (fragPtr < fragRts.Length && fragRts[fragPtr] < refRt - rtTol)
                                fragPtr++;
                            if (fragPtr < fragRts.Length && MathF.Abs(fragRts[fragPtr] - refRt) <= rtTol)
                            {
                                matrix[t * fragmentCount + f] = fragInts[fragPtr];
                                fragPtr++;
                            }
                        }
                    }

                    // ── Peak Group Detection ────────────────────────────────
                    PeakGroup peakGroup = DiaPeakGroupDetector.Detect(
                        matrix, refRts, libIntensities, fragmentCount, timePointCount);
                    result.DetectedPeakGroup = peakGroup;

                    // ── Full-window scoring (backward compatibility) ────────
                    // Apex: find time point with maximum total signal
                    int fullApexIdx = 0;
                    float fullApexSignal = 0f;
                    for (int t = 0; t < timePointCount; t++)
                    {
                        float total = 0f;
                        int rowOffset = t * fragmentCount;
                        for (int f = 0; f < fragmentCount; f++)
                            total += matrix[rowOffset + f];
                        if (total > fullApexSignal) { fullApexSignal = total; fullApexIdx = t; }
                    }

                    result.ApexDotProductScore = CosineActiveFragments(
                        libIntensities, matrix, fullApexIdx * fragmentCount, fragmentCount);
                    result.ApexTimeIndex = fullApexIdx;

                    // Observed apex RT from full window
                    if (fullApexIdx < timePointCount)
                        result.ObservedApexRt = refRts[fullApexIdx];

                    // Full-window temporal cosine
                    ComputeTemporalCosineOnRange(
                        libIntensities, matrix, fragmentCount, timePointCount,
                        0, timePointCount - 1, minActiveFragments: 3,
                        out float fullTemporalScore, out int fullTimePointsUsed);
                    result.TemporalCosineScore = fullTemporalScore;
                    result.TimePointsUsed = fullTimePointsUsed;

                    // Set primary DotProductScore
                    result.DotProductScore = result.TemporalCosineScore;
                    result.RawCosine = result.TemporalCosineScore;

                    // ── Full-window fragment correlations ────────────────────
                    ComputeFragmentCorrelationsOnRange(
                        matrix, fragmentCount, timePointCount, 0, timePointCount - 1,
                        out float fullMeanCorr, out float fullMinCorr);
                    result.MeanFragmentCorrelation = fullMeanCorr;
                    result.MinFragmentCorrelation = fullMinCorr;

                    // ── Peak-boundary-restricted scoring ────────────────────
                    if (peakGroup.IsValid)
                    {
                        int peakLeft = peakGroup.LeftIndex;
                        int peakRight = peakGroup.RightIndex;

                        // Peak apex score
                        result.PeakApexScore = CosineActiveFragments(
                            libIntensities, matrix, peakGroup.ApexIndex * fragmentCount, fragmentCount);

                        // Override observed apex RT with peak-detected apex
                        result.ObservedApexRt = peakGroup.ApexRt;

                        // Peak-restricted temporal cosine
                        ComputeTemporalCosineOnRange(
                            libIntensities, matrix, fragmentCount, timePointCount,
                            peakLeft, peakRight, minActiveFragments: 3,
                            out float peakTemporalScore, out int _);
                        result.PeakTemporalScore = peakTemporalScore;

                        // Peak-restricted fragment correlations
                        ComputeFragmentCorrelationsOnRange(
                            matrix, fragmentCount, timePointCount, peakLeft, peakRight,
                            out float peakMeanCorr, out float peakMinCorr);
                        result.PeakMeanFragCorrelation = peakMeanCorr;
                        result.PeakMinFragCorrelation = peakMinCorr;
                    }
                    else
                    {
                        // No valid peak found — peak features fall back to full-window
                        result.PeakApexScore = result.ApexDotProductScore;
                        result.PeakTemporalScore = result.TemporalCosineScore;
                        result.PeakMeanFragCorrelation = result.MeanFragmentCorrelation;
                        result.PeakMinFragCorrelation = result.MinFragmentCorrelation;
                    }

                    // ── Spectral angle from raw cosine ──────────────────────
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
                finally
                {
                    ArrayPool<float>.Shared.Return(matrixRented);
                }
            }

            return results;
        }

        /// <summary>
        /// Computes temporal cosine score over a restricted range of time points.
        /// This is the core scoring operation shared by full-window and peak-restricted paths.
        /// 
        /// At each time point in [rangeStart, rangeEnd]:
        ///   - Count fragments with nonzero intensity
        ///   - If enough active fragments, compute cosine similarity against library
        ///   - Weight by total intensity at that time point
        ///   - Return the intensity-weighted average cosine
        /// </summary>
        private static void ComputeTemporalCosineOnRange(
            ReadOnlySpan<float> libraryIntensities,
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            int rangeStart,
            int rangeEnd,
            int minActiveFragments,
            out float temporalScore,
            out int validTimePoints)
        {
            float weightedCosineSum = 0f;
            float weightSum = 0f;
            validTimePoints = 0;

            rangeStart = Math.Max(0, rangeStart);
            rangeEnd = Math.Min(timePointCount - 1, rangeEnd);

            for (int t = rangeStart; t <= rangeEnd; t++)
            {
                int rowOffset = t * fragmentCount;
                int activeCount = 0;
                float totalIntensity = 0f;
                for (int f = 0; f < fragmentCount; f++)
                {
                    if (matrix[rowOffset + f] > 0f)
                    {
                        activeCount++;
                        totalIntensity += matrix[rowOffset + f];
                    }
                }

                if (activeCount < minActiveFragments || totalIntensity <= 0f)
                    continue;

                float cos_t = CosineActiveFragments(libraryIntensities, matrix, rowOffset, fragmentCount);
                if (float.IsNaN(cos_t))
                    continue;

                float weight = 1.0f; // uniform weighting
                weightedCosineSum += weight * cos_t;
                weightSum += weight;
                validTimePoints++;
            }

            temporalScore = weightSum > 0f ? weightedCosineSum / weightSum : float.NaN;
        }

        /// <summary>
        /// Computes pairwise Pearson correlation between fragment XICs over a restricted
        /// time range. Used by both full-window and peak-restricted scoring paths.
        /// </summary>
        private static void ComputeFragmentCorrelationsOnRange(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            int rangeStart,
            int rangeEnd,
            out float meanCorrelation,
            out float minCorrelation)
        {
            meanCorrelation = float.NaN;
            minCorrelation = float.NaN;

            rangeStart = Math.Max(0, rangeStart);
            rangeEnd = Math.Min(timePointCount - 1, rangeEnd);
            int rangeLength = rangeEnd - rangeStart + 1;
            if (rangeLength < 3) return;

            // Find fragments with >= 3 nonzero time points within range
            Span<int> detectedFrags = stackalloc int[Math.Min(fragmentCount, 64)];
            int nDetected = 0;

            for (int f = 0; f < fragmentCount && nDetected < detectedFrags.Length; f++)
            {
                int nonzero = 0;
                for (int t = rangeStart; t <= rangeEnd; t++)
                {
                    if (matrix[t * fragmentCount + f] > 0f) nonzero++;
                }
                if (nonzero >= 3) detectedFrags[nDetected++] = f;
            }

            if (nDetected < 2) return;

            float sumCorr = 0f;
            float minCorr = float.MaxValue;
            int nPairs = 0;

            for (int a = 0; a < nDetected; a++)
            {
                for (int b = a + 1; b < nDetected; b++)
                {
                    float r = PearsonCorrelationOnRange(
                        matrix, detectedFrags[a], detectedFrags[b],
                        fragmentCount, rangeStart, rangeEnd);

                    if (!float.IsNaN(r))
                    {
                        sumCorr += r;
                        if (r < minCorr) minCorr = r;
                        nPairs++;
                    }
                }
            }

            if (nPairs > 0)
            {
                meanCorrelation = sumCorr / nPairs;
                minCorrelation = minCorr;
            }
        }

        /// <summary>
        /// Pearson correlation between two fragments over a restricted time range.
        /// Uses only time points where BOTH fragments have nonzero intensity.
        /// </summary>
        private static float PearsonCorrelationOnRange(
            ReadOnlySpan<float> matrix, int fragA, int fragB,
            int fragmentCount, int rangeStart, int rangeEnd)
        {
            float sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0;
            int n = 0;

            for (int t = rangeStart; t <= rangeEnd; t++)
            {
                float a = matrix[t * fragmentCount + fragA];
                float b = matrix[t * fragmentCount + fragB];
                if (a <= 0f || b <= 0f) continue;

                sumA += a;
                sumB += b;
                sumAB += a * b;
                sumA2 += a * a;
                sumB2 += b * b;
                n++;
            }

            if (n < 3) return float.NaN;

            float denom = (n * sumA2 - sumA * sumA) * (n * sumB2 - sumB * sumB);
            if (denom <= 0f) return float.NaN;

            float r = (n * sumAB - sumA * sumB) / MathF.Sqrt(denom);
            return Math.Clamp(r, -1f, 1f);
        }

        /// <summary>
        /// Cosine between library and observed, using ONLY fragments where observed > 0.
        /// Shared helper used by both full-window and peak-restricted scoring.
        /// </summary>
        private static float CosineActiveFragments(
            ReadOnlySpan<float> library,
            ReadOnlySpan<float> matrix,
            int rowOffset,
            int fragmentCount)
        {
            float dot = 0f, normLib = 0f, normObs = 0f;

            for (int f = 0; f < fragmentCount; f++)
            {
                float obs = matrix[rowOffset + f];
                if (obs <= 0f) continue;

                float lib = f < library.Length ? library[f] : 0f;
                dot += lib * obs;
                normLib += lib * lib;
                normObs += obs * obs;
            }

            if (normLib <= 0f || normObs <= 0f)
                return float.NaN;

            return Math.Clamp(dot / (MathF.Sqrt(normLib) * MathF.Sqrt(normObs)), 0f, 1f);
        }
    }
}