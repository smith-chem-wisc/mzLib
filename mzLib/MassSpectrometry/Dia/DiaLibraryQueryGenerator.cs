// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;
using System.Collections.Generic;
using MassSpectrometry.Dia.Calibration;

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

                // BUG FIX (2026): Previously only checked p.RetentionTime.HasValue.
                // For iRT libraries, RetentionTime is null and IrtValue is set instead.
                // The old code silently fell through to the full-run fallback for every
                // iRT-annotated precursor, producing 120-minute windows instead of ±5 min.
                // Fix: resolve IrtValue ?? RetentionTime (IrtValue takes priority when both
                // are set, because calibration is fitted on iRT, not on library RT).
                // If this was intentionally using only RetentionTime (e.g. for a library
                // type that stores calibrated RT in RetentionTime and wants IrtValue ignored),
                // revert this line and document the library type here.
                double? resolvedRt = p.IrtValue ?? p.RetentionTime;

                if (resolvedRt.HasValue)
                {
                    float rt = (float)resolvedRt.Value;
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

        // ════════════════════════════════════════════════════════════════
        //  Phase 15, Prompt 4: RT-Adaptive Window Widths
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Generates extraction queries using a calibration model with RT-adaptive window widths.
        /// Each precursor's RT window is sized based on the local calibration σ at that precursor's
        /// predicted RT, rather than using a single global σ.
        ///
        /// For linear models: GetLocalSigma returns the global σ, so this behaves identically
        /// to GenerateCalibrated with a uniform window.
        /// For piecewise/LOWESS models: GetLocalSigma varies by RT region, producing tighter
        /// windows where calibration is better (dense anchor regions) and wider windows where
        /// calibration is uncertain (gradient edges with fewer anchors).
        ///
        /// Applies a minimum window half-width of minWindowHalfWidthMinutes to prevent
        /// excessively tight windows that could miss valid peptides due to local calibration
        /// artifacts or small-number statistics in per-segment σ.
        ///
        /// Uses the same two-pass (count, then fill) pattern as Generate and GenerateCalibrated.
        /// </summary>
        /// <param name="precursors">Library precursor inputs.</param>
        /// <param name="scanIndex">The DIA scan index.</param>
        /// <param name="parameters">Search parameters (PpmTolerance, RtToleranceMinutes fallback).</param>
        /// <param name="calibration">
        /// The sealed RtCalibrationModel used for RT prediction. For local σ, the method
        /// uses the detailedModel parameter.
        /// </param>
        /// <param name="detailedModel">
        /// The IRtCalibrationModel that provides GetLocalSigma(). If null, falls back to
        /// a uniform window using calibration.SigmaMinutes.
        /// </param>
        /// <param name="sigmaMultiplier">
        /// Window = predicted ± max(sigmaMultiplier × localSigma, minWindowHalfWidthMinutes).
        /// Default 3.0.
        /// </param>
        /// <param name="minWindowHalfWidthMinutes">
        /// Minimum window half-width in minutes. Prevents excessively narrow windows
        /// in regions with very low local σ. Default 0.3 min.
        /// </param>
        /// <returns>GenerationResult with per-precursor RT-adaptive windows.</returns>
        public static GenerationResult GenerateCalibratedAdaptive(
    IList<LibraryPrecursorInput> precursors,
    DiaScanIndex scanIndex,
    DiaSearchParameters parameters,
    RtCalibrationModel calibration,
    IRtCalibrationModel detailedModel = null,
    double sigmaMultiplier = 3.0,
    double minWindowHalfWidthMinutes = 0.3)
        {
            if (precursors == null) throw new ArgumentNullException(nameof(precursors));
            if (scanIndex == null) throw new ArgumentNullException(nameof(scanIndex));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (calibration == null) throw new ArgumentNullException(nameof(calibration));

            float ppmTolerance = parameters.PpmTolerance;

            float globalRtMin = scanIndex.GetGlobalRtMin();
            float globalRtMax = scanIndex.GetGlobalRtMax();

            double uniformHalfWidth = Math.Max(
                calibration.SigmaMinutes * sigmaMultiplier,
                minWindowHalfWidthMinutes);

            // Find first decoy index for paired-target window placement
            int firstDecoyIdx = -1;
            for (int i = 0; i < precursors.Count; i++)
            {
                if (precursors[i].IsDecoy) { firstDecoyIdx = i; break; }
            }

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

            // Second pass: fill with RT-adaptive windows
            int queryIndex = 0;
            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];
                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0 || p.FragmentCount == 0)
                    continue;

                float rtMin, rtMax;

                // For window placement, decoys use their paired target's iRT so they
                // compete in the same local signal environment as their target partner.
                // Assembly still uses the decoy's own irtValue for LibraryRetentionTime,
                // so RtDeviationMinutes is computed correctly relative to the decoy's
                // own predicted RT.
                double? extractionIrt;
                if (p.IsDecoy && firstDecoyIdx >= 0)
                {
                    int pairedTargetIdx = i - firstDecoyIdx;
                    extractionIrt = pairedTargetIdx >= 0 && pairedTargetIdx < firstDecoyIdx
                        ? precursors[pairedTargetIdx].IrtValue ?? precursors[pairedTargetIdx].RetentionTime
                        : p.IrtValue ?? p.RetentionTime;
                }
                else
                {
                    extractionIrt = p.IrtValue ?? p.RetentionTime;
                }

                if (extractionIrt.HasValue)
                {
                    double libRt = extractionIrt.Value;

                    float predictedRt = (float)calibration.ToMinutes(libRt);

                    double localSigma;
                    if (detailedModel != null)
                    {
                        localSigma = detailedModel.GetLocalSigma(libRt);
                        if (double.IsNaN(localSigma) || localSigma <= 0)
                            localSigma = calibration.SigmaMinutes;
                    }
                    else
                    {
                        localSigma = calibration.SigmaMinutes;
                    }

                    double halfWidth = Math.Max(sigmaMultiplier * localSigma, minWindowHalfWidthMinutes);

                    rtMin = (float)(predictedRt - halfWidth);
                    rtMax = (float)(predictedRt + halfWidth);
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

                // BUG FIX (2026): Previously passed input.RetentionTime as libraryRetentionTime.
                // For iRT libraries, RetentionTime is null → LibraryRetentionTime was null on
                // every result → RecalibrateRtDeviations silently skipped all results.
                // Fix: resolve IrtValue ?? RetentionTime so iRT-annotated precursors propagate
                // their iRT value through to the result for RT calibration downstream.
                // If this was intentionally null for iRT libraries (e.g. to suppress RT
                // deviation computation entirely for that library type), revert and document here.
                var result = new DiaSearchResult(
                    sequence: input.Sequence,
                    chargeState: input.ChargeState,
                    precursorMz: input.PrecursorMz,
                    windowId: group.WindowId,
                    isDecoy: input.IsDecoy,
                    fragmentsQueried: group.QueryCount,
                    libraryRetentionTime: input.IrtValue ?? input.RetentionTime,
                    rtWindowStart: group.RtMin,
                    rtWindowEnd: group.RtMax
                );

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
            DiaSearchParameters parameters,
            DiaScanIndex index = null)
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

                // BUG FIX (2026): Same fix as AssembleResults() above — see comment there.
                // Previously passed input.RetentionTime, which is null for iRT libraries,
                // causing LibraryRetentionTime = null on every result and silently disabling
                // RT recalibration for the entire iRT library path.
                var result = new DiaSearchResult(
                    sequence: input.Sequence,
                    chargeState: input.ChargeState,
                    precursorMz: input.PrecursorMz,
                    windowId: group.WindowId,
                    isDecoy: input.IsDecoy,
                    fragmentsQueried: group.QueryCount,
                    libraryRetentionTime: input.IrtValue ?? input.RetentionTime,
                    rtWindowStart: group.RtMin,
                    rtWindowEnd: group.RtMax
                );

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
                    result.ApexScore = float.NaN;
                    result.TemporalScore = float.NaN;
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
                        matrix, refRts, libIntensities, fragmentCount, timePointCount, minCandidateFraction: 0.01f);

                    result.DetectedPeakGroup = peakGroup;

                    // ── Full-window scoring (backward compatibility) ────────
                    // Apex: find time point with maximum total signal
                    // NEW: cosine-weighted apex — picks the time point with best spectral match
                    // weighted by total intensity, so it favors the true peptide peak over
                    // high-intensity interfering signal at the wrong RT.
                    int fullApexIdx = 0;
                    float fullApexScore = -1f;
                    for (int t = 0; t < timePointCount; t++)
                    {
                        int rowOffset = t * fragmentCount;

                        // Count active fragments and total intensity at this time point
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

                        if (activeCount < parameters.MinFragmentsRequired || totalIntensity <= 0f)
                            continue;

                        float cosine = CosineActiveFragments(libIntensities, matrix, rowOffset, fragmentCount);
                        if (float.IsNaN(cosine)) continue;

                        // Score = cosine × log(intensity) — rewards spectral match scaled by signal strength
                        float score = cosine * MathF.Log(totalIntensity + 1f);
                        if (score > fullApexScore)
                        {
                            fullApexScore = score;
                            fullApexIdx = t;
                        }
                    }

                    result.ApexScore = CosineActiveFragments(
                        libIntensities, matrix, fullApexIdx * fragmentCount, fragmentCount);

                    // Full-window temporal cosine
                    ComputeTemporalCosineOnRange(
                        libIntensities, matrix, fragmentCount, timePointCount,
                        0, timePointCount - 1, minActiveFragments: 3,
                        out float fullTemporalScore, out int fullTimePointsUsed);
                    result.TemporalScore = fullTemporalScore;
                    result.TimePointsUsed = fullTimePointsUsed;

                    // Set primary DotProductScore
                    result.DotProductScore = result.TemporalScore;

                    // ── Full-window fragment correlations ────────────────────
                    ComputeFragmentCorrelationsOnRange(
                        matrix, fragmentCount, timePointCount, 0, timePointCount - 1,
                        out float fullMeanCorr, out float fullMinCorr);
                    result.MeanFragCorr = fullMeanCorr;
                    result.MinFragCorr = fullMinCorr;

                    // ── Peak-boundary-restricted scoring ────────────────────
                    if (peakGroup.IsValid)
                    {
                        int peakLeft = peakGroup.LeftIndex;
                        int peakRight = peakGroup.RightIndex;

                        // Peak apex score
                        result.PeakApexScore = CosineActiveFragments(
                            libIntensities, matrix, peakGroup.ApexIndex * fragmentCount, fragmentCount);

                        // Peak-restricted fragment correlations
                        ComputeFragmentCorrelationsOnRange(
                            matrix, fragmentCount, timePointCount, peakLeft, peakRight,
                            out float peakMeanCorr, out float peakMinCorr);
                        result.PeakMeanFragCorr = peakMeanCorr;
                        result.PeakMinFragCorr = peakMinCorr;
                    }
                    else
                    {
                        // No valid peak found — peak features fall back to full-window
                        result.PeakApexScore = result.ApexScore;
                        result.PeakMeanFragCorr = result.MeanFragCorr;
                        result.PeakMinFragCorr = result.MinFragCorr;
                    }
                    // ══════════════════════════════════════════════════════════
                    //  Phase 13: Advanced discriminative features
                    //  Computed from the existing matrix — no extra file reads.
                    // ══════════════════════════════════════════════════════════
                    int apexLocalIdx = peakGroup.IsValid ? peakGroup.ApexIndex : fullApexIdx;
                    result.ObservedApexRt = refRts[apexLocalIdx];

                    // ── Count detected fragments (≥1 nonzero across all scans) ──
                    int detectedFragmentCount = 0;
                    for (int f = 0; f < fragmentCount; f++)
                    {
                        for (int s = 0; s < timePointCount; s++)
                        {
                            if (matrix[s * fragmentCount + f] > 0)
                            {
                                detectedFragmentCount++;
                                break;
                            }
                        }
                    }

                    // 1. Mass accuracy at apex (needs global scan index for m/z lookup)
                    //    Resolve apex RT → nearest scan in this precursor's DIA window
                    //    Uses binary search (O(log N)) since scans are sorted by RT.
                    if (index != null)
                    {
                        float apexRt = refRts[apexLocalIdx];
                        int apexScanGlobalIndex = -1;
                        if (index.TryGetScanRangeForWindow(group.WindowId, out int winStart, out int winCount))
                        {
                            apexScanGlobalIndex = FindClosestScanByRt(index, winStart, winCount, apexRt);
                        }
                        if (apexScanGlobalIndex >= 0)
                        {
                            DiaMassAccuracyHelper.ComputeMassAccuracyAtApex(
                                result, input, index, apexScanGlobalIndex, parameters.PpmTolerance);
                        }
                    }

                    // 2. Signal ratio deviation (needs apex intensities + library intensities)
                    //    Uses ArrayPool to avoid per-precursor heap allocation.
                    float[] apexIntensities = ArrayPool<float>.Shared.Rent(fragmentCount);
                    try
                    {
                        Array.Clear(apexIntensities, 0, fragmentCount);
                        for (int f = 0; f < fragmentCount; f++)
                            apexIntensities[f] = matrix[apexLocalIdx * fragmentCount + f];

                        DiaSignalRatioHelper.ComputeSignalRatioFeatures(
                            apexIntensities, input.FragmentIntensities, fragmentCount, result);
                    }
                    finally
                    {
                        ArrayPool<float>.Shared.Return(apexIntensities);
                    }

                    // 3. Best-fragment reference curve (needs intensity matrix)
                    DiaBestFragmentHelper.ComputeBestFragmentFeatures(
                        matrix, fragmentCount, timePointCount, detectedFragmentCount, result,
                        input.FragmentIntensities.AsSpan(), apexLocalIdx * fragmentCount);

                    // 4. Peak shape ratio features (boundary signal + apex prominence)
                    if (peakGroup.IsValid)
                    {
                        DiaPeakShapeHelper.ComputePeakShapeRatios(
                            matrix, fragmentCount,
                            peakGroup.LeftIndex, peakGroup.RightIndex, peakGroup.ApexIndex,
                            result);
                    }

                    // 5. Smoothed fragment correlations (works on internal pooled copy)
                    DiaSmoothedFeatureHelper.ComputeSmoothedCorrelationFeatures(
                        matrix, fragmentCount, timePointCount, result);

                    // 6. Signal-to-noise (uses apex local index, not global)
                    DiaSmoothedFeatureHelper.ComputeSignalToNoise(
                        matrix, apexLocalIdx, fragmentCount, timePointCount, result);

                    // ── Set peak shape metadata from peak detection ──────────
                    if (peakGroup.IsValid)
                    {
                        result.CandidateCount = peakGroup.CandidateCount;
                        result.PeakWidth = (peakGroup.RightIndex < timePointCount && peakGroup.LeftIndex >= 0)
                            ? refRts[peakGroup.RightIndex] - refRts[peakGroup.LeftIndex]
                            : 0f;
                    }
                    // ── Spectral angle from temporal cosine ─────────────────
                    if (!float.IsNaN(result.TemporalScore))
                    {
                        float clampedCosine = Math.Clamp(result.TemporalScore, 0f, 1f);
                        result.SpectralAngleScore = 1.0f - (2.0f / MathF.PI) * MathF.Acos(clampedCosine);
                    }

                    // 7. Library coverage fraction (intensity-weighted detected fraction)
                    // Numerator: sum of library intensities for fragments with XIC data (XicPointCounts > 0).
                    // Denominator: sum of all library intensities.
                    // Targets: most high-intensity fragments detected → score near 1.0.
                    // Decoys: random m/z values → sparse detection → score near 0.
                    // Computed here (not in DiaFeatureExtractor) because input.FragmentIntensities
                    // is only available during assembly.
                    {
                        float libTotalWeight = 0f;
                        float libDetectedWeight = 0f;
                        for (int f = 0; f < fragmentCount; f++)
                        {
                            float libW = input.FragmentIntensities[f];
                            if (libW <= 0f) continue;
                            libTotalWeight += libW;
                            if (result.XicPointCounts[f] > 0)
                                libDetectedWeight += libW;
                        }
                        result.LibraryCoverageFraction = libTotalWeight > 0f
                            ? libDetectedWeight / libTotalWeight
                            : float.NaN;
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
        /// Binary search for the scan with RT closest to targetRt within [winStart, winStart+winCount).
        /// Scans must be sorted by RT within this range (guaranteed by DiaScanIndex).
        /// Returns the global scan index of the closest scan. O(log N) instead of O(N).
        /// </summary>
        private static int FindClosestScanByRt(DiaScanIndex index, int winStart, int winCount, float targetRt)
        {
            int lo = winStart;
            int hi = winStart + winCount - 1;

            while (lo < hi)
            {
                int mid = lo + (hi - lo) / 2;
                float midRt = index.GetScanRt(mid);
                if (midRt < targetRt)
                    lo = mid + 1;
                else
                    hi = mid;
            }

            // lo is now the first scan with RT >= targetRt (or the last scan if all are <).
            // The closest scan is either lo or lo-1; compare both.
            int best = lo;
            float bestDiff = MathF.Abs(index.GetScanRt(lo) - targetRt);

            if (lo > winStart)
            {
                float prevDiff = MathF.Abs(index.GetScanRt(lo - 1) - targetRt);
                if (prevDiff < bestDiff)
                {
                    best = lo - 1;
                }
            }

            return best;
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

        // ════════════════════════════════════════════════════════════════
        //  Chimeric Score
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Computes a chimeric interference score for each result.
        ///
        /// For each result R, the score is the fraction of R's total extracted fragment
        /// intensity that comes from fragments NOT shared (within ppm) with any other
        /// precursor in the same DIA isolation window:
        ///
        ///   ChimericScore = uncontested_intensity / total_intensity
        ///
        /// A fragment of R is "contested" if any fragment of any other precursor Q in the
        /// same window falls within <paramref name="ppmTolerance"/> ppm of that fragment's m/z.
        ///
        /// Interpretation:
        ///   1.0 — no fragment m/z overlaps with any co-eluting precursor (clean precursor)
        ///   0.0 — every fragment is shared with at least one co-eluting precursor
        ///   x   — fraction x of the total signal comes from uncontested fragments
        ///
        /// When total extracted intensity is zero, ChimericScore is set to 1.0 (no evidence
        /// of contamination).
        ///
        /// Must be called after assembly (ExtractedIntensities populated) and before FDR.
        /// </summary>
        /// <param name="precursors">Library precursor inputs carrying fragment m/z values.</param>
        /// <param name="results">Assembled results. ChimericScore is set in place.</param>
        /// <param name="ppmTolerance">
        /// PPM tolerance for fragment m/z overlap detection. A fragment is contested if
        /// any other precursor in the same window has a fragment within this tolerance.
        /// </param>
        public static void ComputeChimericScores(
            IList<LibraryPrecursorInput> precursors,
            List<DiaSearchResult> results,
            float ppmTolerance)
        {
            if (results == null || results.Count == 0)
                return;

            // Build a lookup: (Sequence, ChargeState, IsDecoy) → precursor index.
            var precursorLookup = new Dictionary<(string, int, bool), int>(results.Count);
            for (int i = 0; i < precursors.Count; i++)
            {
                var key = (precursors[i].Sequence, precursors[i].ChargeState, precursors[i].IsDecoy);
                precursorLookup.TryAdd(key, i);
            }

            // Group result indices by window ID.
            var windowGroups = new Dictionary<int, List<int>>();
            for (int i = 0; i < results.Count; i++)
            {
                int wid = results[i].WindowId;
                if (!windowGroups.TryGetValue(wid, out var list))
                {
                    list = new List<int>();
                    windowGroups[wid] = list;
                }
                list.Add(i);
            }

            foreach (var (_, indices) in windowGroups)
            {
                int n = indices.Count;

                // Cache per-result data needed in the hot loops to avoid repeated
                // property access through the result list.
                var precursorIndices = new int[n];
                var rtStart = new float[n];
                var rtEnd = new float[n];

                for (int a = 0; a < n; a++)
                {
                    var r = results[indices[a]];
                    precursorIndices[a] = precursorLookup.TryGetValue(
                        (r.Sequence, r.ChargeState, r.IsDecoy), out int pi) ? pi : -1;
                    rtStart[a] = r.RtWindowStart;
                    rtEnd[a] = r.RtWindowEnd;
                }

                for (int a = 0; a < n; a++)
                {
                    var ra = results[indices[a]];
                    int piA = precursorIndices[a];
                    if (piA < 0)
                    {
                        ra.ChimericScore = 1.0f;
                        continue;
                    }

                    var precA = precursors[piA];
                    float totalIntensity = 0f;
                    float uncontestedIntensity = 0f;

                    for (int f = 0; f < precA.FragmentCount; f++)
                    {
                        float intensity = ra.ExtractedIntensities[f];
                        totalIntensity += intensity;

                        float mzF = precA.FragmentMzs[f];
                        bool contested = false;

                        for (int b = 0; b < n && !contested; b++)
                        {
                            if (b == a) continue;
                            int piB = precursorIndices[b];
                            if (piB < 0) continue;

                            // Only consider co-eluters: precursors whose RT extraction
                            // window overlaps with result A's window. This prevents the
                            // entire window's ~2800 precursors (spanning the full run)
                            // from contesting every fragment, which would drive all scores
                            // to zero. Two windows overlap when neither ends before the
                            // other starts.
                            if (rtEnd[b] <= rtStart[a] || rtStart[b] >= rtEnd[a])
                                continue;

                            var precB = precursors[piB];
                            for (int g = 0; g < precB.FragmentCount && !contested; g++)
                            {
                                float mzG = precB.FragmentMzs[g];
                                float ppmDiff = MathF.Abs(mzF - mzG) / mzF * 1e6f;
                                if (ppmDiff <= ppmTolerance)
                                    contested = true;
                            }
                        }

                        if (!contested)
                            uncontestedIntensity += intensity;
                    }

                    ra.ChimericScore = totalIntensity > 0f
                        ? uncontestedIntensity / totalIntensity
                        : 1.0f;
                }
            }
        }

    }
}