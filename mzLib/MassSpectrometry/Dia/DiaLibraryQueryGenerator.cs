// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry.Dia.Calibration;
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
        // Add this method to DiaLibraryQueryGenerator (MassSpectrometry/Dia/DiaLibraryQueryGenerator.cs)
        // This is the only query generation method used in the bootstrap calibration path.
        // All other Generate* methods (Generate, GenerateCalibrated, GenerateCalibratedAdaptive)
        // are candidates for deletion once this path is validated end-to-end.

        /// <summary>
        /// Generates extraction queries using the LOWESS calibration model from
        /// <see cref="DiaBootstrapCalibrator.RunBidirectional"/>.
        ///
        /// Window per target precursor:
        ///   predictedRt = lowess.ToMinutes(precursor.RetentionTime)
        ///   hw          = Math.Min(lowess.GetLocalSigma(precursor.RetentionTime), localSigmaCap) * 2.0
        ///   window      = [predictedRt - hw, predictedRt + hw]
        ///
        /// Each decoy receives the same window as its paired target. Pairing is by position:
        /// decoy at index i pairs with target at index (i - firstDecoyIndex). This places
        /// the decoy in the same local RT signal environment as its target for fair
        /// target-decoy competition.
        ///
        /// Precursors with no RetentionTime are skipped and counted in SkippedNoFragments.
        ///
        /// NOTE: Skipping no-RT precursors is a deliberate simplification. In the current
        /// dataset, 93.2% of precursors have library RT so the impact is small. Once
        /// end-to-end results are validated, revisit whether a full-run fallback window
        /// is worth adding for the ~7% without RT. For now, skipping keeps the search
        /// clean and avoids flooding with uncalibrated full-run queries. See Bootstrap
        /// RT Calibration Writeup §7 for context.
        /// </summary>
        /// <param name="precursors">Target + decoy library precursors, targets first.</param>
        /// <param name="scanIndex">The DIA scan index.</param>
        /// <param name="ppmTolerance">Fragment m/z tolerance in ppm.</param>
        /// <param name="lowess">LOWESS calibration model from DiaBootstrapCalibrator.RunBidirectional().</param>
        /// <param name="localSigmaCap">
        /// Cap on local σ before doubling to half-width. Default 1.5 min.
        /// Prevents the σ spike at libRT ~105 (observed σ = 3.16 min) from producing
        /// 6+ min windows. At 1.5 min cap, spike region gets hw = 3.0 min which is
        /// conservative but not catastrophic. See Bootstrap RT Calibration Writeup §7.1.
        /// </param>
        public static GenerationResult GenerateFromLowess(
            IReadOnlyList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            float ppmTolerance,
            LowessRtModel lowess,
            double localSigmaCap = 0.5)
        {
            if (precursors == null) throw new ArgumentNullException(nameof(precursors));
            if (scanIndex == null) throw new ArgumentNullException(nameof(scanIndex));
            if (lowess == null) throw new ArgumentNullException(nameof(lowess));

            // ── Find the first decoy index ───────────────────────────────────────
            int firstDecoyIdx = -1;
            for (int i = 0; i < precursors.Count; i++)
                if (precursors[i].IsDecoy) { firstDecoyIdx = i; break; }

            int targetCount = firstDecoyIdx >= 0 ? firstDecoyIdx : precursors.Count;

            // ── Pre-compute target windows ───────────────────────────────────────
            // Stored so decoys can look up their paired target window by index.
            // Decoys use the same window as their positionally-paired target, placing
            // them in the same local RT signal environment for fair T/D competition.
            var targetRtMin = new float[targetCount];
            var targetRtMax = new float[targetCount];
            var targetHasWindow = new bool[targetCount];

            for (int i = 0; i < targetCount; i++)
            {
                var p = precursors[i];
                if (!p.RetentionTime.HasValue)
                    continue; // NOTE: target has no library RT — skipped (see method summary).

                double libRt = p.RetentionTime.Value;
                double predictedRt = lowess.ToMinutes(libRt);
                double localSigma = lowess.GetLocalSigma(libRt);
                double hw = Math.Min(localSigma, localSigmaCap) * 2.0;

                targetRtMin[i] = (float)(predictedRt - hw);
                targetRtMax[i] = (float)(predictedRt + hw);
                targetHasWindow[i] = true;
            }

            // ── First pass: count valid queries ──────────────────────────────────
            int totalQueryCount = 0;
            int skippedNoWindow = 0;
            int skippedNoFragments = 0;

            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];

                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0) { skippedNoWindow++; continue; }

                if (p.FragmentCount == 0) { skippedNoFragments++; continue; }

                bool hasWindow;
                if (!p.IsDecoy)
                {
                    hasWindow = i < targetCount && targetHasWindow[i];
                }
                else
                {
                    int pairedIdx = firstDecoyIdx >= 0 ? i - firstDecoyIdx : -1;
                    hasWindow = pairedIdx >= 0 && pairedIdx < targetCount && targetHasWindow[pairedIdx];
                }

                if (!hasWindow) { skippedNoFragments++; continue; }

                totalQueryCount += p.FragmentCount;
            }

            // ── Allocate ──────────────────────────────────────────────────────────
            var queries = new FragmentQuery[totalQueryCount];
            var groups = new List<PrecursorQueryGroup>(
                precursors.Count - skippedNoWindow - skippedNoFragments);

            // ── Second pass: fill ─────────────────────────────────────────────────
            int queryIndex = 0;
            for (int i = 0; i < precursors.Count; i++)
            {
                var p = precursors[i];

                int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (windowId < 0 || p.FragmentCount == 0) continue;

                float rtMin, rtMax;

                if (!p.IsDecoy)
                {
                    if (i >= targetCount || !targetHasWindow[i]) continue;
                    rtMin = targetRtMin[i];
                    rtMax = targetRtMax[i];
                }
                else
                {
                    int pairedIdx = firstDecoyIdx >= 0 ? i - firstDecoyIdx : -1;
                    if (pairedIdx < 0 || pairedIdx >= targetCount || !targetHasWindow[pairedIdx]) continue;

                    // Narrow decoy extraction to ±2σ around the paired target's predicted RT.
                    // The full target window (±hw) spans ~1-2 min and contains many other
                    // peptides whose peaks the decoy detector would otherwise find.
                    // Constraining to ±2σ forces decoys to extract data only where their
                    // paired target is expected to elute, so a decoy with no matching
                    // signal gets a flat XIC rather than borrowing a real peptide peak.
                    float targetCenter = (targetRtMin[pairedIdx] + targetRtMax[pairedIdx]) * 0.5f;
                    float halfWidth = (float)(2.0 * lowess.SigmaMinutes);
                    rtMin = targetCenter - halfWidth;
                    rtMax = targetCenter + halfWidth;
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
                        queryId: queryIndex);
                    queryIndex++;
                }

                groups.Add(new PrecursorQueryGroup(
                    inputIndex: i,
                    queryOffset: groupOffset,
                    queryCount: p.FragmentCount,
                    windowId: windowId,
                    rtMin: rtMin,
                    rtMax: rtMax));
            }

            //// ── Second pass: fill ─────────────────────────────────────────────────
            //int queryIndex = 0;
            //for (int i = 0; i < precursors.Count; i++)
            //{
            //    var p = precursors[i];

            //    int windowId = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
            //    if (windowId < 0 || p.FragmentCount == 0) continue;

            //    float rtMin, rtMax;

            //    if (!p.IsDecoy)
            //    {
            //        if (i >= targetCount || !targetHasWindow[i]) continue;
            //        rtMin = targetRtMin[i];
            //        rtMax = targetRtMax[i];
            //    }
            //    else
            //    {
            //        int pairedIdx = firstDecoyIdx >= 0 ? i - firstDecoyIdx : -1;
            //        if (pairedIdx < 0 || pairedIdx >= targetCount || !targetHasWindow[pairedIdx]) continue;
            //        rtMin = targetRtMin[pairedIdx];
            //        rtMax = targetRtMax[pairedIdx];
            //    }

            //    int groupOffset = queryIndex;

            //    for (int f = 0; f < p.FragmentCount; f++)
            //    {
            //        queries[queryIndex] = new FragmentQuery(
            //            targetMz: p.FragmentMzs[f],
            //            tolerancePpm: ppmTolerance,
            //            rtMin: rtMin,
            //            rtMax: rtMax,
            //            windowId: windowId,
            //            queryId: queryIndex);
            //        queryIndex++;
            //    }

            //    groups.Add(new PrecursorQueryGroup(
            //        inputIndex: i,
            //        queryOffset: groupOffset,
            //        queryCount: p.FragmentCount,
            //        windowId: windowId,
            //        rtMin: rtMin,
            //        rtMax: rtMax));
            //}

            return new GenerationResult(queries, groups.ToArray(), skippedNoWindow, skippedNoFragments);
        }


        /// <summary>
        /// Generates FragmentQuery[] from library precursor inputs.
        /// Two-pass: first counts for a single allocation, then fills.
        /// </summary>
        public static GenerationResult Generate(
            IReadOnlyList<LibraryPrecursorInput> precursors,
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
        public static List<DiaSearchResult> AssembleResults(
            IReadOnlyList<LibraryPrecursorInput> precursors,
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
                result.PrecursorIndex = group.InputIndex;

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
        /// Assembles DIA search results from extracted XICs, computing all features
        /// needed for the 38-feature classifier vector.
        ///
        /// For each precursor:
        ///   1. Builds the time × fragment intensity matrix from ExtractionResult
        ///   2. Runs peak group detection (SelectBest) using LOWESS-calibrated predicted RT
        ///   3. Computes all scores and features within peak boundaries
        ///   4. Sets RtDeviationMinutes and RtDeviationSquared against calibrated predicted RT
        ///
        /// The LowessRtModel is the only calibration input. RT windows in GenerationResult
        /// are already set correctly by GenerateFromLowess — this method does not adjust them.
        /// </summary>
        public static List<DiaSearchResult> AssembleResultsWithTemporalScoring(
            IReadOnlyList<LibraryPrecursorInput> precursors,
            GenerationResult generationResult,
            ExtractionResult extractionResult,
            DiaSearchParameters parameters,
            DiaScanIndex index = null,
            LowessRtModel lowess = null)   // NEW
        {
            if (precursors == null) throw new ArgumentNullException(nameof(precursors));
            if (extractionResult == null) throw new ArgumentNullException(nameof(extractionResult));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));

            var results = new List<DiaSearchResult>(generationResult.PrecursorGroups.Length);

            ReadOnlySpan<float> rtBuffer = extractionResult.RtBuffer.AsSpan();
            ReadOnlySpan<float> intensityBuffer = extractionResult.IntensityBuffer.AsSpan();

            // ── Find firstDecoyIdx for pairing ──────────────────────────────────────
            int firstDecoyIdx = -1;
            for (int i = 0; i < precursors.Count; i++)
                if (precursors[i].IsDecoy) { firstDecoyIdx = i; break; }

            // ── Target apex record — populated in pass 1, consumed in pass 2 ────────
            // Keyed by precursor input index (= group.InputIndex for targets).
            // Stores the matrix-local apex scan index and the validated peak group
            // so decoys can use the SAME reference point for ALL feature computation.
            // apexIdx is in [0, timePointCount) of the target's own matrix, not a
            // global scan index — it is only valid for locating an RT which is then
            // matched into the decoy's matrix in pass 2.
            var targetApexRt = new Dictionary<int, float>();   // inputIdx → apexRt (minutes)
            var targetPeakLeft = new Dictionary<int, float>();   // inputIdx → leftRt (minutes)
            var targetPeakRight = new Dictionary<int, float>();   // inputIdx → rightRt (minutes)

            // ── Pass 1: targets only — build matrices, score, record apex ───────────
            for (int g = 0; g < generationResult.PrecursorGroups.Length; g++)
            {
                var group = generationResult.PrecursorGroups[g];
                var input = precursors[group.InputIndex];
                if (input.IsDecoy) continue;

                // Build matrix for this target
                int refFragIdx = -1, maxPts = 0;
                for (int f = 0; f < group.QueryCount; f++)
                {
                    int pts = extractionResult.Results[group.QueryOffset + f].DataPointCount;
                    if (pts > maxPts) { maxPts = pts; refFragIdx = f; }
                }
                if (maxPts == 0) continue;

                int refQi = group.QueryOffset + refFragIdx;
                var refRes = extractionResult.Results[refQi];
                ReadOnlySpan<float> refRts = rtBuffer.Slice(refRes.RtBufferOffset, refRes.DataPointCount);
                int timePointCount = refRts.Length;
                int fragmentCount = group.QueryCount;

                float[] matrixRented = ArrayPool<float>.Shared.Rent(timePointCount * fragmentCount);
                Span<float> matrix = matrixRented.AsSpan(0, timePointCount * fragmentCount);
                matrix.Clear();

                try
                {
                    const float rtTol = 0.01f;
                    ReadOnlySpan<FragmentResult> fragmentResults =
                        extractionResult.Results.AsSpan(group.QueryOffset, group.QueryCount);
                    ReadOnlySpan<float> libIntensities = input.FragmentIntensities.AsSpan();

                    for (int f = 0; f < fragmentCount; f++)
                    {
                        var fr = fragmentResults[f];
                        if (fr.DataPointCount == 0) continue;
                        ReadOnlySpan<float> fragRts = rtBuffer.Slice(fr.RtBufferOffset, fr.DataPointCount);
                        ReadOnlySpan<float> fragInts = intensityBuffer.Slice(fr.IntensityBufferOffset, fr.DataPointCount);
                        int fragPtr = 0;
                        for (int t = 0; t < timePointCount && fragPtr < fragRts.Length; t++)
                        {
                            float refRt = refRts[t];
                            while (fragPtr < fragRts.Length && fragRts[fragPtr] < refRt - rtTol) fragPtr++;
                            if (fragPtr < fragRts.Length && MathF.Abs(fragRts[fragPtr] - refRt) <= rtTol)
                            { matrix[t * fragmentCount + f] = fragInts[fragPtr]; fragPtr++; }
                        }
                    }

                    // Detect peak group for target
                    float? predictedRt = (lowess != null && input.RetentionTime.HasValue)
                        ? (float?)lowess.ToMinutes(input.RetentionTime.Value) : null;

                    float peakSearchRtMin, peakSearchRtMax;
                    if (input.IsDecoy && lowess != null)
                    {
                        float center = (group.RtMin + group.RtMax) * 0.5f;
                        float halfWidth = (float)(2.0 * lowess.SigmaMinutes);
                        peakSearchRtMin = center - halfWidth;
                        peakSearchRtMax = center + halfWidth;
                    }
                    else
                    {
                        peakSearchRtMin = group.RtMin;
                        peakSearchRtMax = group.RtMax;
                    }

                    PeakGroup pg = DiaPeakGroupDetector.SelectBest(
                        matrix, refRts, libIntensities, fragmentCount, timePointCount,
                        out float _,
                        predictedRt: predictedRt,
                        rtWindowHalfWidth: (group.RtMax - group.RtMin) * 0.5f,
                        ms1Index: index,
                        precursorMz: (float)input.PrecursorMz,
                        ppmTolerance: parameters.PpmTolerance,
                        rtMin: peakSearchRtMin,
                        rtMax: peakSearchRtMax);

                    // Find apex index within peak boundaries
                    int searchStart = pg.IsValid ? pg.LeftIndex : 0;
                    int searchEnd = pg.IsValid ? pg.RightIndex : timePointCount - 1;
                    int apexIdx = searchStart;
                    float bestScore = -1f;
                    for (int t = searchStart; t <= searchEnd; t++)
                    {
                        int ro = t * fragmentCount;
                        int active = 0; float tot = 0f;
                        for (int f = 0; f < fragmentCount; f++)
                            if (matrix[ro + f] > 0f) { active++; tot += matrix[ro + f]; }
                        if (active < parameters.MinFragmentsRequired || tot <= 0f) continue;
                        float cos = CosineActiveFragments(libIntensities, matrix, ro, fragmentCount);
                        if (float.IsNaN(cos)) continue;
                        float sc = cos * MathF.Log(tot + 1f);
                        if (sc > bestScore) { bestScore = sc; apexIdx = t; }
                    }

                    // Record apex RT and peak RT bounds for paired decoy use
                    float apexRt = refRts[apexIdx];
                    targetApexRt[group.InputIndex] = apexRt;

                    if (pg.IsValid)
                    {
                        targetPeakLeft[group.InputIndex] = pg.LeftRt;
                        targetPeakRight[group.InputIndex] = pg.RightRt;
                    }
                    else
                    {
                        // No valid peak — record apex ± half FWHM as synthetic bounds
                        float halfFwhm = 0.069f * 0.5f;  // ~half of typical chromatographic FWHM
                        targetPeakLeft[group.InputIndex] = apexRt - halfFwhm;
                        targetPeakRight[group.InputIndex] = apexRt + halfFwhm;
                    }
                }
                finally
                {
                    ArrayPool<float>.Shared.Return(matrixRented);
                }
            }

            // ── Pass 2: all precursors — full scoring ────────────────────────────────
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
                    libraryRetentionTime: input.IrtValue ?? input.RetentionTime,
                    rtWindowStart: group.RtMin,
                    rtWindowEnd: group.RtMax
                );
                result.PrecursorIndex = group.InputIndex;

                int detected = 0;
                for (int f = 0; f < group.QueryCount; f++)
                {
                    int qi = group.QueryOffset + f;
                    var fr = extractionResult.Results[qi];
                    result.ExtractedIntensities[f] = fr.TotalIntensity;
                    result.XicPointCounts[f] = fr.DataPointCount;
                    if (fr.DataPointCount > 0) detected++;
                }
                result.FragmentsDetected = detected;

                if (!result.MeetsMinFragments(parameters.MinFragmentsRequired))
                    continue;

                ReadOnlySpan<FragmentResult> fragmentResults =
                    extractionResult.Results.AsSpan(group.QueryOffset, group.QueryCount);
                ReadOnlySpan<float> libIntensities = input.FragmentIntensities.AsSpan();

                int refFragIdx = -1, maxPts = 0;
                for (int f = 0; f < group.QueryCount; f++)
                {
                    int pts = extractionResult.Results[group.QueryOffset + f].DataPointCount;
                    if (pts > maxPts) { maxPts = pts; refFragIdx = f; }
                }

                if (maxPts == 0)
                {
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

                float[] matrixRented = ArrayPool<float>.Shared.Rent(timePointCount * fragmentCount);
                Span<float> matrix = matrixRented.AsSpan(0, timePointCount * fragmentCount);
                matrix.Clear();

                try
                {
                    // ── Build matrix ────────────────────────────────────────────
                    const float rtTol = 0.01f;
                    for (int f = 0; f < fragmentCount; f++)
                    {
                        var fr = fragmentResults[f];
                        if (fr.DataPointCount == 0) continue;
                        ReadOnlySpan<float> fragRts = rtBuffer.Slice(fr.RtBufferOffset, fr.DataPointCount);
                        ReadOnlySpan<float> fragInts = intensityBuffer.Slice(fr.IntensityBufferOffset, fr.DataPointCount);
                        int fragPtr = 0;
                        for (int t = 0; t < timePointCount && fragPtr < fragRts.Length; t++)
                        {
                            float refRt = refRts[t];
                            while (fragPtr < fragRts.Length && fragRts[fragPtr] < refRt - rtTol) fragPtr++;
                            if (fragPtr < fragRts.Length && MathF.Abs(fragRts[fragPtr] - refRt) <= rtTol)
                            { matrix[t * fragmentCount + f] = fragInts[fragPtr]; fragPtr++; }
                        }
                    }

                    // ── Peak group detection ─────────────────────────────────────
                    float? predictedRt = (lowess != null && input.RetentionTime.HasValue)
                        ? (float?)lowess.ToMinutes(input.RetentionTime.Value) : null;

                    float peakSearchRtMin, peakSearchRtMax;
                    if (input.IsDecoy && lowess != null)
                    {
                        float center = (group.RtMin + group.RtMax) * 0.5f;
                        float halfWidth = (float)(2.0 * lowess.SigmaMinutes);
                        peakSearchRtMin = center - halfWidth;
                        peakSearchRtMax = center + halfWidth;
                    }
                    else
                    {
                        peakSearchRtMin = group.RtMin;
                        peakSearchRtMax = group.RtMax;
                    }

                    PeakGroup peakGroup = DiaPeakGroupDetector.SelectBest(
                        matrix, refRts, libIntensities, fragmentCount, timePointCount,
                        out float ms1Factor,
                        predictedRt: predictedRt,
                        rtWindowHalfWidth: (group.RtMax - group.RtMin) * 0.5f,
                        ms1Index: index,
                        precursorMz: (float)input.PrecursorMz,
                        ppmTolerance: parameters.PpmTolerance,
                        rtMin: peakSearchRtMin,
                        rtMax: peakSearchRtMax);

                    result.DetectedPeakGroup = peakGroup;

                    if (peakGroup.IsValid)
                    {
                        result.CoElutionStd = peakGroup.CoElutionStd;
                        result.CandidateScoreGap = peakGroup.SelectionScore - peakGroup.SecondBestScore;
                        result.Ms1ApexConfirmationScore = ms1Factor;
                    }

                    // ── Resolve scoring anchor: apex index and peak bounds ───────
                    // For targets: use the detected peak group apex and bounds.
                    // For decoys:  pin to the paired target's observed apex RT and
                    //              peak bounds — every feature is then evaluated at
                    //              the SAME RT as the target, not the decoy's best scan.
                    int scoringApexIdx;
                    int scoringLeft, scoringRight;

                    if (input.IsDecoy && firstDecoyIdx >= 0)
                    {
                        int pairedTargetInputIdx = group.InputIndex - firstDecoyIdx;

                        if (pairedTargetInputIdx >= 0
                            && targetApexRt.TryGetValue(pairedTargetInputIdx, out float tApexRt)
                            && targetPeakLeft.TryGetValue(pairedTargetInputIdx, out float tLeftRt)
                            && targetPeakRight.TryGetValue(pairedTargetInputIdx, out float tRightRt))
                        {
                            // Find the decoy matrix scan closest to the target's apex RT
                            scoringApexIdx = 0;
                            float minDist = float.MaxValue;
                            for (int t = 0; t < timePointCount; t++)
                            {
                                float dist = MathF.Abs(refRts[t] - tApexRt);
                                if (dist < minDist) { minDist = dist; scoringApexIdx = t; }
                            }

                            // Find decoy matrix scan indices closest to target peak bounds
                            scoringLeft = 0;
                            scoringRight = timePointCount - 1;
                            float minDistL = float.MaxValue, minDistR = float.MaxValue;
                            for (int t = 0; t < timePointCount; t++)
                            {
                                float distL = MathF.Abs(refRts[t] - tLeftRt);
                                float distR = MathF.Abs(refRts[t] - tRightRt);
                                if (distL < minDistL) { minDistL = distL; scoringLeft = t; }
                                if (distR < minDistR) { minDistR = distR; scoringRight = t; }
                            }
                            // Ensure left <= apex <= right
                            if (scoringLeft > scoringApexIdx) scoringLeft = scoringApexIdx;
                            if (scoringRight < scoringApexIdx) scoringRight = scoringApexIdx;
                        }
                        else
                        {
                            // No paired target found — fall back to own peak
                            scoringApexIdx = peakGroup.IsValid ? peakGroup.ApexIndex : timePointCount / 2;
                            scoringLeft = peakGroup.IsValid ? peakGroup.LeftIndex : 0;
                            scoringRight = peakGroup.IsValid ? peakGroup.RightIndex : timePointCount - 1;
                        }
                    }
                    else
                    {
                        // Target: use peak-restricted apex, fall back to full window
                        int searchStart = peakGroup.IsValid ? peakGroup.LeftIndex : 0;
                        int searchEnd = peakGroup.IsValid ? peakGroup.RightIndex : timePointCount - 1;
                        scoringApexIdx = searchStart;
                        float bestScore = -1f;
                        for (int t = searchStart; t <= searchEnd; t++)
                        {
                            int ro = t * fragmentCount;
                            int active = 0; float tot = 0f;
                            for (int f = 0; f < fragmentCount; f++)
                                if (matrix[ro + f] > 0f) { active++; tot += matrix[ro + f]; }
                            if (active < parameters.MinFragmentsRequired || tot <= 0f) continue;
                            float cos = CosineActiveFragments(libIntensities, matrix, ro, fragmentCount);
                            if (float.IsNaN(cos)) continue;
                            float sc = cos * MathF.Log(tot + 1f);
                            if (sc > bestScore) { bestScore = sc; scoringApexIdx = t; }
                        }
                        scoringLeft = peakGroup.IsValid ? peakGroup.LeftIndex : 0;
                        scoringRight = peakGroup.IsValid ? peakGroup.RightIndex : timePointCount - 1;
                    }

                    // ── All features now use scoringApexIdx / scoringLeft / scoringRight ──

                    result.ApexScore = CosineActiveFragments(
                        libIntensities, matrix, scoringApexIdx * fragmentCount, fragmentCount);
                    result.ObservedApexRt = refRts[scoringApexIdx];

                    // Temporal cosine over scoring peak window
                    ComputeTemporalCosineOnRange(
                        libIntensities, matrix, fragmentCount, timePointCount,
                        scoringLeft, scoringRight, minActiveFragments: 3,
                        out float peakTemporalScore, out int peakTimePointsUsed);
                    result.TemporalScore = peakTemporalScore;
                    result.TimePointsUsed = peakTimePointsUsed;
                    result.DotProductScore = result.TemporalScore;

                    // Spectral angle from peak-restricted temporal cosine
                    if (!float.IsNaN(result.TemporalScore))
                    {
                        float clampedCosine = Math.Clamp(result.TemporalScore, 0f, 1f);
                        result.SpectralAngleScore = 1.0f - (2.0f / MathF.PI) * MathF.Acos(clampedCosine);
                        result.SpectralAngle = result.SpectralAngleScore;
                    }

                    // Full-window fragment correlations (window-wide, independent of apex)
                    ComputeFragmentCorrelationsOnRange(
                        matrix, fragmentCount, timePointCount, 0, timePointCount - 1,
                        out float fullMeanCorr, out float fullMinCorr);
                    result.MeanFragCorr = fullMeanCorr;
                    result.MinFragCorr = fullMinCorr;

                    // Peak-boundary fragment correlations
                    result.PeakApexScore = result.ApexScore;
                    float peakWidthMinutes = (scoringRight < timePointCount && scoringLeft >= 0)
                        ? refRts[scoringRight] - refRts[scoringLeft] : 0f;

                    ComputeApexCenteredFragmentCorrelations(
                        matrix, fragmentCount, timePointCount,
                        scoringApexIdx, peakWidthMinutes, refRts,
                        out float peakMedCorr, out float peakMinCorr);
                    result.PeakMeanFragCorr = peakMedCorr;
                    result.PeakMinFragCorr = peakMinCorr;

                    // RT deviation against LOWESS-predicted RT
                    if (predictedRt.HasValue && !float.IsNaN(result.ObservedApexRt))
                    {
                        float deviation = result.ObservedApexRt - predictedRt.Value;
                        result.RtDeviationMinutes = deviation;
                        result.RtDeviationSquared = deviation * deviation;
                    }

                    // Count detected fragments across all scans
                    int detectedFragmentCount = 0;
                    for (int f = 0; f < fragmentCount; f++)
                        for (int s = 0; s < timePointCount; s++)
                            if (matrix[s * fragmentCount + f] > 0) { detectedFragmentCount++; break; }

                    // 1. Mass accuracy at apex
                    if (index != null)
                    {
                        int apexScanGlobalIndex = -1;
                        if (index.TryGetScanRangeForWindow(group.WindowId, out int winStart, out int winCount))
                            apexScanGlobalIndex = FindClosestScanByRt(
                                index, winStart, winCount, refRts[scoringApexIdx]);
                        if (apexScanGlobalIndex >= 0)
                            DiaMassAccuracyHelper.ComputeMassAccuracyAtApex(
                                result, input, index, apexScanGlobalIndex, parameters.PpmTolerance);
                    }

                    // 2. Signal ratio deviation at apex
                    float[] apexIntensities = ArrayPool<float>.Shared.Rent(fragmentCount);
                    try
                    {
                        Array.Clear(apexIntensities, 0, fragmentCount);
                        for (int f = 0; f < fragmentCount; f++)
                            apexIntensities[f] = matrix[scoringApexIdx * fragmentCount + f];
                        DiaSignalRatioHelper.ComputeSignalRatioFeatures(
                            apexIntensities, input.FragmentIntensities, fragmentCount, result);
                    }
                    finally
                    {
                        ArrayPool<float>.Shared.Return(apexIntensities);
                    }

                    // 3. Best-fragment reference curve
                    float peakWidthForBestFrag = (scoringRight < timePointCount && scoringLeft >= 0)
                        ? refRts[scoringRight] - refRts[scoringLeft] : 0f;
                    DiaBestFragmentHelper.ComputeBestFragmentFeatures(
                        matrix, fragmentCount, timePointCount, detectedFragmentCount, result,
                        input.FragmentIntensities.AsSpan(), scoringApexIdx * fragmentCount,
                        scoringApexIdx, peakWidthForBestFrag, refRts);

                    // 4. Peak shape ratios — use scoring bounds
                    DiaPeakShapeHelper.ComputePeakShapeRatios(
                        matrix, fragmentCount,
                        scoringLeft, scoringRight, scoringApexIdx,
                        result);

                    // 5. Smoothed fragment correlations (full window)
                    DiaSmoothedFeatureHelper.ComputeSmoothedCorrelationFeatures(
                        matrix, fragmentCount, timePointCount, result);

                    // 6. Signal-to-noise around scoring apex
                    DiaSmoothedFeatureHelper.ComputeSignalToNoise(
                        matrix, scoringApexIdx, fragmentCount, timePointCount, result);

                    // Peak shape metadata
                    if (peakGroup.IsValid)
                    {
                        result.CandidateCount = peakGroup.CandidateCount;
                        result.PeakWidth = peakWidthMinutes;
                    }

                    // 7. Library coverage fraction (window-wide)
                    {
                        float libTotalWeight = 0f, libDetectedWeight = 0f;
                        for (int f = 0; f < fragmentCount; f++)
                        {
                            float libW = input.FragmentIntensities[f];
                            if (libW <= 0f) continue;
                            libTotalWeight += libW;
                            if (result.XicPointCounts[f] > 0) libDetectedWeight += libW;
                        }
                        result.LibraryCoverageFraction = libTotalWeight > 0f
                            ? libDetectedWeight / libTotalWeight : float.NaN;
                    }

                    if (!float.IsNaN(result.DotProductScore)
                        && result.DotProductScore < parameters.MinScoreThreshold)
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
        /// Computes pairwise Pearson correlations between fragment XICs in a window
        /// centred on the chromatographic apex, sized to one measured peak half-width
        /// on each side.
        ///
        /// Window sizing:
        ///   scanInterval   = median inter-scan gap in refRts (robust to outliers)
        ///   halfWidthScans = max(MinApexHalfWidthScans,
        ///                        round(peakWidthMinutes / 2 / scanInterval))
        ///   range          = [apexIdx - halfWidthScans, apexIdx + halfWidthScans]
        ///
        /// This gives ~5-7 scans for our typical FWHM of 0.069 min at ~33 scans/min,
        /// replacing the 2-3 scan peak-boundary window (too narrow for reliableComputeFragmentCorrelationsOnRange
        /// Pearson) and the ~50 scan full-window (dominated by baseline noise).
        ///
        /// Summary statistic: MEDIAN pairwise correlation.
        ///   Mean is sensitive to one noisy fragment pair dragging it down.
        ///   Median is robust — one bad pair out of 28 has negligible effect.
        ///
        /// Fragment eligibility: >= 2 nonzero time points within the window.
        ///   (Threshold relaxed from 3 because the window is narrow.)
        ///
        /// Outputs PeakMeanFragCorr (now median) and PeakMinFragCorr unchanged
        /// so no downstream code needs to change.
        /// </summary>
        /// <param name="matrix">Flat row-major [timePoints × fragments] matrix.</param>
        /// <param name="fragmentCount">Number of fragment columns.</param>
        /// <param name="timePointCount">Number of time point rows.</param>
        /// <param name="apexIdx">Apex scan index within the local refRts array.</param>
        /// <param name="peakWidthMinutes">
        /// Measured peak width in minutes (refRts[RightIndex] - refRts[LeftIndex]).
        /// Used to derive halfWidthScans. Pass 0f to use MinApexHalfWidthScans only.
        /// </param>
        /// <param name="refRts">RT timestamps of the reference fragment (used to
        /// compute scan interval). Length must equal timePointCount.</param>
        /// <param name="medianCorrelation">Median pairwise Pearson r. NaN if < 2
        /// eligible fragments or < 1 valid pair.</param>
        /// <param name="minCorrelation">Minimum pairwise Pearson r across all pairs.
        /// NaN if < 1 valid pair.</param>
        private static void ComputeApexCenteredFragmentCorrelations(
            ReadOnlySpan<float> matrix,
            int fragmentCount,
            int timePointCount,
            int apexIdx,
            float peakWidthMinutes,
            ReadOnlySpan<float> refRts,
            out float medianCorrelation,
            out float minCorrelation)
        {
            // Minimum half-width: always use at least this many scans on each side
            // of the apex, regardless of measured peak width.
            const int MinApexHalfWidthScans = 3;

            medianCorrelation = float.NaN;
            minCorrelation = float.NaN;

            if (timePointCount < 2 || fragmentCount < 2) return;

            // ── Compute scan interval from refRts ────────────────────────────
            // Use median of consecutive gaps — robust to any duplicated RT values.
            float scanInterval;
            if (timePointCount >= 2)
            {
                // Collect gaps, take median
                Span<float> gaps = timePointCount <= 65
                    ? stackalloc float[timePointCount - 1]
                    : new float[timePointCount - 1];
                for (int i = 0; i < timePointCount - 1; i++)
                    gaps[i] = refRts[i + 1] - refRts[i];
                // Partial sort to find median (insertion sort on small span)
                var gapsCopy = gaps.ToArray();
                Array.Sort(gapsCopy);
                scanInterval = gapsCopy[gapsCopy.Length / 2];
                if (scanInterval <= 0f)
                    scanInterval = (refRts[timePointCount - 1] - refRts[0])
                                   / Math.Max(timePointCount - 1, 1);
            }
            else
            {
                scanInterval = 0.03f; // ~30ms fallback
            }

            if (scanInterval <= 0f) scanInterval = 0.03f;

            // ── Compute dynamic half-width in scans ──────────────────────────
            int halfWidthScans = MinApexHalfWidthScans;
            if (peakWidthMinutes > 0f)
            {
                int dynamicHalf = (int)Math.Round(peakWidthMinutes / 2f / scanInterval);
                halfWidthScans = Math.Max(MinApexHalfWidthScans, dynamicHalf);
            }

            int rangeStart = Math.Max(0, apexIdx - halfWidthScans);
            int rangeEnd = Math.Min(timePointCount - 1, apexIdx + halfWidthScans);
            int rangeLength = rangeEnd - rangeStart + 1;
            if (rangeLength < 2) return;

            // ── Find fragments with >= 2 nonzero points in range ─────────────
            // Threshold is 2 (not 3) because the window is narrow by design.
            Span<int> detectedFrags = stackalloc int[Math.Min(fragmentCount, 64)];
            int nDetected = 0;

            for (int f = 0; f < fragmentCount && nDetected < detectedFrags.Length; f++)
            {
                int nonzero = 0;
                for (int t = rangeStart; t <= rangeEnd; t++)
                    if (matrix[t * fragmentCount + f] > 0f) nonzero++;
                if (nonzero >= 2)
                    detectedFrags[nDetected++] = f;
            }

            if (nDetected < 2) return;

            // ── Compute all pairwise Pearson correlations ─────────────────────
            int maxPairs = nDetected * (nDetected - 1) / 2;
            float[] pairCorrs = new float[maxPairs];
            int nPairs = 0;
            float minCorr = float.MaxValue;

            for (int a = 0; a < nDetected; a++)
            {
                for (int b = a + 1; b < nDetected; b++)
                {
                    float r = PearsonCorrelationOnRange(
                        matrix, detectedFrags[a], detectedFrags[b],
                        fragmentCount, rangeStart, rangeEnd);

                    if (!float.IsNaN(r))
                    {
                        pairCorrs[nPairs++] = r;
                        if (r < minCorr) minCorr = r;
                    }
                }
            }

            if (nPairs == 0) return;

            // ── Median of pairwise correlations ───────────────────────────────
            // Sort only the populated slice.
            var populated = new float[nPairs];
            Array.Copy(pairCorrs, populated, nPairs);
            Array.Sort(populated);

            medianCorrelation = nPairs % 2 == 1
                ? populated[nPairs / 2]
                : (populated[nPairs / 2 - 1] + populated[nPairs / 2]) * 0.5f;

            minCorrelation = minCorr;
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
        /// precursor in the same DIA isolation window whose chromatographic peak overlaps
        /// with R's peak.
        ///
        ///   ChimericScore = uncontested_intensity / total_intensity
        ///
        /// A fragment of R is "contested" if any fragment of any other precursor Q in the
        /// same window falls within <paramref name="ppmTolerance"/> ppm of that fragment's m/z,
        /// AND Q's peak (apex ± halfFwhm) overlaps with R's peak.
        ///
        /// Co-eluter overlap uses ObservedApexRt ± PeakWidth/2 (falling back to a fixed
        /// 0.075 min half-width when no peak was detected). This is much narrower than the
        /// extraction window (±1 min) and correctly restricts interference detection to
        /// precursors that are physically co-eluting at the chromatographic peak level.
        ///
        /// Interpretation:
        ///   1.0 — no fragment m/z overlaps with any co-eluting precursor (clean precursor)
        ///   0.0 — every fragment is shared with at least one co-eluting precursor
        ///   x   — fraction x of the total signal comes from uncontested fragments
        ///
        /// When total extracted intensity is zero, ChimericScore is set to 1.0 (no evidence
        /// of contamination).
        ///
        /// Must be called after assembly (ExtractedIntensities and ObservedApexRt populated)
        /// and before FDR.
        /// </summary>
        /// <param name="precursors">Library precursor inputs carrying fragment m/z values.</param>
        /// <param name="results">Assembled results. ChimericScore is set in place.</param>
        /// <param name="ppmTolerance">
        /// PPM tolerance for fragment m/z overlap detection. A fragment is contested if
        /// any other co-eluting precursor in the same window has a fragment within this tolerance.
        /// </param>
        public static void ComputeChimericScores(
            IReadOnlyList<LibraryPrecursorInput> precursors,
            List<DiaSearchResult> results,
            float ppmTolerance)
        {
            if (results == null || results.Count == 0)
                return;

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
                var apexRt = new float[n];
                var halfFwhm = new float[n];

                const float fallbackHalfFwhm = 0.075f; // 2× median observed peak FWHM (~0.069 min)

                for (int a = 0; a < n; a++)
                {
                    var r = results[indices[a]];
                    precursorIndices[a] = r.PrecursorIndex;
                    apexRt[a] = r.ObservedApexRt;
                    halfFwhm[a] = (!float.IsNaN(r.PeakWidth) && r.PeakWidth > 0f)
                        ? r.PeakWidth * 0.5f
                        : fallbackHalfFwhm;
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

                    // If result A has no detected apex, fall back to extraction window overlap.
                    // This avoids silently assigning 1.0 to unscored results.
                    bool aHasApex = !float.IsNaN(apexRt[a]);
                    float aLo = aHasApex ? apexRt[a] - halfFwhm[a] : ra.RtWindowStart;
                    float aHi = aHasApex ? apexRt[a] + halfFwhm[a] : ra.RtWindowEnd;

                    var precA = precursors[piA];
                    float totalIntensity = 0f;
                    float uncontestedIntensity = 0f;

                    int fragCount = Math.Min(precA.FragmentCount, ra.ExtractedIntensities?.Length ?? 0);
                    for (int f = 0; f < fragCount; f++)
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

                            // Only consider co-eluters: precursors whose chromatographic peak
                            // overlaps with result A's peak. Uses apex ± halfFwhm, which is
                            // ~0.15 min wide — much narrower than the ±1 min extraction window.
                            // Falls back to extraction window bounds when no apex was detected.
                            bool bHasApex = !float.IsNaN(apexRt[b]);
                            float bLo = bHasApex ? apexRt[b] - halfFwhm[b] : results[indices[b]].RtWindowStart;
                            float bHi = bHasApex ? apexRt[b] + halfFwhm[b] : results[indices[b]].RtWindowEnd;

                            if (aHi <= bLo || bHi <= aLo)
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