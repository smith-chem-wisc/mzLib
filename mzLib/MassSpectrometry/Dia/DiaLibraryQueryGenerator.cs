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
        /// Fragment ion m/z values from the library spectrum.
        /// Length determines the number of queries generated for this precursor.
        /// </summary>
        public readonly float[] FragmentMzs;

        /// <summary>
        /// Fragment ion intensities from the library spectrum (parallel to FragmentMzs).
        /// Used downstream for scoring (library vs extracted intensity comparison).
        /// </summary>
        public readonly float[] FragmentIntensities;

        /// <summary>
        /// Indexed retention time (iRT) value from the library, used for RT calibration.
        /// Null if the library does not provide iRT values. When present, this is the
        /// normalized/indexed RT (e.g., Biognosys iRT scale) as opposed to RetentionTime
        /// which is in run-specific minutes. The calibration system maps iRT → run RT.
        /// </summary>
        public readonly double? IrtValue;

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
            FragmentMzs = fragmentMzs ?? throw new ArgumentNullException(nameof(fragmentMzs));
            FragmentIntensities = fragmentIntensities ?? throw new ArgumentNullException(nameof(fragmentIntensities));
            IrtValue = irtValue;

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
        /// Generates FragmentQuery[] using an iRT calibration model to compute per-precursor
        /// RT windows instead of a fixed ±tolerance.
        /// 
        /// For each precursor with an IrtValue:
        ///   1. Converts iRT → predicted run RT via calibration.Predict(iRT)
        ///   2. Computes the RT window as predicted ± (σ × CalibratedWindowSigmaMultiplier)
        ///      where σ = calibration.SigmaMinutes (residual standard deviation of the fit)
        /// 
        /// Precursors without IrtValue fall back to RetentionTime ± RtToleranceMinutes,
        /// same as the uncalibrated Generate() method.
        /// 
        /// This typically narrows RT windows from ±5 min to ±0.9 min (a ~5–20× reduction),
        /// dramatically improving both extraction speed and scoring quality.
        /// </summary>
        /// <param name="precursors">Library precursor inputs (must have IrtValue set for calibration to apply).</param>
        /// <param name="scanIndex">The DIA scan index for window lookups.</param>
        /// <param name="parameters">Search parameters. Uses PpmTolerance and RtToleranceMinutes (as fallback).</param>
        /// <param name="calibration">
        /// The fitted iRT → RT calibration model. Provides Predict(irt) and SigmaMinutes.
        /// If null, behaves identically to Generate().
        /// </param>
        public static GenerationResult GenerateCalibrated(
            IList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaSearchParameters parameters,
            RtCalibrationModel calibration)
        {
            if (precursors == null) throw new ArgumentNullException(nameof(precursors));
            if (scanIndex == null) throw new ArgumentNullException(nameof(scanIndex));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));

            // If no calibration model, fall back to fixed-window generation
            if (calibration == null)
                return Generate(precursors, scanIndex, parameters);

            float ppmTolerance = parameters.PpmTolerance;
            float fixedRtTolerance = parameters.RtToleranceMinutes;
            float globalRtMin = scanIndex.GetGlobalRtMin();
            float globalRtMax = scanIndex.GetGlobalRtMax();

            // Calibrated window half-width: σ × multiplier (default 3σ ≈ 99.7% coverage)
            double sigmaMultiplier = parameters.CalibratedWindowSigmaMultiplier > 0
                ? parameters.CalibratedWindowSigmaMultiplier
                : 3.0;
            float calibratedHalfWidth = (float)calibration.GetMinutesWindowHalfWidth(sigmaMultiplier);

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
                    // Use calibration model: iRT → predicted RT, then ± calibrated half-width
                    float predictedRt = (float)calibration.ToMinutes(p.IrtValue.Value);
                    rtMin = predictedRt - calibratedHalfWidth;
                    rtMax = predictedRt + calibratedHalfWidth;
                }
                else if (p.RetentionTime.HasValue)
                {
                    // Fallback: use raw RT with fixed tolerance
                    float rt = (float)p.RetentionTime.Value;
                    rtMin = rt - fixedRtTolerance;
                    rtMax = rt + fixedRtTolerance;
                }
                else
                {
                    // No RT info at all: search entire run
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
        /// Assembles DiaSearchResult objects from extraction results using summed TotalIntensity.
        /// This is the original scoring mode — summing all XIC data points per fragment.
        /// 
        /// For better scoring accuracy, prefer the overload that accepts intensityBuffer
        /// to use apex (max single-scan) intensity instead.
        /// </summary>
        public static List<DiaSearchResult> AssembleResults(
            IList<LibraryPrecursorInput> precursors,
            GenerationResult generationResult,
            FragmentResult[] extractionResults,
            DiaSearchParameters parameters,
            IScorer dotProductScorer = null,
            IScorer spectralAngleScorer = null)
        {
            // Delegate to the apex-capable overload with null intensityBuffer,
            // which falls back to TotalIntensity scoring.
            return AssembleResults(
                precursors, generationResult, extractionResults,
                parameters, intensityBuffer: null,
                dotProductScorer, spectralAngleScorer);
        }

        /// <summary>
        /// Assembles DiaSearchResult objects using consensus-apex intensity scoring.
        /// 
        /// When rtBuffer and intensityBuffer are both provided, scoring uses the
        /// "consensus apex" approach:
        ///   1. For each precursor, build a per-RT intensity profile by summing all
        ///      fragment intensities at each unique RT point
        ///   2. Find the RT with the highest total fragment intensity (the consensus apex)
        ///   3. Extract each fragment's intensity at that specific RT
        /// 
        /// This is superior to per-fragment-max-apex because:
        ///   - All fragments are sampled from the same scan (coelution point)
        ///   - The relative intensities match what a single-scan library spectrum represents
        ///   - Interference in one fragment doesn't pull the apex to the wrong scan
        /// 
        /// When buffers are null, falls back to TotalIntensity (summed) scoring.
        /// </summary>
        /// <param name="precursors">Original library precursor inputs.</param>
        /// <param name="generationResult">Query generation metadata (from Generate()).</param>
        /// <param name="extractionResults">Per-query extraction results.</param>
        /// <param name="parameters">Search parameters.</param>
        /// <param name="intensityBuffer">
        /// XIC intensity buffer from ExtractionResult.IntensityBuffer. Pass null for summed scoring.
        /// </param>
        /// <param name="dotProductScorer">Optional normalized dot product scorer.</param>
        /// <param name="spectralAngleScorer">Optional spectral angle scorer.</param>
        /// <param name="rtBuffer">
        /// XIC RT buffer from ExtractionResult.RtBuffer. Required for consensus-apex scoring.
        /// If null but intensityBuffer is provided, falls back to per-fragment max apex.
        /// </param>
        public static List<DiaSearchResult> AssembleResults(
            IList<LibraryPrecursorInput> precursors,
            GenerationResult generationResult,
            FragmentResult[] extractionResults,
            DiaSearchParameters parameters,
            float[] intensityBuffer,
            IScorer dotProductScorer = null,
            IScorer spectralAngleScorer = null,
            float[] rtBuffer = null)
        {
            bool useConsensusApex = intensityBuffer != null && rtBuffer != null;
            bool usePerFragApex = intensityBuffer != null && rtBuffer == null;
            var results = new List<DiaSearchResult>(generationResult.PrecursorGroups.Length);

            // Reusable lookup structures for consensus apex (avoid per-precursor allocation)
            // Key: RT value (float bits as int for exact matching), Value: summed intensity
            var rtToSumIntensity = new Dictionary<int, float>();
            var rtToFragIntensities = new Dictionary<int, float[]>();

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

                int detected = 0;

                if (useConsensusApex)
                {
                    // ── Consensus apex scoring ──────────────────────────────
                    // Pass 1: Build RT → per-fragment intensity map
                    rtToSumIntensity.Clear();
                    rtToFragIntensities.Clear();

                    for (int f = 0; f < group.QueryCount; f++)
                    {
                        int qi = group.QueryOffset + f;
                        var fr = extractionResults[qi];
                        if (fr.DataPointCount == 0) continue;

                        for (int p = 0; p < fr.DataPointCount; p++)
                        {
                            float rt = rtBuffer[fr.RtBufferOffset + p];
                            float intensity = intensityBuffer[fr.IntensityBufferOffset + p];
                            int rtKey = BitConverter.SingleToInt32Bits(rt);

                            if (!rtToSumIntensity.ContainsKey(rtKey))
                            {
                                rtToSumIntensity[rtKey] = 0f;
                                rtToFragIntensities[rtKey] = new float[group.QueryCount];
                            }

                            rtToSumIntensity[rtKey] += intensity;
                            rtToFragIntensities[rtKey][f] = intensity;
                        }
                    }

                    // Pass 2: Find the RT with the highest sum (consensus apex)
                    int bestRtKey = 0;
                    float bestSum = -1f;
                    foreach (var kvp in rtToSumIntensity)
                    {
                        if (kvp.Value > bestSum)
                        {
                            bestSum = kvp.Value;
                            bestRtKey = kvp.Key;
                        }
                    }

                    // Pass 3: Extract each fragment's intensity at the consensus apex RT
                    if (bestSum > 0f && rtToFragIntensities.TryGetValue(bestRtKey, out var apexIntensities))
                    {
                        for (int f = 0; f < group.QueryCount; f++)
                        {
                            result.ExtractedIntensities[f] = apexIntensities[f];
                        }
                    }

                    // Count detected fragments and XIC point counts
                    for (int f = 0; f < group.QueryCount; f++)
                    {
                        int qi = group.QueryOffset + f;
                        var fr = extractionResults[qi];
                        result.XicPointCounts[f] = fr.DataPointCount;
                        if (fr.DataPointCount > 0)
                            detected++;
                    }
                }
                else
                {
                    // ── Per-fragment max apex or summed scoring ──────────────
                    for (int f = 0; f < group.QueryCount; f++)
                    {
                        int qi = group.QueryOffset + f;
                        var fr = extractionResults[qi];

                        if (usePerFragApex && fr.DataPointCount > 0)
                        {
                            float maxIntensity = 0f;
                            int offset = fr.IntensityBufferOffset;
                            for (int p = 0; p < fr.DataPointCount; p++)
                            {
                                float val = intensityBuffer[offset + p];
                                if (val > maxIntensity)
                                    maxIntensity = val;
                            }
                            result.ExtractedIntensities[f] = maxIntensity;
                        }
                        else
                        {
                            result.ExtractedIntensities[f] = fr.TotalIntensity;
                        }

                        result.XicPointCounts[f] = fr.DataPointCount;
                        if (fr.DataPointCount > 0)
                            detected++;
                    }
                }

                result.FragmentsDetected = detected;

                if (!result.MeetsMinFragments(parameters.MinFragmentsRequired))
                    continue;

                if (dotProductScorer != null)
                    result.DotProductScore = dotProductScorer.Score(
                        input.FragmentIntensities, result.ExtractedIntensities);

                if (spectralAngleScorer != null)
                    result.SpectralAngleScore = spectralAngleScorer.Score(
                        input.FragmentIntensities, result.ExtractedIntensities);

                if (!float.IsNaN(result.DotProductScore) && result.DotProductScore < parameters.MinScoreThreshold)
                    continue;

                results.Add(result);
            }

            return results;
        }
    }
}