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

        /// <summary>
        /// Library iRT value (dimensionless indexed retention time scale).
        /// Predicted by Koina/Prosit or loaded from spectral library.
        /// Null if no iRT is available, in which case RetentionTime (minutes) is used.
        /// When iRT calibration is enabled, this is the primary RT coordinate for
        /// candidate selection and RT scoring.
        /// </summary>
        public readonly double? IrtValue;

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
            IrtValue = irtValue;
            IsDecoy = isDecoy;
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
    /// Supports two modes:
    /// 
    /// 1. **Fixed RT window** (original behavior, UseIrtCalibration = false):
    ///    Uses RetentionTime ± RtToleranceMinutes as extraction window.
    /// 
    /// 2. **iRT-calibrated** (UseIrtCalibration = true):
    ///    Phase 1 (broad): provisional iRT mapping + wide window for anchor collection.
    ///    Phase 2+ (refined): calibrated iRT window = k * σ_iRT centered on predicted RT.
    /// 
    /// Also produces PrecursorQueryGroup[] metadata that records which queries belong to
    /// which precursor, enabling downstream result assembly and RT scoring.
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

            /// <summary>RT window lower bound used for queries (minutes)</summary>
            public readonly float RtMin;

            /// <summary>RT window upper bound used for queries (minutes)</summary>
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
        /// Generates FragmentQuery[] from library precursor inputs using fixed RT windows.
        /// This is the original behavior (no iRT calibration).
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
        /// Generates FragmentQuery[] using iRT calibration for RT windowing.
        /// 
        /// The calibration model transforms each precursor's iRT value into an
        /// observed RT window in minutes:
        ///   center_minutes = calibration.ToMinutes(precursor.IrtValue)
        ///   halfWidth_minutes = calibration.GetMinutesWindowHalfWidth(k)
        /// 
        /// Precursors without iRT values fall back to RetentionTime ± RtToleranceMinutes.
        /// </summary>
        /// <param name="precursors">Library precursors to generate queries for.</param>
        /// <param name="scanIndex">DIA scan index (for window ID mapping and RT bounds).</param>
        /// <param name="parameters">Search parameters.</param>
        /// <param name="calibration">
        /// The iRT calibration model. If null or unreliable, falls back to fixed RT windows.
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

            // If no calibration or unreliable, fall back to original method
            if (calibration == null || !calibration.IsReliable)
                return Generate(precursors, scanIndex, parameters);

            float ppmTolerance = parameters.PpmTolerance;
            double k = parameters.CalibratedWindowSigmaMultiplier;
            float fallbackRtTolerance = parameters.RtToleranceMinutes;

            // Calibrated window half-width in minutes
            float windowHalfMinutes = (float)calibration.GetMinutesWindowHalfWidth(k);

            // Fallback RT bounds for precursors without any RT info
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

            // Second pass: fill with calibrated windows
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
                    // Use calibration to convert iRT → minutes center, then apply calibrated window
                    float centerMinutes = (float)calibration.ToMinutes(p.IrtValue.Value);
                    rtMin = centerMinutes - windowHalfMinutes;
                    rtMax = centerMinutes + windowHalfMinutes;
                }
                else if (p.RetentionTime.HasValue)
                {
                    // Fallback: use minutes-based RT with fixed tolerance
                    float rt = (float)p.RetentionTime.Value;
                    rtMin = rt - fallbackRtTolerance;
                    rtMax = rt + fallbackRtTolerance;
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
        /// Assembles DiaSearchResult objects from extraction results.
        /// Call after ExtractBatch() to produce the final result list.
        /// 
        /// When a calibration model is provided, computes RT scores and combined scores.
        /// </summary>
        public static List<DiaSearchResult> AssembleResults(
            IList<LibraryPrecursorInput> precursors,
            GenerationResult generationResult,
            FragmentResult[] extractionResults,
            DiaSearchParameters parameters,
            IScorer dotProductScorer = null,
            IScorer spectralAngleScorer = null,
            RtCalibrationModel calibration = null)
        {
            var results = new List<DiaSearchResult>(generationResult.PrecursorGroups.Length);
            double lambda = parameters.RtScoreLambda;
            bool hasCalibration = calibration != null && calibration.IsReliable;

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
                    rtWindowEnd: group.RtMax,
                    libraryIrt: input.IrtValue
                );

                // Populate fragment extraction data
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

                // Spectral scoring
                if (dotProductScorer != null)
                    result.DotProductScore = dotProductScorer.Score(
                        input.FragmentIntensities, result.ExtractedIntensities);

                if (spectralAngleScorer != null)
                    result.SpectralAngleScore = spectralAngleScorer.Score(
                        input.FragmentIntensities, result.ExtractedIntensities);

                // RT scoring (when calibration is available and precursor has iRT)
                if (hasCalibration && input.IrtValue.HasValue)
                {
                    // Compute calibrated iRT at the center of the extraction window
                    double extractionCenterMinutes = (group.RtMin + group.RtMax) / 2.0;
                    double calibratedIrt = calibration.ToIrt(extractionCenterMinutes);
                    double residualIrt = calibratedIrt - input.IrtValue.Value;

                    result.CalibratedIrt = calibratedIrt;
                    result.IrtResidual = residualIrt;
                    result.RtScore = (float)calibration.ComputeRtScore(residualIrt);
                }

                // Combined score
                float spectralScore = !float.IsNaN(result.SpectralAngleScore)
                    ? result.SpectralAngleScore
                    : result.DotProductScore;

                if (!float.IsNaN(spectralScore))
                {
                    if (!float.IsNaN(result.RtScore) && lambda > 0)
                    {
                        result.CombinedScore = (float)(spectralScore + lambda * result.RtScore);
                    }
                    else
                    {
                        result.CombinedScore = spectralScore;
                    }
                }

                // Score threshold filtering
                float filterScore = !float.IsNaN(result.CombinedScore)
                    ? result.CombinedScore
                    : (!float.IsNaN(result.DotProductScore) ? result.DotProductScore : float.NaN);

                if (!float.IsNaN(filterScore) && filterScore < parameters.MinScoreThreshold)
                    continue;

                results.Add(result);
            }

            return results;
        }

        /// <summary>
        /// Extracts anchor matches from Phase 1 results for calibration fitting.
        /// 
        /// Selects the best-scoring non-decoy result per unique sequence, then returns
        /// parallel arrays of (library_iRT, observed_RT_minutes) for calibration fitting.
        /// 
        /// Only results with both a library iRT and a dot product score above the
        /// anchor threshold are considered.
        /// </summary>
        /// <param name="results">Phase 1 search results.</param>
        /// <param name="precursors">Original precursor inputs (for iRT values).</param>
        /// <param name="generationResult">Query generation metadata (for RT window centers).</param>
        /// <param name="minScore">Minimum spectral score for anchor eligibility.</param>
        /// <param name="libraryIrts">Output: library iRT values of anchors.</param>
        /// <param name="observedRtMinutes">Output: observed RT centers (minutes) of anchors.</param>
        public static void ExtractAnchors(
            IList<DiaSearchResult> results,
            float minScore,
            out double[] libraryIrts,
            out double[] observedRtMinutes)
        {
            // Best score per sequence (unique peptidoform)
            var bestBySequence = new Dictionary<string, DiaSearchResult>();

            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];

                // Only targets with iRT and decent score
                if (r.IsDecoy) continue;
                if (!r.LibraryIrt.HasValue) continue;

                float score = !float.IsNaN(r.SpectralAngleScore)
                    ? r.SpectralAngleScore
                    : r.DotProductScore;
                if (float.IsNaN(score) || score < minScore) continue;

                string key = r.Sequence + "/" + r.ChargeState;
                if (!bestBySequence.TryGetValue(key, out var existing))
                {
                    bestBySequence[key] = r;
                }
                else
                {
                    float existingScore = !float.IsNaN(existing.SpectralAngleScore)
                        ? existing.SpectralAngleScore
                        : existing.DotProductScore;
                    if (score > existingScore)
                        bestBySequence[key] = r;
                }
            }

            // Build anchor arrays
            var irtList = new List<double>(bestBySequence.Count);
            var rtList = new List<double>(bestBySequence.Count);

            foreach (var kvp in bestBySequence)
            {
                var r = kvp.Value;
                // Observed RT = center of the extraction window
                double observedRt = (r.RtWindowStart + r.RtWindowEnd) / 2.0;
                irtList.Add(r.LibraryIrt.Value);
                rtList.Add(observedRt);
            }

            libraryIrts = irtList.ToArray();
            observedRtMinutes = rtList.ToArray();
        }
    }
}
