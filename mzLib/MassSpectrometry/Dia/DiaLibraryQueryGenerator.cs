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
        /// Optional: iRT value for iRT-calibration-based workflows.
        /// Null if the library uses minutes-scale RT instead of iRT.
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
        /// Generates FragmentQuery[] using iRT calibration for per-precursor RT windows.
        /// 
        /// Instead of a fixed ±RtToleranceMinutes window, each precursor's RT window is
        /// computed from its iRT value through the calibration model:
        ///   predicted_RT = model.ToMinutes(precursor.IrtValue)
        ///   window = predicted_RT ± k * σ_minutes
        /// 
        /// where k = parameters.CalibratedWindowSigmaMultiplier (default 3.0)
        /// and σ_minutes = model.SigmaMinutes.
        /// 
        /// For precursors without an iRT value, falls back to RetentionTime ± RtToleranceMinutes.
        /// For precursors with neither iRT nor RT, uses the full run RT range.
        /// </summary>
        /// <param name="precursors">Library precursor inputs with optional IrtValue.</param>
        /// <param name="scanIndex">The DIA scan index for window lookup and RT bounds.</param>
        /// <param name="parameters">Search parameters including PpmTolerance and CalibratedWindowSigmaMultiplier.</param>
        /// <param name="calibration">
        /// The fitted iRT → RT calibration model. If null, falls back to Generate() behavior.
        /// </param>
        public static GenerationResult GenerateCalibrated(
            IList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaSearchParameters parameters,
            RtCalibrationModel calibration)
        {
            if (calibration == null)
                return Generate(precursors, scanIndex, parameters);

            if (precursors == null) throw new ArgumentNullException(nameof(precursors));
            if (scanIndex == null) throw new ArgumentNullException(nameof(scanIndex));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));

            float ppmTolerance = parameters.PpmTolerance;
            float fixedRtTolerance = parameters.RtToleranceMinutes;
            double k = parameters.CalibratedWindowSigmaMultiplier;
            float calibratedHalfWidth = (float)calibration.GetMinutesWindowHalfWidth(k);

            // Fallback RT bounds for precursors with neither iRT nor RT
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

            // Second pass: fill
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
                    // Use calibrated iRT → RT conversion with data-driven window
                    float predictedRt = (float)calibration.ToMinutes(p.IrtValue.Value);
                    rtMin = predictedRt - calibratedHalfWidth;
                    rtMax = predictedRt + calibratedHalfWidth;
                }
                else if (p.RetentionTime.HasValue)
                {
                    // Fallback: fixed window around library RT
                    float rt = (float)p.RetentionTime.Value;
                    rtMin = rt - fixedRtTolerance;
                    rtMax = rt + fixedRtTolerance;
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
        /// This overload uses TotalIntensity from FragmentResult for scoring.
        /// 
        /// For better scoring accuracy, prefer the overload that accepts XIC buffers
        /// which can compute apex (max single-scan) intensity per fragment.
        /// </summary>
        public static List<DiaSearchResult> AssembleResults(
            IList<LibraryPrecursorInput> precursors,
            GenerationResult generationResult,
            FragmentResult[] extractionResults,
            DiaSearchParameters parameters,
            IScorer dotProductScorer = null,
            IScorer spectralAngleScorer = null)
        {
            return AssembleResults(
                precursors, generationResult, extractionResults, parameters,
                intensityBuffer: null, dotProductScorer, spectralAngleScorer);
        }

        /// <summary>
        /// Assembles DiaSearchResult objects from extraction results with XIC buffer access.
        /// 
        /// When parameters.UseApexIntensityForScoring is true and intensityBuffer is provided,
        /// uses the maximum single-scan intensity per fragment for scoring instead of the 
        /// summed TotalIntensity. This produces much better library-vs-extracted correlation
        /// because:
        ///   - Library relative intensities represent the fragment pattern at peak apex
        ///   - Summing over the entire RT window accumulates noise and co-eluting interference
        ///   - Apex intensity captures the true abundance ratio at the chromatographic peak
        /// 
        /// The ExtractedIntensities[] in DiaSearchResult will contain either apex or total
        /// intensity depending on the UseApexIntensityForScoring parameter.
        /// </summary>
        /// <param name="precursors">Original precursor inputs (for library intensities and metadata).</param>
        /// <param name="generationResult">Query generation output mapping precursors to queries.</param>
        /// <param name="extractionResults">Per-query extraction results from ExtractBatch.</param>
        /// <param name="parameters">Search parameters including UseApexIntensityForScoring.</param>
        /// <param name="intensityBuffer">
        /// XIC intensity buffer from extraction (parallel to RT buffer).
        /// If null, TotalIntensity is used regardless of UseApexIntensityForScoring.
        /// </param>
        /// <param name="dotProductScorer">Optional dot product scorer.</param>
        /// <param name="spectralAngleScorer">Optional spectral angle scorer.</param>
        public static List<DiaSearchResult> AssembleResults(
            IList<LibraryPrecursorInput> precursors,
            GenerationResult generationResult,
            FragmentResult[] extractionResults,
            DiaSearchParameters parameters,
            float[] intensityBuffer,
            IScorer dotProductScorer = null,
            IScorer spectralAngleScorer = null)
        {
            bool useApex = parameters.UseApexIntensityForScoring && intensityBuffer != null;

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

                int detected = 0;
                for (int f = 0; f < group.QueryCount; f++)
                {
                    int qi = group.QueryOffset + f;
                    var fr = extractionResults[qi];
                    result.XicPointCounts[f] = fr.DataPointCount;

                    if (fr.DataPointCount > 0)
                    {
                        detected++;

                        if (useApex)
                        {
                            // Compute apex (max) intensity from the XIC buffer
                            float maxIntensity = 0f;
                            int start = fr.IntensityBufferOffset;
                            int end = start + fr.DataPointCount;
                            for (int p = start; p < end; p++)
                            {
                                if (intensityBuffer[p] > maxIntensity)
                                    maxIntensity = intensityBuffer[p];
                            }
                            result.ExtractedIntensities[f] = maxIntensity;
                        }
                        else
                        {
                            result.ExtractedIntensities[f] = fr.TotalIntensity;
                        }
                    }
                    // else: ExtractedIntensities[f] remains 0
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