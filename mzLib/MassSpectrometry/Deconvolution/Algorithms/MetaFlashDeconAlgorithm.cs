// MetaFlashDeconAlgorithm.cs
//
// Status: STEPS 1–5 COMPLETE — full spectral deconvolution implemented
//
// Step 1 — log transform + universal/harmonic pattern construction
// Step 2 — binned candidate mass finding (single-round mass bin computation)
// Step 3 — isotope peak recruitment per candidate per charge
//          + per-charge cosine and noise power collection for Qscoring
// Step 4 — cosine scoring against Averagine; isotope offset optimisation
// Step 5 — IsotopicEnvelope output, deduplication, Qscoring
//
// Reference (paper):
//   Jeong et al. (2020) Cell Systems 10(2):213-218.e6  DOI: 10.1016/j.cels.2020.01.003
//
// Reference (implementation — for faithful reproduction + later fidelity evaluation):
//   OpenMS FLASHDeconv per-spectrum source @ E:\GitClones\OpenMS:
//     src/openms/source/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.cpp  (main pipeline)
//     src/openms/source/ANALYSIS/TOPDOWN/PeakGroup.cpp             (scoring / SNR)
//     src/openms/source/ANALYSIS/TOPDOWN/Qscore.cpp                (logistic Qscore)
//     src/topp/FLASHDeconv.cpp                                     (TOPP defaults)
//   Faithful mapping doc (mzLib step -> OpenMS file:line crosswalk, with the cosine /
//   tol_div_factor / overlap-dedup tolerances we reproduce):
//     E:\CodeReview\MetaFlashDecon\deliverables\reference_flashdecon_faithful_mapping.md
//
// Faithful tolerances (scoped to MetaFlashDecon only; see MetaFlashDeconParameters):
//   - MinCosineScore        = 0.85       (OpenMS min_isotope_cosine MS1; FLASHDeconv.cpp:193)
//   - TolDivFactor          = 2.5        (OpenMS tol_div_factor; FLASHDeconvAlgorithm.cpp:22)
//   - OverlapDedupTolFactor = 1.5        (OpenMS final dedup = 1.5 x ppm; :1223)
//   Candidate voting bins at ppm/2.5; survivors de-duplicated at 1.5 x ppm.

using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;                           // Constants.ProtonMass, C13MinusC12, ToMz/ToMass
using MassSpectrometry.MzSpectra;          // SpectralSimilarity.CosineOfAlignedVectors
using MzLibUtil;                           // PpmTolerance

namespace MassSpectrometry
{
    internal class MetaFlashDeconAlgorithm : DeconvolutionAlgorithm
    {
        // ── Algorithm constants ───────────────────────────────────────────────

        /// <summary>Harmonic denominators from the paper (½, ⅓, ⅖).</summary>
        private static readonly int[] HarmonicDenominators = { 2, 3, 5 };

        /// <summary>
        /// Minimum continuous support peaks to mark a mass bin as a candidate.
        /// Note: because the first peak of a run is credited retroactively when
        /// the second arrives, reaching spc==3 means 3 real peaks confirmed.
        /// </summary>
        private const int MinSupportPeakCount = 3;

        /// <summary>
        /// Maximum intensity ratio between adjacent charge states before the
        /// charge run is considered broken and the support count resets.
        /// </summary>
        private const float MaxChargeIntensityRatio = 10.0f;

        // Candidate-mass voting uses bins FINER than the requested tolerance, exactly as
        // OpenMS does: it divides the ppm tolerance by tol_div_factor (=2.5) before computing
        // the bin multiplication factor (FLASHDeconvAlgorithm.cpp:22; updateMembers_:171:
        //   j *= 1e-6; j /= tol_div_factor; bin_mul_factors_.push_back(1.0 / j);).
        // The wider input tolerance is re-applied only at the final overlap dedup
        // (removeOverlappingPeakGroups_, FLASHDeconvAlgorithm.cpp:1223).
        // The single-arg overload (tolDivFactor implicitly 1.0) is retained for callers
        // and tests that compute the bin factor directly at the raw ppm tolerance.
        private static double ComputeBinMulFactor(double tolerancePpm)
            => ComputeBinMulFactor(tolerancePpm, 1.0);

        private static double ComputeBinMulFactor(double tolerancePpm, double tolDivFactor)
        {
            double effectiveTol = tolerancePpm * 1e-6;
            if (tolDivFactor > 0.0) effectiveTol /= tolDivFactor;
            return 1.0 / effectiveTol;
        }

        // ── Constructor ───────────────────────────────────────────────────────
        internal MetaFlashDeconAlgorithm(DeconvolutionParameters deconParameters)
            : base(deconParameters)
        {
        }

        // ══════════════════════════════════════════════════════════════════════
        // ENTRY POINT
        // ══════════════════════════════════════════════════════════════════════

        protected internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var p = DeconvolutionParameters as MetaFlashDeconParameters
                ?? throw new MzLibException("Deconvolution params and algorithm do not match");

            range ??= spectrum.Range;
            if (spectrum.Size == 0) return Enumerable.Empty<IsotopicEnvelope>();

            // ── Step 1 ────────────────────────────────────────────────────────
            var logPeaks = BuildLogMzPeaks(spectrum, range, p.Polarity);
            if (logPeaks.Count == 0) return Enumerable.Empty<IsotopicEnvelope>();

            int minAbsCharge = Math.Abs(p.MinAssumedChargeState);
            int maxAbsCharge = Math.Abs(p.MaxAssumedChargeState);
            int chargeRange = maxAbsCharge - minAbsCharge + 1;

            double[] universalPattern = BuildUniversalPattern(minAbsCharge, chargeRange);
            double[][] harmonicPatterns = BuildHarmonicPatterns(minAbsCharge, chargeRange);
            // Faithful FLASHDeconv: bin the log-m/z axis at ppm / tol_div_factor (finer than
            // the requested tolerance) for the candidate voting (FLASHDeconvAlgorithm.cpp:22, :171).
            double binMulFactor = ComputeBinMulFactor(p.DeconvolutionTolerancePpm, p.TolDivFactor);

            // ── Step 2 ────────────────────────────────────────────────────────
            var candidates = FindCandidateMasses(
                logPeaks, universalPattern, harmonicPatterns, binMulFactor, p);
            if (candidates.Count == 0) return Enumerable.Empty<IsotopicEnvelope>();

            // ── Steps 3–5: recruit → score → dedup → Qscore ──────────────────
            // ScoreAndBuildEnvelopes collects:
            //   • exact ppm error (isotope index n known during recruitment)
            //   • per-charge cosine for the representative charge (f[1])
            //   • per-charge noise power for the representative charge (f[2]/f[3])
            //   • per-charge summed signal squared for the SNR denominator
            var scoringData = ScoreAndBuildEnvelopes(spectrum, candidates, p).ToList();

            // Deduplicate on envelope identity; propagate the full EnvelopeScoringData
            // through so the scorer receives all per-charge fields for the winner.
            // Faithful FLASHDeconv re-applies the WIDE tolerance only here: the final
            // overlap-dedup window is OverlapDedupTolFactor × the input ppm
            // (= 1.5 × ppm, FLASHDeconvAlgorithm.cpp:1223), in contrast to the finer
            // ppm / tol_div_factor used for candidate voting above.
            double dedupTolPpm = p.DeconvolutionTolerancePpm * p.OverlapDedupTolFactor;
            var dedupedData = MetaFlashDeconDeduplicator.Deduplicate(scoringData, dedupTolPpm);

            return MetaFlashDeconScorer.AssignQscores(dedupedData);
        }

        // ══════════════════════════════════════════════════════════════════════
        // STEP 3 — ISOTOPE RECRUITMENT
        // STEP 4 — COSINE SCORING
        // STEP 5 — IsotopicEnvelope OUTPUT
        // ══════════════════════════════════════════════════════════════════════

        /// <summary>
        /// For each candidate mass, recruits isotope peaks from the raw spectrum
        /// for each charge state in the candidate's range, scores the resulting
        /// isotope intensity distribution against Averagine, and returns accepted
        /// envelopes as <see cref="MetaFlashDeconScorer.EnvelopeScoringData"/>.
        ///
        /// During recruitment the following are accumulated simultaneously:
        /// <list type="bullet">
        ///   <item>Global per-isotope intensity vector (all charges combined) for cosine scoring.</item>
        ///   <item>Per-charge isotope intensity vectors for the representative charge.</item>
        ///   <item>Noise peaks — raw spectrum peaks in the same m/z window as each expected
        ///     isotope position but outside the ±ppm tolerance band — for the representative charge.</item>
        /// </list>
        /// </summary>
        private IEnumerable<MetaFlashDeconScorer.EnvelopeScoringData> ScoreAndBuildEnvelopes(
            MzSpectrum spectrum,
            List<CandidateMass> candidates,
            MetaFlashDeconParameters p)
        {
            var results = new List<MetaFlashDeconScorer.EnvelopeScoringData>();
            var tolerance = new PpmTolerance(p.DeconvolutionTolerancePpm);
            int polSign = Math.Sign((int)p.Polarity); // +1 positive, -1 negative

            foreach (var candidate in candidates)
            {
                // ── Step 3a: get Averagine shape for this mass ─────────────────
                // GetAllTheoreticalMasses/Intensities returns arrays sorted
                // INTENSITY-DESCENDING (index 0 = most abundant / apex isotope).
                // We re-sort to MASS-ASCENDING for use in the per-isotope cosine vector.
                //
                // DiffToMonoisotopic[i] = MostIntenseMasses[i] - monoisotopicMass
                //                       = apexMass - monoMass (already computed, exact)
                int avgIdx = p.AverageResidueModel.GetMostIntenseMassIndex(candidate.Mass);
                double apexDaFromMono = p.AverageResidueModel.GetDiffToMonoisotopic(avgIdx);

                double[] avgMasses = p.AverageResidueModel.GetAllTheoreticalMasses(avgIdx);
                double[] avgIntensRaw = p.AverageResidueModel.GetAllTheoreticalIntensities(avgIdx);

                // Sort to mass-ascending order for cosine scoring vector
                var sortedAvg = avgMasses.Zip(avgIntensRaw).OrderBy(pair => pair.First).ToArray();
                double[] avgIntens = sortedAvg.Select(pair => pair.Second).ToArray();
                int nIso = avgIntens.Length;

                // apexIdx in the mass-ascending array — centers the perIsotopeObs vector
                int apexIdx = Array.IndexOf(avgIntens, avgIntens.Max());

                // Representative charge: midpoint of the candidate's charge range.
                // Used for per-charge cosine and noise collection (f[1]–f[3]).
                int repAbsCharge = (candidate.MinAbsCharge + candidate.MaxAbsCharge) / 2;

                // ── Step 3b: recruit peaks per charge state ────────────────────
                int isoWindowSize = nIso + apexIdx;
                var perIsotopeObs = new double[isoWindowSize]; // global (all charges)
                var repChargeObs = new double[isoWindowSize]; // representative charge only

                // recruitedPeaks: (mz, intensity, absCharge, isotopeIndex n from monoisotopic)
                var recruitedPeaks = new List<(double mz, double intensity, int z, int n)>();

                // Per-charge noise and signal accumulators (representative charge only).
                // noisePower: summed squared intensity of peaks in the isotope window
                //             but OUTSIDE the ±ppm tolerance band.
                // sumSignalSq: summed squared intensity of recruited signal peaks.
                double repChargeNoisePower = 0.0;
                double repChargeSumSignalSq = 0.0;

                for (int z = candidate.MinAbsCharge; z <= candidate.MaxAbsCharge; z++)
                {
                    double isoDelta = Constants.C13MinusC12 / z;
                    double baseMz = candidate.Mass.ToMz(polSign * z);
                    bool isRepZ = (z == repAbsCharge);

                    // Search forward (n=0 is monoisotopic, n>0 are heavier isotopes)
                    for (int n = 0; n < nIso; n++)
                    {
                        double targetMz = baseMz + n * isoDelta;
                        double windowHalfDa = targetMz * p.DeconvolutionTolerancePpm * 1e-6;
                        double windowLo = targetMz - windowHalfDa;
                        double windowHi = targetMz + windowHalfDa;

                        var signalIndices = spectrum.GetPeakIndicesWithinTolerance(targetMz, tolerance);
                        if (signalIndices.Count == 0)
                        {
                            if (isRepZ)
                                CollectNoisePeaks(spectrum, windowLo, windowHi,
                                    signalIndices, ref repChargeNoisePower);
                            break;
                        }

                        int bestIdx = signalIndices.OrderByDescending(i => spectrum.YArray[i]).First();
                        double obsIntensity = spectrum.YArray[bestIdx];
                        double obsMz = spectrum.XArray[bestIdx];

                        int isoSlot = apexIdx + n;
                        if (isoSlot < isoWindowSize)
                        {
                            perIsotopeObs[isoSlot] += obsIntensity;
                            if (isRepZ) repChargeObs[isoSlot] += obsIntensity;
                        }

                        recruitedPeaks.Add((obsMz, obsIntensity, z, n));

                        if (isRepZ)
                        {
                            repChargeSumSignalSq += obsIntensity * obsIntensity;
                            CollectNoisePeaks(spectrum, windowLo, windowHi,
                                signalIndices, ref repChargeNoisePower);
                        }
                    }

                    // Search backward (lighter than monoisotopic, n=−1,−2,…)
                    for (int n = 1; n <= apexIdx; n++)
                    {
                        double targetMz = baseMz - n * isoDelta;
                        if (targetMz <= 0) break;

                        double windowHalfDa = targetMz * p.DeconvolutionTolerancePpm * 1e-6;
                        double windowLo = targetMz - windowHalfDa;
                        double windowHi = targetMz + windowHalfDa;

                        var signalIndices = spectrum.GetPeakIndicesWithinTolerance(targetMz, tolerance);
                        if (signalIndices.Count == 0)
                        {
                            if (isRepZ)
                                CollectNoisePeaks(spectrum, windowLo, windowHi,
                                    signalIndices, ref repChargeNoisePower);
                            break;
                        }

                        int bestIdx = signalIndices.OrderByDescending(i => spectrum.YArray[i]).First();
                        double obsIntensity = spectrum.YArray[bestIdx];
                        double obsMz = spectrum.XArray[bestIdx];

                        int isoSlot = apexIdx - n;
                        if (isoSlot >= 0 && isoSlot < isoWindowSize)
                        {
                            perIsotopeObs[isoSlot] += obsIntensity;
                            if (isRepZ) repChargeObs[isoSlot] += obsIntensity;
                        }

                        recruitedPeaks.Add((obsMz, obsIntensity, z, -n));

                        if (isRepZ)
                        {
                            repChargeSumSignalSq += obsIntensity * obsIntensity;
                            CollectNoisePeaks(spectrum, windowLo, windowHi,
                                signalIndices, ref repChargeNoisePower);
                        }
                    }
                }

                if (recruitedPeaks.Count < p.MinIsotopicPeakCount)
                    continue;

                // ── Step 4a: find best isotope offset ─────────────────────────
                double bestCosine = -1.0;
                int bestOffset = 0;

                for (int offset = -apexIdx; offset <= apexIdx; offset++)
                {
                    var obsSlice = new double[nIso];
                    for (int i = 0; i < nIso; i++)
                    {
                        int slot = i + offset;
                        if (slot >= 0 && slot < isoWindowSize)
                            obsSlice[i] = perIsotopeObs[slot];
                    }

                    double cos = SpectralSimilarity.CosineOfAlignedVectors(obsSlice, avgIntens);
                    if (cos > bestCosine)
                    {
                        bestCosine = cos;
                        bestOffset = offset;
                    }
                }

                // ── Step 4b: per-charge cosine for representative charge ────────
                // Apply the same offset alignment to the representative-charge-only
                // isotope vector to get chargeIsotopeCosine (used in f[1]).
                var repObsSlice = new double[nIso];
                for (int i = 0; i < nIso; i++)
                {
                    int slot = i + bestOffset;
                    if (slot >= 0 && slot < isoWindowSize)
                        repObsSlice[i] = repChargeObs[slot];
                }
                double repChargeCosine = SpectralSimilarity.CosineOfAlignedVectors(repObsSlice, avgIntens);

                // ── Step 4c: apply filters ─────────────────────────────────────
                if (bestCosine < p.MinCosineScore) continue;

                // Deduplicate: same raw spectrum peak recruited from multiple charges.
                var uniquePeaks = recruitedPeaks
                    .GroupBy(pk => spectrum.GetClosestPeakIndex(pk.mz))
                    .Select(g => g.OrderByDescending(pk => pk.intensity).First())
                    .ToList();

                if (uniquePeaks.Count < p.MinIsotopicPeakCount) continue;

                // monoMass = candidate.Mass − apexDaFromMono + bestOffset × C13MinusC12
                double monoMass = candidate.Mass - apexDaFromMono + bestOffset * Constants.C13MinusC12;
                if (monoMass < p.MinMassRange || monoMass > p.MaxMassRange) continue;

                int signedCharge = polSign * repAbsCharge;
                double totalIntensity = uniquePeaks.Sum(pk => pk.intensity);
                var outputPeaks = uniquePeaks.Select(pk => (pk.mz, pk.intensity)).ToList();

                // ── Compute avg ppm error (isotope index n still known) ────────
                double totalPpmError = 0.0;
                foreach (var (obsMz, _, z, n) in uniquePeaks)
                {
                    double isoDelta = Constants.C13MinusC12 / z;
                    double theorMz = candidate.Mass.ToMz(polSign * z) + n * isoDelta;
                    totalPpmError += Math.Abs(obsMz - theorMz) / theorMz * 1e6;
                }
                double avgPpmError = uniquePeaks.Count > 0
                    ? totalPpmError / uniquePeaks.Count
                    : 0.0;

                var envelope = new IsotopicEnvelope(
                    id: 0,
                    peaks: outputPeaks,
                    monoisotopicmass: monoMass,
                    chargestate: signedCharge,
                    intensity: totalIntensity,
                    score: bestCosine);   // Score = global cosine; replaced by Qscore after scoring

                results.Add(new MetaFlashDeconScorer.EnvelopeScoringData(
                    envelope: envelope,
                    avgPpmError: avgPpmError,
                    repChargeIsotopeCosine: repChargeCosine,
                    repChargeNoisePower: repChargeNoisePower,
                    repChargeSumSignalSquared: repChargeSumSignalSq));
            }

            return results;
        }

        // ── Noise collection helper ───────────────────────────────────────────

        /// <summary>
        /// Adds to <paramref name="noisePower"/> the summed squared intensity of
        /// all raw spectrum peaks within [windowLo, windowHi] that are NOT in
        /// <paramref name="signalIndices"/> (i.e. peaks in the isotope window but
        /// outside the ±ppm tolerance band around the expected isotope position).
        /// </summary>
        private static void CollectNoisePeaks(
            MzSpectrum spectrum,
            double windowLo,
            double windowHi,
            IReadOnlyCollection<int> signalIndices,
            ref double noisePower)
        {
            // Binary-search to the first index at or above windowLo.
            // MzSpectrum.XArray is sorted ascending.
            int lo = spectrum.GetClosestPeakIndex(windowLo);
            // Walk forward through all peaks in the window.
            for (int i = lo; i < spectrum.Size; i++)
            {
                double mz = spectrum.XArray[i];
                if (mz > windowHi) break;
                if (mz < windowLo) continue;
                if (signalIndices.Contains(i)) continue; // signal peak, skip
                double intensity = spectrum.YArray[i];
                noisePower += intensity * intensity;
            }
        }

        // ══════════════════════════════════════════════════════════════════════
        // STEP 2 — CANDIDATE MASS FINDING (unchanged)
        // ══════════════════════════════════════════════════════════════════════

        internal sealed class CandidateMass
        {
            public double Mass { get; init; }
            public double LogMass { get; init; }
            public int MinAbsCharge { get; init; }
            public int MaxAbsCharge { get; init; }
            public float SupportIntensity { get; init; }
        }

        internal static List<CandidateMass> FindCandidateMasses(
            List<LogMzPeak> logPeaks,
            double[] universalPattern,
            double[][] harmonicPatterns,
            double binMulFactor,
            MetaFlashDeconParameters deconParams)
        {
            if (logPeaks.Count == 0) return new List<CandidateMass>();

            int chargeRange = universalPattern.Length;
            int minAbsCharge = Math.Abs(deconParams.MinAssumedChargeState);

            double mzBinMinValue = logPeaks[0].LogMz;
            double mzBinMaxValue = logPeaks[^1].LogMz;
            double massBinMinValue = Math.Log(Math.Max(1.0, deconParams.MinMassRange - 50.0));
            double massBinMaxValue = Math.Log(deconParams.MaxMassRange + 200.0);

            int mzBinCount = BinNumber(mzBinMaxValue, mzBinMinValue, binMulFactor) + 1;
            int massBinCount = BinNumber(massBinMaxValue, massBinMinValue, binMulFactor) + 1;

            // Pre-compute integer harmonic bin offsets in mz-bin space.
            // The harmonic check compares mz-bin positions of the current peak against
            // positions that would be occupied by harmonic species at the same mass bin.
            var binHarmonicPatterns = new int[harmonicPatterns.Length][];
            for (int k = 0; k < harmonicPatterns.Length; k++)
            {
                binHarmonicPatterns[k] = new int[chargeRange];
                for (int j = 0; j < chargeRange; j++)
                    binHarmonicPatterns[k][j] = (int)Math.Round(
                        (mzBinMinValue - harmonicPatterns[k][j] - massBinMinValue) * binMulFactor);
            }

            // mz-bin occupancy map — used by the harmonic check
            var mzBinOccupied = new bool[mzBinCount];
            var mzBinIntensity = new float[mzBinCount];
            var peakBinIndex = new int[logPeaks.Count];

            for (int i = 0; i < logPeaks.Count; i++)
            {
                int b = BinNumber(logPeaks[i].LogMz, mzBinMinValue, binMulFactor);
                if (b < 0 || b >= mzBinCount) continue;
                peakBinIndex[i] = b;
                mzBinOccupied[b] = true;
                mzBinIntensity[b] += (float)logPeaks[i].Intensity;
            }

            var massBinOccupied = new bool[massBinCount];
            var massBinIntensity = new float[massBinCount];
            var massBinSupportCount = new int[massBinCount];
            var massBinMinChargeIdx = new int[massBinCount];
            var massBinMaxChargeIdx = new int[massBinCount];
            for (int i = 0; i < massBinCount; i++)
            {
                massBinMinChargeIdx[i] = int.MaxValue;
                massBinMaxChargeIdx[i] = int.MinValue;
            }

            var prevChargeIdx = new int[massBinCount];
            var prevIntensity = new float[massBinCount];
            for (int i = 0; i < massBinCount; i++) prevChargeIdx[i] = chargeRange + 2;

            // Main scan: right-to-left (high logMz → low logMz).
            // Scanning in this direction means charge index j increases as we move
            // to higher logMz, so consecutive charges for the same mass appear as
            // j, j-1 — i.e. prevChargeIdx[massBin] - j == -1 when continuous.
            for (int i = logPeaks.Count - 1; i >= 0; i--)
            {
                int mzBin = peakBinIndex[i];
                double logMz = logPeaks[i].LogMz;
                float intensity = (float)logPeaks[i].Intensity;

                for (int j = 0; j < chargeRange; j++)
                {
                    int massBin = BinNumber(logMz - universalPattern[j], massBinMinValue, binMulFactor);
                    if (massBin < 0 || massBin >= massBinCount) continue;

                    bool chargeIsContinuous = (prevChargeIdx[massBin] - j) == -1
                                             && prevChargeIdx[massBin] <= chargeRange;
                    float ratio = prevIntensity[massBin] <= 0f
                        ? MaxChargeIntensityRatio + 1f
                        : intensity / prevIntensity[massBin];
                    if (ratio < 1f) ratio = 1f / ratio;

                    bool passCheck = chargeIsContinuous && ratio <= MaxChargeIntensityRatio;

                    if (!passCheck)
                    {
                        massBinSupportCount[massBin] = 0;
                    }
                    else
                    {
                        bool isHarmonic = false;
                        float maxHarmonicIntensity = 0f;

                        for (int k = 0; k < binHarmonicPatterns.Length && !isHarmonic; k++)
                        {
                            int hBin = mzBin - binHarmonicPatterns[k][j];
                            if (hBin > 0 && hBin < mzBinCount && hBin != mzBin
                                && mzBinOccupied[hBin])
                            {
                                float hIntensity = mzBinIntensity[hBin];
                                float hRatio = hIntensity > intensity
                                    ? hIntensity / intensity
                                    : intensity / hIntensity;
                                if (hRatio <= MaxChargeIntensityRatio)
                                {
                                    isHarmonic = true;
                                    maxHarmonicIntensity = Math.Max(maxHarmonicIntensity, hIntensity);
                                }
                            }
                        }

                        if (!isHarmonic)
                        {
                            if (massBinSupportCount[massBin] == 0)
                            {
                                // Retroactive credit: this is the second confirming peak,
                                // so we count both this and the previous peak → support = 2.
                                massBinSupportCount[massBin] = 2;
                                massBinIntensity[massBin] += prevIntensity[massBin] + intensity;
                                int prevJ = prevChargeIdx[massBin];
                                if (prevJ < chargeRange)
                                {
                                    if (prevJ < massBinMinChargeIdx[massBin]) massBinMinChargeIdx[massBin] = prevJ;
                                    if (prevJ > massBinMaxChargeIdx[massBin]) massBinMaxChargeIdx[massBin] = prevJ;
                                }
                            }
                            else
                            {
                                massBinSupportCount[massBin]++;
                                massBinIntensity[massBin] += intensity;
                            }

                            if (!massBinOccupied[massBin]
                                && massBinSupportCount[massBin] >= MinSupportPeakCount)
                                massBinOccupied[massBin] = true;

                            if (j < massBinMinChargeIdx[massBin]) massBinMinChargeIdx[massBin] = j;
                            if (j > massBinMaxChargeIdx[massBin]) massBinMaxChargeIdx[massBin] = j;
                        }
                        else
                        {
                            massBinIntensity[massBin] -= maxHarmonicIntensity;
                            if (massBinSupportCount[massBin] > 0) massBinSupportCount[massBin]--;
                        }
                    }

                    prevIntensity[massBin] = intensity;
                    prevChargeIdx[massBin] = j;
                }
            }

            // Conflict resolution — faithful to OpenMS FLASHDeconv filterMassBins_
            // (FLASHDeconvAlgorithm.cpp:543, select_top_N = 3): each m/z peak matches one
            // candidate mass per charge, so it is allowed to "claim" only its TOP-3
            // most-intense occupied mass bins (top-N=3 "to consider frequent coelution"),
            // not just its single best. This bounds how many masses one peak can seed.
            const int SelectTopN = 3;
            var finalMassBins = new bool[massBinCount];
            Span<long> topBins = stackalloc long[SelectTopN];
            Span<float> topInts = stackalloc float[SelectTopN];
            for (int i = 0; i < logPeaks.Count; i++)
            {
                double logMz = logPeaks[i].LogMz;
                for (int t = 0; t < SelectTopN; t++) { topBins[t] = -1; topInts[t] = -1f; }

                for (int j = 0; j < chargeRange; j++)
                {
                    int massBin = BinNumber(logMz - universalPattern[j], massBinMinValue, binMulFactor);
                    if (massBin < 0 || massBin >= massBinCount) continue;
                    if (!massBinOccupied[massBin]) continue;

                    float mbi = massBinIntensity[massBin];
                    if (mbi <= topInts[SelectTopN - 1]) continue;

                    // Insert into the descending-by-intensity top-N (rolling shift).
                    int pos = SelectTopN - 1;
                    while (pos > 0 && topInts[pos - 1] < mbi)
                    {
                        topInts[pos] = topInts[pos - 1];
                        topBins[pos] = topBins[pos - 1];
                        pos--;
                    }
                    topInts[pos] = mbi;
                    topBins[pos] = massBin;
                }

                for (int t = 0; t < SelectTopN; t++)
                    if (topBins[t] >= 0) finalMassBins[topBins[t]] = true;
            }

            var candidates = new List<CandidateMass>();
            for (int b = 0; b < massBinCount; b++)
            {
                if (!finalMassBins[b]) continue;
                if (massBinMinChargeIdx[b] == int.MaxValue || massBinMaxChargeIdx[b] == int.MinValue) continue;

                double logMass = BinValue(b, massBinMinValue, binMulFactor);
                double mass = Math.Exp(logMass);

                if (mass < deconParams.MinMassRange || mass > deconParams.MaxMassRange) continue;

                candidates.Add(new CandidateMass
                {
                    Mass = mass,
                    LogMass = logMass,
                    MinAbsCharge = minAbsCharge + massBinMinChargeIdx[b],
                    MaxAbsCharge = minAbsCharge + massBinMaxChargeIdx[b],
                    SupportIntensity = massBinIntensity[b]
                });
            }

            return candidates;
        }

        // ══════════════════════════════════════════════════════════════════════
        // BIN HELPERS
        // ══════════════════════════════════════════════════════════════════════

        internal static int BinNumber(double value, double minValue, double binMulFactor)
        {
            double raw = (value - minValue) * binMulFactor;
            return raw < 0.0 ? 0 : (int)Math.Round(raw);
        }

        internal static double BinValue(int bin, double minValue, double binMulFactor)
            => minValue + bin / binMulFactor;

        // ══════════════════════════════════════════════════════════════════════
        // STEP 1 HELPERS
        // ══════════════════════════════════════════════════════════════════════

        internal static double GetLogMz(double mz, Polarity polarity)
        {
            double uncharged = mz - Math.Sign((int)polarity) * Constants.ProtonMass;
            return uncharged > 0.0 ? Math.Log(uncharged) : double.NegativeInfinity;
        }

        internal static List<LogMzPeak> BuildLogMzPeaks(
            MzSpectrum spectrum, MzRange range, Polarity polarity)
        {
            var peaks = new List<LogMzPeak>(spectrum.Size);
            for (int i = 0; i < spectrum.Size; i++)
            {
                double mz = spectrum.XArray[i];
                double intensity = spectrum.YArray[i];
                if (mz < range.Minimum || mz > range.Maximum) continue;
                if (intensity <= 0.0) continue;
                double logMz = GetLogMz(mz, polarity);
                if (!double.IsFinite(logMz)) continue;
                peaks.Add(new LogMzPeak(mz, intensity, logMz));
            }
            peaks.Sort((a, b) => a.LogMz.CompareTo(b.LogMz));
            return peaks;
        }

        internal static double[] BuildUniversalPattern(int minAbsCharge, int chargeRange)
        {
            if (chargeRange <= 0) throw new ArgumentOutOfRangeException(nameof(chargeRange), "Must be > 0");
            if (minAbsCharge <= 0) throw new ArgumentOutOfRangeException(nameof(minAbsCharge), "Must be > 0");
            var p = new double[chargeRange];
            for (int j = 0; j < chargeRange; j++) p[j] = -Math.Log(minAbsCharge + j);
            return p;
        }

        internal static double[][] BuildHarmonicPatterns(int minAbsCharge, int chargeRange)
        {
            double[] up = BuildUniversalPattern(minAbsCharge, chargeRange);
            var h = new double[HarmonicDenominators.Length][];
            for (int k = 0; k < HarmonicDenominators.Length; k++)
            {
                int d = HarmonicDenominators[k];
                int n = d / 2;
                h[k] = new double[chargeRange];
                for (int j = 0; j < chargeRange; j++)
                {
                    double b = Math.Exp(-up[j]);
                    double a = j > 0 ? Math.Exp(-up[j - 1]) : 0.0;
                    double hMz = b - (b - a) * (double)n / d;
                    h[k][j] = hMz > 0.0 ? -Math.Log(hMz) : up[j];
                }
            }
            return h;
        }

        // ══════════════════════════════════════════════════════════════════════
        // LogMzPeak STRUCT
        // ══════════════════════════════════════════════════════════════════════

        internal readonly struct LogMzPeak
        {
            public readonly double Mz;
            public readonly double Intensity;
            public readonly double LogMz;
            public readonly int AbsCharge;
            public readonly int IsotopeIndex;

            internal LogMzPeak(double mz, double intensity, double logMz,
                               int absCharge = 0, int isotopeIndex = 0)
            {
                Mz = mz;
                Intensity = intensity;
                LogMz = logMz;
                AbsCharge = absCharge;
                IsotopeIndex = isotopeIndex;
            }

            public LogMzPeak WithAssignment(int absCharge, int isotopeIndex)
                => new LogMzPeak(Mz, Intensity, LogMz, absCharge, isotopeIndex);

            public override string ToString()
                => $"mz={Mz:F4} logMz={LogMz:F6} I={Intensity:G4} z={AbsCharge} iso={IsotopeIndex}";
        }
    }
}