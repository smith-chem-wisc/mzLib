// FLASHDeconvolutionAlgorithm.cs
//
// Status: STEPS 1–5 COMPLETE — full spectral deconvolution implemented
//
// Step 1 — log transform + universal/harmonic pattern construction
// Step 2 — binned candidate mass finding (single-round mass bin computation)
// Step 3 — isotope peak recruitment per candidate per charge
// Step 4 — cosine scoring against Averagine; isotope offset optimisation
// Step 5 — IsotopicEnvelope output
//
// Reference:
//   Jeong et al. (2020) Cell Systems 10(2):213-218.e6
//   DOI: 10.1016/j.cels.2020.01.003
//   OpenMS source: src/openms/source/ANALYSIS/TOPDOWN/SpectralDeconvolution.cpp

using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;                           // Constants.ProtonMass, C13MinusC12, ToMz/ToMass
using MassSpectrometry.MzSpectra;          // SpectralSimilarity.CosineOfAlignedVectors
using MzLibUtil;                           // PpmTolerance

namespace MassSpectrometry
{
    internal class FLASHDeconvolutionAlgorithm : DeconvolutionAlgorithm
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

        private static double ComputeBinMulFactor(double tolerancePpm)
            => 1.0 / (tolerancePpm * 1e-6);

        // ── Constructor ───────────────────────────────────────────────────────
        internal FLASHDeconvolutionAlgorithm(DeconvolutionParameters deconParameters)
            : base(deconParameters)
        {
        }

        // ══════════════════════════════════════════════════════════════════════
        // ENTRY POINT
        // ══════════════════════════════════════════════════════════════════════

        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var p = DeconvolutionParameters as FLASHDeconvolutionParameters
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
            double binMulFactor = ComputeBinMulFactor(p.DeconvolutionTolerancePpm);

            // ── Step 2 ────────────────────────────────────────────────────────
            var candidates = FindCandidateMasses(
                logPeaks, universalPattern, harmonicPatterns, binMulFactor, p);
            if (candidates.Count == 0) return Enumerable.Empty<IsotopicEnvelope>();

            // ── Steps 3–5 ────────────────────────────────────────────────────
            double medianIntensity = FLASHDeconvScorer.ComputeMedianIntensity(spectrum);
            return FLASHDeconvScorer.AssignQscores(
                FLASHDeconvDeduplicator.Deduplicate(ScoreAndBuildEnvelopes(spectrum, candidates, p),
                                                    p.DeconvolutionTolerancePpm),
                medianIntensity,
                p.DeconvolutionTolerancePpm);
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
        /// envelopes.
        /// </summary>
        private IEnumerable<IsotopicEnvelope> ScoreAndBuildEnvelopes(
            MzSpectrum spectrum,
            List<CandidateMass> candidates,
            FLASHDeconvolutionParameters p)
        {
            var results = new List<IsotopicEnvelope>();
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
                // This is the distance in Da from monoisotopic to apex — exactly what we need
                // for the monoisotopic mass calculation, and avoids all fine-resolution array indexing.
                int avgIdx = p.AverageResidueModel.GetMostIntenseMassIndex(candidate.Mass);
                double apexDaFromMono = p.AverageResidueModel.GetDiffToMonoisotopic(avgIdx);

                double[] avgMasses = p.AverageResidueModel.GetAllTheoreticalMasses(avgIdx);
                double[] avgIntensRaw = p.AverageResidueModel.GetAllTheoreticalIntensities(avgIdx);

                // Sort to mass-ascending order for cosine scoring vector
                var sortedAvg = avgMasses.Zip(avgIntensRaw)
                    .OrderBy(pair => pair.First)
                    .ToArray();
                double[] avgIntens = sortedAvg.Select(pair => pair.Second).ToArray();
                int nIso = avgIntens.Length;

                // apexIdx in the mass-ascending array — used to set the isoWindowSize
                // and to center the perIsotopeObs vector.
                int apexIdx = Array.IndexOf(avgIntens, avgIntens.Max());

                // ── Step 3b: recruit peaks per charge state ────────────────────
                // Per-isotope intensity vector accumulated across all charges.
                // Slot 0 = monoisotopic (candidate.Mass), apexIdx = most-abundant isotope.
                // Track (mz, intensity, z, n) so we can compute monoMass correctly per peak.
                int isoWindowSize = nIso + apexIdx;
                var perIsotopeObs = new double[isoWindowSize];
                // recruitedPeaks: (mz, intensity, absCharge, isotopeIndex n from monoisotopic)
                var recruitedPeaks = new List<(double mz, double intensity, int z, int n)>();

                for (int z = candidate.MinAbsCharge; z <= candidate.MaxAbsCharge; z++)
                {
                    double isoDelta = Constants.C13MinusC12 / z;
                    double baseMz = candidate.Mass.ToMz(polSign * z);

                    // Search forward (n=0 is monoisotopic, n>0 are heavier)
                    for (int n = 0; n < nIso; n++)
                    {
                        double targetMz = baseMz + n * isoDelta;
                        var indices = spectrum.GetPeakIndicesWithinTolerance(targetMz, tolerance);
                        if (indices.Count == 0) break;

                        int bestIdx = indices.OrderByDescending(i => spectrum.YArray[i]).First();
                        double obsIntensity = spectrum.YArray[bestIdx];
                        double obsMz = spectrum.XArray[bestIdx];

                        int isoSlot = apexIdx + n;
                        if (isoSlot < isoWindowSize)
                            perIsotopeObs[isoSlot] += obsIntensity;

                        recruitedPeaks.Add((obsMz, obsIntensity, z, n));
                    }

                    // Search backward (lighter than monoisotopic, n=−1,−2,…)
                    for (int n = 1; n <= apexIdx; n++)
                    {
                        double targetMz = baseMz - n * isoDelta;
                        if (targetMz <= 0) break;

                        var indices = spectrum.GetPeakIndicesWithinTolerance(targetMz, tolerance);
                        if (indices.Count == 0) break;

                        int bestIdx = indices.OrderByDescending(i => spectrum.YArray[i]).First();
                        double obsIntensity = spectrum.YArray[bestIdx];
                        double obsMz = spectrum.XArray[bestIdx];

                        int isoSlot = apexIdx - n;
                        if (isoSlot >= 0 && isoSlot < isoWindowSize)
                            perIsotopeObs[isoSlot] += obsIntensity;

                        recruitedPeaks.Add((obsMz, obsIntensity, z, -n));
                    }
                }

                if (recruitedPeaks.Count < p.MinIsotopicPeakCount)
                    continue;

                // ── Step 4a: find best isotope offset ─────────────────────────
                // The candidate mass from Step 2 is approximate. Shift the observed
                // intensity vector relative to the Averagine template to find the
                // offset that maximises cosine similarity. Offset search range: ±apexIdx.
                double bestCosine = -1.0;
                int bestOffset = 0;

                for (int offset = -apexIdx; offset <= apexIdx; offset++)
                {
                    // Align: avgIntens[i] ↔ perIsotopeObs[i + apexIdx + offset]
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

                // ── Step 4b: apply filters ─────────────────────────────────────
                if (bestCosine < p.MinCosineScore) continue;

                // Deduplicate: same raw spectrum peak recruited from multiple charges.
                var uniquePeaks = recruitedPeaks
                    .GroupBy(pk => spectrum.GetClosestPeakIndex(pk.mz))
                    .Select(g => g.OrderByDescending(pk => pk.intensity).First())
                    .ToList();

                if (uniquePeaks.Count < p.MinIsotopicPeakCount) continue;

                // candidate.Mass ≈ apex mass (Step 2 bins on whatever peaks triggered the run,
                // which are near the apex of the isotope distribution).
                // The Averagine apex is apexDaFromMono Da above monoisotopic.
                // bestOffset shifts the alignment frame by that many fine-resolution bins,
                // which corresponds to bestOffset * C13MinusC12 in mass space.
                //
                // monoMass = candidate.Mass - apexDaFromMono + bestOffset * C13MinusC12
                //
                // This is exact regardless of charge and independent of the fine-resolution
                // array indexing, because apexDaFromMono is in Da directly from Averagine.
                double monoMass = candidate.Mass - apexDaFromMono + bestOffset * Constants.C13MinusC12;
                if (monoMass < p.MinMassRange || monoMass > p.MaxMassRange) continue;

                // Representative charge for the signed output
                int repCharge = (candidate.MinAbsCharge + candidate.MaxAbsCharge) / 2;
                int signedCharge = polSign * repCharge;

                double totalIntensity = uniquePeaks.Sum(pk => pk.intensity);
                var outputPeaks = uniquePeaks
                    .Select(pk => (pk.mz, pk.intensity))
                    .ToList();

                results.Add(new IsotopicEnvelope(
                    id: 0,
                    peaks: outputPeaks,
                    monoisotopicmass: monoMass,
                    chargestate: signedCharge,
                    intensity: totalIntensity,
                    score: bestCosine));
            }

            return results;
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
            FLASHDeconvolutionParameters deconParams)
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

            // Pre-compute integer harmonic bin offsets for mz-space lookups
            var binHarmonicPatterns = new int[harmonicPatterns.Length][];
            for (int k = 0; k < harmonicPatterns.Length; k++)
            {
                binHarmonicPatterns[k] = new int[chargeRange];
                for (int j = 0; j < chargeRange; j++)
                    binHarmonicPatterns[k][j] = (int)Math.Round(
                        (mzBinMinValue - harmonicPatterns[k][j] - massBinMinValue) * binMulFactor);
            }

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

            // Main scan: right-to-left (high logMz → low logMz)
            for (int i = logPeaks.Count - 1; i >= 0; i--)
            {
                int mzBin = peakBinIndex[i];
                double logMz = logPeaks[i].LogMz;
                float intensity = (float)logPeaks[i].Intensity;

                for (int j = 0; j < chargeRange; j++)
                {
                    // Single-round mass bin: avoids double-rounding that scatters
                    // peaks from the same mass into adjacent bins.
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

            // Conflict resolution: assign each mz-peak to the highest-intensity mass bin
            var finalMassBins = new bool[massBinCount];
            for (int i = 0; i < logPeaks.Count; i++)
            {
                double logMz = logPeaks[i].LogMz;
                long bestMassBin = -1;
                float bestIntensity = -1f;

                for (int j = 0; j < chargeRange; j++)
                {
                    int massBin = BinNumber(logMz - universalPattern[j], massBinMinValue, binMulFactor);
                    if (massBin < 0 || massBin >= massBinCount) continue;
                    if (!massBinOccupied[massBin]) continue;
                    if (massBinIntensity[massBin] > bestIntensity)
                    {
                        bestIntensity = massBinIntensity[massBin];
                        bestMassBin = massBin;
                    }
                }

                if (bestMassBin >= 0) finalMassBins[bestMassBin] = true;
            }

            var candidates = new List<CandidateMass>();
            for (int b = 0; b < massBinCount; b++)
            {
                if (!finalMassBins[b]) continue;
                if (massBinMinChargeIdx[b] > massBinMaxChargeIdx[b]) continue;

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