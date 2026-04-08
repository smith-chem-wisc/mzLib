// FLASHDeconvolutionAlgorithm.cs
//
// Status: STEP 2 COMPLETE — candidate mass finding
//
// Step 1 (log transform + pattern construction): complete
// Step 2 (binned universal pattern scanning → candidate masses): complete
// Steps 3–5 (isotope clustering, scoring, output): TODO skeleton remains
//
// Reference:
//   Jeong et al. (2020) Cell Systems 10(2):213-218.e6
//   DOI: 10.1016/j.cels.2020.01.003
//   OpenMS source: src/openms/source/ANALYSIS/TOPDOWN/SpectralDeconvolution.cpp
//     → generatePeakGroupsFromSpectrum_(), updateMassBins_(),
//       updateCandidateMassBins_(), filterMassBins_(), getCandidatePeakGroups_()

using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;                           // Constants.ProtonMass, Constants.C13MinusC12
using MassSpectrometry.MzSpectra;          // SpectralSimilarity.CosineOfAlignedVectors
using MzLibUtil;

namespace MassSpectrometry
{
    internal class FLASHDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        // ── Harmonic denominators (½, ⅓, ⅖ from the paper) ──────────────────
        private static readonly int[] HarmonicDenominators = { 2, 3, 5 };

        // ── Minimum number of continuous support peaks required to mark a
        //    mass bin as a candidate (mirrors min_support_peak_count_ in OpenMS,
        //    which defaults to 3 for MS1 and 2 for MS2). Fixed at 3 here; could
        //    be promoted to FLASHDeconvolutionParameters later.
        private const int MinSupportPeakCount = 3;

        // ── Intensity-ratio threshold between consecutive charge states.
        //    If the ratio exceeds this factor the charge run is considered broken
        //    and the support count resets. OpenMS uses 10 for charges ≤ 3 and a
        //    sliding value for higher charges; we use the conservative fixed value.
        private const float MaxChargeIntensityRatio = 10.0f;

        // ── Bin multiplication factor: bins per log unit.
        //    At tolerance T ppm the bin width in log space is approximately T×1e-6,
        //    so bin_mul_factor ≈ 1/(T×1e-6). We compute it from the parameter.
        private static double ComputeBinMulFactor(double tolerancePpm)
            => 1.0 / (tolerancePpm * 1e-6);

        internal FLASHDeconvolutionAlgorithm(DeconvolutionParameters deconParameters)
            : base(deconParameters)
        {
        }

        // ══════════════════════════════════════════════════════════════════════
        // PUBLIC ENTRY POINT
        // ══════════════════════════════════════════════════════════════════════

        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var deconParams = DeconvolutionParameters as FLASHDeconvolutionParameters
                ?? throw new MzLibException("Deconvolution params and algorithm do not match");

            range ??= spectrum.Range;

            if (spectrum.Size == 0)
                return Enumerable.Empty<IsotopicEnvelope>();

            // ── Step 1: log transform + patterns ─────────────────────────────
            var logPeaks = BuildLogMzPeaks(spectrum, range, deconParams.Polarity);
            if (logPeaks.Count == 0)
                return Enumerable.Empty<IsotopicEnvelope>();

            int minAbsCharge = Math.Abs(deconParams.MinAssumedChargeState);
            int maxAbsCharge = Math.Abs(deconParams.MaxAssumedChargeState);
            int chargeRange = maxAbsCharge - minAbsCharge + 1;

            double[] universalPattern = BuildUniversalPattern(minAbsCharge, chargeRange);
            double[][] harmonicPatterns = BuildHarmonicPatterns(minAbsCharge, chargeRange);

            // ── Step 2: candidate mass finding ────────────────────────────────
            double binMulFactor = ComputeBinMulFactor(deconParams.DeconvolutionTolerancePpm);

            var candidateMasses = FindCandidateMasses(
                logPeaks,
                universalPattern,
                harmonicPatterns,
                binMulFactor,
                deconParams);

            if (candidateMasses.Count == 0)
                return Enumerable.Empty<IsotopicEnvelope>();

            // ── Steps 3–5: isotope clustering, scoring, output ────────────────
            //
            // STEP 3 — ISOTOPE ENVELOPE CLUSTERING (to be implemented)
            //   For each CandidateMass, scan the spectrum around the expected m/z
            //   positions for each charge in its charge range. For each charge c,
            //   the most intense peak within mz_delta of (mass + n * C13MinusC12)/c
            //   for n = 0, ±1, ±2, … is recruited. Isotope index is assigned by
            //   round((mz_observed - mz_base) / (C13MinusC12 / c)), then shifted
            //   so that the apex of the Averagine distribution sits at the right index.
            //   Peaks outside [−leftCount, rightCount] from the Averagine apex are excluded.
            //   Use AverageResidueModel.GetAllTheoreticalMasses/Intensities for the lookup.
            //
            // STEP 4 — SCORING AND FILTRATION (to be implemented)
            //   a) Build per-isotope intensity vector from all recruited peaks.
            //   b) Try isotope offsets from −apexIndex to +apexIndex; pick the offset
            //      that maximises SpectralSimilarity.CosineOfAlignedVectors(observed, theoretical).
            //   c) Reject if best cosine < deconParams.MinCosineScore.
            //   d) Determine monoisotopic mass from the winning offset.
            //   e) Reject if monoisotopic mass outside [MinMassRange, MaxMassRange].
            //   f) Reject if recruited peak count < MinIsotopicPeakCount.
            //
            // STEP 5 — OUTPUT (to be implemented)
            //   new IsotopicEnvelope(id: 0, peaks, monoMass, signedCharge, totalIntensity, cosineScore)
            //   Sign the charge: positive for Polarity.Positive, negative for Polarity.Negative.

            return Enumerable.Empty<IsotopicEnvelope>();
        }

        // ══════════════════════════════════════════════════════════════════════
        // STEP 2 — CANDIDATE MASS FINDING
        // ══════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Intermediate result returned by <see cref="FindCandidateMasses"/>.
        /// Represents a candidate mass identified by universal-pattern scanning,
        /// together with the charge range that contributed to it.
        /// Fields are populated incrementally and are readable by Step 3.
        /// </summary>
        internal sealed class CandidateMass
        {
            /// <summary>Candidate neutral monoisotopic mass (exp of the mass-bin log value).</summary>
            public double Mass { get; init; }

            /// <summary>Log of the candidate mass (the mass-bin centre value).</summary>
            public double LogMass { get; init; }

            /// <summary>
            /// Minimum absolute charge state that contributed support peaks
            /// for this candidate (index into the universal pattern + minAbsCharge).
            /// </summary>
            public int MinAbsCharge { get; init; }

            /// <summary>Maximum absolute charge state that contributed support peaks.</summary>
            public int MaxAbsCharge { get; init; }

            /// <summary>
            /// Summed intensity of all non-harmonic support peaks for this candidate.
            /// Used to resolve conflicts when a peak belongs to multiple candidates.
            /// </summary>
            public float SupportIntensity { get; init; }
        }

        /// <summary>
        /// Scans the log-transformed spectrum using the universal charge pattern and
        /// returns a list of candidate masses sorted ascending by mass.
        /// <para>
        /// This is the C# equivalent of OpenMS
        /// <c>updateCandidateMassBins_() + filterMassBins_()</c>. The binning
        /// approach is identical: log-m/z and log-mass values are mapped to integer
        /// bin indices so that the universal-pattern offset becomes a constant integer
        /// offset per charge, enabling O(peaks × charges) scanning.
        /// </para>
        /// </summary>
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

            // ── Set up bin space ──────────────────────────────────────────────
            // mz bin space: logMz values of all peaks → integer bins
            double mzBinMinValue = logPeaks[0].LogMz;       // smallest logMz in spectrum
            double mzBinMaxValue = logPeaks[^1].LogMz;      // largest  logMz in spectrum

            // mass bin space: spans from log(MinMassRange) to
            // log(MaxMassRange + some isotope headroom).  We add headroom so that
            // the most abundant isotope of heavy proteins is still captured.
            double massBinMinValue = Math.Log(Math.Max(1.0, deconParams.MinMassRange - 50.0));
            double massBinMaxValue = Math.Log(deconParams.MaxMassRange + 200.0);

            int mzBinCount = BinNumber(mzBinMaxValue, mzBinMinValue, binMulFactor) + 1;
            int massBinCount = BinNumber(massBinMaxValue, massBinMinValue, binMulFactor) + 1;

            // ── Pre-compute integer universal-pattern offsets ─────────────────
            // Each element j: the integer offset from an mz-bin to the mass-bin
            // corresponding to charge (minAbsCharge + j).
            // offset[j] = round( (mzBinMin − universalPattern[j] − massBinMin) × binMulFactor )
            // This converts "shift by −log(c) in log space" into an integer bin shift.
            var binUniversalPattern = new int[chargeRange];
            for (int j = 0; j < chargeRange; j++)
            {
                binUniversalPattern[j] = (int)Math.Round(
                    (mzBinMinValue - universalPattern[j] - massBinMinValue) * binMulFactor);
            }

            // Harmonic pattern bin offsets [harmonicDenominatorIndex][chargeIndex]
            var binHarmonicPatterns = new int[harmonicPatterns.Length][];
            for (int k = 0; k < harmonicPatterns.Length; k++)
            {
                binHarmonicPatterns[k] = new int[chargeRange];
                for (int j = 0; j < chargeRange; j++)
                {
                    binHarmonicPatterns[k][j] = (int)Math.Round(
                        (mzBinMinValue - harmonicPatterns[k][j] - massBinMinValue) * binMulFactor);
                }
            }

            // ── Build mz-bin arrays from log peaks ────────────────────────────
            // mzBinOccupied[b] = true  when any peak maps to bin b
            // mzBinIntensity[b] = intensity of the peak at bin b (summed if multiple)
            var mzBinOccupied = new bool[mzBinCount];
            var mzBinIntensity = new float[mzBinCount];
            var peakBinIndex = new int[logPeaks.Count];   // bin index per peak

            for (int i = 0; i < logPeaks.Count; i++)
            {
                int b = BinNumber(logPeaks[i].LogMz, mzBinMinValue, binMulFactor);
                if (b < 0 || b >= mzBinCount) continue;
                peakBinIndex[i] = b;
                mzBinOccupied[b] = true;
                mzBinIntensity[b] += (float)logPeaks[i].Intensity;
            }

            // ── Per-mass-bin accumulators ─────────────────────────────────────
            var massBinOccupied = new bool[massBinCount];  // candidate bit
            var massBinIntensity = new float[massBinCount]; // support intensity
            var massBinSupportCount = new int[massBinCount];   // # support peaks
            var massBinMinChargeIdx = new int[massBinCount];   // lowest charge index j
            var massBinMaxChargeIdx = new int[massBinCount];   // highest charge index j

            // Initialise charge range tracking
            for (int i = 0; i < massBinCount; i++)
            {
                massBinMinChargeIdx[i] = int.MaxValue;
                massBinMaxChargeIdx[i] = int.MinValue;
            }

            // Per-mass-bin tracking of previous charge and intensity
            // (for continuity / intensity-ratio checks)
            var prevChargeIdx = new int[massBinCount];
            var prevIntensity = new float[massBinCount];
            for (int i = 0; i < massBinCount; i++) prevChargeIdx[i] = chargeRange + 2; // sentinel: "none"

            // ── Main scan: iterate peaks from HIGH to LOW mz (right to left) ──
            // OpenMS iterates in reverse so that charge states increment (small
            // to large) as we walk from high-mz to low-mz for a given mass.
            // This makes the "continuous charge" check meaningful.
            for (int i = logPeaks.Count - 1; i >= 0; i--)
            {
                int mzBin = peakBinIndex[i];
                double logMz = logPeaks[i].LogMz;
                float intensity = (float)logPeaks[i].Intensity;

                // For each charge state j (low abs charge to high abs charge):
                for (int j = 0; j < chargeRange; j++)
                {
                    // ── SINGLE-ROUND mass bin computation ─────────────────────
                    // Computing massBin as mzBin + binUniversalPattern[j] uses
                    // TWO separate round() calls (one for mzBin, one for the
                    // pre-computed integer offset), which can place peaks from
                    // different charge states of the same mass into adjacent bins
                    // (off by ±1). Instead, we compute the mass bin in one
                    // combined round(), which is mathematically equivalent to
                    // round(log(mass) × binMulFactor) for all charges:
                    //   logMz - universalPattern[j] = log(mass/z) - (-log(z)) = log(mass)
                    int massBin = BinNumber(logMz - universalPattern[j], massBinMinValue, binMulFactor);
                    if (massBin < 0 || massBin >= massBinCount) continue;

                    int absCharge = minAbsCharge + j;

                    // ── Intensity-ratio / continuity check ────────────────────
                    bool chargeIsContinuous = (prevChargeIdx[massBin] - j) == -1
                                             && prevChargeIdx[massBin] <= chargeRange;
                    float ratio = prevIntensity[massBin] <= 0f
                        ? MaxChargeIntensityRatio + 1f
                        : intensity / prevIntensity[massBin];
                    if (ratio < 1f) ratio = 1f / ratio;

                    bool passIntensityCheck = chargeIsContinuous && ratio <= MaxChargeIntensityRatio;

                    if (!passIntensityCheck)
                    {
                        // Broken charge run: reset support count.
                        massBinSupportCount[massBin] = 0;
                    }
                    else
                    {
                        // ── Harmonic check ────────────────────────────────────
                        // A peak is "harmonic" if an interleaved harmonic-pattern
                        // peak is present in the occupied mz-bin set at the same
                        // intensity order of magnitude.
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
                            // When the second peak of a new continuous run arrives it
                            // confirms the first peak retroactively — credit both.
                            // (The first peak of a run cannot credit itself because it
                            // has no prior peak to be "continuous" with; OpenMS handles
                            // this via support_peak_intensity carrying the previous
                            // peak's intensity forward.)
                            if (massBinSupportCount[massBin] == 0)
                            {
                                // This is the second peak confirming the first:
                                // count the predecessor too.
                                massBinSupportCount[massBin] = 2;
                                massBinIntensity[massBin] += prevIntensity[massBin] + intensity;
                                // Also credit the first peak's charge index in the range
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

                            // Mark as candidate once MinSupportPeakCount - 1 continuous
                            // support peaks have been confirmed. The -1 accounts for the
                            // first peak of the run being credited retroactively above:
                            // after the second peak arrives, spc == 2 and we have already
                            // counted the full run of 2 confirmed peaks (the retroactive
                            // first + the current second). With MinSupportPeakCount = 3,
                            // we mark after the third peak sets spc = 3.
                            if (!massBinOccupied[massBin]
                                && massBinSupportCount[massBin] >= MinSupportPeakCount)
                            {
                                massBinOccupied[massBin] = true;
                            }

                            // Track charge range for this mass bin
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

            // ── Conflict resolution: for each mz-peak assign to highest-intensity mass ─
            // When a peak maps to multiple candidate mass bins (different charges),
            // keep only the assignment to the highest-intensity mass bin.
            // Uses the same single-round mass bin computation as the main scan loop.
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

            // ── Collect candidates ────────────────────────────────────────────
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

            return candidates; // already sorted ascending by mass (bins are ordered)
        }

        // ══════════════════════════════════════════════════════════════════════
        // BIN ARITHMETIC HELPERS
        // ══════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Converts a log-space value to an integer bin index.
        /// <c>bin = round((value − minValue) × binMulFactor)</c>
        /// </summary>
        internal static int BinNumber(double value, double minValue, double binMulFactor)
        {
            double raw = (value - minValue) * binMulFactor;
            return raw < 0.0 ? 0 : (int)Math.Round(raw);
        }

        /// <summary>
        /// Converts an integer bin index back to the corresponding log-space value.
        /// <c>value = minValue + bin / binMulFactor</c>
        /// </summary>
        internal static double BinValue(int bin, double minValue, double binMulFactor)
            => minValue + bin / binMulFactor;

        // ══════════════════════════════════════════════════════════════════════
        // STEP 1 METHODS (unchanged from Step 1 — reproduced for completeness)
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
            var pattern = new double[chargeRange];
            for (int j = 0; j < chargeRange; j++)
                pattern[j] = -Math.Log(minAbsCharge + j);
            return pattern;
        }

        internal static double[][] BuildHarmonicPatterns(int minAbsCharge, int chargeRange)
        {
            double[] up = BuildUniversalPattern(minAbsCharge, chargeRange);
            var harmonics = new double[HarmonicDenominators.Length][];
            for (int k = 0; k < HarmonicDenominators.Length; k++)
            {
                int d = HarmonicDenominators[k];
                int n = d / 2;
                harmonics[k] = new double[chargeRange];
                for (int j = 0; j < chargeRange; j++)
                {
                    double b = Math.Exp(-up[j]);
                    double a = j > 0 ? Math.Exp(-up[j - 1]) : 0.0;
                    double harmonicMz = b - (b - a) * (double)n / d;
                    harmonics[k][j] = harmonicMz > 0.0 ? -Math.Log(harmonicMz) : up[j];
                }
            }
            return harmonics;
        }

        // ══════════════════════════════════════════════════════════════════════
        // LogMzPeak STRUCT (unchanged from Step 1)
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