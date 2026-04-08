// FLASHDeconvScorer.cs
//
// Computes a Qscore for IsotopicEnvelope objects produced by
// FLASHDeconvolutionAlgorithm. The model ports the logistic regression in
// OpenMS PeakGroupScoring.cpp (Kyowon Jeong, BSD-3-Clause).
//
// ── OpenMS feature vector (toFeatureVector_) ─────────────────────────────────
//   f[0] = isotopeCosine
//   f[1] = isotopeCosine − chargeIsotopeCosine(repZ)   inter-charge cosine drop
//   f[2] = log2(1 + chargeSNR(repZ))                   per-charge SNR, log-scaled
//   f[3] = log2(1 + chargeSNR(repZ)) − log2(1 + SNR)   per-charge vs global SNR drop
//   f[4] = avgPPMError
//
// ── OpenMS SNR formula (updateSNR_ in PeakGroup.cpp) ─────────────────────────
// For each charge c:
//   cos²      = per_charge_cos[c]²
//   sig_pwr   = per_charge_sum_signal_squared[c] × cos²
//   noise_pwr = per_charge_noise_pwr[c]     ← peaks NOT matching the isotope pattern
//   chargeSNR = (ε + mul_factor × sig_pwr)
//             / (ε + noise_pwr + (1 − cos²) × sum_signal_squared[c])
//
// ── SNR approximation (self-noise) ───────────────────────────────────────────
// We do NOT run a separate noise-peak recruitment pass. OpenMS collects peaks
// in the same m/z window that fall between isotope positions to form noise_pwr.
// Replicating this requires a second scan of the raw spectrum per envelope.
//
// CURRENT APPROXIMATION: use the cosine mismatch fraction as self-noise.
//   chargeSNR ≈ cos² / (1 − cos² + ε)
//
// This is equivalent to the OpenMS formula when external noise_pwr = 0.
// It correctly ranks high-cosine over low-cosine envelopes but cannot
// distinguish two envelopes with the same cosine and different background noise.
//
// TODO: proper noise estimation
// The right fix is to collect non-matching peaks during Step 3 (isotope
// recruitment). After recruiting signal peaks within ±tolerancePpm of each
// expected isotope position, also collect raw spectrum peaks in the same m/z
// window that fall outside the tolerance band. Sum their squared intensities
// per charge to get per_charge_noise_pwr, pass it into AssignQscores alongside
// the envelope, and use it in the SNR formula. This would match the OpenMS
// model exactly for f[2] and f[3] and is the highest-priority scoring
// improvement after the current round of validation.
//
// ── PPM error: exact vs approximate ──────────────────────────────────────────
// f[4] is the average absolute ppm error of recruited peaks vs the theoretical
// Averagine isotope positions. This is most accurately computed during Step 3
// while the isotope index n for each recruited peak is still known:
//   theorMz = monoMass.ToMz(charge) + n * C13/z
//   error   = |obsMz − theorMz| / theorMz × 1e6
//
// The preferred call path is:
//   ScoreAndBuildEnvelopesWithPpmError → pairs of (envelope, avgPpmError)
//   → FLASHDeconvDeduplicator.Deduplicate (with ppm error propagated)
//   → FLASHDeconvScorer.AssignQscores(pairs, tolerancePpm)
//
// A fallback overload AssignQscores(envelopes, tolerancePpm) recomputes ppm
// error post-hoc from the peak list, which is less accurate because the isotope
// index n is reconstructed by rounding rather than known exactly.

using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Logistic-regression Qscore for <see cref="IsotopicEnvelope"/> objects produced
    /// by <see cref="FLASHDeconvolutionAlgorithm"/>. Ported from OpenMS
    /// <c>PeakGroupScoring.cpp</c> (Kyowon Jeong, BSD-3-Clause).
    /// </summary>
    internal static class FLASHDeconvScorer
    {
        // ── Logistic regression weights ───────────────────────────────────────
        // From OpenMS PeakGroupScoring.cpp weight_ vector (stored negated).
        // score = intercept + sum(w[i] * f[i])
        // Qscore = 1 / (1 + exp(score))
        //
        //   w[0] = -21.0476  isotopeCosine       high cosine → low score → high Qscore
        //   w[1] =   1.5045  inter-charge drop   (set to 0, single-charge envelope)
        //   w[2] =  -0.1303  log2(1+chargeSNR)   high SNR → lower score → higher Qscore
        //   w[3] =   0.183   SNR drop            (set to 0, single-charge envelope)
        //   w[4] =   0.1834  avgPPMError         high error → higher score → lower Qscore
        private static readonly double[] Weights = { -21.0476, 1.5045, -0.1303, 0.183, 0.1834 };
        private static readonly double Intercept = 17.804 + 0.5; // +0.5 matches OpenMS getQscore()

        // ── Primary path: pre-computed ppm error ──────────────────────────────

        /// <summary>
        /// Assigns Qscores using the exact ppm error computed during Step 3 while
        /// the isotope index was still known. Preferred over the envelope-only overload.
        /// </summary>
        /// <param name="pairs">
        /// (envelope with Score=cosine, avgPpmError in ppm) from
        /// <see cref="FLASHDeconvolutionAlgorithm.ScoreAndBuildEnvelopesWithPpmError"/>,
        /// after deduplication with ppm error propagated.
        /// </param>
        internal static IEnumerable<IsotopicEnvelope> AssignQscores(
            IEnumerable<(IsotopicEnvelope Envelope, double AvgPpmError)> pairs)
        {
            foreach (var (env, ppmError) in pairs)
            {
                double qscore = ComputeQscore(env.Score, ppmError);
                yield return new IsotopicEnvelope(
                    id: env.PrecursorId,
                    peaks: env.Peaks,
                    monoisotopicmass: env.MonoisotopicMass,
                    chargestate: env.Charge,
                    intensity: env.TotalIntensity,
                    score: qscore);
            }
        }

        // ── Fallback path: envelope only (ppm error reconstructed post-hoc) ──

        /// <summary>
        /// Assigns Qscores using only the envelope. Ppm error is reconstructed
        /// from the (mz, intensity) peak list by rounding to the nearest isotope
        /// index — less accurate than the primary path but usable when the
        /// pre-computed error is not available.
        /// </summary>
        internal static IEnumerable<IsotopicEnvelope> AssignQscores(
            IEnumerable<IsotopicEnvelope> envelopes,
            double tolerancePpm = 10.0)
        {
            foreach (var env in envelopes)
            {
                double ppmError = RecomputeAvgPpmError(env, tolerancePpm);
                double qscore = ComputeQscore(env.Score, ppmError);
                yield return new IsotopicEnvelope(
                    id: env.PrecursorId,
                    peaks: env.Peaks,
                    monoisotopicmass: env.MonoisotopicMass,
                    chargestate: env.Charge,
                    intensity: env.TotalIntensity,
                    score: qscore);
            }
        }

        // ── Core computation ──────────────────────────────────────────────────

        /// <summary>
        /// Computes the Qscore from a cosine and a known ppm error.
        /// </summary>
        internal static double ComputeQscore(double cosine, double avgPpmError)
        {
            // f[0]: isotopeCosine
            double f0 = cosine;

            // f[1]: inter-charge cosine drop — zero (single-charge envelope)
            double f1 = 0.0;

            // f[2]: log2(1 + chargeSNR) via self-noise approximation
            // SNR = cos² / (1 − cos² + ε)
            // See file header for the TODO on proper external noise collection.
            double cos2 = cosine * cosine;
            double selfNoiseSNR = cos2 / (1.0 - cos2 + 1e-6);
            double f2 = Math.Log(1.0 + selfNoiseSNR, 2.0);

            // f[3]: per-charge vs global SNR drop — zero (same reason as f[1])
            double f3 = 0.0;

            // f[4]: avg ppm error — passed in directly
            double f4 = avgPpmError;

            double linearScore = Intercept
                + Weights[0] * f0
                + Weights[1] * f1
                + Weights[2] * f2
                + Weights[3] * f3
                + Weights[4] * f4;

            return 1.0 / (1.0 + Math.Exp(linearScore));
        }

        // ── Post-hoc ppm error reconstruction (fallback) ─────────────────────

        private static double RecomputeAvgPpmError(IsotopicEnvelope env, double tolerancePpm)
        {
            if (env.Peaks.Count == 0 || env.Charge == 0) return tolerancePpm;
            int absCharge = Math.Abs(env.Charge);
            double isotopeStepMz = Constants.C13MinusC12 / absCharge;
            double apexMz = env.Peaks.MaxBy(p => p.intensity).mz;
            double totalPpm = 0.0;
            foreach (var (mz, _) in env.Peaks)
            {
                int n = (int)Math.Round((mz - apexMz) / isotopeStepMz);
                double theorMz = apexMz + n * isotopeStepMz;
                totalPpm += Math.Abs(mz - theorMz) / theorMz * 1e6;
            }
            return totalPpm / env.Peaks.Count;
        }

        // ── Utility ───────────────────────────────────────────────────────────

        /// <summary>
        /// Median intensity of the spectrum. Kept for callers that may want it
        /// for other purposes; no longer used internally by <see cref="AssignQscores"/>.
        /// </summary>
        internal static double ComputeMedianIntensity(MzSpectrum spectrum)
        {
            if (spectrum == null || spectrum.Size == 0) return 1.0;
            var intensities = spectrum.YArray.ToList();
            intensities.Sort();
            int mid = intensities.Count / 2;
            return intensities.Count % 2 == 0
                ? (intensities[mid - 1] + intensities[mid]) / 2.0
                : intensities[mid];
        }
    }
}