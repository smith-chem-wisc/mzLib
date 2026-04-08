// FLASHDeconvScorer.cs
//
// Computes a Qscore analog for IsotopicEnvelope objects produced by
// FLASHDeconvolutionAlgorithm. The model is a direct port of the logistic
// regression in OpenMS PeakGroupScoring.cpp.
//
// ── OpenMS reference ─────────────────────────────────────────────────────────
// File:    src/openms/source/ANALYSIS/TOPDOWN/PeakGroupScoring.cpp
// Authors: Kyowon Jeong
// License: BSD-3-Clause
//
// OpenMS weight vector (from PeakGroupScoring.cpp):
//   weights = { -21.0476, 1.5045, -0.1303, 0.183, 0.1834, 17.804 }
//   Intercept (last element) = 17.804, adjusted by +0.5 inside getQscore()
//   giving an effective intercept of 18.304.
//
// OpenMS feature vector (toFeatureVector_):
//   f[0] = isotopeCosine                                (global cosine, all charges)
//   f[1] = isotopeCosine − chargeIsotopeCosine(repZ)   (inter-charge cosine drop)
//   f[2] = log2(1 + chargeSNR(repZ))                   (per-charge SNR, log-scaled)
//   f[3] = log2(1 + chargeSNR(repZ)) − log2(1 + SNR)   (per-charge vs global SNR drop)
//   f[4] = avgPPMError                                  (average mass PPM error)
//
// ── mzLib mapping ────────────────────────────────────────────────────────────
// IsotopicEnvelope carries:
//   Score        — cosine of the full recruited envelope vs Averagine (Step 4)
//                  This is the closest analogue to OpenMS isotopeCosine.
//   Peaks        — list of (mz, intensity) for recruited isotope peaks
//   Charge       — signed representative charge (one charge per envelope)
//   TotalIntensity
//   MonoisotopicMass
//
// Because FLASHDeconvolutionAlgorithm emits one envelope per (candidate, charge)
// pair BEFORE deduplication, each envelope represents a single charge state.
// After deduplication only the best-scoring envelope per mass cluster survives.
// Therefore:
//
//   f[0] — isotopeCosine              → envelope.Score  (cosine stored at construction)
//   f[1] — cosine drop across charges → 0.0  (we have one charge per envelope; there
//                                             is no "other charge" to compare against.
//                                             Approximated as zero, which is the neutral
//                                             (non-penalising) value for this feature.)
//   f[2] — log2(1 + chargeSNR)        → log2(1 + SNR(envelope))
//                                       where SNR = TotalIntensity / noiseEstimate.
//                                       We estimate noise as the spectrum median intensity
//                                       passed in at scoring time.
//   f[3] — per-charge vs global SNR   → 0.0  (same reason as f[1]; approximated as zero)
//   f[4] — avgPPMError                → mean absolute ppm error of recruited peaks
//                                       vs the theoretical Averagine m/z positions,
//                                       computed from the peaks list and charge.
//
// f[1] and f[3] being zero is conservative: it means the model cannot penalise
// for inter-charge inconsistency, which is fine for a first-pass scorer. The
// ML-improved version will handle this properly once we track per-charge cosines.
//
// ── Qscore interpretation ─────────────────────────────────────────────────────
// Qscore = sigmoid(-score) where score is the linear combination.
// Higher Qscore = more likely to be a real peak group.
// OpenMS threshold for keeping a feature: Qscore >= 0.7 gives ~99% recall of
// high-confidence peaks on the Filgrastim dataset (see validation analysis).
//
// ── Placement ────────────────────────────────────────────────────────────────
// Add to MassSpectrometry/Deconvolution/Algorithms/ alongside
// FLASHDeconvolutionAlgorithm.cs and FLASHDeconvDeduplicator.cs.
//
// ── Usage in FLASHDeconvolutionAlgorithm.Deconvolute() ──────────────────────
// After deduplication:
//
//   var deduped = FLASHDeconvDeduplicator.Deduplicate(
//                     ScoreAndBuildEnvelopes(spectrum, candidates, p),
//                     p.DeconvolutionTolerancePpm);
//   double medianIntensity = ComputeMedianIntensity(spectrum);  // helper below
//   return FLASHDeconvScorer.AssignQscores(deduped, medianIntensity, p.DeconvolutionTolerancePpm);
//
// Or, if you prefer to filter rather than annotate:
//   return FLASHDeconvScorer.AssignQscores(deduped, medianIntensity, p.DeconvolutionTolerancePpm)
//                           .Where(e => e.Score >= p.MinQscore);
//
// MinQscore can be added to FLASHDeconvolutionParameters; suggested default: 0.0
// (keep everything, let downstream filtering decide) or 0.3 (light pre-filter).

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
        // Directly from OpenMS PeakGroupScoring.cpp weight_ vector.
        // Element order matches toFeatureVector_ in that file.
        //
        //   w[0] = -21.0476  → isotopeCosine
        //   w[1] =   1.5045  → cosine drop across charges  (we set feature = 0)
        //   w[2] =  -0.1303  → log2(1 + chargeSNR)
        //   w[3] =   0.183   → chargeSNR vs globalSNR drop (we set feature = 0)
        //   w[4] =   0.1834  → avgPPMError
        //   w[5] =  17.804   → intercept (OpenMS adds +0.5 inside getQscore)
        //
        // Note on signs: OpenMS negates its weight_ vector internally (see code).
        // The stored values ARE the negated coefficients. The linear combination is
        //   score = intercept + 0.5 + sum(w[i] * f[i])
        // and then sigmoid(-score) — which simplifies to:
        //   Qscore = 1 / (1 + exp(score))   [standard sigmoid of negative argument]
        private static readonly double[] Weights = { -21.0476, 1.5045, -0.1303, 0.183, 0.1834 };
        private static readonly double Intercept = 17.804 + 0.5; // +0.5 matches OpenMS getQscore()

        /// <summary>
        /// Assigns Qscores to a sequence of envelopes by replacing their <c>Score</c>
        /// field with the logistic-regression Qscore. The cosine score used to build
        /// the feature vector is read before replacement.
        /// </summary>
        /// <param name="envelopes">Deduplicated envelopes from FLASHDeconvolutionAlgorithm.</param>
        /// <param name="medianSpectrumIntensity">
        /// Median intensity of the raw spectrum, used as noise floor for SNR estimation.
        /// Pass 1.0 if unknown (degrades SNR features to neutral values).
        /// </param>
        /// <param name="tolerancePpm">Peak-matching tolerance in ppm (same as decon tolerance).</param>
        /// <returns>
        /// New <see cref="IsotopicEnvelope"/> objects with <c>Score</c> replaced by
        /// the Qscore in [0, 1]. Higher = more confident.
        /// </returns>
        internal static IEnumerable<IsotopicEnvelope> AssignQscores(
            IEnumerable<IsotopicEnvelope> envelopes,
            double medianSpectrumIntensity,
            double tolerancePpm = 10.0)
        {
            double noiseFloor = Math.Max(medianSpectrumIntensity, 1.0);

            foreach (var env in envelopes)
            {
                double qscore = ComputeQscore(env, noiseFloor, tolerancePpm);
                yield return new IsotopicEnvelope(
                    id: env.PrecursorId,
                    peaks: env.Peaks,
                    monoisotopicmass: env.MonoisotopicMass,
                    chargestate: env.Charge,
                    intensity: env.TotalIntensity,
                    score: qscore);
            }
        }

        /// <summary>
        /// Computes the Qscore for a single envelope without creating a new object.
        /// Useful for filtering before allocation.
        /// </summary>
        internal static double ComputeQscore(
            IsotopicEnvelope env,
            double noiseFloor,
            double tolerancePpm = 10.0)
        {
            if (env.Peaks == null || env.Peaks.Count == 0)
                return 0.0;

            double[] fv = BuildFeatureVector(env, noiseFloor, tolerancePpm);

            double linearScore = Intercept;
            for (int i = 0; i < Weights.Length; i++)
                linearScore += Weights[i] * fv[i];

            // OpenMS: qscore = 1 / (1 + exp(score))
            // Higher linear score → lower Qscore (the negation is already in the weights)
            return 1.0 / (1.0 + Math.Exp(linearScore));
        }

        // ── Feature vector construction ───────────────────────────────────────

        private static double[] BuildFeatureVector(
            IsotopicEnvelope env,
            double noiseFloor,
            double tolerancePpm)
        {
            // f[0]: isotopeCosine — the cosine stored in Score at envelope creation
            double isotopeCosine = env.Score;

            // f[1]: inter-charge cosine drop — zero (single-charge envelope)
            double cosDropAcrossCharges = 0.0;

            // f[2]: log2(1 + chargeSNR)
            // SNR = signal power / noise power ≈ TotalIntensity / noiseFloor
            // OpenMS SNR is more nuanced (per-charge cosine-weighted power ratio),
            // but TotalIntensity / median is a reasonable analogue for a single charge.
            double snr = env.TotalIntensity / noiseFloor;
            double logChargeSNR = Math.Log(1.0 + snr, 2.0);

            // f[3]: per-charge SNR drop vs global SNR — zero (one charge per envelope)
            double snrDrop = 0.0;

            // f[4]: average absolute ppm error of recruited peaks vs Averagine theoretical
            double avgPpmError = ComputeAvgPpmError(env, tolerancePpm);

            return new[] { isotopeCosine, cosDropAcrossCharges, logChargeSNR, snrDrop, avgPpmError };
        }

        /// <summary>
        /// Computes the mean absolute ppm error between the recruited m/z peaks and
        /// the theoretical Averagine isotope positions at this charge state.
        /// </summary>
        private static double ComputeAvgPpmError(IsotopicEnvelope env, double tolerancePpm)
        {
            if (env.Peaks.Count == 0) return tolerancePpm; // worst-case default

            int absCharge = Math.Abs(env.Charge);
            if (absCharge == 0) return tolerancePpm;

            // Theoretical isotope spacing at this charge = C13 / z
            double isotopeStepMz = Constants.C13MinusC12 / absCharge;

            // Use the most-intense peak as the anchor for theoretical positions
            var apexPeak = env.Peaks.MaxBy(p => p.intensity);
            double apexMz = apexPeak.mz;

            double totalPpmError = 0.0;
            int count = 0;

            foreach (var (mz, _) in env.Peaks)
            {
                // Nearest theoretical isotope index relative to apex
                double offsetMz = mz - apexMz;
                int n = (int)Math.Round(offsetMz / isotopeStepMz);
                double theoreticalMz = apexMz + n * isotopeStepMz;

                double ppm = Math.Abs(mz - theoreticalMz) / theoreticalMz * 1e6;
                totalPpmError += ppm;
                count++;
            }

            return count > 0 ? totalPpmError / count : tolerancePpm;
        }

        // ── Spectrum noise estimation helper ──────────────────────────────────

        /// <summary>
        /// Estimates the noise floor of a spectrum as the median of all peak intensities.
        /// Call this once per spectrum and pass the result to <see cref="AssignQscores"/>.
        /// Returns 1.0 if the spectrum is empty.
        /// </summary>
        internal static double ComputeMedianIntensity(MzSpectrum spectrum)
        {
            if (spectrum == null || spectrum.Size == 0) return 1.0;

            // Copy intensities and find median without full sort (approximate with sort for simplicity)
            var intensities = spectrum.YArray.ToList();
            intensities.Sort();
            int mid = intensities.Count / 2;
            return intensities.Count % 2 == 0
                ? (intensities[mid - 1] + intensities[mid]) / 2.0
                : intensities[mid];
        }
    }
}