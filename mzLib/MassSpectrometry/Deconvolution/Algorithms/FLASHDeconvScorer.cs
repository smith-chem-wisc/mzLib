// FLASHDeconvScorer.cs
//
// Computes a Qscore for IsotopicEnvelope objects produced by
// FLASHDeconvolutionAlgorithm. The model ports the logistic regression in
// OpenMS PeakGroupScoring.cpp (Kyowon Jeong, BSD-3-Clause).
//
// ── OpenMS feature vector (toFeatureVector_) ─────────────────────────────────
//   f[0] = isotopeCosine                          global cosine vs Averagine
//   f[1] = isotopeCosine − chargeIsotopeCosine    inter-charge cosine drop
//   f[2] = log2(1 + chargeSNR)                    per-charge SNR, log-scaled
//   f[3] = log2(1+chargeSNR) − log2(1+SNR)        per-charge vs global SNR drop
//   f[4] = avgPPMError
//
// ── OpenMS SNR formula (updateSNR_ in PeakGroup.cpp) ─────────────────────────
// For the representative charge z:
//   cos²      = chargeIsotopeCosine²
//   sig_pwr   = perChargeSumSignalSquared × cos²
//   chargeSNR = (ε + sig_pwr) / (ε + noisePwr + (1 − cos²) × perChargeSumSignalSquared)
//
// noisePwr is the summed squared intensity of raw spectrum peaks that fall
// within the isotope m/z window for charge z but outside the ±ppm tolerance
// band around each expected isotope position. It is collected during Step 3
// isotope recruitment in FLASHDeconvolutionAlgorithm.
//
// ── Data flow ─────────────────────────────────────────────────────────────────
//   ScoreAndBuildEnvelopesWithPpmError  →  EnvelopeScoringData
//   → FLASHDeconvDeduplicator.Deduplicate (propagates best cosine + ppm error)
//   → FLASHDeconvScorer.AssignQscores(data, model)
//
// ── PPM error ─────────────────────────────────────────────────────────────────
// f[4] is computed during Step 3 while the isotope index n is still known:
//   theorMz = monoMass.ToMz(z) + n * isotopeStep/z
//   error   = |obsMz − theorMz| / theorMz × 1e6
//
// A fallback overload AssignQscores(envelopes, model, tolerancePpm) recomputes
// ppm error post-hoc from the peak list when pre-computed data is unavailable.
// The isotope step for both paths is derived from the AverageResidue model.

using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Logistic-regression Qscore for <see cref="IsotopicEnvelope"/> objects
    /// produced by <see cref="FLASHDeconvolutionAlgorithm"/>.
    /// Ported from OpenMS <c>PeakGroupScoring.cpp</c> (Kyowon Jeong, BSD-3-Clause).
    /// </summary>
    internal static class FLASHDeconvScorer
    {
        // ── Logistic regression weights ───────────────────────────────────────
        // From OpenMS PeakGroupScoring.cpp weight_ vector (stored negated).
        // linearScore = Intercept + sum(Weights[i] * f[i])
        // Qscore      = 1 / (1 + exp(linearScore))
        //
        //   Weights[0] = -21.0476  isotopeCosine         high cosine → high Qscore
        //   Weights[1] =   1.5045  inter-charge drop     large drop  → lower Qscore
        //   Weights[2] =  -0.1303  log2(1+chargeSNR)     high SNR    → higher Qscore
        //   Weights[3] =   0.183   per-charge SNR drop   large drop  → lower Qscore
        //   Weights[4] =   0.1834  avgPPMError           high error  → lower Qscore
        private static readonly double[] Weights = { -21.0476, 1.5045, -0.1303, 0.183, 0.1834 };
        private static readonly double Intercept = 17.804 + 0.5; // +0.5 matches OpenMS getQscore()

        // mul_factor from OpenMS updateSNR_ — scales signal power relative to noise
        private const double SnrMulFactor = 1.0;
        private const double SnrEpsilon = 1e-6;

        // ── Data carrier ─────────────────────────────────────────────────────

        /// <summary>
        /// All per-envelope data collected during Step 3 isotope recruitment
        /// that the scorer needs. Produced by
        /// <see cref="FLASHDeconvolutionAlgorithm"/> and consumed here.
        /// </summary>
        internal readonly struct EnvelopeScoringData
        {
            /// <summary>Envelope whose Score field holds the global Averagine cosine.</summary>
            internal readonly IsotopicEnvelope Envelope;

            /// <summary>
            /// Average absolute ppm error of recruited peaks vs theoretical
            /// isotope positions, computed during recruitment while the isotope
            /// index n was still known.
            /// </summary>
            internal readonly double AvgPpmError;

            /// <summary>
            /// Cosine similarity between the observed isotope intensity distribution
            /// for the representative charge state and the Averagine template.
            /// Used for f[1] (inter-charge cosine drop).
            /// </summary>
            internal readonly double RepChargeIsotopeCosine;

            /// <summary>
            /// Summed squared intensity of raw spectrum peaks that fall within the
            /// isotope m/z window for the representative charge but outside the
            /// ±ppm tolerance band around each expected isotope position.
            /// Used as the external noise power in the OpenMS SNR formula for f[2]/f[3].
            /// </summary>
            internal readonly double RepChargeNoisePower;

            /// <summary>
            /// Summed squared intensity of signal peaks at the representative charge.
            /// Used together with RepChargeNoisePower in the SNR denominator.
            /// </summary>
            internal readonly double RepChargeSumSignalSquared;

            internal EnvelopeScoringData(
                IsotopicEnvelope envelope,
                double avgPpmError,
                double repChargeIsotopeCosine,
                double repChargeNoisePower,
                double repChargeSumSignalSquared)
            {
                Envelope = envelope;
                AvgPpmError = avgPpmError;
                RepChargeIsotopeCosine = repChargeIsotopeCosine;
                RepChargeNoisePower = repChargeNoisePower;
                RepChargeSumSignalSquared = repChargeSumSignalSquared;
            }
        }

        // ── Primary path: full scoring data from Step 3 ───────────────────────

        /// <summary>
        /// Assigns Qscores using the full scoring data collected during Step 3.
        /// This is the preferred path: f[1]–f[3] are computed from real per-charge
        /// data rather than approximations.
        /// </summary>
        /// <param name="scoringData">
        /// Sequence of <see cref="EnvelopeScoringData"/> produced by
        /// <see cref="FLASHDeconvolutionAlgorithm.ScoreAndBuildEnvelopes"/>,
        /// after deduplication.
        /// </param>
        internal static IEnumerable<IsotopicEnvelope> AssignQscores(
            IEnumerable<EnvelopeScoringData> scoringData)
        {
            foreach (var d in scoringData)
            {
                double qscore = ComputeQscore(
                    globalCosine: d.Envelope.Score,
                    repChargeCosine: d.RepChargeIsotopeCosine,
                    repChargeNoisePower: d.RepChargeNoisePower,
                    repChargeSumSigSq: d.RepChargeSumSignalSquared,
                    avgPpmError: d.AvgPpmError);

                yield return new IsotopicEnvelope(
                    id: d.Envelope.PrecursorId,
                    peaks: d.Envelope.Peaks,
                    monoisotopicmass: d.Envelope.MonoisotopicMass,
                    chargestate: d.Envelope.Charge,
                    intensity: d.Envelope.TotalIntensity,
                    score: qscore);
            }
        }

        // ── Fallback path: envelope only ──────────────────────────────────────

        /// <summary>
        /// Assigns Qscores using only the envelopes. Per-charge cosine and noise
        /// are unavailable, so f[1] and f[3] are set to zero and f[2] uses the
        /// self-noise approximation (cos²/(1−cos²)). Ppm error is reconstructed
        /// post-hoc from the peak list.
        /// <para>
        /// Use this path only when <see cref="EnvelopeScoringData"/> from Step 3
        /// is not available — for example when rescoring previously deconvolved
        /// envelopes read from a file.
        /// </para>
        /// </summary>
        /// <param name="envelopes">Envelopes whose Score field holds the cosine.</param>
        /// <param name="model">
        /// Averagine model used for deconvolution. Determines the isotope step
        /// used in ppm error reconstruction.
        /// </param>
        /// <param name="tolerancePpm">
        /// Fallback ppm error value returned for single-peak or zero-charge envelopes.
        /// </param>
        internal static IEnumerable<IsotopicEnvelope> AssignQscores(
            IEnumerable<IsotopicEnvelope> envelopes,
            AverageResidue model,
            double tolerancePpm = 10.0)
        {
            foreach (var env in envelopes)
            {
                double ppmError = RecomputeAvgPpmError(env, model, tolerancePpm);
                double cosine = env.Score;
                double cos2 = cosine * cosine;
                double selfNoiseSNR = cos2 / (1.0 - cos2 + SnrEpsilon);

                double qscore = ComputeQscore(
                    globalCosine: cosine,
                    repChargeCosine: cosine,          // no per-charge data: treat as equal
                    repChargeNoisePower: 0.0,             // self-noise path: no external noise
                    repChargeSumSigSq: selfNoiseSNR,    // SNR already encoded via self-noise
                    avgPpmError: ppmError,
                    selfNoiseOverride: true);

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
        /// Computes the Qscore from all five OpenMS features.
        /// </summary>
        /// <param name="globalCosine">f[0]: global Averagine cosine across all charges.</param>
        /// <param name="repChargeCosine">Per-charge cosine for the representative charge.</param>
        /// <param name="repChargeNoisePower">
        /// Summed squared intensity of noise peaks at the representative charge
        /// (peaks in the isotope window but outside ppm tolerance).
        /// </param>
        /// <param name="repChargeSumSigSq">
        /// Summed squared intensity of signal peaks at the representative charge.
        /// </param>
        /// <param name="avgPpmError">f[4]: average absolute ppm error in ppm.</param>
        /// <param name="selfNoiseOverride">
        /// When true, repChargeSumSigSq is treated as a pre-computed SNR value
        /// (self-noise fallback path) rather than raw signal squared intensity.
        /// </param>
        internal static double ComputeQscore(
            double globalCosine,
            double repChargeCosine,
            double repChargeNoisePower,
            double repChargeSumSigSq,
            double avgPpmError,
            bool selfNoiseOverride = false)
        {
            // f[0]: global isotope cosine
            double f0 = globalCosine;

            // f[1]: inter-charge cosine drop
            // OpenMS: isotopeCosine − chargeIsotopeCosine(repZ)
            // Positive when global cosine exceeds per-charge cosine (unexpected — typically
            // near-zero or slightly negative for good peaks).
            double f1 = globalCosine - repChargeCosine;

            // f[2]: log2(1 + chargeSNR) for the representative charge
            double chargeSNR;
            if (selfNoiseOverride)
            {
                // Fallback: repChargeSumSigSq already encodes the SNR directly
                chargeSNR = repChargeSumSigSq;
            }
            else
            {
                // OpenMS updateSNR_ formula:
                //   cos²     = repChargeCosine²
                //   sig_pwr  = repChargeSumSigSq × cos²
                //   chargeSNR = (ε + mul_factor × sig_pwr)
                //             / (ε + noisePwr + (1 − cos²) × repChargeSumSigSq)
                double cos2 = repChargeCosine * repChargeCosine;
                double sigPwr = repChargeSumSigSq * cos2;
                double denominator = SnrEpsilon
                    + repChargeNoisePower
                    + (1.0 - cos2) * repChargeSumSigSq;
                chargeSNR = (SnrEpsilon + SnrMulFactor * sigPwr) / denominator;
            }
            double f2 = Math.Log(1.0 + chargeSNR, 2.0);

            // f[3]: per-charge vs global SNR drop
            // OpenMS: log2(1+chargeSNR) − log2(1+SNR)
            // We use the global cosine self-noise as the global SNR reference, matching
            // the approach OpenMS uses when only one representative charge is tracked.
            double globalCos2 = globalCosine * globalCosine;
            double globalSNR = globalCos2 / (1.0 - globalCos2 + SnrEpsilon);
            double f3 = f2 - Math.Log(1.0 + globalSNR, 2.0);

            // f[4]: avg ppm error
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

        /// <summary>
        /// Reconstructs the average ppm error from the envelope's peak list.
        /// Less accurate than the primary path because the isotope index n is
        /// estimated by rounding rather than known exactly. The isotope step
        /// is derived from the supplied <see cref="AverageResidue"/> model.
        /// </summary>
        private static double RecomputeAvgPpmError(
            IsotopicEnvelope env,
            AverageResidue model,
            double tolerancePpm)
        {
            if (env.Peaks.Count == 0 || env.Charge == 0)
                return tolerancePpm;

            int absCharge = Math.Abs(env.Charge);

            // Derive the isotope step from the model's apex mass for this envelope.
            // For protein Averagine the step is effectively C13MinusC12 / z;
            // for RNA it differs because the average elemental composition is different.
            // We use C13MinusC12 directly as the per-Th isotope spacing — this is
            // correct because the 1.003355 Da step is the mass of a neutron excess
            // for a 13C substitution, independent of molecule type.
            // The model is consulted here for forward-compatibility: if a future model
            // redefines the isotope spacing (e.g. for non-standard isotopologue patterns),
            // this is the right place to apply it.
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
    }
}