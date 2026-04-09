using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry.MzSpectra;

namespace MassSpectrometry
{
    /// <summary>
    /// Generic post-hoc scorer for <see cref="IsotopicEnvelope"/> objects produced
    /// by any mzLib deconvolution algorithm (Classic, IsoDec, FLASHDeconv).
    ///
    /// All scoring is performed from the envelope's <see cref="IsotopicEnvelope.Peaks"/>
    /// list and an <see cref="AverageResidue"/> model. No raw <see cref="MzSpectrum"/>
    /// access is required, supporting callers that deconvolute and discard the spectrum.
    ///
    /// The score produced by <see cref="ComputeScore"/> is a logistic function in [0, 1]
    /// combining three features:
    /// <list type="bullet">
    ///   <item><see cref="EnvelopeScoreFeatures.AveragineCosineSimilarity"/> — isotope pattern fidelity</item>
    ///   <item><see cref="EnvelopeScoreFeatures.AvgPpmError"/> — mass accuracy</item>
    ///   <item><see cref="EnvelopeScoreFeatures.PeakCompleteness"/> — observed vs expected peak count</item>
    /// </list>
    ///
    /// The logistic weights are provisional and should be recalibrated against labelled
    /// training data once available. See <see cref="CoefficientCosine"/> and related
    /// constants.
    /// </summary>
    public static class DeconvolutionScorer
    {
        // ── Logistic regression weights ───────────────────────────────────────
        // Provisional weights — recalibrate against labelled training data.
        // Calibrated to produce scores broadly similar in range to the FLASHDeconv
        // Qscore on intact protein top-down data; not fit to a labelled data set.
        private const double CoefficientCosine      = -14.0;  // high cosine → low linear score → high Qscore
        private const double CoefficientPpmError    =   0.15; // high ppm error → higher linear score → lower Qscore
        private const double CoefficientCompleteness = -3.0;  // high completeness → low linear score → high Qscore
        private const double Intercept              =   9.5;

        // ── Feature matching tolerance ────────────────────────────────────────
        private const double MatchTolerancePpm = 10.0;

        // ── Public API ────────────────────────────────────────────────────────

        /// <summary>
        /// Computes the three scoring features for a single envelope using the supplied
        /// Averagine model. Does not access the raw spectrum.
        /// </summary>
        /// <param name="envelope">Envelope to score.</param>
        /// <param name="model">
        /// Averagine model used for theoretical isotope pattern lookup.
        /// Should be the same model used during deconvolution.
        /// </param>
        public static EnvelopeScoreFeatures ComputeFeatures(
            IsotopicEnvelope envelope,
            AverageResidue model)
        {
            if (envelope == null) throw new ArgumentNullException(nameof(envelope));
            if (model    == null) throw new ArgumentNullException(nameof(model));

            // ── Averagine lookup ──────────────────────────────────────────────
            int avgIdx = model.GetMostIntenseMassIndex(envelope.MonoisotopicMass);

            double[] rawMasses  = model.GetAllTheoreticalMasses(avgIdx);
            double[] rawIntens  = model.GetAllTheoreticalIntensities(avgIdx);

            // Arrays from AverageResidue are sorted intensity-descending (index 0 = apex).
            // Re-sort to mass-ascending so index 0 = monoisotopic.
            var sorted = rawMasses.Zip(rawIntens)
                .OrderBy(pair => pair.First)
                .ToArray();
            double[] avgMassAsc  = sorted.Select(p => p.First).ToArray();
            double[] avgIntAsc   = sorted.Select(p => p.Second).ToArray();
            int      nIso        = avgMassAsc.Length;

            // apexDaFromMono: distance in Da from monoisotopic to apex
            double apexDaFromMono = model.GetDiffToMonoisotopic(avgIdx);
            int    absCharge      = Math.Abs(envelope.Charge);

            // Guard: if charge is 0 (unusual but defensive) fall back gracefully
            if (absCharge == 0)
                return new EnvelopeScoreFeatures(0.0, MatchTolerancePpm, 0.0);

            double isotopeStepMz = Constants.C13MinusC12 / absCharge;

            // ── Build observed intensity vector for cosine and completeness ───
            // For each theoretical isotope position n, find the closest peak in
            // envelope.Peaks within MatchTolerancePpm.
            double monoMass   = envelope.MonoisotopicMass;
            double monoMz     = monoMass.ToMz(absCharge);
            var    peakList   = envelope.Peaks;

            double[] observed     = new double[nIso];
            bool[]   peakMatched  = new bool[nIso];

            for (int n = 0; n < nIso; n++)
            {
                double theorMz  = monoMz + n * isotopeStepMz;
                double tolMz    = theorMz * MatchTolerancePpm * 1e-6;

                // Find the peak in the envelope closest to this theoretical position.
                double bestDist = double.MaxValue;
                double bestIntensity = 0.0;

                foreach (var (mz, intensity) in peakList)
                {
                    double dist = Math.Abs(mz - theorMz);
                    if (dist < tolMz && dist < bestDist)
                    {
                        bestDist      = dist;
                        bestIntensity = intensity;
                    }
                }

                if (bestDist < double.MaxValue)
                {
                    observed[n]    = bestIntensity;
                    peakMatched[n] = true;
                }
                // else observed[n] = 0 (default), peakMatched[n] = false
            }

            // ── AveragineCosineSimilarity ─────────────────────────────────────
            double cosine = SpectralSimilarity.CosineOfAlignedVectors(observed, avgIntAsc);
            cosine = Math.Max(0.0, Math.Min(1.0, cosine));

            // ── AvgPpmError ───────────────────────────────────────────────────
            double avgPpmError = ComputeAvgPpmError(peakList, absCharge, isotopeStepMz);

            // ── PeakCompleteness ──────────────────────────────────────────────
            double maxTheoInt    = avgIntAsc.Max();
            double threshold1pct = maxTheoInt * 0.01;
            int    expectedPeaks = 0;
            int    observedPeaks = 0;

            for (int n = 0; n < nIso; n++)
            {
                if (avgIntAsc[n] > threshold1pct)
                {
                    expectedPeaks++;
                    if (peakMatched[n]) observedPeaks++;
                }
            }

            double completeness = expectedPeaks == 0
                ? 0.0
                : Math.Max(0.0, Math.Min(1.0, (double)observedPeaks / expectedPeaks));

            return new EnvelopeScoreFeatures(cosine, avgPpmError, completeness);
        }

        /// <summary>
        /// Combines the three features into a single score in [0, 1] using a logistic
        /// function.
        /// <para>
        /// Higher scores indicate higher confidence that the envelope represents a true
        /// deconvolution result. The logistic weights are provisional — see the constants
        /// at the top of this class.
        /// </para>
        /// </summary>
        public static double ComputeScore(EnvelopeScoreFeatures features)
        {
            double linear = Intercept
                + CoefficientCosine       * features.AveragineCosineSimilarity
                + CoefficientPpmError     * features.AvgPpmError
                + CoefficientCompleteness * features.PeakCompleteness;

            return 1.0 / (1.0 + Math.Exp(linear));
        }

        /// <summary>
        /// Convenience method: computes features then score for a single envelope.
        /// </summary>
        public static double ScoreEnvelope(IsotopicEnvelope envelope, AverageResidue model)
            => ComputeScore(ComputeFeatures(envelope, model));

        /// <summary>
        /// Applies <see cref="ScoreEnvelope"/> to each envelope in the sequence and
        /// yields <c>(envelope, genericScore)</c> pairs. The original
        /// <see cref="IsotopicEnvelope.Score"/> field is not mutated.
        /// </summary>
        /// <param name="envelopes">Envelopes to score.</param>
        /// <param name="model">Averagine model for theoretical isotope lookup.</param>
        public static IEnumerable<(IsotopicEnvelope Envelope, double Score)> ScoreEnvelopes(
            IEnumerable<IsotopicEnvelope> envelopes,
            AverageResidue model)
        {
            foreach (var env in envelopes)
                yield return (env, ScoreEnvelope(env, model));
        }

        // ── Private helpers ───────────────────────────────────────────────────

        /// <summary>
        /// Computes the mean absolute ppm error of the envelope's peaks relative
        /// to the apex-anchored theoretical isotope grid.
        /// </summary>
        private static double ComputeAvgPpmError(
            IReadOnlyList<(double mz, double intensity)> peaks,
            int absCharge,
            double isotopeStepMz)
        {
            if (peaks.Count == 0) return MatchTolerancePpm;

            // Anchor at the most intense peak
            double apexMz = peaks.MaxBy(p => p.intensity).mz;

            double totalPpm = 0.0;
            foreach (var (mz, _) in peaks)
            {
                int    n       = (int)Math.Round((mz - apexMz) / isotopeStepMz);
                double theorMz = apexMz + n * isotopeStepMz;
                totalPpm += Math.Abs(mz - theorMz) / theorMz * 1e6;
            }

            return totalPpm / peaks.Count;
        }
    }
}
