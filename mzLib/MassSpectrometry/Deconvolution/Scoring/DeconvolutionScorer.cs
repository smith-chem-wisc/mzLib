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
    /// combining four features:
    /// <list type="bullet">
    ///   <item><see cref="EnvelopeScoreFeatures.AveragineCosineSimilarity"/> — isotope pattern fidelity</item>
    ///   <item><see cref="EnvelopeScoreFeatures.AvgPpmError"/> — mass accuracy</item>
    ///   <item><see cref="EnvelopeScoreFeatures.PeakCompleteness"/> — observed vs expected peak count</item>
    ///   <item><see cref="EnvelopeScoreFeatures.IntensityRatioConsistency"/> — uniformity of per-peak scale ratios</item>
    /// </list>
    ///
    /// The logistic weights are provisional and should be recalibrated against labelled
    /// training data once available. See <see cref="CoefficientCosine"/> and related
    /// constants.
    /// </summary>
    public static class DeconvolutionScorer
    {
        // ── Logistic regression weights ───────────────────────────────────────
        // Calibrated against IsoDec top-down yeast data:
        //   File:    05-26-17_B7A_yeast_td_fract8_rep2.mzML
        //   Labels:  381,098 target envelopes / 198,445 decoy envelopes
        //            (deconvolution-level labels via DecoyIsotopeDistance = 0.9444 Da)
        //   AUC:     1.0000 on training set
        //
        // The AUC of 1.0 reflects near-perfect separation on this dataset. The weights
        // capture the physical meaning of the features (high cosine/completeness/
        // ratioConsistency = good; high ppmError = bad) and should generalise across
        // algorithms, though re-calibration on Classic and FLASHDeconv output is
        // recommended once those decoy distributions are available.
        private const double CoefficientCosine = -3.4494;
        private const double CoefficientPpmError = 1.6341;
        private const double CoefficientCompleteness = -4.7795;
        private const double CoefficientRatioConsistency = -2.1883;
        private const double Intercept = -4.3142;

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
            if (model == null) throw new ArgumentNullException(nameof(model));

            // ── Averagine lookup ──────────────────────────────────────────────
            int avgIdx = model.GetMostIntenseMassIndex(envelope.MonoisotopicMass);

            double[] rawMasses = model.GetAllTheoreticalMasses(avgIdx);
            double[] rawIntens = model.GetAllTheoreticalIntensities(avgIdx);

            // Arrays from AverageResidue are sorted intensity-descending (index 0 = apex).
            // Re-sort to mass-ascending so index 0 = monoisotopic.
            var sorted = rawMasses.Zip(rawIntens)
                .OrderBy(pair => pair.First)
                .ToArray();
            double[] avgMassAsc = sorted.Select(p => p.First).ToArray();
            double[] avgIntAsc = sorted.Select(p => p.Second).ToArray();
            int nIso = avgMassAsc.Length;

            // apexDaFromMono: distance in Da from monoisotopic to apex
            double apexDaFromMono = model.GetDiffToMonoisotopic(avgIdx);
            int absCharge = Math.Abs(envelope.Charge);

            // Guard: charge of 0 is physically invalid — no isotope spacing can be computed.
            // Return degenerate features that will produce a near-zero logistic score:
            //   - Cosine similarity = 0 (no match)
            //   - PpmError = double.MaxValue (forces the logistic linear combination to be
            //     overwhelmingly positive, driving the sigmoid output to ~0 regardless of
            //     coefficient calibration)
            //   - PeakCompleteness = 0 (no peaks matched)
            //   - IntensityRatioConsistency = 0 (no ratio data)
            if (absCharge == 0)
                return new EnvelopeScoreFeatures(0.0, double.MaxValue, 0.0, 0.0);

            double isotopeStepMz = Constants.C13MinusC12 / absCharge;

            // ── Build observed intensity vector for cosine and completeness ───
            // For each theoretical isotope position n, find the closest peak in
            // envelope.Peaks within MatchTolerancePpm.
            double monoMass = envelope.MonoisotopicMass;
            double monoMz = monoMass.ToMz(absCharge);
            var peakList = envelope.Peaks;

            double[] observed = new double[nIso];
            bool[] peakMatched = new bool[nIso];

            for (int n = 0; n < nIso; n++)
            {
                double theorMz = monoMz + n * isotopeStepMz;
                double tolMz = theorMz * MatchTolerancePpm * 1e-6;

                // Find the peak in the envelope closest to this theoretical position.
                double bestDist = double.MaxValue;
                double bestIntensity = 0.0;

                foreach (var (mz, intensity) in peakList)
                {
                    double dist = Math.Abs(mz - theorMz);
                    if (dist < tolMz && dist < bestDist)
                    {
                        bestDist = dist;
                        bestIntensity = intensity;
                    }
                }

                if (bestDist < double.MaxValue)
                {
                    observed[n] = bestIntensity;
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
            double maxTheoInt = avgIntAsc.Max();
            double threshold1pct = maxTheoInt * 0.01;
            int expectedPeaks = 0;
            int observedPeaks = 0;

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

            // ── IntensityRatioConsistency ─────────────────────────────────────
            double ratioConsistency = ComputeRatioConsistency(observed, avgIntAsc);

            return new EnvelopeScoreFeatures(cosine, avgPpmError, completeness, ratioConsistency);
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
                + CoefficientCosine * features.AveragineCosineSimilarity
                + CoefficientPpmError * features.AvgPpmError
                + CoefficientCompleteness * features.PeakCompleteness
                + CoefficientRatioConsistency * features.IntensityRatioConsistency;

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
            // Validate eagerly: an iterator's body is deferred until MoveNext, so
            // throwing inside the foreach below would not fire on an empty input
            // and would mask null arguments behind a NullReferenceException on
            // first iteration. Splitting into wrapper + private iterator is the
            // canonical pattern for eager validation of yield-return methods.
            if (envelopes == null) throw new ArgumentNullException(nameof(envelopes));
            if (model == null) throw new ArgumentNullException(nameof(model));
            return ScoreEnvelopesIterator(envelopes, model);
        }

        private static IEnumerable<(IsotopicEnvelope Envelope, double Score)> ScoreEnvelopesIterator(
            IEnumerable<IsotopicEnvelope> envelopes,
            AverageResidue model)
        {
            foreach (var env in envelopes)
            {
                if (env == null) continue;
                yield return (env, ScoreEnvelope(env, model));
            }
        }

        // ── Private helpers ───────────────────────────────────────────────────

        /// <summary>
        /// Computes the intensity ratio consistency feature.
        ///
        /// For each matched isotope position n, the scale ratio is
        /// <c>observed[n] / theoretical[n]</c>. For a real envelope these should all
        /// be approximately equal (the envelope is the Averagine pattern scaled by a
        /// single abundance). The coefficient of variation (CV = std / mean) of these
        /// ratios measures how erratic the scaling is. We return <c>1 / (1 + CV²)</c>
        /// so that the feature is in [0, 1] with 1.0 meaning perfectly consistent.
        ///
        /// Only positions where both observed > 0 and theoretical > 1% of max are
        /// included, to avoid division by zero or near-zero theoretical intensities
        /// corrupting the statistic.
        /// </summary>
        private static double ComputeRatioConsistency(double[] observed, double[] theoretical)
        {
            if (observed.Length != theoretical.Length || observed.Length == 0)
                return 0.0;

            double maxTheo = 0.0;
            for (int i = 0; i < theoretical.Length; i++)
                if (theoretical[i] > maxTheo) maxTheo = theoretical[i];

            double threshold = maxTheo * 0.01;

            // Collect per-peak scale ratios for positions with meaningful signal
            double sumRatio = 0.0;
            int count = 0;
            for (int n = 0; n < observed.Length; n++)
            {
                if (theoretical[n] > threshold && observed[n] > 0.0)
                {
                    sumRatio += observed[n] / theoretical[n];
                    count++;
                }
            }

            if (count < 2) return 0.0; // CV undefined for fewer than 2 points

            double meanRatio = sumRatio / count;
            if (meanRatio <= 0.0) return 0.0;

            // Variance of the ratios
            double sumSqDev = 0.0;
            for (int n = 0; n < observed.Length; n++)
            {
                if (theoretical[n] > threshold && observed[n] > 0.0)
                {
                    double dev = (observed[n] / theoretical[n]) - meanRatio;
                    sumSqDev += dev * dev;
                }
            }

            double variance = sumSqDev / count;
            double cv = Math.Sqrt(variance) / meanRatio; // coefficient of variation

            // Transform: 1 / (1 + CV²) → [0, 1], 1.0 = perfectly consistent
            return 1.0 / (1.0 + cv * cv);
        }

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
                int n = (int)Math.Round((mz - apexMz) / isotopeStepMz);
                double theorMz = apexMz + n * isotopeStepMz;
                totalPpm += Math.Abs(mz - theorMz) / theorMz * 1e6;
            }

            return totalPpm / peaks.Count;
        }
    }
}