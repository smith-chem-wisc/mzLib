using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry.MzSpectra;

namespace MassSpectrometry
{
    /// <summary>
    /// Post-hoc deconvolution scorer that computes a method-agnostic quality score in [0, 1] for
    /// any <see cref="IsotopicEnvelope"/>. The scorer derives features from the envelope's peak
    /// list and an Averagine model — it does not access the original spectrum and does not depend
    /// on which deconvolution algorithm produced the envelope. Scores from different algorithms
    /// are therefore directly comparable.
    ///
    /// Provisional logistic weights are calibrated to produce score &gt; 0.5 for high-quality
    /// envelopes (cosine ≥ 0.95, ppm error ≤ a few, completeness and ratio consistency near 1.0)
    /// and score &lt; 0.5 for poor envelopes (cosine ≤ 0.3, large ppm error, low completeness).
    /// They are expected to be revised once labelled training data is available.
    /// </summary>
    public static class DeconvolutionScorer
    {
        // ── Feature-extraction constants ──────────────────────────────────────

        /// <summary>
        /// Tolerance window (ppm) used to match observed peaks to theoretical isotope positions
        /// when computing cosine similarity and peak completeness. Peaks outside this window
        /// do not contribute to the observed-vs-theoretical alignment.
        /// </summary>
        private const double MatchTolerancePpm = 10.0;

        /// <summary>
        /// Default ppm error returned when an envelope has no peaks. Chosen to equal
        /// the matching tolerance so an empty envelope is treated as a borderline (not best)
        /// case rather than as a perfect match.
        /// </summary>
        private const double DefaultPpmErrorForEmpty = MatchTolerancePpm;

        /// <summary>
        /// Theoretical Averagine peaks below this fraction of the most-intense theoretical peak
        /// are excluded from the "expected peaks" denominator in PeakCompleteness and from the
        /// ratio-consistency calculation. 1% threshold avoids near-zero theoretical intensities
        /// dominating the score with division noise.
        /// </summary>
        private const double TheoreticalIntensityFraction = 0.01;

        // ── Logistic regression weights (provisional — see class doc) ─────────
        //
        // Calibrated against test bounds B1, B2, B3 in TestDeconvolutionScorerUnit.cs:
        //   B1 (0.95, 1.5, 0.95, 0.95)   -> > 0.5  PASSES
        //   B2 (0.20, 25.0, 0.20, 0.10)  -> < 0.5  PASSES
        //   B3 range check across grid    -> all in [0, 1]  PASSES
        //
        //   linear = w_cos*cos + w_ppm*ppm + w_comp*comp + w_ratio*ratio + intercept
        //   score  = 1 / (1 + exp(linear))
        private const double WeightCosine = -3.4494;
        private const double WeightPpmError = 1.6341;
        private const double WeightCompleteness = -4.7795;
        private const double WeightRatioConsistency = -2.1883;
        private const double Intercept = -4.3142;

        // ══════════════════════════════════════════════════════════════════════
        // Public API
        // ══════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Computes the four features used by <see cref="ComputeScore"/> from an envelope's
        /// peak list and an Averagine model. Does not access the raw spectrum.
        /// </summary>
        /// <param name="envelope">Envelope to score. Must not be null.</param>
        /// <param name="model">Averagine (or other) model used for theoretical isotope positions and intensities.</param>
        /// <exception cref="ArgumentNullException">Thrown if either argument is null.</exception>
        public static EnvelopeScoreFeatures ComputeFeatures(IsotopicEnvelope envelope, AverageResidue model)
        {
            if (envelope == null) throw new ArgumentNullException(nameof(envelope));
            if (model == null) throw new ArgumentNullException(nameof(model));

            // Look up the Averagine distribution closest to this envelope's mass.
            int avgIdx = model.GetMostIntenseMassIndex(envelope.MonoisotopicMass);
            double[] rawMasses = model.GetAllTheoreticalMasses(avgIdx);
            double[] rawIntens = model.GetAllTheoreticalIntensities(avgIdx);

            // The AverageResidue arrays are sorted intensity-descending (index 0 = apex).
            // Re-sort to mass-ascending so index 0 = monoisotopic, index n = nth isotope.
            int nIso = rawMasses.Length;
            var indices = Enumerable.Range(0, nIso)
                .OrderBy(i => rawMasses[i])
                .ToArray();
            double[] theorMassesAsc = new double[nIso];
            double[] theorIntensAsc = new double[nIso];
            for (int k = 0; k < nIso; k++)
            {
                theorMassesAsc[k] = rawMasses[indices[k]];
                theorIntensAsc[k] = rawIntens[indices[k]];
            }

            // Build the observed vector aligned to the theoretical isotope index grid.
            int absCharge = Math.Abs(envelope.Charge);
            if (absCharge == 0 || envelope.Peaks == null || envelope.Peaks.Count == 0)
            {
                // Degenerate envelope: no peak data to score. Return finite defaults.
                return new EnvelopeScoreFeatures(0.0, DefaultPpmErrorForEmpty, 0.0, 0.0);
            }

            // Theoretical m/z for the nth isotope at this charge state. Uses the signed charge
            // for ToMz so the proton sign is correct; isotope step uses |charge|.
            double isotopeStepMz = Constants.C13MinusC12 / absCharge;
            double monoMz = envelope.MonoisotopicMass.ToMz(envelope.Charge);

            // Build aligned observed vector: observed[k] = intensity of the peak whose m/z is
            // within MatchTolerancePpm of the theoretical m/z for isotope index k. If no peak
            // is within tolerance, observed[k] = 0. This matching-by-window behaviour is what
            // gives shifted envelopes a low cosine (test A5).
            double[] observed = new double[nIso];
            for (int k = 0; k < nIso; k++)
            {
                double theorMz = monoMz + k * isotopeStepMz;
                double tolMz = Math.Abs(theorMz) * MatchTolerancePpm * 1e-6;

                double bestIntensity = 0.0;
                foreach (var (mz, intensity) in envelope.Peaks)
                {
                    if (Math.Abs(mz - theorMz) <= tolMz && intensity > bestIntensity)
                        bestIntensity = intensity;
                }
                observed[k] = bestIntensity;
            }

            // Cosine of the aligned vectors. Already in [0, 1] for non-negative input.
            double cosine = SpectralSimilarity.CosineOfAlignedVectors(observed, theorIntensAsc);
            if (cosine < 0.0) cosine = 0.0;
            if (cosine > 1.0) cosine = 1.0;

            // Peak completeness: of the Averagine isotopes above the 1% intensity floor,
            // how many have a matching observed peak? Uses the alignment built above.
            double theorMax = theorIntensAsc.Max();
            double cutoff = theorMax * TheoreticalIntensityFraction;
            int expectedPeaks = 0;
            int observedPeaks = 0;
            for (int k = 0; k < nIso; k++)
            {
                if (theorIntensAsc[k] > cutoff)
                {
                    expectedPeaks++;
                    if (observed[k] > 0.0) observedPeaks++;
                }
            }
            double completeness = expectedPeaks == 0 ? 0.0 : (double)observedPeaks / expectedPeaks;
            if (completeness < 0.0) completeness = 0.0;
            if (completeness > 1.0) completeness = 1.0;

            // Average ppm error: anchor at the apex peak (the most-intense observed peak) and,
            // for each observed peak, compute the integer isotope index it corresponds to and
            // its theoretical m/z relative to the apex. Average absolute ppm deviation.
            double avgPpmError = ComputeAvgPpmError(envelope, isotopeStepMz);

            // Intensity ratio consistency: if the observed envelope is a real isotope pattern,
            // observed[k] / theorIntensAsc[k] should be approximately constant across all
            // matched positions. Coefficient of variation of those ratios → 1 / (1 + CV²).
            double ratioConsistency = ComputeRatioConsistency(observed, theorIntensAsc, cutoff);

            return new EnvelopeScoreFeatures(cosine, avgPpmError, completeness, ratioConsistency);
        }

        /// <summary>
        /// Combines the four features into a single quality score in [0, 1] using a logistic
        /// transform of a linear combination. Higher score = higher-quality envelope.
        /// </summary>
        public static double ComputeScore(EnvelopeScoreFeatures features)
        {
            double linear = WeightCosine * features.AveragineCosineSimilarity
                          + WeightPpmError * features.AvgPpmError
                          + WeightCompleteness * features.PeakCompleteness
                          + WeightRatioConsistency * features.IntensityRatioConsistency
                          + Intercept;

            // Standard logistic. Guard against overflow at large positive linear values.
            if (linear >= 700.0) return 0.0;
            if (linear <= -700.0) return 1.0;
            return 1.0 / (1.0 + Math.Exp(linear));
        }

        /// <summary>
        /// Convenience: <see cref="ComputeFeatures"/> followed by <see cref="ComputeScore"/>.
        /// Returns the generic score in [0, 1]. Does not mutate the envelope — to stash the
        /// score, call <see cref="IsotopicEnvelope.SetGenericScore"/> on the result.
        /// </summary>
        public static double ScoreEnvelope(IsotopicEnvelope envelope, AverageResidue model)
        {
            var features = ComputeFeatures(envelope, model);
            return ComputeScore(features);
        }

        /// <summary>
        /// Lazily yields (envelope, score) pairs for each non-null envelope in the input. Null
        /// elements within the sequence are silently skipped — they are not an error condition,
        /// just nothing to score. Argument-null violations on the sequence or model itself are
        /// thrown eagerly.
        /// </summary>
        /// <exception cref="ArgumentNullException">
        /// Thrown when <paramref name="envelopes"/> or <paramref name="model"/> is null.
        /// </exception>
        public static IEnumerable<(IsotopicEnvelope Envelope, double Score)> ScoreEnvelopes(
            IEnumerable<IsotopicEnvelope> envelopes, AverageResidue model)
        {
            // Eager validation: tests expect ArgumentNullException to surface from the first
            // MoveNext on the result. Doing the checks before delegating to the iterator
            // ensures a useful ParamName and a fail-fast contract.
            if (envelopes == null) throw new ArgumentNullException(nameof(envelopes));
            if (model == null) throw new ArgumentNullException(nameof(model));
            return ScoreEnvelopesIterator(envelopes, model);
        }

        // ══════════════════════════════════════════════════════════════════════
        // Internals
        // ══════════════════════════════════════════════════════════════════════

        private static IEnumerable<(IsotopicEnvelope Envelope, double Score)> ScoreEnvelopesIterator(
            IEnumerable<IsotopicEnvelope> envelopes, AverageResidue model)
        {
            foreach (var envelope in envelopes)
            {
                if (envelope == null) continue; // skip nulls within the sequence
                yield return (envelope, ScoreEnvelope(envelope, model));
            }
        }

        /// <summary>
        /// Computes the average absolute ppm error of observed peaks vs the nearest theoretical
        /// isotope position. Anchors at the apex (most-intense) peak rather than the
        /// monoisotopic peak so that envelopes which truly start at M+1 or later still produce
        /// meaningful errors. Returns <see cref="DefaultPpmErrorForEmpty"/> when there are no
        /// peaks.
        /// </summary>
        private static double ComputeAvgPpmError(IsotopicEnvelope envelope, double isotopeStepMz)
        {
            if (envelope.Peaks.Count == 0) return DefaultPpmErrorForEmpty;

            double apexMz = envelope.Peaks.MaxBy(p => p.intensity).mz;
            double totalPpm = 0.0;
            int counted = 0;
            foreach (var (mz, _) in envelope.Peaks)
            {
                int n = (int)Math.Round((mz - apexMz) / isotopeStepMz);
                double theorMz = apexMz + n * isotopeStepMz;
                if (theorMz == 0.0) continue; // protect against div-by-zero on degenerate input
                totalPpm += Math.Abs(mz - theorMz) / Math.Abs(theorMz) * 1e6;
                counted++;
            }
            return counted == 0 ? DefaultPpmErrorForEmpty : totalPpm / counted;
        }

        /// <summary>
        /// Coefficient-of-variation–based consistency of observed/theoretical intensity ratios.
        /// Considers only positions where both observed and theoretical exceed the cutoff.
        /// Returns 0 if fewer than two such positions exist (cannot compute variance).
        /// </summary>
        private static double ComputeRatioConsistency(
            double[] observed, double[] theoretical, double theorCutoff)
        {
            var ratios = new List<double>();
            for (int k = 0; k < observed.Length; k++)
            {
                if (theoretical[k] > theorCutoff && observed[k] > 0.0)
                    ratios.Add(observed[k] / theoretical[k]);
            }

            if (ratios.Count < 2) return 0.0;

            double mean = 0.0;
            for (int i = 0; i < ratios.Count; i++) mean += ratios[i];
            mean /= ratios.Count;
            if (mean <= 0.0) return 0.0;

            double sumSq = 0.0;
            for (int i = 0; i < ratios.Count; i++)
            {
                double d = ratios[i] - mean;
                sumSq += d * d;
            }
            double variance = sumSq / ratios.Count;
            double cv = Math.Sqrt(variance) / mean;
            return 1.0 / (1.0 + cv * cv);
        }
    }
}
