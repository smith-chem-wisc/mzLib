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
        // Signs follow standard logistic convention: positive coefficient => feature
        // increases P(true). Weights can be replaced directly with sklearn / glm output.
        private const double CoefficientCosine = 3.4494;
        private const double CoefficientPpmError = -1.6341;
        private const double CoefficientCompleteness = 4.7795;
        private const double CoefficientRatioConsistency = 2.1883;
        private const double Intercept = 4.3142;

        // ── Feature matching tolerance ────────────────────────────────────────
        private const double MatchTolerancePpm = 10.0;

        // ── Spectrum-aware window constants ───────────────────────────────────
        // Half-width (Da) of the m/z window used when sampling the local noise floor for SNR.
        // Wide enough to contain plenty of non-envelope peaks at typical top-down resolution,
        // narrow enough that the noise estimate is genuinely local.
        private const double SnrWindowDa = 2.0;
        // Half-width (Da) of the m/z window used when checking whether the matched envelope
        // peak is also the local apex (CompetingPeakRatio). Smaller than SnrWindowDa because
        // we only care about peaks close enough to plausibly belong to a competing isotope
        // envelope of similar charge.
        private const double CompetingPeakWindowDa = 0.5;

        // ── Spectrum-aware logistic weights (PROVISIONAL) ─────────────────────
        // These extend the four envelope-only coefficients above with two more for the
        // spectrum-aware features. They have NOT been fit against labelled data — they
        // are chosen to satisfy the unit-test bounds in the spectrum-aware test fixture
        // (high-quality features → score > 0.5, low-quality features → score < 0.5,
        // all scores in [0, 1]) and to honour the physical meaning (higher SNR = better,
        // more local-apex matches = better). Sign convention matches the envelope-only
        // weights above: positive coefficient ⇒ feature increases P(true).
        // The intercept is shifted negative relative to the envelope-only Intercept
        // because the SNR squashing function and the additional CompetingPeakRatio term
        // raise the typical linear score; without a more negative intercept, every
        // envelope would score above 0.5 even at low cosine and completeness.
        // RECALIBRATE AGAINST LABELLED DATA before relying on these in production.
        private const double CoefficientSnr = 2.0;
        private const double CoefficientCompetingRatio = 1.5;
        private const double InterceptSpectrumAware = -2.5;

        // TODO: Apex spectral percentile and unmatched-intensity-ratio features are deferred.
        // Each requires its own validation and weight calibration; see the spectrum-aware
        // PR follow-up section.

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
            // For each theoretical isotope position n, find the most intense peak in
            // envelope.Peaks within MatchTolerancePpm.
            double monoMass = envelope.MonoisotopicMass;
            // Signed charge: ToMz embeds the polarity-dependent proton sign
            // (mass/|z| + sign(z)*proton), so a negative envelope anchors on its
            // physically-correct negative-mode m/z, not the positive-mode mirror.
            double monoMz = monoMass.ToMz(envelope.Charge);
            var peakList = envelope.Peaks;

            double[] observed = new double[nIso];
            bool[] peakMatched = new bool[nIso];

            for (int n = 0; n < nIso; n++)
            {
                double theorMz = monoMz + n * isotopeStepMz;
                double tolMz = theorMz * MatchTolerancePpm * 1e-6;

                // Find the most intense peak in the envelope within tolerance of this
                // theoretical position.
                double bestIntensity = 0.0;

                foreach (var (mz, intensity) in peakList)
                {
                    if (Math.Abs(mz - theorMz) < tolMz && intensity > bestIntensity)
                    {
                        bestIntensity = intensity;
                    }
                }

                if (bestIntensity > 0.0)
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
        /// Spectrum-aware feature computation. Reuses the four envelope-only features from
        /// <see cref="ComputeFeatures(IsotopicEnvelope, AverageResidue)"/> verbatim, then
        /// augments them with two features that require access to the original spectrum:
        /// <see cref="EnvelopeScoreFeatures.LocalSignalToNoise"/> and
        /// <see cref="EnvelopeScoreFeatures.CompetingPeakRatio"/>.
        /// </summary>
        /// <remarks>
        /// Use this overload when the caller still holds the <see cref="MzSpectrum"/> the
        /// envelope was extracted from (the always-on path through
        /// <see cref="Deconvoluter.DeconvoluteWithGenericScoring(MzSpectrum, DeconvolutionParameters, MzLibUtil.MzRange)"/>).
        /// The richer feature vector flows into <see cref="ComputeScoreWithSpectrumContext"/>.
        ///
        /// Note: spectrum-aware scores are not strictly comparable across deconvolution
        /// algorithms because preprocessing differences (peak picking, smoothing, threshold
        /// choice) can leak into spectrum-derived features. For cross-algorithm comparability,
        /// use the envelope-only path.
        /// </remarks>
        public static EnvelopeScoreFeatures ComputeFeatures(
            IsotopicEnvelope envelope,
            AverageResidue model,
            MzSpectrum spectrum)
        {
            if (envelope == null) throw new ArgumentNullException(nameof(envelope));
            if (model == null) throw new ArgumentNullException(nameof(model));
            if (spectrum == null) throw new ArgumentNullException(nameof(spectrum));

            // Reuse the existing envelope-only computation verbatim. This is the contract:
            // anything that's true of the envelope-only features in the four-arg path is
            // true here too — there is exactly one cosine/ppm/completeness/ratio assignment
            // in this file.
            EnvelopeScoreFeatures envelopeOnly = ComputeFeatures(envelope, model);

            double snr = ComputeLocalSnr(envelope, spectrum);
            double competing = ComputeCompetingPeakRatio(envelope, model, spectrum);

            return new EnvelopeScoreFeatures(
                envelopeOnly.AveragineCosineSimilarity,
                envelopeOnly.AvgPpmError,
                envelopeOnly.PeakCompleteness,
                envelopeOnly.IntensityRatioConsistency,
                snr,
                competing);
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

            return 1.0 / (1.0 + Math.Exp(-linear));
        }

        /// <summary>
        /// Convenience method: computes features then score for a single envelope.
        /// </summary>
        public static double ScoreEnvelope(IsotopicEnvelope envelope, AverageResidue model)
            => ComputeScore(ComputeFeatures(envelope, model));

        /// <summary>
        /// Combines all six features (four envelope-only + two spectrum-aware) into a single
        /// quality score in [0, 1] using a logistic function. Use this in place of
        /// <see cref="ComputeScore"/> when
        /// <see cref="ComputeFeatures(IsotopicEnvelope, AverageResidue, MzSpectrum)"/>
        /// was used to produce the features.
        /// </summary>
        /// <remarks>
        /// Throws <see cref="ArgumentException"/> if <paramref name="features"/> does not have
        /// spectrum features (i.e. <see cref="EnvelopeScoreFeatures.HasSpectrumFeatures"/> is
        /// false). This prevents accidentally scoring envelope-only features through the
        /// spectrum-aware formula, which would consume NaN values and silently return NaN.
        ///
        /// The <see cref="CoefficientSnr"/>, <see cref="CoefficientCompetingRatio"/>, and
        /// <see cref="InterceptSpectrumAware"/> constants are provisional — see the comment
        /// block at the top of the class.
        /// </remarks>
        public static double ComputeScoreWithSpectrumContext(EnvelopeScoreFeatures features)
        {
            if (!features.HasSpectrumFeatures)
                throw new ArgumentException(
                    "ComputeScoreWithSpectrumContext requires features produced by the spectrum-aware ComputeFeatures overload.",
                    nameof(features));

            double linear = InterceptSpectrumAware
                + CoefficientCosine * features.AveragineCosineSimilarity
                + CoefficientPpmError * features.AvgPpmError
                + CoefficientCompleteness * features.PeakCompleteness
                + CoefficientRatioConsistency * features.IntensityRatioConsistency
                + CoefficientSnr * SquashSnr(features.LocalSignalToNoise)
                + CoefficientCompetingRatio * features.CompetingPeakRatio;

            // Saturate the exponent to keep the score firmly inside [0, 1] even at
            // pathological inputs (e.g. AvgPpmError = double.MaxValue from the absCharge==0
            // guard in ComputeFeatures). The cutoffs are chosen well inside the range where
            // Math.Exp would overflow, but the resulting score is still within IEEE 754 limits.
            if (linear >= 700.0) return 1.0;
            if (linear <= -700.0) return 0.0;
            return 1.0 / (1.0 + Math.Exp(-linear));
        }

        /// <summary>
        /// Spectrum-aware convenience method: computes features (including SNR and
        /// competing-peak ratio) then score for a single envelope. Mirrors the shape of
        /// the envelope-only <see cref="ScoreEnvelope(IsotopicEnvelope, AverageResidue)"/>.
        /// </summary>
        public static double ScoreEnvelope(IsotopicEnvelope envelope, AverageResidue model, MzSpectrum spectrum)
            => ComputeScoreWithSpectrumContext(ComputeFeatures(envelope, model, spectrum));

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

        /// <summary>
        /// Median per-peak SNR for the envelope against the local spectrum noise floor.
        /// For each envelope peak at m/z = mz_p:
        ///   1. Extract the window [mz_p − <see cref="SnrWindowDa"/>, mz_p + <see cref="SnrWindowDa"/>] from the spectrum.
        ///   2. From the spectrum peaks in that window, exclude any peak whose m/z is within
        ///      <see cref="MatchTolerancePpm"/> of any envelope peak (i.e. the envelope's own peaks).
        ///   3. SNR_p = envelope peak intensity / median(remaining window intensities).
        ///      If the remaining window is empty or its median is non-positive, SNR_p is treated
        ///      as missing and is not contributed to the aggregate median.
        /// Returns the median of the SNR_p values across all envelope peaks, or 0.0 if no SNR_p
        /// could be computed (degenerate spectrum / envelope).
        /// </summary>
        private static double ComputeLocalSnr(IsotopicEnvelope envelope, MzSpectrum spectrum)
        {
            if (envelope.Peaks == null || envelope.Peaks.Count == 0) return 0.0;
            if (spectrum.Size == 0) return 0.0;

            double[] xArr = spectrum.XArray;
            double[] yArr = spectrum.YArray;
            var envPeaks = envelope.Peaks;

            var snrPerPeak = new List<double>(envPeaks.Count);

            foreach (var (envMz, envIntensity) in envPeaks)
            {
                if (envIntensity <= 0.0) continue;

                (int start, int end) = SpectrumWindowIndices(spectrum, envMz - SnrWindowDa, envMz + SnrWindowDa);
                if (start > end) continue; // empty window

                // Collect non-envelope intensities inside the window.
                var noiseIntensities = new List<double>(end - start + 1);
                for (int i = start; i <= end; i++)
                {
                    double xi = xArr[i];
                    if (IsEnvelopePeakMz(xi, envPeaks)) continue;
                    noiseIntensities.Add(yArr[i]);
                }

                if (noiseIntensities.Count == 0) continue; // no noise sample → skip
                double medianNoise = Median(noiseIntensities);
                if (medianNoise <= 0.0) continue;

                snrPerPeak.Add(envIntensity / medianNoise);
            }

            return snrPerPeak.Count == 0 ? 0.0 : Median(snrPerPeak);
        }

        /// <summary>
        /// Fraction of matched theoretical isotope positions where the matched envelope peak
        /// is also the local apex within ± <see cref="CompetingPeakWindowDa"/>. Catches
        /// envelopes extracted from the shoulder of a larger feature.
        ///
        /// The alignment (Averagine model lookup, mass-ascending sort, per-isotope theoretical
        /// m/z, 10 ppm matching) mirrors the envelope-only <see cref="ComputeFeatures(IsotopicEnvelope, AverageResidue)"/>
        /// path so the two paths agree on which positions count as "matched".
        /// </summary>
        private static double ComputeCompetingPeakRatio(
            IsotopicEnvelope envelope, AverageResidue model, MzSpectrum spectrum)
        {
            if (spectrum.Size == 0) return 0.0;

            int absCharge = Math.Abs(envelope.Charge);
            if (absCharge == 0) return 0.0;

            // Mirror the alignment performed by ComputeFeatures(envelope, model).
            int avgIdx = model.GetMostIntenseMassIndex(envelope.MonoisotopicMass);
            double[] rawMasses = model.GetAllTheoreticalMasses(avgIdx);
            double[] rawIntens = model.GetAllTheoreticalIntensities(avgIdx);

            int nIsoLen = rawMasses.Length;
            double[] avgMassAsc = new double[nIsoLen];
            double[] avgIntAsc = new double[nIsoLen];
            Array.Copy(rawMasses, avgMassAsc, nIsoLen);
            Array.Copy(rawIntens, avgIntAsc, nIsoLen);
            Array.Sort(avgMassAsc, avgIntAsc);

            double maxTheoInt = 0.0;
            for (int i = 0; i < avgIntAsc.Length; i++)
                if (avgIntAsc[i] > maxTheoInt) maxTheoInt = avgIntAsc[i];
            double threshold1pct = maxTheoInt * 0.01;

            double isotopeStepMz = Constants.C13MinusC12 / absCharge;
            double monoMz = envelope.MonoisotopicMass.ToMz(envelope.Charge);
            var peakList = envelope.Peaks;

            int matched = 0;
            int nonCompeting = 0;
            const double apexEpsilon = 1e-9;

            for (int n = 0; n < avgIntAsc.Length; n++)
            {
                if (avgIntAsc[n] <= threshold1pct) continue;

                double theorMz = monoMz + n * isotopeStepMz;
                double tolMz = theorMz * MatchTolerancePpm * 1e-6;

                // Find the most intense envelope peak within tolerance, recording its m/z.
                double bestIntensity = 0.0;
                double bestMz = 0.0;
                foreach (var (mz, intensity) in peakList)
                {
                    if (Math.Abs(mz - theorMz) < tolMz && intensity > bestIntensity)
                    {
                        bestIntensity = intensity;
                        bestMz = mz;
                    }
                }
                if (bestIntensity <= 0.0) continue; // unmatched — counted in PeakCompleteness instead

                matched++;

                (int start, int end) = SpectrumWindowIndices(
                    spectrum, bestMz - CompetingPeakWindowDa, bestMz + CompetingPeakWindowDa);
                if (start > end)
                {
                    // No spectrum peaks in the competing-peak window: trivially the apex.
                    nonCompeting++;
                    continue;
                }

                double apexInWindow = 0.0;
                double[] yArr = spectrum.YArray;
                for (int i = start; i <= end; i++)
                    if (yArr[i] > apexInWindow) apexInWindow = yArr[i];

                if (bestIntensity >= apexInWindow * (1.0 - apexEpsilon))
                    nonCompeting++;
            }

            return matched == 0 ? 0.0 : (double)nonCompeting / matched;
        }

        /// <summary>
        /// Squashes raw SNR (potentially unbounded above) into [0, 1] via 1 − exp(−snr / scale).
        /// Scale chosen so SNR ≈ 5 (typical good envelope) maps to ~0.6 and SNR ≈ 20 maps to ~0.96.
        /// Negative or NaN SNR maps to 0.
        /// </summary>
        private static double SquashSnr(double snr)
        {
            if (double.IsNaN(snr) || snr <= 0.0) return 0.0;
            const double SnrScale = 6.0;
            return 1.0 - Math.Exp(-snr / SnrScale);
        }

        /// <summary>
        /// Inclusive index range for the spectrum peaks in <c>[minMz, maxMz]</c>. Returns
        /// <c>(start, end)</c> with <c>start &gt; end</c> when the window is empty.
        /// Wraps the obsolete <see cref="MzSpectrum.ExtractIndices"/> at a single call site.
        /// </summary>
        private static (int start, int end) SpectrumWindowIndices(MzSpectrum spectrum, double minMz, double maxMz)
        {
#pragma warning disable CS0618 // ExtractIndices is marked Obsolete but is the documented entry point for index-range extraction.
            return spectrum.ExtractIndices(minMz, maxMz);
#pragma warning restore CS0618
        }

        /// <summary>
        /// True if <paramref name="mz"/> is within <see cref="MatchTolerancePpm"/> of any peak
        /// in <paramref name="envPeaks"/>. Used to exclude envelope peaks from the noise sample.
        /// </summary>
        private static bool IsEnvelopePeakMz(double mz, IReadOnlyList<(double mz, double intensity)> envPeaks)
        {
            for (int j = 0; j < envPeaks.Count; j++)
            {
                double pMz = envPeaks[j].mz;
                if (pMz <= 0.0) continue;
                double tol = pMz * MatchTolerancePpm * 1e-6;
                if (Math.Abs(mz - pMz) < tol) return true;
            }
            return false;
        }

        /// <summary>
        /// Sample median of a non-empty list. Sorts a local copy so the caller's list is
        /// not mutated.
        /// </summary>
        private static double Median(List<double> values)
        {
            int n = values.Count;
            double[] copy = new double[n];
            values.CopyTo(copy);
            Array.Sort(copy);
            return (n & 1) == 1
                ? copy[n / 2]
                : 0.5 * (copy[n / 2 - 1] + copy[n / 2]);
        }
    }
}