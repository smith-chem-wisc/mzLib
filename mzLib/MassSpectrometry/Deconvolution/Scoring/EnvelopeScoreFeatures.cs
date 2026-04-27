namespace MassSpectrometry
{
    /// <summary>
    /// The computed feature vector for a single <see cref="IsotopicEnvelope"/>,
    /// produced by <see cref="DeconvolutionScorer.ComputeFeatures"/> and consumed
    /// by <see cref="DeconvolutionScorer.ComputeScore"/>.
    ///
    /// All four features are computable from the envelope's peak list and an
    /// <see cref="AverageResidue"/> model alone — no raw spectrum access is required.
    /// This allows the scorer to operate after the spectrum has been discarded
    /// (deconvolute-and-discard callers).
    /// </summary>
    public readonly struct EnvelopeScoreFeatures
    {
        /// <summary>
        /// Cosine similarity between the observed isotope intensity distribution
        /// and the theoretical Averagine distribution at the envelope's monoisotopic
        /// mass. Computed by aligning observed peak intensities to Averagine isotope
        /// positions at 10 ppm tolerance.
        /// Range: [0, 1]. Higher values indicate better agreement with the expected
        /// isotope pattern.
        /// </summary>
        public readonly double AveragineCosineSimilarity;

        /// <summary>
        /// Mean absolute ppm error of the observed peaks relative to their
        /// theoretical m/z positions, anchored at the apex (most intense) peak.
        /// Units: ppm. Lower values indicate better mass accuracy.
        /// </summary>
        public readonly double AvgPpmError;

        /// <summary>
        /// Fraction of expected Averagine isotope peaks (those with theoretical
        /// intensity above 1% of the theoretical maximum) that were actually
        /// observed in the envelope within 10 ppm.
        /// Range: [0, 1]. A value of 1.0 means all significant isotope peaks were
        /// observed; 0.5 means half were observed.
        /// </summary>
        public readonly double PeakCompleteness;

        /// <summary>
        /// Consistency of the observed-to-theoretical intensity ratios across all
        /// matched isotope peaks, expressed as <c>1 / (1 + CV²)</c> where CV is the
        /// coefficient of variation (std / mean) of the per-peak scale ratios
        /// <c>observed[n] / theoretical[n]</c>.
        ///
        /// Physical interpretation: a real isotope envelope is the Averagine pattern
        /// scaled by a single abundance factor, so every observed[n] / theoretical[n]
        /// ratio should be approximately the same constant (low CV → value near 1.0).
        /// A noise envelope has erratic ratios with no physical coherence
        /// (high CV → value near 0.0).
        ///
        /// This feature is complementary to <see cref="AveragineCosineSimilarity"/>:
        /// cosine captures global shape agreement; ratio consistency captures the
        /// uniformity of the per-peak scale errors. A noisy envelope can achieve
        /// moderate cosine similarity while having highly inconsistent ratios.
        ///
        /// Range: [0, 1]. Returns 0.0 when fewer than 2 peaks are matched (CV
        /// undefined) or when the mean ratio is zero.
        /// </summary>
        public readonly double IntensityRatioConsistency;

        /// <summary>
        /// Median signal-to-noise ratio across the envelope's matched isotope peaks. For each
        /// matched peak, SNR is the peak's intensity divided by the median intensity of all
        /// spectrum peaks in a wider m/z window around it that excludes the envelope's own
        /// peaks. NaN when the scorer was called without a spectrum.
        /// </summary>
        public readonly double LocalSignalToNoise;

        /// <summary>
        /// Fraction of theoretical isotope positions (above the 1% intensity floor) where the
        /// matched observed peak is the tallest peak within a small m/z window — i.e. the
        /// envelope is not extracted from the shoulder of a larger feature. Range [0, 1].
        /// NaN when the scorer was called without a spectrum.
        /// </summary>
        public readonly double CompetingPeakRatio;

        /// <summary>
        /// True if spectrum-aware features were computed (i.e. a spectrum was supplied to
        /// <see cref="DeconvolutionScorer.ComputeFeatures(IsotopicEnvelope, AverageResidue, MzSpectrum)"/>).
        /// Equivalent to <c>!double.IsNaN(LocalSignalToNoise)</c>.
        /// </summary>
        public bool HasSpectrumFeatures => !double.IsNaN(LocalSignalToNoise);

        /// <summary>
        /// Constructs an <see cref="EnvelopeScoreFeatures"/> with all four feature values.
        /// The <paramref name="intensityRatioConsistency"/> parameter defaults to 0.0
        /// for backward compatibility with any call sites that pass only three values.
        /// Spectrum-aware fields are set to <see cref="double.NaN"/>;
        /// <see cref="HasSpectrumFeatures"/> will be <c>false</c>.
        /// </summary>
        public EnvelopeScoreFeatures(
            double averagineCosineSimilarity,
            double avgPpmError,
            double peakCompleteness,
            double intensityRatioConsistency = 0.0)
            : this(averagineCosineSimilarity, avgPpmError, peakCompleteness,
                   intensityRatioConsistency, double.NaN, double.NaN)
        { }

        /// <summary>
        /// Constructs an <see cref="EnvelopeScoreFeatures"/> with all six feature values,
        /// including the spectrum-aware <paramref name="localSignalToNoise"/> and
        /// <paramref name="competingPeakRatio"/>. Use this constructor (or, more commonly,
        /// let <see cref="DeconvolutionScorer.ComputeFeatures(IsotopicEnvelope, AverageResidue, MzSpectrum)"/>
        /// build the struct) when both spectrum-aware features have been measured.
        /// </summary>
        public EnvelopeScoreFeatures(
            double averagineCosineSimilarity,
            double avgPpmError,
            double peakCompleteness,
            double intensityRatioConsistency,
            double localSignalToNoise,
            double competingPeakRatio)
        {
            AveragineCosineSimilarity = averagineCosineSimilarity;
            AvgPpmError = avgPpmError;
            PeakCompleteness = peakCompleteness;
            IntensityRatioConsistency = intensityRatioConsistency;
            LocalSignalToNoise = localSignalToNoise;
            CompetingPeakRatio = competingPeakRatio;
        }

        /// <inheritdoc/>
        public override string ToString()
        {
            string envOnly =
                $"Cosine={AveragineCosineSimilarity:F4}  AvgPpmErr={AvgPpmError:F2}  " +
                $"Completeness={PeakCompleteness:F4}  RatioConsistency={IntensityRatioConsistency:F4}";
            return HasSpectrumFeatures
                ? envOnly + $"  SNR={LocalSignalToNoise:F2}  CompetingPeakRatio={CompetingPeakRatio:F4}"
                : envOnly;
        }
    }
}