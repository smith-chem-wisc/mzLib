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
    /// <remarks>
    /// Kept as a readonly struct (4 doubles, 32 bytes) rather than a class or record
    /// class because <see cref="DeconvolutionScorer.ComputeFeatures"/> is invoked once
    /// per envelope on the batch-scoring hot path (potentially thousands per scan), and
    /// stack-local value semantics avoid the per-envelope heap allocation a reference
    /// type would incur. The 32-byte copy cost on parameter passes is acceptable
    /// because the type only crosses the ComputeFeatures -> ComputeScore boundary once.
    /// </remarks>
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
        /// Constructs an <see cref="EnvelopeScoreFeatures"/> with all four feature values.
        /// The <paramref name="intensityRatioConsistency"/> parameter defaults to 0.0
        /// for backward compatibility with any call sites that pass only three values.
        /// </summary>
        public EnvelopeScoreFeatures(
            double averagineCosineSimilarity,
            double avgPpmError,
            double peakCompleteness,
            double intensityRatioConsistency = 0.0)
        {
            AveragineCosineSimilarity = averagineCosineSimilarity;
            AvgPpmError = avgPpmError;
            PeakCompleteness = peakCompleteness;
            IntensityRatioConsistency = intensityRatioConsistency;
        }

        /// <inheritdoc/>
        public override string ToString()
            => $"Cosine={AveragineCosineSimilarity:F4}  AvgPpmErr={AvgPpmError:F2}  " +
               $"Completeness={PeakCompleteness:F4}  RatioConsistency={IntensityRatioConsistency:F4}";
    }
}