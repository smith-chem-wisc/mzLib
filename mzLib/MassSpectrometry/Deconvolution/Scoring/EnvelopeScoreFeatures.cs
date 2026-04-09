namespace MassSpectrometry
{
    /// <summary>
    /// The computed feature vector for a single <see cref="IsotopicEnvelope"/>,
    /// produced by <see cref="DeconvolutionScorer.ComputeFeatures"/> and consumed
    /// by <see cref="DeconvolutionScorer.ComputeScore"/>.
    ///
    /// All three features are computable from the envelope's peak list and an
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
        /// Constructs an <see cref="EnvelopeScoreFeatures"/> with all three feature values.
        /// </summary>
        public EnvelopeScoreFeatures(
            double averagineCosineSimilarity,
            double avgPpmError,
            double peakCompleteness)
        {
            AveragineCosineSimilarity = averagineCosineSimilarity;
            AvgPpmError               = avgPpmError;
            PeakCompleteness          = peakCompleteness;
        }

        /// <inheritdoc/>
        public override string ToString()
            => $"Cosine={AveragineCosineSimilarity:F4}  AvgPpmErr={AvgPpmError:F2}  Completeness={PeakCompleteness:F4}";
    }
}
