using System.Globalization;

namespace MassSpectrometry
{
    /// <summary>
    /// Computed feature vector for a single <see cref="IsotopicEnvelope"/>. Produced by
    /// <see cref="DeconvolutionScorer.ComputeFeatures"/> and consumed by
    /// <see cref="DeconvolutionScorer.ComputeScore"/>. Pure data carrier — no scoring logic.
    /// </summary>
    public readonly struct EnvelopeScoreFeatures
    {
        /// <summary>
        /// Cosine similarity between the observed isotope intensity vector and the theoretical
        /// Averagine distribution at the envelope's mass. Range [0, 1]. Higher = better match.
        /// </summary>
        public readonly double AveragineCosineSimilarity;

        /// <summary>
        /// Mean absolute ppm error of observed peaks vs their theoretical m/z positions, anchored
        /// at the apex peak. Lower = better mass accuracy.
        /// </summary>
        public readonly double AvgPpmError;

        /// <summary>
        /// Fraction of expected Averagine isotope peaks (above 1% of theoretical max intensity)
        /// that were actually observed. Range [0, 1]. Higher = more complete envelope.
        /// </summary>
        public readonly double PeakCompleteness;

        /// <summary>
        /// Uniformity of the per-peak observed/theoretical intensity ratios, expressed as
        /// 1 / (1 + CV²) where CV is the coefficient of variation of those ratios. Range [0, 1].
        /// A real isotope envelope is the Averagine pattern scaled by a single abundance factor,
        /// so observed/theoretical should be a constant across peaks (low CV → value near 1.0).
        /// A noise envelope has erratic ratios (high CV → value near 0.0).
        /// </summary>
        public readonly double IntensityRatioConsistency;

        /// <summary>
        /// Constructs an <see cref="EnvelopeScoreFeatures"/> with explicit values. Intended for
        /// test code; production code should call <see cref="DeconvolutionScorer.ComputeFeatures"/>.
        /// </summary>
        public EnvelopeScoreFeatures(
            double averagineCosineSimilarity,
            double avgPpmError,
            double peakCompleteness,
            double intensityRatioConsistency)
        {
            AveragineCosineSimilarity = averagineCosineSimilarity;
            AvgPpmError = avgPpmError;
            PeakCompleteness = peakCompleteness;
            IntensityRatioConsistency = intensityRatioConsistency;
        }

        public override string ToString()
        {
            return string.Format(
                CultureInfo.InvariantCulture,
                "Cosine={0:F4}, PpmError={1:F4}, Completeness={2:F4}, RatioConsistency={3:F4}",
                AveragineCosineSimilarity, AvgPpmError, PeakCompleteness, IntensityRatioConsistency);
        }
    }
}