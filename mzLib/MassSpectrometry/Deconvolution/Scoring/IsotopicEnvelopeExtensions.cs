namespace MassSpectrometry
{
    /// <summary>
    /// Convenience extensions for <see cref="IsotopicEnvelope"/> that compute and stash the
    /// generic deconvolution score on the envelope. Designed for downstream callers
    /// (e.g. MetaMorpheus) that hold an envelope and want a one-line way to get a
    /// comparable [0, 1] score for it.
    /// </summary>
    public static class IsotopicEnvelopeExtensions
    {
        /// <summary>
        /// Returns the envelope's <see cref="IsotopicEnvelope.GenericScore"/>, computing and
        /// caching it on first access. The score is in [0, 1] and is comparable across
        /// deconvolution algorithms because it is recomputed from the envelope's peak list
        /// using the supplied <paramref name="model"/>.
        /// </summary>
        /// <remarks>
        /// Idempotent: if <see cref="IsotopicEnvelope.GenericScore"/> has already been set
        /// (e.g. by <see cref="Deconvoluter.DeconvoluteWithGenericScoring(MzSpectrum, DeconvolutionParameters, MzLibUtil.MzRange)"/>),
        /// the cached value is returned without recomputation.
        /// </remarks>
        public static double GetOrComputeGenericScore(this IsotopicEnvelope envelope, AverageResidue model)
        {
            if (envelope == null) throw new System.ArgumentNullException(nameof(envelope));
            if (model == null) throw new System.ArgumentNullException(nameof(model));

            if (envelope.GenericScore.HasValue) return envelope.GenericScore.Value;

            double score = DeconvolutionScorer.ScoreEnvelope(envelope, model);
            envelope.SetGenericScore(score);
            return score;
        }

        /// <summary>
        /// Convenience overload that pulls the <see cref="AverageResidue"/> from the supplied
        /// <paramref name="parameters"/>. Useful when the caller already holds a parameters
        /// object (the common case in MetaMorpheus task code).
        /// </summary>
        public static double GetOrComputeGenericScore(
            this IsotopicEnvelope envelope, DeconvolutionParameters parameters)
        {
            if (parameters == null) throw new System.ArgumentNullException(nameof(parameters));
            return GetOrComputeGenericScore(envelope, parameters.AverageResidueModel);
        }
    }
}