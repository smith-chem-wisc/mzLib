using Chemistry;

namespace MassSpectrometry
{
    /// <summary>
    /// Creates decoy <see cref="DeconvolutionParameters"/> instances for target-decoy
    /// FDR estimation.
    ///
    /// The decoy strategy runs the same deconvolution algorithm on the same spectrum
    /// with an isotope spacing that is physically impossible for any real molecule.
    /// The canonical decoy distance is 0.9444 Da — taken directly from OpenMS
    /// FLASHDeconv <c>SpectralDeconvolution.cpp</c> (<c>noise_iso_delta_ = 0.9444</c>).
    /// This value is close to the true C13-C12 spacing (1.003355 Da) but does not
    /// correspond to any real charge-state isotope series, so decoy envelopes that
    /// survive deconvolution represent false-positive deconvolution noise.
    ///
    /// The decoy parameters are a shallow clone of the target parameters with only
    /// <see cref="DeconvolutionParameters.DecoyIsotopeDistance"/> and
    /// <see cref="DeconvolutionParameters.IsDecoyRun"/> changed. All other settings
    /// (charge range, mass range, tolerance, Averagine model) are identical, so decoy
    /// envelopes face exactly the same quality filters as targets.
    /// </summary>
    internal static class DeconvolutionDecoyGenerator
    {
        /// <summary>
        /// The shifted isotope spacing used for decoy deconvolution runs (Da).
        /// Value from OpenMS FLASHDeconv <c>SpectralDeconvolution.cpp</c>:
        /// <c>noise_iso_delta_ = 0.9444</c>.
        /// </summary>
        internal const double DecoyIsotopeDistance = 0.9444;

        /// <summary>
        /// Creates a shallow clone of <paramref name="parameters"/> suitable for a
        /// decoy deconvolution pass.
        ///
        /// The returned object is the same concrete subtype as the input, so the
        /// <see cref="Deconvoluter"/> factory will route to the same algorithm.
        /// Only <see cref="DeconvolutionParameters.DecoyIsotopeDistance"/> and
        /// <see cref="DeconvolutionParameters.IsDecoyRun"/> differ from the original.
        /// </summary>
        /// <param name="parameters">Target parameters to clone.</param>
        /// <returns>
        /// A new parameters instance of the same concrete type with decoy settings applied.
        /// </returns>
        internal static DeconvolutionParameters MakeDecoyParameters(DeconvolutionParameters parameters)
        {
            // ShallowClone is virtual — subclasses that cache derived state (e.g.
            // IsoDecDeconvolutionParameters caches IsoSettings) override it to
            // invalidate that cache on the clone.
            DeconvolutionParameters clone = parameters.ShallowClone();
            clone.DecoyIsotopeDistance = DecoyIsotopeDistance;
            clone.IsDecoyRun           = true;
            return clone;
        }
    }
}
