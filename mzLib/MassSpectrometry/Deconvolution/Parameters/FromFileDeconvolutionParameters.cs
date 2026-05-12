#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    /// <summary>
    /// Deconvolution parameters that wrap a pre-computed list of MS1 features —
    /// typically loaded by a Readers consumer from a FlashDeconv / TopFD
    /// <c>.ms1.feature</c> file or produced by mzLib's own whole-file
    /// deconvolution. Pairing them with MS2 scans is delegated to
    /// <see cref="FromFileDeconvolutionDeconvolutionAlgorithm"/>, which the
    /// <see cref="Deconvoluter"/> factory selects when this parameter type is supplied.
    /// </summary>
    /// <remarks>
    /// The caller is responsible for loading the features. This keeps the
    /// <c>MassSpectrometry → Readers</c> dependency arrow out of the picture —
    /// callers in Readers (or downstream consumers) hand in
    /// <see cref="ISingleChargeMs1Feature"/> instances; this class simply indexes them
    /// in memory for the algorithm to filter by m/z + RT range.
    /// </remarks>
    public class FromFileDeconvolutionParameters : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; } = DeconvolutionType.FromFile;

        /// <summary>
        /// The pre-loaded per-charge features the algorithm will filter against.
        /// </summary>
        public IReadOnlyList<ISingleChargeMs1Feature> Features { get; }

        public FromFileDeconvolutionParameters(
            IEnumerable<ISingleChargeMs1Feature> features,
            int minCharge,
            int maxCharge,
            Polarity polarity = Polarity.Positive)
            : base(minCharge, maxCharge, polarity)
        {
            if (features is null) throw new ArgumentNullException(nameof(features));
            Features = features.ToList();
        }

        // No decoy support for from-file deconvolution: the masses are taken as
        // authoritative from the producer, so flipping the isotope spacing has no
        // physical meaning here.
        public override DeconvolutionParameters? ToDecoyParameters() => null;
    }
}
