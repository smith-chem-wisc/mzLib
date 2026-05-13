using System.Collections.Generic;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// Deconvolution algorithm that yields envelopes for pre-deconvoluted MS1 features
    /// loaded by <see cref="FromFileDeconvolutionParameters"/> from a feature file
    /// (FlashDeconv / TopFD <c>.ms1.feature</c>, Dinosaur <c>.feature.tsv</c>, ...).
    /// </summary>
    /// <remarks>
    /// Pure join. For each <see cref="ISingleChargeMs1Feature"/> in the parameters,
    /// yields a synthetic <see cref="IsotopicEnvelope"/> when:
    /// <list type="bullet">
    /// <item><description>The feature's m/z falls inside <c>[range.Minimum, range.Maximum]</c>.</description></item>
    /// <item><description>The feature's RT window
    /// <c>[RetentionTimeStart, RetentionTimeEnd]</c> overlaps
    /// <c>[range.MinimumRt, range.MaximumRt]</c>.</description></item>
    /// <item><description>The feature's <see cref="ISingleChargeMs1Feature.Charge"/>
    /// falls inside <c>[MinAssumedChargeState, MaxAssumedChargeState]</c>.</description></item>
    /// </list>
    /// The <see cref="MzSpectrum"/> argument is ignored — the masses come from the
    /// pre-deconvoluted features, not from the spectrum.
    /// </remarks>
    internal class FromFileDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        internal FromFileDeconvolutionAlgorithm(DeconvolutionParameters deconParameters)
            : base(deconParameters)
        {
        }

        protected override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var fromFileParams = DeconvolutionParameters as FromFileDeconvolutionParameters
                ?? throw new MzLibException(
                    "FromFileDeconvolutionAlgorithm requires FromFileDeconvolutionParameters");

            var rangeWithRt = range as MzRtRange
                ?? throw new MzLibException(
                    "FromFileDeconvolutionAlgorithm requires an MzRtRange (m/z range + RT bounds)");

            // Features are sorted by m/z. Binary-search for the lower bound, iterate
            // forward, and break as soon as we pass the upper bound — O(log N + k).
            var features = fromFileParams.Features;
            int start = fromFileParams.FindFirstIndexAtOrAbove(rangeWithRt.Minimum);
            for (int i = start; i < features.Count; i++)
            {
                var feat = features[i];
                if (feat.Mz > rangeWithRt.Maximum) yield break;

                if (feat.Charge < fromFileParams.MinAssumedChargeState) continue;
                if (feat.Charge > fromFileParams.MaxAssumedChargeState) continue;
                if (feat.RetentionTimeEnd < rangeWithRt.MinimumRt) continue;
                if (feat.RetentionTimeStart > rangeWithRt.MaximumRt) continue;

                yield return new IsotopicEnvelope(
                    feat.Mz.ToMass(feat.Charge),
                    feat.Intensity,
                    feat.Charge);
            }
        }
    }
}
