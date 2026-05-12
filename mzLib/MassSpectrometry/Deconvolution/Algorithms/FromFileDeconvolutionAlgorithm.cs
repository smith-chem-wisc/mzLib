using System.Collections.Generic;
using Chemistry;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Deconvolution algorithm that yields envelopes for pre-deconvoluted MS1 features
    /// (typically loaded from a FlashDeconv / TopFD <c>.ms1.feature</c> file via the
    /// Readers project, then wrapped in <see cref="FromFileDeconvolutionParameters"/>).
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
    /// pre-deconvoluted features, not from the spectrum. The returned envelope is
    /// built via <see cref="IsotopicEnvelope(double, double, int)"/>; its <c>Peaks</c>
    /// list contains a single synthetic entry at the feature's m/z. Per-peak data
    /// from the producer is not surfaced because external <c>.ms1.feature</c> formats
    /// do not carry it.
    /// </remarks>
    internal class FromFileDeconvolutionAlgorithm : DeconvolutionAlgorithm
    {
        internal FromFileDeconvolutionAlgorithm(DeconvolutionParameters deconParameters)
            : base(deconParameters)
        {
        }

        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var fromFileParams = DeconvolutionParameters as FromFileDeconvolutionParameters
                ?? throw new MzLibException(
                    "FromFileDeconvolutionAlgorithm requires FromFileDeconvolutionParameters");

            var rangeWithRt = range as MzRtRange
                ?? throw new MzLibException(
                    "FromFileDeconvolutionAlgorithm requires an MzRtRange (m/z range + RT bounds)");

            foreach (var feat in fromFileParams.Features)
            {
                if (feat.Charge < fromFileParams.MinAssumedChargeState) continue;
                if (feat.Charge > fromFileParams.MaxAssumedChargeState) continue;
                if (feat.Mz < rangeWithRt.Minimum) continue;
                if (feat.Mz > rangeWithRt.Maximum) continue;
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
