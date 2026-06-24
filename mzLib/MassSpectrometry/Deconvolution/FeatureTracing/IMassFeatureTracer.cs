using System.Collections.Generic;
using MassSpectrometry.Deconvolution.Consensus;

namespace MassSpectrometry.Deconvolution.FeatureTracing
{
    /// <summary>
    /// Aggregates per-MS1-scan deconvoluted isotopic envelopes (from any deconvolution
    /// algorithm) into whole-file, cross-scan <see cref="MassFeature"/>s.
    ///
    /// <para>
    /// This is the pluggable boundary that lets MetaFlashDecon emit an <c>_ms1.feature</c>
    /// file: the per-spectrum deconvolution step is decoupled from the feature-tracing
    /// strategy, so different tracers can be swapped and compared on the <i>same</i> input
    /// envelopes. Two strategies are intended:
    /// </para>
    /// <list type="bullet">
    ///   <item><see cref="ConsensusMassFeatureTracer"/> — charge-locked traces stitched
    ///   across charge (the consensus mass-tracing pipeline).</item>
    ///   <item>a FLASHDeconv-style tracer (Phase 3) — neutral-mass mass-trace detection with
    ///   charge re-attached afterward, faithful to OpenMS FLASHDeconv.</item>
    /// </list>
    /// Both produce <see cref="MassFeature"/>s, which a caller writes to an
    /// <c>_ms1.feature</c> file via <c>Ms1FeatureFile.FromMassFeatures(...)</c>.
    /// </summary>
    public interface IMassFeatureTracer
    {
        /// <summary>
        /// Traces whole-file mass features from per-scan deconvolution output.
        /// </summary>
        /// <param name="ms1Scans">MS1 scans, in acquisition (scan-number) order.</param>
        /// <param name="perScanEnvelopes">
        /// Deconvoluted envelopes for each scan, aligned 1:1 and in the same order as
        /// <paramref name="ms1Scans"/>.
        /// </param>
        /// <returns>The cross-scan mass features detected over the whole file.</returns>
        List<MassFeature> TraceFeatures(
            IReadOnlyList<MsDataScan> ms1Scans,
            IReadOnlyList<IReadOnlyList<IsotopicEnvelope>> perScanEnvelopes);
    }
}
