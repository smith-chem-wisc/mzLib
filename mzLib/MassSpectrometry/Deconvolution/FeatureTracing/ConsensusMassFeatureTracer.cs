using System.Collections.Generic;
using System.Linq;
using MassSpectrometry.Deconvolution.Consensus;

namespace MassSpectrometry.Deconvolution.FeatureTracing
{
    /// <summary>
    /// <see cref="IMassFeatureTracer"/> backed by the consensus mass-tracing pipeline:
    /// <see cref="MassTraceBuilder.BuildTraces"/> (charge-locked traces) →
    /// <see cref="TraceCorrector.Correct"/> (resolution-aware off-by-one correction) →
    /// <see cref="MassFeatureBuilder.BuildFeatures"/> (cross-charge stitching).
    /// <para>
    /// Lets MetaFlashDecon's per-scan envelopes be traced by the existing consensus grouping
    /// so it can be compared side-by-side against the FLASHDeconv-native tracer on identical
    /// input. Pure delegation: the algorithm-specific work lives in the consensus types.
    /// </para>
    /// </summary>
    public sealed class ConsensusMassFeatureTracer : IMassFeatureTracer
    {
        /// <summary>Default Pass-B trace-grouping tolerance (Da).</summary>
        public const double DefaultTraceToleranceDa = 1.5;
        /// <summary>Default maximum consecutive missing scans tolerated within a trace.</summary>
        public const int DefaultMaxGap = 1;
        /// <summary>Default cross-charge stitching mass tolerance (ppm).</summary>
        public const double DefaultCrossChargeMassPpm = 10.0;

        private readonly double _traceToleranceDa;
        private readonly int _maxGap;
        private readonly double _crossChargeMassPpm;

        /// <param name="traceToleranceDa">Pass-B trace-grouping tolerance (Da).</param>
        /// <param name="maxGap">Maximum consecutive missing scans tolerated within a trace.</param>
        /// <param name="crossChargeMassPpm">Cross-charge stitching mass tolerance (ppm).</param>
        public ConsensusMassFeatureTracer(
            double traceToleranceDa = DefaultTraceToleranceDa,
            int maxGap = DefaultMaxGap,
            double crossChargeMassPpm = DefaultCrossChargeMassPpm)
        {
            _traceToleranceDa = traceToleranceDa;
            _maxGap = maxGap;
            _crossChargeMassPpm = crossChargeMassPpm;
        }

        /// <inheritdoc />
        public List<MassFeature> TraceFeatures(
            IReadOnlyList<MsDataScan> ms1Scans,
            IReadOnlyList<IReadOnlyList<IsotopicEnvelope>> perScanEnvelopes)
        {
            var traces = MassTraceBuilder.BuildTraces(ms1Scans, perScanEnvelopes, _traceToleranceDa, _maxGap);
            // Correct() may split one trace into multiple corrected traces (e.g. a resolvable
            // deamidation co-group), so flatten with SelectMany.
            var corrected = traces.SelectMany(t => TraceCorrector.Correct(t)).ToList();
            return MassFeatureBuilder.BuildFeatures(corrected, _crossChargeMassPpm);
        }
    }
}
