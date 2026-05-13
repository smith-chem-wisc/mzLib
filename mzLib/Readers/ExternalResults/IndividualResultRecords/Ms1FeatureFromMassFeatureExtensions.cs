using System.Linq;
using MassSpectrometry.Deconvolution.Consensus;

namespace Readers
{
    /// <summary>
    /// Convert a consensus-tracer <see cref="MassFeature"/> into an
    /// <see cref="Ms1Feature"/> record for round-trip through
    /// <see cref="Ms1FeatureFile"/>. The intended use: write the
    /// consensus pipeline's output as a FLASHDeconv-style
    /// <c>_ms1.feature</c> file so downstream consumers (e.g.,
    /// MetaMorpheus's <c>FromFileDeconvolutionParameters</c>) can
    /// pair the features with MS2 scans without any new wire format.
    /// </summary>
    public static class Ms1FeatureFromMassFeatureExtensions
    {
        /// <summary>
        /// Build an <see cref="Ms1Feature"/> record from a
        /// <see cref="MassFeature"/>. Field mapping follows the
        /// FLASHDeconv canonical schema (writer side uses the first
        /// <c>[Name(...)]</c> alias as the column name; reader side
        /// accepts both FLASHDeconv and newer TopFD aliases).
        /// </summary>
        /// <param name="feature">A finalised cross-charge feature
        /// (caller is responsible for having invoked
        /// <see cref="MassFeature.Finalise"/>).</param>
        /// <param name="sequentialId">Feature ID for the produced row.
        /// The consensus pipeline's internal ID is not used because it
        /// may collide across pipelines; downstream callers usually
        /// just want 1..N for a single file.</param>
        /// <param name="sampleId">Sample identifier. Default 0 fits the
        /// single-mzML case; multi-file producers can pass a fraction
        /// index here.</param>
        /// <param name="fractionId">Fraction identifier (both min and
        /// max are set to this value). Default 0 fits a single fraction.</param>
        public static Ms1Feature ToMs1Feature(
            this MassFeature feature,
            int sequentialId,
            int sampleId = 0,
            int fractionId = 0)
        {
            // Apex = max-intensity envelope of the max-intensity constituent
            // trace. "Apex of the dominant charge state at its most intense
            // scan" matches what FLASHDeconv/TopFD typically report.
            var dominantTrace = feature.Traces
                .OrderByDescending(t => t.TotalIntensity)
                .First();
            var apexEnvelope = dominantTrace.Envelopes
                .OrderByDescending(e => e.Intensity)
                .First();

            return new Ms1Feature
            {
                SampleId = sampleId,
                Id = sequentialId,
                Mass = feature.ConsensusMass,
                Intensity = feature.SummedIntensity,
                RetentionTimeBegin = feature.RTStart,
                RetentionTimeEnd = feature.RTEnd,
                RetentionTimeApex = apexEnvelope.RT,
                IntensityApex = apexEnvelope.Intensity,
                ChargeStateMin = feature.Charges.Min(),
                ChargeStateMax = feature.Charges.Max(),
                FractionIdMin = fractionId,
                FractionIdMax = fractionId,
            };
        }
    }
}
