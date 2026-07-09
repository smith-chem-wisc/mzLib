using System.Collections.Generic;
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
            EnsureFinalised(feature);
            var (apexRt, apexIntensity) = Apex(feature);
            // Single row spanning the full observed charge span. For a gapped charge
            // set this collapses the gaps; callers that need exact fidelity should use
            // ToMs1Features (which FromMassFeatures does).
            return BuildRow(feature, sequentialId, sampleId, fractionId,
                feature.Charges.Min(), feature.Charges.Max(), apexRt, apexIntensity);
        }

        /// <summary>
        /// Emit one <see cref="Ms1Feature"/> row per <em>contiguous</em> run of charge
        /// states. A gapped charge set (e.g. {10, 12, 15}, which arises when an
        /// intermediate charge falls below the deconvolution score cutoff) would
        /// otherwise be written as a single Min=10/Max=15 row that the reader re-expands
        /// to 10..15 inclusive -- fabricating charges 11/13/14. Splitting into contiguous
        /// runs lets the Min/Max-only <c>_ms1.feature</c> schema carry the exact set with
        /// no invented charges. A contiguous set yields exactly one row, identical to
        /// <see cref="ToMs1Feature"/>. Rows are returned without an Id assigned; the
        /// caller numbers them.
        /// </summary>
        public static IReadOnlyList<Ms1Feature> ToMs1Features(
            this MassFeature feature,
            int sampleId = 0,
            int fractionId = 0)
        {
            EnsureFinalised(feature);
            var (apexRt, apexIntensity) = Apex(feature);
            var rows = new List<Ms1Feature>();
            foreach (var (lo, hi) in ContiguousChargeRuns(feature.Charges))
                rows.Add(BuildRow(feature, id: 0, sampleId, fractionId, lo, hi, apexRt, apexIntensity));
            return rows;
        }

        private static void EnsureFinalised(MassFeature feature)
        {
            if (feature.Traces.Count == 0)
                throw new System.ArgumentException("MassFeature has no traces.", nameof(feature));
            if (feature.Charges.Count == 0)
                throw new System.ArgumentException(
                    "MassFeature.Finalise() must be called before ToMs1Feature(...).", nameof(feature));
        }

        // Apex = max-intensity envelope of the max-intensity constituent trace. "Apex of
        // the dominant charge state at its most intense scan" matches what FLASHDeconv/
        // TopFD typically report.
        private static (double apexRt, double apexIntensity) Apex(MassFeature feature)
        {
            var dominantTrace = feature.Traces
                .OrderByDescending(t => t.TotalIntensity)
                .First();
            var apexEnvelope = dominantTrace.Envelopes
                .OrderByDescending(e => e.Intensity)
                .First();
            return (apexEnvelope.RT, apexEnvelope.Intensity);
        }

        /// <summary>Distinct charges sorted ascending, grouped into maximal contiguous (lo, hi) runs.</summary>
        private static IEnumerable<(int lo, int hi)> ContiguousChargeRuns(IEnumerable<int> charges)
        {
            var sorted = charges.Distinct().OrderBy(z => z).ToList();
            int runStart = sorted[0], prev = sorted[0];
            for (int i = 1; i < sorted.Count; i++)
            {
                if (sorted[i] == prev + 1) { prev = sorted[i]; continue; }
                yield return (runStart, prev);
                runStart = prev = sorted[i];
            }
            yield return (runStart, prev);
        }

        private static Ms1Feature BuildRow(MassFeature feature, int id, int sampleId, int fractionId,
            int chargeMin, int chargeMax, double apexRt, double apexIntensity) =>
            new Ms1Feature
            {
                SampleId = sampleId,
                Id = id,
                Mass = feature.ConsensusMass,
                Intensity = feature.SummedIntensity,
                RetentionTimeBegin = feature.RTStart,
                RetentionTimeEnd = feature.RTEnd,
                RetentionTimeApex = apexRt,
                IntensityApex = apexIntensity,
                ChargeStateMin = chargeMin,
                ChargeStateMax = chargeMax,
                FractionIdMin = fractionId,
                FractionIdMax = fractionId,
            };
    }
}
