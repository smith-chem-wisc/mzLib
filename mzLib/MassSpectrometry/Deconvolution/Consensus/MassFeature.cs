using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry.Deconvolution.Consensus
{
    /// <summary>
    /// CONSENSUS PIPELINE OVERVIEW (the type hierarchy, narrowest to widest).
    ///
    /// The consensus mass-tracing pipeline turns raw per-scan deconvolution
    /// output into cross-charge features through four nested levels. Each level
    /// groups the one below it:
    ///
    ///   IsotopicEnvelope      one deconvolved species (mass, intensity, charge)
    ///     |                   reported in a SINGLE MS1 scan. Raw input; not a
    ///     |                   type in this namespace.
    ///     |   grouped by <see cref="MassTraceBuilder"/> (charge-locked,
    ///     |   anchor-mass + scan-adjacency)
    ///     v
    ///   <see cref="MassTrace"/>        the same species followed across adjacent
    ///     |                   scans at ONE charge state. Its envelopes are the
    ///     |                   raw per-scan tuples above.
    ///     |   corrected by <see cref="TraceCorrector"/> (off-by-one rescue +
    ///     |   resolution-aware splitting)
    ///     v
    ///   <see cref="CorrectedTrace"/>   a MassTrace after correction, holding
    ///     |                   <see cref="CorrectedEnvelope"/> entries (original
    ///     |                   + corrected mass per envelope) and a per-trace
    ///     |                   <see cref="CorrectedTrace.ConsensusMass"/>. Still
    ///     |                   one charge state.
    ///     |   stitched by <see cref="MassFeatureBuilder"/> (cross-charge,
    ///     |   ppm mass agreement + RT overlap)
    ///     v
    ///   MassFeature           THIS type: one species across ALL of its charge
    ///                         states. The widest grouping, and the unit the
    ///                         writer turns into an Ms1Feature row.
    ///
    /// So: envelopes nest inside traces, traces (after correction) nest inside
    /// features. Charge is fixed within a trace and varies across a feature.
    ///
    /// A cross-charge-state consensus feature: a group of
    /// <see cref="CorrectedTrace"/> entries (each at one charge state) whose
    /// consensus masses agree within a ppm tolerance and whose RT windows
    /// overlap.
    ///
    /// A real proteoform or peptide produces envelopes at several charge
    /// states (BU peptides commonly +2/+3; TD proteoforms across +10..+15
    /// or wider). The trace builder is charge-locked, so each charge gets
    /// its own trace. <see cref="MassFeatureBuilder"/> stitches those
    /// per-charge traces back together. Charge multiplicity then becomes
    /// a confidence signal: a feature seen at multiple charges is
    /// corroborated by independent charge calculations; a single-charge,
    /// single-envelope feature is more likely to be noise.
    ///
    /// Mutation policy: callers append to <see cref="Traces"/> during
    /// construction, then call <see cref="Finalise"/> exactly once to
    /// derive the aggregate fields. After Finalise, the feature is
    /// treated as read-only by downstream consumers.
    /// </summary>
    public sealed class MassFeature
    {
        public int Id;
        public List<CorrectedTrace> Traces = new();
        public double ConsensusMass;
        public HashSet<int> Charges = new();
        public int ChargeCount => Charges.Count;
        public int MaxTraceLength;
        public double RTStart;
        public double RTEnd;
        public double SummedIntensity;

        /// <summary>
        /// Populate derived fields from the current <see cref="Traces"/>
        /// list. Idempotent; safe to re-run after the trace list changes.
        /// <see cref="ConsensusMass"/> is the intensity-weighted mean of
        /// per-trace consensus masses (heavier traces contribute more).
        /// </summary>
        public void Finalise()
        {
            if (Traces.Count == 0)
                throw new System.InvalidOperationException("Cannot Finalise a MassFeature with no traces.");

            Charges = new HashSet<int>(Traces.Select(t => t.Charge));
            SummedIntensity = Traces.Sum(t => t.TotalIntensity);
            RTStart = Traces.Min(t => t.FirstRT);
            RTEnd = Traces.Max(t => t.LastRT);
            MaxTraceLength = Traces.Max(t => t.Envelopes.Count);

            // Intensity-weighted mean of per-trace consensus masses. Heavier
            // traces (more intense, more confident) contribute more.
            double w = Traces.Sum(t => t.TotalIntensity);
            ConsensusMass = w == 0
                ? Traces.Average(t => t.ConsensusMass)
                : Traces.Sum(t => t.ConsensusMass * t.TotalIntensity) / w;
        }
    }
}
