using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry.Deconvolution.Consensus
{
    /// <summary>
    /// A single mass trace: a sequence of per-scan envelopes that the grouper
    /// believes belong to the same species. Same charge across all entries,
    /// scan-adjacency &lt;= MaxGap, mass within tolerance of <see cref="AnchorMass"/>.
    ///
    /// <see cref="MassTraceBuilder"/> populates these in scan-order. Off-by-one
    /// correction is then handled by <see cref="TraceCorrector"/>, which wraps
    /// each trace in a <see cref="CorrectedTrace"/> rather than mutating it.
    /// Mutation policy on the trace itself: contents are appended during
    /// construction and read-only thereafter.
    /// </summary>
    public sealed class MassTrace
    {
        public int Id;
        public int Charge;

        /// <summary>
        /// First envelope's mass; never updated after trace creation. Keeping
        /// the anchor fixed bounds the trace to AnchorMass +/- tolerance and
        /// prevents drift over long traces.
        /// </summary>
        public double AnchorMass;

        public List<(int ScanIndex, int ScanNumber, double RT, double Mass, double Intensity)> Envelopes
            = new();

        public int LastScanIndex => Envelopes[^1].ScanIndex;
        public double MinMass => Envelopes.Min(e => e.Mass);
        public double MaxMass => Envelopes.Max(e => e.Mass);
        public double MassSpread => MaxMass - MinMass;
    }
}
