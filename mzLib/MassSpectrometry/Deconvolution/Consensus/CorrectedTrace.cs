using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry.Deconvolution.Consensus
{
    /// <summary>
    /// A <see cref="MassTrace"/> after off-by-one correction by
    /// <see cref="TraceCorrector"/>. <see cref="ConsensusMass"/> is the
    /// per-trace (weighted) median; envelopes whose original mass differed
    /// by approximately +/-1.00335 Da from it have been snapped to the
    /// consensus and have <see cref="CorrectedEnvelope.WasCorrected"/> set.
    /// Original masses are preserved on each envelope for diagnostic
    /// inspection.
    /// </summary>
    public sealed class CorrectedTrace
    {
        public int Id;
        public int Charge;
        public double ConsensusMass;
        public List<CorrectedEnvelope> Envelopes = new();
        public double OriginalSpread;
        public double CorrectedSpread;
        public int CorrectionCount => Envelopes.Count(e => e.WasCorrected);

        public double FirstRT => Envelopes[0].RT;
        public double LastRT => Envelopes[^1].RT;
        public double TotalIntensity => Envelopes.Sum(e => e.Intensity);
    }
}
