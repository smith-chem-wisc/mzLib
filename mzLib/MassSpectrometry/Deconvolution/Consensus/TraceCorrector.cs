using System.Linq;

namespace MassSpectrometry.Deconvolution.Consensus
{
    /// <summary>
    /// Off-by-one correction. Computes a per-trace consensus mass (weighted
    /// or unweighted median, caller's choice) and rescues envelopes whose
    /// mass differs from consensus by approximately +/-1.00335 Da.
    ///
    /// Median over weighted mean: a mean can be shifted off the true mass
    /// by a minority of off-by-one outliers, especially if they carry high
    /// weight. A (weighted) median is immune unless &gt; 50% of the trace's
    /// total weight is off-by-one.
    ///
    /// Rescue (not drop): a flagged envelope keeps its place in the trace
    /// with its <see cref="CorrectedEnvelope.CorrectedMass"/> set to
    /// consensus and <see cref="CorrectedEnvelope.WasCorrected"/> = true.
    /// The original mass is preserved for diagnostic inspection.
    ///
    /// Weight schemes the caller can pick from:
    ///   - uniform (every envelope counts equally) -- the default.
    ///   - intensity (TotalIntensity per envelope) -- trusts intense scans.
    ///   - scorer (e.g., DeconvolutionScorer.ScoreEnvelope) -- trusts
    ///     well-shaped envelopes regardless of intensity.
    /// </summary>
    public static class TraceCorrector
    {
        /// <summary>
        /// Average mass spacing between adjacent C12/C13 isotope peaks. An
        /// off-by-one anchor error shifts the reported monoisotope by this.
        /// </summary>
        public const double IsotopeSpacingDa = 1.00335;

        /// <summary>
        /// Half-window for matching |delta - +/-isotope| &lt;= this to declare
        /// an envelope an off-by-one outlier.
        /// </summary>
        public const double DefaultOffByOneWindowDa = 0.05;

        /// <summary>
        /// Type of the per-envelope weight function. Inputs match the named
        /// fields of <see cref="MassTrace.Envelopes"/> tuple entries; output
        /// is a non-negative weight. A function returning 1.0 for every
        /// envelope is the unweighted-median case.
        /// </summary>
        public delegate double EnvelopeWeight(int scanIndex, int scanNumber, double rt, double mass, int charge, double intensity);

        public static readonly EnvelopeWeight UniformWeight = (_, _, _, _, _, _) => 1.0;
        public static readonly EnvelopeWeight IntensityWeight = (_, _, _, _, _, intensity) => intensity;

        public static CorrectedTrace Correct(MassTrace t)
            => Correct(t, UniformWeight, DefaultOffByOneWindowDa);

        public static CorrectedTrace Correct(MassTrace t, EnvelopeWeight weight, double offByOneWindowDa = DefaultOffByOneWindowDa)
        {
            // Weighted median: sort envelopes by mass, take the value at which
            // cumulative weight first reaches half the total. Falls back to
            // ordinary median when all weights are equal.
            var items = t.Envelopes
                .Select(e => (Mass: e.Mass, Weight: System.Math.Max(0, weight(e.ScanIndex, e.ScanNumber, e.RT, e.Mass, t.Charge, e.Intensity))))
                .OrderBy(x => x.Mass)
                .ToArray();
            double totalWeight = items.Sum(x => x.Weight);
            double consensus;
            if (totalWeight <= 0)
            {
                // Degenerate (all weights zero) -- fall back to unweighted median.
                var masses = t.Envelopes.Select(e => e.Mass).OrderBy(m => m).ToArray();
                consensus = masses[masses.Length / 2];
            }
            else
            {
                double half = totalWeight / 2.0;
                double cum = 0;
                consensus = items[^1].Mass;
                foreach (var (mass, w) in items)
                {
                    cum += w;
                    if (cum >= half) { consensus = mass; break; }
                }
            }

            var result = new CorrectedTrace
            {
                Id = t.Id,
                Charge = t.Charge,
                ConsensusMass = consensus,
            };

            double originalMin = double.PositiveInfinity, originalMax = double.NegativeInfinity;
            double correctedMin = double.PositiveInfinity, correctedMax = double.NegativeInfinity;

            foreach (var e in t.Envelopes)
            {
                double delta = e.Mass - consensus;
                bool isOffByOne =
                    System.Math.Abs(delta - IsotopeSpacingDa) <= offByOneWindowDa ||
                    System.Math.Abs(delta + IsotopeSpacingDa) <= offByOneWindowDa;

                double correctedMass = isOffByOne ? consensus : e.Mass;

                result.Envelopes.Add(new CorrectedEnvelope
                {
                    ScanIndex = e.ScanIndex,
                    ScanNumber = e.ScanNumber,
                    RT = e.RT,
                    OriginalMass = e.Mass,
                    CorrectedMass = correctedMass,
                    Charge = t.Charge,
                    Intensity = e.Intensity,
                    WasCorrected = isOffByOne,
                });

                if (e.Mass < originalMin) originalMin = e.Mass;
                if (e.Mass > originalMax) originalMax = e.Mass;
                if (correctedMass < correctedMin) correctedMin = correctedMass;
                if (correctedMass > correctedMax) correctedMax = correctedMass;
            }

            result.OriginalSpread = originalMax - originalMin;
            result.CorrectedSpread = correctedMax - correctedMin;
            return result;
        }
    }
}
