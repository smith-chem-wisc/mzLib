using System.Collections.Generic;
using System.Linq;
using Chemistry;

namespace MassSpectrometry.Deconvolution.Consensus
{
    /// <summary>
    /// Off-by-one correction with resolution-aware splitting.
    ///
    /// Computes a per-trace consensus mass (weighted or unweighted median, caller's
    /// choice). Envelopes whose mass differs from consensus by approximately an integer
    /// multiple of the C12/C13 spacing (1.00335 Da) are off-by-one anchor errors and are
    /// rescued (snapped to consensus). Envelopes that sit at a resolvable, NON-isotope
    /// offset (e.g. a deamidated form at +0.98402 Da) are a distinct co-grouped species:
    /// they are SPLIT OUT into their own <see cref="CorrectedTrace"/> rather than erased.
    ///
    /// The discriminating window is derived from the data, not hardwired: an off-by-one
    /// match (and a cluster boundary) is declared within k*sigma of the integer isotope
    /// offset, where sigma is the per-trace mass scatter estimated by the median absolute
    /// deviation of the envelopes about consensus. On high-resolution data sigma is small,
    /// so the 0.0193 Da gap between the isotope spacing (1.00335) and deamidation (0.98402)
    /// is resolved and the deamidated form surfaces as its own feature. On low-resolution
    /// data sigma is large, the two cannot be told apart, and the conservative fallback is
    /// to treat the offset as an isotope error and merge it. The only constants are
    /// <see cref="DefaultSigmaMultiple"/>, the cap <see cref="MaxWindowDa"/>, and the
    /// physics constant <see cref="IsotopeSpacingDa"/>.
    ///
    /// Median over weighted mean: a mean can be shifted off the true mass by a minority of
    /// off-by-one outliers, especially if they carry high weight. A (weighted) median is
    /// immune unless &gt; 50% of the trace's total weight is off-by-one.
    ///
    /// Rescue (not drop): a flagged envelope keeps its place in the trace with its
    /// <see cref="CorrectedEnvelope.CorrectedMass"/> set to consensus and
    /// <see cref="CorrectedEnvelope.WasCorrected"/> = true; the original mass is preserved.
    ///
    /// Weight schemes the caller can pick from:
    ///   - uniform (every envelope counts equally) -- the default.
    ///   - intensity (TotalIntensity per envelope) -- trusts intense scans.
    ///   - scorer (e.g., DeconvolutionScorer.ScoreEnvelope) -- trusts well-shaped envelopes.
    /// </summary>
    public static class TraceCorrector
    {
        /// <summary>
        /// Average mass spacing between adjacent C12/C13 isotope peaks. An off-by-one
        /// anchor error shifts the reported monoisotope by an integer multiple of this.
        /// Aliases the shared <see cref="Constants.C13MinusC12"/> rather than redefining it.
        /// </summary>
        public const double IsotopeSpacingDa = Constants.C13MinusC12;

        /// <summary>
        /// Half-window, in units of estimated mass-scatter sigma, for matching an envelope
        /// to an integer isotope offset (off-by-one) and for delimiting distinct mass
        /// clusters. 3 sigma is the conventional choice.
        /// </summary>
        public const double DefaultSigmaMultiple = 3.0;

        /// <summary>
        /// Numerical floor on the per-trace sigma estimate (Da). Stops a zero-scatter
        /// cluster (identical or synthetic masses) from yielding a zero-width window;
        /// real measurement scatter normally dominates this.
        /// </summary>
        public const double MinSigmaDa = 1e-4;

        /// <summary>
        /// Upper bound on the k*sigma window (Da). Held below half the isotope spacing so
        /// the core (n=0) and adjacent off-by-one (n=+/-1) bands stay disjoint. Only engages
        /// on pathologically scattered traces, where it stops the window from swallowing an
        /// implausibly wide mass range.
        /// </summary>
        public const double MaxWindowDa = IsotopeSpacingDa * 0.45;

        /// <summary>
        /// Type of the per-envelope weight function. Inputs match the named fields of
        /// <see cref="MassTrace.Envelopes"/> tuple entries; output is a non-negative weight.
        /// A function returning 1.0 for every envelope is the unweighted-median case.
        /// </summary>
        public delegate double EnvelopeWeight(int scanIndex, int scanNumber, double rt, double mass, int charge, double intensity);

        public static readonly EnvelopeWeight UniformWeight = (_, _, _, _, _, _) => 1.0;
        public static readonly EnvelopeWeight IntensityWeight = (_, _, _, _, _, intensity) => intensity;

        public static List<CorrectedTrace> Correct(MassTrace t)
            => Correct(t, UniformWeight, DefaultSigmaMultiple);

        /// <summary>
        /// Correct a trace, returning one <see cref="CorrectedTrace"/> for the base species
        /// plus one for each resolvable distinct co-grouped species split out of it.
        /// </summary>
        public static List<CorrectedTrace> Correct(MassTrace t, EnvelopeWeight weight, double sigmaMultiple = DefaultSigmaMultiple)
        {
            if (t.Envelopes.Count == 0)
                throw new System.ArgumentException("Cannot correct a MassTrace with no envelopes.", nameof(t));

            var results = new List<CorrectedTrace>();
            CorrectGroup(t.Id, t.Charge, t.Envelopes.ToList(), weight, sigmaMultiple, results);
            return results;
        }

        /// <summary>
        /// Median-anchor one species: keep the consensus cluster and its integer-isotope
        /// shadows, correct the shadows, and recurse on whatever envelopes sit at a
        /// resolvable non-isotope offset (a distinct species).
        /// </summary>
        private static void CorrectGroup(
            int id, int charge,
            List<(int ScanIndex, int ScanNumber, double RT, double Mass, double Intensity)> envs,
            EnvelopeWeight weight, double sigmaMultiple, List<CorrectedTrace> results)
        {
            double consensus = WeightedMedian(envs, weight, charge);
            double window = System.Math.Min(
                System.Math.Max(MinSigmaDa, EstimateSigma(envs, consensus)) * sigmaMultiple,
                MaxWindowDa);

            var belongs = new List<(int ScanIndex, int ScanNumber, double RT, double Mass, double Intensity)>();
            var other = new List<(int ScanIndex, int ScanNumber, double RT, double Mass, double Intensity)>();
            foreach (var e in envs)
            {
                double delta = e.Mass - consensus;
                bool isCore = System.Math.Abs(delta) <= window;
                long n = (long)System.Math.Round(delta / IsotopeSpacingDa);
                bool isShadow = n != 0 && System.Math.Abs(delta - n * IsotopeSpacingDa) <= window;
                (isCore || isShadow ? belongs : other).Add(e);
            }

            var result = new CorrectedTrace { Id = id, Charge = charge, ConsensusMass = consensus };
            double oMin = double.PositiveInfinity, oMax = double.NegativeInfinity;
            double cMin = double.PositiveInfinity, cMax = double.NegativeInfinity;
            foreach (var e in belongs)
            {
                bool wasCorrected = System.Math.Abs(e.Mass - consensus) > window; // a shadow, not core
                double correctedMass = wasCorrected ? consensus : e.Mass;

                result.Envelopes.Add(new CorrectedEnvelope
                {
                    ScanIndex = e.ScanIndex,
                    ScanNumber = e.ScanNumber,
                    RT = e.RT,
                    OriginalMass = e.Mass,
                    CorrectedMass = correctedMass,
                    Charge = charge,
                    Intensity = e.Intensity,
                    WasCorrected = wasCorrected,
                });

                if (e.Mass < oMin) oMin = e.Mass;
                if (e.Mass > oMax) oMax = e.Mass;
                if (correctedMass < cMin) cMin = correctedMass;
                if (correctedMass > cMax) cMax = correctedMass;
            }
            result.OriginalSpread = oMax - oMin;
            result.CorrectedSpread = cMax - cMin;
            results.Add(result);

            // No minimum-support gate before splitting: a distinct cluster becomes its own
            // trace even if it has only a few envelopes. A count threshold here cannot tell a
            // faint real species (e.g. a low-abundance deamidated form, often just a scan or
            // two) from a noise envelope -- both are few-envelope -- so it would discard low
            // signal, not just noise. Quality-based culling belongs downstream where the
            // information lives: per-envelope deconvolution scoring / ML rescoring, and
            // cross-charge corroboration at the MassFeature level (ChargeCount / multiplicity).
            //
            // The split-off species is itself a trace: re-anchor and correct it (it may carry
            // its own off-by-one shadows). Each pass strips >=1 envelope (consensus is always
            // an actual envelope mass, hence always core), so the recursion terminates.
            if (other.Count > 0)
                CorrectGroup(id, charge, other, weight, sigmaMultiple, results);
        }

        private static double WeightedMedian(
            List<(int ScanIndex, int ScanNumber, double RT, double Mass, double Intensity)> envs,
            EnvelopeWeight weight, int charge)
        {
            var items = envs
                .Select(e => (e.Mass, Weight: System.Math.Max(0, weight(e.ScanIndex, e.ScanNumber, e.RT, e.Mass, charge, e.Intensity))))
                .OrderBy(x => x.Mass)
                .ToArray();
            double totalWeight = items.Sum(x => x.Weight);
            if (totalWeight <= 0)
            {
                // Degenerate (all weights zero) -- fall back to unweighted median.
                var masses = envs.Select(e => e.Mass).OrderBy(m => m).ToArray();
                return masses[masses.Length / 2];
            }
            double half = totalWeight / 2.0, cum = 0, consensus = items[^1].Mass;
            foreach (var (mass, w) in items)
            {
                cum += w;
                if (cum >= half) { consensus = mass; break; }
            }
            return consensus;
        }

        /// <summary>
        /// Robust per-trace mass scatter via the median absolute deviation about consensus
        /// (scaled to a Gaussian sigma). MAD ignores the minority off-isotope/distinct
        /// clusters, so it reads the precision of the dominant species.
        ///
        /// Limitation: this holds only while the dominant cluster is a clear majority. Near
        /// a 50/50 split the median deviation can land in the gap and inflate sigma, which
        /// suppresses the split. Acceptable for the minority-PTM case this targets; a
        /// tightest-cluster estimator would be the hardening if a balanced-split regime ever
        /// matters.
        /// </summary>
        private static double EstimateSigma(
            List<(int ScanIndex, int ScanNumber, double RT, double Mass, double Intensity)> envs,
            double consensus)
        {
            var devs = envs.Select(e => System.Math.Abs(e.Mass - consensus)).OrderBy(d => d).ToArray();
            double mad = devs[devs.Length / 2];
            return mad / 0.6744897501960817; // MAD -> sigma for a normal distribution
        }
    }
}
