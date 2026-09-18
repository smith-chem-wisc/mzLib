using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry.Deconvolution.Consensus
{
    /// <summary>
    /// Assigns a false-discovery-rate q-value to each <see cref="MassFeature"/> by competing the
    /// real features against decoy features built from them.
    ///
    /// <para><b>Why features need their own null</b></para>
    /// <para>
    /// mzLib already estimates FDR for individual isotopic envelopes during deconvolution
    /// (<see cref="DeconvolutionQValueCalculator"/> against <see cref="DecoyAveragine"/>). That
    /// operates one scan at a time and says nothing about an assembled feature, which is what a
    /// consumer of an <c>_ms1.feature</c> file actually receives. A feature list with no error
    /// rate can only be used whole or thresholded by guesswork.
    /// </para>
    ///
    /// <para><b>How the decoys are built</b></para>
    /// <para>
    /// Each decoy keeps a real feature's charge and retention time and displaces its mass by a
    /// half-integer multiple of the <sup>13</sup>C spacing, 5 to 50 Da away. Keeping charge and
    /// retention time means the decoy probes the same spectrum at the same density as its target,
    /// so the comparison isolates the one thing that makes a feature real: whether a complete
    /// isotope envelope actually sits at that mass. The half-integer offset places the decoy comb
    /// between the real teeth rather than on them.
    /// </para>
    /// <para>
    /// Displacing the mass, rather than perturbing the isotope spacing as the envelope-level decoy
    /// does, matters at top-down charge states. A spacing perturbation is a fixed offset in neutral
    /// mass, so the matcher, which works in m/z, sees it divided by charge; past roughly z = 4 at a
    /// 20 ppm tolerance the perturbed comb lands back on the real peaks and the decoy stops being
    /// wrong. A whole-feature displacement of 5 to 50 Da does not shrink with charge.
    /// </para>
    ///
    /// <para><b>What the q-value means, and does not</b></para>
    /// <para>
    /// It estimates the fraction of surviving features that are not real isotope envelopes. It does
    /// NOT estimate the fraction that fail to correspond to an identifiable molecule: most true
    /// features are never identified by any search, and a feature can be a perfectly real envelope
    /// of a species no database contains. Read it as "is this a real envelope", not "is this a real
    /// proteoform".
    /// </para>
    /// </summary>
    public static class ConsensusFeatureFdr
    {
        /// <summary>Default peak-matching tolerance when probing a hypothesis.</summary>
        public const double DefaultTolerancePpm = 20.0;

        /// <summary>Smallest and largest decoy mass displacement, in daltons.</summary>
        public const double MinDecoyOffsetDa = 5.0;
        public const double MaxDecoyOffsetDa = 50.0;

        /// <summary>
        /// Scores every feature and its decoy at the feature's apex scan, then assigns
        /// <see cref="MassFeature.QValue"/>.
        /// </summary>
        /// <param name="features">Finalised features. Mutated in place.</param>
        /// <param name="ms1Scans">The MS1 scans the features were traced from, ascending in retention time.</param>
        /// <param name="model">Isotope model; use the same one the features were deconvolved with.</param>
        /// <param name="tolerancePpm">Peak-matching tolerance for the probe.</param>
        /// <param name="seed">Fixed so a rerun reproduces the same decoys. Vary it to confirm a result does not depend on one draw.</param>
        /// <returns>The number of features that received a finite q-value.</returns>
        /// <remarks>
        /// Target and decoy are both scored by probing the apex scan, rather than reusing
        /// <see cref="MassFeature.QualityScore"/>, because a decoy has no trace to aggregate over.
        /// Comparing a trace-averaged target score against a single-scan decoy score would make the
        /// target look better for a reason that has nothing to do with being real. The trace
        /// aggregate remains the descriptive per-feature number; this is the matched one.
        /// </remarks>
        public static int AssignQValues(
            IReadOnlyList<MassFeature> features,
            IReadOnlyList<MsDataScan> ms1Scans,
            AverageResidue model,
            double tolerancePpm = DefaultTolerancePpm,
            int seed = 42)
        {
            if (features == null) throw new ArgumentNullException(nameof(features));
            if (ms1Scans == null) throw new ArgumentNullException(nameof(ms1Scans));
            if (model == null) throw new ArgumentNullException(nameof(model));
            if (features.Count == 0) return 0;
            if (ms1Scans.Count == 0) throw new ArgumentException("no MS1 scans supplied", nameof(ms1Scans));

            var rng = new Random(seed);
            var scored = new List<(MassFeature Feature, double Target)>();
            var decoyScores = new List<double>();

            foreach (var f in features)
            {
                MsDataScan apex = NearestScan(ms1Scans, ApexRetentionTime(f));
                int charge = RepresentativeCharge(f);

                double target = IsotopicEnvelopeProbe.Score(
                    apex.MassSpectrum, f.ConsensusMass, charge, model, tolerancePpm);

                double decoyMass = f.ConsensusMass + DecoyOffset(rng);
                double decoy = IsotopicEnvelopeProbe.Score(
                    apex.MassSpectrum, decoyMass, charge, model, tolerancePpm);

                // A hypothesis the spectrum cannot support at all scores NaN. For the target that
                // means no q-value can be assigned. For the decoy it means the null found nothing,
                // which is a legitimate outcome and enters the null distribution as the lowest
                // possible score rather than being dropped -- dropping it would bias the null
                // upward and make the FDR look better than it is.
                if (!double.IsNaN(target))
                    scored.Add((f, target));
                decoyScores.Add(double.IsNaN(decoy) ? 0.0 : decoy);

                f.QValue = double.NaN;
            }

            if (scored.Count == 0) return 0;

            double[] q = DeconvolutionQValueCalculator.AssignQValues(
                scored.Select(s => s.Target).ToList(), decoyScores);

            for (int i = 0; i < scored.Count; i++)
                scored[i].Feature.QValue = q[i];

            return scored.Count;
        }

        /// <summary>Half-integer multiple of the 13C spacing, 5 to 50 Da, either side.</summary>
        private static double DecoyOffset(Random rng)
        {
            double magnitude = MinDecoyOffsetDa + rng.NextDouble() * (MaxDecoyOffsetDa - MinDecoyOffsetDa);
            double sign = rng.Next(2) == 0 ? -1.0 : 1.0;
            double k = Math.Round(magnitude / Constants.C13MinusC12);
            return sign * (k + 0.5) * Constants.C13MinusC12;
        }

        /// <summary>Retention time of the most intense envelope in the feature.</summary>
        private static double ApexRetentionTime(MassFeature f)
        {
            var envelopes = f.Traces.SelectMany(t => t.Envelopes).ToList();
            if (envelopes.Count == 0) return (f.RTStart + f.RTEnd) / 2;
            double best = envelopes.Max(e => e.Intensity);
            return envelopes.First(e => e.Intensity == best).RT;
        }

        /// <summary>
        /// Charge of the feature's most intense trace. The abundant charge is the one most likely
        /// to have a complete envelope at the apex; the midpoint of the charge range can fall on a
        /// charge state carrying almost no signal, which would make a real feature look unsupported.
        /// </summary>
        private static int RepresentativeCharge(MassFeature f)
        {
            if (f.Traces.Count == 0) return f.Charges.Count > 0 ? f.Charges.Min() : 1;
            return f.Traces.OrderByDescending(t => t.TotalIntensity).First().Charge;
        }

        private static MsDataScan NearestScan(IReadOnlyList<MsDataScan> scans, double rt)
        {
            MsDataScan best = scans[0];
            double bestDelta = Math.Abs(scans[0].RetentionTime - rt);
            for (int i = 1; i < scans.Count; i++)
            {
                double d = Math.Abs(scans[i].RetentionTime - rt);
                if (d < bestDelta) { bestDelta = d; best = scans[i]; }
            }
            return best;
        }
    }
}
