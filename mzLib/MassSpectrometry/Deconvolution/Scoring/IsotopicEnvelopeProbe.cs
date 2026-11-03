using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Gathers the observed peaks that a species of a given monoisotopic mass and charge would
    /// occupy in a spectrum, and returns them as an <see cref="IsotopicEnvelope"/> that can be
    /// handed to <see cref="DeconvolutionScorer"/>.
    ///
    /// <para>
    /// Deconvolution answers "what is in this spectrum". This answers the narrower question
    /// "how well does this spectrum support THIS hypothesis", which deconvolution cannot be asked
    /// directly because it only ever reports the hypotheses it chose. Scoring a mass and charge
    /// that deconvolution did not report is exactly what a decoy needs: a decoy feature has no
    /// peaks of its own, so it cannot be scored by any path that starts from an observed envelope.
    /// </para>
    ///
    /// <para>
    /// The gathered envelope may be empty or nearly so, and that is the informative case. A
    /// hypothesis the spectrum does not support returns few peaks and scores badly, which is the
    /// behaviour a null model depends on.
    /// </para>
    /// </summary>
    public static class IsotopicEnvelopeProbe
    {
        /// <summary>
        /// Collects the peaks in <paramref name="spectrum"/> consistent with a species of
        /// <paramref name="monoisotopicMass"/> at <paramref name="charge"/>.
        /// </summary>
        /// <param name="spectrum">The spectrum to probe.</param>
        /// <param name="monoisotopicMass">Hypothesised monoisotopic mass.</param>
        /// <param name="charge">Hypothesised charge; must be non-zero.</param>
        /// <param name="model">Isotope model supplying the theoretical envelope.</param>
        /// <param name="tolerancePpm">Peak-matching tolerance.</param>
        /// <returns>
        /// An envelope holding every theoretical isotope position that found a peak within
        /// tolerance, or <c>null</c> when the spectrum is empty or not one peak matched. Null
        /// rather than an empty envelope because <see cref="IsotopicEnvelope"/> requires at least
        /// one peak to describe anything, and callers must decide what an unsupported hypothesis
        /// means for them rather than receive a zero that looks like a measurement.
        /// </returns>
        public static IsotopicEnvelope Gather(
            MzSpectrum spectrum,
            double monoisotopicMass,
            int charge,
            AverageResidue model,
            double tolerancePpm)
        {
            if (spectrum == null) throw new ArgumentNullException(nameof(spectrum));
            if (model == null) throw new ArgumentNullException(nameof(model));
            if (charge == 0) throw new ArgumentOutOfRangeException(nameof(charge), "charge must be non-zero");
            if (spectrum.Size == 0) return null;

            int massIndex = model.GetMostIntenseMassIndex(monoisotopicMass);
            double[] theorMasses = model.GetAllTheoreticalMasses(massIndex);
            double[] theorIntensities = model.GetAllTheoreticalIntensities(massIndex);
            if (theorMasses.Length == 0) return null;

            // The model's arrays are indexed from the most intense peak, and describe an averagine
            // of the table's mass rather than of the exact hypothesised mass. Anchoring on the
            // apex and shifting the whole comb keeps the spacing right while placing it at the
            // mass actually being tested.
            double apexMass = monoisotopicMass + model.GetDiffToMonoisotopic(massIndex);
            double offset = apexMass - theorMasses[0];

            var peaks = new List<(double mz, double intensity)>();
            var ratios = new List<double>();
            double totalIntensity = 0;

            for (int i = 0; i < theorMasses.Length; i++)
            {
                double targetMass = theorMasses[i] + offset;
                double targetMz = targetMass.ToMz(charge);

                int idx = spectrum.GetClosestPeakIndex(targetMz);
                double mz = spectrum.XArray[idx];
                double intensity = spectrum.YArray[idx];
                double observedMass = mz.ToMass(charge);

                if (Math.Abs(observedMass - targetMass) / targetMass * 1e6 > tolerancePpm)
                    continue;
                // A single observed peak must not be claimed by two theoretical positions; at high
                // charge adjacent teeth can round to the same peak.
                if (peaks.Any(p => p.mz == mz))
                    continue;

                peaks.Add((mz, intensity));
                totalIntensity += intensity;
                if (intensity > 0)
                    ratios.Add(theorIntensities[i] / intensity);
            }

            if (peaks.Count == 0) return null;

            return new IsotopicEnvelope(
                peaks.Select(p => (p.mz, p.intensity)).ToList(),
                monoisotopicMass,
                charge,
                totalIntensity,
                ratios.Count > 1 ? ratios.StandardDeviation() : 0.0);
        }

        /// <summary>
        /// Probes for the hypothesis and scores it with the spectrum-aware scorer, returning
        /// <see cref="double.NaN"/> when the spectrum supports no peaks at all. NaN, not 0, so a
        /// hypothesis that could not be evaluated stays distinguishable from one that evaluated
        /// badly; a caller building a score distribution must decide which of those it wants.
        /// </summary>
        public static double Score(
            MzSpectrum spectrum,
            double monoisotopicMass,
            int charge,
            AverageResidue model,
            double tolerancePpm)
        {
            var env = Gather(spectrum, monoisotopicMass, charge, model, tolerancePpm);
            if (env == null) return double.NaN;
            return DeconvolutionScorer.ComputeScoreWithSpectrumContext(
                DeconvolutionScorer.ComputeFeatures(env, model, spectrum));
        }
    }
}
