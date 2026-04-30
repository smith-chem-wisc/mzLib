using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;

namespace Test.MassSpectrometry.Deconvolution
{
    /// <summary>
    /// Shared test helpers for deconvolution scorer tests.
    /// </summary>
    internal static class DeconvolutionTestHelpers
    {
        internal static readonly AverageResidue Model = new Averagine();
        internal const double TestMass = 5000.0;
        internal const int TestCharge = 5;

        /// <summary>
        /// Builds a perfect synthetic envelope: peaks placed at exact theoretical
        /// m/z positions with Averagine-proportional intensities. No noise, no
        /// missing peaks. The peak list is suitable for verifying ideal feature values.
        /// </summary>
        internal static IsotopicEnvelope BuildPerfectEnvelope(
            double monoMass = TestMass,
            int charge = TestCharge,
            double baseIntens = 1e6,
            double score = 0.999)
        {
            int avgIdx = Model.GetMostIntenseMassIndex(monoMass);

            double[] rawMasses = Model.GetAllTheoreticalMasses(avgIdx);
            double[] rawIntens = Model.GetAllTheoreticalIntensities(avgIdx);

            // Mass-ascending sort so index 0 = monoisotopic
            var sorted = rawMasses.Zip(rawIntens)
                .OrderBy(p => p.First)
                .ToArray();

            int absCharge = Math.Abs(charge);
            double isotopeStep = Constants.C13MinusC12 / absCharge;
            double monoMz = monoMass.ToMz(charge);
            var peaks = new List<(double mz, double intensity)>();

            for (int n = 0; n < sorted.Length; n++)
            {
                double intensity = baseIntens * sorted[n].Second;
                if (intensity < baseIntens * 0.001) continue; // skip near-zero peaks
                double mz = monoMz + n * isotopeStep;
                peaks.Add((mz, intensity));
            }

            return new IsotopicEnvelope(
                id: 0,
                peaks: peaks,
                monoisotopicmass: monoMass,
                chargestate: charge,
                intensity: peaks.Sum(p => p.intensity),
                score: score);
        }
    }
}
