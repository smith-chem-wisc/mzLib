using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Unit tests for the resolved "most abundant mass" precursor-selection support (Strategy B).
    ///
    /// Terminology pinned by these tests:
    ///  • most-abundant mass = the neutral mass of the single most intense (tallest) isotopic peak,
    ///    proton-corrected (<see cref="IsotopicEnvelope.MostAbundantObservedNeutralMass"/>);
    ///  • most-abundant offset = GetDiffToMonoisotopic(GetMostIntenseMassIndex(mono)) on the averagine
    ///    model — the gap from the monoisotopic mass to that tallest isotopologue.
    /// The intensity-weighted average (centroid) mass and the isotopically-unresolved path are a
    /// separate change and are tested with that work, not here. All tests build synthetic envelopes
    /// from the Averagine model — no deconvolution is run.
    /// </summary>
    [TestFixture]
    public sealed class TestMostAbundantMass
    {
        private static readonly AverageResidue Model = new Averagine();

        // Most-abundant offset for a monoisotopic mass = the averagine diff-to-monoisotopic of the
        // nearest mass bin (the mass-keyed composition consumers use in place of a dedicated method).
        private static double MostAbundantOffset(double monoMass) => Model.GetDiffToMonoisotopic(Model.GetMostIntenseMassIndex(monoMass));

        /// <summary>
        /// Builds a perfect synthetic envelope: peaks at exact theoretical m/z with
        /// Averagine-proportional intensities. (Same construction as TestDeconvolutionScorerUnit.)
        /// </summary>
        private static List<(double mz, double intensity)> BuildPerfectPeaks(double monoMass, int charge, double baseIntens = 1e6)
        {
            int avgIdx = Model.GetMostIntenseMassIndex(monoMass);
            double[] rawMasses = Model.GetAllTheoreticalMasses(avgIdx);
            double[] rawIntens = Model.GetAllTheoreticalIntensities(avgIdx);

            var sorted = rawMasses.Zip(rawIntens).OrderBy(p => p.First).ToArray();

            double isotopeStep = Constants.C13MinusC12 / charge;
            double monoMz = monoMass.ToMz(charge);
            var peaks = new List<(double mz, double intensity)>();
            for (int n = 0; n < sorted.Length; n++)
            {
                double intensity = baseIntens * sorted[n].Second;
                if (intensity < baseIntens * 0.001) continue;
                peaks.Add((monoMz + n * isotopeStep, intensity));
            }
            return peaks;
        }

        private static IsotopicEnvelope BuildPerfectEnvelope(double monoMass, int charge, double baseIntens = 1e6)
        {
            var peaks = BuildPerfectPeaks(monoMass, charge, baseIntens);
            return new IsotopicEnvelope(0, peaks, monoMass, charge, peaks.Sum(p => p.intensity), 0.999);
        }

        // ── AverageResidue most-abundant offset ─────────────────────────────────

        [Test]
        public void MostAbundantOffset_IsNonNegativeAndNonDecreasingWithMass()
        {
            double[] masses = { 500, 2000, 5000, 10000, 20000, 40000 };
            double prev = double.NegativeInfinity;
            foreach (double m in masses)
            {
                double offset = MostAbundantOffset(m);
                Assert.That(offset, Is.GreaterThanOrEqualTo(-1e-6)); // ~0 at tiny mass (mono IS most abundant)
                Assert.That(offset, Is.GreaterThanOrEqualTo(prev - 1e-6), $"offset decreased at mass {m}");
                prev = offset;
            }
        }

        [Test]
        public void MostAbundantOffset_IsNearZeroForSmallMass()
        {
            // A small peptide's monoisotopic peak is (nearly) the most abundant.
            Assert.That(MostAbundantOffset(500), Is.LessThan(0.5));
        }

        [Test]
        public void MostAbundantOffset_GrowsRoughlyOneNeutronPer1600Da()
        {
            // ~1 13C neutron (~1.00235 Da) per ~1.6 kDa. Assert the offset at 16 kDa is in a
            // physically plausible band (≈ 9–11 Da) rather than an exact value.
            double offset = MostAbundantOffset(16000);
            Assert.That(offset, Is.GreaterThan(8.0).And.LessThan(12.0));
        }

        // ── IsotopicEnvelope most-abundant observed mass ────────────────────────

        [Test]
        public void MostAbundantObservedNeutralMass_IsProtonCorrectedNeutralMass()
        {
            const int charge = 10;
            var env = BuildPerfectEnvelope(15000, charge);

            // The proton-corrected neutral mass equals the most intense peak's m/z .ToMass(charge)...
            double mostIntenseMz = env.Peaks.MaxBy(p => p.intensity).mz;
            Assert.That(env.MostAbundantObservedNeutralMass, Is.EqualTo(mostIntenseMz.ToMass(charge)).Within(1e-6));

            // ...and it differs from the un-proton-corrected mz*|z| field by exactly z proton masses.
            Assert.That(env.MostAbundantObservedIsotopicMass - env.MostAbundantObservedNeutralMass,
                Is.EqualTo(charge * Constants.ProtonMass).Within(1e-6));
        }

        [Test]
        public void MostAbundantObservedNeutralMass_MatchesMonoPlusAveragineOffset()
        {
            // Ties the two features together: the observed most-abundant neutral mass of a perfect
            // envelope ≈ candidate monoisotopic + averagine most-abundant offset.
            const double mono = 12000;
            const int charge = 12;
            var env = BuildPerfectEnvelope(mono, charge);

            double predicted = mono + MostAbundantOffset(mono);
            Assert.That(env.MostAbundantObservedNeutralMass, Is.EqualTo(predicted).Within(0.15));
        }

        [Test]
        public void DeconvolutionConstructor_ComputesMostAbundantObservedNeutralMass()
        {
            // The 5-arg mzLib-deconvolution constructor must compute the most-abundant observed mass.
            const double mono = 12000;
            const int charge = 12;
            var peaks = BuildPerfectPeaks(mono, charge);
            var env = new IsotopicEnvelope(peaks, mono, charge, peaks.Sum(p => p.intensity), 0.5);

            double mostIntenseMz = peaks.MaxBy(p => p.intensity).mz;
            Assert.That(env.MostAbundantObservedNeutralMass, Is.EqualTo(mostIntenseMz.ToMass(charge)).Within(1e-6));
        }

        [Test]
        public void FileReadEnvelope_HasNoMostAbundantPeak_ReturnsSentinel()
        {
            // The file-read constructor carries a neutral mass but no observed isotopic envelope, so the
            // most-abundant observed mass is undefined: both the m/z×|charge| form and the proton-corrected
            // form report the -1 sentinel rather than a synthetic value.
            const double mono = 8000;
            const int charge = 8;
            var env = new IsotopicEnvelope(mono, 1e6, charge);

            Assert.That(env.MostAbundantObservedIsotopicMass, Is.EqualTo(-1));
            Assert.That(env.MostAbundantObservedNeutralMass, Is.EqualTo(-1));
        }
    }
}
