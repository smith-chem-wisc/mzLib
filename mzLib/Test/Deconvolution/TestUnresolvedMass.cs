using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Unit tests for the isotopically-<b>unresolved</b> (high-mass) precursor-selection support: the
    /// average (centroid) mass path that the resolved most-abundant feature defers to.
    ///
    /// Terminology pinned here:
    ///  • average / centroid mass = intensity-weighted mean neutral mass over the whole envelope
    ///    (<see cref="IsotopicEnvelope.AverageObservedMass"/>), the abundance-weighted counterpart to the
    ///    monoisotopic mass — distinct from the most-abundant (tallest single peak) mass;
    ///  • average / centroid offset = <see cref="AverageResidue.GetAverageOffset"/>, the gap from the
    ///    monoisotopic mass to that centroid;
    ///  • <see cref="EnvelopeResolution"/> selects which observed mass candidate selection uses.
    /// All tests build synthetic envelopes from the Averagine model — no deconvolution is run.
    /// </summary>
    [TestFixture]
    public sealed class TestUnresolvedMass
    {
        private static readonly AverageResidue Model = new Averagine();

        private static double MostAbundantOffset(double monoMass) => Model.GetDiffToMonoisotopic(Model.GetMostIntenseMassIndex(monoMass));

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

        // ── Average (centroid) offset ───────────────────────────────────────────

        [Test]
        public void AverageOffset_IsAtLeastMostAbundantOffset()
        {
            // The intensity-weighted mean of a protein isotope envelope sits at or above the
            // most-abundant isotopologue.
            double[] masses = { 2000, 10000, 30000 };
            foreach (double m in masses)
            {
                double avgOffset = Model.GetAverageOffset(m);
                double mostAbundant = MostAbundantOffset(m);
                Assert.That(avgOffset, Is.GreaterThanOrEqualTo(mostAbundant - 1e-6), $"at mass {m}");
                // centroid sits just above the most-abundant peak (≲1 Da here); a sort-order/sign
                // regression would inflate the offset by ~the most-abundant offset itself (many Da).
                Assert.That(avgOffset, Is.LessThanOrEqualTo(mostAbundant + 2.0), $"inflated average offset at mass {m}");
            }
        }

        // ── Average (centroid) observed mass ────────────────────────────────────

        [Test]
        public void AverageObservedMass_MatchesMonoPlusAveragineOffset()
        {
            // The observed centroid neutral mass of a perfect envelope ≈ candidate monoisotopic +
            // averagine average offset.
            const double mono = 12000;
            const int charge = 12;
            var env = BuildPerfectEnvelope(mono, charge);

            double predicted = mono + Model.GetAverageOffset(mono);
            Assert.That(env.AverageObservedMass, Is.EqualTo(predicted).Within(0.15));
        }

        [Test]
        public void AverageObservedMass_LiesWithinTheEnvelope()
        {
            const double mono = 12000;
            const int charge = 12;
            var env = BuildPerfectEnvelope(mono, charge);

            Assert.That(env.AverageObservedMass, Is.GreaterThanOrEqualTo(mono));
            // at/above the most-abundant peak for these right-skewed envelopes...
            Assert.That(env.AverageObservedMass, Is.GreaterThanOrEqualTo(env.MostAbundantObservedNeutralMass - 1e-6));
            // ...and never above the heaviest peak (a charge/weighting regression would inflate it).
            double heaviestNeutralMass = env.Peaks.Max(p => p.mz).ToMass(charge);
            Assert.That(env.AverageObservedMass, Is.LessThanOrEqualTo(heaviestNeutralMass + 1e-6));
        }

        [Test]
        public void AverageObservedMass_ZeroTotalIntensity_FallsBackToMostIntensePeakMass()
        {
            // Degenerate envelope (all intensities zero): AverageObservedMass must take the
            // totalIntensity == 0 fallback (most-intense peak's m/z) and return a real mass, not NaN.
            const int charge = 1;
            var peaks = new List<(double mz, double intensity)> { (1000.0, 0.0), (1000.5, 0.0) };
            var env = new IsotopicEnvelope(0, peaks, 999.0, charge, 0.0, 0.5);

            Assert.That(double.IsNaN(env.AverageObservedMass), Is.False);
            Assert.That(env.AverageObservedMass, Is.EqualTo(env.MostAbundantObservedNeutralMass).Within(1e-6));
        }

        [Test]
        public void FileReadEnvelope_HasNoAverageMass_ReturnsSentinel()
        {
            // A neutral mass read from a file has no observed envelope, so the centroid is undefined and
            // reports the same -1 sentinel as the most-abundant masses.
            var env = new IsotopicEnvelope(8000, 1e6, 8);
            Assert.That(env.AverageObservedMass, Is.EqualTo(-1));
        }

        // ── Resolution state ────────────────────────────────────────────────────

        [Test]
        public void Resolution_DefaultsToResolved_AndCanBeSet()
        {
            var env = BuildPerfectEnvelope(10000, 10);
            Assert.That(env.Resolution, Is.EqualTo(EnvelopeResolution.Resolved));

            env.SetResolution(EnvelopeResolution.Unresolved);
            Assert.That(env.Resolution, Is.EqualTo(EnvelopeResolution.Unresolved));
        }

        [Test]
        public void Resolution_ParticipatesInEquality()
        {
            var a = BuildPerfectEnvelope(10000, 10);
            var b = BuildPerfectEnvelope(10000, 10);
            Assert.That(a, Is.EqualTo(b));

            b.SetResolution(EnvelopeResolution.Unresolved);
            Assert.That(a, Is.Not.EqualTo(b));
        }
    }
}
