using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Unit tests for the "most abundant mass" precursor-selection support (Strategy B):
    /// the averagine apex offset (diff-to-monoisotopic of the nearest mass bin) and
    /// <see cref="AverageResidue.GetAverageOffset"/>, and the proton-corrected
    /// <see cref="IsotopicEnvelope.MostAbundantObservedMass"/> /
    /// <see cref="IsotopicEnvelope.AverageObservedMass"/> / <see cref="IsotopicEnvelope.Resolution"/>.
    /// All tests build synthetic envelopes from the Averagine model — no deconvolution is run.
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

        // ── AverageResidue offsets ──────────────────────────────────────────────

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

        // ── IsotopicEnvelope observed masses ────────────────────────────────────

        [Test]
        public void MostAbundantObservedMass_IsProtonCorrectedNeutralMass()
        {
            const int charge = 10;
            var env = BuildPerfectEnvelope(15000, charge);

            // The proton-corrected neutral mass equals the most intense peak's m/z .ToMass(charge)...
            double mostIntenseMz = env.Peaks.MaxBy(p => p.intensity).mz;
            Assert.That(env.MostAbundantObservedMass, Is.EqualTo(mostIntenseMz.ToMass(charge)).Within(1e-6));

            // ...and it differs from the un-proton-corrected mz*|z| field by exactly z proton masses.
            Assert.That(env.MostAbundantObservedIsotopicMass - env.MostAbundantObservedMass,
                Is.EqualTo(charge * Constants.ProtonMass).Within(1e-6));
        }

        [Test]
        public void MostAbundantObservedMass_MatchesMonoPlusAveragineOffset()
        {
            // Ties the two features together: the observed most-abundant neutral mass of a perfect
            // envelope ≈ candidate monoisotopic + averagine most-abundant offset.
            const double mono = 12000;
            const int charge = 12;
            var env = BuildPerfectEnvelope(mono, charge);

            double predicted = mono + MostAbundantOffset(mono);
            Assert.That(env.MostAbundantObservedMass, Is.EqualTo(predicted).Within(0.15));
        }

        [Test]
        public void AverageObservedMass_MatchesMonoPlusAveragineOffset()
        {
            // Symmetric to MostAbundantObservedMass_MatchesMonoPlusAveragineOffset: the observed centroid
            // neutral mass of a perfect envelope ≈ candidate monoisotopic + averagine average offset.
            const double mono = 12000;
            const int charge = 12;
            var env = BuildPerfectEnvelope(mono, charge);

            double predicted = mono + Model.GetAverageOffset(mono);
            Assert.That(env.AverageObservedMass, Is.EqualTo(predicted).Within(0.15));
        }

        [Test]
        public void AverageObservedMass_IsBetweenMonoAndAboveMostAbundant()
        {
            const double mono = 12000;
            const int charge = 12;
            var env = BuildPerfectEnvelope(mono, charge);

            Assert.That(env.AverageObservedMass, Is.GreaterThanOrEqualTo(mono));
            // centroid sits at/above the most-abundant peak for these right-skewed envelopes
            Assert.That(env.AverageObservedMass, Is.GreaterThanOrEqualTo(env.MostAbundantObservedMass - 1e-6));
            // ...and the "between" claim needs an upper bound: the centroid cannot exceed the heaviest
            // peak in the envelope (a charge/weighting regression would inflate it past the envelope span).
            double heaviestNeutralMass = env.Peaks.Max(p => p.mz).ToMass(charge);
            Assert.That(env.AverageObservedMass, Is.LessThanOrEqualTo(heaviestNeutralMass + 1e-6));
        }

        [Test]
        public void AverageObservedMass_ZeroTotalIntensity_FallsBackToMostIntensePeakMass()
        {
            // Degenerate envelope (all intensities zero): AverageObservedMass must take the
            // totalIntensity == 0 fallback (most-intense peak's m/z) and return a real mass, not NaN
            // from a divide-by-zero. Exercises the fallback branch in the centroid computation.
            const int charge = 1;
            var peaks = new List<(double mz, double intensity)> { (1000.0, 0.0), (1000.5, 0.0) };
            var env = new IsotopicEnvelope(0, peaks, 999.0, charge, 0.0, 0.5);

            Assert.That(double.IsNaN(env.AverageObservedMass), Is.False);
            Assert.That(env.AverageObservedMass, Is.EqualTo(env.MostAbundantObservedMass).Within(1e-6));
        }

        [Test]
        public void DeconvolutionConstructor_ComputesObservedMasses()
        {
            // The 5-arg mzLib-deconvolution constructor must compute the same observed masses as the 6-arg path.
            const double mono = 12000;
            const int charge = 12;
            var peaks = BuildPerfectPeaks(mono, charge);
            var env = new IsotopicEnvelope(peaks, mono, charge, peaks.Sum(p => p.intensity), 0.5);

            double mostIntenseMz = peaks.MaxBy(p => p.intensity).mz;
            Assert.That(env.MostAbundantObservedMass, Is.EqualTo(mostIntenseMz.ToMass(charge)).Within(1e-6));
            Assert.That(env.AverageObservedMass, Is.GreaterThanOrEqualTo(mono));
            Assert.That(env.AverageObservedMass, Is.EqualTo(mono + Model.GetAverageOffset(mono)).Within(0.15));
        }

        [Test]
        public void SinglePeakEnvelope_AllMassesEqualMonoisotopic()
        {
            // The file-read constructor produces a one-peak envelope at the monoisotopic mass.
            const double mono = 8000;
            const int charge = 8;
            var env = new IsotopicEnvelope(mono, 1e6, charge);

            Assert.That(env.MostAbundantObservedMass, Is.EqualTo(mono).Within(1e-6));
            Assert.That(env.AverageObservedMass, Is.EqualTo(mono).Within(1e-6));
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
