using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Issue #245: the isotope walk only ever examined the single nearest peak to the sought isotope, so an
    /// interferer nearer to the theoretical m/z masked a usable isotope sitting just behind it. Because the walk
    /// breaks on the first miss, that did not merely swap one isotope - it truncated the envelope and shifted the
    /// reported monoisotopic mass by a whole isotope.
    /// </summary>
    public class TestDeconvolutionPeakChoiceWithinTolerance
    {
        private const int Charge = 2;
        private const double TolerancePpm = 20;
        private const double IntensityRatioLimit = 3;
        private static readonly MzRange WholeSpectrum = new MzRange(0, 5000);

        /// <summary>
        /// An exact averagine envelope for <paramref name="mostIntenseMass"/>, with the second-most-intense
        /// isotope displaced by <paramref name="realPeakPpmError"/> so that an interferer can be placed nearer
        /// to the theoretical m/z than the real peak is.
        /// </summary>
        private static (MzSpectrum Spectrum, double TargetMz, double ExpectedIntensity) Build(
            double mostIntenseMass, double realPeakPpmError, double interfererPpmError, double interfererFoldIntensity)
        {
            var averagine = new Averagine();
            int massIndex = averagine.GetMostIntenseMassIndex(mostIntenseMass);
            double[] theorMasses = averagine.GetAllTheoreticalMasses(massIndex);
            double[] theorIntensities = averagine.GetAllTheoreticalIntensities(massIndex);
            double shift = mostIntenseMass - theorMasses[0];

            double targetMz = (theorMasses[1] + shift).ToMz(Charge);
            double expectedIntensity = 1e6 * theorIntensities[1] / theorIntensities[0];

            var peaks = new List<(double Mz, double Intensity)>();
            for (int i = 0; i < Math.Min(6, theorMasses.Length); i++)
            {
                double mz = i == 1
                    ? targetMz * (1 + realPeakPpmError / 1e6)
                    : (theorMasses[i] + shift).ToMz(Charge);
                peaks.Add((mz, 1e6 * theorIntensities[i] / theorIntensities[0]));
            }

            if (interfererFoldIntensity > 0)
            {
                peaks.Add((targetMz * (1 + interfererPpmError / 1e6), expectedIntensity * interfererFoldIntensity));
            }

            var ordered = peaks.OrderBy(p => p.Mz).ToList();
            return (new MzSpectrum(ordered.Select(p => p.Mz).ToArray(), ordered.Select(p => p.Intensity).ToArray(), false),
                targetMz, expectedIntensity);
        }

        private static string[] Describe(IEnumerable<IsotopicEnvelope> envelopes) => envelopes
            .Select(e => $"z={e.Charge} mono={e.MonoisotopicMass:F4} peaks={e.Peaks.Count}")
            .OrderBy(s => s)
            .ToArray();

        /// <summary>
        /// An interferer that is nearer to the theoretical m/z but has the wrong intensity must not displace
        /// the real isotope, whether it is far too intense or far too weak.
        /// </summary>
        [TestCase(20.0)]
        [TestCase(1.0 / 20.0)]
        public void AnInterfererNearerToTheoreticalDoesNotDisplaceTheRealIsotope(double interfererFoldIntensity)
        {
            var parameters = new ClassicDeconvolutionParameters(1, 6, TolerancePpm, IntensityRatioLimit);

            var clean = Build(1600.0, realPeakPpmError: 12, interfererPpmError: 0, interfererFoldIntensity: 0);
            var withInterferer = Build(1600.0, realPeakPpmError: 12, interfererPpmError: -4, interfererFoldIntensity);

            string[] expected = Describe(Deconvoluter.Deconvolute(clean.Spectrum, parameters, WholeSpectrum));
            string[] actual = Describe(Deconvoluter.Deconvolute(withInterferer.Spectrum, parameters, WholeSpectrum));

            Assert.That(expected, Is.Not.Empty, "the clean envelope must deconvolute, or the test proves nothing");
            Assert.That(expected.Single(), Does.Contain("mono=1600.0000").And.Contain("peaks=6"));
            Assert.That(actual, Is.EqualTo(expected));
        }

        /// <summary>
        /// The surviving peak feeds the monoisotopic mass estimate, so the envelope must contain the real
        /// isotope's m/z and not the interferer's.
        /// </summary>
        [Test]
        public void TheEnvelopeContainsTheRealIsotopeNotTheNearerInterferer()
        {
            var parameters = new ClassicDeconvolutionParameters(1, 6, TolerancePpm, IntensityRatioLimit);
            var built = Build(1600.0, realPeakPpmError: 12, interfererPpmError: -4, interfererFoldIntensity: 20);

            double realMz = built.TargetMz * (1 + 12e-6);
            double interfererMz = built.TargetMz * (1 - 4e-6);

            var envelope = Deconvoluter.Deconvolute(built.Spectrum, parameters, WholeSpectrum)
                .OrderByDescending(e => e.Peaks.Count)
                .First();

            Assert.That(envelope.Peaks.Select(p => p.mz), Has.Member(realMz));
            Assert.That(envelope.Peaks.Select(p => p.mz), Has.No.Member(interfererMz));
            Assert.That(envelope.MonoisotopicMass, Is.EqualTo(1600.0).Within(0.01));
        }

        /// <summary>
        /// With only one peak in the window the choice is unambiguous, so nothing about the single-candidate
        /// case may change - including when the sole in-window peak fails the intensity ratio and the walk
        /// must still stop there.
        /// </summary>
        [Test]
        public void ASolePeakFailingTheIntensityRatioStillEndsTheWalk()
        {
            var parameters = new ClassicDeconvolutionParameters(1, 6, TolerancePpm, IntensityRatioLimit);
            var built = Build(1600.0, realPeakPpmError: 0, interfererPpmError: 0, interfererFoldIntensity: 0);

            // replace the second-most-intense isotope with a peak 20x too intense; no alternative exists
            var mz = built.Spectrum.XArray.ToArray();
            var intensities = built.Spectrum.YArray.ToArray();
            int targetIndex = Array.FindIndex(mz, x => Math.Abs(x - built.TargetMz) < 1e-6);
            Assert.That(targetIndex, Is.GreaterThanOrEqualTo(0));
            intensities[targetIndex] = built.ExpectedIntensity * 20;

            var envelopes = Deconvoluter
                .Deconvolute(new MzSpectrum(mz, intensities, false), parameters, WholeSpectrum)
                .ToList();

            Assert.That(envelopes.Any(e => Math.Abs(e.MonoisotopicMass - 1600.0) < 0.01
                                           && e.Peaks.Count == 6), Is.False,
                "the full ladder must not be recovered when the only in-tolerance peak has the wrong intensity");
        }
    }
}
