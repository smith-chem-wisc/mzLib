using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Two peaks sharing an m/z make the isotope spacing exactly zero, and both classic charge-state
    /// loops divide by that spacing. These tests pin the guard to skipping the duplicate neighbour
    /// rather than ending the scan: the loop's else branch is a break, so guarding in the outer
    /// condition instead of inside it silently truncates charge-state discovery for that candidate.
    /// </summary>
    public class TestDeconvolutionDuplicateMz
    {
        private static readonly MzRange WholeSpectrum = new MzRange(0, 5000);

        /// <summary>An isotopic ladder with the monoisotopic peak most intense.</summary>
        private static (double[] Mz, double[] Intensities) Ladder(double monoisotopicMass, int charge, int peaks)
        {
            var mz = new List<double>();
            var intensities = new List<double>();
            for (int i = 0; i < peaks; i++)
            {
                mz.Add((monoisotopicMass + i * Constants.C13MinusC12).ToMz(charge));
                intensities.Add(1000.0 / (i + 1));
            }
            return (mz.ToArray(), intensities.ToArray());
        }

        private static MzSpectrum Spectrum((double[] Mz, double[] Intensities) ladder) =>
            new MzSpectrum(ladder.Mz.ToArray(), ladder.Intensities.ToArray(), false);

        private static MzSpectrum SpectrumWithDuplicateAt((double[] Mz, double[] Intensities) ladder, int index)
        {
            var mz = ladder.Mz.ToList();
            var intensities = ladder.Intensities.ToList();
            mz.Insert(index + 1, ladder.Mz[index]);
            intensities.Insert(index + 1, ladder.Intensities[index] * 0.9);
            return new MzSpectrum(mz.ToArray(), intensities.ToArray(), false);
        }

        private static string[] Describe(IEnumerable<IsotopicEnvelope> envelopes) => envelopes
            .Select(e => $"z={e.Charge} mono={e.MonoisotopicMass:F4} peaks={e.Peaks.Count} intensity={e.TotalIntensity:F2}")
            .ToArray();

        /// <summary>
        /// The live path, reached by every real caller through <see cref="Deconvoluter"/>.
        /// </summary>
        [TestCase(0)]
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(3)]
        public void ClassicDeconvolutionIsUnaffectedByADuplicatedMz(int duplicateAt)
        {
            var ladder = Ladder(700.0, 1, 4);
            var parameters = new ClassicDeconvolutionParameters(1, 12, 20, 3);

            string[] clean = Describe(Deconvoluter.Deconvolute(Spectrum(ladder), parameters, WholeSpectrum));
            string[] duplicated = Describe(Deconvoluter.Deconvolute(
                SpectrumWithDuplicateAt(ladder, duplicateAt), parameters, WholeSpectrum));

            Assert.That(clean, Is.Not.Empty, "the ladder itself must deconvolute, or the test proves nothing");
            Assert.That(duplicated, Is.EqualTo(clean));
        }

        /// <summary>
        /// The obsolete copy in <see cref="MzSpectrum"/>, which duplicates the same loop and is where
        /// guarding in the outer condition drops the most intense peak out of the envelope.
        /// </summary>
        [TestCase(0)]
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(3)]
        public void ObsoleteMzSpectrumDeconvolutionIsUnaffectedByADuplicatedMz(int duplicateAt)
        {
            var ladder = Ladder(700.0, 1, 4);

#pragma warning disable CS0618 // deliberately covering the obsolete overload; it clones the live loop
            string[] clean = Describe(Spectrum(ladder).Deconvolute(WholeSpectrum, 1, 12, 20, 3));
            string[] duplicated = Describe(
                SpectrumWithDuplicateAt(ladder, duplicateAt).Deconvolute(WholeSpectrum, 1, 12, 20, 3));
#pragma warning restore CS0618

            Assert.That(clean, Is.Not.Empty, "the ladder itself must deconvolute, or the test proves nothing");
            Assert.That(duplicated, Is.EqualTo(clean));
        }
    }
}
