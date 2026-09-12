using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSpectrumSmoothing
    {
        /// <summary>A unit square wave ten points wide, centred in a fifty-point spectrum.</summary>
        private static double[] SquareWave()
        {
            double[] intensities = new double[50];
            for (int i = 24; i < 34; i++)
            {
                intensities[i] = 1.0;
            }

            return intensities;
        }

        private static MzSpectrum SquareWaveSpectrum()
        {
            double[] intensities = SquareWave();
            double[] mzs = Enumerable.Range(0, intensities.Length).Select(i => (double)i).ToArray();
            return new MzSpectrum(mzs, intensities, true);
        }

        [TestCase(0)]
        [TestCase(1)]
        [TestCase(5)]
        public void IterationsOfZeroLeavesTheDataAlone(int halfWindowSize)
        {
            double[] input = { 5, 4, 3, 2, 1 };

            Assert.That(MzSpectrum.KZ1D(input, halfWindowSize, 0), Is.EqualTo(input),
                "zero passes means no smoothing was asked for, so the input has to come back");
        }

        [TestCase(1)]
        [TestCase(3)]
        public void AHalfWindowOfZeroLeavesTheDataAlone(int iterations)
        {
            double[] input = { 1, 2, 3, 4, 5 };

            Assert.That(MzSpectrum.KZ1D(input, 0, iterations), Is.EqualTo(input),
                "a window of one point is its own average, however many passes are made");
        }

        /// <summary>
        /// The end of the array is where an off-by-one hides: a window clamped to Length - 1 rather than
        /// Length drops the final channel out of every average, including its own.
        /// </summary>
        [Test]
        public void TheFinalChannelIsNotDiscarded()
        {
            double[] allSignalInTheLastChannel = { 0, 0, 0, 0, 100 };

            double[] smoothed = MzSpectrum.KZ1D(allSignalInTheLastChannel, 1, 1);

            Assert.That(smoothed.Sum(), Is.GreaterThan(0.0),
                "the whole spectrum's intensity sat in the final channel and has gone missing");
            Assert.That(smoothed[4], Is.EqualTo(50.0).Within(1e-9), "mean of {0, 100}");
            Assert.That(smoothed[3], Is.EqualTo(100.0 / 3).Within(1e-9), "mean of {0, 0, 100}");
        }

        [Test]
        public void ASingleChannelSpectrumIsReturnedUnchanged()
        {
            Assert.That(MzSpectrum.KZ1D(new double[] { 7 }, 1, 1), Is.EqualTo(new double[] { 7 }));
        }

        /// <summary>
        /// The reference skips non-finite values rather than summing them, so one NaN costs its own
        /// channel and not every channel whose window overlaps it.
        /// </summary>
        [Test]
        public void ANonFiniteValueDoesNotPoisonItsNeighbours()
        {
            double[] withAHole = { 1, 1, double.NaN, 1, 1 };

            double[] smoothed = MzSpectrum.KZ1D(withAHole, 1, 1);

            Assert.That(smoothed, Has.None.Matches<double>(double.IsNaN),
                "a single NaN spread across the whole spectrum");
            Assert.That(smoothed, Is.All.EqualTo(1.0).Within(1e-9));
        }

        /// <summary>The one case where NaN is the right answer: nothing finite to average.</summary>
        [Test]
        public void AWindowWithNoFiniteValueIsNaN()
        {
            double[] nothingFinite = { double.NaN, double.PositiveInfinity, double.NaN };

            Assert.That(MzSpectrum.KZ1D(nothingFinite, 1, 1), Is.All.Matches<double>(double.IsNaN));
        }

        [Test]
        public void NullIntensitiesAreRefused()
        {
            var ex = Assert.Throws<MzLibException>(() => MzSpectrum.KZ1D(null, 1, 1));
            Assert.That(ex!.Message, Does.Contain("was given null"));
        }

        [Test]
        public void ANegativeHalfWindowSizeIsRefused()
        {
            var ex = Assert.Throws<MzLibException>(() => MzSpectrum.KZ1D(new double[] { 1, 2 }, -1, 1));
            Assert.That(ex!.Message, Does.Contain("half-window size"));
            Assert.That(ex.Message, Does.Contain("-1"), "the offending value belongs in the message");
        }

        [Test]
        public void ANegativeIterationCountIsRefused()
        {
            var ex = Assert.Throws<MzLibException>(() => MzSpectrum.KZ1D(new double[] { 1, 2 }, 1, -3));
            Assert.That(ex!.Message, Does.Contain("iteration count"));
            Assert.That(ex.Message, Does.Contain("-3"));
        }

        [Test]
        public void TheInputArrayIsNotModified()
        {
            double[] input = { 1, 9, 1 };
            double[] copy = (double[])input.Clone();

            MzSpectrum.KZ1D(input, 1, 3);

            Assert.That(input, Is.EqualTo(copy));
        }

        /// <summary>
        /// While the smeared signal stays clear of both ends of the array, every point's window holds a
        /// full complement of neighbours and the filter only moves intensity around.
        /// </summary>
        [TestCase(3, 2)]
        [TestCase(5, 2)]
        public void SmoothingConservesTotalIntensityAwayFromTheEnds(int halfWindowSize, int iterations)
        {
            double[] wave = SquareWave();

            double[] smoothed = MzSpectrum.KZ1D(wave, halfWindowSize, iterations);

            Assert.That(smoothed.Sum(), Is.EqualTo(wave.Sum()).Within(1e-9));
        }

        /// <summary>
        /// And once the window is clamped at an end it averages fewer points, so the total drifts. This
        /// is the reference implementation's behaviour, not a defect, but it is the reason the output is
        /// not quantitative near the first or last few points -- so it is pinned rather than left to be
        /// rediscovered.
        /// </summary>
        [Test]
        public void ClampingAtTheEndsDoesNotConserveTotalIntensity()
        {
            double[] wave = SquareWave();

            // Wide enough, and iterated enough, that the smeared wave reaches the ends of the array.
            double reachesTheEnds = MzSpectrum.KZ1D(wave, 5, 4).Sum();

            Assert.That(reachesTheEnds, Is.Not.EqualTo(wave.Sum()).Within(1e-9));
            Assert.That(reachesTheEnds, Is.EqualTo(10.0236).Within(1e-4),
                "drift is deterministic; if this number moves, the edge handling changed");
        }

        /// <summary>
        /// Pins the trade-off the remarks describe, so a later change cannot quietly alter it: this
        /// filter buys smoothness with peak height.
        /// </summary>
        [Test]
        public void AWiderWindowLowersThePeakFurther()
        {
            double[] wave = SquareWave();

            double unsmoothed = wave.Max();
            double narrow = MzSpectrum.KZ1D(wave, 3, 2).Max();
            double wide = MzSpectrum.KZ1D(wave, 5, 2).Max();

            Assert.That(unsmoothed, Is.EqualTo(1.0));
            Assert.That(narrow, Is.LessThan(unsmoothed));
            Assert.That(wide, Is.LessThan(narrow));
        }

        [Test]
        public void AWindowWiderThanTheSpectrumGivesTheMeanOfTheWholeScan()
        {
            double[] wave = SquareWave();
            double mean = wave.Average();

            double[] smoothed = MzSpectrum.KZ1D(wave, wave.Length + 1, 2);

            Assert.That(smoothed, Is.All.EqualTo(mean).Within(1e-9));
        }

        [Test]
        public void SmoothSpectrumKZKeepsTheMzValuesAndSmoothsRowOne()
        {
            MzSpectrum spectrum = SquareWaveSpectrum();

            double[,] result = spectrum.SmoothSpectrumKZ(5, 2);

            Assert.That(result.GetLength(0), Is.EqualTo(2));
            Assert.That(result.GetLength(1), Is.EqualTo(50));

            for (int i = 0; i < 50; i++)
            {
                Assert.That(result[0, i], Is.EqualTo(spectrum.XArray[i]), $"m/z moved at index {i}");
            }

            Assert.That(result[1, 0], Is.EqualTo(0.0).Within(1e-9),
                "far from the wave, so still zero");
            Assert.That(result[1, 23], Is.GreaterThan(0.0),
                "just outside the wave, so smoothing has carried signal into it");
        }
    }
}
