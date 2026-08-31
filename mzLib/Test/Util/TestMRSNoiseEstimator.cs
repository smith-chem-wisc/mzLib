using System;
using System.Diagnostics.CodeAnalysis;
using MzLibUtil.NoiseEstimation;
using NUnit.Framework;

namespace Test.Util
{
    /// <summary>
    /// MRSNoiseEstimator had no direct test. Its only caller is SpectralWeighting, so every mutant in
    /// the wavelet path had to be detected -- if at all -- through a spectral averaging assertion several
    /// layers up. These are characterization tests: they pin what master computes today on a fixed,
    /// deterministic input, so that a change to the numerics has to be deliberate.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class TestMRSNoiseEstimator
    {
        /// <summary>
        /// Measured on master. Pinned, not derived -- see the class comment.
        /// </summary>
        private const double ExpectedNoise = 6.8650750184874045d;

        /// <summary>
        /// A deterministic signal: a smooth baseline plus a repeating sawtooth ripple. No RNG, so the
        /// expected value below is stable across platforms and runtimes.
        /// </summary>
        private static double[] BuildSignal(int length)
        {
            double[] signal = new double[length];
            for (int i = 0; i < length; i++)
            {
                signal[i] = 10d * Math.Sin(i * 0.05d) + ((i * 37 % 17) - 8) * 0.1d;
            }
            return signal;
        }

        [Test]
        public static void TestNoiseEstimateOfDeterministicSignal()
        {
            double[] signal = BuildSignal(512);

            bool converged = MRSNoiseEstimator.MRSNoiseEstimation(signal, 0.01, out double noiseEstimate);

            Assert.That(converged, Is.True);
            Assert.That(noiseEstimate, Is.EqualTo(ExpectedNoise).Within(1e-6));
        }

        [Test]
        public static void TestNoiseEstimateScalesWithSignalAmplitude()
        {
            double[] signal = BuildSignal(512);
            double[] amplified = new double[signal.Length];
            for (int i = 0; i < signal.Length; i++)
            {
                amplified[i] = signal[i] * 3d;
            }

            MRSNoiseEstimator.MRSNoiseEstimation(signal, 0.01, out double noiseEstimate);
            MRSNoiseEstimator.MRSNoiseEstimation(amplified, 0.01, out double amplifiedNoiseEstimate);

            // the estimator is linear in the signal, so tripling the signal triples the noise estimate
            Assert.That(amplifiedNoiseEstimate, Is.EqualTo(3d * noiseEstimate).Within(1e-6));
        }
    }
}
