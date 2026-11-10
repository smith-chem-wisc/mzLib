using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;

namespace Test.Omics.Modifications
{
    /// <summary>
    /// <see cref="IsobaricMassTag.GetReporterIonSignalToNoise"/>: each reporter ion's intensity, as
    /// <see cref="IsobaricMassTag.GetReporterIonIntensities"/> finds it, over the scan's noise there.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class IsobaricMassTagSignalToNoiseTests
    {
        // Noise 1100 and baseline 100 at m/z 100, rising linearly to 2100 and 200 at m/z 200.
        private static readonly double[,] NoiseTable = MsDataScan.CreateNoiseData(new[] { 100.0, 200.0 }, new[] { 1100.0, 2100.0 }, new[] { 100.0, 200.0 });

        /// <summary>The signal-to-noise the table above gives, written out independently of the code under test.</summary>
        private static double Expected(double mz, double intensity)
        {
            double noise = 1100 + (mz - 100) * 10;
            double baseline = 100 + (mz - 100);
            return (intensity - baseline) / (noise - baseline);
        }

        private static MsDataScan Scan(double[] mzs, double[] intensities, double[,] noiseData) =>
            new MsDataScan(new MzSpectrum(mzs, intensities, false), 2, 3, true, Polarity.Positive, 1, new MzRange(100, 500),
                "FTMS + c NSI sps ms3", MZAnalyzerType.Orbitrap, intensities.Sum(), null, noiseData, "scan=2");

        private static IsobaricMassTag Tmt10()
        {
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(IsobaricMassTagType.TMT10, out IsobaricMassTag tag), Is.True);
            return tag;
        }

        [Test]
        public void EachReporterIonsSignalToNoiseIsItsIntensityOverTheNoiseThere()
        {
            IsobaricMassTag tag = Tmt10();
            double[] intensities = Enumerable.Range(1, tag.ReporterIonMzs.Length).Select(k => 10000.0 * k).ToArray();
            MsDataScan scan = Scan(tag.ReporterIonMzs.ToArray(), intensities, NoiseTable);

            double[] signalToNoise = tag.GetReporterIonSignalToNoise(scan);

            Assert.That(signalToNoise, Has.Length.EqualTo(tag.ReporterIonMzs.Length));
            for (int i = 0; i < signalToNoise.Length; i++)
            {
                Assert.That(signalToNoise[i], Is.EqualTo(Expected(tag.ReporterIonMzs[i], intensities[i])).Within(1e-12).Percent, tag.ChannelLabels[i]);
            }
        }

        [Test]
        public void TheIntensityIsTheOneGetReporterIonIntensitiesFindsAndAMissingReporterIonIsZero()
        {
            IsobaricMassTag tag = Tmt10();
            double mz126 = tag.ReporterIonMzs[0];
            double mz127N = tag.ReporterIonMzs[1];
            // Two peaks within tolerance of 126: the more intense one counts, as for the intensities. 127N has no peak.
            MsDataScan scan = Scan(new[] { mz126 - 0.001, mz126 + 0.002, mz127N + 0.01 }, new[] { 30000.0, 50000.0, 90000.0 }, NoiseTable);

            double[] intensities = tag.GetReporterIonIntensities(scan.MassSpectrum);
            double[] signalToNoise = tag.GetReporterIonSignalToNoise(scan);

            Assert.That(intensities[0], Is.EqualTo(50000.0));
            Assert.That(signalToNoise[0], Is.EqualTo(Expected(mz126, 50000.0)).Within(1e-12).Percent);
            Assert.That(intensities[1], Is.EqualTo(0.0));
            Assert.That(signalToNoise[1], Is.EqualTo(0.0));
            Assert.That(signalToNoise.Skip(2), Is.All.EqualTo(0.0));
        }

        [Test]
        public void NoScanNoNoiseTableOrNoSpectrumGivesNull()
        {
            IsobaricMassTag tag = Tmt10();
            double[] mzs = tag.ReporterIonMzs.ToArray();
            double[] intensities = mzs.Select(_ => 1000.0).ToArray();

            Assert.That(tag.GetReporterIonSignalToNoise(null), Is.Null);
            Assert.That(tag.GetReporterIonSignalToNoise(Scan(mzs, intensities, null)), Is.Null, "an ion-trap scan has no noise table");
            Assert.That(tag.GetReporterIonSignalToNoise(Scan(new double[0], new double[0], NoiseTable)), Is.Null);
        }

        [Test]
        public void AReporterIonWhereNoiseDoesNotExceedBaselineIsNaN()
        {
            IsobaricMassTag tag = Tmt10();
            double[,] flat = MsDataScan.CreateNoiseData(new[] { 100.0, 200.0 }, new[] { 100.0, 100.0 }, new[] { 100.0, 100.0 });
            MsDataScan scan = Scan(new[] { tag.ReporterIonMzs[0] }, new[] { 5000.0 }, flat);

            double[] signalToNoise = tag.GetReporterIonSignalToNoise(scan);

            Assert.That(signalToNoise[0], Is.NaN);
            Assert.That(signalToNoise.Skip(1), Is.All.EqualTo(0.0), "no peak is zero, as its intensity is");
        }
    }
}
