using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Readers;

namespace Test.Omics.Modifications
{
    /// <summary>
    /// <see cref="IsobaricMassTag.GetReporterIonReadings"/>: each reporter ion as a value, as below the scan's detection
    /// floor, or as not recorded, with the floor where one is known.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class IsobaricMassTagReadingTests
    {
        // Noise 1100 at m/z 100 rising linearly to 2100 at m/z 200 (baseline 100 to 200), as in IsobaricMassTagSignalToNoiseTests.
        private static readonly double[,] NoiseTable = MsDataScan.CreateNoiseData(new[] { 100.0, 200.0 }, new[] { 1100.0, 2100.0 }, new[] { 100.0, 200.0 });

        /// <summary>The noise the table above gives at <paramref name="mz"/>, written out independently of the code under test.</summary>
        private static double ExpectedNoise(double mz) => 1100 + (mz - 100) * 10;

        private static MsDataScan Scan(double[] mzs, double[] intensities, MzRange window, double[,] noiseData = null,
            MZAnalyzerType analyzer = MZAnalyzerType.Orbitrap) =>
            new MsDataScan(new MzSpectrum(mzs, intensities, false), 2, 3, true, Polarity.Positive, 1, window,
                "FTMS + c NSI sps ms3", analyzer, intensities.Sum(), null, noiseData, "scan=2");

        private static IsobaricMassTag Tag(IsobaricMassTagType type)
        {
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(type, out IsobaricMassTag tag), Is.True);
            return tag;
        }

        private static void AssertNotRecorded(ReporterIonReading reading, string label)
        {
            Assert.That(reading.Kind, Is.EqualTo(ReporterIonReadingKind.NotRecorded), label);
            Assert.That(reading.Intensity, Is.NaN, label);
            Assert.That(reading.Floor, Is.Null, label);
            Assert.That(reading.FloorSource, Is.EqualTo(DetectionFloorSource.None), label);
        }

        [Test]
        public void AReporterIonBelowTheScanWindowIsNotRecordedAndOneInsideWithoutAPeakIsBelowTheFloor()
        {
            IsobaricMassTag tag = Tag(IsobaricMassTagType.TMT10);
            double[] mz = tag.ReporterIonMzs;
            // The window starts at 128, between 127C (127.131) and 128N (128.128). Peaks at 128N and 129C, and one far away.
            MsDataScan scan = Scan(new[] { mz[3], mz[6], 300.0 }, new[] { 20000.0, 40000.0, 1000.0 }, new MzRange(128, 500));

            ReporterIonReading[] readings = tag.GetReporterIonReadings(scan);
            double[] intensities = tag.GetReporterIonIntensities(scan.MassSpectrum);

            Assert.That(readings, Has.Length.EqualTo(10));
            for (int i = 0; i < 3; i++)
            {
                AssertNotRecorded(readings[i], tag.ChannelLabels[i]);
            }
            for (int i = 3; i < 10; i++)
            {
                string label = tag.ChannelLabels[i];
                // No noise table: the floor is the weakest peak among the reporter ions, 20000 (the 1000 at m/z 300 is outside).
                Assert.That(readings[i].Floor, Is.EqualTo(20000.0), label);
                Assert.That(readings[i].FloorSource, Is.EqualTo(DetectionFloorSource.ReporterRegionMinimum), label);
                if (i == 3 || i == 6)
                {
                    Assert.That(readings[i].Kind, Is.EqualTo(ReporterIonReadingKind.Value), label);
                    Assert.That(readings[i].Intensity, Is.EqualTo(intensities[i]), label);
                }
                else
                {
                    Assert.That(readings[i].Kind, Is.EqualTo(ReporterIonReadingKind.BelowFloor), label);
                    Assert.That(readings[i].Intensity, Is.NaN, label + ": no value is never 0");
                    Assert.That(intensities[i], Is.EqualTo(0.0), label + ": the intensities still write 0, unchanged");
                }
            }
            Assert.That(readings[3].Intensity, Is.EqualTo(20000.0));
            Assert.That(readings[6].Intensity, Is.EqualTo(40000.0));
        }

        [Test]
        public void TheNoiseTableGivesTheFloorAtEachReporterIon()
        {
            IsobaricMassTag tag = Tag(IsobaricMassTagType.TMT10);
            double[] mz = tag.ReporterIonMzs;
            MsDataScan scan = Scan(new[] { mz[0] }, new[] { 50000.0 }, new MzRange(100, 500), NoiseTable);

            ReporterIonReading[] readings = tag.GetReporterIonReadings(scan);

            Assert.That(readings[0].Kind, Is.EqualTo(ReporterIonReadingKind.Value));
            Assert.That(readings[0].Intensity, Is.EqualTo(50000.0));
            Assert.That(readings.Skip(1).Select(r => r.Kind), Is.All.EqualTo(ReporterIonReadingKind.BelowFloor));
            for (int i = 0; i < readings.Length; i++)
            {
                Assert.That(readings[i].Floor, Is.EqualTo(ExpectedNoise(mz[i])).Within(1e-9), tag.ChannelLabels[i]);
                Assert.That(readings[i].FloorSource, Is.EqualTo(DetectionFloorSource.NoiseTable), tag.ChannelLabels[i]);
            }
        }

        [Test]
        public void WithoutANoiseTableTheFloorIsTheWeakestReporterRegionPeakElseTheWeakestPeakInTheScan()
        {
            IsobaricMassTag tag = Tag(IsobaricMassTagType.TMT10);
            double[] mz = tag.ReporterIonMzs;

            // 129.5 lies between the reporter ions and matches none of them: it counts toward the region's floor. The weaker
            // peaks at 110 and 300 are outside the region and do not.
            MsDataScan withRegionPeak = Scan(new[] { 110.0, mz[0], 129.5, 300.0 }, new[] { 500.0, 9000.0, 1200.0, 100.0 }, new MzRange(100, 500));
            ReporterIonReading[] readings = tag.GetReporterIonReadings(withRegionPeak);
            Assert.That(readings.Select(r => r.Floor), Is.All.EqualTo(1200.0));
            Assert.That(readings.Select(r => r.FloorSource), Is.All.EqualTo(DetectionFloorSource.ReporterRegionMinimum));
            Assert.That(readings[0].Kind, Is.EqualTo(ReporterIonReadingKind.Value));
            Assert.That(readings.Skip(1).Select(r => r.Kind), Is.All.EqualTo(ReporterIonReadingKind.BelowFloor));

            // No peak between the reporter ions: the weakest peak in the scan, 800.
            MsDataScan withoutRegionPeak = Scan(new[] { 110.0, 300.0 }, new[] { 800.0, 5000.0 }, new MzRange(100, 500));
            readings = tag.GetReporterIonReadings(withoutRegionPeak);
            Assert.That(readings.Select(r => r.Kind), Is.All.EqualTo(ReporterIonReadingKind.BelowFloor));
            Assert.That(readings.Select(r => r.Floor), Is.All.EqualTo(800.0));
            Assert.That(readings.Select(r => r.FloorSource), Is.All.EqualTo(DetectionFloorSource.ScanMinimum));
        }

        [Test]
        public void AReporterIonAboveTheScanWindowIsNotRecorded()
        {
            IsobaricMassTag tag = Tag(IsobaricMassTagType.TMT10);
            // The window ends at 127.128, between 127N (127.125) and 127C (127.131). Its only peak, at 110, is outside the
            // reporter region, so the floor is the weakest peak in the scan.
            MsDataScan scan = Scan(new[] { 110.0 }, new[] { 600.0 }, new MzRange(100, 127.128));

            ReporterIonReading[] readings = tag.GetReporterIonReadings(scan);

            for (int i = 0; i < 2; i++)
            {
                Assert.That(readings[i].Kind, Is.EqualTo(ReporterIonReadingKind.BelowFloor), tag.ChannelLabels[i]);
                Assert.That(readings[i].Floor, Is.EqualTo(600.0), tag.ChannelLabels[i]);
                Assert.That(readings[i].FloorSource, Is.EqualTo(DetectionFloorSource.ScanMinimum), tag.ChannelLabels[i]);
            }
            for (int i = 2; i < readings.Length; i++)
            {
                AssertNotRecorded(readings[i], tag.ChannelLabels[i]);
            }
        }

        [Test]
        public void ANoiseTableWithoutAPositiveNoiseLevelGivesWayToTheReporterRegionPeaks()
        {
            IsobaricMassTag tag = Tag(IsobaricMassTagType.TMT10);
            double[,] zeroNoise = MsDataScan.CreateNoiseData(new[] { 100.0, 200.0 }, new[] { 0.0, 0.0 }, new[] { 0.0, 0.0 });
            MsDataScan scan = Scan(new[] { tag.ReporterIonMzs[0] }, new[] { 9000.0 }, new MzRange(100, 500), zeroNoise);

            ReporterIonReading[] readings = tag.GetReporterIonReadings(scan);

            Assert.That(readings.Select(r => r.Floor), Is.All.EqualTo(9000.0));
            Assert.That(readings.Select(r => r.FloorSource), Is.All.EqualTo(DetectionFloorSource.ReporterRegionMinimum));
        }

        [Test]
        public void AScanWhosePeaksAreAllZeroIsBelowTheFloorWithNoFloorToGive()
        {
            IsobaricMassTag tag = Tag(IsobaricMassTagType.TMT10);
            MsDataScan scan = Scan(new[] { tag.ReporterIonMzs[0], 300.0 }, new[] { 0.0, 0.0 }, new MzRange(100, 500));

            ReporterIonReading[] readings = tag.GetReporterIonReadings(scan);

            Assert.That(readings.Select(r => r.Kind), Is.All.EqualTo(ReporterIonReadingKind.BelowFloor));
            Assert.That(readings.Select(r => r.Intensity), Is.All.NaN);
            Assert.That(readings.Select(r => r.Floor), Is.All.Null);
            Assert.That(readings.Select(r => r.FloorSource), Is.All.EqualTo(DetectionFloorSource.None));
        }

        [TestCase(MZAnalyzerType.IonTrap2D)]
        [TestCase(MZAnalyzerType.IonTrap3D)]
        public void AnIonTrapScanCannotSayAReporterIonWasAbsentButAPeakItFindsIsStillAValue(MZAnalyzerType ionTrap)
        {
            IsobaricMassTag tag = Tag(IsobaricMassTagType.TMT10);
            double[] mz = tag.ReporterIonMzs;
            MsDataScan scan = Scan(new[] { mz[0], 300.0 }, new[] { 7000.0, 2000.0 }, new MzRange(100, 2000), analyzer: ionTrap);

            ReporterIonReading[] readings = tag.GetReporterIonReadings(scan);

            Assert.That(readings[0].Kind, Is.EqualTo(ReporterIonReadingKind.Value));
            Assert.That(readings[0].Intensity, Is.EqualTo(7000.0));
            for (int i = 1; i < readings.Length; i++)
            {
                AssertNotRecorded(readings[i], tag.ChannelLabels[i]);
            }
        }

        [Test]
        public void WithoutAScanWindowTheSpanOfThePeaksIsWhatWasRecorded()
        {
            IsobaricMassTag tag = Tag(IsobaricMassTagType.TMT10);
            // The lowest peak is at 127.5, above 126, 127N (127.125) and 127C (127.131).
            MsDataScan scan = Scan(new[] { 127.5, 400.0 }, new[] { 1000.0, 2000.0 }, null);

            ReporterIonReading[] readings = tag.GetReporterIonReadings(scan);

            for (int i = 0; i < 3; i++)
            {
                AssertNotRecorded(readings[i], tag.ChannelLabels[i]);
            }
            Assert.That(readings.Skip(3).Select(r => r.Kind), Is.All.EqualTo(ReporterIonReadingKind.BelowFloor));
            Assert.That(readings.Skip(3).Select(r => r.Floor), Is.All.EqualTo(1000.0));
        }

        [Test]
        public void AnEmptySpectrumRecordedNothingAndNoScanGivesNull()
        {
            IsobaricMassTag tag = Tag(IsobaricMassTagType.TMT10);

            ReporterIonReading[] readings = tag.GetReporterIonReadings(Scan(new double[0], new double[0], new MzRange(100, 500), NoiseTable));

            Assert.That(readings, Has.Length.EqualTo(10));
            for (int i = 0; i < readings.Length; i++)
            {
                AssertNotRecorded(readings[i], tag.ChannelLabels[i]);
            }
            Assert.That(tag.GetReporterIonReadings(null), Is.Null);
        }

        /// <summary>
        /// A real TMT11 SPS-MS3 run. Its 37 MS2 scans are ion-trap scans (ITMS): none can say a reporter ion was absent,
        /// and 33 of them start above the reporter ions anyway. Its 37 MS3 scans are Orbitrap scans from m/z 110. The
        /// counts were derived independently of mzLib, by decoding the mzML in Python and matching the published TMT11
        /// reporter m/z within 0.003 Da: 405 values and 2 reporter ions without a peak, both in scan 28298 (130C and 131C),
        /// whose weakest reporter-region peak is 3743.227.
        /// </summary>
        [Test]
        public void ARealSpsMs3RunsIonTrapMs2sRecordNoReporterIonAndItsMs3sAreValuesOrBelowTheFloor()
        {
            IsobaricMassTag tag = Tag(IsobaricMassTagType.TMT11);
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "MS3_TMT11_Mouse_snip.mzML");
            var file = MsDataFileReader.GetDataFile(path);
            file.LoadAllStaticData();
            var scans = file.GetAllScansList();

            var ms2 = scans.Where(s => s.MsnOrder == 2).ToList();
            var ms2Readings = ms2.SelectMany(s => tag.GetReporterIonReadings(s)).ToList();
            Assert.That(ms2, Has.Count.EqualTo(37));
            Assert.That(ms2.Select(s => s.MzAnalyzer), Is.All.EqualTo(MZAnalyzerType.IonTrap2D));
            Assert.That(ms2Readings.Select(r => r.Kind), Is.All.EqualTo(ReporterIonReadingKind.NotRecorded));
            Assert.That(ms2Readings, Has.Count.EqualTo(37 * 11));

            var ms3 = scans.Where(s => s.MsnOrder == 3).ToList();
            var ms3Readings = ms3.SelectMany(s => tag.GetReporterIonReadings(s)).ToList();
            Assert.That(ms3, Has.Count.EqualTo(37));
            Assert.That(ms3Readings.Count(r => r.Kind == ReporterIonReadingKind.Value), Is.EqualTo(405));
            Assert.That(ms3Readings.Count(r => r.Kind == ReporterIonReadingKind.BelowFloor), Is.EqualTo(2));
            Assert.That(ms3Readings.Count(r => r.Kind == ReporterIonReadingKind.NotRecorded), Is.EqualTo(0));
            Assert.That(ms3Readings.Select(r => r.FloorSource), Is.All.EqualTo(DetectionFloorSource.ReporterRegionMinimum), "the mzML has no noise table");

            // The snip renumbers its scans by position; the instrument's scan number is in the native ID.
            MsDataScan scan28298 = ms3.Single(s => s.NativeId.EndsWith("scan=28298"));
            ReporterIonReading[] readings = tag.GetReporterIonReadings(scan28298);
            int[] belowFloor = Enumerable.Range(0, 11).Where(i => readings[i].Kind == ReporterIonReadingKind.BelowFloor).ToArray();
            Assert.That(belowFloor.Select(i => tag.ChannelLabels[i]), Is.EqualTo(new[] { "130C", "131C" }));
            Assert.That(readings[belowFloor[0]].Floor, Is.EqualTo(3743.227).Within(0.001));
        }
    }
}
