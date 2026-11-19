using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.SpectraFileReading
{
    /// <summary>
    /// The sampled noise table (<see cref="MsDataScan.NoiseData"/>): read from Thermo files, written to and read from
    /// mzML with the PSI-MS sampled noise arrays, read from the older non-standard arrays this library used to write,
    /// and the noise, baseline and signal-to-noise read off it.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestNoiseData
    {
        private static string DataFile(string name) => Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", name);

        /// <summary>
        /// Values read with Thermo's own RawFileReader (8.0.37, the copy in Readers/Thermo) from sliced_ethcd.raw scan
        /// 1: the first and last points of its noise table (single precision there), and three centroids with the
        /// signal-to-noise Thermo reports for them (CentroidStream.GetLabelPeaks()[i].SignalToNoise).
        /// </summary>
        private const float Scan1FirstMz = 83.10098f, Scan1FirstNoise = 224846.47f, Scan1FirstBaseline = -32542.05f;
        private const float Scan1LastMz = 4467.6797f, Scan1LastNoise = 701565.56f, Scan1LastBaseline = 35161.848f;

        private static readonly (double Mz, double Intensity, double ThermoSignalToNoise)[] Scan1Centroids =
        {
            (504.5805969238281, 6057006.5, 11.761663269613294),
            (815.9088134765625, 3846208.25, 6.546732416183238),
            (1798.9915771484375, 1278683.25, 1.8968680092394143),
        };

        private static MsDataScan ScanWithNoise(double[,] noiseData, int scanNumber = 1) =>
            new MsDataScan(new MzSpectrum(new[] { 126.0, 500.0 }, new[] { 1000.0, 2000.0 }, false), scanNumber, 1, true,
                Polarity.Positive, scanNumber, new MzRange(100, 2000), "FTMS + p NSI Full ms", MZAnalyzerType.Orbitrap,
                3000, null, noiseData, "scan=" + scanNumber);

        #region the table and what is read off it

        [Test]
        public void CreateNoiseDataOrdersThePointsByMzAndKeepsEachPointsValuesTogether()
        {
            double[,] table = MsDataScan.CreateNoiseData(new[] { 300.0, 100.0, 200.0 }, new[] { 30.0, 10.0, 20.0 }, new[] { 3.0, 1.0, 2.0 });

            Assert.That(table.GetLength(0), Is.EqualTo(3));
            Assert.That(Enumerable.Range(0, 3).Select(i => (table[0, i], table[1, i], table[2, i])),
                Is.EqualTo(new[] { (100.0, 10.0, 1.0), (200.0, 20.0, 2.0), (300.0, 30.0, 3.0) }));
        }

        [Test]
        public void CreateNoiseDataRefusesMismatchedArraysAndGivesNullForNoPoints()
        {
            Assert.That(MsDataScan.CreateNoiseData(new double[0], new double[0], new double[0]), Is.Null);
            Assert.Throws<ArgumentException>(() => MsDataScan.CreateNoiseData(new[] { 1.0 }, new[] { 1.0, 2.0 }, new[] { 1.0 }));
            Assert.Throws<ArgumentException>(() => MsDataScan.CreateNoiseData(new[] { 1.0 }, new[] { 1.0 }, new double[0]));
            Assert.Throws<ArgumentNullException>(() => MsDataScan.CreateNoiseData(null, new[] { 1.0 }, new[] { 1.0 }));
        }

        [Test]
        public void NoiseAndBaselineAreInterpolatedBetweenPointsAndHeldBeyondTheEnds()
        {
            var scan = ScanWithNoise(MsDataScan.CreateNoiseData(new[] { 100.0, 200.0, 400.0 }, new[] { 1000.0, 2000.0, 6000.0 }, new[] { -10.0, 30.0, 70.0 }));

            Assert.That(scan.TryGetNoiseAndBaseline(150, out double noise, out double baseline), Is.True);
            Assert.That((noise, baseline), Is.EqualTo((1500.0, 10.0)));
            Assert.That(scan.TryGetNoiseAndBaseline(300, out noise, out baseline), Is.True);
            Assert.That((noise, baseline), Is.EqualTo((4000.0, 50.0)));
            Assert.That(scan.TryGetNoiseAndBaseline(200, out noise, out baseline), Is.True);
            Assert.That((noise, baseline), Is.EqualTo((2000.0, 30.0)));
            Assert.That(scan.TryGetNoiseAndBaseline(50, out noise, out baseline), Is.True);
            Assert.That((noise, baseline), Is.EqualTo((1000.0, -10.0)));
            Assert.That(scan.TryGetNoiseAndBaseline(900, out noise, out baseline), Is.True);
            Assert.That((noise, baseline), Is.EqualTo((6000.0, 70.0)));
        }

        [Test]
        public void ASingleNoisePointHoldsEverywhere()
        {
            var scan = ScanWithNoise(MsDataScan.CreateNoiseData(new[] { 300.0 }, new[] { 50.0 }, new[] { 5.0 }));

            Assert.That(scan.TryGetNoiseAndBaseline(100, out double noise, out double baseline), Is.True);
            Assert.That((noise, baseline), Is.EqualTo((50.0, 5.0)));
            Assert.That(scan.TryGetNoiseAndBaseline(1000, out noise, out baseline), Is.True);
            Assert.That((noise, baseline), Is.EqualTo((50.0, 5.0)));
        }

        [Test]
        public void AScanWithoutAUsableNoiseTableGivesNoNoiseAndNoSignalToNoise()
        {
            // Null, as an ion-trap scan has; and the [2,1] table the Bruker reader builds, which is not three rows.
            foreach (double[,] noiseData in new[] { null, new double[,] { { 0 }, { 0 } } })
            {
                var scan = ScanWithNoise(noiseData);
                Assert.That(scan.TryGetNoiseAndBaseline(300, out double noise, out double baseline), Is.False);
                Assert.That(double.IsNaN(noise) && double.IsNaN(baseline), Is.True);
                Assert.That(scan.GetSignalToNoise(300, 1000), Is.NaN);
            }
        }

        [Test]
        public void SignalToNoiseIsIntensityAboveBaselineOverNoiseAboveBaseline()
        {
            var scan = ScanWithNoise(MsDataScan.CreateNoiseData(new[] { 100.0, 300.0 }, new[] { 1100.0, 1100.0 }, new[] { 100.0, 100.0 }));

            Assert.That(scan.GetSignalToNoise(200, 5100), Is.EqualTo(5.0));
            Assert.That(scan.GetSignalToNoise(200, 100), Is.EqualTo(0.0));
        }

        [Test]
        public void SignalToNoiseIsNaNWhereNoiseDoesNotExceedBaseline()
        {
            var scan = ScanWithNoise(MsDataScan.CreateNoiseData(new[] { 100.0, 300.0 }, new[] { 100.0, 50.0 }, new[] { 100.0, 100.0 }));

            Assert.That(scan.GetSignalToNoise(100, 5000), Is.NaN);
            Assert.That(scan.GetSignalToNoise(300, 5000), Is.NaN);
        }

        #endregion

        #region Thermo

        [Test]
        public void ThermoScansCarryTheNoiseTableThermoRecords()
        {
            var file = MsDataFileReader.GetDataFile(DataFile("sliced_ethcd.raw"));
            file.LoadAllStaticData();
            MsDataScan[] scans = file.GetAllScansList().ToArray();

            Assert.That(scans, Has.Length.EqualTo(6));
            foreach (MsDataScan scan in scans)
            {
                Assert.That(scan.NoiseData, Is.Not.Null, $"scan {scan.OneBasedScanNumber}");
                Assert.That(scan.NoiseData.GetLength(0), Is.EqualTo(3));
                Assert.That(scan.NoiseData.GetLength(1), Is.EqualTo(40), $"scan {scan.OneBasedScanNumber}");
                double[] mzs = Enumerable.Range(0, 40).Select(i => scan.NoiseData[0, i]).ToArray();
                Assert.That(mzs, Is.Ordered.Ascending);
            }

            double[,] table = scans[0].NoiseData;
            Assert.That((table[0, 0], table[1, 0], table[2, 0]), Is.EqualTo(((double)Scan1FirstMz, (double)Scan1FirstNoise, (double)Scan1FirstBaseline)));
            Assert.That((table[0, 39], table[1, 39], table[2, 39]), Is.EqualTo(((double)Scan1LastMz, (double)Scan1LastNoise, (double)Scan1LastBaseline)));
        }

        [Test]
        public void ThermoSignalToNoiseIsReproducedFromTheNoiseTable()
        {
            var file = MsDataFileReader.GetDataFile(DataFile("sliced_ethcd.raw"));
            file.LoadAllStaticData();
            MsDataScan scan = file.GetOneBasedScan(1);

            foreach (var (mz, intensity, thermo) in Scan1Centroids)
            {
                int peak = scan.MassSpectrum.GetClosestPeakIndex(mz);
                Assert.That(scan.MassSpectrum.XArray[peak], Is.EqualTo(mz));
                Assert.That(scan.MassSpectrum.YArray[peak], Is.EqualTo(intensity));
                Assert.That(scan.GetSignalToNoise(mz, intensity), Is.EqualTo(thermo).Within(1e-4).Percent, $"m/z {mz}");
            }
        }

        [Test]
        public void ThermoDynamicReadingGivesTheSameNoiseTable()
        {
            var staticFile = MsDataFileReader.GetDataFile(DataFile("sliced_ethcd.raw"));
            staticFile.LoadAllStaticData();
            var dynamicFile = MsDataFileReader.GetDataFile(DataFile("sliced_ethcd.raw"));
            dynamicFile.InitiateDynamicConnection();
            try
            {
                for (int n = 1; n <= 6; n++)
                {
                    Assert.That(dynamicFile.GetOneBasedScanFromDynamicConnection(n).NoiseData, Is.EqualTo(staticFile.GetOneBasedScan(n).NoiseData), $"scan {n}");
                }
            }
            finally
            {
                dynamicFile.CloseDynamicConnection();
            }
        }

        #endregion

        #region mzML

        /// <summary>The scan's noise table read back, statically and dynamically.</summary>
        private static IEnumerable<(string How, MsDataScan Scan)> ReadBack(string path, int scanNumber)
        {
            var staticFile = MsDataFileReader.GetDataFile(path);
            staticFile.LoadAllStaticData();
            yield return ("static", staticFile.GetOneBasedScan(scanNumber));

            var dynamicFile = MsDataFileReader.GetDataFile(path);
            dynamicFile.InitiateDynamicConnection();
            try
            {
                yield return ("dynamic", dynamicFile.GetOneBasedScanFromDynamicConnection(scanNumber));
            }
            finally
            {
                dynamicFile.CloseDynamicConnection();
            }
        }

        private static double[] DecodeArray(string spectrumXml, string accession)
        {
            Match array = Regex.Match(spectrumXml, "<binaryDataArray [^>]*>(?:(?!</binaryDataArray>).)*?accession=\"" + accession + "\".*?</binaryDataArray>", RegexOptions.Singleline);
            Assert.That(array.Success, Is.True, accession + " array not written");
            byte[] bytes = Convert.FromBase64String(Regex.Match(array.Value, "<binary>(.*?)</binary>", RegexOptions.Singleline).Groups[1].Value);
            return Enumerable.Range(0, bytes.Length / 8).Select(k => BitConverter.ToDouble(bytes, k * 8)).ToArray();
        }

        [Test]
        public void MzmlWritesTheNoiseTableAsPsiMsSampledNoiseArraysAndReadsItBack()
        {
            double[,] noiseData = MsDataScan.CreateNoiseData(new[] { 100.0, 250.0, 900.0 }, new[] { 1000.0, 1500.0, 4000.0 }, new[] { -50.0, 20.0, 300.0 });
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "noise_psi_roundtrip.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(new[] { ScanWithNoise(noiseData) }), path, false);
            try
            {
                string xml = File.ReadAllText(path);
                Assert.That(xml, Does.Not.Contain("kelleherCustomType"));
                Assert.That(xml, Does.Not.Contain("MS:1000786"));
                Assert.That(DecodeArray(xml, "MS:1002743"), Is.EqualTo(new[] { 100.0, 250.0, 900.0 }));
                // PSI-MS: the sampled noise intensity is the noise above the baseline.
                Assert.That(DecodeArray(xml, "MS:1002744"), Is.EqualTo(new[] { 1050.0, 1480.0, 3700.0 }));
                Assert.That(DecodeArray(xml, "MS:1002745"), Is.EqualTo(new[] { -50.0, 20.0, 300.0 }));

                foreach (var (how, scan) in ReadBack(path, 1))
                {
                    Assert.That(scan.NoiseData, Is.EqualTo(noiseData), how);
                }
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void ThermoSignalToNoiseSurvivesAnMzmlRoundTrip()
        {
            // As MetaMorpheus's Calibrate -> Search does: read the raw file, write mzML, search the mzML.
            var raw = MsDataFileReader.GetDataFile(DataFile("sliced_ethcd.raw"));
            raw.LoadAllStaticData();
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "noise_thermo_roundtrip.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(raw, path, false);
            try
            {
                foreach (var (how, scan) in ReadBack(path, 1))
                {
                    Assert.That(scan.NoiseData, Is.EqualTo(raw.GetOneBasedScan(1).NoiseData).Within(1e-9).Percent, how);
                    foreach (var (mz, intensity, thermo) in Scan1Centroids)
                    {
                        Assert.That(scan.GetSignalToNoise(mz, intensity), Is.EqualTo(thermo).Within(1e-4).Percent, $"{how} m/z {mz}");
                    }
                }
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void MzmlWrittenBeforeThePsiMsTermsIsReadByItsOwnLabels()
        {
            // SmallCalibratibleYeast.mzml carries the non-standard arrays this library wrote before: MS:1000786 with value
            // "mass", "SignalToNoise" (the noise level itself) and "baseline". Scan 1's first and last points, decoded
            // from the file.
            foreach (var (how, scan) in ReadBack(DataFile("SmallCalibratibleYeast.mzml"), 1))
            {
                double[,] table = scan.NoiseData;
                Assert.That(table, Is.Not.Null, how);
                Assert.That((table.GetLength(0), table.GetLength(1)), Is.EqualTo((3, 10)), how);
                Assert.That((table[0, 0], table[1, 0], table[2, 0]), Is.EqualTo((311.45068359375, 27591.466796875, 188.302734375)), how);
                Assert.That((table[0, 9], table[1, 9], table[2, 9]), Is.EqualTo((1256.1822509765625, 35649.3046875, 4908.7109375)), how);
            }
        }

        [Test]
        public void MzmlWithoutNoiseArraysReadsWithoutANoiseTable()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "noise_none_roundtrip.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(new[] { ScanWithNoise(null) }), path, false);
            try
            {
                Assert.That(File.ReadAllText(path), Does.Not.Contain("MS:100274"));
                foreach (var (how, scan) in ReadBack(path, 1))
                {
                    Assert.That(scan.NoiseData, Is.Null, how);
                }
            }
            finally
            {
                File.Delete(path);
            }
        }

        #endregion
    }
}
