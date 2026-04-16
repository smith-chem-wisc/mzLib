using System;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test.FileReadingTests
{
    [TestFixture]
    internal class ImspExportServiceTests
    {
        [Test]
        public static void WriteBytesProducesExpectedTinyImsp()
        {
            var scans = MakeTinyScans();
            byte[] bytes = new ImspExportService().WriteBytes(scans, intensityThreshold: 0);

            using (var reader = new BinaryReader(new MemoryStream(bytes)))
            {
                Assert.That(reader.ReadBytes(4), Is.EqualTo(new[] { (byte)'I', (byte)'M', (byte)'S', (byte)'P' }));
                Assert.That(reader.ReadUInt32(), Is.EqualTo(1));
                Assert.That(reader.ReadUInt32(), Is.EqualTo(100));
                Assert.That(reader.ReadUInt32(), Is.EqualTo(6));
                Assert.That(reader.ReadUInt32(), Is.EqualTo(9));
                Assert.That(reader.ReadUInt32(), Is.EqualTo(3));

                AssertScan(reader, 1, 1.0, 225000);
                AssertScan(reader, 2, 2.0, 250000);
                AssertScan(reader, 3, 3.0, 210000);

                AssertBin(reader, 50012, 0, 2);
                AssertBin(reader, 60023, 2, 1);
                AssertBin(reader, 75057, 3, 2);
                AssertBin(reader, 85068, 5, 1);
                AssertBin(reader, 100090, 6, 2);
                AssertBin(reader, 125035, 8, 1);

                AssertPeak(reader, 5001234, 50000, 0);
                AssertPeak(reader, 5001234, 60000, 1);
                AssertPeak(reader, 6002345, 45000, 2);
                AssertPeak(reader, 7505678, 100000, 0);
                AssertPeak(reader, 7505678, 110000, 1);
                AssertPeak(reader, 8506789, 95000, 2);
                AssertPeak(reader, 10009012, 75000, 0);
                AssertPeak(reader, 10009012, 70000, 2);
                AssertPeak(reader, 12503456, 80000, 1);

                Assert.That(reader.BaseStream.Position, Is.EqualTo(reader.BaseStream.Length));
            }
        }

        [Test]
        public static void WriteFileReturnsOutputPathCompatibleBytes()
        {
            string outputPath = Path.Combine(TestContext.CurrentContext.WorkDirectory, "tiny-known.imsp");
            var service = new ImspExportService();

            int peakCount = service.WriteFile(MakeTinyScans(), outputPath, intensityThreshold: 0);

            Assert.That(peakCount, Is.EqualTo(9));
            Assert.That(File.Exists(outputPath), Is.True);
            Assert.That(File.ReadAllBytes(outputPath), Is.EqualTo(service.WriteBytes(MakeTinyScans(), intensityThreshold: 0)));
        }

        private static MsDataScan[] MakeTinyScans()
        {
            return new[]
            {
                MakeScan(1, 1.0, new[] { (500.1234, 50000.0), (750.5678, 100000.0), (1000.9012, 75000.0) }),
                MakeScan(2, 2.0, new[] { (500.1234, 60000.0), (750.5678, 110000.0), (1250.3456, 80000.0) }),
                MakeScan(3, 3.0, new[] { (600.2345, 45000.0), (850.6789, 95000.0), (1000.9012, 70000.0) }),
            };
        }

        private static MsDataScan MakeScan(int oneBasedScanNumber, double retentionTime, (double mz, double intensity)[] peaks)
        {
            double[] mzArray = peaks.Select(p => p.mz).ToArray();
            double[] intensityArray = peaks.Select(p => p.intensity).ToArray();

            return new MsDataScan(
                massSpectrum: new MzSpectrum(mzArray, intensityArray, false),
                oneBasedScanNumber: oneBasedScanNumber,
                msnOrder: 1,
                isCentroid: true,
                polarity: Polarity.Positive,
                retentionTime: retentionTime,
                scanWindowRange: new MzRange(mzArray.Min(), mzArray.Max()),
                scanFilter: "FTMS + p NSI Full ms",
                mzAnalyzer: MZAnalyzerType.Orbitrap,
                totalIonCurrent: intensityArray.Sum(),
                injectionTime: 1.0,
                noiseData: null,
                nativeId: "scan=" + oneBasedScanNumber);
        }

        private static void AssertScan(BinaryReader reader, uint scanNumber, double retentionTime, float tic)
        {
            Assert.That(reader.ReadUInt32(), Is.EqualTo(scanNumber));
            Assert.That(reader.ReadDouble(), Is.EqualTo(retentionTime));
            Assert.That(reader.ReadSingle(), Is.EqualTo(tic));
        }

        private static void AssertBin(BinaryReader reader, uint binIndex, uint peakOffset, uint peakCount)
        {
            Assert.That(reader.ReadUInt32(), Is.EqualTo(binIndex));
            Assert.That(reader.ReadUInt32(), Is.EqualTo(peakOffset));
            Assert.That(reader.ReadUInt32(), Is.EqualTo(peakCount));
        }

        private static void AssertPeak(BinaryReader reader, uint mzTenThousandths, float intensity, uint scanIndex)
        {
            Assert.That(reader.ReadUInt32(), Is.EqualTo(mzTenThousandths));
            Assert.That(reader.ReadSingle(), Is.EqualTo(intensity));
            Assert.That(reader.ReadUInt32(), Is.EqualTo(scanIndex));
        }
    }
}
