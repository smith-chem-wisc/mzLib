using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.FlashLFQ
{
    /// <summary>
    /// IndexPeaks now counts the peaks in each bin before filling it, so each bin's list is allocated once at its final
    /// size instead of growing by doubling. The index must be exactly what the previous code built: the same bins, the
    /// same peaks in the same order, and the same scan table.
    /// </summary>
    [TestFixture]
    public static class TestIndexPeaksExactCapacity
    {
        [Test]
        public static void IndexMatchesThePreviousCodeOnARealFile()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData", "20100614_Velos1_TaGe_SA_K562_3.mzML");
            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData();
            MsDataScan[] scans = reader.GetMS1Scans().Where(s => s != null && s.MsnOrder == 1).OrderBy(s => s.OneBasedScanNumber).ToArray();

            AssertSameIndex(scans, expectedPeaks: 90_000);
        }

        [Test]
        public static void IndexMatchesThePreviousCodeOnDenseRandomScans()
        {
            var random = new Random(23);
            for (int trial = 0; trial < 20; trial++)
            {
                var scans = new MsDataScan[random.Next(1, 30)];
                for (int s = 0; s < scans.Length; s++)
                {
                    // Near-duplicate m/z values put many peaks into the same bins, including bins at both rounding edges
                    double[] mz = Enumerable.Range(0, random.Next(1, 2000))
                        .Select(_ => Math.Round(100 + random.NextDouble() * 50, random.Next(2, 6))).Order().ToArray();
                    double[] intensity = mz.Select(_ => random.NextDouble() * 1e6).ToArray();
                    scans[s] = new MsDataScan(new MzSpectrum(mz, intensity, false), s + 1, 1, true, Polarity.Positive, 1.0 + s / 10.0,
                        new MzRange(100, 200), "f", MZAnalyzerType.Orbitrap, intensity.Sum(), 1.0, null, "scan=" + (s + 1));
                }

                AssertSameIndex(scans, expectedPeaks: 1);
            }
        }

        private static void AssertSameIndex(MsDataScan[] scans, int expectedPeaks)
        {
            var previous = new IndexProbe();
            var current = new IndexProbe();
            Assert.That(current.IndexPeaks(scans), Is.EqualTo(previous.PreviousIndexPeaks(scans)));

            Assert.That(current.Bins.Length, Is.EqualTo(previous.Bins.Length));
            int peaks = 0;
            for (int bin = 0; bin < previous.Bins.Length; bin++)
            {
                var expected = previous.Bins[bin];
                var actual = current.Bins[bin];
                if (expected == null)
                {
                    Assert.That(actual, Is.Null, $"bin {bin}");
                    continue;
                }

                Assert.That(actual.Count, Is.EqualTo(expected.Count), $"bin {bin}");
                Assert.That(actual.Capacity, Is.EqualTo(actual.Count), $"bin {bin} should be allocated at its final size");
                for (int i = 0; i < expected.Count; i++)
                {
                    Assert.That(actual[i].ZeroBasedScanIndex, Is.EqualTo(expected[i].ZeroBasedScanIndex));
                    Assert.That(BitConverter.SingleToInt32Bits(actual[i].Mz), Is.EqualTo(BitConverter.SingleToInt32Bits(expected[i].Mz)));
                    Assert.That(BitConverter.SingleToInt32Bits(actual[i].Intensity), Is.EqualTo(BitConverter.SingleToInt32Bits(expected[i].Intensity)));
                    Assert.That(BitConverter.SingleToInt32Bits(actual[i].RetentionTime), Is.EqualTo(BitConverter.SingleToInt32Bits(expected[i].RetentionTime)));
                }
                peaks += expected.Count;
            }
            Assert.That(peaks, Is.GreaterThanOrEqualTo(expectedPeaks));

            Assert.That(current.ScanInfoArray.Length, Is.EqualTo(previous.ScanInfoArray.Length));
            for (int s = 0; s < previous.ScanInfoArray.Length; s++)
            {
                Assert.That(current.ScanInfoArray[s].OneBasedScanNumber, Is.EqualTo(previous.ScanInfoArray[s].OneBasedScanNumber));
                Assert.That(current.ScanInfoArray[s].ZeroBasedScanIndex, Is.EqualTo(previous.ScanInfoArray[s].ZeroBasedScanIndex));
                Assert.That(current.ScanInfoArray[s].RetentionTime, Is.EqualTo(previous.ScanInfoArray[s].RetentionTime));
            }
        }

        private class IndexProbe : IndexingEngine<IndexedMassSpectralPeak>
        {
            public List<IndexedMassSpectralPeak>[] Bins => IndexedPeaks;

            /// <summary>
            /// Verbatim copy of IndexingEngine.IndexPeaks before it counted each bin first
            /// </summary>
            public bool PreviousIndexPeaks(MsDataScan[] scanArray)
            {
                if (scanArray.IsNullOrEmpty() || scanArray.All(p => p == null))
                    return false;

                IndexedPeaks = new List<IndexedMassSpectralPeak>[(int)Math.Ceiling(scanArray.Where(p => p != null
                    && p.MassSpectrum.LastX != null).Max(p => p.MassSpectrum.LastX.Value) * BinsPerDalton) + 1];
                ScanInfoArray = new ScanInfo[scanArray.Length];

                for (int scanIndex = 0; scanIndex < scanArray.Length; scanIndex++)
                {
                    ScanInfoArray[scanIndex] = new ScanInfo(scanArray[scanIndex].OneBasedScanNumber, scanIndex, scanArray[scanIndex].RetentionTime, scanArray[scanIndex].MsnOrder);

                    for (int j = 0; j < scanArray[scanIndex].MassSpectrum.XArray.Length; j++)
                    {
                        int roundedMz = (int)Math.Round(scanArray[scanIndex].MassSpectrum.XArray[j] * BinsPerDalton, 0);
                        IndexedPeaks[roundedMz] ??= new List<IndexedMassSpectralPeak>();
                        IndexedPeaks[roundedMz].Add((IndexedMassSpectralPeak)(IIndexedPeak)
                            new IndexedMassSpectralPeak(
                                scanArray[scanIndex].MassSpectrum.XArray[j],
                                scanArray[scanIndex].MassSpectrum.YArray[j],
                                scanIndex,
                                scanArray[scanIndex].RetentionTime));
                    }
                }
                if (IndexedPeaks == null || IndexedPeaks.Length == 0)
                    return false;
                else
                    return true;
            }
        }
    }
}
