using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using SpectralAveraging;
using SpectralAveragingExtensions;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestAveragingSpectraFileHandling
    {
        public static List<MsDataScan> DummyAllMs1Scans { get; set; }
        public static List<MsDataScan> DummyDDAScansInOrder { get; set; }
        public static List<MsDataScan> DummyDDAScansOutOfOrder { get; set; }
        public static List<MsDataScan> ActualScans { get; set; }
        public static MzLibSpectralAveragingOptions MzLibSpectralAveragingOptions { get; set; }

        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            ActualScans = SpectraFileHandler.LoadAllScansFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"AveragingTestData\TDYeastFractionMS1.mzML"));
            MzLibSpectralAveragingOptions = new();
            MzLibSpectralAveragingOptions.SpectralAveragingOptions.RejectionType = RejectionType.NoRejection;
            MzLibSpectralAveragingOptions.SpectralAveragingOptions.WeightingType = WeightingType.NoWeight;
            MzLibSpectralAveragingOptions.SpectralAveragingOptions.PerformNormalization = false;

            double[] xArray = new double[] { 100.1453781, 200, 300, 400, 500, 600, 700, 800, 900.4123745 };
            double[] yArray1 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0, 0 };
            double[] yArray2 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0, 0 };
            double[] yArray3 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0, 0 };
            double[] yArray4 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0, 0 };
            double[] yArray5 = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            List<MzSpectrum> mzSpectra = new();
            mzSpectra.Add(new(xArray, yArray1, true));
            mzSpectra.Add(new(xArray, yArray2, true));
            mzSpectra.Add(new(xArray, yArray3, true));
            mzSpectra.Add(new(xArray, yArray4, true));
            mzSpectra.Add(new(xArray, yArray5, true));
            mzSpectra.Add(new(xArray, yArray1, true));
            mzSpectra.Add(new(xArray, yArray2, true));
            mzSpectra.Add(new(xArray, yArray3, true));
            mzSpectra.Add(new(xArray, yArray4, true));
            mzSpectra.Add(new(xArray, yArray5, true));

            DummyAllMs1Scans = new();
            for (int i = 0; i < mzSpectra.Count; i++)
            {
                MsDataScan dummyScan = new(mzSpectra[i], i + 1, 1, true, Polarity.Positive, i,
                    mzSpectra[i].Range, null, MZAnalyzerType.Orbitrap, 10,
                    1, null, ActualScans.First().NativeId);
                DummyAllMs1Scans.Add(dummyScan);
            }

            // DDA creation
            DummyDDAScansInOrder = new();
            DummyDDAScansOutOfOrder = new();
            List<MsDataScan> dummyDDAMs1Scans = new();
            List<MsDataScan> dummyDDAMs2Scans = new();
            List<MzSpectrum> spectraForDummyDDAMs1Scans = mzSpectra.Take(5).ToList();
            List<MzSpectrum> spectraForDummyDDAMs2Scans = mzSpectra;
            spectraForDummyDDAMs2Scans.AddRange(mzSpectra.Take(5));

            // create Ms1's
            for (int i = 0; i < spectraForDummyDDAMs1Scans.Count; i++)
            {
                MsDataScan dummyScan = new(mzSpectra[i], i + 1, 1, true, Polarity.Positive, i,
                    mzSpectra[i].Range, null, MZAnalyzerType.Orbitrap, 10,
                    1, null, ActualScans.First().NativeId);
                dummyDDAMs1Scans.Add(dummyScan);
            }

            // create Ms2's
            for (int i = 0; i < spectraForDummyDDAMs2Scans.Count; i++)
            {
                MsDataScan dummyScan = new(mzSpectra[i], i + 1, 2, true, Polarity.Positive, i,
                    mzSpectra[i].Range, null, MZAnalyzerType.Orbitrap, 10,
                    1, null, ActualScans.First().NativeId);
                dummyDDAMs2Scans.Add(dummyScan);
            }

            // create in order dda
            int scanIndex = 1;
            for (int i = 0; i < dummyDDAMs1Scans.Count; i++)
            {
                dummyDDAMs1Scans[i].SetOneBasedScanNumber(scanIndex);
                DummyDDAScansInOrder.Add(dummyDDAMs1Scans[i]);
                int precursorIndex = dummyDDAMs1Scans[i].OneBasedScanNumber;
                scanIndex++;
                // add 3 ms2's to the ms1
                for (int j = 0; j < 3; j++)
                {
                    dummyDDAMs2Scans[3 * i + j].SetOneBasedScanNumber(scanIndex);
                    dummyDDAMs2Scans[3 * i + j].SetOneBasedPrecursorScanNumber(precursorIndex);
                    DummyDDAScansInOrder.Add(dummyDDAMs2Scans[3 * i + j]);
                    scanIndex++;
                }
            }

            // changed the dda scans so their retention times are in the correct order, could be a more elegent way to do it, but I am tired and lazy
            int retentionTimeIndex = 1;
            List<MsDataScan> formattedScans = new List<MsDataScan>();
            foreach (var dataScan in DummyDDAScansInOrder.OrderBy(p => p.OneBasedScanNumber))
            {
                MsDataScan newScan = new(dataScan.MassSpectrum, dataScan.OneBasedScanNumber, dataScan.MsnOrder,
                    dataScan.IsCentroid, dataScan.Polarity, retentionTimeIndex,
                    dataScan.ScanWindowRange, dataScan.ScanFilter, dataScan.MzAnalyzer, dataScan.TotalIonCurrent,
                    dataScan.InjectionTime, dataScan.NoiseData, dataScan.NativeId,
                    oneBasedPrecursorScanNumber: dataScan.OneBasedPrecursorScanNumber);
                retentionTimeIndex++;
                formattedScans.Add(newScan);
            }
            DummyDDAScansInOrder = formattedScans;

            // create out of order dda
            DummyDDAScansOutOfOrder.AddRange(DummyDDAScansInOrder.Take(5));
            Stack<MsDataScan> dummyMs1ScanStack = new Stack<MsDataScan>(dummyDDAMs1Scans.GetRange(2, 3));
            Stack<MsDataScan> dummyMs2ScanStack = new Stack<MsDataScan>(dummyDDAMs2Scans.GetRange(3, 12));
            MsDataScan scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(6);
            scan.SetOneBasedPrecursorScanNumber(5);
            DummyDDAScansOutOfOrder.Add(scan);
            scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(7);
            scan.SetOneBasedPrecursorScanNumber(5);
            DummyDDAScansOutOfOrder.Add(scan);

            scan = dummyMs1ScanStack.Pop();
            scan.SetOneBasedScanNumber(8);
            DummyDDAScansOutOfOrder.Add(scan);
            scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(9);
            scan.SetOneBasedPrecursorScanNumber(5);
            DummyDDAScansOutOfOrder.Add(scan);
            scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(10);
            scan.SetOneBasedPrecursorScanNumber(8);
            DummyDDAScansOutOfOrder.Add(scan);
            scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(11);
            scan.SetOneBasedPrecursorScanNumber(8);
            DummyDDAScansOutOfOrder.Add(scan);
            scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(12);
            scan.SetOneBasedPrecursorScanNumber(8);
            DummyDDAScansOutOfOrder.Add(scan);

            scan = dummyMs1ScanStack.Pop();
            scan.SetOneBasedScanNumber(13);
            DummyDDAScansOutOfOrder.Add(scan);
            scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(14);
            scan.SetOneBasedPrecursorScanNumber(13);
            DummyDDAScansOutOfOrder.Add(scan);

            scan = dummyMs1ScanStack.Pop();
            scan.SetOneBasedScanNumber(15);
            DummyDDAScansOutOfOrder.Add(scan);
            scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(16);
            scan.SetOneBasedPrecursorScanNumber(15);
            DummyDDAScansOutOfOrder.Add(scan);
            scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(17);
            scan.SetOneBasedPrecursorScanNumber(15);
            DummyDDAScansOutOfOrder.Add(scan);
            scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(18);
            scan.SetOneBasedPrecursorScanNumber(15);
            DummyDDAScansOutOfOrder.Add(scan);
            scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(19);
            scan.SetOneBasedPrecursorScanNumber(13);
            DummyDDAScansOutOfOrder.Add(scan);
            scan = dummyMs2ScanStack.Pop();
            scan.SetOneBasedScanNumber(20);
            scan.SetOneBasedPrecursorScanNumber(15);
            DummyDDAScansOutOfOrder.Add(scan);

            formattedScans = new();
            retentionTimeIndex = 1;
            foreach (var dataScan in DummyDDAScansOutOfOrder.OrderBy(p => p.OneBasedScanNumber))
            {
                MsDataScan newScan = new(dataScan.MassSpectrum, dataScan.OneBasedScanNumber, dataScan.MsnOrder,
                    dataScan.IsCentroid, dataScan.Polarity, retentionTimeIndex,
                    dataScan.ScanWindowRange, dataScan.ScanFilter, dataScan.MzAnalyzer, dataScan.TotalIonCurrent,
                    dataScan.InjectionTime, dataScan.NoiseData, dataScan.NativeId,
                    oneBasedPrecursorScanNumber: dataScan.OneBasedPrecursorScanNumber);
                retentionTimeIndex++;
                formattedScans.Add(newScan);
            }
            DummyDDAScansOutOfOrder = formattedScans;
        }

        [Test]
        public static void TestAverageAll()
        {
            MzLibSpectralAveragingOptions.SpectraFileProcessingType = SpectraFileProcessingType.AverageAll;
            MsDataScan[] averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            double[] expected = new double[] { 0, 4, 0, 0, 0, 0, 0, 8, 0 };
            Assert.That(averagedScans.Length == 1);
            Assert.That(expected.SequenceEqual(averagedScans.First().MassSpectrum.YArray));

            averagedScans = SpectraFileProcessing.ProcessSpectra(ActualScans, MzLibSpectralAveragingOptions);
            Assert.That(averagedScans.Length == 1);
        }

        [Test]
        public static void TestAverageEveryNScansWithoutOverlap()
        {
            MzLibSpectralAveragingOptions.SpectraFileProcessingType = SpectraFileProcessingType.AverageEverynScans;
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 5;
            MsDataScan[] averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            double[] expected = new double[] { 0, 4, 0, 0, 0, 0, 0, 8, 0 };
            Assert.That(averagedScans.Length == 2);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            Assert.That(averagedScans[1].MassSpectrum.YArray.SequenceEqual(expected));

            MzLibSpectralAveragingOptions.ScanOverlap = 4;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            expected = new double[] { 0, 4, 0, 0, 0, 0, 0, 8, 0 };
            Assert.That(averagedScans.Length == 2);
            Assert.That(MzLibSpectralAveragingOptions.ScanOverlap == 0);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            Assert.That(averagedScans[1].MassSpectrum.YArray.SequenceEqual(expected));

            averagedScans = SpectraFileProcessing.ProcessSpectra(ActualScans, MzLibSpectralAveragingOptions);
            Assert.That(averagedScans.Length == 91);
        }

        [Test]
        public static void TestAverageEveryNScansWithOverlap()
        {
            MzLibSpectralAveragingOptions.SpectraFileProcessingType =
                SpectraFileProcessingType.AverageEverynScansWithOverlap;

            MzLibSpectralAveragingOptions.ScanOverlap = 1;
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 2;
            MsDataScan[] averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            double[] expected = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
            Assert.That(averagedScans.Length == 9);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 3;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            Assert.That(averagedScans.Length == 4);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 4;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            Assert.That(averagedScans.Length == 3);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 5;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            expected = new double[] { 0, 4, 0, 0, 0, 0, 0, 8, 0 };
            Assert.That(averagedScans.Length == 2);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));


            MzLibSpectralAveragingOptions.ScanOverlap = 2;
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 3;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            expected = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
            Assert.That(averagedScans.Length == 8);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 4;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            Assert.That(averagedScans.Length == 4);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 5;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            expected = new double[] { 0, 4, 0, 0, 0, 0, 0, 8, 0 };
            Assert.That(averagedScans.Length == 2);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));

            MzLibSpectralAveragingOptions.ScanOverlap = 3;
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 4;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            expected = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
            Assert.That(averagedScans.Length == 7);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 5;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyAllMs1Scans, MzLibSpectralAveragingOptions);
            expected = new double[] { 0, 4, 0, 0, 0, 0, 0, 8, 0 };
            Assert.That(averagedScans.Length == 3);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
        }

        [Test]
        public static void TestAverageDDAScansWithoutOverlapInOrder()
        {
            MzLibSpectralAveragingOptions.SpectraFileProcessingType = SpectraFileProcessingType.AverageDDAScans;
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 2;
            MsDataScan[] averagedScans = SpectraFileProcessing.ProcessSpectra(DummyDDAScansInOrder, MzLibSpectralAveragingOptions);
            double[] expected = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 2);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 12);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            int?[] expectedNullable = new int?[] { null, 1, 1, 1, 1, 1, 1, null, 8, 8, 8, 8, 8, 8 };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));
        }

        [Test]
        public static void TestAverageDDAScansWithoutOverlapOutOfOrder()
        {
            MzLibSpectralAveragingOptions.SpectraFileProcessingType = SpectraFileProcessingType.AverageDDAScans;
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 2;
            MsDataScan[] averagedScans = SpectraFileProcessing.ProcessSpectra(DummyDDAScansOutOfOrder, MzLibSpectralAveragingOptions);
            double[] expected = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 2);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 11);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            int?[] expectedNullable = new int?[] { null, 1, 1, 1, 1, 1, 1, null, 8, 8, 8, 8, 8 };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));
        }

        [Test]
        public static void TestAverageDDAScansWithOverlapInOrder()
        {
            MzLibSpectralAveragingOptions.SpectraFileProcessingType = SpectraFileProcessingType.AverageDDAScansWithOverlap;
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 2;
            MzLibSpectralAveragingOptions.ScanOverlap = 1;
            MsDataScan[] averagedScans = SpectraFileProcessing.ProcessSpectra(DummyDDAScansInOrder, MzLibSpectralAveragingOptions);
            double[] expected = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 4);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 15);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            expected = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            Assert.AreEqual(expected, averagedScans.Last(p => p.MsnOrder == 1).MassSpectrum.YArray);
            int?[] expectedNullable = new int?[]
            {
                null, 1, 1, 1, null, 5, 5, 5, null, 9, 9, 9, null, 13, 13, 13, 13, 13, 13
            };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));

            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 3;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyDDAScansInOrder, MzLibSpectralAveragingOptions);
            expected = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 2);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 15);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            expected = new double[] { 0, 5.0 * (2.0 / 3.0), 0, 0, 0, 0, 0, 10.0 * (2.0 / 3.0), 0 };
            Assert.AreEqual(expected.Select(p => Math.Round(p, 4)), averagedScans.Last(p => p.MsnOrder == 1).MassSpectrum.YArray.Select(p => Math.Round(p, 4)));
            expectedNullable = new int?[]
            {
                null, 1, 1, 1, 1, 1, 1, null, 8, 8, 8, 8, 8, 8, 8, 8, 8
            };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));

            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 3;
            MzLibSpectralAveragingOptions.ScanOverlap = 2;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyDDAScansInOrder, MzLibSpectralAveragingOptions);
            expected = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 3);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 15);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            expected = new double[] { 0, 5.0 * (2.0 / 3.0), 0, 0, 0, 0, 0, 10.0 * (2.0 / 3.0), 0 };
            Assert.AreEqual(expected.Select(p => Math.Round(p, 4)), averagedScans.Last(p => p.MsnOrder == 1).MassSpectrum.YArray.Select(p => Math.Round(p, 4)));
            expectedNullable = new int?[]
            {
                null, 1, 1, 1, null, 5, 5, 5, null, 9, 9, 9, 9, 9, 9, 9, 9, 9
            };
            Assert.AreEqual(expectedNullable, averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray());
        }
        [Test]
        public static void TestAverageDDAScansWithOverlapOutOfOrder()
        {
            MzLibSpectralAveragingOptions.SpectraFileProcessingType = SpectraFileProcessingType.AverageDDAScansWithOverlap;
            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 2;
            MzLibSpectralAveragingOptions.ScanOverlap = 1;
            MsDataScan[] averagedScans = SpectraFileProcessing.ProcessSpectra(DummyDDAScansOutOfOrder, MzLibSpectralAveragingOptions);
            double[] expected = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 4);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 15);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            Assert.AreEqual(expected, averagedScans.Last(p => p.MsnOrder == 1).MassSpectrum.YArray);
            int?[] expectedNullable = new int?[]
            {
                null, 1, 1, 1, null, 5, 5, 5, null, 9, 9, 9, null, 13, 13, 13, 13, 13, 13
            };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));

            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 3;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyDDAScansOutOfOrder, MzLibSpectralAveragingOptions);
            expected = new double[] { 0, 5.0 * (2.0 / 3.0), 0, 0, 0, 0, 0, 10.0 * (2.0 / 3.0), 0 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 2);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 15);
            Assert.AreEqual(expected.Select(p => Math.Round(p, 4)), averagedScans[0].MassSpectrum.YArray.Select(p => Math.Round(p, 4)));
            Assert.AreEqual(expected.Select(p => Math.Round(p, 4)), averagedScans.Last(p => p.MsnOrder == 1).MassSpectrum.YArray.Select(p => Math.Round(p, 4)));
            expectedNullable = new int?[]
            {
                null, 1, 1, 1, 1, 1, 1, null, 8, 8, 8, 8, 8, 8, 8, 8, 8
            };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));

            MzLibSpectralAveragingOptions.NumberOfScansToAverage = 3;
            MzLibSpectralAveragingOptions.ScanOverlap = 2;
            averagedScans = SpectraFileProcessing.ProcessSpectra(DummyDDAScansOutOfOrder, MzLibSpectralAveragingOptions);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 3);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 15);
            Assert.AreEqual(expected.Select(p => Math.Round(p, 4)), averagedScans[0].MassSpectrum.YArray.Select(p => Math.Round(p, 4)));
            Assert.AreEqual(expected.Select(p => Math.Round(p, 4)), averagedScans.Last(p => p.MsnOrder == 1).MassSpectrum.YArray.Select(p => Math.Round(p, 4)));
            expectedNullable = new int?[]
            {
                null, 1, 1, 1, null, 5, 5, 5, null, 9, 9, 9, 9, 9, 9, 9, 9, 9
            };
            Assert.AreEqual(expectedNullable, averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray());
        }

        [Test]
        public static void TestLoadingRawFilesAndSourceFiles()
        {
            string spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles/small.raw");
            List<MsDataScan> scans = SpectraFileHandler.LoadAllScansFromFile(spectraPath);
            Assert.That(scans.Count == 48);

            Readers.SourceFile file = SpectraFileHandler.GetSourceFile(spectraPath);
            Assert.That(file.NativeIdFormat == "Thermo nativeID format");

            string badPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles/small.toml");
            try
            {
                SpectraFileHandler.LoadAllScansFromFile(badPath);
            }
            catch (MzLibException e)
            {
                Assert.That(e.Message == "File extension not supported.");
            }
            catch (Exception e)
            {
                Assert.That(false);
            }

            try
            {
                SpectraFileHandler.GetSourceFile(badPath);
            }
            catch (MzLibException e)
            {
                Assert.That(e.Message == "File extension not supported.");
            }
            catch (Exception e)
            {
                Assert.That(false);
            }
        }
    }
}
