using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using Readers;
using SpectralAveraging;

namespace Test.AveragingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSpectraFileAveraging
    {
        #region TestCase Setup

        private static double[] xArray = new double[] { 100.1453781, 200, 300, 400, 500, 600, 700, 800, 900.4123745 };
        private static double[] yArray1 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
        private static double[] yArray2 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
        private static double[] yArray3 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
        private static double[] yArray4 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0 };
        private static double[] yArray5 = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        private static List<MzSpectrum> mzSpectra = new()
        {
            new(xArray, yArray1, true),
            new(xArray, yArray2, true),
            new(xArray, yArray3, true),
            new(xArray, yArray4, true),
            new(xArray, yArray5, true),
            new(xArray, yArray1, true),
            new(xArray, yArray2, true),
            new(xArray, yArray3, true),
            new(xArray, yArray4, true),
            new(xArray, yArray5, true)
        };

        public static List<MsDataScan> DummyAllMs1Scans
        {
            get
            {
                var dummyScans = new List<MsDataScan>();

                MzSpectrum[] mzSpectraArr = new MzSpectrum[10];
                TestSpectraFileAveraging.mzSpectra.CopyTo(mzSpectraArr);
                var mzSpectra = mzSpectraArr.ToList();

                for (int i = 0; i < mzSpectra.Count; i++)
                {
                    MsDataScan dummyScan = new(mzSpectra[i], i + 1, 1, true, Polarity.Positive, i,
                        mzSpectra[i].Range, null, MZAnalyzerType.Orbitrap, 10,
                        1, null, NativeId);
                    dummyScans.Add(dummyScan);
                }
                return dummyScans;
            }
        }

        public static List<MsDataScan> DummyDDAScansInOrder
        {
            get
            {
                MzSpectrum[] mzSpectraArr = new MzSpectrum[10];
                TestSpectraFileAveraging.mzSpectra.CopyTo(mzSpectraArr);
                var mzSpectra = mzSpectraArr.ToList();

                var dummyDDAScansInOrder = new List<MsDataScan>();
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
                        1, null, NativeId);
                    dummyDDAMs1Scans.Add(dummyScan);
                }

                // create Ms2's
                for (int i = 0; i < spectraForDummyDDAMs2Scans.Count; i++)
                {
                    MsDataScan dummyScan = new(mzSpectra[i], i + 1, 2, true, Polarity.Positive, i,
                        mzSpectra[i].Range, null, MZAnalyzerType.Orbitrap, 10,
                        1, null, NativeId);
                    dummyDDAMs2Scans.Add(dummyScan);
                }

                // create in order dda
                int scanIndex = 1;
                for (int i = 0; i < dummyDDAMs1Scans.Count; i++)
                {
                    dummyDDAMs1Scans[i].SetOneBasedScanNumber(scanIndex);
                    dummyDDAScansInOrder.Add(dummyDDAMs1Scans[i]);
                    int precursorIndex = dummyDDAMs1Scans[i].OneBasedScanNumber;
                    scanIndex++;
                    // add 3 ms2's to the ms1
                    for (int j = 0; j < 3; j++)
                    {
                        dummyDDAMs2Scans[3 * i + j].SetOneBasedScanNumber(scanIndex);
                        dummyDDAMs2Scans[3 * i + j].SetOneBasedPrecursorScanNumber(precursorIndex);
                        dummyDDAScansInOrder.Add(dummyDDAMs2Scans[3 * i + j]);
                        scanIndex++;
                    }
                }

                // changed the dda scans so their retention times are in the correct order, could be a more elegent way to do it, but I am tired and lazy
                int retentionTimeIndex = 1;
                List<MsDataScan> formattedScans = new List<MsDataScan>();
                foreach (var dataScan in dummyDDAScansInOrder.OrderBy(p => p.OneBasedScanNumber))
                {
                    MsDataScan newScan = new(dataScan.MassSpectrum, dataScan.OneBasedScanNumber, dataScan.MsnOrder,
                        dataScan.IsCentroid, dataScan.Polarity, retentionTimeIndex,
                        dataScan.ScanWindowRange, dataScan.ScanFilter, dataScan.MzAnalyzer, dataScan.TotalIonCurrent,
                        dataScan.InjectionTime, dataScan.NoiseData, dataScan.NativeId,
                        oneBasedPrecursorScanNumber: dataScan.OneBasedPrecursorScanNumber);
                    retentionTimeIndex++;
                    formattedScans.Add(newScan);
                }
                return formattedScans;
            }
        }

        public static List<MsDataScan> DummyDDAScansOutOfOrder
        {
            get
            {
                MzSpectrum[] mzSpectraArr = new MzSpectrum[10];
                TestSpectraFileAveraging.mzSpectra.CopyTo(mzSpectraArr);
                var mzSpectra = mzSpectraArr.ToList();

                // DDA creation
                var dummyDDAScansOutOfOrder = new List<MsDataScan> ();
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
                        1, null, NativeId);
                    dummyDDAMs1Scans.Add(dummyScan);
                }

                // create Ms2's
                for (int i = 0; i < spectraForDummyDDAMs2Scans.Count; i++)
                {
                    MsDataScan dummyScan = new(mzSpectra[i], i + 1, 2, true, Polarity.Positive, i,
                        mzSpectra[i].Range, null, MZAnalyzerType.Orbitrap, 10,
                        1, null, NativeId);
                    dummyDDAMs2Scans.Add(dummyScan);
                }

                // create out of order dda
                dummyDDAScansOutOfOrder.AddRange(DummyDDAScansInOrder.Take(5));
                Stack<MsDataScan> dummyMs1ScanStack = new Stack<MsDataScan>(dummyDDAMs1Scans.GetRange(2, 3));
                Stack<MsDataScan> dummyMs2ScanStack = new Stack<MsDataScan>(dummyDDAMs2Scans.GetRange(3, 12));
                MsDataScan scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(6);
                scan.SetOneBasedPrecursorScanNumber(5);
                dummyDDAScansOutOfOrder.Add(scan);
                scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(7);
                scan.SetOneBasedPrecursorScanNumber(5);
                dummyDDAScansOutOfOrder.Add(scan);

                scan = dummyMs1ScanStack.Pop();
                scan.SetOneBasedScanNumber(8);
                dummyDDAScansOutOfOrder.Add(scan);
                scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(9);
                scan.SetOneBasedPrecursorScanNumber(5);
                dummyDDAScansOutOfOrder.Add(scan);
                scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(10);
                scan.SetOneBasedPrecursorScanNumber(8);
                dummyDDAScansOutOfOrder.Add(scan);
                scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(11);
                scan.SetOneBasedPrecursorScanNumber(8);
                dummyDDAScansOutOfOrder.Add(scan);
                scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(12);
                scan.SetOneBasedPrecursorScanNumber(8);
                dummyDDAScansOutOfOrder.Add(scan);

                scan = dummyMs1ScanStack.Pop();
                scan.SetOneBasedScanNumber(13);
                dummyDDAScansOutOfOrder.Add(scan);
                scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(14);
                scan.SetOneBasedPrecursorScanNumber(13);
                dummyDDAScansOutOfOrder.Add(scan);

                scan = dummyMs1ScanStack.Pop();
                scan.SetOneBasedScanNumber(15);
                dummyDDAScansOutOfOrder.Add(scan);
                scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(16);
                scan.SetOneBasedPrecursorScanNumber(15);
                dummyDDAScansOutOfOrder.Add(scan);
                scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(17);
                scan.SetOneBasedPrecursorScanNumber(15);
                dummyDDAScansOutOfOrder.Add(scan);
                scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(18);
                scan.SetOneBasedPrecursorScanNumber(15);
                dummyDDAScansOutOfOrder.Add(scan);
                scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(19);
                scan.SetOneBasedPrecursorScanNumber(13);
                dummyDDAScansOutOfOrder.Add(scan);
                scan = dummyMs2ScanStack.Pop();
                scan.SetOneBasedScanNumber(20);
                scan.SetOneBasedPrecursorScanNumber(15);
                dummyDDAScansOutOfOrder.Add(scan);

                List<MsDataScan> formattedScans = new();
                var retentionTimeIndex = 1;
                foreach (var dataScan in dummyDDAScansOutOfOrder.OrderBy(p => p.OneBasedScanNumber))
                {
                    MsDataScan newScan = new(dataScan.MassSpectrum, dataScan.OneBasedScanNumber, dataScan.MsnOrder,
                        dataScan.IsCentroid, dataScan.Polarity, retentionTimeIndex,
                        dataScan.ScanWindowRange, dataScan.ScanFilter, dataScan.MzAnalyzer, dataScan.TotalIonCurrent,
                        dataScan.InjectionTime, dataScan.NoiseData, dataScan.NativeId,
                        oneBasedPrecursorScanNumber: dataScan.OneBasedPrecursorScanNumber);
                    retentionTimeIndex++;
                    formattedScans.Add(newScan);
                }
                return formattedScans;
            }
        }

        #endregion

        public static List<MsDataScan> ActualScans => MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"AveragingTestData\TDYeastFractionMS1.mzML")).GetAllScansList().Take(50).ToList();

        public static string NativeId;

        public static SpectralAveragingParameters SpectralAveragingParameters { get; set; }

        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            SpectralAveragingParameters = new();
            SpectralAveragingParameters.OutlierRejectionType = OutlierRejectionType.NoRejection;
            SpectralAveragingParameters.SpectralWeightingType = SpectraWeightingType.WeightEvenly;
            SpectralAveragingParameters.NormalizationType = NormalizationType.NoNormalization;
            SpectralAveragingParameters.MaxThreadsToUsePerFile = 12;
            NativeId = ActualScans.First().NativeId;
        }


        [Test]
        public static void TestAverageAll()
        {
            SpectralAveragingParameters.SpectraFileAveragingType = SpectraFileAveragingType.AverageAll;
            MsDataScan[] averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            double[] expected = new double[] { 4, 8 };
            Assert.That(averagedScans.Length == 1);
            Assert.That(expected.SequenceEqual(averagedScans.First().MassSpectrum.YArray));

            averagedScans = SpectraFileAveraging.AverageSpectraFile(ActualScans.Take(30).ToList(), SpectralAveragingParameters);
            Assert.That(averagedScans.Length == 1);
        }

        [Test]
        public static void TestAverageEveryNScansWithoutOverlap()
        {
            SpectralAveragingParameters.SpectraFileAveragingType = SpectraFileAveragingType.AverageEverynScans;
            SpectralAveragingParameters.NumberOfScansToAverage = 5;
            MsDataScan[] averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            double[] expected = new double[] {4, 8};
            Assert.That(averagedScans.Length == 2);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            Assert.That(averagedScans[1].MassSpectrum.YArray.SequenceEqual(expected));

            SpectralAveragingParameters.ScanOverlap = 4;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            expected = new double[] { 4, 8 };
            Assert.That(averagedScans.Length == 2);
            Assert.That(SpectralAveragingParameters.ScanOverlap == 0);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            Assert.That(averagedScans[1].MassSpectrum.YArray.SequenceEqual(expected));

            averagedScans = SpectraFileAveraging.AverageSpectraFile(ActualScans, SpectralAveragingParameters);
            Assert.That(averagedScans.Length == 10);
        }

        [Test]
        public static void TestAverageEveryNScansWithOverlap()
        {
            SpectralAveragingParameters.SpectraFileAveragingType =
                SpectraFileAveragingType.AverageEverynScansWithOverlap;

            SpectralAveragingParameters.ScanOverlap = 1;
            SpectralAveragingParameters.NumberOfScansToAverage = 2;
            MsDataScan[] averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            double[] expected = new double[] { 5, 10 };
            Assert.That(averagedScans.Length == 9);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            SpectralAveragingParameters.NumberOfScansToAverage = 3;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            Assert.That(averagedScans.Length == 4);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            SpectralAveragingParameters.NumberOfScansToAverage = 4;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            Assert.That(averagedScans.Length == 3);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            SpectralAveragingParameters.NumberOfScansToAverage = 5;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            expected = new double[] { 4, 8 };
            Assert.That(averagedScans.Length == 2);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));


            SpectralAveragingParameters.ScanOverlap = 2;
            SpectralAveragingParameters.NumberOfScansToAverage = 3;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            expected = new double[] { 5, 10 };
            Assert.That(averagedScans.Length == 8);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            SpectralAveragingParameters.NumberOfScansToAverage = 4;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            Assert.That(averagedScans.Length == 4);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            SpectralAveragingParameters.NumberOfScansToAverage = 5;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            expected = new double[] { 4, 8 };
            Assert.That(averagedScans.Length == 2);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));

            SpectralAveragingParameters.ScanOverlap = 3;
            SpectralAveragingParameters.NumberOfScansToAverage = 4;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            expected = new double[] { 5, 10 };
            Assert.That(averagedScans.Length == 7);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
            SpectralAveragingParameters.NumberOfScansToAverage = 5;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyAllMs1Scans, SpectralAveragingParameters);
            expected = new double[] { 4, 8 };
            Assert.That(averagedScans.Length == 3);
            Assert.That(averagedScans.First().MassSpectrum.YArray.SequenceEqual(expected));
        }

        [Test]
        public static void TestAverageDDAScansWithoutOverlapInOrder()
        {
            SpectralAveragingParameters.SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans;
            SpectralAveragingParameters.NumberOfScansToAverage = 2;
            MsDataScan[] averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyDDAScansInOrder, SpectralAveragingParameters);
            double[] expected = new double[] { 5, 10 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 2);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 12);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            int?[] expectedNullable = new int?[] { null, 1, 1, 1, 1, 1, 1, null, 8, 8, 8, 8, 8, 8 };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));
        }

        [Test]
        public static void TestAverageDDAScansWithoutOverlapOutOfOrder()
        {
            SpectralAveragingParameters.SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans;
            SpectralAveragingParameters.NumberOfScansToAverage = 2;
            MsDataScan[] averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyDDAScansOutOfOrder, SpectralAveragingParameters);
            double[] expected = new double[] { 5, 10 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 2);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 11);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            int?[] expectedNullable = new int?[] { null, 1, 1, 1, 1, 1, 1, null, 8, 8, 8, 8, 8 };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));
        }

        [Test]
        public static void TestAverageDDAScansWithOverlapInOrder()
        {
            SpectralAveragingParameters.SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScansWithOverlap;
            SpectralAveragingParameters.NumberOfScansToAverage = 2;
            SpectralAveragingParameters.ScanOverlap = 1;
            MsDataScan[] averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyDDAScansInOrder, SpectralAveragingParameters);
            double[] expected = new double[] { 5, 10 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 4);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 15);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            expected = new double[] { 2.5, 5 };
            //Assert.AreEqual(expected, averagedScans.Last(p => p.MsnOrder == 1).MassSpectrum.YArray);
            int?[] expectedNullable = new int?[]
            {
                null, 1, 1, 1, null, 5, 5, 5, null, 9, 9, 9, null, 13, 13, 13, 13, 13, 13
            };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));

            SpectralAveragingParameters.NumberOfScansToAverage = 3;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyDDAScansInOrder, SpectralAveragingParameters);
            expected = new double[] { 5, 10 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 2);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 15);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            expected = new double[] { 5.0 * (2.0 / 3.0), 10.0 * (2.0 / 3.0) };
            Assert.AreEqual(expected.Select(p => Math.Round(p, 4)), averagedScans.Last(p => p.MsnOrder == 1).MassSpectrum.YArray.Select(p => Math.Round(p, 4)));
            expectedNullable = new int?[]
            {
                null, 1, 1, 1, 1, 1, 1, null, 8, 8, 8, 8, 8, 8, 8, 8, 8
            };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));

            SpectralAveragingParameters.NumberOfScansToAverage = 3;
            SpectralAveragingParameters.ScanOverlap = 2;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyDDAScansInOrder, SpectralAveragingParameters);
            expected = new double[] { 5, 10 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 3);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 15);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            expected = new double[] { 5.0 * (2.0 / 3.0), 10.0 * (2.0 / 3.0) };
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
            SpectralAveragingParameters.SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScansWithOverlap;
            SpectralAveragingParameters.NumberOfScansToAverage = 2;
            SpectralAveragingParameters.ScanOverlap = 1;
            MsDataScan[] averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyDDAScansOutOfOrder, SpectralAveragingParameters);
            double[] expected = new double[] { 5, 10 };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 4);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 15);
            Assert.That(averagedScans[0].MassSpectrum.YArray.SequenceEqual(expected));
            Assert.AreEqual(expected, averagedScans.Last(p => p.MsnOrder == 1).MassSpectrum.YArray);
            int?[] expectedNullable = new int?[]
            {
                null, 1, 1, 1, null, 5, 5, 5, null, 9, 9, 9, null, 13, 13, 13, 13, 13, 13
            };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));

            SpectralAveragingParameters.NumberOfScansToAverage = 3;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyDDAScansOutOfOrder, SpectralAveragingParameters);
            expected = new double[] { 5.0 * (2.0 / 3.0), 10.0 * (2.0 / 3.0) };
            Assert.That(averagedScans.Count(p => p.MsnOrder == 1) == 2);
            Assert.That(averagedScans.Count(p => p.MsnOrder == 2) == 15);
            Assert.AreEqual(expected.Select(p => Math.Round(p, 4)), averagedScans[0].MassSpectrum.YArray.Select(p => Math.Round(p, 4)));
            Assert.AreEqual(expected.Select(p => Math.Round(p, 4)), averagedScans.Last(p => p.MsnOrder == 1).MassSpectrum.YArray.Select(p => Math.Round(p, 4)));
            expectedNullable = new int?[]
            {
                null, 1, 1, 1, 1, 1, 1, null, 8, 8, 8, 8, 8, 8, 8, 8, 8
            };
            Assert.That(averagedScans.Select(p => p.OneBasedPrecursorScanNumber).ToArray().SequenceEqual(expectedNullable));

            SpectralAveragingParameters.NumberOfScansToAverage = 3;
            SpectralAveragingParameters.ScanOverlap = 2;
            averagedScans = SpectraFileAveraging.AverageSpectraFile(DummyDDAScansOutOfOrder, SpectralAveragingParameters);
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
    }
}
