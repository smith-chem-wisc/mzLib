using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Windows.Automation.Peers;
using Easy.Common.Extensions;
using FlashLFQ.Alex_project;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using OxyPlot.Axes;
using Readers;
using Plotly.NET;
using Plotly.NET.CSharp;
using Chart = Plotly.NET.CSharp.Chart;
using GenericChartExtensions = Plotly.NET.CSharp.GenericChartExtensions;
using Plotly.NET.LayoutObjects;
using Plotly.NET.TraceObjects;
using System.ComponentModel;
using System.Drawing;
using Color = Plotly.NET.Color;
using System.Collections;
using System.Security.Policy;



namespace Test.FileReadingTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestRawFileReader
    {
        [Test]
        public void TestFileDoesntExist()
        {
            string fakePath = "fakePath.raw";
            var reader = MsDataFileReader.GetDataFile(fakePath);
            Assert.Throws<FileNotFoundException>(() =>
            {
                reader.InitiateDynamicConnection();
            });
        }

        #region Testing Exceptions

        [Test]
        public void TestRawFileReaderFileNotFoundException()
        {
            var fakeRawFile = "asdasd.raw";

            var ex = Assert.Throws<FileNotFoundException>(() => MsDataFileReader.GetDataFile(fakeRawFile).LoadAllStaticData());

            Assert.That(ex.Message, Is.EqualTo(new FileNotFoundException().Message));
        }

        #endregion


        [Test]
        public void TestScanDescription()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "ScanDescriptionTestData.raw");
            var scans = MsDataFileReader.GetDataFile(filePath).GetAllScansList();
            var ms1Scans = scans.Where(x => x.MsnOrder == 1).ToList();
            var ms2Scans = scans.Where(x => x.MsnOrder == 2).ToList();

            ms1Scans.ForEach(x => Assert.That(x.ScanDescription, Is.EqualTo(null)));
            ms2Scans.ForEach(x => Assert.That(x.ScanDescription, Is.EqualTo("Testing2")));
        }

        /// <summary>
        /// Tests LoadAllStaticData for ThermoRawFileReader
        /// </summary>
        /// <param name="infile"></param>
        /// <param name="outfile1"></param>
        /// <param name="outfile2"></param>
        [Test]
        [TestCase("testFileWMS2.raw", "a.mzML", "aa.mzML")]
        [TestCase("small.raw", "a.mzML", "aa.mzML")]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw", "a.mzML", "aa.mzML")]
        public static void TestLoadAllStaticDataRawFileReader(string infile, string outfile1, string outfile2)
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", infile);
            outfile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", outfile1);
            outfile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", outfile2);

            string dummyPath = "aaljienxmbelsiemxmbmeba.raw";
            Assert.Throws<FileNotFoundException>(() =>
            {
                var dummyReader = MsDataFileReader.GetDataFile(dummyPath);
                dummyReader.LoadAllStaticData();
            });

            // testing load with multiple threads 
            var parallelReader = MsDataFileReader.GetDataFile(path);
            parallelReader.LoadAllStaticData(null, maxThreads: 4);

            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData(null, maxThreads: 1);
            reader.ExportAsMzML(outfile1, false);
            reader.LoadAllStaticData();
            reader.ExportAsMzML(outfile2, true);
            var readerMzml = MsDataFileReader.GetDataFile(outfile2);
            readerMzml.LoadAllStaticData();
            Console.WriteLine($"Analysis time for TestLoadAllStaticDataRawFileReader({infile}): {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
        }

        // for one ms1 scan pairing
        public static List<(double,double)> CombinedMassSpectrumWithTolerance(double[] FirstMz, double[] SecondMz, double[]FirstIntensity, double[]SecondIntensity, double ppmTolerance)
        {
            int i = 0;
            int j = 0;
            List<(double,double)> mzIntensityPairs = new List<(double, double)>();
            while (i < FirstMz.Length && j < SecondMz.Length)
            {
                double diff = Math.Abs((FirstMz[i] - SecondMz[j]) / SecondMz[j] *1000000);
                if (diff <= ppmTolerance)
                {
                    double newMZ = (FirstMz[i]* FirstIntensity[i] + SecondMz[j]* SecondIntensity[j]) / (FirstIntensity[i] + SecondIntensity[j]);
                    mzIntensityPairs.Add((newMZ, FirstIntensity[i] + SecondIntensity[j]));
                    i++;
                    j++;
                }
                else if (FirstMz[i] < SecondMz[j])
                {
                    mzIntensityPairs.Add(((FirstMz[i]) / 1, (FirstIntensity[i]) / 1));
                    i++;
                }
                else
                {
                    mzIntensityPairs.Add(((SecondMz[j]) / 1, (SecondIntensity[j]) / 1));
                    j++;
                }
            }
            return mzIntensityPairs;
        }

        public static List<MsDataScan> combinedMS1(List<MsDataScan> first, List<MsDataScan> second) 
        {
            List<MsDataScan> combinedScans = new List<MsDataScan>();
            for (int i = 0; i < first.Count; i++)
            {
                bool added = false;
                for (int j = 0; j < second.Count; j++)
                {
                    if (first[i].OneBasedScanNumber == second[j].OneBasedScanNumber)
                    {
                        var one = first[i];
                        var two = second[j];
                        var newSpectrum = CombinedMassSpectrumWithTolerance(one.MassSpectrum.XArray, two.MassSpectrum.XArray, one.MassSpectrum.YArray, two.MassSpectrum.YArray, 10);
                        double[] mz = newSpectrum.Select(x => x.Item1).ToArray();
                        double[] intensities = newSpectrum.Select(x => x.Item2).ToArray();
                        bool shouldCopy = true;
                        var combinedMassSpectrum = new MzSpectrum(mz, intensities, shouldCopy);
                        var combinedMs1 = new MsDataScan(combinedMassSpectrum, one.OneBasedScanNumber, one.MsnOrder, one.IsCentroid, one.Polarity, one.RetentionTime, one.ScanWindowRange,
                            one.ScanFilter, one.MzAnalyzer, one.TotalIonCurrent, one.InjectionTime, one.NoiseData, one.NativeId, one.SelectedIonMZ, one.SelectedIonChargeStateGuess,
                            one.SelectedIonIntensity, one.IsolationMz, one.IsolationWidth, one.DissociationType, one.OneBasedPrecursorScanNumber, one.SelectedIonMonoisotopicGuessMz, one.HcdEnergy, one.ScanDescription);
                        combinedScans.Add(combinedMs1);
                        added = true;
                    }
                }

                if (!added) 
                {
                    combinedScans.Add(first[i]);
                }
            }
            return combinedScans;
        }

        public static List<MsDataScan> buildScans(List<MsDataScan> Ms1scans, List<MsDataScan> Ms2scans, double Shifttime = 0)
        {
            List<MsDataScan> scansForTheNewFile = new List<MsDataScan>();


            foreach (var scan in Ms1scans)
            {
                scansForTheNewFile.Add(scan);
            }
            foreach (var scan in Ms2scans)
            {
                scansForTheNewFile.Add(scan);
            }
            scansForTheNewFile = scansForTheNewFile.OrderBy(x => x.OneBasedScanNumber).ToList();

            if (Shifttime != 0)
            {
                scansForTheNewFile = SwitchRetentionTime(scansForTheNewFile, Shifttime);
            }

            return scansForTheNewFile;


        }


        public static List<MsDataScan> SwitchRetentionTime(List<MsDataScan> MsScans, double timeShift)
        {
            List<MsDataScan> output = new List<MsDataScan>();

            foreach (var one in MsScans)
            {
                var modifiedMs1 = new MsDataScan(one.MassSpectrum, one.OneBasedScanNumber, one.MsnOrder, one.IsCentroid, one.Polarity, one.RetentionTime + timeShift, one.ScanWindowRange,
                one.ScanFilter, one.MzAnalyzer, one.TotalIonCurrent, one.InjectionTime, one.NoiseData, one.NativeId, one.SelectedIonMZ, one.SelectedIonChargeStateGuess,
                one.SelectedIonIntensity, one.IsolationMz, one.IsolationWidth, one.DissociationType, one.OneBasedPrecursorScanNumber, one.SelectedIonMonoisotopicGuessMz, one.HcdEnergy, one.ScanDescription);
                output.Add(modifiedMs1);
            }

            return output;
        }
        public static List<MsDataScan> ShiftMs1(List<MsDataScan> Ms1Scans,int crossNumber)
        {
            List<MsDataScan> shiftedScans = new List<MsDataScan>();
            for (int i = 0; i < Ms1Scans.Count; i++)
            {
                var one = Ms1Scans[i];
                if (i + crossNumber >= Ms1Scans.Count || i + crossNumber < 0)
                {
                    continue;
                }
                int scanNumber_crossed = Ms1Scans[i + crossNumber].OneBasedScanNumber;
                double retentionTime_crossed = Ms1Scans[i + crossNumber].RetentionTime;
                var newSpectrum = new MzSpectrum(one.MassSpectrum.XArray, one.MassSpectrum.YArray, true);
                var shiftedMs1 = new MsDataScan(newSpectrum, scanNumber_crossed, one.MsnOrder, one.IsCentroid, one.Polarity, retentionTime_crossed, one.ScanWindowRange,
                one.ScanFilter, one.MzAnalyzer, one.TotalIonCurrent, one.InjectionTime, one.NoiseData, one.NativeId, one.SelectedIonMZ, one.SelectedIonChargeStateGuess,
                one.SelectedIonIntensity, one.IsolationMz, one.IsolationWidth, one.DissociationType, one.OneBasedPrecursorScanNumber, one.SelectedIonMonoisotopicGuessMz, one.HcdEnergy, one.ScanDescription);
                shiftedScans.Add(shiftedMs1);
            }
            return shiftedScans;
        }

        public static List<(double, double)> XIC_convert(List<MsDataScan> Ms1Scan, double mass)
        {
            List<(double, double)> XIC = new List<(double, double)>();

            for (int i = 0; i < Ms1Scan.Count; i++)
            {
                var one = Ms1Scan[i];
                var mz = one.MassSpectrum.XArray;
                var intensity = one.MassSpectrum.YArray;
                for (int j = 0; j < mz.Length; j++)
                {
                    if (Math.Abs(mz[j] - mass) < 0.2)
                    {
                        XIC.Add((one.RetentionTime, intensity[j]));
                    }
                }
            }
            return XIC;
        }


        [Test]
        public void test()
        {
            string filePath = "E:\\MBR\\testRawFile\\20100604_Velos1_TaGe_SA_A549_3.raw";
            var outfile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "out1.mzml");
            var outfile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "out2.mzml");
            var reader = MsDataFileReader.GetDataFile(filePath);
            reader.LoadAllStaticData(null, maxThreads: 1);

            List<MsDataScan> ms1scans_OriginPeak = reader.GetAllScansList().Where(x => x.MsnOrder == 1).ToList();
            List<MsDataScan> ms1scans_BackPeak = ShiftMs1(ms1scans_OriginPeak, 50);
            List<MsDataScan> ms1scans_FrontPeak = ShiftMs1(ms1scans_OriginPeak, -50);

            var ms2scans = reader.GetAllScansList().Where(x => x.MsnOrder == 2).ToList();

            var ms1scans_combined = combinedMS1(ms1scans_OriginPeak, ms1scans_BackPeak);
            var ms1scans_combined2 = combinedMS1(ms1scans_OriginPeak, ms1scans_FrontPeak);
            var ms1scans_combined2_forDraw = SwitchRetentionTime(ms1scans_combined2, 3);

           var scansForTheNewFile = buildScans(ms1scans_combined, ms2scans);
            var scansForTheNewFile2 = buildScans(ms1scans_combined2, ms2scans,3);

            string outPath = filePath.Replace(".raw", "_first.mzML").ToString();
            string outPath_2 = filePath.Replace(".raw", "_second.mzML").ToString();

            SourceFile sourceFile = new SourceFile(reader.SourceFile.NativeIdFormat,
                reader.SourceFile.MassSpectrometerFileFormat, reader.SourceFile.CheckSum, reader.SourceFile.FileChecksumType, reader.SourceFile.Uri, reader.SourceFile.Id, reader.SourceFile.FileName);


            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(scansForTheNewFile.ToArray(), sourceFile), outPath, false);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(scansForTheNewFile2.ToArray(), sourceFile), outPath_2, false);











            // make a chart
            var XIC_firstRun = XIC_convert(ms1scans_OriginPeak, 1200.0);
            var XIC_BackPeak = XIC_convert(ms1scans_BackPeak, 1200.0);
            var XIC_FrontPeak = XIC_convert(ms1scans_FrontPeak, 1200.0);
            var XIC_combined = XIC_convert(ms1scans_combined, 1200.0);
            var XIC_combined2 = XIC_convert(ms1scans_combined2_forDraw, 1200.0);



            var individualPlot = Chart.Combine(new[]
            {
                Chart.Scatter<double, double, string>(XIC_firstRun.Select(p=>p.Item1), XIC_firstRun.Select(p=>p.Item2),
                StyleParam.Mode.Lines_Markers, MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(ColorKeyword.Blue), Name: "original"),
                Chart.Scatter<double, double, string>(XIC_BackPeak.Select(p=>p.Item1), XIC_BackPeak.Select(p=>p.Item2),
                StyleParam.Mode.Lines_Markers, MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(ColorKeyword.Red), Name: "original_Back"),
                Chart.Scatter<double, double, string>(XIC_FrontPeak.Select(p=>p.Item1), XIC_FrontPeak.Select(p=>p.Item2),
                StyleParam.Mode.Lines_Markers, MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(ColorKeyword.Green), Name: "original_Front"),
            })
                .WithTitle("Data XIC")
                .WithXAxisStyle<double, double, string>(Title: Title.init("Times"))
                .WithYAxisStyle<double, double, string>(Title: Title.init("Intensities"))
            .WithSize(1600, 400);


            var combinedPlot = Chart.Combine(new[]
            {
                Chart.Scatter<double, double, string>(XIC_combined.Select(p=>p.Item1), XIC_combined.Select(p=>p.Item2),
                StyleParam.Mode.Lines_Markers, MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(ColorKeyword.Black), Name: "original_combined"),
            })
                .WithTitle("Data XIC")
                .WithXAxisStyle<double, double, string>(Title: Title.init("Times"))
                .WithYAxisStyle<double, double, string>(Title: Title.init("Intensities"))
            .WithSize(1600, 400);


            var combinedPlot2 = Chart.Combine(new[]
            {
                Chart.Scatter<double, double, string>(XIC_combined2.Select(p=>p.Item1), XIC_combined2.Select(p=>p.Item2),
                StyleParam.Mode.Lines_Markers, MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(ColorKeyword.Black), Name: "original_combined2"),
            })
                .WithTitle("Data XIC")
                .WithXAxisStyle<double, double, string>(Title: Title.init("Times"))
                .WithYAxisStyle<double, double, string>(Title: Title.init("Intensities"))
            .WithSize(1600, 400);


            var stacked = Chart.Grid(new[] { individualPlot, combinedPlot, combinedPlot2 }, 3, 1)
                .WithSize(1600, 1200);
            GenericChartExtensions.Show(stacked);


            //reader.ExportAsMzML(outfile1, false);
            //reader.LoadAllStaticData();
            //reader.ExportAsMzML(outfile2, true);
            //var readerMzml = MsDataFileReader.GetDataFile(outfile2);
            //readerMzml.LoadAllStaticData();
        }


        [Test]
        public void TestThermoGetSourceFile()
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "small.raw");
            var reader = MsDataFileReader.GetDataFile(path);
            SourceFile sf = reader.GetSourceFile();
            Assert.That(sf.NativeIdFormat, Is.EqualTo(@"Thermo nativeID format"));
        }

        /// <summary>
        /// Tests the dynamic connection for thermorawfilereader
        /// </summary>
        [Test]
        public static void TestDynamicConnectionRawFileReader()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            var path1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "small.raw");
            var thermoDynamic1 = MsDataFileReader.GetDataFile(path1);
            thermoDynamic1.InitiateDynamicConnection();

            var path2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "testFileWMS2.raw");
            var thermoDynamic2 = MsDataFileReader.GetDataFile(path2);
            thermoDynamic2.InitiateDynamicConnection();

            var msOrders = thermoDynamic1.GetMsOrderByScanInDynamicConnection();
            Assert.That(msOrders != null && msOrders.Length > 0);

            var a = thermoDynamic1.GetOneBasedScanFromDynamicConnection(1);
            Assert.That(a != null);

            var b = thermoDynamic2.GetOneBasedScanFromDynamicConnection(1);
            Assert.That(b != null);

            Assert.That(a.MassSpectrum.XArray.Length != b.MassSpectrum.XArray.Length);

            a = thermoDynamic1.GetOneBasedScanFromDynamicConnection(10000);
            thermoDynamic1.CloseDynamicConnection();
            thermoDynamic2.CloseDynamicConnection();

            Console.WriteLine($"Analysis time for TestDynamicConnectionRawFileReader: {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
            Assert.That(a == null);
        }

        /// <summary>
        /// Tests peak filtering for ThermoRawFileReader
        /// </summary>
        /// <param name="infile"></param>
        [Test]
        [TestCase("testFileWMS2.raw")]
        [TestCase("small.raw")]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw")]
        public static void TestPeakFilteringRawFileReader(string infile)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            var filterParams = new FilteringParams(200, 0.01, 0, 1, false, true, true);

            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", infile);
            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData(filterParams, maxThreads: 1);
            var rawScans = reader.GetAllScansList();
            foreach (var scan in rawScans)
            {
                Assert.That(scan.MassSpectrum.XArray.Length <= 200);
            }

            string outfile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", Path.GetFileNameWithoutExtension(infile) + ".mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(reader, outfile1, false);
            var mzml = MsDataFileReader.GetDataFile(outfile1);
            mzml.LoadAllStaticData(filterParams, maxThreads: 1);

            var mzmlScans = mzml.GetAllScansList();
            for (int i = 0; i < mzmlScans.Count; i++)
            {
                var mzmlScan = mzmlScans[i];
                var rawScan = rawScans[i];

                for (int j = 0; j < mzmlScan.MassSpectrum.XArray.Length; j++)
                {
                    double roundedMzmlMz = Math.Round(mzmlScan.MassSpectrum.XArray[j], 2);
                    double roundedRawMz = Math.Round(rawScan.MassSpectrum.XArray[j], 2);

                    Assert.AreEqual(roundedMzmlMz, roundedRawMz);

                    double roundedMzmlIntensity = Math.Round(mzmlScan.MassSpectrum.XArray[j], 0);
                    double roundedRawIntensity = Math.Round(rawScan.MassSpectrum.XArray[j], 0);

                    Assert.AreEqual(roundedMzmlIntensity, roundedRawIntensity);
                }
            }

            Console.WriteLine($"Analysis time for TestPeakFilteringRawFileReader: {stopwatch.Elapsed.Hours}h " +
                $"{stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
        }

        /// <summary>
        /// Test Thermo License for ThermoRawFileReader
        /// </summary>
        [Test]
        public static void TestThermoLicence()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            var licence = ThermoRawFileReaderLicence.ThermoLicenceText;
            Assert.That(licence.Length > 100);

            Console.WriteLine($"Analysis time for TestThermoLicence: {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
        }

        /// <summary>
        /// Test that raw files can be opened dynamically in ThermoRawFileReader
        /// </summary>
        /// <param name="fileName"></param>
        [Test]
        [TestCase("small.RAW")]
        [TestCase("testFileWMS2.raw")]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw")]
        public static void TestDynamicRaw(string fileName)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);

            var staticRaw = MsDataFileReader.GetDataFile(filePath);
            staticRaw.LoadAllStaticData();
            staticRaw.InitiateDynamicConnection();

            foreach (MsDataScan staticScan in staticRaw.GetAllScansList())
            {
                MsDataScan dynamicScan = staticRaw.GetOneBasedScanFromDynamicConnection(staticScan.OneBasedScanNumber);

                Assert.IsFalse(staticScan.MassSpectrum.YArray.Contains(0));
                Assert.IsFalse(dynamicScan.MassSpectrum.YArray.Contains(0));
                Assert.That(dynamicScan.OneBasedScanNumber == staticScan.OneBasedScanNumber);
                Assert.That(dynamicScan.MsnOrder == staticScan.MsnOrder);
                Assert.That(dynamicScan.RetentionTime == staticScan.RetentionTime);
                Assert.That(dynamicScan.Polarity == staticScan.Polarity);
                Assert.That(dynamicScan.ScanWindowRange.Minimum == staticScan.ScanWindowRange.Minimum);
                Assert.That(dynamicScan.ScanWindowRange.Maximum == staticScan.ScanWindowRange.Maximum);
                Assert.That(dynamicScan.ScanFilter == staticScan.ScanFilter);
                Assert.That(dynamicScan.NativeId == staticScan.NativeId);
                Assert.That(dynamicScan.IsCentroid == staticScan.IsCentroid);
                Assert.That(dynamicScan.IsCentroid == staticScan.IsCentroid);
                Assert.That(dynamicScan.InjectionTime == staticScan.InjectionTime);
                Assert.That(dynamicScan.NoiseData == staticScan.NoiseData);

                Assert.That(dynamicScan.IsolationMz == staticScan.IsolationMz);
                Assert.That(dynamicScan.SelectedIonChargeStateGuess == staticScan.SelectedIonChargeStateGuess);
                Assert.That(dynamicScan.SelectedIonIntensity == staticScan.SelectedIonIntensity);
                Assert.That(dynamicScan.SelectedIonMZ == staticScan.SelectedIonMZ);
                Assert.That(dynamicScan.DissociationType == staticScan.DissociationType);
                Assert.That(dynamicScan.IsolationWidth == staticScan.IsolationWidth);
                Assert.That(dynamicScan.OneBasedPrecursorScanNumber == staticScan.OneBasedPrecursorScanNumber);
                Assert.That(dynamicScan.SelectedIonMonoisotopicGuessIntensity == staticScan.SelectedIonMonoisotopicGuessIntensity);
                Assert.That(dynamicScan.SelectedIonMonoisotopicGuessMz == staticScan.SelectedIonMonoisotopicGuessMz);

                if (dynamicScan.IsolationRange != null || staticScan.IsolationRange != null)
                {
                    Assert.That(dynamicScan.IsolationRange.Minimum == staticScan.IsolationRange.Minimum);
                    Assert.That(dynamicScan.IsolationRange.Maximum == staticScan.IsolationRange.Maximum);
                }

                Assert.That(dynamicScan.MassSpectrum.XArray.Length == staticScan.MassSpectrum.XArray.Length);
                Assert.That(dynamicScan.MassSpectrum.YArray.Length == staticScan.MassSpectrum.YArray.Length);

                for (int i = 0; i < staticScan.MassSpectrum.XArray.Length; i++)
                {
                    double staticMz = staticScan.MassSpectrum.XArray[i];
                    double staticIntensity = staticScan.MassSpectrum.YArray[i];

                    double dynamicMz = dynamicScan.MassSpectrum.XArray[i];
                    double dynamicIntensity = dynamicScan.MassSpectrum.YArray[i];

                    Assert.That(dynamicMz == staticMz);
                    Assert.That(dynamicIntensity == staticIntensity);
                }
            }
        }

        /// <summary>
        /// Tests that you can read EtHCD files in ThermoRawFileReader
        /// </summary>
        [Test]
        public static void TestEthcdReading()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "sliced_ethcd.raw");
            var spectra = MsDataFileReader.GetDataFile(filePath);
            spectra.LoadAllStaticData(null, 1);
            var hcdScan = spectra.GetOneBasedScan(5);
            Assert.That(hcdScan.DissociationType == DissociationType.HCD);
            var ethcdScan = spectra.GetOneBasedScan(6);
            Assert.That(ethcdScan.DissociationType == DissociationType.EThcD);
        }
    }
}