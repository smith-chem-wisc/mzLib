using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Windows.Automation.Peers;
using Easy.Common.Extensions;
using FlashLFQ.Alex_project;
using MassSpectrometry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
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
using System.Collections.Generic;
using System.Xml.Schema;
using FlashLFQ;
using Microsoft.FSharp.Collections;
using OxyPlot.Wpf;
using mzPlot;
using Readers.Generated;
using System.Reflection;



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


        /// <summary>
        /// Add the same MS1 scans to the existed MS1 scans with the specific time shift. 
        /// Ex. PeakA in 3 min, and we add the same PeakA in 3 min later, the output will be PeakA in 3 min and PeakA in 6 min
        /// </summary>
        /// <param name="msScans"> The input MS1 scans</param>
        /// <param name="timeShift">The time shift to inserted the MS1 scans</param>
        /// <returns></returns>
        public static List<MsDataScan> AddingPeak(List<MsDataScan> msScans, double timeShift)
        {
            List<MsDataScan> output = new List<MsDataScan>();
            int InsertedCounter = 0;
            int currentCounter = 0;
            int Ms1TracingNumber = 0;

            foreach (var ms1Scan in msScans.Where(p => p.MsnOrder == 1))
            {

                if (ms1Scan.RetentionTime + timeShift > msScans.Last().RetentionTime)
                {
                    break;
                }

                foreach (var scan in msScans.Where(p => p.OneBasedScanNumber > currentCounter))
                {
                    if (scan.MsnOrder == 1)
                    {
                        Ms1TracingNumber = scan.OneBasedScanNumber + InsertedCounter;
                    }

                    if (scan.RetentionTime < ms1Scan.RetentionTime + timeShift)
                    {
                        MsDataScan newScan = new MsDataScan(scan.MassSpectrum, scan.OneBasedScanNumber + InsertedCounter, scan.MsnOrder, scan.IsCentroid, scan.Polarity, scan.RetentionTime, scan.ScanWindowRange,
                            scan.ScanFilter, scan.MzAnalyzer, scan.TotalIonCurrent, scan.InjectionTime, scan.NoiseData, scan.NativeId, scan.SelectedIonMZ, scan.SelectedIonChargeStateGuess,
                            scan.SelectedIonIntensity, scan.IsolationMz, scan.IsolationWidth, scan.DissociationType, Ms1TracingNumber, scan.SelectedIonMonoisotopicGuessMz, scan.HcdEnergy, scan.ScanDescription);
                        output.Add(newScan);

                        currentCounter = scan.OneBasedScanNumber;
                    }
                    else
                    {
                        MsDataScan insertedScan = new MsDataScan(ms1Scan.MassSpectrum, currentCounter + InsertedCounter + 1, ms1Scan.MsnOrder, ms1Scan.IsCentroid, ms1Scan.Polarity, ms1Scan.RetentionTime + timeShift, ms1Scan.ScanWindowRange,
                            ms1Scan.ScanFilter, ms1Scan.MzAnalyzer, ms1Scan.TotalIonCurrent, ms1Scan.InjectionTime, ms1Scan.NoiseData, ms1Scan.NativeId, ms1Scan.SelectedIonMZ, ms1Scan.SelectedIonChargeStateGuess,
                            ms1Scan.SelectedIonIntensity, ms1Scan.IsolationMz, ms1Scan.IsolationWidth, ms1Scan.DissociationType, null, ms1Scan.SelectedIonMonoisotopicGuessMz, ms1Scan.HcdEnergy, ms1Scan.ScanDescription);
                        output.Add(insertedScan);
                        InsertedCounter++;
                        break;
                    }

                }

            }

            return output;
        }

        /// <summary>
        /// Shift the retention time of the whole MS scans (MS1 and MS2) with the specific time.
        /// </summary>
        /// <param name="MsScans"> The MS scans need to be shifted</param>
        /// <param name="timeShift">The shift time</param>
        /// <returns></returns>
        public static List<MsDataScan> SwitchRetentionTime(List<MsDataScan> MsScans, double timeShift)
        {
            List<MsDataScan> output = new List<MsDataScan>();

            foreach (var one in MsScans)
            {
                var modifiedMsScan = new MsDataScan(one.MassSpectrum, one.OneBasedScanNumber, one.MsnOrder, one.IsCentroid, one.Polarity, one.RetentionTime + timeShift, one.ScanWindowRange,
                one.ScanFilter, one.MzAnalyzer, one.TotalIonCurrent, one.InjectionTime, one.NoiseData, one.NativeId, one.SelectedIonMZ, one.SelectedIonChargeStateGuess,
                one.SelectedIonIntensity, one.IsolationMz, one.IsolationWidth, one.DissociationType, one.OneBasedPrecursorScanNumber, one.SelectedIonMonoisotopicGuessMz, one.HcdEnergy, one.ScanDescription);
                output.Add(modifiedMsScan);
            }

            return output;
        }

        public static XIC XIC_convert(List<MsDataScan> Ms2Scan, double mass, bool isReference, int tolerance = 20)
        {
            List<(double, double)> XIC = new List<(double, double)>();
            List<XIC> xICs = new List<XIC>();

            double xIC_tolerance = (double)tolerance/1000;

            for (int i = 0; i < Ms2Scan.Count; i++)
            {
                var scan = Ms2Scan[i];
                var mz = scan.MassSpectrum.XArray;
                var intensity = scan.MassSpectrum.YArray;
                Dictionary<double, double> peaks = new Dictionary<double, double>();
                for (int j = 0; j < mz.Length; j++)
                {
                    if (Math.Abs(mz[j] - mass) < xIC_tolerance)
                    {
                        peaks.Add(mz[j], intensity[j]);
                    }                   
                }

                if (peaks.Count()!= 0)
                {
                    var peak_intensity = peaks.OrderBy(p => p.Key - mass).First().Value; // the closest peak to the mass
                    XIC.Add((scan.RetentionTime, peak_intensity));
                }
            }

            List<IndexedMassSpectralPeak> peakList = new List<IndexedMassSpectralPeak>();
            for (int i = 0; i < XIC.Count(); i++) 
            {
                peakList.Add(new IndexedMassSpectralPeak(mass, XIC[i].Item2, 0, XIC[i].Item1));
            }

            
            return new XIC(peakList,mass, null, isReference);
        }

        /// <summary>
        /// Testing method and will not use in the metamorpheus. The function will load two raw files and generate the XIC. Then try to group them.
        /// Is good for examine the grouping quality.
        /// </summary>
        /// <param name="peakOne"></param>
        /// <param name="peakTwo"></param>
        /// <param name="mass"> The mass that used for generating the XIC</param>
        /// <returns></returns>
        public static XICGroups XICGrouping(List<List<MsDataScan>> MS1List, double mass, double countNum = 0.55, double resoultion = 0.03)
        {
            bool isFirst = true;
            List<XIC> xICs = new List<XIC>();
            foreach (var ms1Scan in MS1List)
            {
                if (isFirst) 
                {
                    xICs.Add(XIC_convert(ms1Scan, mass, true)); 
                    isFirst = false;
                }
                else
                {
                    xICs.Add(XIC_convert(ms1Scan, mass, false));
                }
            }


            XICGroups xICGroups = null;
            if (xICs.Count()>1 && !xICs.Any(p => p == null))
            {
                xICGroups = new XICGroups(xICs, countNum, resoultion);
            }

            return xICGroups;

        }


        public static void DrawPlot(List<List<MsDataScan>> peakGroup, double mass, bool sharedExteMode = false, bool smooth = false, List<string> nameList = null)
        {
            //XIC generating
            List<XIC> xICs = new List<XIC>();
            bool isFirst = true;
            foreach (var peak in peakGroup)
            {
                if(isFirst)
                {
                    xICs.Add(XIC_convert(peak, mass, true));
                    isFirst = false;
                }
                else
                {
                    xICs.Add(XIC_convert(peak, mass, false));
                }
            }

            if (nameList == null)
            {
                nameList = new List<string>();
                for (int i = 0; i < peakGroup.Count; i++)
                {
                    nameList.Add("Run" + i);
                }
            }

            //Plotting
            List<GenericChart> plots = new List<GenericChart>();
            List<ColorKeyword> colorKeywords = new List<ColorKeyword> { ColorKeyword.Red, ColorKeyword.Blue, ColorKeyword.Green, ColorKeyword.Orange, ColorKeyword.Purple, ColorKeyword.Pink, ColorKeyword.Brown, ColorKeyword.Cyan, ColorKeyword.Magenta, ColorKeyword.Gray, ColorKeyword.Black };
            int index = 0;

            if (!smooth && !sharedExteMode)
            {
                foreach (var xic in xICs)
                {
                    var PlotRef = Chart.Combine(new[]
                    {
                    Chart.Line<double, double, string>(xic.Ms1Peaks.Select(p=>p.RetentionTime), xic.Ms1Peaks.Select(p=>p.Intensity), Name: nameList[index],
                     MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[index])),
                })
                        .WithTitle("Data XIC")
                        .WithXAxisStyle<double, double, string>(Title: Title.init("Times"))
                        .WithYAxisStyle<double, double, string>(Title: Title.init("Intensities"))
                    .WithSize(1600, 200);

                    plots.Add(PlotRef);
                    index++;
                }
            }

            else if (sharedExteMode) 
            {
                var xicGorup = XICGrouping(peakGroup, mass, 0.3, 0.05);
                var PlotRef = Chart.Combine(new[]
                    {
                    Chart.Line<double, double, string>(xICs[0].SmoothedRetentionTime, xICs[0].SmoothedIntensity, Name: nameList[index],
                     MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[index])),
                    Chart.Point<double, double, string>(xicGorup.ExtremaInRef.Select(p=>p.Key), xicGorup.ExtremaInRef.Select(p=>p.Value), Name: nameList[index],
                                        MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[index+1])),
                }
                )
                        .WithTitle("Data XIC")
                        .WithXAxisStyle<double, double, string>(Title: Title.init("Times"))
                        .WithYAxisStyle<double, double, string>(Title: Title.init("Intensities"))
                    .WithSize(1600, 800);
                GenericChartExtensions.Show(PlotRef);
                int iii = 0;
            }

            


            else
            {
                foreach (var xic in xICs)
                {
                    var PlotRef = Chart.Combine(new[]
                    {
                    Chart.Line<double, double, string>(xic.SmoothedRetentionTime, xic.SmoothedIntensity, Name: nameList[index],
                     MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[index])),
                })
                        .WithTitle("Data XIC")
                        .WithXAxisStyle<double, double, string>(Title: Title.init("Times"))
                        .WithYAxisStyle<double, double, string>(Title: Title.init("Intensities"))
                    .WithSize(1600, 200);

                    plots.Add(PlotRef);
                    index++;
                }
            }
            

            var stacked = Chart.Grid(plots, plots.Count(), 1)
                .WithSize(1600, 200 * plots.Count());
            GenericChartExtensions.Show(stacked);

        }

        [Test]
        public void XICSmoothViewer() 
        {
            string file9 = "FS1_14373_3.raw";
            string file10 = "FS1_MS_14373_6_05082021.raw";
            string filePath9 = "E:\\GitClones\\mzLib\\mzLib\\TestFlashLFQ\\TestData\\isobaricFolder\\" + file9;
            string filePath10 = "E:\\GitClones\\mzLib\\mzLib\\TestFlashLFQ\\TestData\\isobaricFolder\\" + file10;
            var reader9 = MsDataFileReader.GetDataFile(filePath9);
            var reader10 = MsDataFileReader.GetDataFile(filePath10);
            reader9.LoadAllStaticData(null, maxThreads: 1);
            reader10.LoadAllStaticData(null, maxThreads: 1);
            List<MsDataScan> ms1scans9 = reader9.GetAllScansList().ToList();
            List<MsDataScan> ms1scans10 = reader10.GetAllScansList().ToList();
            var run9 = ms1scans9.Where(p => p.MsnOrder > 1).ToList();
            var run10 = ms1scans10.Where(p => p.MsnOrder > 1).ToList();
            //DrawPlot(new List<List<MsDataScan>> { run9,run10}, 518.2038, false, false, new List<string>(){file9,file10});
            DrawPlot(new List<List<MsDataScan>> { run9, run10 }, 1078.482, false ,false, new List<string>() { file9, file10 });
            var xicGroup = XICGrouping( new List<List<MsDataScan>>{ run9,run10 }, 1078.482, 0.55, 0.2);
            int iii = 0;
        }

        [Test]
        public void IsotopicViewer() 
        {
            string file = "FS1_14373_18_20210729230343.raw";
            string filePath = "E:\\GitClones\\mzLib\\mzLib\\TestFlashLFQ\\TestData\\isobaricFolder\\" + file;
            var reader = MsDataFileReader.GetDataFile(filePath);
            reader.LoadAllStaticData(null, maxThreads: 1);
            List<MsDataScan> ms1scans = reader.GetAllScansList().ToList().Where(p => p.MsnOrder == 1).ToList();
            List<double> isotopicevenlope = new List<double>() { 1078.48, 1078.73, 1078.98, 1079.23, 1079.49, 1079.74, 1079.99, 1080.24 };
            List<string> nameList = new List<string>() { "1078.48", "1078.73", "1078.98", "1079.23", "1079.49", "1079.74", "1079.99", "1080.24" };
            List<ColorKeyword> colorKeywords = new List<ColorKeyword> { ColorKeyword.Red, ColorKeyword.Blue, ColorKeyword.Green, ColorKeyword.Orange, 
                ColorKeyword.Purple, ColorKeyword.Pink, ColorKeyword.Brown, ColorKeyword.Cyan, ColorKeyword.Magenta, ColorKeyword.Gray, ColorKeyword.Black };
            
            List<XIC> xICs = new List<XIC>();
            foreach (var isopeak in isotopicevenlope) 
            {
                xICs.Add(XIC_convert(ms1scans, isopeak, false));
            }
         
            var plot = Chart.Combine(new[]
            {
                    Chart.Line<double, double, string>(xICs[0].Ms1Peaks.Select(p=>p.RetentionTime), xICs[0].Ms1Peaks.Select(p=>p.Intensity), Name: nameList[0],
                     MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[0])),
                    Chart.Line<double, double, string>(xICs[1].Ms1Peaks.Select(p=>p.RetentionTime), xICs[1].Ms1Peaks.Select(p=>p.Intensity), Name: nameList[1],
                     MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[1])),
                    Chart.Line<double, double, string>(xICs[2].Ms1Peaks.Select(p=>p.RetentionTime), xICs[2].Ms1Peaks.Select(p=>p.Intensity), Name: nameList[2],
                     MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[2])),
                    Chart.Line<double, double, string>(xICs[3].Ms1Peaks.Select(p=>p.RetentionTime), xICs[3].Ms1Peaks.Select(p=>p.Intensity), Name: nameList[3],
                     MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[3])),
                    Chart.Line<double, double, string>(xICs[4].Ms1Peaks.Select(p=>p.RetentionTime), xICs[4].Ms1Peaks.Select(p=>p.Intensity), Name: nameList[4],
                     MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[4])),
                    Chart.Line<double, double, string>(xICs[5].Ms1Peaks.Select(p=>p.RetentionTime), xICs[5].Ms1Peaks.Select(p=>p.Intensity), Name: nameList[5],
                     MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[5])),
                    Chart.Line<double, double, string>(xICs[6].Ms1Peaks.Select(p=>p.RetentionTime), xICs[6].Ms1Peaks.Select(p=>p.Intensity), Name: nameList[6],
                     MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[6])),
                    Chart.Line<double, double, string>(xICs[7].Ms1Peaks.Select(p=>p.RetentionTime), xICs[7].Ms1Peaks.Select(p=>p.Intensity), Name: nameList[7],
                     MarkerSymbol: StyleParam.MarkerSymbol.Circle, MarkerColor: Color.fromKeyword(colorKeywords[7])),
                })
                       .WithTitle("Data XIC")
                       .WithXAxisStyle<double, double, string>(Title: Title.init("Times"))
                       .WithYAxisStyle<double, double, string>(Title: Title.init("Intensities"))
                   .WithSize(1600, 800);
            GenericChartExtensions.Show(plot);
        }

        [Test]
        public void XicViewer() 
        {
            string file1 = "FS1_14373_3.raw";
            string file2 = "FS1_MS_14373_6_05082021.raw";
            string file3 = "FS1_MS_14373_9_05082021.raw";
            string file4 = "FS1_14373_12_20210729014106.raw";
            string file5 = "FS1_14373_15.raw";
            string file6 = "FS1_14373_18_20210729230343.raw";
            string file7 = "FS1_MS_14373_21_05082021.raw";
            string file8 = "FS1_MS_14373_24_05082021.raw";
            string file9 = "FS1_14373_27.raw";
            string file10 = "FS1_14373_30.raw";
            List<string> nameList = new List<string>(){ file1, file2, file3, file4, file5, file6, file7, file8, file9, file10};
            //List<string> nameList = new List<string>() {  file2, file3, file4, file8,  file10 };
            List<string> filePathList = new List<string>();
            List<MsDataFile> readerList = new List<MsDataFile>();
            List<List<MsDataScan>> ms1ScansList = new List<List<MsDataScan>>();


            foreach (string name in nameList) 
            {
                filePathList.Add("E:\\GitClones\\mzLib\\mzLib\\TestFlashLFQ\\TestData\\isobaricFolder\\" + name);
            }

            foreach (string path in filePathList)
            {
                var reader = MsDataFileReader.GetDataFile(path);
                reader.LoadAllStaticData(null, maxThreads: 1);
                readerList.Add(reader);
            }
            
            foreach (var reader in readerList)
            {
                ms1ScansList.Add(reader.GetAllScansList().Where(p => p.MsnOrder == 1).ToList());
            }

            var xicGroup = XICGrouping(ms1ScansList, 1078.4825, 0.3, 0.05);
            DrawPlot(ms1ScansList, 1078.4825, false, true, nameList);
            int iiii = 0;
        }

        [Test]
        public void test_generateIsobaricData_noRTDiff()
        {
            string filePath = "E:\\MBR\\testRawFile\\20100604_Velos1_TaGe_SA_A549_3.raw";
            var reader = MsDataFileReader.GetDataFile(filePath);
            reader.LoadAllStaticData(null, maxThreads: 1);

            List<MsDataScan> ms1scans_OriginPeak = reader.GetAllScansList().ToList();

            List<MsDataScan> File1 = AddingPeak(ms1scans_OriginPeak, 3);
            List<MsDataScan> File2 = AddingPeak(ms1scans_OriginPeak, -3);
            List<MsDataScan> File2_shifted = SwitchRetentionTime(File2, +3);
            
            var peakRef = ms1scans_OriginPeak.Where(p=>p.MsnOrder == 1).ToList();
            var peakOne = File1.Where(p=>p.MsnOrder == 1).ToList();
            var peakTwo = File2.Where(p=>p.MsnOrder == 1).ToList();
            var peakTwo_shifted = File2_shifted.Where(p => p.MsnOrder == 1).ToList();



            //var xicGroups = XICGrouping(peakOne, peakTwo_shifted, 1199.88);
            //DrawPlot(peakRef,peakOne, peakTwo_shifted, 1020.192, true);
            //DrawPlot(peakRef, peakOne, peakTwo_shifted, 1020.192, false);

            string outPath = filePath.Replace(".raw", "_first.mzML").ToString();
            string outPath_2 = filePath.Replace(".raw", "_second.mzML").ToString();

            SourceFile sourceFile = new SourceFile(reader.SourceFile.NativeIdFormat,
                reader.SourceFile.MassSpectrometerFileFormat, reader.SourceFile.CheckSum, reader.SourceFile.FileChecksumType, reader.SourceFile.Uri, reader.SourceFile.Id, reader.SourceFile.FileName);


            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(File1.ToArray(), sourceFile), outPath, false);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(File2_shifted.ToArray(), sourceFile), outPath_2, false);


            //reader.ExportAsMzML(outfile1, false);
            //reader.LoadAllStaticData();
            //reader.ExportAsMzML(outfile2, true);
            //var readerMzml = MsDataFileReader.GetDataFile(outfile2);
            //readerMzml.LoadAllStaticData();
        }

        [Test]
        public void test_generateIsobaricData_withRTDiff()
        {
            string filePath = "E:\\MBR\\testRawFile\\20100604_Velos1_TaGe_SA_A549_3.raw";
            var reader = MsDataFileReader.GetDataFile(filePath);
            reader.LoadAllStaticData(null, maxThreads: 1);

            List<MsDataScan> ms1scans_OriginPeak = reader.GetAllScansList().ToList();

            List<MsDataScan> File1 = AddingPeak(ms1scans_OriginPeak, 3);
            List<MsDataScan> File2 = AddingPeak(ms1scans_OriginPeak, -3);
            List<MsDataScan> File3 = SwitchRetentionTime(File2, +1);
            List<MsDataScan> File4 = SwitchRetentionTime(File2, +1.5);
            List<MsDataScan> File5 = SwitchRetentionTime(File2, +2);

            List<List<MsDataScan>> MSPeaks = new List<List<MsDataScan>> { File1, File2, File3, File4, File5};
            var Ms1Peaks = MSPeaks.Select(p => p.Where(p => p.MsnOrder == 1).ToList()).ToList();


            //int Good = 0;
            //int Bad = 0;
            //Dictionary<double, double> goodIndex = new Dictionary<double, double>();
            //Dictionary<double, double> badIndex = new Dictionary<double, double>();

            //for (double i = 300.0; i < 400.0; i = i + 2.0)
            //{
            //    var groups = XICGrouping(peakOne, peakTwo, i);

            //    if (groups != null && groups.RTDict[1] > 2)
            //    {
            //        Good++;
            //        goodIndex.Add(i, groups.RTDict[1]);
            //    }
            //    if (groups != null && groups.RTDict[1] < 0.5)
            //    {
            //        Bad++;
            //        badIndex.Add(i, groups.RTDict[1]);
            //    }
            //}

            //var xicGroups = XICGrouping(peakOne, peakTwo, 818.3857);

            DrawPlot(Ms1Peaks, 387.2589);

            string outPath = filePath.Replace(".raw", "_first_shiftcase.mzML").ToString();
            string outPath_2 = filePath.Replace(".raw", "_second_shiftcase.mzML").ToString();
            string outPath_3 = filePath.Replace(".raw", "_third_shiftcase.mzML").ToString();
            string outPath_4 = filePath.Replace(".raw", "_fourth_shiftcase.mzML").ToString();
            string outPath_5 = filePath.Replace(".raw", "_fifth_shiftcase.mzML").ToString();

            SourceFile sourceFile = new SourceFile(reader.SourceFile.NativeIdFormat,
                reader.SourceFile.MassSpectrometerFileFormat, reader.SourceFile.CheckSum, reader.SourceFile.FileChecksumType, reader.SourceFile.Uri, reader.SourceFile.Id, reader.SourceFile.FileName);


            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(File1.ToArray(), sourceFile), outPath, false);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(File2.ToArray(), sourceFile), outPath_2, false);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(File3.ToArray(), sourceFile), outPath_3, false);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(File4.ToArray(), sourceFile), outPath_4, false);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(File5.ToArray(), sourceFile), outPath_5, false);


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