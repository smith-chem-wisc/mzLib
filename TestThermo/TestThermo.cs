using Chemistry;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.IO;
using System.Linq;

namespace TestThermo
{
    [TestFixture]
    public sealed class TestThermo
    {
        #region Public Methods

        [Test]
        public static void ReadWriteReadEtc()
        {
            {
                ThermoStaticData a = ThermoStaticData.LoadAllStaticData(@"testFileWMS2.raw");

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(a, "a.mzML", false);

                var aa = Mzml.LoadAllStaticData("a.mzML");

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(aa, "aa.mzML", true);

                Mzml.LoadAllStaticData("aa.mzML");
            }
            {
                ThermoStaticData a = ThermoStaticData.LoadAllStaticData(@"small.raw");

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(a, "a.mzML", false);

                var aa = Mzml.LoadAllStaticData("a.mzML");

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(aa, "aa.mzML", true);

                Mzml.LoadAllStaticData("aa.mzML");
            }
            {
                ThermoStaticData a = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw");

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(a, "a.mzML", false);

                var aa = Mzml.LoadAllStaticData("a.mzML");

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(aa, "aa.mzML", true);

                Mzml.LoadAllStaticData("aa.mzML");
            }
        }

        [Test]
        public static void ThermoLoadError()
        {
            Assert.Throws<MzLibException>(() => ThermoStaticData.LoadAllStaticData(@"aaa.RAW"));
        }

        [Test]
        public static void ThermoReaderNotInstalled()
        {
            bool check = ThermoFile.CheckForMsFileReader();
            Assert.IsTrue(check);
        }

        [Test]
        public static void LoadCompressedMzml()
        {
            ThermoStaticData a = ThermoStaticData.LoadAllStaticData(@"small.RAW");

            Mzml b = Mzml.LoadAllStaticData(@"smallCentroid.mzML");

            Assert.AreEqual(a.NumSpectra, b.NumSpectra);

            Assert.AreEqual(a.GetOneBasedScan(1).MassSpectrum.XofPeakWithHighestY, b.GetOneBasedScan(1).MassSpectrum.XofPeakWithHighestY, 1e-8);
            Assert.IsTrue(Math.Abs((a.GetOneBasedScan(1).MassSpectrum.YofPeakWithHighestY - b.GetOneBasedScan(1).MassSpectrum.YofPeakWithHighestY) / b.GetOneBasedScan(1).MassSpectrum.YofPeakWithHighestY) < 1e-8);

            Assert.AreEqual(a.GetOneBasedScan(2).MassSpectrum.XofPeakWithHighestY, b.GetOneBasedScan(2).MassSpectrum.XofPeakWithHighestY, 1e-8);
            Assert.IsTrue(Math.Abs((a.GetOneBasedScan(2).MassSpectrum.YofPeakWithHighestY - b.GetOneBasedScan(2).MassSpectrum.YofPeakWithHighestY) / b.GetOneBasedScan(1).MassSpectrum.YofPeakWithHighestY) < 1e-8);

            Assert.AreEqual(a.GetOneBasedScan(3).MassSpectrum.XofPeakWithHighestY, b.GetOneBasedScan(3).MassSpectrum.XofPeakWithHighestY, 1e-8);
            Assert.IsTrue(Math.Abs((a.GetOneBasedScan(3).MassSpectrum.YofPeakWithHighestY - b.GetOneBasedScan(3).MassSpectrum.YofPeakWithHighestY) / b.GetOneBasedScan(1).MassSpectrum.YofPeakWithHighestY) < 1e-8);
        }

        [Test]
        public static void LoadThermoTest2()
        {
            ThermoStaticData a = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw");
            Assert.AreEqual(360, a.NumSpectra);
            Assert.GreaterOrEqual(1000, a.GetOneBasedScan(1).MassSpectrum.Extract(0, 500).Last().X);
            Assert.AreEqual(2, a.GetOneBasedScan(1).MassSpectrum.FilterByY(5e6, double.MaxValue).Count());
            var ye = a.GetOneBasedScan(1).MassSpectrum.CopyTo2DArray();
            Assert.AreEqual(77561752, a.GetOneBasedScan(1).TotalIonCurrent);
            Assert.AreEqual(144, a.GetClosestOneBasedSpectrumNumber(2));

            var newSpectrum = new ThermoSpectrum(a.GetOneBasedScan(51).MassSpectrum);

            Assert.AreEqual(1120, a.GetOneBasedScan(1).MassSpectrum.Size);

            var newDeconvolution = a.GetOneBasedScan(1).MassSpectrum.Deconvolute(new MzRange(double.MinValue, double.MaxValue), 1, 10, 1, 4).ToList();

            Assert.IsTrue(newDeconvolution.Any(b => Math.Abs(b.peaks.First().Item1.ToMass(b.charge) - 523.257) < 0.001));

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(a, Path.Combine(TestContext.CurrentContext.TestDirectory, "convertedThermo.mzML"), false);

            var sdafaf = a.Deconvolute(null, null, 1, 30, 10, 3, 10, b => true).OrderByDescending(b => b.NumPeaks).First();

            Assert.IsTrue(Math.Abs(262.64 - sdafaf.Mass.ToMz(2)) <= 0.01);

            using (ThermoDynamicData dynamicThermo = ThermoDynamicData.InitiateDynamicConnection(@"05-13-16_cali_MS_60K-res_MS.raw"))
            {
                Assert.AreEqual(136, dynamicThermo.GetClosestOneBasedSpectrumNumber(1.89));
                dynamicThermo.ClearCachedScans();
            }

            Mzml readCovertedMzmlFile = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "convertedThermo.mzML"));

            Assert.AreEqual(a.First().Polarity, readCovertedMzmlFile.First().Polarity);
        }

        [Test]
        public static void WindowFilteringStaticTest()
        {
            //test window number of 1
            ThermoStaticData a_w = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw");
            ThermoStaticData b_w = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(top: 400, num:1), trimMs1Peaks: true);
            ThermoStaticData c_w = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(ratio: 0.001, num: 1), trimMs1Peaks: true);
            ThermoStaticData d_w = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(ratio: 0.001, top: 400, num: 1), trimMs1Peaks: true);

            var aLen = a_w.GetOneBasedScan(1).MassSpectrum.Size;
            var bLen = b_w.GetOneBasedScan(1).MassSpectrum.Size;
            var cLen = c_w.GetOneBasedScan(1).MassSpectrum.Size;
            var dLen = d_w.GetOneBasedScan(1).MassSpectrum.Size;

            Assert.AreEqual(Math.Min(bLen, cLen), dLen);

            var aLen2 = a_w.GetOneBasedScan(2).MassSpectrum.Size;
            var bLen2 = b_w.GetOneBasedScan(2).MassSpectrum.Size;
            var cLen2 = c_w.GetOneBasedScan(2).MassSpectrum.Size;
            var dLen2 = d_w.GetOneBasedScan(2).MassSpectrum.Size;

            Assert.AreEqual(Math.Min(bLen2, cLen2), dLen2);

            var aLen3 = a_w.GetOneBasedScan(3).MassSpectrum.Size;
            var bLen3 = b_w.GetOneBasedScan(3).MassSpectrum.Size;
            var cLen3 = c_w.GetOneBasedScan(3).MassSpectrum.Size;
            var dLen3 = d_w.GetOneBasedScan(3).MassSpectrum.Size;

            Assert.AreEqual(Math.Min(bLen3, cLen3), dLen3);
        }

        [Test]
        public static void MultiWindowFiltering()
        {
            //tests for filtering with window
            ThermoStaticData a_w = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(num: 1), trimMs1Peaks: true);
            Assert.AreEqual(1120, a_w.GetOneBasedScan(1).MassSpectrum.Size);
            //number of 2
            ThermoStaticData b_w = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(top: 200,num: 3), trimMs1Peaks: true);
            Assert.AreEqual(600, b_w.GetOneBasedScan(1).MassSpectrum.Size);
            //number of 4
            ThermoStaticData c_w = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(top: 200, num: 4), trimMs1Peaks: true);
            Assert.AreEqual(800, c_w.GetOneBasedScan(1).MassSpectrum.Size);
            //number of 6, which doesn't divide 1120
            ThermoStaticData d_w = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(top: 150, num: 6), trimMs1Peaks: true);
            Assert.AreEqual(900, d_w.GetOneBasedScan(1).MassSpectrum.Size);
        }

        [Test]
        public static void LoadThermoFiltered()
        {
            ThermoStaticData a = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw");
            ThermoStaticData b = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(top: 400), trimMs1Peaks: true);
            ThermoStaticData c = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(ratio: 0.001), trimMs1Peaks: true);
            ThermoStaticData d = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(ratio: 0.001, top: 400), trimMs1Peaks: true);

            var aLen = a.GetOneBasedScan(1).MassSpectrum.Size;
            var bLen = b.GetOneBasedScan(1).MassSpectrum.Size;
            var cLen = c.GetOneBasedScan(1).MassSpectrum.Size;
            var dLen = d.GetOneBasedScan(1).MassSpectrum.Size;

            Assert.AreEqual(Math.Min(bLen, cLen), dLen);
        }

        [Test]
        public static void LoadThermoFiltered2()
        {
            ThermoStaticData a = ThermoStaticData.LoadAllStaticData(@"small.raw");
            ThermoStaticData b = ThermoStaticData.LoadAllStaticData(@"small.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(top: 40), trimMs1Peaks: true, trimMsMsPeaks: true);
            ThermoStaticData c = ThermoStaticData.LoadAllStaticData(@"small.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(ratio: 0.1), trimMs1Peaks: true, trimMsMsPeaks: true);
            ThermoStaticData d = ThermoStaticData.LoadAllStaticData(@"small.raw", thermoParams: new MsDataFile<IThermoScan>.FilteringParams(ratio: 0.1, top: 40), trimMs1Peaks: true, trimMsMsPeaks: true);

            var aLen = a.GetOneBasedScan(1).MassSpectrum.Size;
            var bLen = b.GetOneBasedScan(1).MassSpectrum.Size;
            var cLen = c.GetOneBasedScan(1).MassSpectrum.Size;
            var dLen = d.GetOneBasedScan(1).MassSpectrum.Size;

            Assert.AreEqual(Math.Min(bLen, cLen), dLen);

            var aLen2 = a.GetOneBasedScan(2).MassSpectrum.Size;
            var bLen2 = b.GetOneBasedScan(2).MassSpectrum.Size;
            var cLen2 = c.GetOneBasedScan(2).MassSpectrum.Size;
            var dLen2 = d.GetOneBasedScan(2).MassSpectrum.Size;

            Assert.AreEqual(Math.Min(bLen2, cLen2), dLen2);

            var aLen3 = a.GetOneBasedScan(3).MassSpectrum.Size;
            var bLen3 = b.GetOneBasedScan(3).MassSpectrum.Size;
            var cLen3 = c.GetOneBasedScan(3).MassSpectrum.Size;
            var dLen3 = d.GetOneBasedScan(3).MassSpectrum.Size;

            Assert.AreEqual(Math.Min(bLen3, cLen3), dLen3);
        }

        [Test]
        public static void LoadThermoTest3()
        {
            ThermoStaticData a = ThermoStaticData.LoadAllStaticData(@"small.RAW");

            Assert.IsTrue(a.Where(eb => eb.MsnOrder > 1).Count() > 0);

            Assert.IsTrue(a.Where(eb => eb.MsnOrder == 1).Count() > 0);

            Assert.IsFalse(a.ThermoGlobalParams.MonoisotopicselectionEnabled);

            var hehe = a.First(b => b.MsnOrder > 1) as ThermoScanWithPrecursor;

            var prec = a.GetOneBasedScan(hehe.OneBasedPrecursorScanNumber.Value);

            Assert.IsNull(hehe.SelectedIonChargeStateGuess);

            Assert.IsNull(hehe.SelectedIonIntensity);

            hehe.ComputeSelectedPeakIntensity(prec.MassSpectrum);

            Assert.AreEqual(1017759, hehe.SelectedIonIntensity, 1);

            Assert.IsNull(hehe.SelectedIonMonoisotopicGuessIntensity);

            hehe.ComputeMonoisotopicPeakIntensity(prec.MassSpectrum);

            Assert.AreEqual(1017759, hehe.SelectedIonMonoisotopicGuessIntensity, 1);
        }

        [Test]
        public static void ThermoSpectrumTest()
        {
            double[] mz = new double[] { 1 };
            double[] intensity = new double[] { 1 };
            ThermoSpectrum s1 = new ThermoSpectrum(mz, intensity, false);
            ThermoSpectrum s2 = new ThermoSpectrum(mz, intensity, false);
            s1.ReplaceXbyApplyingFunction((a) => 4);
            Assert.AreEqual(4, s2.XArray[0]);
        }

        [Test]
        public static void ThermoDynamicTest()
        {
            ThermoDynamicData dynamicThermo = ThermoDynamicData.InitiateDynamicConnection(@"testFileWMS2.raw");
            var ms1scan = dynamicThermo.GetOneBasedScan(1);
            ThermoScanWithPrecursor ms2scan = dynamicThermo.GetOneBasedScan(651) as ThermoScanWithPrecursor;
            Assert.That(ms1scan.OneBasedScanNumber == 1);
            Assert.That(ms2scan.OneBasedScanNumber == 651);
            Assert.That(Math.Round(ms2scan.RetentionTime, 2) == 12.16);
            Assert.That(ms2scan.OneBasedPrecursorScanNumber == 650);
            Assert.That(ms2scan.SelectedIonMZ == 442.67);
            var t = dynamicThermo.ThermoGlobalParams.msOrderByScan;
            Assert.That(t[0] == 1);
            Assert.That(t[5] == 1);
            Assert.That(t[649] == 1);
            Assert.That(t[650] == 2);
            Assert.That(!t.Where(v => v == 0).Any());
        }

        [Test]
        public static void TestSummedMsDataFile()
        {
            ThermoStaticData rawFile = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw");

            // 3 scans

            SummedMsDataFile summed3 = new SummedMsDataFile(rawFile, 3, 10);

            Assert.AreEqual(rawFile.NumSpectra - 2, summed3.NumSpectra);

            var resultingTic = summed3.GetOneBasedScan(1).TotalIonCurrent;
            var mySummedTic = rawFile.GetOneBasedScan(1).MassSpectrum.SumOfAllY + rawFile.GetOneBasedScan(2).MassSpectrum.SumOfAllY + rawFile.GetOneBasedScan(3).MassSpectrum.SumOfAllY;
            var instrumentSummedTic = rawFile.GetOneBasedScan(1).TotalIonCurrent + rawFile.GetOneBasedScan(2).TotalIonCurrent + rawFile.GetOneBasedScan(3).TotalIonCurrent;

            // Tics are approximately what they should be
            Assert.IsTrue(Math.Abs(resultingTic - mySummedTic) / mySummedTic < 1e-4);
            Assert.IsTrue(Math.Abs(resultingTic - instrumentSummedTic) / instrumentSummedTic < 1e-1);

            // Equal to representative
            Assert.AreEqual(summed3.GetOneBasedScan(1).RetentionTime, rawFile.GetOneBasedScan(2).RetentionTime);

            Assert.IsTrue(summed3.GetOneBasedScan(1).MassSpectrum.Size <= rawFile.GetOneBasedScan(1).MassSpectrum.Size + rawFile.GetOneBasedScan(2).MassSpectrum.Size + rawFile.GetOneBasedScan(3).MassSpectrum.Size);
            Assert.IsTrue(summed3.GetOneBasedScan(1).MassSpectrum.Size >= rawFile.GetOneBasedScan(1).MassSpectrum.Size);
            Assert.IsTrue(summed3.GetOneBasedScan(1).MassSpectrum.Size >= rawFile.GetOneBasedScan(2).MassSpectrum.Size);
            Assert.IsTrue(summed3.GetOneBasedScan(1).MassSpectrum.Size >= rawFile.GetOneBasedScan(3).MassSpectrum.Size);

            Assert.IsTrue(summed3.GetOneBasedScan(1).MassSpectrum.YofPeakWithHighestY == rawFile.GetOneBasedScan(1).MassSpectrum.YofPeakWithHighestY + rawFile.GetOneBasedScan(2).MassSpectrum.YofPeakWithHighestY + rawFile.GetOneBasedScan(3).MassSpectrum.YofPeakWithHighestY);

            // Interval of 893-899 mz

            Assert.AreEqual(2, rawFile.GetOneBasedScan(1).MassSpectrum.NumPeaksWithinRange(893, 899));
            Assert.AreEqual(2, rawFile.GetOneBasedScan(2).MassSpectrum.NumPeaksWithinRange(893, 899));
            Assert.AreEqual(1, rawFile.GetOneBasedScan(3).MassSpectrum.NumPeaksWithinRange(893, 899));

            // One peak persists across the three scans! So instead of 5 see three peaks in summed
            Assert.AreEqual(3, summed3.GetOneBasedScan(1).MassSpectrum.NumPeaksWithinRange(893, 899));

            Assert.AreEqual(summed3.GetOneBasedScan(1).MassSpectrum.FirstX, Math.Min(Math.Min(rawFile.GetOneBasedScan(1).MassSpectrum.FirstX, rawFile.GetOneBasedScan(2).MassSpectrum.FirstX), rawFile.GetOneBasedScan(3).MassSpectrum.FirstX));

            Assert.AreEqual(summed3.GetOneBasedScan(1).MassSpectrum.LastX, Math.Max(Math.Max(rawFile.GetOneBasedScan(1).MassSpectrum.LastX, rawFile.GetOneBasedScan(2).MassSpectrum.LastX), rawFile.GetOneBasedScan(3).MassSpectrum.LastX));

            // 5 scans
            SummedMsDataFile summed5 = new SummedMsDataFile(rawFile, 5, 10);

            Assert.AreEqual(rawFile.NumSpectra - 4, summed5.NumSpectra);

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(summed5, "testSummed.mzML", false);

            var ok = Mzml.LoadAllStaticData("testSummed.mzML");

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(ok, "testSummed2.mzML", false);

            Mzml.LoadAllStaticData("testSummed2.mzML");
        }

        [Test]
        public static void WriteIndexedMzmlFromThermoTest()
        {
            var smallThermo = ThermoStaticData.LoadAllStaticData(@"small.raw");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(smallThermo, Path.Combine(TestContext.CurrentContext.TestDirectory, "Hi.mzML"), true);
            var smallMzml = Mzml.LoadAllStaticData(@"Hi.mzML");
            Assert.AreEqual(smallMzml.NumSpectra, 48);
            Assert.AreEqual(smallMzml.GetOneBasedScan(8).OneBasedScanNumber, 8);
            Assert.AreEqual(smallThermo.GetOneBasedScan(5).RetentionTime, smallMzml.GetOneBasedScan(5).RetentionTime);
        }

        [OneTimeSetUp]
        public void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
        }

        #endregion Public Methods
    }
}