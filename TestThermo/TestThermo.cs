using Chemistry;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Linq;

namespace TestThermo
{
    [TestFixture]
    public sealed class TestThermo
    {

        #region Public Methods

        [OneTimeSetUp]
        public void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
        }

        [Test]
        public void ThermoLoadError()
        {
            Assert.Throws<MzLibException>(() => ThermoStaticData.LoadAllStaticData(@"aaa.RAW"));
        }

        [Test]
        public void LoadThermoTest2()
        {
            ThermoStaticData a = ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw");
            Assert.AreEqual(360, a.NumSpectra);
            var ok = a.GetOneBasedScan(1).MassSpectrum.GetNoises();
            Assert.AreEqual(2401.57, ok[0], 0.01);
            Assert.GreaterOrEqual(1000, a.GetOneBasedScan(1).MassSpectrum.Extract(0, 500).Last().X);
            Assert.AreEqual(2, a.GetOneBasedScan(1).MassSpectrum.FilterByY(5e6, double.MaxValue).Count());
            var ye = a.GetOneBasedScan(1).MassSpectrum.CopyTo2DArray();
            Assert.AreEqual(1, ye[4, 1119]);
            Assert.AreEqual("(195.0874,1.021401E+07) z = +1 SN = 4170.38", a.GetOneBasedScan(1).MassSpectrum.PeakWithHighestY.ToString());
            Assert.AreEqual(77561752, a.GetOneBasedScan(1).TotalIonCurrent);
            Assert.AreEqual(144, a.GetClosestOneBasedSpectrumNumber(2));

            var newSpectrum = new ThermoSpectrum(a.GetOneBasedScan(51).MassSpectrum);
            Assert.AreEqual(22246 / 5574.8, newSpectrum.GetSignalToNoise(1), 0.01);

            Assert.AreEqual(1, newSpectrum.GetCharges()[1]);
            Assert.AreEqual(102604, newSpectrum.GetResolutions()[1]);

            Assert.AreEqual(1120, a.GetOneBasedScan(1).MassSpectrum.Size);

            var newDeconvolution = a.GetOneBasedScan(1).MassSpectrum.Deconvolute(new MzRange(double.MinValue, double.MaxValue), 10, 1, 4, b => true).ToList();

            Assert.IsTrue(newDeconvolution.Any(b => Math.Abs(b.peaks.First().Mz.ToMass(b.charge) - 523.257) < 0.001));

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(a, "convertedThermo.mzML", false);

            var sdafaf = a.Deconvolute(null, null, 30, 10, 3, b => true, 10, b => true).OrderByDescending(b => b.NumPeaks).First();

            Assert.IsTrue(Math.Abs(262.64 - sdafaf.Mass.ToMz(2)) <= 0.01);

            using (ThermoDynamicData dynamicThermo = ThermoDynamicData.InitiateDynamicConnection(@"05-13-16_cali_MS_60K-res_MS.raw"))
            {
                Assert.AreEqual(136, dynamicThermo.GetClosestOneBasedSpectrumNumber(1.89));
                dynamicThermo.ClearCachedScans();
            }
        }

        [Test]
        public void LoadThermoTest3()
        {
            ThermoStaticData a = ThermoStaticData.LoadAllStaticData(@"small.RAW");

            Assert.IsTrue(a.Where(eb => eb.MsnOrder > 1).Count() > 0);

            Assert.IsTrue(a.Where(eb => eb.MsnOrder == 1).Count() > 0);

            Assert.IsFalse(a.ThermoGlobalParams.MonoisotopicselectionEnabled);

            var hehe = a.First(b => b.MsnOrder > 1) as ThermoScanWithPrecursor;

            var prec = a.GetOneBasedScan(hehe.OneBasedPrecursorScanNumber);

            Assert.IsNull(hehe.SelectedIonChargeStateGuess);

            Assert.IsNull(hehe.SelectedIonIntensity);

            hehe.ComputeSelectedPeakIntensity(prec.MassSpectrum);

            Assert.AreEqual(1017759, hehe.SelectedIonIntensity, 1);

            Assert.IsNull(hehe.SelectedIonMonoisotopicGuessIntensity);

            hehe.ComputeMonoisotopicPeakIntensity(prec.MassSpectrum);

            Assert.AreEqual(1017759, hehe.SelectedIonMonoisotopicGuessIntensity, 1);
        }

        [Test]
        public void ThermoSpectrumTest()
        {
            double[] resolutions = new double[] { 1 };
            int[] charge = new int[] { 1 };
            double[] mz = new double[] { 1 };
            double[] intensity = new double[] { 1 };
            double[] noise = new double[] { 1 };
            ThermoSpectrum s1 = new ThermoSpectrum(mz, intensity, noise, charge, resolutions, false);
            ThermoSpectrum s2 = new ThermoSpectrum(mz, intensity, noise, charge, resolutions, false);
            s1.ReplaceXbyApplyingFunction((a) => 4);
            Assert.AreEqual(4, s2[0].Mz);
        }

        [Test]
        public void ThermoDynamicTest()
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
        public void TestSummedMsDataFile()
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

            Assert.IsTrue(summed3.GetOneBasedScan(1).MassSpectrum.PeakWithHighestY.Intensity == rawFile.GetOneBasedScan(1).MassSpectrum.PeakWithHighestY.Intensity + rawFile.GetOneBasedScan(2).MassSpectrum.PeakWithHighestY.Intensity + rawFile.GetOneBasedScan(3).MassSpectrum.PeakWithHighestY.Intensity);

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
        }

        #endregion Public Methods

    }
}