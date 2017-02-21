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
        public void setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
        }

        [Test]
        public void LoadThermoTest()
        {
            ThermoDynamicData a = ThermoDynamicData.InitiateDynamicConnection(@"Shew_246a_LCQa_15Oct04_Andro_0904-2_4-20.RAW");
            Assert.AreEqual(3316, a.NumSpectra);
            Assert.AreEqual(3316, a.NumSpectra);

            var scan = a.GetOneBasedScan(53);
            Assert.AreEqual(1.2623333333333333, scan.RetentionTime);
            Assert.AreEqual(1, scan.MsnOrder);
            Assert.AreEqual("+ c ESI Full ms [400.00-2000.00]", scan.ScanFilter);

            var spectrum = a.GetOneBasedScan(53).MassSpectrum;

            var peak = spectrum.PeakWithHighestY;
            Assert.AreEqual(75501, peak.Intensity);

            Assert.AreEqual(1, spectrum.FilterByY(7.5e4, double.MaxValue).Count());
            Assert.AreEqual(2, spectrum.Extract(new DoubleRange(923, 928)).Count());

            Assert.AreEqual(double.NaN, spectrum.GetSignalToNoise(1));

            Assert.AreEqual("1.3", a.thermoGlobalParams.pbstrInstSoftwareVersion);
            var ms2scan = a.GetOneBasedScan(948) as IMsDataScanWithPrecursor<ThermoSpectrum>;
            Assert.IsNull(ms2scan.SelectedIonGuessChargeStateGuess);
            var precursorScan = a.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber);
            ms2scan.RecomputeChargeState(precursorScan.MassSpectrum.Extract(ms2scan.IsolationMz - 2.1, ms2scan.IsolationMz + 2.1).ToList(), 0.1, 4);
            Assert.AreEqual(1, ms2scan.SelectedIonGuessChargeStateGuess);
            Assert.IsNull(ms2scan.SelectedIonGuessIntensity);
            ms2scan.RecomputeSelectedPeak(precursorScan.MassSpectrum);
            Assert.AreEqual(4125760, ms2scan.SelectedIonGuessIntensity);

            Assert.AreEqual("LCQ", a.thermoGlobalParams.pbstrInstName);
            Assert.AreEqual("LCQ", a.thermoGlobalParams.pbstrInstModel);

            Assert.AreEqual(false, a.thermoGlobalParams.MonoisotopicselectionEnabled);
        }

        [Test]
        public void ThermoLoadError()
        {
            Assert.Throws<ThermoReadException>(() => ThermoStaticData.LoadAllStaticData(@"aaa.RAW"));
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

            //Assert.AreEqual(0.98, a.GetElapsedScanTime(100), 0.01);

            //var cromatogram = a.GetTicChroma();

            //Assert.AreEqual(360, cromatogram.Size);
            //Assert.AreEqual(0.01, cromatogram.FirstTime, 0.002);
            //Assert.AreEqual(2.788433333, cromatogram.PeakWithHighestY.Time, 0.0001);
            //Assert.AreEqual(2.788433333, cromatogram.GetApex(0, 5).Time, 0.0001);

            var newSpectrum = new ThermoSpectrum(a.GetOneBasedScan(51).MassSpectrum);
            Assert.AreEqual(22246 / 5574.8, newSpectrum.GetSignalToNoise(1), 0.01);

            Assert.AreEqual(1, newSpectrum.GetCharges()[1]);
            Assert.AreEqual(102604, newSpectrum.GetResolutions()[1]);

            //Assert.AreEqual(181, newSpectrum.NewSpectrumExtract(500, 1000).Size);

            //Assert.AreEqual(0, newSpectrum.NewSpectrumExtract(-3, -2).Size);

            //var hm = newSpectrum.NewSpectrumExtract(501, 502);

            //Assert.AreEqual(0, hm.Size);

            Assert.AreEqual(1120, a.GetOneBasedScan(1).MassSpectrum.Size);

            //var b = new ThermoStaticData.LoadAllStaticData(@"05-13-16_cali_MS_60K-res_MS.raw", 400);

            //Assert.AreEqual(400, b.GetOneBasedScan(1).MassSpectrum.Size);

            //Assert.AreEqual(0, b.Where(eb => eb.MsnOrder > 1).Count());

            //Assert.AreEqual(false, b.MonoisotopicPrecursorSelectionEnabled);

            //IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> uu = b;

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

            Assert.IsFalse(a.thermoGlobalParams.MonoisotopicselectionEnabled);
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
            Assert.AreEqual(4, s2.First().Mz);
        }

        #endregion Public Methods

    }
}