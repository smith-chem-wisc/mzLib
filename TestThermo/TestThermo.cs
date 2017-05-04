using Chemistry;
using IO.MzML;
using IO.Thermo;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
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

            var newSpectrum = new ThermoSpectrum(a.GetOneBasedScan(51).MassSpectrum);
            Assert.AreEqual(22246 / 5574.8, newSpectrum.GetSignalToNoise(1), 0.01);

            Assert.AreEqual(1, newSpectrum.GetCharges()[1]);
            Assert.AreEqual(102604, newSpectrum.GetResolutions()[1]);

            Assert.AreEqual(1120, a.GetOneBasedScan(1).MassSpectrum.Size);

            List<IO.Thermo.Deconvolution.PossibleProteoform> cool = a.GetOneBasedScan(1).MassSpectrum.SpecialThermoDeconvolution(0.1).ToList();

            Assert.AreEqual(523.257, cool[0].GetMonoisotopicMass(), 0.001);

            var newDeconvolution = a.GetOneBasedScan(1).MassSpectrum.Deconvolute(new MzRange(double.MinValue, double.MaxValue), 10, new Tolerance("1 PPM"), 4).ToList();

            Assert.IsTrue(newDeconvolution.Any(b => Math.Abs(b.Item1.First().Mz.ToMass(b.Item2) - 523.257) < 0.001));

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(a, "convertedThermo.mzML", false);

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
            Assert.AreEqual(4, s2.First().Mz);
        }

        #endregion Public Methods

    }
}