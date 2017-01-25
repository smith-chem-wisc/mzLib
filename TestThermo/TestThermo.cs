using IO.Thermo;
using NUnit.Framework;
using Spectra;
using System;
using System.IO;
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
            ThermoRawFile a = new ThermoRawFile(@"Shew_246a_LCQa_15Oct04_Andro_0904-2_4-20.RAW");
            a.Open();
            a.Open();
            Assert.AreEqual(3316, a.NumSpectra);
            Assert.AreEqual(3316, a.NumSpectra);

            var scan = a.GetOneBasedScan(53);
            Assert.AreEqual(1.2623333333333333, scan.RetentionTime);
            Assert.AreEqual(1, scan.MsnOrder);
            Assert.AreEqual("controllerType=0 controllerNumber=1 scan=53", scan.Id);
            Assert.AreEqual("+ c ESI Full ms [400.00-2000.00]", scan.ScanFilter);

            var spectrum = a.GetOneBasedScan(53).MassSpectrum;

            var peak = spectrum.PeakWithHighestY;
            Assert.AreEqual(75501, peak.Intensity);

			Assert.AreEqual(1, spectrum.NewSpectrumFilterByY(7.5e4, double.MaxValue).Count);
            Assert.AreEqual(2, spectrum.NewSpectrumExtract(new DoubleRange(923, 928)).Count);

            Assert.AreEqual(double.NaN, spectrum.GetSignalToNoise(1));

            Assert.AreEqual("1.3", a.GetSofwareVersion());
            double ya;
            a.GetOneBasedScan(948).TryGetSelectedIonGuessIntensity(out ya);
            Assert.AreEqual(4125760, ya);

            Assert.AreEqual("LCQ", a.GetInstrumentName());
            Assert.AreEqual("LCQ", a.GetInstrumentModel());

            Assert.AreEqual(0, a.GetMsxPrecursors(1289).Count);
            Assert.AreEqual(1, a.GetMsxPrecursors(1290).Count);
            Assert.AreEqual(1194.53, a.GetMsxPrecursors(1290).First());

			Assert.AreEqual(false, a.MonoisotopicPrecursorSelectionEnabled);
        }

        [Test]
        public void ThermoLoadError()
        {
            ThermoRawFile a = new ThermoRawFile(@"aaa.RAW");
            Assert.Throws<IOException>(() =>  a.Open());
        }

        [Test]
        public void LoadThermoTest2()
        {
            ThermoRawFile a = new ThermoRawFile(@"05-13-16_cali_MS_60K-res_MS.raw");
            a.Open();
            Assert.AreEqual(360, a.NumSpectra);
            var ok = a.GetOneBasedScan(1).MassSpectrum.GetNoises();
            Assert.AreEqual(2401.57, ok[0], 0.01);
            ThermoSpectrum ok2 = a.GetOneBasedScan(1).MassSpectrum.NewSpectrumExtract(0, 500);
            Assert.GreaterOrEqual(1000, a.GetOneBasedScan(1).MassSpectrum.NewSpectrumExtract(0, 500).LastX);
			Assert.AreEqual(2, a.GetOneBasedScan(1).MassSpectrum.NewSpectrumFilterByY(5e6, double.MaxValue).Count);
            var ye = a.GetOneBasedScan(1).MassSpectrum.CopyTo2DArray();
            Assert.AreEqual(1, ye[4, 1119]);
            Assert.AreEqual("(195.0874,1.0214E+07) z = +1 SN = 4170.38", a.GetOneBasedScan(1).MassSpectrum.PeakWithHighestY.ToString());
            Assert.AreEqual(77561752, a.GetTic(1));
            Assert.AreEqual(144, a.GetClosestOneBasedSpectrumNumber(2));

            Assert.AreEqual(0.98, a.GetElapsedScanTime(100), 0.01);

            var cromatogram = a.GetTicChroma();

            Assert.AreEqual(360, cromatogram.Count);
            Assert.AreEqual(0.01, cromatogram.FirstTime, 0.002);
            Assert.AreEqual(2.788433333, cromatogram.PeakWithHighestY.Time, 0.0001);
            Assert.AreEqual(2.788433333, cromatogram.GetApex(0, 5).Time, 0.0001);

            var newSpectrum = new ThermoSpectrum(a.GetOneBasedScan(51).MassSpectrum);
            Assert.AreEqual(22246 / 5574.8, newSpectrum.GetSignalToNoise(1), 0.01);

            Assert.AreEqual(1, newSpectrum.GetCharges()[1]);
            Assert.AreEqual(102604, newSpectrum.GetResolutions()[1]);

            Assert.AreEqual(181, newSpectrum.NewSpectrumExtract(500, 1000).Count);

            Assert.AreEqual(0, newSpectrum.NewSpectrumExtract(-3, -2).Count);

            var hm = newSpectrum.NewSpectrumExtract(501, 502);

            Assert.AreEqual(0, hm.Count);

            Assert.AreEqual(1120, a.GetOneBasedScan(1).MassSpectrum.Count);

            a.Close();

            var b = new ThermoRawFile(@"05-13-16_cali_MS_60K-res_MS.raw", 400);
            b.Open();

            Assert.AreEqual(400, b.GetOneBasedScan(1).MassSpectrum.Count);


			Assert.AreEqual(0, b.Where(eb=>eb.MsnOrder>1).Count());

			Assert.AreEqual(false, b.MonoisotopicPrecursorSelectionEnabled);

        }

        [Test]
		public void LoadThermoTest3()
		{
			ThermoRawFile a = new ThermoRawFile(@"small.RAW");
			a.Open();

			Assert.AreEqual(0, a.Where(eb => eb.MsnOrder > 1).Count());

			Assert.AreEqual(0, a.Where(eb => eb.MsnOrder == 1).Count());

			a.Close();

		}

		#endregion Public Methods

	}
}