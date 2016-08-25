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
            Assert.AreEqual(1, a.FirstSpectrumNumber);
            Assert.AreEqual(3316, a.LastSpectrumNumber);
            Assert.AreEqual(3316, a.LastSpectrumNumber);

            var scan = a.GetScan(53);
            Assert.AreEqual(1.2623333333333333, scan.RetentionTime);
            Assert.AreEqual(1, scan.MsnOrder);
            Assert.AreEqual("controllerType=0 controllerNumber=1 scan=53", scan.id);
            Assert.AreEqual("+ c ESI Full ms [400.00-2000.00]", scan.ScanFilter);


            var spectrum = a.GetScan(53).MassSpectrum;

            var peak = spectrum.PeakWithHighestY;
            Assert.AreEqual(75501, peak.Intensity);

            Assert.AreEqual(1, spectrum.newSpectrumFilterByY(7.5e4).Count);
            Assert.AreEqual(2, spectrum.newSpectrumExtract(new DoubleRange(923, 928)).Count);


            Assert.AreEqual(double.NaN, spectrum.GetSignalToNoise(1));

            Assert.AreEqual("1.3", a.GetSofwareVersion());
            double ya;
            a.GetScan(948).TryGetSelectedIonGuessIntensity(out ya);
            Assert.AreEqual(4125760, ya);

            Assert.AreEqual("LCQ", a.GetInstrumentName());
            Assert.AreEqual("LCQ", a.GetInstrumentModel());

            Assert.AreEqual(0, a.GetMSXPrecursors(1289).Count);
            Assert.AreEqual(1, a.GetMSXPrecursors(1290).Count);
            Assert.AreEqual(1194.53, a.GetMSXPrecursors(1290).First());

        }


        [Test]
        public void ThermoLoadError()
        {
            ThermoRawFile a = new ThermoRawFile(@"aaa.RAW");
            Assert.Throws<IOException>(() => { a.Open(); });
        }
        [Test]
        public void LoadThermoTest2()
        {
            ThermoRawFile a = new ThermoRawFile(@"05-13-16_cali_MS_60K-res_MS.raw");
            a.Open();
            Assert.AreEqual(360, a.LastSpectrumNumber);
            var ok = a.GetScan(1).MassSpectrum.GetNoises();
            Assert.AreEqual(2401.57, ok[0], 0.01);
            ThermoSpectrum ok2 = a.GetScan(1).MassSpectrum.newSpectrumExtract(0, 500);
            Assert.GreaterOrEqual(1000, a.GetScan(1).MassSpectrum.newSpectrumExtract(0, 500).LastX);
            Assert.AreEqual(2, a.GetScan(1).MassSpectrum.newSpectrumFilterByY(5e6).Count);
            var ye = a.GetScan(1).MassSpectrum.CopyTo2DArray();
            Assert.AreEqual(1, ye[4, 1119]);
            Assert.AreEqual("(195.0874,1.0214E+07) z = +1 SN = 4170.38", a.GetScan(1).MassSpectrum.PeakWithHighestY.ToString());
            Assert.AreEqual(77561752, a.GetTIC(1));
            Assert.AreEqual(144, a.GetSpectrumNumber(2));



            Assert.AreEqual(0.98, a.GetElapsedScanTime(100), 0.01);

            var cromatogram = a.GetTICChroma();

            Assert.AreEqual(360, cromatogram.Count);
            Assert.AreEqual(0.01, cromatogram.FirstTime, 0.002);
            Assert.AreEqual(2.788433333, cromatogram.PeakWithHighestY.Time, 0.0001);
            Assert.AreEqual(2.788433333, cromatogram.GetApex(0, 5).Time, 0.0001);

            var newSpectrum = new ThermoSpectrum(a.GetScan(51).MassSpectrum);
            Assert.AreEqual(22246 / 5574.8, newSpectrum.GetSignalToNoise(1), 0.01);


            Assert.AreEqual(1, newSpectrum.GetCharges()[1]);
            Assert.AreEqual(102604, newSpectrum.GetResolutions()[1]);

            Assert.AreEqual(181, newSpectrum.newSpectrumExtract(500, 1000).Count);

            Assert.AreEqual(0, newSpectrum.newSpectrumExtract(-3, -2).Count);

            var hm = newSpectrum.newSpectrumExtract(501, 502);

            Assert.AreEqual(0, hm.Count);
        }
    }
}