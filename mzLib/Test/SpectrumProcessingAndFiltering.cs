using Chemistry;
using IO.MzML;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class SpectrumProcessingAndFiltering
    {
        [Test]
        [TestCase(10, null, 100, false)]
        [TestCase(10, 1900, 100, false)]
        [TestCase(10, 400, 100, false)]
        [TestCase(10, 400, 100, true)]
        public static void TestFilteringPeaksTopN_MultipleWindows(int peaksToKeep, int? nominalWindowWidthDaltons, int peakCount, bool normalize)
        {
            double[] mzArray = new double[100];
            double[] intArray = new double[100];

            for (int i = 0; i < 100; i++)
            {
                mzArray[i] = 100d + ((double)i / 100d) * 1900;
                intArray[i] = i;
            }

            FilteringParams f = new FilteringParams(peaksToKeep, null, nominalWindowWidthDaltons, null, normalize, false, false);

            MsDataFile.WindowModeHelper(ref intArray, ref mzArray, f, 100, 2000, false);

            if (nominalWindowWidthDaltons.HasValue)
            {
                Assert.LessOrEqual((decimal)mzArray.Count(), (decimal)peaksToKeep * (decimal)nominalWindowWidthDaltons);
            }
            else
            {
                Assert.LessOrEqual((decimal)mzArray.Count(), (decimal)peaksToKeep * (decimal)1.0);
            }
            if (normalize)
            {
                Assert.That(50, Is.EqualTo(intArray.Max()).Within(0.1));
            }
        }

        [Test]
        public static void TestFilterLowIntensity()
        {
            double[] mzArray = new double[100];
            double[] intArray = new double[100];

            for (int i = 0; i < 100; i++)
            {
                mzArray[i] = (double)i;
                intArray[i] = (double)i;
            }

            FilteringParams f = new FilteringParams(null, 0.05, null, null, false, false, false);

            MsDataFile.WindowModeHelper(ref intArray, ref mzArray, f, mzArray.Min(), mzArray.Max());

            //The first five intensities are below 5% and therefore removed.
            Assert.AreEqual(95, intArray.Count());
            Assert.AreEqual(95, mzArray.Count());
        }

        [Test]
        public static void TestFilterNumberPeaksPerWindowOneWindow()
        {
            double[] mzArray = new double[200];
            double[] intArray = new double[200];

            for (int i = 0; i < 200; i++)
            {
                mzArray[i] = (double)(i + 1);
                intArray[i] = (double)(i + 1);
            }

            FilteringParams f = new FilteringParams(100, null, null, null, false, false, false);

            MsDataFile.WindowModeHelper(ref intArray, ref mzArray, f, mzArray.Min(), mzArray.Max());

            Assert.AreEqual(100, intArray.Count());
            Assert.AreEqual(100, mzArray.Count());
            Assert.AreEqual(101, intArray.Min());
            Assert.AreEqual(200, intArray.Max());
            Assert.AreEqual(101, mzArray.Min());
            Assert.AreEqual(200, mzArray.Max());
        }

        [Test]
        public static void TestFilterNumberPeaksPerWindowTenWindows()
        {
            double[] mzArray = new double[200];
            double[] intArray = new double[200];

            for (int i = 0; i < 200; i++)
            {
                mzArray[i] = (double)(i + 1);
                intArray[i] = (double)(i + 1);
            }

            FilteringParams f = new FilteringParams(10, null, 20, 10, false, false, false);

            MsDataFile.WindowModeHelper(ref intArray, ref mzArray, f, mzArray.Min(), mzArray.Max());

            Assert.AreEqual(100, intArray.Count());
            Assert.AreEqual(100, mzArray.Count());
            Assert.AreEqual(11, intArray.Min());
            Assert.AreEqual(200, intArray.Max());
            Assert.AreEqual(11, mzArray.Min());
            Assert.AreEqual(200, mzArray.Max());
        }

        [Test]
        public static void TestFilterKeepsPeaksWithHighestIntensity()
        {
            double[] mzArray = new double[200];
            double[] intArray = new double[200];

            Random r = new Random();

            for (int i = 0; i < 200; i++)
            {
                mzArray[i] = (double)(i + 1);
                intArray[i] = (i * Math.Abs(r.Next(1, 100)) + 1d);
            }

            List<double> l = intArray.ToList();

            l.Sort((x, y) => y.CompareTo(x));
            l = l.Take(100).ToList();

            FilteringParams f = new FilteringParams(100, null, null, null, false, false, false);

            MsDataFile.WindowModeHelper(ref intArray, ref mzArray, f, mzArray.Min(), mzArray.Max());

            List<double> myOut = intArray.ToList();
            myOut.Sort((x, y) => y.CompareTo(x));
            Assert.IsTrue(l.SequenceEqual(myOut));
        }

        [Test]
        public static void TestXcorrFilteringPeaksTopN_MultipleWindows()
        {
            List<double> masses = new List<double>();
            List<double> intensities = new List<double>();

            int startMass = 1;
            int incrementMass = 5;

            double startIntensity = 0.5;
            double incrementIntensity = 5;

            while (startMass < 1969)
            {
                masses.Add(startMass * 1.0005079);
                intensities.Add(startIntensity);
                startMass = startMass + incrementMass;
                startIntensity = startIntensity + incrementIntensity * startMass;
            }

            double[] mzArray = masses.ToArray();
            double[] intArray = intensities.ToArray();
            Array.Sort(mzArray, intArray);

            var spectrum = new MzSpectrum(mzArray, intArray, false);
            spectrum.XCorrPrePreprocessing(mzArray.Min(), mzArray.Max(), 241.122);

            //first mz rounded to nearest discrete mass bin 1.0005079
            Assert.AreEqual(Math.Round(96.0487584, 5), Math.Round(spectrum.XArray.Min(), 5));

            //last mz rounded to nearest discrete mass bin 1.0005079
            Assert.AreEqual(Math.Round(1966.998531, 5), Math.Round(spectrum.XArray.Max(), 5));

            //peaks within 1.5 thomson of precursor 241.122 are absent
            double precursorIntensity = 0;
            for (int i = 0; i < spectrum.XArray.Length; i++)
            {
                if (spectrum.XArray[i] > (241.122 - 1.5) && spectrum.XArray[i] < (241.122 + 1.5))
                {
                    precursorIntensity += spectrum.YArray[i];
                }
            }
            Assert.AreEqual(0, precursorIntensity);

            //not zero intensities
            Assert.AreEqual(374, spectrum.YArray.Where(i => i > 0).ToList().Count);

            //zero intensities. there shouldn't be any.
            Assert.AreEqual(0, spectrum.YArray.Where(i => i == 0).ToList().Count);

            //first peak with intensity
            Assert.AreEqual(Math.Round(21.170981474538145, 5), Math.Round(spectrum.YArray[0], 5));

            //last peak with intensity
            Assert.AreEqual(Math.Round(39.674211517419245, 5), Math.Round(spectrum.YArray[373], 5));
        }

        [Test]
        public static void ProcessXcorrInMzSpectrum()
        {
            Dictionary<string, MsDataFile> MyMsDataFiles = new Dictionary<string, MsDataFile>();
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "BinGenerationTest.mzML");
            FilteringParams filter = new FilteringParams(200, 0.01, null, 1, false, false, true);

            MyMsDataFiles[origDataFile] = Mzml.LoadAllStaticData(origDataFile, filter, 1);

            var scans = MyMsDataFiles[origDataFile].GetAllScansList();

            foreach (MsDataScan scan in scans.Where(s => s.MsnOrder > 1))
            {
                scan.MassSpectrum.XCorrPrePreprocessing(0, 1968 * 1.0005079, scan.IsolationMz.Value);
            }

            Assert.AreEqual(6, scans[0].MassSpectrum.XArray.Count());
            Assert.AreEqual(20, scans[1].MassSpectrum.XArray.Count());
        }

        [Test]
        public static void ProcessXcorrInB6MzSpectrum()
        {
            Dictionary<string, MsDataFile> MyMsDataFiles = new Dictionary<string, MsDataFile>();
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\sliced_b6.mzML");
            FilteringParams filter = new FilteringParams(200, 0.01, null, 1, false, false, false);

            string expectedResultFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\Working_86.tsv");

            List<string> expectedResults = File.ReadAllLines(expectedResultFile, Encoding.UTF8).ToList();

            MyMsDataFiles[origDataFile] = Mzml.LoadAllStaticData(origDataFile, filter, 1);

            var scans = MyMsDataFiles[origDataFile].GetAllScansList();

            List<double> xArrayProcessed = new List<double>();
            foreach (MsDataScan scan in scans.Where(s => s.MsnOrder > 1))
            {
                if (scan.OneBasedScanNumber == 86)
                {
                    scan.MassSpectrum.XCorrPrePreprocessing(0, 1969, scan.IsolationMz.Value);
                    xArrayProcessed = scan.MassSpectrum.XArray.ToList();
                }
            }

            for (int i = 0; i < expectedResults.Count; i++)
            {
                Assert.That(double.Parse(expectedResults[i]), Is.EqualTo(xArrayProcessed[i]).Within(0.001));
            }
        }

        [Test]
        public static void XcorrTestBorrowedFromMM()
        {
            double[] mzs = new double[] { 130.0499, 148.0604, 199.1077, 209.0921, 227.1026, 245.0768, 263.0874, 296.1605, 306.1448, 324.1554, 358.1609, 376.1714, 397.2082, 407.1925, 425.2031, 459.2086, 477.2191, 510.2922, 520.2766, 538.2871, 556.2613, 574.2719, 625.3192, 635.3035, 653.3141, 685.3039, 703.3145, 782.3567, 800.3672 };
            double[] intensities = new double[] { 20, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 };

            double precursorMass = 799.359968;

            var ms2 = new MzSpectrum(mzs, intensities, false);

            double rt = 1;
            int precursorZ = 1;

            var myScan = new MsDataScan(ms2, 1, 2, true, Polarity.Positive, rt + 0.01, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=2", precursorMass.ToMz(precursorZ), precursorZ, 1, precursorMass.ToMz(precursorZ), 1.0, DissociationType.HCD, 1, precursorMass.ToMz(precursorZ));

            ms2.XCorrPrePreprocessing(0, 1969, precursorMass.ToMz(1));

            List<double> X = new List<double>() { 130.07, 148.08, 199.1, 209.11, 227.12, 245.12, 263.13, 296.15, 306.16, 324.16, 358.18, 376.19, 397.2, 407.21, 425.22, 459.23, 477.24, 510.26, 520.26, 538.27, 556.28, 574.29, 625.32, 635.32, 653.33, 685.35, 703.36, 782.4 };
            List<double> Y = new List<double>() { 49.10, 34.13, 47.78, 48.11, 48.01, 47.68, 47.35, 47.68, 47.68, 47.68, 47.35, 47.68, 47.68, 47.68, 47.68, 47.68, 47.68, 47.68, 47.68, 48.01, 48.01, 47.68, 48.01, 48.01, 48.34, 48.34, 48.68, 49.67, 49.67 };

            for (int i = 0; i < X.Count; i++)
            {
                Assert.AreEqual(X[i], Math.Round(ms2.XArray[i], 2));
                Assert.AreEqual(Y[i], Math.Round(ms2.YArray[i], 2));
            }
        }

        [Test]
        public static void TestMSDataScan()
        {
            Dictionary<string, MsDataFile> MyMsDataFiles = new Dictionary<string, MsDataFile>();
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\sliced_b6.mzML");
            FilteringParams filter = new FilteringParams(200, 0.01, null, 1, false, false, false);

            MyMsDataFiles[origDataFile] = Mzml.LoadAllStaticData(origDataFile, filter, 1);

            var scans = MyMsDataFiles[origDataFile].GetAllScansList();

            MsDataScan[] ms1scans = scans.Where(s => s.MsnOrder == 1).ToArray();
            MsDataScan[] ms2scans = scans.Where(s => s.MsnOrder == 2).ToArray();

            MsDataScan ms1scanExample = ms1scans.Where(s => s.OneBasedScanNumber == 148).First();
            MsDataScan ms2scanExample = ms2scans.Where(s => s.OneBasedPrecursorScanNumber == 148).First();

            Assert.IsNull(ms1scanExample.DissociationType);
            Assert.IsNull(ms1scanExample.HcdEnergy);
            Assert.That(2.35, Is.EqualTo(ms1scanExample.InjectionTime).Within(0.01));
            Assert.IsTrue(ms1scanExample.IsCentroid);
            Assert.IsNull(ms1scanExample.IsolationMz);
            Assert.IsNull(ms1scanExample.IsolationRange);
            Assert.IsNull(ms1scanExample.IsolationWidth);
            Assert.AreEqual(2951, ms1scanExample.MassSpectrum.XArray.Length);
            Assert.AreEqual(1, ms1scanExample.MsnOrder);
            Assert.AreEqual(MZAnalyzerType.Quadrupole, ms1scanExample.MzAnalyzer);
            Assert.AreEqual("controllerType=0 controllerNumber=1 scan=148", ms1scanExample.NativeId);
            Assert.IsNull(ms1scanExample.NoiseData);
            Assert.IsNull(ms1scanExample.OneBasedPrecursorScanNumber);
            Assert.AreEqual(148, ms1scanExample.OneBasedScanNumber);
            Assert.AreEqual(Polarity.Positive, ms1scanExample.Polarity);
            Assert.That(77.08, Is.EqualTo(ms1scanExample.RetentionTime).Within(0.01));
            Assert.AreEqual("FTMS + p NSI Full ms [350.0000-1350.0000]", ms1scanExample.ScanFilter);
            Assert.AreEqual(350, ms1scanExample.ScanWindowRange.Minimum);
            Assert.AreEqual(1350, ms1scanExample.ScanWindowRange.Maximum);
            Assert.IsNull(ms1scanExample.SelectedIonChargeStateGuess);
            Assert.IsNull(ms1scanExample.SelectedIonIntensity);
            Assert.IsNull(ms1scanExample.SelectedIonMonoisotopicGuessIntensity);
            Assert.IsNull(ms1scanExample.SelectedIonMonoisotopicGuessMz);
            Assert.IsNull(ms1scanExample.SelectedIonMZ);
            Assert.AreEqual(1222352000, ms1scanExample.TotalIonCurrent);

            Assert.AreEqual(DissociationType.HCD, ms2scanExample.DissociationType);
            Assert.IsNull(ms2scanExample.HcdEnergy);
            Assert.That(15, Is.EqualTo(ms2scanExample.InjectionTime).Within(0.01));
            Assert.IsTrue(ms2scanExample.IsCentroid);
            Assert.That(518.23, Is.EqualTo(ms2scanExample.IsolationMz).Within(0.01));
            Assert.That(517.88, Is.EqualTo(ms2scanExample.IsolationRange.Minimum).Within(0.01));
            Assert.That(518.58, Is.EqualTo(ms2scanExample.IsolationRange.Maximum).Within(0.01));
            Assert.That(0.70, Is.EqualTo(ms2scanExample.IsolationWidth).Within(0.01));
            Assert.AreEqual(871, ms2scanExample.MassSpectrum.XArray.Length);
            Assert.AreEqual(2, ms2scanExample.MsnOrder);
            Assert.AreEqual(MZAnalyzerType.Quadrupole, ms2scanExample.MzAnalyzer);
            Assert.AreEqual("controllerType=0 controllerNumber=1 scan=149", ms2scanExample.NativeId);
            Assert.IsNull(ms2scanExample.NoiseData);
            Assert.AreEqual(148, ms2scanExample.OneBasedPrecursorScanNumber);
            Assert.AreEqual(149, ms2scanExample.OneBasedScanNumber);
            Assert.AreEqual(Polarity.Positive, ms2scanExample.Polarity);
            Assert.That(77.08, Is.EqualTo(ms2scanExample.RetentionTime).Within(0.01));
            Assert.AreEqual("ITMS + c NSI t d Full ms2 518.2317@hcd25.00 [200.0000-1200.0000]", ms2scanExample.ScanFilter);
            Assert.AreEqual(200, ms2scanExample.ScanWindowRange.Minimum);
            Assert.AreEqual(1200, ms2scanExample.ScanWindowRange.Maximum);
            Assert.AreEqual(3, ms2scanExample.SelectedIonChargeStateGuess);
            Assert.AreEqual(8984495, ms2scanExample.SelectedIonIntensity);
            Assert.IsNull(ms2scanExample.SelectedIonMonoisotopicGuessIntensity);
            Assert.That(518.23, Is.EqualTo(ms2scanExample.SelectedIonMonoisotopicGuessMz).Within(0.01));
            Assert.That(518.23, Is.EqualTo(ms2scanExample.SelectedIonMZ).Within(0.01));
            Assert.That(3748516.16, Is.EqualTo(ms2scanExample.TotalIonCurrent).Within(0.01));
        }
    }
}