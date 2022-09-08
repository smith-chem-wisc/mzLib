using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using NUnit.Framework;
using SpectralAveraging;
using SpectralAveragingExtensions;
using SpectralAveragingExtensions.Util;

namespace Test
{
    public class SpectralAveragingExtensionsTest
    {
        public static string OutputDirectory;
        public static string SpectraPath;
        public static SpectralAveragingOptions Options;
        public static List<MsDataScan> Scans;

        [OneTimeSetUp]
        public static void OneTimeSetUp()
        {
            Options = new SpectralAveragingOptions();
            OutputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"AveragingTestData");
            SpectraPath = Path.Combine(OutputDirectory, "TDYeastFractionMS1.mzML");
            Scans = SpectraFileHandler.LoadAllScansFromFile(SpectraPath);
        }

        [Test]
        public static void TestCombineSpectra()
        {
            MzSpectrum compositeMsDataScanSpectrum = Scans.Take(5).ToList().CombineSpectra(Options);
            MzSpectrum compositeMzSpectrumSpectrum = Scans.Select(p => p.MassSpectrum).Take(5).ToList()
                .CombineSpectra(Options);
            Assert.That(compositeMzSpectrumSpectrum.Equals(compositeMsDataScanSpectrum));
        }

        [Test]
        public static void TestSingleObjectNormalization()
        {
            MsDataScan scan = Scans.First();
            MzSpectrum spectrum = new(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray, true);
            double[] expectedY = scan.MassSpectrum.YArray.Select(p => p / scan.MassSpectrum.SumOfAllY).ToArray();

            scan.NormalizeSpectrumToTic();
            spectrum.NormalizeSpectrumToTic();
            Assert.That(spectrum.Equals(scan.MassSpectrum));

            Assert.That(scan.MassSpectrum.YArray.SequenceEqual(expectedY));
            Assert.That(spectrum.YArray.SequenceEqual(expectedY));
        }

        [Test]
        public static void TestMultipleObjectAbsoluteNormalization()
        {
            List<MsDataScan> scans = Scans.Take(5).ToList();
            List<MzSpectrum> spectra = new();
            foreach (var scan in scans)
            {
                spectra.Add(new(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray, true));
            }
            List<double[]> expectedYs = scans
                .Select(p => p.MassSpectrum.YArray.Select(m => m / p.MassSpectrum.SumOfAllY).ToArray()).ToList();

            scans.NormalizeSpectrumToTic(false);
            spectra.NormalizeSpectrumToTic(false);

            for (var i = 0; i < spectra.Count; i++)
            {
                Assert.That(scans[i].MassSpectrum.Equals(spectra[i]));
                Assert.That(scans[i].MassSpectrum.YArray.SequenceEqual(expectedYs[i]));
                Assert.That(spectra[i].YArray.SequenceEqual(expectedYs[i]));
            }
        }

        [Test]
        public static void TestMultipleObjectRelativeNormalization()
        {
            List<MsDataScan> scans = Scans.Take(5).ToList();
            List<MzSpectrum> spectra = new();
            foreach (var scan in scans)
            {
                spectra.Add(new(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray, true));
            }
            List<double[]> expectedYs = scans
                .Select(p => p.MassSpectrum.YArray.Select(m => m / p.MassSpectrum.SumOfAllY * spectra.Average(p => p.SumOfAllY)).ToArray()).ToList();

            scans.NormalizeSpectrumToTic(true);
            spectra.NormalizeSpectrumToTic(true);

            for (var i = 0; i < spectra.Count; i++)
            {
                Assert.That(scans[i].MassSpectrum.Equals(spectra[i]));
                Assert.That(scans[i].MassSpectrum.YArray.SequenceEqual(expectedYs[i]));
                Assert.That(spectra[i].YArray.SequenceEqual(expectedYs[i]));
            }
        }
    }
}
