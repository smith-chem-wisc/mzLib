using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MzLibSpectralAveraging;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSpectralAveragingExtensions
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
            List<MsDataScan> scans = SpectraFileHandler.LoadAllScansFromFile(SpectraPath).Take(5).ToList();
            List<MzSpectrum> spectra = new();
            foreach (var scan in scans)
            {
                spectra.Add(new(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray, true));
            }

            MzSpectrum compositeMsDataScanSpectrum = scans.AverageSpectra(Options);
            MzSpectrum compositeMzSpectrumSpectrum = spectra.AverageSpectra(Options);
            Assert.That(compositeMzSpectrumSpectrum.Equals(compositeMsDataScanSpectrum));
        }
    }
}
