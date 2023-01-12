using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IO.MzML;
using MassSpectrometry;
using MzLibSpectralAveraging;
using NUnit.Framework;
using MzLibUtil;

namespace Test.AveragingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class TestAveragingExtensions
    {
        public static List<MzSpectrum> DummyMzSpectra { get; set; }
        public static List<MsDataScan> ActualScans { get; set; }
        public static List<MzSpectrum> DummyMzCopy
        {
            get
            {
                List<MzSpectrum> newList = new();
                foreach (var spec in DummyMzSpectra)
                {
                    newList.Add(new MzSpectrum(spec.XArray, spec.YArray, true));
                }
                return newList;
            }
        }

        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            ActualScans = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"AveragingTestData\TDYeastFractionMS1.mzML")).GetAllScansList();
            double[] xArray = new double[] { 100.1453781, 200, 300, 400, 500, 600, 700, 800, 900.4123745 };
            double[] yArray1 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0, 0 };
            double[] yArray2 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0, 0 };
            double[] yArray3 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0, 0 };
            double[] yArray4 = new double[] { 0, 5, 0, 0, 0, 0, 0, 10, 0, 0 };
            double[] yArray5 = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            List<MzSpectrum> mzSpectra = new();
            mzSpectra.Add(new(xArray, yArray1, true));
            mzSpectra.Add(new(xArray, yArray2, true));
            mzSpectra.Add(new(xArray, yArray3, true));
            mzSpectra.Add(new(xArray, yArray4, true));
            mzSpectra.Add(new(xArray, yArray5, true));
            mzSpectra.Add(new(xArray, yArray1, true));
            mzSpectra.Add(new(xArray, yArray2, true));
            mzSpectra.Add(new(xArray, yArray3, true));
            mzSpectra.Add(new(xArray, yArray4, true));
            mzSpectra.Add(new(xArray, yArray5, true));

            DummyMzSpectra = mzSpectra;
        }

        [Test]
        [TestCase(5)]
        [TestCase(10)]
        public static void TestAverageExtensions(int scansToTake)
        {
            SpectralAveragingParameters parameters = new();
            List<MsDataScan> scans = ActualScans.Take(scansToTake).ToList();
            List<MzSpectrum> spectra = new();

            foreach (var scan in scans)
            {
                var spec = new MzSpectrum(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray, true);
                spectra.Add(spec);
            }

            var xArrays = spectra.Select(p => p.XArray).ToArray();
            var yArrays = spectra.Select(p => p.YArray.SubArray(0, p.YArray.Length)).ToArray();
            var xyJagged = SpectralAveraging.AverageSpectra(xArrays, yArrays, parameters);

            MzSpectrum compositeArrayMzSpectrum = new MzSpectrum(xyJagged[0], xyJagged[1], true);
            MzSpectrum compositeMsDataScanSpectrum = scans.AverageSpectra(parameters);
            MzSpectrum compositeMzSpectrumSpectrum = spectra.AverageSpectra(parameters);

            Assert.That(compositeMzSpectrumSpectrum.Equals(compositeMsDataScanSpectrum));
            Assert.That(compositeMzSpectrumSpectrum.Equals(compositeArrayMzSpectrum));
            Assert.That(compositeMsDataScanSpectrum.Equals(compositeArrayMzSpectrum));
        }

        [Test]
        [TestCase(5, NormalizationType.NoNormalization)]
        [TestCase(5, NormalizationType.AbsoluteToTic)]
        [TestCase(5, NormalizationType.RelativeToTics)]
        [TestCase(10, NormalizationType.NoNormalization)]
        [TestCase(10, NormalizationType.AbsoluteToTic)]
        [TestCase(10, NormalizationType.RelativeToTics)]
        public static void TestNormalizationExtensions(int scansToTake, NormalizationType type)
        {
            // setup values
            SpectralAveragingParameters parameters = new() {NormalizationType = type};
            List<MsDataScan> scans = ActualScans.Take(scansToTake).ToList();
            List<MzSpectrum> spectra = new();
            foreach (var scan in scans)
            {
                var spec = new MzSpectrum(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray, true);
                spectra.Add(spec);
            }
            var yArrays = spectra.Select(p => p.YArray.SubArray(0, p.YArray.Length)).ToArray();

            // normalize
            scans.NormalizeSpectra(type);
            spectra.NormalizeSpectra(type);
            SpectraNormalization.NormalizeSpectra(yArrays, type);

            var msDataScanyArrays = scans.Select(p => p.MassSpectrum.YArray).ToArray();
            var mzSpectrumyArrays = spectra.Select(p => p.YArray).ToArray();

            for (int i = 0; i < yArrays.Length; i++)
            {
                Assert.That(yArrays[i].SequenceEqual(msDataScanyArrays[i]));
                Assert.That(yArrays[i].SequenceEqual(mzSpectrumyArrays[i]));
            }
        }
    }
}
