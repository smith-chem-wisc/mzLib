using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using IO.MzML;
using MassSpectrometry;
using MzLibSpectralAveraging;
using NUnit.Framework;

namespace Test.AveragingTests
{
    [ExcludeFromCodeCoverage]
    public static class TestMerging
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
        public static void TestMzBinning()
        {
            SpectralAveragingParameters parameters = new();
            MzSpectrum[] mzSpectras = new MzSpectrum[DummyMzSpectra.Count];
            DummyMzCopy.CopyTo(mzSpectras);
            var compositeSpectra = mzSpectras.AverageSpectra(parameters);
                

            double[] expected = new double[] { 0, 3.2, 0, 0, 0, 0, 0, 6.4, 0 };
            Assert.That(compositeSpectra.XArray.Length == compositeSpectra.YArray.Length);
            Assert.That(expected.SequenceEqual(compositeSpectra.YArray));

            parameters.PerformNormalization = false;
            DummyMzCopy.CopyTo(mzSpectras);
            compositeSpectra = mzSpectras.AverageSpectra(parameters);
            expected = new double[] { 0, 4, 0, 0, 0, 0, 0, 8, 0 };
            Assert.That(compositeSpectra.XArray.Length == compositeSpectra.YArray.Length);
            Assert.That(expected.SequenceEqual(compositeSpectra.YArray));
        }

        [Test]
        public static void TestAverageSpectraExtensions()
        {
            SpectralAveragingParameters parameters = new();
            List<MsDataScan> scans = ActualScans.Take(5).ToList();
            List<MzSpectrum> spectra = new();
            foreach (var scan in scans)
            {
                spectra.Add(new(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray, true));
            }

            MzSpectrum compositeMsDataScanSpectrum = scans.AverageSpectra(parameters);
            MzSpectrum compositeMzSpectrumSpectrum = spectra.AverageSpectra(parameters);
            Assert.That(compositeMzSpectrumSpectrum.Equals(compositeMsDataScanSpectrum));
        }


        [Test]
        public static void TestAverageSpectraError()
        {
            SpectralAveragingParameters parameters = new SpectralAveragingParameters();
            parameters.SpectraMergingType = (SpectraMergingType)(-1);
            MzSpectrum[] mzSpectras = new MzSpectrum[DummyMzSpectra.Count];
            DummyMzCopy.CopyTo(mzSpectras);
            try
            {
                var compositeSpectraValues = mzSpectras.AverageSpectra(parameters);
                Assert.That(false);
            }
            catch (NotImplementedException)
            {

            }
            catch (Exception)
            {
                Assert.That(false);
            }
        }

    }
}
