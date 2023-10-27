using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using SpectralAveraging;


namespace Test.AveragingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class TestAveraging
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

        private static double[][] xArrays;
        private static double[][] yArrays;

        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            ActualScans = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory,
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
            xArrays = new[]
            {
                new double[] { 0, 1, 2, 3, 3.49, 4 },
                new double[] { 0, 1, 2, 3, 4 },
                new double[] { 0.1, 1.1, 2.1, 3.1, 4.1}
            };
            yArrays = new[]
            {
                new double[] { 10, 11, 12, 12, 13, 14 },
                new double[] { 11, 12, 13, 14, 15 },
                new double[] { 20, 25, 30, 35, 40 }
            };
        }

        [Test]
        public static void TestMzBinning()
        {
            SpectralAveragingParameters parameters = new();
            MzSpectrum[] mzSpectras = new MzSpectrum[DummyMzSpectra.Count];
            DummyMzCopy.CopyTo(mzSpectras);
            var compositeSpectra = mzSpectras.AverageSpectra(parameters);

            double[] expected = new[] { 3.2, 6.4};
            Assert.That(compositeSpectra.XArray.Length == compositeSpectra.YArray.Length);
            Assert.That(expected.SequenceEqual(compositeSpectra.YArray));

            parameters.NormalizationType = NormalizationType.NoNormalization;
            DummyMzCopy.CopyTo(mzSpectras);
            compositeSpectra = mzSpectras.AverageSpectra(parameters);
            expected = new[] { 4.0, 8.0};
            Assert.That(compositeSpectra.XArray.Length == compositeSpectra.YArray.Length);
            Assert.That(expected.SequenceEqual(compositeSpectra.YArray));
        }

        [Test]
        public static void TestWhenAllValuesGetRejected()
        {
            SpectralAveragingParameters parameters = new()
            {
                SpectraFileAveragingType = SpectraFileAveragingType.AverageAll,
                OutlierRejectionType = OutlierRejectionType.MinMaxClipping,
                SpectralWeightingType = SpectraWeightingType.WeightEvenly,
                NormalizationType = NormalizationType.NoNormalization,
                BinSize = 1,
            };

            List<MzSpectrum> spectra = new()
            {
                new MzSpectrum(new[] { 2.0, 3.0, 3.1, 4.0, 4.1 }, new[] { 1.0, 1.1, 1.2, 1.5, 1.6 }, false),
                new MzSpectrum(new[] { 2.0, 3.0, 3.1, 4.0, 4.1 }, new[] { 1.0, 1.3, 1.4, 1.7, 1.8 }, false),
            };

            var averagedSpectra = spectra.AverageSpectra(parameters);
            Assert.That(averagedSpectra.XArray.SequenceEqual(new [] {3.05, 4.05}));
            Assert.That(averagedSpectra.YArray.SequenceEqual(new[]{1.25, 1.65}));
        }

        [Test]
        public static void TestBinningMethod()
        {
            SpectralAveragingParameters parameters = new() { BinSize = 1 };

            var methodInfo = typeof(SpectraAveraging).GetMethod("GetBins", BindingFlags.NonPublic | BindingFlags.Static);
            var bins = (Dictionary<int, List<BinnedPeak>>)methodInfo.Invoke(null, new object?[] { xArrays, yArrays, parameters.BinSize });
            Assert.That(bins != null);
            Assert.That(bins.Count == 5);

            var bin = bins[0];
            Assert.That(bin.Select(p => p.Mz).SequenceEqual(new [] {0, 0, 0.1}));
            Assert.That(bin.Select(p => p.Intensity).SequenceEqual(new [] {10.0, 11, 20}));

            bin = bins[1];
            Assert.That(bin.Select(p => p.Mz).SequenceEqual(new[] { 1, 1, 1.1 }));
            Assert.That(bin.Select(p => p.Intensity).SequenceEqual(new[] { 11.0, 12, 25 }));

            bin = bins[2];
            Assert.That(bin.Select(p => p.Mz).SequenceEqual(new[] { 2, 2, 2.1 }));
            Assert.That(bin.Select(p => p.Intensity).SequenceEqual(new[] { 12, 13, 30.0 }));

            bin = bins[3];
            Assert.That(bin.Select(p => p.Mz).SequenceEqual(new[] { 3, 3.49, 3, 3.1 }));
            Assert.That(bin.Select(p => p.Intensity).SequenceEqual(new[] { 12.0, 13, 14, 35 }));

            bin = bins[4];
            Assert.That(bin.Select(p => p.Mz).SequenceEqual(new[] { 4, 4, 4.1 }));
            Assert.That(bin.Select(p => p.Intensity).SequenceEqual(new[] { 14.0, 15, 40 }));
        }

        [Test]
        public static void TestAverageSpectraError()
        {
            SpectralAveragingParameters parameters = new SpectralAveragingParameters();
            parameters.SpectralAveragingType = (SpectralAveragingType)(-1);

            var exception = Assert.Throws<MzLibException>(() =>
            {
                DummyMzSpectra.AverageSpectra(parameters);
            });
            Assert.That(exception.Message == "Spectrum Averaging Type Not Yet Implemented");
        }

        [Test]
        public static void TestBinnedPeakStruct()
        {
            BinnedPeak peak = new(1, 20.0, 25.5, 1);
            Assert.That(peak.Bin, Is.EqualTo(1));
            Assert.That(peak.Mz, Is.EqualTo(20.0));
            Assert.That(peak.Intensity, Is.EqualTo(25.5));
            Assert.That(peak.SpectraId, Is.EqualTo(1));
            Assert.That(peak.ToString(), Is.EqualTo("20 : 25.5 : 1"));
        }

    }
}
