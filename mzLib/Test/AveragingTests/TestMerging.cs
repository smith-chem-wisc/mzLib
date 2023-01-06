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
        public static void TestSpectralBinning()
        {
            SpectralAveragingOptions options = new();
            MzSpectrum[] mzSpectras = new MzSpectrum[DummyMzSpectra.Count];
            DummyMzCopy.CopyTo(mzSpectras);
            double[][] compositeSpectraValues = SpectralMerging.CombineSpectra(mzSpectras.Select(p => p.XArray).ToArray(),
                mzSpectras.Select(p => p.YArray).ToArray(),
                mzSpectras.Select(p => p.SumOfAllY).ToArray(), mzSpectras.Count(), options);
            double[] expected = new double[] { 0, 3.2, 0, 0, 0, 0, 0, 6.4, 0 };
            Assert.That(compositeSpectraValues[0].Length == compositeSpectraValues[1].Length);
            Assert.That(expected.SequenceEqual(compositeSpectraValues[1]));

            options.PerformNormalization = false;
            DummyMzCopy.CopyTo(mzSpectras);
            compositeSpectraValues = SpectralMerging.CombineSpectra(mzSpectras.Select(p => p.XArray).ToArray(),
                mzSpectras.Select(p => p.YArray).ToArray(),
                mzSpectras.Select(p => p.SumOfAllY).ToArray(), mzSpectras.Count(), options);
            expected = new double[] { 0, 4, 0, 0, 0, 0, 0, 8, 0 };
            Assert.That(compositeSpectraValues[0].Length == compositeSpectraValues[1].Length);
            Assert.That(expected.SequenceEqual(compositeSpectraValues[1]));
        }

        [Test]
        public static void TestMostSimilarSpectrum()
        {
            SpectralAveragingOptions options = new() {SpectrumMergingType = SpectrumMergingType.MostSimilarSpectrum};
            MzSpectrum[] mzSpectras = new MzSpectrum[DummyMzSpectra.Count];
            DummyMzCopy.CopyTo(mzSpectras);
            try
            {
                double[][] compositeSpectraValues = SpectralMerging.CombineSpectra(
                    mzSpectras.Select(p => p.XArray).ToArray(),
                    mzSpectras.Select(p => p.YArray).ToArray(),
                    mzSpectras.Select(p => p.SumOfAllY).ToArray(), 
                    mzSpectras.Count(), options);
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

        [Test]
        public static void TestCombineSpectraError()
        {
            SpectralAveragingOptions options = new SpectralAveragingOptions();
            options.SpectrumMergingType = (SpectrumMergingType)(-1);
            MzSpectrum[] mzSpectras = new MzSpectrum[DummyMzSpectra.Count];
            DummyMzCopy.CopyTo(mzSpectras);
            try
            {
                double[][] compositeSpectraValues = SpectralMerging.CombineSpectra(
                    mzSpectras.Select(p => p.XArray).ToArray(),
                    mzSpectras.Select(p => p.YArray).ToArray(),
                    mzSpectras.Select(p => p.SumOfAllY).ToArray(),
                    mzSpectras.Count(), options);
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

        [Test]
        public static void TestProcessSingleMzArray()
        {
            double[] arr = new double[] { 1, 2, 3, 4, 5 };
            double[] arr1 = new double[] { 1, 0, 0, 0, 0 };
            double[] arr2 = new double[] { 0, 0, 0, 0, 0};
            double[] arr3 = new double[] { 100, 10, 1};
            SpectralAveragingOptions options = new SpectralAveragingOptions();
            options.RejectionType = RejectionType.NoRejection;
            options.PerformNormalization = false;
            options.RejectionType = RejectionType.NoRejection;

            var normalResult = SpectralMerging.ProcessSingleMzArray(arr, options);
            Assert.That(normalResult == 3);

            var allZerosResult = SpectralMerging.ProcessSingleMzArray(arr2, options);
            Assert.That(allZerosResult == 0);

            var allButOneZeroResult = SpectralMerging.ProcessSingleMzArray(arr1, options);
            Assert.That(allButOneZeroResult == 0);

            options.RejectionType = RejectionType.MinMaxClipping;
            var oneLeftResult = SpectralMerging.ProcessSingleMzArray(arr3, options);
            Assert.That(oneLeftResult == 0);
        }
    }
}
