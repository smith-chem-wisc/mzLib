using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using MassSpectrometry;
using MzLibSpectralAveraging;
using MzLibUtil;
using NUnit.Framework;

namespace Test.AveragingTests
{
    [ExcludeFromCodeCoverage]
    public static class GeneralMiscTests
    {
        [Test]
        public static void TestOptionToString()
        {
            SpectralAveragingOptions options1 = new SpectralAveragingOptions();
            SpectralAveragingOptions options2 = new SpectralAveragingOptions(){Percentile = 0.2, OutlierRejectionType = OutlierRejectionType.PercentileClipping};
            SpectralAveragingOptions options3 = new SpectralAveragingOptions()
                { OutlierRejectionType = OutlierRejectionType.SigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1 };
            SpectralAveragingOptions options4 = new SpectralAveragingOptions()
                { OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1 };
            SpectralAveragingOptions options5 = new SpectralAveragingOptions()
                { OutlierRejectionType = OutlierRejectionType.WinsorizedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1 };
            SpectralAveragingOptions options6 = new SpectralAveragingOptions()
            {
                OutlierRejectionType = OutlierRejectionType.WinsorizedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1,
                SpectralWeightingType = SpectraWeightingType.WeightEvenly
            };
            SpectralAveragingOptions options7 = new SpectralAveragingOptions()
            {
                OutlierRejectionType = OutlierRejectionType.WinsorizedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1,
                SpectralWeightingType = SpectraWeightingType.TicValue, PerformNormalization = false
            };

            string opt1String = options1.ToString();
            string opt2String = options2.ToString();
            string opt3String = options3.ToString();
            string opt4String = options4.ToString();
            string opt5String = options5.ToString();
            string opt6String = options6.ToString();
            string opt7String = options7.ToString();

            Assert.That(opt1String.Equals("NoRejection_None_Normalized_BinSize-0.01"));
            Assert.That(opt2String.Equals("PercentileClipping_None_Normalized_Percentile-0.2_BinSize-0.01"));
            Assert.That(opt3String.Equals("SigmaClipping_None_Normalized_MinSigma-1_MaxSigma-1_BinSize-0.01"));
            Assert.That(opt4String.Equals("AveragedSigmaClipping_None_Normalized_MinSigma-1_MaxSigma-1_BinSize-0.01"));
            Assert.That(opt5String.Equals("WinsorizedSigmaClipping_None_Normalized_MinSigma-1_MaxSigma-1_BinSize-0.01"));
            Assert.That(opt6String.Equals("WinsorizedSigmaClipping_None_Normalized_MinSigma-1_MaxSigma-1_BinSize-0.01"));
            Assert.That(opt7String.Equals("WinsorizedSigmaClipping_TicValue_MinSigma-1_MaxSigma-1_BinSize-0.01"));
        }

        [Test]
        public static void TestMzLibSpectralAveragingOptionsConstructors()
        {
            SpectralAveragingOptions defaultOptions = new();
            defaultOptions.SetDefaultValues();
            SpectralAveragingOptions options = new();
            // testing if values are equal to one another
            foreach (var property in defaultOptions.GetType().GetProperties())
            {
                var propertyType = property.PropertyType;
                var defaultProperty = Convert.ChangeType(defaultOptions.GetType().GetProperty(property.Name)?.GetValue(defaultOptions), propertyType);
                var testProperty = Convert.ChangeType(options.GetType().GetProperty(property.Name)?.GetValue(options),
                    propertyType);
                if (propertyType == typeof(SpectralAveragingOptions))
                {
                    foreach (var averagingProperty in ((SpectralAveragingOptions)defaultProperty)?.GetType().GetProperties())
                    {
                        // ensure wrapped property is equal to base property
                        defaultProperty =
                            Convert.ChangeType(
                                defaultOptions.GetType().GetProperty(averagingProperty.Name)?.GetValue(defaultOptions),
                                averagingProperty.PropertyType);
                        var defaultAvgProperty =
                            Convert.ChangeType(
                                defaultOptions.GetType().GetProperty(averagingProperty.Name)
                                    ?.GetValue(defaultOptions), averagingProperty.PropertyType);
                        Assert.That(defaultProperty?.ToString() == defaultAvgProperty?.ToString());

                        // ensure wrapped property is equal to base property
                        testProperty =
                            Convert.ChangeType(
                                options.GetType().GetProperty(averagingProperty.Name)
                                    ?.GetValue(options), averagingProperty.PropertyType);
                        var testAvgProperty =
                            Convert.ChangeType(
                                options.GetType().GetProperty(averagingProperty.Name)
                                    ?.GetValue(options),
                                averagingProperty.PropertyType);
                        Assert.That(testProperty?.ToString() == testAvgProperty?.ToString());
                        Assert.That(testAvgProperty?.ToString() == defaultAvgProperty?.ToString());
                    }
                }
                else
                {
                    Assert.That(defaultProperty?.ToString() == testProperty?.ToString());
                }
            }

            options.OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping;
            options.SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation;
            options.PerformNormalization = false;
            options.Percentile = 2;
            options.MinSigmaValue = 2;
            options.MaxSigmaValue = 2;
            options.BinSize = 2;

            Assert.That(options.OutlierRejectionType == OutlierRejectionType.AveragedSigmaClipping);
            Assert.That(options.SpectralWeightingType == SpectraWeightingType.MrsNoiseEstimation);
            Assert.That(options.PerformNormalization == false);
            Assert.That(Math.Abs(options.Percentile - 2) < 0.001);
            Assert.That(Math.Abs(options.MinSigmaValue - 2) < 0.001);
            Assert.That(Math.Abs(options.MaxSigmaValue - 2) < 0.001);
            Assert.That(Math.Abs(options.BinSize - 2) < 0.001);


        }

        [Test]
        [TestCase("DataFiles/small.RAW", 48, "Thermo nativeID format")]
        [TestCase("DataFiles/sliced_ethcd.raw", 6, "Thermo nativeID format")]
        [TestCase("DataFiles/SmallCalibratibleYeast.mzml",142, "Thermo nativeID format")]
        [TestCase("DataFiles/tester.mzML", 7, null)]
        public static void TestLoadingRawFilesAndSourceFiles(string filePath, int expectedScanCount, string sourceFormat)
        {
            string spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, filePath);
            List<MsDataScan> scans = SpectraFileHandler.LoadAllScansFromFile(spectraPath);
            Assert.That(scans.Count == expectedScanCount);

            SourceFile file = SpectraFileHandler.GetSourceFile(spectraPath);
            Assert.That(file.NativeIdFormat == sourceFormat);
        }

        [Test]
        public static void TestLoadingFileErrors()
        {
            string badPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles/small.toml");

            var exception = Assert.Throws<MzLibException>(() =>
            {
                SpectraFileHandler.LoadAllScansFromFile(badPath);
            });
            Assert.That(exception.Message == "Cannot load spectra");

            exception = Assert.Throws<MzLibException>(() =>
            {
                SpectraFileHandler.GetSourceFile(badPath);
            });
            Assert.That(exception.Message == "Cannot access SourceFile");
        }
    }
}
