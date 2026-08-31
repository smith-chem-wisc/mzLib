using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using SpectralAveraging;

namespace Test.AveragingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class TestParametersAndFileLoading
    {
        [Test]
        public static void TestOptionToString()
        {
            SpectralAveragingParameters options1 = new SpectralAveragingParameters();
            SpectralAveragingParameters options2 = new SpectralAveragingParameters(){Percentile = 0.2, OutlierRejectionType = OutlierRejectionType.PercentileClipping};
            SpectralAveragingParameters options3 = new SpectralAveragingParameters()
                { OutlierRejectionType = OutlierRejectionType.SigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1 };
            SpectralAveragingParameters options4 = new SpectralAveragingParameters()
                { OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1 };
            SpectralAveragingParameters options5 = new SpectralAveragingParameters()
                { OutlierRejectionType = OutlierRejectionType.WinsorizedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1 };
            SpectralAveragingParameters options6 = new SpectralAveragingParameters()
            {
                OutlierRejectionType = OutlierRejectionType.WinsorizedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1,
                SpectralWeightingType = SpectraWeightingType.WeightEvenly
            };
            SpectralAveragingParameters options7 = new SpectralAveragingParameters()
            {
                OutlierRejectionType = OutlierRejectionType.WinsorizedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1,
                SpectralWeightingType = SpectraWeightingType.TicValue, NormalizationType = NormalizationType.NoNormalization
            };

            string opt1String = options1.ToString();
            string opt2String = options2.ToString();
            string opt3String = options3.ToString();
            string opt4String = options4.ToString();
            string opt5String = options5.ToString();
            string opt6String = options6.ToString();
            string opt7String = options7.ToString();

            Assert.That(opt1String.Equals("NoRejection_WeightEvenly_RelativeToTics_BinSize-0.01"));
            Assert.That(opt2String.Equals("PercentileClipping_WeightEvenly_RelativeToTics_BinSize-0.01_Percentile-0.2"));
            Assert.That(opt3String.Equals("SigmaClipping_WeightEvenly_RelativeToTics_BinSize-0.01_MinSigma-1_MaxSigma-1"));
            Assert.That(opt4String.Equals("AveragedSigmaClipping_WeightEvenly_RelativeToTics_BinSize-0.01_MinSigma-1_MaxSigma-1"));
            Assert.That(opt5String.Equals("WinsorizedSigmaClipping_WeightEvenly_RelativeToTics_BinSize-0.01_MinSigma-1_MaxSigma-1"));
            Assert.That(opt6String.Equals("WinsorizedSigmaClipping_WeightEvenly_RelativeToTics_BinSize-0.01_MinSigma-1_MaxSigma-1"));
            Assert.That(opt7String.Equals("WinsorizedSigmaClipping_TicValue_NoNormalization_BinSize-0.01_MinSigma-1_MaxSigma-1"));
        }

        [Test]
        public static void SpectralAveragingParametersConstructor()
        {
            SpectralAveragingParameters defaultParameters = new();
            defaultParameters.SetDefaultValues();
            SpectralAveragingParameters parameters = new();
            // testing if values are equal to one another
            foreach (var property in defaultParameters.GetType().GetProperties())
            {
                var propertyType = property.PropertyType;
                var defaultProperty = Convert.ChangeType(defaultParameters.GetType().GetProperty(property.Name)?.GetValue(defaultParameters), propertyType);
                var testProperty = Convert.ChangeType(parameters.GetType().GetProperty(property.Name)?.GetValue(parameters),
                    propertyType);
                if (propertyType == typeof(SpectralAveragingParameters))
                {
                    foreach (var averagingProperty in ((SpectralAveragingParameters)defaultProperty)?.GetType().GetProperties())
                    {
                        // ensure wrapped property is equal to base property
                        defaultProperty =
                            Convert.ChangeType(
                                defaultParameters.GetType().GetProperty(averagingProperty.Name)?.GetValue(defaultParameters),
                                averagingProperty.PropertyType);
                        var defaultAvgProperty =
                            Convert.ChangeType(
                                defaultParameters.GetType().GetProperty(averagingProperty.Name)
                                    ?.GetValue(defaultParameters), averagingProperty.PropertyType);
                        Assert.That(defaultProperty?.ToString() == defaultAvgProperty?.ToString());

                        // ensure wrapped property is equal to base property
                        testProperty =
                            Convert.ChangeType(
                                parameters.GetType().GetProperty(averagingProperty.Name)
                                    ?.GetValue(parameters), averagingProperty.PropertyType);
                        var testAvgProperty =
                            Convert.ChangeType(
                                parameters.GetType().GetProperty(averagingProperty.Name)
                                    ?.GetValue(parameters),
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

            parameters.OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping;
            parameters.SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation;
            parameters.NormalizationType = NormalizationType.NoNormalization;
            parameters.Percentile = 2;
            parameters.MinSigmaValue = 2;
            parameters.MaxSigmaValue = 2;
            parameters.BinSize = 2;

            Assert.That(parameters.OutlierRejectionType == OutlierRejectionType.AveragedSigmaClipping);
            Assert.That(parameters.SpectralWeightingType == SpectraWeightingType.MrsNoiseEstimation);
            Assert.That(parameters.NormalizationType == NormalizationType.NoNormalization);
            Assert.That(Math.Abs(parameters.Percentile - 2) < 0.001);
            Assert.That(Math.Abs(parameters.MinSigmaValue - 2) < 0.001);
            Assert.That(Math.Abs(parameters.MaxSigmaValue - 2) < 0.001);
            Assert.That(Math.Abs(parameters.BinSize - 2) < 0.001);


        }

        [Test]
        [TestCase(new[] { 0.1 }, new[] { 5 }, new[] { 0 }, new[] { 1.3 }, new[] { 0.8 }, 
            new[] {SpectraWeightingType.MrsNoiseEstimation, SpectraWeightingType.WeightEvenly, SpectraWeightingType.TicValue}, 
            new[] {OutlierRejectionType.PercentileClipping, OutlierRejectionType.MinMaxClipping, OutlierRejectionType.SigmaClipping, 
                OutlierRejectionType.WinsorizedSigmaClipping, OutlierRejectionType.AveragedSigmaClipping, OutlierRejectionType.BelowThresholdRejection},
            new[] {NormalizationType.RelativeToTics, NormalizationType.NoNormalization})]
        [TestCase(new[] { 0.01, 0.1 }, new[] { 5, 10, 20 }, new[] { 2, 3, 4 }, new[] { 1.0, 2.0, 3.0 }, 
            new[] { 0.5, 0.7, 0.9 }, new[] { SpectraWeightingType.MrsNoiseEstimation, SpectraWeightingType.WeightEvenly, SpectraWeightingType.TicValue },
            new[] {OutlierRejectionType.PercentileClipping, OutlierRejectionType.MinMaxClipping, OutlierRejectionType.SigmaClipping,
                OutlierRejectionType.WinsorizedSigmaClipping, OutlierRejectionType.AveragedSigmaClipping, OutlierRejectionType.BelowThresholdRejection},
            new[] { NormalizationType.RelativeToTics, NormalizationType.NoNormalization })]
        [TestCase(new[] { 0.01, 0.1 }, new[] { 5, 10, 20 }, new[] { 2, 3, 4, 10, 15 }, new[] { 1.0, 2.0, 3.0, 4.0 },
            new[] { 0.5, 0.7, 0.9 }, new[] { SpectraWeightingType.MrsNoiseEstimation, SpectraWeightingType.WeightEvenly, SpectraWeightingType.TicValue },
            new[] {OutlierRejectionType.PercentileClipping, OutlierRejectionType.MinMaxClipping, OutlierRejectionType.SigmaClipping,
                OutlierRejectionType.WinsorizedSigmaClipping, OutlierRejectionType.AveragedSigmaClipping, OutlierRejectionType.BelowThresholdRejection},
            new[] { NormalizationType.RelativeToTics, NormalizationType.NoNormalization })]
        public static void TestGenerateSpectralAveragingParameters(double[] binSizes,
            int[] numberOfScansToAverage, int[] scanOverlap, double[] sigmas, double[] percentiles, 
            SpectraWeightingType[] weightingTypes, OutlierRejectionType[] outlierRejectionTypes, 
            NormalizationType[] normalizationTypes)
        {
            int rejectionTypes = outlierRejectionTypes
                .Count(p => !p.ToString().Contains("Sigma") && !p.ToString().Contains("Percent"));

            int sigmaTypes = (int)Math.Pow(sigmas.Length, 2) *
                             outlierRejectionTypes.Count(p => p.ToString().Contains("Sigma"));

            int scanToAverageCount = 0;
            foreach (var scanCount in numberOfScansToAverage)
            {
                foreach (var overlap in scanOverlap)
                {
                    if (overlap < scanCount)
                        scanToAverageCount++;
                }
            }

            var averagingParamCount =
                weightingTypes.Length * normalizationTypes.Length * (sigmaTypes + percentiles.Length + rejectionTypes) * binSizes.Length * scanToAverageCount;


            var result = SpectralAveragingParameters.GenerateSpectralAveragingParameters(binSizes,
                numberOfScansToAverage, scanOverlap, sigmas, percentiles, weightingTypes, outlierRejectionTypes,
                normalizationTypes);
            Assert.That(result.Count, Is.EqualTo(averagingParamCount));
        }
    }
}
