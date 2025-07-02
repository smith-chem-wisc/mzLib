using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Algorithms;
using MassSpectrometry.Deconvolution.Parameters;
using MzLibUtil;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class TestFlashDeconv2Deconvolution
    {
        private FlashDeconv2 _flashDeconv2;

        [SetUp]
        public void Setup()
        {
            var deconParams = new FlashDeconvDeconvolutionParameters(1,60);
            _flashDeconv2 = new FlashDeconv2(deconParams);
        }

        [Test]
        // Tests that LogTransformSpectrum filters out low-intensity peaks and applies log transform to X values.
        public void LogTransformSpectrum_FiltersAndTransformsCorrectly()
        {
            var x = new[] { 100.0, 200.0, 300.0 };
            var y = new[] { 0.005, 0.02, 0.5 };
            var spectrum = new MzSpectrum(x, y, false);

            var result = _flashDeconv2.GetType()
                .GetMethod("LogTransformSpectrum", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(_flashDeconv2, new object[] { spectrum, 0.01 }) as MzSpectrum;

            Assert.That(result.XArray.Length, Is.EqualTo(2));
            Assert.That(result.YArray.Length, Is.EqualTo(2));
            Assert.That(result.XArray, Is.All.GreaterThan(0));
        }

        [Test]
        // Tests that AllAcceptibleLogMzDifferencesForAdjacentValues returns the correct number of differences and correct values.
        public void AllAcceptibleLogMzDifferencesForAdjacentValues_ReturnsExpectedCount()
        {
            var result = _flashDeconv2.GetType()
                .GetMethod("AllAcceptibleLogMzDifferencesForAdjacentValues", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(_flashDeconv2, new object[] { 1, 5 }) as List<double>;

            Assert.That(result.Count, Is.EqualTo(5));
            Assert.That(result[0], Is.EqualTo(Math.Log(2) - Math.Log(1)).Within(1e-10));
        }

        [Test]
        // Tests that FindMatchingGroups finds at least one group with more than one peak when differences match.
        public void FindMatchingGroups_FindsExpectedGroups()
        {
            double[] logX = { Math.Log(100), Math.Log(200), Math.Log(300) };
            double[] y = { 10, 20, 30 };
            var diffs = new List<double> { Math.Log(200) - Math.Log(100), Math.Log(300) - Math.Log(200) };

            var result = typeof(FlashDeconv2)
                .GetMethod("FindMatchingGroups", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static)
                .Invoke(null, new object[] { logX, y, diffs }) as List<(double[] X, double[] Y, int[] ChargeState)>;

            Assert.That(result, Is.Not.Null);
            Assert.That(result.Count, Is.GreaterThanOrEqualTo(1));
            Assert.That(result[0].X.Length, Is.GreaterThan(1));
        }

        [Test]
        // Tests that RemoveSubsetGroups removes groups that are subsets of larger groups.
        public void RemoveSubsetGroups_RemovesSubsets()
        {
            var groups = new List<(double[] X, double[] Y)>
            {
                (new double[] { 1, 2 }, new double[] { 10, 20 }),
                (new double[] { 1, 2, 3 }, new double[] { 10, 20, 30 })
            };

            var result = typeof(FlashDeconv2)
                .GetMethod("RemoveSubsetGroups", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static)
                .Invoke(null, new object[] { groups }) as List<(double[] X, double[] Y)>;

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].X.Length, Is.EqualTo(3));
        }

        [Test]
        // Tests that TransformGroupsToExpX correctly exponentiates the X values in each group.
        public void TransformGroupsToExpX_TransformsCorrectly()
        {
            var groups = new List<(double[] X, double[] Y)>
            {
                (new double[] { Math.Log(100), Math.Log(200) }, new double[] { 10, 20 })
            };

            var result = typeof(FlashDeconv2)
                .GetMethod("TransformGroupsToExpX", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static)
                .Invoke(null, new object[] { groups }) as List<(double[] X, double[] Y)>;

            Assert.That(result[0].X[0], Is.EqualTo(100).Within(1e-10));
            Assert.That(result[0].X[1], Is.EqualTo(200).Within(1e-10));
        }

        [Test]
        // Tests that FilterMassIntensityGroupsByPpmTolerance separates groups into likely correct and incorrect based on ppm tolerance.
        public void FilterMassIntensityGroupsByPpmTolerance_FiltersCorrectly()
        {
            var groups = new List<(double[] neutralMass, double[] intensity)>
            {
                (new double[] { 1000, 1000.01 }, new double[] { 10, 20 }),
                (new double[] { 1000, 2000 }, new double[] { 10, 20 })
            };

            FlashDeconv2.FilterMassIntensityGroupsByPpmTolerance(
                groups,
                out var likelyCorrect,
                out var likelyIncorrect,
                correctPpmTolerance: 20_000, // Large tolerance to force first group as correct
                incorrectPpmTolerance: 10);

            Assert.That(likelyCorrect.Count, Is.EqualTo(1));
            Assert.That(likelyIncorrect.Count, Is.EqualTo(1));
        }

        [Test]
        // Tests that GetMostCommonNeutralMassAndSummedIntensity finds the mode cluster and sums the correct intensities.
        public void GetMostCommonNeutralMassAndSummedIntensity_ReturnsExpected()
        {
            var groups = new List<(double[] neutralMass, double[] intensity)>
            {
                (new double[] { 1000, 1000.01, 2000 }, new double[] { 10, 20, 30 })
            };

            var result = FlashDeconv2.GetMostCommonNeutralMassAndSummedIntensity(groups, ppmTolerance: 20);

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].summedIntensity, Is.EqualTo(30)); // Only 2000 is outside 20 ppm, so 10+20=30
        }

        [Test]
        // Tests that CreateNeutralMassIntensityGroups creates a group with the expected number of neutral masses.
        public void CreateNeutralMassIntensityGroups_CreatesExpectedGroups()
        {
            var groups = new List<(double[] X, double[] Y, int[] ChargeState)>
            {
                (new double[] { Math.Log(100) }, new double[] { 10 }, new int[] { 1 })
            };

            var result = typeof(FlashDeconv2)
                .GetMethod("CreateNeutralMassIntensityGroups", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static)
                .Invoke(null, new object[] { groups }) as List<(double[] neutralMass, double[] intensity)>;

            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].neutralMass.Length, Is.EqualTo(1));
        }

        [Test]
        // Tests that LogMzDependentTolerance returns a positive tolerance value for a given log(m/z).
        public void LogMzDependentTolerance_ReturnsExpectedTolerance()
        {
            double logMz = Math.Log(1000);
            var result = typeof(FlashDeconv2)
                .GetMethod("LogMzDependentTolerance", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static)
                .Invoke(null, new object[] { logMz, 250.0 });

            Assert.That(result, Is.TypeOf<double>());
            Assert.That((double)result, Is.GreaterThan(0));
        }

        [Test]
        // Tests that NeutralMassFromLogMz returns a positive neutral mass for a given log(m/z) and charge.
        public void NeutralMassFromLogMz_ReturnsExpectedMass()
        {
            double logMz = Math.Log(1000);
            int charge = 2;
            var result = typeof(FlashDeconv2)
                .GetMethod("NeutralMassFromLogMz", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static)
                .Invoke(null, new object[] { logMz, charge });

            Assert.That(result, Is.TypeOf<double>());
            Assert.That((double)result, Is.GreaterThan(0));
        }

        [Test]
        // Tests that Deconvolute returns an empty result when given an empty spectrum.
        public void Deconvolute_ReturnsEmptyListForEmptySpectrum()
        {
            var spectrum = new MzSpectrum(new double[0], new double[0], false);
            var range = new MzRange(0, 1000);

            var result = _flashDeconv2.Deconvolute(spectrum, range);

            Assert.That(result, Is.Empty);
        }
    }
}
