using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Nett;
using Test.FileReadingTests;
using UsefulProteomicsDatabases;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using Stopwatch = System.Diagnostics.Stopwatch;
using TopDownProteomics;

namespace Test.FlashLFQ
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestFlashLfqParameters
    {
        [Test]
        public static void TestTomlReadAndWrite()
        {
            FlashLfqParameters flashParams = new FlashLfqParameters() { 
                Normalize = true,
                MatchBetweenRuns = true,
                PpmTolerance = 7,
                IsotopePpmTolerance = 6,
                IsoTracker = true,
                DonorCriterion = DonorCriterion.Intensity,
                DonorQValueThreshold = 0.005,
                BayesianProteinQuant = true
            };

            // write the parameters to a toml file
            string tomlFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "FlashLfqParams.toml");
            Toml.WriteFile(flashParams, tomlFilePath);

            var readToml = Toml.ReadFile<FlashLfqParameters>(tomlFilePath);

            // Check non-default values specified in the original parameters object
            Assert.That(readToml.Normalize);
            Assert.That(readToml.MatchBetweenRuns);
            Assert.That(readToml.BayesianProteinQuant);
            Assert.That(readToml.IsoTracker);
            Assert.That(readToml.IsotopePpmTolerance, Is.EqualTo(6));
            Assert.That(readToml.PpmTolerance, Is.EqualTo(7));
            Assert.That(readToml.DonorCriterion, Is.EqualTo(DonorCriterion.Intensity));
            Assert.That(readToml.DonorQValueThreshold, Is.EqualTo(0.005));

            // Check default values
            Assert.That(readToml.NumIsotopesRequired, Is.EqualTo(2));
            Assert.That(readToml.IdSpecificChargeState, Is.False);
            Assert.That(readToml.UseSharedPeptidesForProteinQuant, Is.False);
            Assert.That(readToml.QuantifyAmbiguousPeptides, Is.False);
            Assert.That(readToml.Silent, Is.False);
            Assert.That(readToml.MaxThreads, Is.EqualTo(-1));
            Assert.That(readToml.Integrate, Is.False);
            Assert.That(readToml.RequireMsmsIdInCondition, Is.False);
            Assert.That(readToml.MbrPpmTolerance, Is.EqualTo(10));
            Assert.That(readToml.MaxMbrRtWindow, Is.EqualTo(1.0));
            Assert.That(readToml.MbrQValueThreshold, Is.EqualTo(0.05));
            Assert.That(readToml.ProteinQuantFoldChangeCutoff, Is.EqualTo(0.1));
            Assert.That(readToml.ProteinQuantBaseCondition, Is.Null);
            Assert.That(readToml.McmcSteps, Is.EqualTo(3000));
            Assert.That(readToml.McmcBurninSteps, Is.EqualTo(1000));
            Assert.That(readToml.PairedSamples, Is.False);
            Assert.That(readToml.RandomSeed, Is.Null);
        }   
    }
}
