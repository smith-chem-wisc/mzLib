using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FlyabilityModels;

namespace Test.KoinaTests.BenchmarkTests
{
    [TestFixture]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class DetectabilityBenchmarkTests
    {
        private static HashSet<string> GenerateUniquePeptides(int count, int length)
        {
            var aminoacids = "ACDEFGHIKLMNPQRSTVWY".ToArray();
            var peptides = new HashSet<string>();
            while (peptides.Count < count)
            {
                peptides.Add(new string(Random.Shared.GetItems(aminoacids, length)));
            }
            return peptides;
        }

        [Test]
        [Explicit("Massive test, takes a long time to run")]
        [Category("Performance Benchmark")]
        /// <summary>
        /// Performance benchmark for pfly_2024_fine_tuned.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 128 peptides, Predict() takes about 0.75 seconds to run.
        ///  - With 500,000 unique peptides of length 40, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 4 million peptides without hitting client issues.
        /// </summary>
        public static void BenchmarkPFly2024FineTuned()
        {
            var peptides = GenerateUniquePeptides(500000, 40);
            var modelInputs = peptides.Select(p => new DetectabilityPredictionInput(p)).ToList();
            var model = new PFly2024FineTuned();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(500000));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict 500,000 peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }
    }
}
