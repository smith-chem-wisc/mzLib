using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.CCSModels;

namespace Test.KoinaTests.BenchmarkTests
{
    [TestFixture]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class CCSBenchmarkTests
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
        /// Performance benchmark for IM2Deep.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkIM2Deep()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new CCSPredictionInput(p, 2)).ToList();
            var model = new IM2Deep();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(500000));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict 500,000 peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }

        [Test]
        [Explicit("Massive test, takes a long time to run")]
        [Category("Performance Benchmark")]
        /// <summary>
        /// Performance benchmark for AlphaPeptDeep_ccs_generic.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkAlphaPeptDeepCCSGeneric()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new CCSPredictionInput(p, 2)).ToList();
            var model = new AlphaPeptDeepCCSGeneric();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(500000));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict 500,000 peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }
    }
}
