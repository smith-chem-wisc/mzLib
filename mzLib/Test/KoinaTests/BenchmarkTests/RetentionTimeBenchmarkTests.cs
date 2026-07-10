using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;

namespace Test.KoinaTests.BenchmarkTests
{
    [TestFixture]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class RetentionTimeBenchmarkTests
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
        /// Performance benchmark for Prosit_2019_irt.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes under 1.5 seconds to run.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkProsit2019iRT()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput(p)).ToList();
            var model = new Prosit2019iRT();
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
        /// Performance benchmark for Prosit_2020_irt_TMT.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes under 1.5 seconds to run.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        ///  - All peptides must include a required N-terminal modification (TMT6plex used here).
        /// </summary>
        public static void BenchmarkProsit2020iRTTMT()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput("[Common Fixed:TMT6plex on N-terminus]" + p)).ToList();
            var model = new Prosit2020iRTTMT();
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
        /// Performance benchmark for Deeplc_hela_hf.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkDeeplcHelaHf()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput(p)).ToList();
            var model = new DeeplcHelaHf();
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
        /// Performance benchmark for AlphaPeptDeep_rt_generic.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkAlphaPeptDeepRTGeneric()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput(p)).ToList();
            var model = new AlphaPeptDeepRTGeneric();
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
        /// Performance benchmark for Chronologer_RT.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkChronologerRT()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput(p)).ToList();
            var model = new ChronologerRT();
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
        /// Performance benchmark for Prosit_2024_irt_cit.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkProsit2024iRTCit()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput(p)).ToList();
            var model = new Prosit2024iRTCit();
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
        /// Performance benchmark for Prosit_2024_irt_PTMs_gl.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkProsit2024iRTPTMsGl()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput(p)).ToList();
            var model = new Prosit2024iRTPTMsGl();
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
        /// Performance benchmark for Prosit_2025_irt_22PTM.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkProsit2025iRT22PTM()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput(p)).ToList();
            var model = new Prosit2025iRT22PTM();
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
        /// Performance benchmark for Prosit_2025_irt_40PTM.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkProsit2025iRT40PTM()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput(p)).ToList();
            var model = new Prosit2025iRT40PTM();
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
        /// Performance benchmark for Prosit_2025_irt_lac.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 30, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 500, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkProsit2025iRTLac()
        {
            var peptides = GenerateUniquePeptides(500000, 30);
            var modelInputs = peptides.Select(p => new RetentionTimePredictionInput(p)).ToList();
            var model = new Prosit2025iRTLac();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(500000));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict 500,000 peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }
    }
}
