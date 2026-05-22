using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;

namespace Test.KoinaTests.BenchmarkTests
{
    [TestFixture]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class FragmentIntensityBenchmarkTests
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

        // ═══════════════════════════════════════════════════════════════════════════
        // Prosit Family
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        [Explicit("Massive test, takes a long time to run")]
        [Category("Performance Benchmark")]
        /// <summary>
        /// Performance benchmark for Prosit_2019_intensity.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkProsit2019Intensity()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, null)).ToList();
            var model = new Prosit2019Intensity();
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
        /// Performance benchmark for Prosit_2020_intensity_HCD.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes under 10 seconds to run.
        ///  - With 1,000 unique peptides of length 30, the test takes approximately 36 minutes to run.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - For anything above 500k peptides, the MaxNumberOfBatchesPerRequest should be less than 500.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 4 million peptides without hitting client issues.
        /// </summary>
        public static void BenchmarkProsit2020IntensityHCD()
        {
            var peptides = GenerateUniquePeptides(1000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, null)).ToList();
            var model = new Prosit2020IntensityHCD();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(1000));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict 1,000 peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }

        [Test]
        [Explicit("Massive test, takes a long time to run")]
        [Category("Performance Benchmark")]
        /// <summary>
        /// Performance benchmark for Prosit_2020_intensity_CID.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkProsit2020IntensityCID()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, "CID")).ToList();
            var model = new Prosit2020IntensityCID();
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
        /// Performance benchmark for Prosit_2020_intensity_TMT.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - Requires N-terminal TMT/iTRAQ labeling on all peptides.
        /// </summary>
        public static void BenchmarkProsit2020IntensityTMT()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput("[Common Fixed:TMT6plex on N-terminus]" + p, 2, 35, null, "HCD")).ToList();
            var model = new Prosit2020IntensityTMT();
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
        /// Performance benchmark for Prosit_2023_intensity_timsTOF.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkProsit2023IntensityTimsTOF()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, null)).ToList();
            var model = new Prosit2023IntensityTimsTOF();
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
        /// Performance benchmark for Prosit_2024_intensity_cit.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - Requires fragmentation type input.
        /// </summary>
        public static void BenchmarkProsit2024IntensityCit()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, "HCD")).ToList();
            var model = new Prosit2024IntensityCit();
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
        /// Performance benchmark for Prosit_2024_intensity_PTMs_gl.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - Requires fragmentation type input.
        /// </summary>
        public static void BenchmarkProsit2024IntensityPTMsGl()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, "HCD")).ToList();
            var model = new Prosit2024IntensityPTMsGl();
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
        /// Performance benchmark for Prosit_2025_intensity_22PTM.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - Requires fragmentation type input.
        /// </summary>
        public static void BenchmarkProsit2025Intensity22PTM()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, "HCD")).ToList();
            var model = new Prosit2025Intensity22PTM();
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
        /// Performance benchmark for Prosit_2025_intensity_40PTM.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - Requires fragmentation type input.
        /// </summary>
        public static void BenchmarkProsit2025Intensity40PTM()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, "HCD")).ToList();
            var model = new Prosit2025Intensity40PTM();
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
        /// Performance benchmark for Prosit_2025_intensity_lac.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - Requires fragmentation type and instrument type inputs.
        /// </summary>
        public static void BenchmarkProsit2025IntensityLac()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, "QE", "HCD")).ToList();
            var model = new Prosit2025IntensityLac();
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
        /// Performance benchmark for Prosit_2025_intensity_ptm2.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - Requires fragmentation type input.
        /// </summary>
        public static void BenchmarkProsit2025IntensityPtm2()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, "HCD")).ToList();
            var model = new Prosit2025IntensityPtm2();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(500000));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict 500,000 peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // MS2PIP Family
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        [Explicit("Massive test, takes a long time to run")]
        [Category("Performance Benchmark")]
        /// <summary>
        /// Performance benchmark for ms2pip_HCD2021.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkMs2PipHCD2021()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, null)).ToList();
            var model = new Ms2PipHCD2021();
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
        /// Performance benchmark for ms2pip_CID_TMT.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - Requires N-terminal TMT/iTRAQ labeling on all peptides.
        /// </summary>
        public static void BenchmarkMs2PipCIDTMT()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput("[Common Fixed:TMT6plex on N-terminus]" + p, 2, 35, null, null)).ToList();
            var model = new Ms2PipCIDTMT();
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
        /// Performance benchmark for ms2pip_Immuno_HCD.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 12 (max 15 for this model), test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - MaxPeptideLength is 15 for this model (shorter than other MS2PIP models).
        /// </summary>
        public static void BenchmarkMs2PipImmunoHCD()
        {
            var peptides = GenerateUniquePeptides(500000, 12);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, null)).ToList();
            var model = new Ms2PipImmunoHCD();
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
        /// Performance benchmark for ms2pip_TTOF5600.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkMs2PipTTOF5600()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, null)).ToList();
            var model = new Ms2PipTTOF5600();
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
        /// Performance benchmark for ms2pip_iTRAQphospho.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - Requires iTRAQ/phospho modifications on peptides.
        /// </summary>
        public static void BenchmarkMs2PipITRAQPhospho()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, null)).ToList();
            var model = new Ms2PipITRAQPhospho();
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
        /// Performance benchmark for ms2pip_timsTOF2023.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkMs2PipTimsTOF2023()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, null)).ToList();
            var model = new Ms2PipTimsTOF2023();
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
        /// Performance benchmark for ms2pip_timsTOF2024.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkMs2PipTimsTOF2024()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, null)).ToList();
            var model = new Ms2PipTimsTOF2024();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(500000));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict 500,000 peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }

        // ═══════════════════════════════════════════════════════════════════════════
        // Other Fragment Intensity Models
        // ═══════════════════════════════════════════════════════════════════════════

        [Test]
        [Explicit("Massive test, takes a long time to run")]
        [Category("Performance Benchmark")]
        /// <summary>
        /// Performance benchmark for AlphaPeptDeep_ms2_generic.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - Requires instrument type input.
        /// </summary>
        public static void BenchmarkAlphaPeptDeepMs2Generic()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, "QE", null)).ToList();
            var model = new AlphaPeptDeepMs2Generic();
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
        /// Performance benchmark for Altimeter_2024_intensities.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        /// </summary>
        public static void BenchmarkAltimeter2024Intensities()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, null, null)).ToList();
            var model = new Altimeter2024Intensities();
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
        /// Performance benchmark for UniSpec.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500,000 unique peptides of length 20, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500,000 peptides.
        ///  - Requires instrument type input.
        /// </summary>
        public static void BenchmarkUniSpec()
        {
            var peptides = GenerateUniquePeptides(500000, 20);
            var modelInputs = peptides.Select(p => new FragmentIntensityPredictionInput(p, 2, 35, "QE", null)).ToList();
            var model = new UniSpec();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(500000));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict 500,000 peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }
    }
}
