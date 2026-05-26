using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.CrosslinkIntensityModels;

namespace Test.KoinaTests.BenchmarkTests
{
    [TestFixture]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class CrosslinkBenchmarkTests
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

        /// <summary>
        /// Converts a bare peptide into one with a crosslink annotation at a fixed position (index 4).
        /// The Koina helpers for CMS2/NMS2 require K[UNIMOD:NNNN] in the alpha sequence to determine
        /// the crosslink position for fragment annotation generation.
        /// </summary>
        private static string AnnotateCrosslink(string peptide, string unimodId)
        {
            if (peptide.Length < 6)
                return peptide[..1] + "K[UNIMOD:" + unimodId + "]" + peptide[2..];
            return peptide[..4] + "K[UNIMOD:" + unimodId + "]" + peptide[5..];
        }

        [Test]
        [Explicit("Massive test, takes a long time to run")]
        [Category("Performance Benchmark")]
        /// <summary>
        /// Performance benchmark for Prosit_2023_intensity_XL_CMS2.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500 unique crosslinked peptide pairs of length 15, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500 peptides.
        ///  - Beta peptides are reversed alpha peptides.
        ///  - Alpha peptides include K[UNIMOD:1896] for the CMS2 helper to locate crosslink.
        /// </summary>
        public static void BenchmarkProsit2023IntensityXLCMS2()
        {
            var alphaPeptides = GenerateUniquePeptides(500, 15).Select(p => AnnotateCrosslink(p, "1896")).ToList();
            var betaPeptides = GenerateUniquePeptides(500, 15).Select(p => AnnotateCrosslink(p, "1896")).ToList();
            var modelInputs = alphaPeptides.Zip(betaPeptides, (alpha, beta) => new CrosslinkIntensityPredictionInput(alpha, beta, 2, 35)).ToList();
            var model = new Prosit2023IntensityXLCMS2();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(500));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict 500 crosslinked peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }

        [Test]
        [Explicit("Massive test, takes a long time to run")]
        [Category("Performance Benchmark")]
        /// <summary>
        /// Performance benchmark for Prosit_2023_intensity_XL_CMS3.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500 unique crosslinked peptide pairs of length 15, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500 peptides.
        ///  - Beta peptides are reversed alpha peptides.
        ///  - Uses fixed NCE of 35 (no collision energy input required).
        ///  - CMS3 helper does not require the crosslink annotation, but annotating is harmless.
        /// </summary>
        public static void BenchmarkProsit2023IntensityXLCMS3()
        {
            var alphaPeptides = GenerateUniquePeptides(500, 15).Select(p => AnnotateCrosslink(p, "1881")).ToList();
            var betaPeptides = GenerateUniquePeptides(500, 15).Select(p => AnnotateCrosslink(p, "1881")).ToList();
            var modelInputs = alphaPeptides.Zip(betaPeptides, (alpha, beta) => new CrosslinkIntensityPredictionInput(alpha, beta, 2, 35)).ToList();
            var model = new Prosit2023IntensityXLCMS3();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(500));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict 500 crosslinked peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }

        [Test]
        [Explicit("Massive test, takes a long time to run")]
        [Category("Performance Benchmark")]
        /// <summary>
        /// Performance benchmark for Prosit_2024_intensity_XL_NMS2.
        /// 
        /// Notes: 
        ///  - MaxBatchSize: 1000 peptides, Predict() takes ____ per batch.
        ///  - With 500 unique crosslinked peptide pairs of length 15, test takes approximately ____.
        ///  - Timing is linearly proportional to peptide count due to throttled approach.
        ///  - Default MaxNumberOfBatchesPerRequest: 250, tested up to 500 peptides.
        ///  - Beta peptides are reversed alpha peptides.
        ///  - Alpha peptides include K[UNIMOD:1898] for the NMS2 helper to locate crosslink.
        /// </summary>
        public static void BenchmarkProsit2024IntensityXLNMS2()
        {
            var alphaPeptides = GenerateUniquePeptides(500, 15).Select(p => AnnotateCrosslink(p, "1898")).ToList();
            var betaPeptides = GenerateUniquePeptides(500, 15).Select(p => AnnotateCrosslink(p, "1898")).ToList();
            var modelInputs = alphaPeptides.Zip(betaPeptides, (alpha, beta) => new CrosslinkIntensityPredictionInput(alpha, beta, 2, 35)).ToList();
            var model = new Prosit2024IntensityXLNMS2();
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.DoesNotThrow(() => model.Predict(modelInputs));
            watch.Stop();
            Assert.That(model.Predictions.Count, Is.EqualTo(500));
            Assert.That(model.ValidInputsMask, Is.All.True);
            Console.WriteLine($"Time taken to predict 500 crosslinked peptides: {watch.Elapsed.Minutes}min {watch.Elapsed.Seconds}s {watch.Elapsed.Milliseconds}ms");
        }
    }
}
