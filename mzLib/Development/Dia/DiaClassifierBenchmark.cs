// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Phase 14 Benchmark: Gradient-Boosted Tree vs Linear Discriminant
// Location: mzLib/Development/Dia/FDR/DiaClassifierBenchmark.cs
//
// Run from Development/Program.cs:
//   DiaClassifierBenchmark.RunAll();

using System;
using System.Buffers;
using System.Diagnostics;

namespace MassSpectrometry.Dia.Benchmarks
{
    public static class DiaClassifierBenchmark
    {
        public static void RunAll()
        {
            Console.WriteLine("=" + new string('=', 69));
            Console.WriteLine("Phase 14 Benchmark: GBT vs LDA Classifier Comparison");
            Console.WriteLine("=" + new string('=', 69));
            Console.WriteLine();

            BenchmarkLinearProblem();
            BenchmarkXorProblem();
            BenchmarkConcentricCircles();
            BenchmarkRealisticDiaFeatures();
            BenchmarkTrainingThroughput();
            BenchmarkPredictionThroughput();
            BenchmarkScalingBehavior();

            Console.WriteLine();
            Console.WriteLine("All benchmarks complete.");
        }

        // ════════════════════════════════════════════════════════════════
        //  AUC Comparison Benchmarks
        // ════════════════════════════════════════════════════════════════

        private static void BenchmarkLinearProblem()
        {
            Console.WriteLine("--- Linear Separation (sanity check) ---");
            var rng = new Random(42);
            int n = 5000;
            int fc = 10;
            var (X, y) = GenerateLinearProblem(rng, n, fc, separation: 3f);
            ShuffleData(X, y, n, fc, rng);

            float gbtAuc = TrainAndEvaluateGbt(X, y, n, fc);
            float ldaAuc = TrainAndEvaluateLda(X, y, n, fc);

            Console.WriteLine($"  GBT AUC: {gbtAuc:F4}");
            Console.WriteLine($"  LDA AUC: {ldaAuc:F4}");
            Console.WriteLine($"  Delta:   {gbtAuc - ldaAuc:+0.0000;-0.0000}");
            Console.WriteLine();
        }

        private static void BenchmarkXorProblem()
        {
            Console.WriteLine("--- XOR Problem (non-linear, 2 features) ---");
            var rng = new Random(42);
            int n = 4000;
            int fc = 2;
            var X = new float[n * fc];
            var y = new float[n];

            for (int i = 0; i < n; i++)
            {
                float x1 = (float)(rng.NextDouble() * 4 - 2);
                float x2 = (float)(rng.NextDouble() * 4 - 2);
                X[i * 2] = x1;
                X[i * 2 + 1] = x2;
                y[i] = (x1 > 0) != (x2 > 0) ? 1f : 0f;
            }
            // XOR data is already interleaved by construction — no shuffle needed

            float gbtAuc = TrainAndEvaluateGbt(X, y, n, fc);
            float ldaAuc = TrainAndEvaluateLda(X, y, n, fc);

            Console.WriteLine($"  GBT AUC: {gbtAuc:F4}");
            Console.WriteLine($"  LDA AUC: {ldaAuc:F4}");
            Console.WriteLine($"  Delta:   {gbtAuc - ldaAuc:+0.0000;-0.0000}  ← GBT advantage on non-linear");
            Console.WriteLine();
        }

        private static void BenchmarkConcentricCircles()
        {
            Console.WriteLine("--- Concentric Circles (non-linear, 2 features) ---");
            var rng = new Random(42);
            int n = 4000;
            int fc = 2;
            var X = new float[n * fc];
            var y = new float[n];

            for (int i = 0; i < n; i++)
            {
                float angle = (float)(rng.NextDouble() * 2 * Math.PI);
                bool isTarget = i < n / 2;
                float radius = isTarget ? 1f : 3f;
                radius += (float)(rng.NextDouble() * 0.5 - 0.25);
                X[i * 2] = radius * MathF.Cos(angle);
                X[i * 2 + 1] = radius * MathF.Sin(angle);
                y[i] = isTarget ? 1f : 0f;
            }
            ShuffleData(X, y, n, fc, rng);

            float gbtAuc = TrainAndEvaluateGbt(X, y, n, fc);
            float ldaAuc = TrainAndEvaluateLda(X, y, n, fc);

            Console.WriteLine($"  GBT AUC: {gbtAuc:F4}");
            Console.WriteLine($"  LDA AUC: {ldaAuc:F4}");
            Console.WriteLine($"  Delta:   {gbtAuc - ldaAuc:+0.0000;-0.0000}");
            Console.WriteLine();
        }

        private static void BenchmarkRealisticDiaFeatures()
        {
            Console.WriteLine("--- Realistic DIA Features (28 features, non-linear interactions) ---");
            var rng = new Random(42);
            int n = 50000;
            int fc = 28;
            int nTargets = 39000;
            var X = new float[n * fc];
            var y = new float[n];

            for (int i = 0; i < n; i++)
            {
                bool isTarget = i < nTargets;
                y[i] = isTarget ? 1f : 0f;
                int offset = i * fc;

                bool isGoodTarget = isTarget && rng.NextDouble() < 0.6;
                float baseQ = isGoodTarget ? 0.8f + Noise(rng, 0.1f) :
                              isTarget ? 0.4f + Noise(rng, 0.2f) :
                              Noise(rng, 0.2f) + 0.15f;

                X[offset + 0] = baseQ * 0.85f + Noise(rng, 0.15f);
                X[offset + 1] = X[offset] * 0.98f + Noise(rng, 0.02f);
                X[offset + 2] = X[offset] * 0.97f + Noise(rng, 0.03f);
                X[offset + 3] = MathF.Min(baseQ + 0.2f + Noise(rng, 0.1f), 1f);
                X[offset + 4] = (isGoodTarget ? 8f : isTarget ? 6f : 4f) + Noise(rng, 1.5f);
                X[offset + 5] = baseQ * 0.9f + Noise(rng, 0.15f);
                X[offset + 6] = baseQ * 0.6f + Noise(rng, 0.2f);
                X[offset + 7] = (isTarget ? 0.3f : 0.6f) + Noise(rng, 0.2f);
                X[offset + 8] = (isTarget ? 15f : 8f) + Noise(rng, 5f);
                X[offset + 9] = (isTarget ? 30f : 15f) + Noise(rng, 10f);
                X[offset + 10] = 0.3f + Noise(rng, 0.15f);

                float rtDev = isGoodTarget ? MathF.Abs(Noise(rng, 0.3f)) :
                              isTarget ? MathF.Abs(Noise(rng, 1.0f)) :
                              MathF.Abs(Noise(rng, 2.0f));
                X[offset + 11] = rtDev;
                X[offset + 12] = rtDev * rtDev;

                X[offset + 13] = baseQ * 0.8f + Noise(rng, 0.15f);
                X[offset + 14] = baseQ * 0.7f + Noise(rng, 0.15f);
                X[offset + 15] = (isTarget ? 1.2f : 2.5f) + Noise(rng, 0.8f);
                X[offset + 16] = baseQ * 5f + Noise(rng, 0.5f);
                X[offset + 17] = baseQ * 0.85f + Noise(rng, 0.15f);
                X[offset + 18] = baseQ * 0.5f + Noise(rng, 0.2f);
                X[offset + 19] = (isTarget ? 0.15f : 0.35f) + Noise(rng, 0.1f);
                X[offset + 20] = baseQ * 0.8f + Noise(rng, 0.1f);
                X[offset + 21] = (isTarget ? 0.3f : 0.8f) + Noise(rng, 0.2f);
                X[offset + 22] = (isTarget ? 0.6f : 1.5f) + Noise(rng, 0.3f);
                X[offset + 23] = (isTarget ? 0.2f : 0.5f) + Noise(rng, 0.15f);
                X[offset + 24] = baseQ * 0.75f + Noise(rng, 0.15f);
                X[offset + 25] = (isTarget ? 3f : 1f) + Noise(rng, 1.5f);
                X[offset + 26] = (isTarget ? 0.15f : 0.25f) + Noise(rng, 0.1f);
                X[offset + 27] = (isGoodTarget ? 2.5f : 1.5f) + Noise(rng, 0.5f);
            }
            ShuffleData(X, y, n, fc, rng);

            var sw = Stopwatch.StartNew();
            float gbtAuc = TrainAndEvaluateGbt(X, y, n, fc);
            long gbtMs = sw.ElapsedMilliseconds;

            sw.Restart();
            float ldaAuc = TrainAndEvaluateLda(X, y, n, fc);
            long ldaMs = sw.ElapsedMilliseconds;

            Console.WriteLine($"  Samples:    {n:N0} ({nTargets:N0} targets, {n - nTargets:N0} decoys)");
            Console.WriteLine($"  Features:   {fc}");
            Console.WriteLine($"  GBT AUC:    {gbtAuc:F4}  ({gbtMs}ms train+eval)");
            Console.WriteLine($"  LDA AUC:    {ldaAuc:F4}  ({ldaMs}ms train+eval)");
            Console.WriteLine($"  Delta:      {gbtAuc - ldaAuc:+0.0000;-0.0000}  ← expected improvement");
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  Throughput Benchmarks
        // ════════════════════════════════════════════════════════════════

        private static void BenchmarkTrainingThroughput()
        {
            Console.WriteLine("--- Training Throughput ---");
            int n = 50000;
            int fc = 28;
            var rng = new Random(42);
            var (X, y) = GenerateLinearProblem(rng, n, fc, separation: 2f);
            ShuffleData(X, y, n, fc, rng);

            // Warm up
            var warmup = new DiaGradientBoostedClassifier(fc, numTrees: 5, maxDepth: 2);
            warmup.Train(X, y, n);

            int runs = 3;
            long gbtTotalMs = 0;
            for (int r = 0; r < runs; r++)
            {
                var gbt = new DiaGradientBoostedClassifier(fc, numTrees: 100, maxDepth: 4,
                    learningRate: 0.1f, minSamplesLeaf: 20, numBins: 64, subsampleFraction: 0.8f);
                var sw = Stopwatch.StartNew();
                gbt.Train(X, y, n);
                gbtTotalMs += sw.ElapsedMilliseconds;
            }

            float gbtAvg = gbtTotalMs / (float)runs;
            Console.WriteLine($"  Dataset: {n:N0} samples × {fc} features");
            Console.WriteLine($"  GBT (100 trees, depth 4): {gbtAvg:F1}ms  ({n / Math.Max(gbtAvg / 1000, 0.001):N0} samples/sec)");
            Console.WriteLine();
        }

        private static void BenchmarkPredictionThroughput()
        {
            Console.WriteLine("--- Prediction Throughput ---");
            int n = 50000;
            int fc = 28;
            var rng = new Random(42);
            var (X, y) = GenerateLinearProblem(rng, n, fc, separation: 2f);
            ShuffleData(X, y, n, fc, rng);

            var gbt = new DiaGradientBoostedClassifier(fc, numTrees: 100, maxDepth: 4);
            gbt.Train(X, y, n);

            var scores = new float[n];
            gbt.PredictBatch(X, n, scores);

            int runs = 10;
            GC.Collect();
            GC.WaitForPendingFinalizers();

            var sw = Stopwatch.StartNew();
            for (int r = 0; r < runs; r++)
                gbt.PredictBatch(X, n, scores);
            sw.Stop();
            float gbtPredMs = sw.ElapsedMilliseconds / (float)runs;

            Console.WriteLine($"  Dataset: {n:N0} samples × {fc} features");
            Console.WriteLine($"  GBT predict: {gbtPredMs:F2}ms  ({n / Math.Max(gbtPredMs / 1000, 0.001):N0} samples/sec)");
            Console.WriteLine($"  Context: FDR runs ~5-7 iterations × {n:N0} predictions");
            Console.WriteLine($"  Est. total FDR prediction time (7 iters): {gbtPredMs * 7:F0}ms");
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  Scaling Benchmark
        // ════════════════════════════════════════════════════════════════

        private static void BenchmarkScalingBehavior()
        {
            Console.WriteLine("--- Scaling: Training Time vs Sample Count ---");
            int fc = 28;
            int[] sampleCounts = { 5000, 10000, 25000, 50000, 100000 };

            Console.WriteLine($"  {"Samples",-12} {"GBT (ms)",-12} {"GBT AUC",-10}");

            foreach (int n in sampleCounts)
            {
                var rng = new Random(42);
                var X = new float[n * fc];
                var y = new float[n];
                for (int i = 0; i < n; i++)
                {
                    y[i] = i < n / 2 ? 1f : 0f;
                    for (int f = 0; f < fc; f++)
                    {
                        float signal = y[i] > 0.5f ? 1f : -1f;
                        X[i * fc + f] = (f < 2)
                            ? signal + (float)(rng.NextDouble() * 2 - 1)
                            : (float)(rng.NextDouble() * 2 - 1);
                    }
                    // Non-linear interaction
                    if (y[i] > 0.5f && X[i * fc] * X[i * fc + 1] < 0)
                        y[i] = rng.NextDouble() < 0.3 ? 0f : 1f;
                }
                ShuffleData(X, y, n, fc, rng);

                var gbt = new DiaGradientBoostedClassifier(fc, numTrees: 100, maxDepth: 4,
                    minSamplesLeaf: 20, numBins: 64);
                var sw = Stopwatch.StartNew();
                gbt.Train(X, y, n);
                long gbtMs = sw.ElapsedMilliseconds;

                // Evaluate on held-out 20%
                int trainN = (int)(n * 0.8);
                int testN = n - trainN;
                // Since data is shuffled, we can just split
                float gbtAuc = gbt.ComputeAuc(
                    X.AsSpan(trainN * fc, testN * fc),
                    y.AsSpan(trainN, testN), testN);

                Console.WriteLine($"  {n,-12:N0} {gbtMs,-12} {gbtAuc,-10:F4}");
            }

            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  Helpers
        // ════════════════════════════════════════════════════════════════

        private static (float[] X, float[] y) GenerateLinearProblem(
            Random rng, int n, int fc, float separation)
        {
            var X = new float[n * fc];
            var y = new float[n];
            for (int i = 0; i < n; i++)
            {
                y[i] = i < n / 2 ? 1f : 0f;
                for (int f = 0; f < fc; f++)
                {
                    float center = y[i] > 0.5f ? separation / 2f : -separation / 2f;
                    X[i * fc + f] = center + (float)(rng.NextDouble() * 2 - 1);
                }
            }
            return (X, y);
        }

        /// <summary>
        /// Fisher-Yates shuffle of rows in the feature matrix + label array.
        /// Critical: without this, an 80/20 split on ordered data gives a test set
        /// with only one class, yielding AUC = 0.5 (undefined).
        /// </summary>
        private static void ShuffleData(float[] X, float[] y, int n, int fc, Random rng)
        {
            for (int i = n - 1; i > 0; i--)
            {
                int j = rng.Next(i + 1);
                if (i == j) continue;

                // Swap labels
                (y[i], y[j]) = (y[j], y[i]);

                // Swap feature rows
                int rowI = i * fc;
                int rowJ = j * fc;
                for (int f = 0; f < fc; f++)
                    (X[rowI + f], X[rowJ + f]) = (X[rowJ + f], X[rowI + f]);
            }
        }

        /// <summary>
        /// Trains GBT on first 80% of (shuffled) data, evaluates AUC on remaining 20%.
        /// </summary>
        private static float TrainAndEvaluateGbt(float[] X, float[] y, int n, int fc)
        {
            int trainN = (int)(n * 0.8);
            int testN = n - trainN;

            var gbt = new DiaGradientBoostedClassifier(fc, numTrees: 100, maxDepth: 4,
                learningRate: 0.1f, minSamplesLeaf: 20, numBins: 64, subsampleFraction: 0.8f);
            gbt.Train(X.AsSpan(0, trainN * fc), y.AsSpan(0, trainN), trainN);
            return gbt.ComputeAuc(
                X.AsSpan(trainN * fc, testN * fc),
                y.AsSpan(trainN, testN), testN);
        }

        /// <summary>
        /// Trains a diagonal Fisher LDA on first 80%, evaluates on remaining 20%.
        /// </summary>
        private static float TrainAndEvaluateLda(float[] X, float[] y, int n, int fc)
        {
            int trainN = (int)(n * 0.8);
            int testN = n - trainN;

            var mean0 = new float[fc];
            var mean1 = new float[fc];
            int count0 = 0, count1 = 0;

            for (int i = 0; i < trainN; i++)
            {
                int offset = i * fc;
                if (y[i] > 0.5f)
                {
                    count1++;
                    for (int f = 0; f < fc; f++) mean1[f] += X[offset + f];
                }
                else
                {
                    count0++;
                    for (int f = 0; f < fc; f++) mean0[f] += X[offset + f];
                }
            }

            if (count0 == 0 || count1 == 0)
                return 0.5f; // degenerate

            for (int f = 0; f < fc; f++) { mean0[f] /= count0; mean1[f] /= count1; }

            var variance = new float[fc];
            for (int i = 0; i < trainN; i++)
            {
                int offset = i * fc;
                var mean = y[i] > 0.5f ? mean1 : mean0;
                for (int f = 0; f < fc; f++)
                {
                    float diff = X[offset + f] - mean[f];
                    variance[f] += diff * diff;
                }
            }

            var weights = new float[fc];
            float bias = 0f;
            for (int f = 0; f < fc; f++)
            {
                variance[f] /= trainN;
                variance[f] = MathF.Max(variance[f], 1e-8f);
                weights[f] = (mean1[f] - mean0[f]) / variance[f];
                bias -= 0.5f * weights[f] * (mean1[f] + mean0[f]);
            }

            var scores = new float[testN];
            for (int i = 0; i < testN; i++)
            {
                int offset = (trainN + i) * fc;
                float logit = bias;
                for (int f = 0; f < fc; f++) logit += weights[f] * X[offset + f];
                scores[i] = 1f / (1f + MathF.Exp(-logit));
            }

            return DiaGradientBoostedClassifier.CalculateAuc(
                            scores.AsSpan(), y.AsSpan(trainN, testN), testN);
        }

        private static float Noise(Random rng, float scale) =>
            (float)(rng.NextDouble() * 2 - 1) * scale;
    }
}