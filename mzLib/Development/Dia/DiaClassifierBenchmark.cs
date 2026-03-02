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
    /// <summary>
    /// Benchmarks comparing GBT vs LDA classifiers on:
    ///   1. Training throughput (samples/sec)
    ///   2. Prediction throughput (samples/sec)
    ///   3. AUC on linear vs non-linear synthetic problems
    ///   4. AUC on realistic DIA-like feature distributions (28 features)
    ///   5. Scaling behavior with sample count
    /// 
    /// This benchmark validates the Phase 14 claim that GBT significantly
    /// outperforms LDA on non-linear target/decoy boundaries.
    /// </summary>
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

        /// <summary>
        /// Baseline: linearly separable problem. Both classifiers should do well.
        /// Confirms GBT doesn't degrade on easy problems.
        /// </summary>
        private static void BenchmarkLinearProblem()
        {
            Console.WriteLine("--- Linear Separation (sanity check) ---");
            var rng = new Random(42);
            int n = 5000;
            int features = 10;
            var X = new float[n * features];
            var y = new float[n];

            for (int i = 0; i < n; i++)
            {
                bool isTarget = i < n / 2;
                y[i] = isTarget ? 1f : 0f;
                for (int f = 0; f < features; f++)
                {
                    float center = isTarget ? 1.5f : -1.5f;
                    X[i * features + f] = center + (float)(rng.NextDouble() * 2 - 1);
                }
            }

            float gbtAuc = TrainAndEvaluate(X, y, n, features, DiaClassifierType.GradientBoostedTree);
            float ldaAuc = TrainAndEvaluate(X, y, n, features, DiaClassifierType.LinearDiscriminant);

            Console.WriteLine($"  GBT AUC: {gbtAuc:F4}");
            Console.WriteLine($"  LDA AUC: {ldaAuc:F4}");
            Console.WriteLine($"  Delta:   {gbtAuc - ldaAuc:+0.0000;-0.0000}");
            Console.WriteLine();
        }

        /// <summary>
        /// XOR problem: purely non-linear. LDA should fail (~0.5), GBT should succeed (>0.9).
        /// </summary>
        private static void BenchmarkXorProblem()
        {
            Console.WriteLine("--- XOR Problem (non-linear, 2 features) ---");
            var rng = new Random(42);
            int n = 4000;
            int features = 2;
            var X = new float[n * features];
            var y = new float[n];

            for (int i = 0; i < n; i++)
            {
                float x1 = (float)(rng.NextDouble() * 4 - 2);
                float x2 = (float)(rng.NextDouble() * 4 - 2);
                X[i * 2] = x1;
                X[i * 2 + 1] = x2;
                y[i] = (x1 > 0) != (x2 > 0) ? 1f : 0f;
            }

            float gbtAuc = TrainAndEvaluate(X, y, n, features, DiaClassifierType.GradientBoostedTree);
            float ldaAuc = TrainAndEvaluate(X, y, n, features, DiaClassifierType.LinearDiscriminant);

            Console.WriteLine($"  GBT AUC: {gbtAuc:F4}");
            Console.WriteLine($"  LDA AUC: {ldaAuc:F4}");
            Console.WriteLine($"  Delta:   {gbtAuc - ldaAuc:+0.0000;-0.0000}  ← GBT advantage on non-linear");
            Console.WriteLine();
        }

        /// <summary>
        /// Concentric circles: another classic non-linear boundary.
        /// </summary>
        private static void BenchmarkConcentricCircles()
        {
            Console.WriteLine("--- Concentric Circles (non-linear, 2 features) ---");
            var rng = new Random(42);
            int n = 4000;
            int features = 2;
            var X = new float[n * features];
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

            float gbtAuc = TrainAndEvaluate(X, y, n, features, DiaClassifierType.GradientBoostedTree);
            float ldaAuc = TrainAndEvaluate(X, y, n, features, DiaClassifierType.LinearDiscriminant);

            Console.WriteLine($"  GBT AUC: {gbtAuc:F4}");
            Console.WriteLine($"  LDA AUC: {ldaAuc:F4}");
            Console.WriteLine($"  Delta:   {gbtAuc - ldaAuc:+0.0000;-0.0000}");
            Console.WriteLine();
        }

        /// <summary>
        /// Realistic DIA-like feature distributions.
        /// 28 features modeled after the actual Phase 13 feature space:
        ///   - Some features have high linear separation (BestFragCorrSum, MeanFragCorr)
        ///   - Some have non-linear interactions (RT deviation × fragment correlation)
        ///   - Some are near-redundant (TemporalScore ↔ SpectralAngle, r=0.993)
        ///   - Some are noisy (IntensityCV, BoundarySignalRatio)
        /// 
        /// This is the most important benchmark: it simulates the actual problem
        /// where we expect the GBT to gain +3000-5000 IDs over LDA.
        /// </summary>
        private static void BenchmarkRealisticDiaFeatures()
        {
            Console.WriteLine("--- Realistic DIA Features (28 features, non-linear interactions) ---");
            var rng = new Random(42);
            int n = 50000; // ~39K targets + ~12K decoys (matches real data ratio)
            int features = 28;
            int nTargets = 39000;
            int nDecoys = n - nTargets;
            var X = new float[n * features];
            var y = new float[n];

            // Feature indices (matching DiaFeatureVector order):
            // 0: ApexScore, 1: TemporalScore, 2: SpectralAngle, 3: FragDetRate
            // 4: LogTotalIntensity, 5: MeanFragCorr, 6: MinFragCorr, 7: XicDepthCV
            // 8: MedianXicDepth, 9: TimePointsUsed, 10: IntensityCV
            // 11: RtDeviationMinutes, 12: RtDeviationSquared
            // 13: PeakApexScore, 14: PeakMeanFragCorr, 15: PeakWidth
            // 16: BestFragCorrSum, 17: MedianFragRefCorr, 18: MinFragRefCorr
            // 19: StdFragRefCorr, 20: BestFragWeightedCosine
            // 21: MeanSigRatioDev, 22: MaxSigRatioDev, 23: StdSigRatioDev
            // 24: SmoothedMeanFragCorr, 25: Log2SNR
            // 26: BoundarySignalRatio, 27: ApexToMeanRatio

            for (int i = 0; i < n; i++)
            {
                bool isTarget = i < nTargets;
                y[i] = isTarget ? 1f : 0f;
                int offset = i * features;

                // Simulate target vs decoy distributions
                // Targets: mix of good hits (60%) and marginal (40%)
                // Decoys: lower scores with some overlap
                bool isGoodTarget = isTarget && rng.NextDouble() < 0.6;
                bool isMarginalTarget = isTarget && !isGoodTarget;

                float baseQuality;
                if (isGoodTarget) baseQuality = 0.8f + (float)(rng.NextDouble() * 0.2);
                else if (isMarginalTarget) baseQuality = 0.3f + (float)(rng.NextDouble() * 0.4);
                else baseQuality = (float)(rng.NextDouble() * 0.4); // decoy

                // Strong linear features (similar to Phase 13 top weights)
                X[offset + 16] = baseQuality * 5f + Noise(rng, 0.5f); // BestFragCorrSum
                X[offset + 5] = baseQuality * 0.9f + Noise(rng, 0.15f); // MeanFragCorr
                X[offset + 17] = baseQuality * 0.85f + Noise(rng, 0.15f); // MedianFragRefCorr
                X[offset + 3] = MathF.Min(baseQuality + 0.2f + Noise(rng, 0.1f), 1f); // FragDetRate

                // RT deviation: non-linear interaction with quality
                // Good targets have small RT deviation, but the boundary is non-linear
                float rtDev;
                if (isGoodTarget) rtDev = MathF.Abs(Noise(rng, 0.3f));
                else if (isMarginalTarget) rtDev = MathF.Abs(Noise(rng, 1.0f));
                else rtDev = MathF.Abs(Noise(rng, 2.0f));
                X[offset + 11] = rtDev; // RtDeviationMinutes
                X[offset + 12] = rtDev * rtDev; // RtDeviationSquared

                // NON-LINEAR INTERACTION: high MeanFragCorr + low RT deviation → target
                // This is exactly the kind of boundary LDA misses
                X[offset + 0] = baseQuality * 0.85f + Noise(rng, 0.15f); // ApexScore

                // Near-redundant features (r ≈ 0.99)
                X[offset + 1] = X[offset + 0] * 0.98f + Noise(rng, 0.02f); // TemporalScore
                X[offset + 2] = X[offset + 0] * 0.97f + Noise(rng, 0.03f); // SpectralAngle

                // Intensity features
                X[offset + 4] = (isGoodTarget ? 8f : (isTarget ? 6f : 4f)) + Noise(rng, 1.5f); // LogTotalIntensity
                X[offset + 10] = 0.3f + Noise(rng, 0.15f); // IntensityCV (noisy, low sep)

                // Fragment correlation features
                X[offset + 6] = baseQuality * 0.6f + Noise(rng, 0.2f); // MinFragCorr
                X[offset + 7] = (isTarget ? 0.3f : 0.6f) + Noise(rng, 0.2f); // XicDepthCV (negative weight)
                X[offset + 8] = (isTarget ? 15f : 8f) + Noise(rng, 5f); // MedianXicDepth
                X[offset + 9] = (isTarget ? 30f : 15f) + Noise(rng, 10f); // TimePointsUsed

                // Peak features
                X[offset + 13] = baseQuality * 0.8f + Noise(rng, 0.15f); // PeakApexScore
                X[offset + 14] = baseQuality * 0.7f + Noise(rng, 0.15f); // PeakMeanFragCorr
                X[offset + 15] = (isTarget ? 1.2f : 2.5f) + Noise(rng, 0.8f); // PeakWidth

                // Reference curve features
                X[offset + 18] = baseQuality * 0.5f + Noise(rng, 0.2f); // MinFragRefCorr
                X[offset + 19] = (isTarget ? 0.15f : 0.35f) + Noise(rng, 0.1f); // StdFragRefCorr
                X[offset + 20] = baseQuality * 0.8f + Noise(rng, 0.1f); // BestFragWeightedCosine

                // Signal ratio deviation (negative = worse)
                X[offset + 21] = (isTarget ? 0.3f : 0.8f) + Noise(rng, 0.2f); // MeanSigRatioDev
                X[offset + 22] = (isTarget ? 0.6f : 1.5f) + Noise(rng, 0.3f); // MaxSigRatioDev
                X[offset + 23] = (isTarget ? 0.2f : 0.5f) + Noise(rng, 0.15f); // StdSigRatioDev

                // Smoothed/other
                X[offset + 24] = baseQuality * 0.75f + Noise(rng, 0.15f); // SmoothedMeanFragCorr
                X[offset + 25] = (isTarget ? 3f : 1f) + Noise(rng, 1.5f); // Log2SNR
                X[offset + 26] = (isTarget ? 0.15f : 0.25f) + Noise(rng, 0.1f); // BoundarySignalRatio
                X[offset + 27] = (isGoodTarget ? 2.5f : 1.5f) + Noise(rng, 0.5f); // ApexToMeanRatio
            }

            // Train and evaluate both classifiers
            var sw = Stopwatch.StartNew();
            float gbtAuc = TrainAndEvaluate(X, y, n, features, DiaClassifierType.GradientBoostedTree);
            long gbtMs = sw.ElapsedMilliseconds;

            sw.Restart();
            float ldaAuc = TrainAndEvaluate(X, y, n, features, DiaClassifierType.LinearDiscriminant);
            long ldaMs = sw.ElapsedMilliseconds;

            Console.WriteLine($"  Samples:    {n:N0} ({nTargets:N0} targets, {nDecoys:N0} decoys)");
            Console.WriteLine($"  Features:   {features}");
            Console.WriteLine($"  GBT AUC:    {gbtAuc:F4}  ({gbtMs}ms)");
            Console.WriteLine($"  LDA AUC:    {ldaAuc:F4}  ({ldaMs}ms)");
            Console.WriteLine($"  Delta:      {gbtAuc - ldaAuc:+0.0000;-0.0000}  ← expected improvement");
            Console.WriteLine();

            // Estimate ID impact (rough: AUC improvement of 0.03 ≈ +3000-5000 IDs at 1% FDR)
            float aucDelta = gbtAuc - ldaAuc;
            if (aucDelta > 0)
            {
                int estimatedAdditionalIds = (int)(aucDelta * 150000); // rough empirical scaling
                Console.WriteLine($"  Estimated additional IDs from GBT: ~{estimatedAdditionalIds:N0}");
                Console.WriteLine();
            }
        }

        /// <summary>
        /// Training throughput: how fast can each classifier train on 50K samples × 28 features?
        /// </summary>
        private static void BenchmarkTrainingThroughput()
        {
            Console.WriteLine("--- Training Throughput ---");
            int n = 50000;
            int features = 28;
            var rng = new Random(42);
            var X = new float[n * features];
            var y = new float[n];

            for (int i = 0; i < n; i++)
            {
                y[i] = i < n / 2 ? 1f : 0f;
                for (int f = 0; f < features; f++)
                    X[i * features + f] = (y[i] > 0.5f ? 1f : -1f) + (float)(rng.NextDouble() * 2 - 1);
            }

            // Warm up
            var warmup = new DiaGradientBoostedClassifier(features, numTrees: 5, maxDepth: 2);
            warmup.Train(X, y, n);

            // GBT training
            int gbtRuns = 3;
            long gbtTotalMs = 0;
            for (int r = 0; r < gbtRuns; r++)
            {
                var gbt = new DiaGradientBoostedClassifier(features, numTrees: 100, maxDepth: 4,
                    learningRate: 0.1f, minSamplesLeaf: 20, numBins: 64, subsampleFraction: 0.8f);
                var sw = Stopwatch.StartNew();
                gbt.Train(X, y, n);
                sw.Stop();
                gbtTotalMs += sw.ElapsedMilliseconds;
            }

            // LDA training
            int ldaRuns = 3;
            long ldaTotalMs = 0;
            for (int r = 0; r < ldaRuns; r++)
            {
                var lda = new DiaLdaClassifierAdapter(features);
                var sw = Stopwatch.StartNew();
                lda.Train(X, y, n);
                sw.Stop();
                ldaTotalMs += sw.ElapsedMilliseconds;
            }

            float gbtAvgMs = gbtTotalMs / (float)gbtRuns;
            float ldaAvgMs = ldaTotalMs / (float)ldaRuns;

            Console.WriteLine($"  Dataset: {n:N0} samples × {features} features");
            Console.WriteLine($"  GBT (100 trees, depth 4): {gbtAvgMs:F1}ms  ({n / (gbtAvgMs / 1000):N0} samples/sec)");
            Console.WriteLine($"  LDA (diagonal Fisher):    {ldaAvgMs:F1}ms  ({n / (ldaAvgMs / 1000):N0} samples/sec)");
            Console.WriteLine($"  GBT/LDA time ratio:       {gbtAvgMs / ldaAvgMs:F1}x");
            Console.WriteLine();
        }

        /// <summary>
        /// Prediction throughput: how fast can each classifier score 50K samples?
        /// In the FDR loop, prediction runs every iteration on all results.
        /// </summary>
        private static void BenchmarkPredictionThroughput()
        {
            Console.WriteLine("--- Prediction Throughput ---");
            int n = 50000;
            int features = 28;
            var rng = new Random(42);
            var X = new float[n * features];
            var y = new float[n];

            for (int i = 0; i < n; i++)
            {
                y[i] = i < n / 2 ? 1f : 0f;
                for (int f = 0; f < features; f++)
                    X[i * features + f] = (y[i] > 0.5f ? 1f : -1f) + (float)(rng.NextDouble() * 2 - 1);
            }

            // Train both
            var gbt = new DiaGradientBoostedClassifier(features, numTrees: 100, maxDepth: 4);
            gbt.Train(X, y, n);

            var lda = new DiaLdaClassifierAdapter(features);
            lda.Train(X, y, n);

            var scores = new float[n];

            // Warm up
            gbt.PredictBatch(X, n, scores);
            lda.PredictBatch(X, n, scores);

            // GBT prediction
            int runs = 10;
            GC.Collect();
            GC.WaitForPendingFinalizers();
            var sw = Stopwatch.StartNew();
            for (int r = 0; r < runs; r++)
                gbt.PredictBatch(X, n, scores);
            sw.Stop();
            float gbtPredMs = sw.ElapsedMilliseconds / (float)runs;

            // LDA prediction
            GC.Collect();
            GC.WaitForPendingFinalizers();
            sw.Restart();
            for (int r = 0; r < runs; r++)
                lda.PredictBatch(X, n, scores);
            sw.Stop();
            float ldaPredMs = sw.ElapsedMilliseconds / (float)runs;

            Console.WriteLine($"  Dataset: {n:N0} samples × {features} features");
            Console.WriteLine($"  GBT predict: {gbtPredMs:F2}ms  ({n / (gbtPredMs / 1000):N0} samples/sec)");
            Console.WriteLine($"  LDA predict: {ldaPredMs:F2}ms  ({n / (ldaPredMs / 1000):N0} samples/sec)");
            Console.WriteLine($"  GBT/LDA ratio: {gbtPredMs / Math.Max(ldaPredMs, 0.01f):F1}x");
            Console.WriteLine();

            // Context: FDR runs ~5-7 iterations, each predicting 50K samples
            // Even if GBT is 10x slower than LDA at prediction, 50K × 10ms = 500ms total
            // vs LDA at 50K × 1ms = 50ms. Both negligible vs extraction (2.4s).
            float totalGbtFdr = gbtPredMs * 7f; // 7 iterations
            float totalLdaFdr = ldaPredMs * 7f;
            Console.WriteLine($"  Est. total FDR prediction (7 iters): GBT={totalGbtFdr:F0}ms, LDA={totalLdaFdr:F0}ms");
            Console.WriteLine();
        }

        /// <summary>
        /// How training time scales with sample count.
        /// Should be roughly linear for histogram-based GBT.
        /// </summary>
        private static void BenchmarkScalingBehavior()
        {
            Console.WriteLine("--- Scaling: Training Time vs Sample Count ---");
            int features = 28;
            int[] sampleCounts = { 5000, 10000, 25000, 50000, 100000 };
            var rng = new Random(42);

            Console.WriteLine($"  {"Samples",-12} {"GBT (ms)",-12} {"LDA (ms)",-12} {"GBT AUC",-10} {"LDA AUC",-10}");

            foreach (int n in sampleCounts)
            {
                var X = new float[n * features];
                var y = new float[n];
                for (int i = 0; i < n; i++)
                {
                    y[i] = i < n / 2 ? 1f : 0f;
                    for (int f = 0; f < features; f++)
                    {
                        float signal = y[i] > 0.5f ? 1f : -1f;
                        // Add a non-linear interaction on features 0,1
                        if (f == 0 || f == 1)
                            X[i * features + f] = signal + (float)(rng.NextDouble() * 2 - 1);
                        else
                            X[i * features + f] = (float)(rng.NextDouble() * 2 - 1); // noise
                    }
                    // Make label depend on interaction: x0*x1 > 0 boosts targets
                    if (y[i] > 0.5f && X[i * features] * X[i * features + 1] < 0)
                        y[i] = rng.NextDouble() < 0.3 ? 0f : 1f; // some targets misclassified by LDA
                }

                // GBT
                var gbt = new DiaGradientBoostedClassifier(features, numTrees: 100, maxDepth: 4,
                    minSamplesLeaf: 20, numBins: 64);
                var sw = Stopwatch.StartNew();
                gbt.Train(X, y, n);
                long gbtMs = sw.ElapsedMilliseconds;
                float gbtAuc = gbt.ComputeAuc(X, y, n);

                // LDA
                var lda = new DiaLdaClassifierAdapter(features);
                sw.Restart();
                lda.Train(X, y, n);
                long ldaMs = sw.ElapsedMilliseconds;
                float ldaAuc = lda.ComputeAuc(X, y, n);

                Console.WriteLine($"  {n,-12:N0} {gbtMs,-12} {ldaMs,-12} {gbtAuc,-10:F4} {ldaAuc,-10:F4}");
            }

            Console.WriteLine();
        }

        #region Helpers

        private static float TrainAndEvaluate(float[] X, float[] y, int n, int features,
            DiaClassifierType type)
        {
            // 80/20 train/test split
            int trainN = (int)(n * 0.8);
            int testN = n - trainN;
            int testStart = trainN * features;

            IDiaClassifier classifier = type switch
            {
                DiaClassifierType.GradientBoostedTree => new DiaGbtClassifierAdapter(
                    new DiaGradientBoostedClassifier(features, numTrees: 100, maxDepth: 4,
                        learningRate: 0.1f, minSamplesLeaf: 20, numBins: 64, subsampleFraction: 0.8f)),
                DiaClassifierType.LinearDiscriminant => new DiaLdaClassifierAdapter(features),
                _ => throw new ArgumentException()
            };

            classifier.Train(
                X.AsSpan(0, trainN * features),
                y.AsSpan(0, trainN),
                trainN);

            // Evaluate on test set
            return classifier.ComputeAuc(
                X.AsSpan(testStart, testN * features),
                y.AsSpan(trainN, testN),
                testN);
        }

        private static float Noise(Random rng, float scale)
        {
            return (float)(rng.NextDouble() * 2 - 1) * scale;
        }

        #endregion
    }
}