// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: Development/Dia/Phase10_5_FeatureRefinementBenchmark.cs

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using MassSpectrometry.Dia;

namespace Development.Dia
{
    /// <summary>
    /// Phase 10.5 Feature Refinement Benchmark.
    /// 
    /// Validates changes against the PXD005573 HeLa dataset:
    ///   10.5a — Feature Pruning: RawCosine and RtWindowHalfWidth removed
    ///   10.5b — RtDeviationSquared added as quadratic penalty
    ///   10.5c — Regularization sweep: λ ∈ {1e-4..1e-1} with 3-fold CV
    ///   10.5d — Interaction feature: ApexScore × MeanFragCorr (keep if AUC > 0.001 improvement)
    /// 
    /// Entry point:
    ///   var results = pipeline.GetScoredResults(); // from Phase 10 pipeline
    ///   Phase10_5_FeatureRefinementBenchmark.RunAll(results);
    /// </summary>
    public static class Phase10_5_FeatureRefinementBenchmark
    {
        public static void RunAll(List<DiaSearchResult> results, string outputTsvPath = null)
        {
            Console.WriteLine("╔══════════════════════════════════════════════════════╗");
            Console.WriteLine("║       Phase 10.5 — Feature Refinement Benchmark     ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════╝");
            Console.WriteLine();

            // ── 10.5a: Validate feature pruning ─────────────────────────
            ValidateFeaturePruning();

            // ── Compute feature vectors ─────────────────────────────────
            Console.WriteLine("Computing feature vectors...");
            var sw = Stopwatch.StartNew();
            var vectors = new DiaFeatureVector[results.Count];
            for (int i = 0; i < results.Count; i++)
                vectors[i] = DiaFeatureExtractor.ComputeFeatures(results[i], i);
            sw.Stop();
            Console.WriteLine($"  {vectors.Length:N0} vectors in {sw.ElapsedMilliseconds} ms");
            Console.WriteLine();

            // ── 10.5b: Validate RtDeviationSquared ──────────────────────
            ValidateRtDeviationSquared(vectors);

            // ── Feature separation analysis (updated 13 features) ───────
            AnalyzeFeatureSeparation(vectors);

            // ── Prepare training data (40/40 quantile split) ────────────
            SplitByApexQuantile(vectors,
                out DiaFeatureVector[] positives, out DiaFeatureVector[] negatives);

            Console.WriteLine($"Training set: {positives.Length} positive (top 40%), {negatives.Length} negative (bottom 40%)");
            Console.WriteLine();

            // ── 10.5c: Regularization sweep ─────────────────────────────
            float bestLambda = RunRegularizationSweep(positives, negatives);

            // ── 10.5d: Interaction feature test ─────────────────────────
            TestInteractionFeature(positives, negatives, bestLambda);

            // ── Feature correlation matrix (for diagnostics) ────────────
            PrintFeatureCorrelations(vectors);

            // ── Write updated TSV ───────────────────────────────────────
            if (outputTsvPath != null)
                WriteFeatureTsv(results, vectors, outputTsvPath);

            Console.WriteLine("Phase 10.5 validation complete.");
        }

        // ════════════════════════════════════════════════════════════════
        //  10.5a: Feature Pruning Validation
        // ════════════════════════════════════════════════════════════════

        private static void ValidateFeaturePruning()
        {
            Console.WriteLine("─── 10.5a: Feature Pruning Validation ───");

            bool ok = true;

            // Verify feature count
            Console.Write($"  ClassifierFeatureCount = {DiaFeatureVector.ClassifierFeatureCount}  ");
            if (DiaFeatureVector.ClassifierFeatureCount == 13)
                Console.WriteLine("PASS ✓");
            else { Console.WriteLine("FAIL ✗ (expected 13)"); ok = false; }

            // Verify RawCosine removed
            bool noRawCosine = !Array.Exists(DiaFeatureVector.FeatureNames, n => n == "RawCosine");
            Console.WriteLine($"  RawCosine removed: {(noRawCosine ? "PASS ✓" : "FAIL ✗")}");
            ok &= noRawCosine;

            // Verify RtWindowHalfWidth removed
            bool noRtWindow = !Array.Exists(DiaFeatureVector.FeatureNames, n => n == "RtWindowHalfWidth");
            Console.WriteLine($"  RtWindowHalfWidth removed: {(noRtWindow ? "PASS ✓" : "FAIL ✗")}");
            ok &= noRtWindow;

            // Verify RtDeviationSquared added
            bool hasRtSq = Array.Exists(DiaFeatureVector.FeatureNames, n => n == "RtDeviationSquared");
            Console.WriteLine($"  RtDeviationSquared present: {(hasRtSq ? "PASS ✓" : "FAIL ✗")}");
            ok &= hasRtSq;

            // Verify FeatureNames.Length == ClassifierFeatureCount
            bool namesMatch = DiaFeatureVector.FeatureNames.Length == DiaFeatureVector.ClassifierFeatureCount;
            Console.WriteLine($"  FeatureNames.Length matches: {(namesMatch ? "PASS ✓" : "FAIL ✗")}");
            ok &= namesMatch;

            // WriteTo/ReadFrom roundtrip
            var fv = new DiaFeatureVector
            {
                ApexScore = 0.95f,
                TemporalScore = 0.78f,
                SpectralAngle = 0.55f,
                MeanFragmentCorrelation = 0.65f,
                MinFragmentCorrelation = -0.2f,
                FragmentDetectionRate = 0.92f,
                LogTotalIntensity = 5.5f,
                IntensityCV = 1.1f,
                MedianXicDepth = 25f,
                XicDepthCV = 0.3f,
                TimePointsUsed = 42,
                RtDeviationMinutes = 0.12f,
                RtDeviationSquared = 0.0144f
            };
            Span<float> buf = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];
            fv.WriteTo(buf);
            bool rtSqMatch = MathF.Abs(buf[12] - 0.0144f) < 1e-6f;
            Console.WriteLine($"  WriteTo roundtrip (RtDevSq): {(rtSqMatch ? "PASS ✓" : "FAIL ✗")}");
            ok &= rtSqMatch;

            Console.WriteLine($"\n  Feature names: {string.Join(", ", DiaFeatureVector.FeatureNames)}");
            Console.WriteLine($"\n  Overall: {(ok ? "ALL CHECKS PASSED ✓" : "SOME CHECKS FAILED ✗")}");
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  10.5b: RtDeviationSquared Validation
        // ════════════════════════════════════════════════════════════════

        private static void ValidateRtDeviationSquared(DiaFeatureVector[] vectors)
        {
            Console.WriteLine("─── 10.5b: RtDeviationSquared Validation ───");

            int nonZero = 0, errorCount = 0;
            double sumSq = 0, sumLin = 0;

            for (int i = 0; i < vectors.Length; i++)
            {
                float lin = vectors[i].RtDeviationMinutes;
                float sq = vectors[i].RtDeviationSquared;
                float expected = lin * lin;

                if (MathF.Abs(sq - expected) > 1e-4f) errorCount++;
                if (sq > 1e-8f) nonZero++;
                sumSq += sq;
                sumLin += lin;
            }

            Console.WriteLine($"  Vectors: {vectors.Length:N0}");
            Console.WriteLine($"  Non-zero RtDeviationSquared: {nonZero:N0} ({100.0 * nonZero / vectors.Length:F1}%)");
            Console.WriteLine($"  Mean RtDeviation: {sumLin / vectors.Length:F4} min");
            Console.WriteLine($"  Mean RtDeviationSquared: {sumSq / vectors.Length:F6} min²");
            Console.WriteLine($"  Consistency errors (|sq - lin²| > 1e-4): {errorCount}  {(errorCount == 0 ? "PASS ✓" : "WARN")}");
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  Feature Separation Analysis
        // ════════════════════════════════════════════════════════════════

        private static void AnalyzeFeatureSeparation(DiaFeatureVector[] vectors)
        {
            Console.WriteLine("─── Feature Separation (top 30% vs bottom 30% by ApexScore) ───");

            // Sort by ApexScore
            int[] sortedIdx = Enumerable.Range(0, vectors.Length)
                .OrderBy(i => vectors[i].ApexScore).ToArray();

            int q30 = (int)(vectors.Length * 0.30);
            int q70 = (int)(vectors.Length * 0.70);

            int n = DiaFeatureVector.ClassifierFeatureCount;
            Span<float> buf = stackalloc float[n];

            Console.WriteLine($"\n  {"Feature",-22} {"Separation",10} {"HiQ Median",12} {"LoQ Median",12}");
            Console.WriteLine($"  {new string('-', 56)}");

            for (int f = 0; f < n; f++)
            {
                var hiq = new List<float>();
                var loq = new List<float>();

                for (int r = 0; r < q30; r++)
                {
                    vectors[sortedIdx[r]].WriteTo(buf);
                    if (!float.IsNaN(buf[f])) loq.Add(buf[f]);
                }
                for (int r = q70; r < vectors.Length; r++)
                {
                    vectors[sortedIdx[r]].WriteTo(buf);
                    if (!float.IsNaN(buf[f])) hiq.Add(buf[f]);
                }

                float hiqMed = Median(hiq);
                float loqMed = Median(loq);
                float hiqMad = MAD(hiq, hiqMed);
                float loqMad = MAD(loq, loqMed);
                float denom = hiqMad + loqMad;
                float sep = denom > 1e-8f ? MathF.Abs(hiqMed - loqMed) / denom : 0f;

                Console.WriteLine($"  {DiaFeatureVector.FeatureNames[f],-22} {sep,10:F2} {hiqMed,12:F3} {loqMed,12:F3}");
            }
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  10.5c: Regularization Sweep
        // ════════════════════════════════════════════════════════════════

        private static float RunRegularizationSweep(
            DiaFeatureVector[] positives, DiaFeatureVector[] negatives)
        {
            Console.WriteLine("─── 10.5c: Regularization Sweep (3-fold CV, LogReg) ───");

            var results = DiaLinearDiscriminant.SweepRegularization(
                positives.AsSpan(), negatives.AsSpan(),
                lambdas: null, kFolds: 3, learningRate: 0.01f,
                maxEpochs: 300, batchSize: 256, useInteraction: false);

            Console.WriteLine($"\n  {"Lambda",10} {"AUC",8} {"TP@1%FP",8}");
            Console.WriteLine($"  {new string('-', 26)}");

            double bestAuc = 0;
            float bestLambda = 1e-3f;

            foreach (var (lambda, auc, tp1fp) in results)
            {
                string marker = "";
                if (auc > bestAuc) { bestAuc = auc; bestLambda = lambda; marker = " ◄"; }
                Console.WriteLine($"  {lambda,10:E1} {auc,8:F4} {tp1fp,8:F4}{marker}");
            }

            Console.WriteLine($"\n  Best λ = {bestLambda:E1} (AUC = {bestAuc:F4})");

            // Train final model with best lambda
            var model = DiaLinearDiscriminant.TrainLogisticRegression(
                positives.AsSpan(), negatives.AsSpan(),
                learningRate: 0.01f, l2Lambda: bestLambda, maxEpochs: 300);

            Console.WriteLine($"\n{model.DescribeWeights()}");
            Console.WriteLine();

            return bestLambda;
        }

        // ════════════════════════════════════════════════════════════════
        //  10.5d: Interaction Feature Test
        // ════════════════════════════════════════════════════════════════

        private static void TestInteractionFeature(
            DiaFeatureVector[] positives, DiaFeatureVector[] negatives, float lambda)
        {
            Console.WriteLine("─── 10.5d: Interaction Feature Test (ApexScore × MeanFragCorr) ───");

            var (aucBase, tp1fpBase) = DiaLinearDiscriminant.CrossValidate(
                positives, negatives, kFolds: 3,
                learningRate: 0.01f, l2Lambda: lambda, maxEpochs: 300, batchSize: 256,
                useInteraction: false);

            var (aucInter, tp1fpInter) = DiaLinearDiscriminant.CrossValidate(
                positives, negatives, kFolds: 3,
                learningRate: 0.01f, l2Lambda: lambda, maxEpochs: 300, batchSize: 256,
                useInteraction: true);

            Console.WriteLine($"\n  Without interaction: AUC={aucBase:F4}  TP@1%FP={tp1fpBase:F4}");
            Console.WriteLine($"  With interaction:    AUC={aucInter:F4}  TP@1%FP={tp1fpInter:F4}");

            double diff = aucInter - aucBase;
            bool keep = diff > 0.001;
            Console.WriteLine($"  ΔAUC = {diff:+0.0000;-0.0000}");
            Console.WriteLine($"  Decision: {(keep ? "KEEP interaction (ΔAUC > 0.001)" : "DROP interaction (ΔAUC ≤ 0.001)")}");

            if (keep)
            {
                var model = DiaLinearDiscriminant.TrainLogisticRegression(
                    positives.AsSpan(), negatives.AsSpan(),
                    learningRate: 0.01f, l2Lambda: lambda, maxEpochs: 300,
                    useInteraction: true);
                Console.WriteLine($"\n{model.DescribeWeights()}");
            }
            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  Feature Correlations
        // ════════════════════════════════════════════════════════════════

        private static void PrintFeatureCorrelations(DiaFeatureVector[] vectors)
        {
            Console.WriteLine("─── Feature Correlation Matrix (top pairs) ───");

            int n = DiaFeatureVector.ClassifierFeatureCount;
            int m = vectors.Length;

            // Extract all features
            float[][] feats = new float[m][];
            Span<float> buf = stackalloc float[n];
            for (int i = 0; i < m; i++)
            {
                vectors[i].WriteTo(buf);
                feats[i] = new float[n];
                buf.CopyTo(feats[i]);
            }

            // Compute means
            float[] means = new float[n];
            for (int j = 0; j < n; j++)
            {
                double sum = 0;
                for (int i = 0; i < m; i++)
                {
                    float v = float.IsNaN(feats[i][j]) ? 0f : feats[i][j];
                    sum += v;
                }
                means[j] = (float)(sum / m);
            }

            // Compute pairwise Pearson correlations
            var pairs = new List<(string A, string B, float R)>();
            for (int a = 0; a < n; a++)
            {
                for (int b = a + 1; b < n; b++)
                {
                    double sumAB = 0, sumA2 = 0, sumB2 = 0;
                    for (int i = 0; i < m; i++)
                    {
                        float va = (float.IsNaN(feats[i][a]) ? 0f : feats[i][a]) - means[a];
                        float vb = (float.IsNaN(feats[i][b]) ? 0f : feats[i][b]) - means[b];
                        sumAB += va * vb;
                        sumA2 += va * va;
                        sumB2 += vb * vb;
                    }
                    float denom = (float)(Math.Sqrt(sumA2) * Math.Sqrt(sumB2));
                    float r = denom > 1e-8f ? (float)(sumAB / denom) : 0f;
                    if (MathF.Abs(r) > 0.3f)
                        pairs.Add((DiaFeatureVector.FeatureNames[a], DiaFeatureVector.FeatureNames[b], r));
                }
            }

            pairs.Sort((x, y) => MathF.Abs(y.R).CompareTo(MathF.Abs(x.R)));
            Console.WriteLine($"\n  Pairs with |r| > 0.3:");
            foreach (var (a, b, r) in pairs.Take(15))
                Console.WriteLine($"  {a,-22} ↔ {b,-22}  r = {r:F3}");

            Console.WriteLine();
        }

        // ════════════════════════════════════════════════════════════════
        //  Training Data Split
        // ════════════════════════════════════════════════════════════════

        private static void SplitByApexQuantile(
            DiaFeatureVector[] vectors,
            out DiaFeatureVector[] positives,
            out DiaFeatureVector[] negatives)
        {
            int[] sortedIdx = Enumerable.Range(0, vectors.Length)
                .OrderBy(i => vectors[i].ApexScore).ToArray();

            int q40 = (int)(vectors.Length * 0.40);
            int q60 = (int)(vectors.Length * 0.60);

            negatives = new DiaFeatureVector[q40];
            for (int i = 0; i < q40; i++)
                negatives[i] = vectors[sortedIdx[i]];

            int posCount = vectors.Length - q60;
            positives = new DiaFeatureVector[posCount];
            for (int i = 0; i < posCount; i++)
                positives[i] = vectors[sortedIdx[q60 + i]];
        }

        // ════════════════════════════════════════════════════════════════
        //  TSV Output
        // ════════════════════════════════════════════════════════════════

        private static void WriteFeatureTsv(
            List<DiaSearchResult> results, DiaFeatureVector[] vectors, string path)
        {
            Console.WriteLine($"Writing feature TSV: {path}");
            using var writer = new StreamWriter(path);

            // Header
            writer.Write("Sequence\tCharge\tPrecursorMz\tWindowId\tIsDecoy");
            for (int f = 0; f < DiaFeatureVector.ClassifierFeatureCount; f++)
            {
                writer.Write('\t');
                writer.Write(DiaFeatureVector.FeatureNames[f]);
            }
            writer.WriteLine();

            // Data
            int n = DiaFeatureVector.ClassifierFeatureCount;
            Span<float> buf = stackalloc float[n];
            for (int i = 0; i < results.Count; i++)
            {
                var r = results[i];
                writer.Write($"{r.Sequence}\t{r.ChargeState}\t{r.PrecursorMz:F4}\t{r.WindowId}\t{r.IsDecoy}");
                vectors[i].WriteTo(buf);
                for (int f = 0; f < n; f++)
                {
                    writer.Write('\t');
                    writer.Write(buf[f].ToString("G6"));
                }
                writer.WriteLine();
            }
            Console.WriteLine($"  {results.Count:N0} rows × {n + 5} columns");
        }

        // ════════════════════════════════════════════════════════════════
        //  Statistical Helpers
        // ════════════════════════════════════════════════════════════════

        private static float Median(List<float> values)
        {
            if (values.Count == 0) return 0f;
            values.Sort();
            int n = values.Count;
            return n % 2 == 1 ? values[n / 2] : (values[n / 2 - 1] + values[n / 2]) / 2f;
        }

        private static float MAD(List<float> values, float median)
        {
            if (values.Count == 0) return 0f;
            var devs = new List<float>(values.Count);
            for (int i = 0; i < values.Count; i++)
                devs.Add(MathF.Abs(values[i] - median));
            return Median(devs);
        }
    }
}