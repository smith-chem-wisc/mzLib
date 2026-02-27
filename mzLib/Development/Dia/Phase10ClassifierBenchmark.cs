// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry;
using MassSpectrometry.Dia;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace Development.Dia
{
    /// <summary>
    /// Phase 10: Multi-Feature Classifier Benchmark (Revised)
    /// 
    /// Now uses Koina/Prosit-predicted fragment library (.msp) instead of
    /// synthetic fragments. This dramatically improves fragment detection rate,
    /// enabling temporal scoring, coelution features, and RT deviation to work.
    /// </summary>
    public static class Phase10ClassifierBenchmark
    {
        /// <summary>
        /// Run the full Phase 10 benchmark.
        /// </summary>
        /// <param name="rawFilePath">Path to the Thermo .raw (or .mzML) DIA file</param>
        /// <param name="mspLibraryPath">Path to Koina .msp predicted library</param>
        /// <param name="groundTruthTsvPath">Path to DIA-NN ground truth TSV (for RT lookup)</param>
        /// <param name="outputDir">Output directory for TSV export</param>
        public static void RunAll(
            string rawFilePath,
            string mspLibraryPath,
            string groundTruthTsvPath,
            string outputDir = null)
        {
            outputDir ??= Path.GetDirectoryName(rawFilePath);
            Console.WriteLine("================================================================");
            Console.WriteLine("  Phase 10: Multi-Feature Classifier Benchmark (Koina Library)");
            Console.WriteLine("================================================================");
            Console.WriteLine();

            // -- Step 1: Load Koina library with RT from ground truth ----------
            Console.WriteLine("--- Loading Koina .msp library ---------------------------------");
            var sw = Stopwatch.StartNew();

            Dictionary<string, double> rtLookup = null;
            if (!string.IsNullOrEmpty(groundTruthTsvPath) && File.Exists(groundTruthTsvPath))
            {
                rtLookup = KoinaMspParser.BuildRtLookupFromDiannTsv(groundTruthTsvPath);
            }

            var precursors = KoinaMspParser.Parse(mspLibraryPath, rtLookup, minIntensity: 0.05f);
            Console.WriteLine($"  Library load: {sw.ElapsedMilliseconds}ms");
            Console.WriteLine($"  Precursors: {precursors.Count:N0}");

            // Quick stats
            int hasRt = precursors.Count(p => p.RetentionTime.HasValue);
            float avgFrags = precursors.Count > 0 ? (float)precursors.Sum(p => p.FragmentCount) / precursors.Count : 0;
            Console.WriteLine($"  With RT: {hasRt:N0} ({100.0 * hasRt / precursors.Count:F1}%)");
            Console.WriteLine($"  Avg fragments/precursor: {avgFrags:F1}");
            Console.WriteLine();

            // -- Step 2: Load raw file and build index -------------------------
            Console.WriteLine("--- Loading raw file and building index ------------------------");
            sw.Restart();

            MsDataFile msDataFile;
            string ext = Path.GetExtension(rawFilePath).ToLowerInvariant();
            if (ext == ".raw")
                msDataFile = new ThermoRawFileReader(rawFilePath);
            else
                msDataFile = new Mzml(rawFilePath);
            msDataFile.LoadAllStaticData();
            var loadTime = sw.Elapsed;

            sw.Restart();
            var scans = msDataFile.GetAllScansList().ToArray();
            using var index = DiaScanIndexBuilder.Build(scans);
            var buildTime = sw.Elapsed;
            Console.WriteLine($"  Load: {loadTime.TotalSeconds:F1}s | Build: {buildTime.TotalMilliseconds:F0}ms");
            Console.WriteLine($"  Scans: {index.ScanCount:N0} | Windows: {index.WindowCount} | Peaks: {index.TotalPeakCount:N0}");
            Console.WriteLine();

            // -- Step 3: Run extraction + temporal scoring ----------------------
            Console.WriteLine("--- Running extraction pipeline --------------------------------");
            var parameters = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 1.14f,
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0f,
                MaxThreads = -1,
                ScoringStrategy = ScoringStrategy.TemporalCosine,
            };

            sw.Restart();
            var genResult = DiaLibraryQueryGenerator.Generate(precursors, index, parameters);
            Console.WriteLine($"  Query generation: {sw.ElapsedMilliseconds}ms | {genResult.Queries.Length:N0} queries");
            Console.WriteLine($"  Skipped (no window): {genResult.SkippedNoWindow}");

            sw.Restart();
            using var orchestrator = new DiaExtractionOrchestrator(index);
            var extractionResult = orchestrator.ExtractAll(genResult.Queries,
                maxDegreeOfParallelism: parameters.EffectiveMaxThreads);
            Console.WriteLine($"  Extraction: {sw.ElapsedMilliseconds}ms | {extractionResult.TotalDataPoints:N0} data points");

            sw.Restart();
            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursors, genResult, extractionResult, parameters);
            Console.WriteLine($"  Temporal scoring: {sw.ElapsedMilliseconds}ms | {results.Count:N0} results");
            Console.WriteLine();

            // Quick diagnostic: fragment detection rate with real library
            if (results.Count > 0)
            {
                float medDetRate = results.Select(r => r.FragmentDetectionRate).OrderBy(x => x).ElementAt(results.Count / 2);
                float medFragDet = results.Select(r => (float)r.FragmentsDetected).OrderBy(x => x).ElementAt(results.Count / 2);
                int hasTemp = results.Count(r => r.TemporalCosineScore > 0 && !float.IsNaN(r.TemporalCosineScore));
                int hasCorr = results.Count(r => !float.IsNaN(r.MeanFragmentCorrelation));
                int hasApexRt = results.Count(r => !float.IsNaN(r.ObservedApexRt));

                Console.WriteLine("--- Library Quality Diagnostics --------------------------------");
                Console.WriteLine($"  Median fragment detection rate: {medDetRate:F3} ({medFragDet:F0} fragments)");
                Console.WriteLine($"  Temporal score populated: {hasTemp:N0} ({100.0 * hasTemp / results.Count:F1}%)");
                Console.WriteLine($"  Fragment correlation populated: {hasCorr:N0} ({100.0 * hasCorr / results.Count:F1}%)");
                Console.WriteLine($"  Observed apex RT populated: {hasApexRt:N0} ({100.0 * hasApexRt / results.Count:F1}%)");
                Console.WriteLine();
            }

            // -- Step 4: Compute feature vectors -------------------------------
            Console.WriteLine("--- Computing feature vectors ----------------------------------");
            sw.Restart();
            var features = new DiaFeatureVector[results.Count];
            for (int i = 0; i < results.Count; i++)
                features[i] = DiaFeatureExtractor.ComputeFeatures(results[i], i);
            Console.WriteLine($"  Feature computation: {sw.ElapsedMilliseconds}ms | {features.Length:N0} vectors");
            Console.WriteLine();

            // -- Step 5: Analyze feature distributions -------------------------
            AnalyzeFeatureDistributions(features);

            // -- Step 6: Compare classifiers -----------------------------------
            CompareClassifiers(features);

            // -- Step 7: Export TSV --------------------------------------------
            string tsvPath = Path.Combine(outputDir, "phase10_features_koina.tsv");
            ExportFeatureTsv(features, results, tsvPath);
            Console.WriteLine($"  Feature TSV exported: {tsvPath}");
            Console.WriteLine();
        }

        // =================================================================
        //  Feature distribution analysis
        // =================================================================

        private static void AnalyzeFeatureDistributions(DiaFeatureVector[] features)
        {
            Console.WriteLine("--- Feature Distributions --------------------------------------");
            Console.WriteLine();

            int n = features.Length;
            int nF = DiaFeatureVector.ClassifierFeatureCount;

            float[][] allValues = new float[nF][];
            for (int j = 0; j < nF; j++) allValues[j] = new float[n];

            Span<float> buf = stackalloc float[nF];
            for (int i = 0; i < n; i++)
            {
                features[i].WriteTo(buf);
                for (int j = 0; j < nF; j++) allValues[j][i] = buf[j];
            }

            // Quantile-based split: top 30% by ApexScore vs bottom 30%
            var apexSorted = features.Select(f => f.ApexScore).OrderBy(x => x).ToArray();
            float hiThresh = apexSorted[(int)(n * 0.7)];
            float loThresh = apexSorted[(int)(n * 0.3)];

            var highQ = features.Where(f => f.ApexScore >= hiThresh).ToArray();
            var lowQ = features.Where(f => f.ApexScore <= loThresh).ToArray();

            Console.WriteLine($"  Total: {n:N0} | High quality (top 30%): {highQ.Length:N0} | Low quality (bottom 30%): {lowQ.Length:N0}");
            Console.WriteLine($"  Split thresholds: apex >= {hiThresh:F3} (high), apex <= {loThresh:F3} (low)");
            Console.WriteLine();

            Console.WriteLine("  " + "Feature".PadRight(22) +
                "Median".PadLeft(8) + "Mean".PadLeft(8) + "Std".PadLeft(8) +
                " | " + "HiQ Med".PadLeft(8) + "LoQ Med".PadLeft(8) + "Sep".PadLeft(6));
            Console.WriteLine("  " + new string('-', 22) +
                new string('-', 8) + new string('-', 8) + new string('-', 8) +
                " | " + new string('-', 8) + new string('-', 8) + new string('-', 6));

            for (int j = 0; j < nF; j++)
            {
                var sorted = allValues[j].ToArray();
                Array.Sort(sorted);
                float median = sorted[n / 2];
                float mean = sorted.Average();
                float std = MathF.Sqrt(sorted.Average(v => (v - mean) * (v - mean)));

                float hiMed = MedianForGroup(highQ, j);
                float loMed = MedianForGroup(lowQ, j);
                float sep = std > 1e-6f ? MathF.Abs(hiMed - loMed) / std : 0f;

                Console.WriteLine("  " + DiaFeatureVector.FeatureNames[j].PadRight(22) +
                    median.ToString("F4").PadLeft(8) + mean.ToString("F4").PadLeft(8) +
                    std.ToString("F4").PadLeft(8) + " | " +
                    hiMed.ToString("F4").PadLeft(8) + loMed.ToString("F4").PadLeft(8) +
                    sep.ToString("F2").PadLeft(6));
            }
            Console.WriteLine();

            // Apex/Temporal complementarity
            int bothHigh = features.Count(f => f.ApexScore > 0.7f && f.TemporalScore > 0.7f);
            int apexOnly = features.Count(f => f.ApexScore > 0.7f && f.TemporalScore <= 0.7f);
            int tempOnly = features.Count(f => f.ApexScore <= 0.7f && f.TemporalScore > 0.7f);
            int neither = features.Count(f => f.ApexScore <= 0.7f && f.TemporalScore <= 0.7f);
            Console.WriteLine("  Apex/Temporal complementarity at 0.7 threshold:");
            Console.WriteLine($"    Both >0.7:  {bothHigh,6} ({100.0 * bothHigh / n:F1}%)");
            Console.WriteLine($"    Apex only:  {apexOnly,6} ({100.0 * apexOnly / n:F1}%)");
            Console.WriteLine($"    Temp only:  {tempOnly,6} ({100.0 * tempOnly / n:F1}%)");
            Console.WriteLine($"    Neither:    {neither,6} ({100.0 * neither / n:F1}%)");
            Console.WriteLine();
        }

        private static float MedianForGroup(DiaFeatureVector[] group, int featureIndex)
        {
            if (group.Length == 0) return float.NaN;
            float[] vals = new float[group.Length];
            Span<float> buf = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];
            for (int i = 0; i < group.Length; i++)
            {
                group[i].WriteTo(buf);
                vals[i] = buf[featureIndex];
            }
            Array.Sort(vals);
            return vals[vals.Length / 2];
        }

        // =================================================================
        //  Classifier comparison (3-fold CV)
        // =================================================================

        private static void CompareClassifiers(DiaFeatureVector[] features)
        {
            Console.WriteLine("--- Classifier Comparison (3-fold CV) --------------------------");
            Console.WriteLine();

            // Quantile-based split: top 40% positive, bottom 40% negative
            var sorted = features.OrderBy(f => MathF.Max(f.ApexScore, f.TemporalScore)).ToArray();
            int n40 = (int)(sorted.Length * 0.4);
            var negatives = sorted.Take(n40).ToArray();
            var positives = sorted.Skip(sorted.Length - n40).ToArray();

            Console.WriteLine($"  Positives (top 40%): {positives.Length:N0}");
            Console.WriteLine($"  Negatives (bottom 40%): {negatives.Length:N0}");
            Console.WriteLine($"  Excluded (middle 20%): {sorted.Length - positives.Length - negatives.Length:N0}");
            Console.WriteLine();

            if (positives.Length < 100 || negatives.Length < 100)
            {
                Console.WriteLine("  Insufficient data for cross-validation.");
                return;
            }

            int nFolds = 3;
            var rng = new Random(42);
            var posShuf = positives.OrderBy(_ => rng.Next()).ToArray();
            var negShuf = negatives.OrderBy(_ => rng.Next()).ToArray();

            var types = new[]
            {
                ClassifierType.MaxApexTemporal,
                ClassifierType.FixedLinearCombination,
                ClassifierType.LDA,
                ClassifierType.LogisticRegression,
            };

            var allScores = new Dictionary<ClassifierType, List<(float score, bool isPos)>>();
            foreach (var ct in types) allScores[ct] = new List<(float, bool)>();

            for (int fold = 0; fold < nFolds; fold++)
            {
                var (posTrain, posTest) = SplitFold(posShuf, fold, nFolds);
                var (negTrain, negTest) = SplitFold(negShuf, fold, nFolds);

                foreach (var ct in types)
                {
                    DiaLinearDiscriminant classifier;
                    switch (ct)
                    {
                        case ClassifierType.MaxApexTemporal:
                            classifier = DiaLinearDiscriminant.CreateMaxApexTemporal(); break;
                        case ClassifierType.FixedLinearCombination:
                            classifier = DiaLinearDiscriminant.CreateFixedLinear(); break;
                        case ClassifierType.LDA:
                            classifier = DiaLinearDiscriminant.TrainLDA(posTrain, negTrain); break;
                        case ClassifierType.LogisticRegression:
                            classifier = DiaLinearDiscriminant.TrainLogisticRegression(
                                posTrain, negTrain, learningRate: 0.05f, l2Lambda: 0.001f, maxEpochs: 300); break;
                        default: continue;
                    }

                    foreach (var fv in posTest) allScores[ct].Add((classifier.Score(in fv), true));
                    foreach (var fv in negTest) allScores[ct].Add((classifier.Score(in fv), false));
                }
            }

            Console.WriteLine("  " + "Classifier".PadRight(28) +
                "AUC".PadLeft(7) + "Acc".PadLeft(7) + "TP@5%FP".PadLeft(8) + "TP@1%FP".PadLeft(8));
            Console.WriteLine("  " + new string('-', 28) +
                new string('-', 7) + new string('-', 7) + new string('-', 8) + new string('-', 8));

            foreach (var ct in types)
            {
                var m = ComputeMetrics(allScores[ct]);
                Console.WriteLine("  " + ct.ToString().PadRight(28) +
                    m.Auc.ToString("F4").PadLeft(7) + m.BestAcc.ToString("F4").PadLeft(7) +
                    m.TpAt5Fp.ToString("F4").PadLeft(8) + m.TpAt1Fp.ToString("F4").PadLeft(8));
            }
            Console.WriteLine();

            // Final classifiers
            Console.WriteLine("--- Final Classifiers (trained on all data) --------------------");
            Console.WriteLine();

            var ldaFinal = DiaLinearDiscriminant.TrainLDA(positives, negatives);
            Console.WriteLine(ldaFinal.DescribeWeights());

            var lrFinal = DiaLinearDiscriminant.TrainLogisticRegression(
                positives, negatives, learningRate: 0.05f, l2Lambda: 0.001f, maxEpochs: 300);
            Console.WriteLine(lrFinal.DescribeWeights());

            // Score distribution
            Console.WriteLine("--- Score Distributions (all precursors) -----------------------");
            Console.WriteLine();

            var allClassifiers = new (string Name, DiaLinearDiscriminant C)[]
            {
                ("MaxApexTemporal", DiaLinearDiscriminant.CreateMaxApexTemporal()),
                ("FixedLinear", DiaLinearDiscriminant.CreateFixedLinear()),
                ("LDA", ldaFinal),
                ("LogisticRegression", lrFinal),
            };

            Console.WriteLine("  " + "Classifier".PadRight(24) +
                "Median".PadLeft(8) + ">0.5".PadLeft(8) + ">0.7".PadLeft(8) + ">0.8".PadLeft(8));
            Console.WriteLine("  " + new string('-', 24) +
                new string('-', 8) + new string('-', 8) + new string('-', 8) + new string('-', 8));

            foreach (var (name, c) in allClassifiers)
            {
                float[] scores = new float[features.Length];
                for (int i = 0; i < features.Length; i++) scores[i] = c.Score(in features[i]);
                Array.Sort(scores);

                float median = scores[scores.Length / 2];
                Console.WriteLine("  " + name.PadRight(24) +
                    median.ToString("F4").PadLeft(8) +
                    Pct(CountAbove(scores, 0.5f), features.Length).PadLeft(8) +
                    Pct(CountAbove(scores, 0.7f), features.Length).PadLeft(8) +
                    Pct(CountAbove(scores, 0.8f), features.Length).PadLeft(8));
            }
            Console.WriteLine();
        }

        // =================================================================
        //  Metrics
        // =================================================================

        private struct Metrics { public float Auc, BestAcc, TpAt5Fp, TpAt1Fp; }

        private static Metrics ComputeMetrics(List<(float score, bool isPos)> data)
        {
            data.Sort((a, b) => b.score.CompareTo(a.score));
            int totalPos = data.Count(s => s.isPos);
            int totalNeg = data.Count - totalPos;
            float auc = 0, prevFpr = 0, prevTpr = 0, tpAt5 = 0, tpAt1 = 0, bestAcc = 0;
            int tp = 0, fp = 0;

            for (int i = 0; i < data.Count; i++)
            {
                if (data[i].isPos) tp++; else fp++;
                float tpr = (float)tp / totalPos;
                float fpr = totalNeg > 0 ? (float)fp / totalNeg : 0;
                auc += (fpr - prevFpr) * (tpr + prevTpr) / 2f;
                prevFpr = fpr; prevTpr = tpr;
                float acc = (float)(tp + totalNeg - fp) / data.Count;
                if (acc > bestAcc) bestAcc = acc;
                if (fpr <= 0.05f) tpAt5 = tpr;
                if (fpr <= 0.01f) tpAt1 = tpr;
            }
            return new Metrics { Auc = auc, BestAcc = bestAcc, TpAt5Fp = tpAt5, TpAt1Fp = tpAt1 };
        }

        // =================================================================
        //  TSV export
        // =================================================================

        private static void ExportFeatureTsv(
            DiaFeatureVector[] features, List<DiaSearchResult> results, string path)
        {
            using var w = new StreamWriter(path);
            w.Write("Sequence\tCharge\tPrecursorMz\tIsDecoy\tFragDet\tFragQueried");
            w.Write("\tApexTimeIndex\tTimePointsUsed\tScoringStrategy");
            w.Write("\tObservedApexRt\tMeanFragCorr\tMinFragCorr");
            for (int j = 0; j < DiaFeatureVector.ClassifierFeatureCount; j++)
                w.Write("\t" + DiaFeatureVector.FeatureNames[j]);
            w.WriteLine();

            Span<float> buf = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];
            for (int i = 0; i < features.Length; i++)
            {
                var r = results[i];
                w.Write(r.Sequence + "\t" + r.ChargeState + "\t" +
                    r.PrecursorMz.ToString("F4") + "\t" + r.IsDecoy);
                w.Write("\t" + r.FragmentsDetected + "\t" + r.FragmentsQueried);
                w.Write("\t" + r.ApexTimeIndex + "\t" + r.TimePointsUsed + "\t" + r.ScoringStrategyUsed);
                w.Write("\t" + r.ObservedApexRt.ToString("F4") +
                    "\t" + r.MeanFragmentCorrelation.ToString("F4") +
                    "\t" + r.MinFragmentCorrelation.ToString("F4"));

                features[i].WriteTo(buf);
                for (int j = 0; j < DiaFeatureVector.ClassifierFeatureCount; j++)
                    w.Write("\t" + buf[j].ToString("G6"));
                w.WriteLine();
            }
        }

        // =================================================================
        //  Utility
        // =================================================================

        private static (DiaFeatureVector[] train, DiaFeatureVector[] test) SplitFold(
            DiaFeatureVector[] data, int fold, int nFolds)
        {
            int foldSize = data.Length / nFolds;
            int testStart = fold * foldSize;
            int testEnd = fold == nFolds - 1 ? data.Length : testStart + foldSize;
            var train = new List<DiaFeatureVector>(data.Length);
            var test = new List<DiaFeatureVector>(foldSize + 1);
            for (int i = 0; i < data.Length; i++)
            {
                if (i >= testStart && i < testEnd) test.Add(data[i]);
                else train.Add(data[i]);
            }
            return (train.ToArray(), test.ToArray());
        }

        private static int CountAbove(float[] sorted, float threshold)
        {
            int idx = Array.BinarySearch(sorted, threshold);
            if (idx < 0) idx = ~idx;
            return sorted.Length - idx;
        }

        private static string Pct(int count, int total) =>
            (100.0 * count / total).ToString("F1") + "%";
    }
}