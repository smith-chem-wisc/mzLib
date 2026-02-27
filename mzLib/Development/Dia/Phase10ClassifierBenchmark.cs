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
    /// Phase 10: Multi-Feature Classifier Benchmark
    /// 
    /// Runs on the real HeLa DIA dataset (PXD005573) with DIA-NN ground truth.
    /// 
    /// This benchmark:
    ///   1. Loads the full pipeline (raw file -> index -> extraction -> temporal scoring)
    ///   2. Computes feature vectors for all ~39K precursors
    ///   3. Analyzes feature distributions (high-quality vs low-quality)
    ///   4. Trains multiple classifier approaches in cross-validation
    ///   5. Compares discrimination power head-to-head
    ///   6. Exports TSV for external analysis
    /// </summary>
    public static class Phase10ClassifierBenchmark
    {
        public static void RunAll(
            string rawFilePath,
            string groundTruthTsvPath,
            string outputDir = null)
        {
            outputDir ??= Path.GetDirectoryName(rawFilePath);
            Console.WriteLine("================================================================");
            Console.WriteLine("  Phase 10: Multi-Feature Classifier Benchmark");
            Console.WriteLine("================================================================");
            Console.WriteLine();

            // -- Step 1: Load ground truth -----------------------------------------
            Console.WriteLine("--- Loading DIA-NN ground truth --------------------------------");
            var groundTruth = LoadDiannGroundTruth(groundTruthTsvPath);
            Console.WriteLine($"  Precursors loaded: {groundTruth.Count:N0}");
            Console.WriteLine();

            // -- Step 2: Load raw file and build index -----------------------------
            Console.WriteLine("--- Loading raw file and building index ------------------------");
            var sw = Stopwatch.StartNew();

            MsDataFile msDataFile;
            string ext = Path.GetExtension(rawFilePath).ToLowerInvariant();
            if (ext == ".raw")
            {
                msDataFile = new ThermoRawFileReader(rawFilePath);
            }
            else
            {
                msDataFile = new Mzml(rawFilePath);
            }
            msDataFile.LoadAllStaticData();
            var loadTime = sw.Elapsed;

            sw.Restart();
            var scans = msDataFile.GetAllScansList().ToArray();
            using var index = DiaScanIndexBuilder.Build(scans);
            var buildTime = sw.Elapsed;
            Console.WriteLine($"  Load: {loadTime.TotalSeconds:F1}s | Build: {buildTime.TotalMilliseconds:F0}ms");
            Console.WriteLine($"  Scans: {index.ScanCount:N0} | Windows: {index.WindowCount} | Peaks: {index.TotalPeakCount:N0}");
            Console.WriteLine();

            // -- Step 3: Build library from ground truth ---------------------------
            Console.WriteLine("--- Building library inputs ------------------------------------");
            var precursors = BuildPrecursorInputs(groundTruth, index);
            Console.WriteLine($"  Library precursors: {precursors.Count:N0}");
            Console.WriteLine();

            // -- Step 4: Run extraction + temporal scoring pipeline -----------------
            Console.WriteLine("--- Running extraction pipeline (k=1.0 calibrated windows) -----");
            var parameters = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 1.14f, // k=1.0 sigma from Phase 9 optimum
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0f,
                MaxThreads = -1,
                ScoringStrategy = ScoringStrategy.TemporalCosine,
            };

            sw.Restart();
            var genResult = DiaLibraryQueryGenerator.Generate(precursors, index, parameters);
            Console.WriteLine($"  Query generation: {sw.ElapsedMilliseconds}ms | {genResult.Queries.Length:N0} queries");

            sw.Restart();
            using var orchestrator = new DiaExtractionOrchestrator(index);
            var extractionResult = orchestrator.ExtractAll(genResult.Queries,
                maxDegreeOfParallelism: parameters.EffectiveMaxThreads);
            Console.WriteLine($"  Extraction: {sw.ElapsedMilliseconds}ms | {extractionResult.TotalDataPoints:N0} data points");

            // Use AssembleResultsWithTemporalScoring - the Phase 9 path that computes
            // ApexDotProductScore, TemporalCosineScore, and the primary DotProductScore
            sw.Restart();
            var results = DiaLibraryQueryGenerator.AssembleResultsWithTemporalScoring(
                precursors, genResult, extractionResult, parameters);
            Console.WriteLine($"  Temporal scoring: {sw.ElapsedMilliseconds}ms | {results.Count:N0} results");
            Console.WriteLine();

            // -- Step 5: Compute feature vectors -----------------------------------
            Console.WriteLine("--- Computing feature vectors ----------------------------------");
            sw.Restart();
            var features = ComputeAllFeatures(results);
            Console.WriteLine($"  Feature computation: {sw.ElapsedMilliseconds}ms | {features.Length:N0} vectors");
            Console.WriteLine();

            // -- Step 6: Analyze feature distributions -----------------------------
            AnalyzeFeatureDistributions(features);

            // -- Step 7: Compare classifiers via cross-validation ------------------
            CompareClassifiers(features);

            // -- Step 8: Export feature TSV -----------------------------------------
            string tsvPath = Path.Combine(outputDir, "phase10_features.tsv");
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

            // Collect per-feature values
            float[][] allValues = new float[nF][];
            for (int j = 0; j < nF; j++) allValues[j] = new float[n];

            Span<float> buf = stackalloc float[nF];
            for (int i = 0; i < n; i++)
            {
                features[i].WriteTo(buf);
                for (int j = 0; j < nF; j++) allValues[j][i] = buf[j];
            }

            // Split by score quality: max(apex, temporal) > 0.7 vs < 0.3
            var highQ = features.Where(f => MathF.Max(f.ApexScore, f.TemporalScore) > 0.7f).ToArray();
            var lowQ = features.Where(f => MathF.Max(f.ApexScore, f.TemporalScore) < 0.3f).ToArray();

            Console.WriteLine($"  Total: {n:N0} | High quality (>0.7): {highQ.Length:N0} | Low quality (<0.3): {lowQ.Length:N0}");
            Console.WriteLine();

            // Header - use padding with PadLeft/PadRight to avoid interpolation alignment issues
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
            Console.WriteLine("  Sep = |median_high - median_low| / std -- higher = more discriminative");
            Console.WriteLine();

            // Apex/temporal complementarity
            int bothHigh = features.Count(f => f.ApexScore > 0.7f && f.TemporalScore > 0.7f);
            int apexOnly = features.Count(f => f.ApexScore > 0.7f && f.TemporalScore <= 0.7f);
            int tempOnly = features.Count(f => f.ApexScore <= 0.7f && f.TemporalScore > 0.7f);
            int neither = features.Count(f => f.ApexScore <= 0.7f && f.TemporalScore <= 0.7f);
            int union = bothHigh + apexOnly + tempOnly;
            Console.WriteLine("  Apex/Temporal complementarity at 0.7 threshold:");
            Console.WriteLine($"    Both >0.7: {bothHigh,6} ({100.0 * bothHigh / n:F1}%)");
            Console.WriteLine($"    Apex only: {apexOnly,6} ({100.0 * apexOnly / n:F1}%)");
            Console.WriteLine($"    Temp only: {tempOnly,6} ({100.0 * tempOnly / n:F1}%)");
            Console.WriteLine($"    Neither:   {neither,6} ({100.0 * neither / n:F1}%)");
            Console.WriteLine($"    Union:     {union,6} ({100.0 * union / n:F1}%)");
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
        //  Classifier comparison (3-fold cross-validation)
        // =================================================================

        private static void CompareClassifiers(DiaFeatureVector[] features)
        {
            Console.WriteLine("--- Classifier Comparison (3-fold CV) --------------------------");
            Console.WriteLine();

            // Semi-supervised split
            var positives = features.Where(f => MathF.Max(f.ApexScore, f.TemporalScore) > 0.6f).ToArray();
            var negatives = features.Where(f => MathF.Max(f.ApexScore, f.TemporalScore) < 0.3f).ToArray();

            Console.WriteLine($"  Positives (max score > 0.6): {positives.Length:N0}");
            Console.WriteLine($"  Negatives (max score < 0.3): {negatives.Length:N0}");
            Console.WriteLine($"  Excluded (0.3-0.6):          {features.Length - positives.Length - negatives.Length:N0}");
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
                            classifier = DiaLinearDiscriminant.CreateMaxApexTemporal();
                            break;
                        case ClassifierType.FixedLinearCombination:
                            classifier = DiaLinearDiscriminant.CreateFixedLinear();
                            break;
                        case ClassifierType.LDA:
                            classifier = DiaLinearDiscriminant.TrainLDA(posTrain, negTrain);
                            break;
                        case ClassifierType.LogisticRegression:
                            classifier = DiaLinearDiscriminant.TrainLogisticRegression(
                                posTrain, negTrain, learningRate: 0.05f, l2Lambda: 0.001f, maxEpochs: 300);
                            break;
                        default:
                            continue;
                    }

                    foreach (var fv in posTest)
                        allScores[ct].Add((classifier.Score(in fv), true));
                    foreach (var fv in negTest)
                        allScores[ct].Add((classifier.Score(in fv), false));
                }
            }

            // Report header
            Console.WriteLine("  " + "Classifier".PadRight(28) +
                "AUC".PadLeft(7) + "Acc".PadLeft(7) +
                "TP@5%FP".PadLeft(8) + "TP@1%FP".PadLeft(8));
            Console.WriteLine("  " + new string('-', 28) +
                new string('-', 7) + new string('-', 7) +
                new string('-', 8) + new string('-', 8));

            foreach (var ct in types)
            {
                var m = ComputeMetrics(allScores[ct]);
                Console.WriteLine("  " + ct.ToString().PadRight(28) +
                    m.Auc.ToString("F4").PadLeft(7) +
                    m.BestAcc.ToString("F4").PadLeft(7) +
                    m.TpAt5Fp.ToString("F4").PadLeft(8) +
                    m.TpAt1Fp.ToString("F4").PadLeft(8));
            }
            Console.WriteLine();

            // -- Train final classifiers on all data and report weights -----------
            Console.WriteLine("--- Final Classifiers (trained on all data) --------------------");
            Console.WriteLine();

            var ldaFinal = DiaLinearDiscriminant.TrainLDA(positives, negatives);
            Console.WriteLine(ldaFinal.DescribeWeights());

            var lrFinal = DiaLinearDiscriminant.TrainLogisticRegression(
                positives, negatives, learningRate: 0.05f, l2Lambda: 0.001f, maxEpochs: 300);
            Console.WriteLine(lrFinal.DescribeWeights());

            // -- Score distribution on ALL precursors -----------------------------
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
                "Median".PadLeft(8) + ">0.5".PadLeft(8) +
                ">0.7".PadLeft(8) + ">0.8".PadLeft(8));
            Console.WriteLine("  " + new string('-', 24) +
                new string('-', 8) + new string('-', 8) +
                new string('-', 8) + new string('-', 8));

            foreach (var (name, c) in allClassifiers)
            {
                float[] scores = new float[features.Length];
                for (int i = 0; i < features.Length; i++)
                    scores[i] = c.Score(in features[i]);
                Array.Sort(scores);

                float median = scores[scores.Length / 2];
                int a50 = CountAbove(scores, 0.5f);
                int a70 = CountAbove(scores, 0.7f);
                int a80 = CountAbove(scores, 0.8f);

                Console.WriteLine("  " + name.PadRight(24) +
                    median.ToString("F4").PadLeft(8) +
                    Pct(a50, features.Length).PadLeft(8) +
                    Pct(a70, features.Length).PadLeft(8) +
                    Pct(a80, features.Length).PadLeft(8));
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

            float auc = 0, prevFpr = 0, prevTpr = 0;
            float tpAt5 = 0, tpAt1 = 0, bestAcc = 0;
            int tp = 0, fp = 0;

            for (int i = 0; i < data.Count; i++)
            {
                if (data[i].isPos) tp++; else fp++;
                float tpr = (float)tp / totalPos;
                float fpr = totalNeg > 0 ? (float)fp / totalNeg : 0;

                auc += (fpr - prevFpr) * (tpr + prevTpr) / 2f;
                prevFpr = fpr;
                prevTpr = tpr;

                float acc = (float)(tp + totalNeg - fp) / data.Count;
                if (acc > bestAcc) bestAcc = acc;
                if (fpr <= 0.05f) tpAt5 = tpr;
                if (fpr <= 0.01f) tpAt1 = tpr;
            }

            return new Metrics { Auc = auc, BestAcc = bestAcc, TpAt5Fp = tpAt5, TpAt1Fp = tpAt1 };
        }

        // =================================================================
        //  Data loading
        // =================================================================

        private struct DiannEntry
        {
            public string Sequence;
            public int Charge;
            public double PrecursorMz;
            public double RetentionTime;
            public float[] FragMzs;
            public float[] FragIntensities;
        }

        private static Dictionary<string, DiannEntry> LoadDiannGroundTruth(string tsvPath)
        {
            var result = new Dictionary<string, DiannEntry>();
            var lines = File.ReadAllLines(tsvPath);
            if (lines.Length < 2) return result;

            var header = lines[0].Split('\t');
            var cols = new Dictionary<string, int>();
            for (int i = 0; i < header.Length; i++) cols[header[i].Trim()] = i;

            int seqCol = FindCol(cols, "Modified.Sequence", "Stripped.Sequence", "Sequence");
            int chargeCol = FindCol(cols, "Precursor.Charge", "Charge");
            int mzCol = FindCol(cols, "Precursor.Mz", "PrecursorMz");
            int rtCol = FindCol(cols, "RT", "RetentionTime", "iRT");
            int fragMzCol = FindColOpt(cols, "Fragment.Mz.Predicted", "Fragment.Info");
            int fragIntCol = FindColOpt(cols, "Fragment.Quant.Raw", "Fragment.Quant.Corrected");

            for (int row = 1; row < lines.Length; row++)
            {
                var f = lines[row].Split('\t');
                if (f.Length <= Math.Max(seqCol, Math.Max(chargeCol, Math.Max(mzCol, rtCol))))
                    continue;

                string seq = f[seqCol].Trim();
                if (!int.TryParse(f[chargeCol], out int charge)) continue;
                if (!double.TryParse(f[mzCol], out double mz)) continue;
                if (!double.TryParse(f[rtCol], out double rt)) continue;

                string key = seq + "/" + charge;
                if (result.ContainsKey(key)) continue;

                float[] fragMzs = ParseFloats(fragMzCol >= 0 && fragMzCol < f.Length ? f[fragMzCol] : null);
                float[] fragInt = ParseFloats(fragIntCol >= 0 && fragIntCol < f.Length ? f[fragIntCol] : null);

                if (fragMzs == null || fragMzs.Length == 0)
                {
                    fragMzs = SyntheticFragments(mz, charge);
                    fragInt = new float[fragMzs.Length];
                    for (int i = 0; i < fragInt.Length; i++)
                        fragInt[i] = 1000f * (fragMzs.Length - i);
                }
                else if (fragInt == null || fragInt.Length != fragMzs.Length)
                {
                    fragInt = new float[fragMzs.Length];
                    for (int i = 0; i < fragInt.Length; i++) fragInt[i] = 1000f;
                }

                result[key] = new DiannEntry
                {
                    Sequence = seq,
                    Charge = charge,
                    PrecursorMz = mz,
                    RetentionTime = rt,
                    FragMzs = fragMzs,
                    FragIntensities = fragInt
                };
            }
            return result;
        }

        private static List<LibraryPrecursorInput> BuildPrecursorInputs(
            Dictionary<string, DiannEntry> gt, DiaScanIndex index)
        {
            var inputs = new List<LibraryPrecursorInput>(gt.Count);
            int skipped = 0;
            foreach (var e in gt.Values)
            {
                if (index.FindWindowForPrecursorMz(e.PrecursorMz) < 0) { skipped++; continue; }
                inputs.Add(new LibraryPrecursorInput(
                    e.Sequence, e.PrecursorMz, e.Charge, e.RetentionTime,
                    isDecoy: false, e.FragMzs, e.FragIntensities));
            }
            if (skipped > 0)
                Console.WriteLine($"  Skipped (outside windows): {skipped}");
            return inputs;
        }

        // =================================================================
        //  Feature computation
        // =================================================================

        private static DiaFeatureVector[] ComputeAllFeatures(List<DiaSearchResult> results)
        {
            var features = new DiaFeatureVector[results.Count];
            for (int i = 0; i < results.Count; i++)
            {
                features[i] = DiaFeatureExtractor.ComputeFeatures(results[i], i);
            }
            return features;
        }

        // =================================================================
        //  TSV export
        // =================================================================

        private static void ExportFeatureTsv(
            DiaFeatureVector[] features, List<DiaSearchResult> results, string path)
        {
            using var w = new StreamWriter(path);

            // Header
            w.Write("Sequence\tCharge\tPrecursorMz\tIsDecoy\tFragDet\tFragQueried");
            w.Write("\tApexTimeIndex\tTimePointsUsed\tScoringStrategy");
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
                w.Write("\t" + r.ApexTimeIndex + "\t" + r.TimePointsUsed +
                    "\t" + r.ScoringStrategyUsed);

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
                if (i >= testStart && i < testEnd)
                    test.Add(data[i]);
                else
                    train.Add(data[i]);
            }
            return (train.ToArray(), test.ToArray());
        }

        private static int CountAbove(float[] sorted, float threshold)
        {
            int idx = Array.BinarySearch(sorted, threshold);
            if (idx < 0) idx = ~idx;
            return sorted.Length - idx;
        }

        private static string Pct(int count, int total)
        {
            return (100.0 * count / total).ToString("F1") + "%";
        }

        private static int FindCol(Dictionary<string, int> cols, params string[] names)
        {
            foreach (var name in names)
                if (cols.TryGetValue(name, out int idx)) return idx;
            throw new KeyNotFoundException(
                "Column not found. Tried: " + string.Join(", ", names) +
                ". Available: " + string.Join(", ", cols.Keys));
        }

        private static int FindColOpt(Dictionary<string, int> cols, params string[] names)
        {
            foreach (var name in names)
                if (cols.TryGetValue(name, out int idx)) return idx;
            return -1;
        }

        private static float[] ParseFloats(string s)
        {
            if (string.IsNullOrWhiteSpace(s)) return null;
            var parts = s.Split(';', StringSplitOptions.RemoveEmptyEntries);
            var vals = new List<float>(parts.Length);
            foreach (var p in parts)
                if (float.TryParse(p.Trim(), out float v)) vals.Add(v);
            return vals.Count > 0 ? vals.ToArray() : null;
        }

        private static float[] SyntheticFragments(double precursorMz, int charge)
        {
            double mass = (precursorMz - 1.00728) * charge;
            int nFrags = Math.Clamp((int)(mass / 100), 4, 12);
            float[] frags = new float[nFrags];
            double step = (mass * 0.8) / nFrags;
            for (int i = 0; i < nFrags; i++)
                frags[i] = (float)(147.0 + i * step);
            return frags;
        }
    }
}