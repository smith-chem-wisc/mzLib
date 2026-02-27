// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Classifier approaches for combining DIA features into a single discriminant score.
    /// </summary>
    public enum ClassifierType
    {
        MaxApexTemporal,
        FixedLinearCombination,
        LDA,
        LogisticRegression,
    }

    /// <summary>
    /// A trained linear discriminant classifier for DIA precursor scoring.
    /// Immutable after training. Thread-safe for scoring.
    /// 
    /// Phase 10.5 changes:
    ///   - Updated for 13-feature vector (removed RawCosine + RtWindowHalfWidth, added RtDeviationSquared)
    ///   - Added SweepRegularization() for Phase 10.5c lambda optimization
    ///   - Added interaction feature support (ApexScore × MeanFragCorr) for Phase 10.5d
    ///   - Added WeightChangeL2() for convergence detection (Phase 11 prep)
    ///   - Updated FixedLinear weights for new feature layout
    /// </summary>
    public sealed class DiaLinearDiscriminant
    {
        public ClassifierType Type { get; }
        public float[] Weights { get; }
        public float Bias { get; }
        public float[] FeatureMeans { get; }
        public float[] FeatureStds { get; }

        /// <summary>
        /// Whether the ApexScore × MeanFragCorr interaction feature is active.
        /// When true, Weights has length ClassifierFeatureCount + 1.
        /// </summary>
        public bool UseInteractionFeature { get; }

        /// <summary>Total weight count including any interaction features.</summary>
        public int TotalWeightCount => UseInteractionFeature
            ? DiaFeatureVector.ClassifierFeatureCount + 1
            : DiaFeatureVector.ClassifierFeatureCount;

        private DiaLinearDiscriminant(ClassifierType type, float[] weights, float bias,
            float[] featureMeans = null, float[] featureStds = null,
            bool useInteraction = false)
        {
            Type = type;
            Weights = weights;
            Bias = bias;
            FeatureMeans = featureMeans;
            FeatureStds = featureStds;
            UseInteractionFeature = useInteraction;
        }

        // ════════════════════════════════════════════════════════════════
        //  Scoring
        // ════════════════════════════════════════════════════════════════

        public float Score(in DiaFeatureVector fv)
        {
            if (Type == ClassifierType.MaxApexTemporal)
                return MathF.Max(fv.ApexScore, fv.TemporalScore);

            int baseCount = DiaFeatureVector.ClassifierFeatureCount;
            Span<float> features = stackalloc float[baseCount];
            fv.WriteTo(features);

            if (FeatureMeans != null && FeatureStds != null)
            {
                for (int i = 0; i < baseCount; i++)
                {
                    float v = float.IsNaN(features[i]) ? 0f : features[i];
                    float std = FeatureStds[i];
                    features[i] = std > 1e-8f ? (v - FeatureMeans[i]) / std : 0f;
                }
            }

            float score = Bias;
            for (int i = 0; i < baseCount; i++)
                score += Weights[i] * features[i];

            // Interaction term: ApexScore (0) × MeanFragCorr (3), standardized
            if (UseInteractionFeature && Weights.Length > baseCount)
            {
                float interactionVal = features[0] * features[3]; // already standardized
                score += Weights[baseCount] * interactionVal;
            }

            return score;
        }

        public void ScoreBatch(ReadOnlySpan<DiaFeatureVector> features, Span<float> scores)
        {
            for (int i = 0; i < features.Length; i++)
                scores[i] = Score(in features[i]);
        }

        // ════════════════════════════════════════════════════════════════
        //  Factory methods
        // ════════════════════════════════════════════════════════════════

        public static DiaLinearDiscriminant CreateMaxApexTemporal()
        {
            var weights = new float[DiaFeatureVector.ClassifierFeatureCount];
            return new DiaLinearDiscriminant(ClassifierType.MaxApexTemporal, weights, 0f);
        }

        /// <summary>
        /// Hand-tuned weights — updated for 13-feature layout (Phase 10.5).
        /// Feature order: Apex, Temporal, SpectralAngle,
        ///   MeanFragCorr, MinFragCorr, FragDetRate, LogIntensity, IntCV,
        ///   MedianXicDepth, XicDepthCV, TimePointsUsed, RtDev, RtDevSq
        /// </summary>
        public static DiaLinearDiscriminant CreateFixedLinear()
        {
            var w = new float[DiaFeatureVector.ClassifierFeatureCount];
            w[0] = 0.30f;   // ApexScore
            w[1] = 0.20f;   // TemporalScore
            w[2] = 0.0f;    // SpectralAngle (correlated with temporal)
            w[3] = 0.20f;   // MeanFragCorr — high-leverage feature
            w[4] = 0.05f;   // MinFragCorr — weakest-link detection
            w[5] = 0.08f;   // FragmentDetectionRate
            w[6] = 0.03f;   // LogTotalIntensity
            w[7] = -0.04f;  // IntensityCV (lower = better)
            w[8] = 0.02f;   // MedianXicDepth
            w[9] = -0.02f;  // XicDepthCV (lower = better)
            w[10] = 0.02f;  // TimePointsUsed
            w[11] = -0.05f; // RtDeviationMinutes (lower = better)
            w[12] = -0.02f; // RtDeviationSquared (quadratic penalty)
            return new DiaLinearDiscriminant(ClassifierType.FixedLinearCombination, w, 0f);
        }

        // ════════════════════════════════════════════════════════════════
        //  LDA Training
        // ════════════════════════════════════════════════════════════════

        public static DiaLinearDiscriminant TrainLDA(
            ReadOnlySpan<DiaFeatureVector> positives,
            ReadOnlySpan<DiaFeatureVector> negatives)
        {
            int n = DiaFeatureVector.ClassifierFeatureCount;
            float[] meanPos = new float[n], meanNeg = new float[n];
            float[] globalMean = new float[n], globalStd = new float[n];
            ComputeMeanAndStd(positives, negatives, meanPos, meanNeg, globalMean, globalStd);

            float[][] posF = ExtractAndStandardize(positives, globalMean, globalStd);
            float[][] negF = ExtractAndStandardize(negatives, globalMean, globalStd);

            float[] sMeanPos = ComputeMean(posF, n);
            float[] sMeanNeg = ComputeMean(negF, n);

            float[,] sw = new float[n, n];
            AccumulateCovariance(posF, sMeanPos, sw);
            AccumulateCovariance(negF, sMeanNeg, sw);
            float ridge = 1e-4f * (posF.Length + negF.Length);
            for (int i = 0; i < n; i++) sw[i, i] += ridge;

            float[] meanDiff = new float[n];
            for (int j = 0; j < n; j++) meanDiff[j] = sMeanPos[j] - sMeanNeg[j];
            float[] weights = SolveLinearSystem(sw, meanDiff);

            float wNorm = 0;
            for (int j = 0; j < n; j++) wNorm += weights[j] * weights[j];
            wNorm = MathF.Sqrt(wNorm);
            if (wNorm > 1e-8f)
                for (int j = 0; j < n; j++) weights[j] /= wNorm;

            float midpoint = 0;
            for (int j = 0; j < n; j++)
                midpoint += weights[j] * (sMeanPos[j] + sMeanNeg[j]) / 2f;

            return new DiaLinearDiscriminant(
                ClassifierType.LDA, weights, -midpoint, globalMean, globalStd);
        }

        // ════════════════════════════════════════════════════════════════
        //  Logistic Regression Training
        // ════════════════════════════════════════════════════════════════

        public static DiaLinearDiscriminant TrainLogisticRegression(
            ReadOnlySpan<DiaFeatureVector> positives,
            ReadOnlySpan<DiaFeatureVector> negatives,
            float learningRate = 0.01f,
            float l2Lambda = 0.001f,
            int maxEpochs = 200,
            int batchSize = 256,
            bool useInteraction = false)
        {
            int baseN = DiaFeatureVector.ClassifierFeatureCount;
            int totalN = useInteraction ? baseN + 1 : baseN;

            float[] globalMean = new float[baseN], globalStd = new float[baseN];
            ComputeMeanAndStd(positives, negatives, null, null, globalMean, globalStd);

            float[][] posF = ExtractStandardizeAndExpand(positives, globalMean, globalStd, useInteraction);
            float[][] negF = ExtractStandardizeAndExpand(negatives, globalMean, globalStd, useInteraction);

            int total = posF.Length + negF.Length;
            float[][] allF = new float[total][];
            float[] labels = new float[total];
            Array.Copy(posF, 0, allF, 0, posF.Length);
            Array.Copy(negF, 0, allF, posF.Length, negF.Length);
            for (int i = 0; i < posF.Length; i++) labels[i] = 1f;

            float[] weights = new float[totalN];
            float bias = 0f;
            var rng = new Random(42);
            int[] indices = new int[total];
            for (int i = 0; i < total; i++) indices[i] = i;

            for (int epoch = 0; epoch < maxEpochs; epoch++)
            {
                for (int i = total - 1; i > 0; i--)
                {
                    int j = rng.Next(i + 1);
                    (indices[i], indices[j]) = (indices[j], indices[i]);
                }

                for (int batchStart = 0; batchStart < total; batchStart += batchSize)
                {
                    int batchEnd = Math.Min(batchStart + batchSize, total);
                    int batchLen = batchEnd - batchStart;
                    float[] gradW = new float[totalN];
                    float gradB = 0f;

                    for (int b = batchStart; b < batchEnd; b++)
                    {
                        int idx = indices[b];
                        float[] x = allF[idx];
                        float y = labels[idx];
                        float z = bias;
                        for (int j = 0; j < totalN; j++) z += weights[j] * x[j];
                        float p = Sigmoid(z);
                        float error = p - y;
                        for (int j = 0; j < totalN; j++) gradW[j] += error * x[j];
                        gradB += error;
                    }

                    float invBatch = 1f / batchLen;
                    for (int j = 0; j < totalN; j++)
                    {
                        gradW[j] = gradW[j] * invBatch + l2Lambda * weights[j];
                        weights[j] -= learningRate * gradW[j];
                    }
                    bias -= learningRate * gradB * invBatch;
                }
            }

            // Expand globalMean/Std if interaction feature is used
            float[] expandedMean = globalMean;
            float[] expandedStd = globalStd;
            if (useInteraction)
            {
                expandedMean = new float[totalN];
                expandedStd = new float[totalN];
                Array.Copy(globalMean, expandedMean, baseN);
                Array.Copy(globalStd, expandedStd, baseN);
                // Interaction mean/std approximation: compute from training data
                ComputeInteractionStats(allF, baseN, out float intMean, out float intStd);
                expandedMean[baseN] = intMean;
                expandedStd[baseN] = intStd;
            }

            return new DiaLinearDiscriminant(
                ClassifierType.LogisticRegression, weights, bias,
                expandedMean, expandedStd, useInteraction);
        }

        // ════════════════════════════════════════════════════════════════
        //  Phase 10.5c: Regularization Sweep
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Cross-validates logistic regression across a range of L2 regularization strengths.
        /// Returns (lambda, AUC, TP@1%FP) for each value.
        /// 
        /// Uses k-fold CV on the provided positive/negative feature vectors.
        /// </summary>
        public static List<(float Lambda, double Auc, double TpAt1Fp)> SweepRegularization(
            ReadOnlySpan<DiaFeatureVector> positives,
            ReadOnlySpan<DiaFeatureVector> negatives,
            float[] lambdas = null,
            int kFolds = 3,
            float learningRate = 0.01f,
            int maxEpochs = 200,
            int batchSize = 256,
            bool useInteraction = false)
        {
            lambdas ??= new float[] { 1e-4f, 5e-4f, 1e-3f, 5e-3f, 1e-2f, 5e-2f, 1e-1f };

            // Convert to arrays for fold splitting
            var posArray = new DiaFeatureVector[positives.Length];
            positives.CopyTo(posArray);
            var negArray = new DiaFeatureVector[negatives.Length];
            negatives.CopyTo(negArray);

            var results = new List<(float Lambda, double Auc, double TpAt1Fp)>(lambdas.Length);

            foreach (float lambda in lambdas)
            {
                var (auc, tp1fp) = CrossValidate(posArray, negArray, kFolds,
                    learningRate, lambda, maxEpochs, batchSize, useInteraction);
                results.Add((lambda, auc, tp1fp));
            }

            return results;
        }

        /// <summary>
        /// K-fold cross-validation. Returns (AUC, TP@1%FP) computed from held-out scores.
        /// </summary>
        public static (double Auc, double TpAt1Fp) CrossValidate(
            DiaFeatureVector[] positives,
            DiaFeatureVector[] negatives,
            int kFolds,
            float learningRate,
            float l2Lambda,
            int maxEpochs,
            int batchSize,
            bool useInteraction = false)
        {
            int nPos = positives.Length;
            int nNeg = negatives.Length;
            int posFoldSize = nPos / kFolds;
            int negFoldSize = nNeg / kFolds;

            var allScores = new List<(float Score, bool IsPositive)>();

            for (int fold = 0; fold < kFolds; fold++)
            {
                int posTestStart = fold * posFoldSize;
                int posTestEnd = (fold == kFolds - 1) ? nPos : posTestStart + posFoldSize;
                int negTestStart = fold * negFoldSize;
                int negTestEnd = (fold == kFolds - 1) ? nNeg : negTestStart + negFoldSize;

                // Build train sets (excluding test fold)
                var trainPos = new List<DiaFeatureVector>(nPos);
                var trainNeg = new List<DiaFeatureVector>(nNeg);

                for (int i = 0; i < nPos; i++)
                    if (i < posTestStart || i >= posTestEnd)
                        trainPos.Add(positives[i]);
                for (int i = 0; i < nNeg; i++)
                    if (i < negTestStart || i >= negTestEnd)
                        trainNeg.Add(negatives[i]);

                var model = TrainLogisticRegression(
                    trainPos.ToArray().AsSpan(),
                    trainNeg.ToArray().AsSpan(),
                    learningRate, l2Lambda, maxEpochs, batchSize, useInteraction);

                // Score test fold
                for (int i = posTestStart; i < posTestEnd; i++)
                    allScores.Add((model.Score(in positives[i]), true));
                for (int i = negTestStart; i < negTestEnd; i++)
                    allScores.Add((model.Score(in negatives[i]), false));
            }

            return ComputeRocMetrics(allScores);
        }

        // ════════════════════════════════════════════════════════════════
        //  Diagnostics
        // ════════════════════════════════════════════════════════════════

        public string DescribeWeights()
        {
            var sb = new System.Text.StringBuilder();
            sb.AppendLine("Classifier: " + Type);
            sb.AppendLine("Bias: " + Bias.ToString("F6"));
            if (UseInteractionFeature)
                sb.AppendLine("Interaction: ApexScore × MeanFragCorr (enabled)");
            sb.AppendLine("Feature weights (sorted by |weight|):");

            int baseN = DiaFeatureVector.ClassifierFeatureCount;
            int totalN = Weights.Length;
            var indexed = new (string Name, float Weight, float AbsWeight)[totalN];

            for (int i = 0; i < baseN; i++)
                indexed[i] = (DiaFeatureVector.FeatureNames[i], Weights[i], MathF.Abs(Weights[i]));
            if (totalN > baseN)
                indexed[baseN] = ("ApexScore×MeanFragCorr", Weights[baseN], MathF.Abs(Weights[baseN]));

            Array.Sort(indexed, (a, b) => b.AbsWeight.CompareTo(a.AbsWeight));

            for (int i = 0; i < indexed.Length; i++)
            {
                string sign = indexed[i].Weight >= 0 ? "+" : "";
                sb.AppendLine("  " + sign + indexed[i].Weight.ToString("F4").PadLeft(8) + "  " + indexed[i].Name);
            }
            return sb.ToString();
        }

        /// <summary>
        /// Computes L2 norm of the weight vector (excluding bias).
        /// Used for convergence checking in iterative retraining (Phase 11).
        /// </summary>
        public float WeightNormL2()
        {
            float sum = 0f;
            for (int i = 0; i < Weights.Length; i++)
                sum += Weights[i] * Weights[i];
            return MathF.Sqrt(sum);
        }

        /// <summary>
        /// Fractional L2 change: ||w_new - w_old|| / ||w_old||.
        /// Used for convergence detection in iterative retraining (Phase 11).
        /// </summary>
        public float WeightChangeL2(float[] previousWeights)
        {
            if (previousWeights == null || previousWeights.Length != Weights.Length)
                return float.PositiveInfinity;

            float diffSq = 0f, oldSq = 0f;
            for (int i = 0; i < Weights.Length; i++)
            {
                float diff = Weights[i] - previousWeights[i];
                diffSq += diff * diff;
                oldSq += previousWeights[i] * previousWeights[i];
            }
            return oldSq > 1e-12f ? MathF.Sqrt(diffSq) / MathF.Sqrt(oldSq) : float.PositiveInfinity;
        }

        // ════════════════════════════════════════════════════════════════
        //  ROC / AUC Metrics
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Computes AUC and TP@1%FP from a list of scored results.
        /// </summary>
        public static (double Auc, double TpAt1Fp) ComputeRocMetrics(
            List<(float Score, bool IsPositive)> scoredResults)
        {
            scoredResults.Sort((a, b) => b.Score.CompareTo(a.Score));

            int totalPos = 0, totalNeg = 0;
            for (int i = 0; i < scoredResults.Count; i++)
            {
                if (scoredResults[i].IsPositive) totalPos++;
                else totalNeg++;
            }
            if (totalPos == 0 || totalNeg == 0) return (0.5, 0.0);

            double auc = 0;
            int tp = 0, fp = 0;
            double prevFpr = 0, prevTpr = 0;
            double tpAt1Fp = 0;
            bool found = false;

            for (int i = 0; i < scoredResults.Count; i++)
            {
                if (scoredResults[i].IsPositive) tp++;
                else fp++;

                double tpr = (double)tp / totalPos;
                double fpr = (double)fp / totalNeg;

                auc += (fpr - prevFpr) * (tpr + prevTpr) / 2.0;

                if (!found && fpr > 0.01)
                {
                    tpAt1Fp = prevTpr;
                    found = true;
                }

                prevFpr = fpr;
                prevTpr = tpr;
            }
            auc += (1.0 - prevFpr) * (1.0 + prevTpr) / 2.0;
            if (!found) tpAt1Fp = prevTpr;

            return (auc, tpAt1Fp);
        }

        // ════════════════════════════════════════════════════════════════
        //  Internal Helpers
        // ════════════════════════════════════════════════════════════════

        private static void ComputeMeanAndStd(
            ReadOnlySpan<DiaFeatureVector> positives, ReadOnlySpan<DiaFeatureVector> negatives,
            float[] meanPos, float[] meanNeg, float[] globalMean, float[] globalStd)
        {
            int n = DiaFeatureVector.ClassifierFeatureCount;
            int totalCount = positives.Length + negatives.Length;
            Span<float> buf = stackalloc float[n];
            float[] sumAll = new float[n], sumSqAll = new float[n];
            float[] sumPos = meanPos != null ? new float[n] : null;
            float[] sumNeg = meanNeg != null ? new float[n] : null;

            for (int i = 0; i < positives.Length; i++)
            {
                positives[i].WriteTo(buf);
                for (int j = 0; j < n; j++)
                {
                    float v = float.IsNaN(buf[j]) ? 0f : buf[j];
                    sumAll[j] += v; sumSqAll[j] += v * v;
                    if (sumPos != null) sumPos[j] += v;
                }
            }
            for (int i = 0; i < negatives.Length; i++)
            {
                negatives[i].WriteTo(buf);
                for (int j = 0; j < n; j++)
                {
                    float v = float.IsNaN(buf[j]) ? 0f : buf[j];
                    sumAll[j] += v; sumSqAll[j] += v * v;
                    if (sumNeg != null) sumNeg[j] += v;
                }
            }
            for (int j = 0; j < n; j++)
            {
                globalMean[j] = sumAll[j] / totalCount;
                float variance = sumSqAll[j] / totalCount - globalMean[j] * globalMean[j];
                globalStd[j] = variance > 0 ? MathF.Sqrt(variance) : 1e-8f;
                if (meanPos != null) meanPos[j] = positives.Length > 0 ? sumPos[j] / positives.Length : 0;
                if (meanNeg != null) meanNeg[j] = negatives.Length > 0 ? sumNeg[j] / negatives.Length : 0;
            }
        }

        private static float[][] ExtractAndStandardize(
            ReadOnlySpan<DiaFeatureVector> vectors, float[] mean, float[] std)
        {
            int n = DiaFeatureVector.ClassifierFeatureCount;
            float[][] result = new float[vectors.Length][];
            Span<float> buf = stackalloc float[n];
            for (int i = 0; i < vectors.Length; i++)
            {
                vectors[i].WriteTo(buf);
                result[i] = new float[n];
                for (int j = 0; j < n; j++)
                {
                    float v = float.IsNaN(buf[j]) ? 0f : buf[j];
                    result[i][j] = std[j] > 1e-8f ? (v - mean[j]) / std[j] : 0f;
                }
            }
            return result;
        }

        /// <summary>
        /// Extracts, standardizes, and optionally appends interaction feature.
        /// </summary>
        private static float[][] ExtractStandardizeAndExpand(
            ReadOnlySpan<DiaFeatureVector> vectors, float[] mean, float[] std,
            bool useInteraction)
        {
            int baseN = DiaFeatureVector.ClassifierFeatureCount;
            int totalN = useInteraction ? baseN + 1 : baseN;
            float[][] result = new float[vectors.Length][];
            Span<float> buf = stackalloc float[baseN];

            for (int i = 0; i < vectors.Length; i++)
            {
                vectors[i].WriteTo(buf);
                result[i] = new float[totalN];
                for (int j = 0; j < baseN; j++)
                {
                    float v = float.IsNaN(buf[j]) ? 0f : buf[j];
                    result[i][j] = std[j] > 1e-8f ? (v - mean[j]) / std[j] : 0f;
                }
                if (useInteraction)
                {
                    // Interaction: standardized ApexScore (0) × standardized MeanFragCorr (3)
                    result[i][baseN] = result[i][0] * result[i][3];
                }
            }
            return result;
        }

        private static void ComputeInteractionStats(float[][] allF, int baseN,
            out float mean, out float std)
        {
            double sum = 0, sumSq = 0;
            for (int i = 0; i < allF.Length; i++)
            {
                float v = allF[i].Length > baseN ? allF[i][baseN] : 0f;
                sum += v;
                sumSq += v * v;
            }
            mean = (float)(sum / allF.Length);
            float var = (float)(sumSq / allF.Length - (double)mean * mean);
            std = var > 0 ? MathF.Sqrt(var) : 1e-8f;
        }

        private static float[] ComputeMean(float[][] data, int n)
        {
            float[] mean = new float[n];
            for (int i = 0; i < data.Length; i++)
                for (int j = 0; j < n; j++) mean[j] += data[i][j];
            for (int j = 0; j < n; j++) mean[j] /= data.Length;
            return mean;
        }

        private static void AccumulateCovariance(float[][] features, float[] mean, float[,] cov)
        {
            int n = mean.Length;
            for (int i = 0; i < features.Length; i++)
                for (int j = 0; j < n; j++)
                {
                    float dj = features[i][j] - mean[j];
                    for (int k = j; k < n; k++)
                    {
                        float v = dj * (features[i][k] - mean[k]);
                        cov[j, k] += v;
                        if (j != k) cov[k, j] += v;
                    }
                }
        }

        private static float[] SolveLinearSystem(float[,] A, float[] b)
        {
            int n = b.Length;
            float[,] L = new float[n, n];
            bool ok = true;
            for (int i = 0; i < n && ok; i++)
                for (int j = 0; j <= i; j++)
                {
                    float sum = A[i, j];
                    for (int k = 0; k < j; k++) sum -= L[i, k] * L[j, k];
                    if (i == j) { if (sum <= 0) { ok = false; break; } L[i, j] = MathF.Sqrt(sum); }
                    else L[i, j] = sum / L[j, j];
                }
            float[] x = new float[n];
            if (ok)
            {
                float[] y = new float[n];
                for (int i = 0; i < n; i++) { float s = b[i]; for (int j = 0; j < i; j++) s -= L[i, j] * y[j]; y[i] = s / L[i, i]; }
                for (int i = n - 1; i >= 0; i--) { float s = y[i]; for (int j = i + 1; j < n; j++) s -= L[j, i] * x[j]; x[i] = s / L[i, i]; }
            }
            else
                for (int i = 0; i < n; i++) x[i] = A[i, i] > 1e-10f ? b[i] / A[i, i] : 0f;
            return x;
        }

        private static float Sigmoid(float x)
        {
            if (x > 20f) return 1f;
            if (x < -20f) return 0f;
            return 1f / (1f + MathF.Exp(-x));
        }
    }
}