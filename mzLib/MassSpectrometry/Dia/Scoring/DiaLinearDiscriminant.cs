// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Enumerates the classifier approaches available for combining DIA features
    /// into a single discriminant score.
    /// </summary>
    public enum ClassifierType
    {
        /// <summary>max(ApexScore, TemporalScore) — Phase 9 baseline, no training required</summary>
        MaxApexTemporal,

        /// <summary>Weighted linear combination with hand-tuned weights</summary>
        FixedLinearCombination,

        /// <summary>Linear Discriminant Analysis: projects onto the axis maximizing class separation</summary>
        LDA,

        /// <summary>Logistic regression: maximum-likelihood linear boundary via gradient descent</summary>
        LogisticRegression,
    }

    /// <summary>
    /// A trained linear discriminant classifier for DIA precursor scoring.
    /// 
    /// Takes a DiaFeatureVector and produces a single composite score that
    /// combines all features optimally (or near-optimally, depending on approach).
    /// The score can then be used for target-decoy FDR estimation.
    /// 
    /// Immutable after training. Thread-safe for scoring.
    /// </summary>
    public sealed class DiaLinearDiscriminant
    {
        /// <summary>Which approach was used to train this classifier.</summary>
        public ClassifierType Type { get; }

        /// <summary>Feature weights (length = DiaFeatureVector.ClassifierFeatureCount)</summary>
        public float[] Weights { get; }

        /// <summary>Intercept/bias term.</summary>
        public float Bias { get; }

        /// <summary>Feature-wise means used for centering (null if not centering)</summary>
        public float[] FeatureMeans { get; }

        /// <summary>Feature-wise standard deviations used for scaling (null if not scaling)</summary>
        public float[] FeatureStds { get; }

        private DiaLinearDiscriminant(ClassifierType type, float[] weights, float bias,
            float[] featureMeans = null, float[] featureStds = null)
        {
            Type = type;
            Weights = weights;
            Bias = bias;
            FeatureMeans = featureMeans;
            FeatureStds = featureStds;
        }

        /// <summary>
        /// Scores a single feature vector. Higher = more likely true positive.
        /// For MaxApexTemporal, returns max(apex, temporal) instead of a linear combination.
        /// </summary>
        public float Score(in DiaFeatureVector fv)
        {
            if (Type == ClassifierType.MaxApexTemporal)
                return MathF.Max(fv.ApexScore, fv.TemporalScore);

            Span<float> features = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];
            fv.WriteTo(features);

            // Apply standardization if trained with it
            if (FeatureMeans != null && FeatureStds != null)
            {
                for (int i = 0; i < features.Length; i++)
                {
                    float std = FeatureStds[i];
                    features[i] = std > 1e-8f
                        ? (features[i] - FeatureMeans[i]) / std
                        : 0f;
                }
            }

            float score = Bias;
            for (int i = 0; i < features.Length; i++)
                score += Weights[i] * features[i];

            return score;
        }

        /// <summary>
        /// Scores a batch of feature vectors.
        /// </summary>
        public void ScoreBatch(ReadOnlySpan<DiaFeatureVector> features, Span<float> scores)
        {
            for (int i = 0; i < features.Length; i++)
                scores[i] = Score(in features[i]);
        }

        // ════════════════════════════════════════════════════════════════════
        //  Factory methods — each returns an immutable classifier instance
        // ════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Baseline: max(ApexScore, TemporalScore). No training, no weights.
        /// Reference to measure how much the classifier improves things.
        /// </summary>
        public static DiaLinearDiscriminant CreateMaxApexTemporal()
        {
            var weights = new float[DiaFeatureVector.ClassifierFeatureCount];
            return new DiaLinearDiscriminant(ClassifierType.MaxApexTemporal, weights, 0f);
        }

        /// <summary>
        /// Hand-tuned weights based on feature analysis.
        /// "Domain expert" baseline that doesn't require training data.
        /// </summary>
        public static DiaLinearDiscriminant CreateFixedLinear()
        {
            var weights = new float[DiaFeatureVector.ClassifierFeatureCount];
            // Feature order: Apex, Temporal, RawCosine, SpectralAngle,
            //   FragDetRate, LogIntensity, IntCV, MedianXicDepth, XicDepthCV,
            //   TimePointsUsed, RtDeviation, RtWindowHalfWidth
            weights[0] = 0.40f;  // ApexScore — best single predictor from Phase 9
            weights[1] = 0.30f;  // TemporalScore — complementary to apex
            weights[2] = 0.0f;   // RawCosine — correlated with apex/temporal
            weights[3] = 0.0f;   // SpectralAngle — correlated with RawCosine
            weights[4] = 0.10f;  // FragmentDetectionRate — higher is better
            weights[5] = 0.05f;  // LogTotalIntensity — proxy for abundance
            weights[6] = -0.05f; // IntensityCV — lower is better
            weights[7] = 0.02f;  // MedianXicDepth — more points = better
            weights[8] = -0.03f; // XicDepthCV — lower is better
            weights[9] = 0.02f;  // TimePointsUsed — more temporal evidence
            weights[10] = -0.05f; // RtDeviationMinutes — lower is better
            weights[11] = 0.0f;  // RtWindowHalfWidth — not discriminative

            return new DiaLinearDiscriminant(ClassifierType.FixedLinearCombination, weights, 0f);
        }

        /// <summary>
        /// Linear Discriminant Analysis (LDA).
        /// 
        /// Finds the projection w = Sw^{-1} (μ+ − μ−) maximizing between-class / within-class variance.
        /// Features are standardized before training for numerical stability.
        /// </summary>
        public static DiaLinearDiscriminant TrainLDA(
            ReadOnlySpan<DiaFeatureVector> positives,
            ReadOnlySpan<DiaFeatureVector> negatives)
        {
            int n = DiaFeatureVector.ClassifierFeatureCount;

            // Step 1: Compute standardization stats and class means
            float[] meanPos = new float[n], meanNeg = new float[n];
            float[] globalMean = new float[n], globalStd = new float[n];
            ComputeMeanAndStd(positives, negatives, meanPos, meanNeg, globalMean, globalStd);

            // Step 2: Standardize and recompute class means in standardized space
            float[][] posF = ExtractAndStandardize(positives, globalMean, globalStd);
            float[][] negF = ExtractAndStandardize(negatives, globalMean, globalStd);

            float[] sMeanPos = ComputeMean(posF, n);
            float[] sMeanNeg = ComputeMean(negF, n);

            // Step 3: Pooled within-class covariance + ridge regularization
            float[,] sw = new float[n, n];
            AccumulateCovariance(posF, sMeanPos, sw);
            AccumulateCovariance(negF, sMeanNeg, sw);
            float ridge = 1e-4f * (posF.Length + negF.Length);
            for (int i = 0; i < n; i++) sw[i, i] += ridge;

            // Step 4: Solve w = Sw^{-1} * (μ+ − μ-)
            float[] meanDiff = new float[n];
            for (int j = 0; j < n; j++) meanDiff[j] = sMeanPos[j] - sMeanNeg[j];
            float[] weights = SolveLinearSystem(sw, meanDiff);

            // Normalize weights
            float wNorm = 0;
            for (int j = 0; j < n; j++) wNorm += weights[j] * weights[j];
            wNorm = MathF.Sqrt(wNorm);
            if (wNorm > 1e-8f)
                for (int j = 0; j < n; j++) weights[j] /= wNorm;

            // Bias: midpoint between class means maps to 0
            float midpoint = 0;
            for (int j = 0; j < n; j++)
                midpoint += weights[j] * (sMeanPos[j] + sMeanNeg[j]) / 2f;

            return new DiaLinearDiscriminant(
                ClassifierType.LDA, weights, -midpoint, globalMean, globalStd);
        }

        /// <summary>
        /// Logistic regression via mini-batch gradient descent.
        /// Learns P(target | features) = sigmoid(w · x + b).
        /// Minimizes cross-entropy loss with L2 regularization.
        /// </summary>
        public static DiaLinearDiscriminant TrainLogisticRegression(
            ReadOnlySpan<DiaFeatureVector> positives,
            ReadOnlySpan<DiaFeatureVector> negatives,
            float learningRate = 0.01f,
            float l2Lambda = 0.001f,
            int maxEpochs = 200,
            int batchSize = 256)
        {
            int n = DiaFeatureVector.ClassifierFeatureCount;

            float[] globalMean = new float[n], globalStd = new float[n];
            ComputeMeanAndStd(positives, negatives, null, null, globalMean, globalStd);

            float[][] posF = ExtractAndStandardize(positives, globalMean, globalStd);
            float[][] negF = ExtractAndStandardize(negatives, globalMean, globalStd);

            // Combine with labels
            int total = posF.Length + negF.Length;
            float[][] allF = new float[total][];
            float[] labels = new float[total];
            Array.Copy(posF, 0, allF, 0, posF.Length);
            Array.Copy(negF, 0, allF, posF.Length, negF.Length);
            for (int i = 0; i < posF.Length; i++) labels[i] = 1f;

            float[] weights = new float[n];
            float bias = 0f;

            var rng = new Random(42);
            int[] indices = new int[total];
            for (int i = 0; i < total; i++) indices[i] = i;

            for (int epoch = 0; epoch < maxEpochs; epoch++)
            {
                // Shuffle
                for (int i = total - 1; i > 0; i--)
                {
                    int j = rng.Next(i + 1);
                    (indices[i], indices[j]) = (indices[j], indices[i]);
                }

                for (int batchStart = 0; batchStart < total; batchStart += batchSize)
                {
                    int batchEnd = Math.Min(batchStart + batchSize, total);
                    int batchLen = batchEnd - batchStart;

                    float[] gradW = new float[n];
                    float gradB = 0f;

                    for (int b = batchStart; b < batchEnd; b++)
                    {
                        int idx = indices[b];
                        float[] x = allF[idx];
                        float y = labels[idx];

                        float z = bias;
                        for (int j = 0; j < n; j++) z += weights[j] * x[j];
                        float p = Sigmoid(z);

                        float error = p - y;
                        for (int j = 0; j < n; j++)
                            gradW[j] += error * x[j];
                        gradB += error;
                    }

                    float invBatch = 1f / batchLen;
                    for (int j = 0; j < n; j++)
                    {
                        gradW[j] = gradW[j] * invBatch + l2Lambda * weights[j];
                        weights[j] -= learningRate * gradW[j];
                    }
                    bias -= learningRate * gradB * invBatch;
                }
            }

            return new DiaLinearDiscriminant(
                ClassifierType.LogisticRegression, weights, bias, globalMean, globalStd);
        }

        // ════════════════════════════════════════════════════════════════════
        //  Reporting
        // ════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Returns a formatted string describing the classifier weights and feature importance.
        /// </summary>
        public string DescribeWeights()
        {
            var sb = new System.Text.StringBuilder();
            sb.AppendLine($"Classifier: {Type}");
            sb.AppendLine($"Bias: {Bias:F6}");
            sb.AppendLine("Feature weights (sorted by |weight|):");

            var indexed = new (string Name, float Weight, float AbsWeight)[Weights.Length];
            for (int i = 0; i < Weights.Length; i++)
                indexed[i] = (DiaFeatureVector.FeatureNames[i], Weights[i], MathF.Abs(Weights[i]));

            Array.Sort(indexed, (a, b) => b.AbsWeight.CompareTo(a.AbsWeight));

            for (int i = 0; i < indexed.Length; i++)
            {
                string sign = indexed[i].Weight >= 0 ? "+" : "";
                sb.AppendLine($"  {sign}{indexed[i].Weight,8:F4}  {indexed[i].Name}");
            }
            return sb.ToString();
        }

        // ════════════════════════════════════════════════════════════════════
        //  Helper methods
        // ════════════════════════════════════════════════════════════════════

        private static void ComputeMeanAndStd(
            ReadOnlySpan<DiaFeatureVector> positives,
            ReadOnlySpan<DiaFeatureVector> negatives,
            float[] meanPos, float[] meanNeg,
            float[] globalMean, float[] globalStd)
        {
            int n = DiaFeatureVector.ClassifierFeatureCount;
            int totalCount = positives.Length + negatives.Length;

            Span<float> buf = stackalloc float[n];
            float[] sumAll = new float[n];
            float[] sumSqAll = new float[n];
            float[] sumPos = meanPos != null ? new float[n] : null;
            float[] sumNeg = meanNeg != null ? new float[n] : null;

            for (int i = 0; i < positives.Length; i++)
            {
                positives[i].WriteTo(buf);
                for (int j = 0; j < n; j++)
                {
                    float v = float.IsNaN(buf[j]) ? 0f : buf[j];
                    sumAll[j] += v;
                    sumSqAll[j] += v * v;
                    if (sumPos != null) sumPos[j] += v;
                }
            }
            for (int i = 0; i < negatives.Length; i++)
            {
                negatives[i].WriteTo(buf);
                for (int j = 0; j < n; j++)
                {
                    float v = float.IsNaN(buf[j]) ? 0f : buf[j];
                    sumAll[j] += v;
                    sumSqAll[j] += v * v;
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

        private static float[] ComputeMean(float[][] data, int n)
        {
            float[] mean = new float[n];
            for (int i = 0; i < data.Length; i++)
                for (int j = 0; j < n; j++)
                    mean[j] += data[i][j];
            for (int j = 0; j < n; j++) mean[j] /= data.Length;
            return mean;
        }

        private static void AccumulateCovariance(float[][] features, float[] mean, float[,] cov)
        {
            int n = mean.Length;
            for (int i = 0; i < features.Length; i++)
            {
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
        }

        /// <summary>
        /// Solves Ax = b via Cholesky decomposition. Falls back to diagonal if Cholesky fails.
        /// </summary>
        private static float[] SolveLinearSystem(float[,] A, float[] b)
        {
            int n = b.Length;
            float[,] L = new float[n, n];
            bool ok = true;

            for (int i = 0; i < n && ok; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    float sum = A[i, j];
                    for (int k = 0; k < j; k++) sum -= L[i, k] * L[j, k];

                    if (i == j)
                    {
                        if (sum <= 0) { ok = false; break; }
                        L[i, j] = MathF.Sqrt(sum);
                    }
                    else
                    {
                        L[i, j] = sum / L[j, j];
                    }
                }
            }

            float[] x = new float[n];
            if (ok)
            {
                float[] y = new float[n];
                for (int i = 0; i < n; i++)
                {
                    float sum = b[i];
                    for (int j = 0; j < i; j++) sum -= L[i, j] * y[j];
                    y[i] = sum / L[i, i];
                }
                for (int i = n - 1; i >= 0; i--)
                {
                    float sum = y[i];
                    for (int j = i + 1; j < n; j++) sum -= L[j, i] * x[j];
                    x[i] = sum / L[i, i];
                }
            }
            else
            {
                for (int i = 0; i < n; i++)
                    x[i] = A[i, i] > 1e-10f ? b[i] / A[i, i] : 0f;
            }
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