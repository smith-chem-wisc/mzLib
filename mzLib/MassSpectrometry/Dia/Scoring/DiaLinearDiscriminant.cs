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
    /// </summary>
    public sealed class DiaLinearDiscriminant
    {
        public ClassifierType Type { get; }
        public float[] Weights { get; }
        public float Bias { get; }
        public float[] FeatureMeans { get; }
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

        public float Score(in DiaFeatureVector fv)
        {
            if (Type == ClassifierType.MaxApexTemporal)
                return MathF.Max(fv.ApexScore, fv.TemporalScore);

            Span<float> features = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];
            fv.WriteTo(features);

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
        /// Hand-tuned weights — updated for 14 features including coelution.
        /// Feature order: Apex, Temporal, RawCosine, SpectralAngle,
        ///   MeanFragCorr, MinFragCorr, FragDetRate, LogIntensity, IntCV,
        ///   MedianXicDepth, XicDepthCV, TimePointsUsed, RtDev, RtWindowHW
        /// </summary>
        public static DiaLinearDiscriminant CreateFixedLinear()
        {
            var w = new float[DiaFeatureVector.ClassifierFeatureCount];
            w[0] = 0.30f;   // ApexScore
            w[1] = 0.20f;   // TemporalScore
            w[2] = 0.0f;    // RawCosine (correlated with apex/temporal)
            w[3] = 0.0f;    // SpectralAngle (correlated)
            w[4] = 0.20f;   // MeanFragCorr — high-leverage new feature
            w[5] = 0.05f;   // MinFragCorr — weakest-link detection
            w[6] = 0.08f;   // FragmentDetectionRate
            w[7] = 0.03f;   // LogTotalIntensity
            w[8] = -0.04f;  // IntensityCV (lower = better)
            w[9] = 0.02f;   // MedianXicDepth
            w[10] = -0.02f; // XicDepthCV (lower = better)
            w[11] = 0.02f;  // TimePointsUsed
            w[12] = -0.05f; // RtDeviationMinutes (lower = better)
            w[13] = 0.0f;   // RtWindowHalfWidth
            return new DiaLinearDiscriminant(ClassifierType.FixedLinearCombination, w, 0f);
        }

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
                        for (int j = 0; j < n; j++) gradW[j] += error * x[j];
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

        public string DescribeWeights()
        {
            var sb = new System.Text.StringBuilder();
            sb.AppendLine("Classifier: " + Type);
            sb.AppendLine("Bias: " + Bias.ToString("F6"));
            sb.AppendLine("Feature weights (sorted by |weight|):");

            var indexed = new (string Name, float Weight, float AbsWeight)[Weights.Length];
            for (int i = 0; i < Weights.Length; i++)
                indexed[i] = (DiaFeatureVector.FeatureNames[i], Weights[i], MathF.Abs(Weights[i]));
            Array.Sort(indexed, (a, b) => b.AbsWeight.CompareTo(a.AbsWeight));

            for (int i = 0; i < indexed.Length; i++)
            {
                string sign = indexed[i].Weight >= 0 ? "+" : "";
                sb.AppendLine("  " + sign + indexed[i].Weight.ToString("F4").PadLeft(8) + "  " + indexed[i].Name);
            }
            return sb.ToString();
        }

        // ════════════════════════════════════════════════════════════════
        //  Helpers (unchanged from previous version)
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