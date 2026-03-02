// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;
using System.Collections.Generic;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Abstraction for DIA classifiers used in the iterative FDR workflow.
    /// Both DiaLinearDiscriminant and DiaGradientBoostedClassifier implement this.
    /// This allows DiaFdrEngine to be agnostic to the classifier type.
    /// </summary>
    public interface IDiaClassifier
    {
        /// <summary>Number of features expected per sample</summary>
        int FeatureCount { get; }

        /// <summary>Whether the model has been trained</summary>
        bool IsTrained { get; }

        /// <summary>
        /// Trains on row-major feature matrix with binary labels.
        /// Labels: 1.0f = target, 0.0f = decoy.
        /// </summary>
        void Train(ReadOnlySpan<float> features, ReadOnlySpan<float> labels, int numSamples,
            ReadOnlySpan<float> sampleWeights = default);

        /// <summary>
        /// Batch prediction. Writes P(target | features) to scores span.
        /// </summary>
        void PredictBatch(ReadOnlySpan<float> features, int numSamples, Span<float> scores);

        /// <summary>Computes AUC on the given labeled dataset.</summary>
        float ComputeAuc(ReadOnlySpan<float> features, ReadOnlySpan<float> labels, int numSamples);
    }

    /// <summary>
    /// Specifies which classifier to use in the DIA FDR workflow.
    /// </summary>
    public enum DiaClassifierType
    {
        /// <summary>Linear Discriminant Analysis (existing, fast, but limited to linear boundaries)</summary>
        LinearDiscriminant,

        /// <summary>Gradient-Boosted Trees (new, non-linear, captures feature interactions)</summary>
        GradientBoostedTree
    }

    /// <summary>
    /// Iterative semi-supervised FDR estimation for DIA search results.
    /// 
    /// Workflow:
    ///   1. Initial labeling: targets with ApexScore > threshold → positive, decoys → negative
    ///   2. Train classifier (LDA or GBT) on (features, labels)
    ///   3. Score all results with classifier → ClassifierScore
    ///   4. Run target-decoy FDR estimation → q-values
    ///   5. Re-label: confident targets (q ≤ 0.01) → positive
    ///   6. Repeat until convergence or max iterations
    /// 
    /// Phase 14 changes:
    ///   - Added IDiaClassifier interface for classifier abstraction
    ///   - Added DiaGradientBoostedClassifier as non-linear option
    ///   - GBT is now the default classifier (DiaClassifierType.GradientBoostedTree)
    ///   - Added early stopping when ID count stabilizes (±0.5%)
    ///   - Added class balancing via sample weights
    ///   - Training AUC tracked per iteration for convergence monitoring
    /// </summary>
    public sealed class DiaFdrEngine
    {
        private readonly DiaSearchParameters _parameters;

        /// <summary>Results of the most recent FDR estimation run</summary>
        public DiaFdrRunSummary LastRunSummary { get; private set; }

        public DiaFdrEngine(DiaSearchParameters parameters)
        {
            _parameters = parameters ?? throw new ArgumentNullException(nameof(parameters));
        }

        /// <summary>
        /// Runs iterative semi-supervised FDR estimation.
        /// 
        /// Modifies DiaSearchResult.ClassifierScore and DiaSearchResult.FdrInfo in-place.
        /// Returns count of results passing 1% FDR.
        /// </summary>
        /// <param name="results">All search results (targets + decoys). Modified in-place.</param>
        /// <param name="featureVectors">Pre-computed feature vectors (parallel to results).</param>
        /// <param name="classifierType">Which classifier to use. Default: GBT.</param>
        /// <param name="maxIterations">Maximum FDR iterations. Default: 10.</param>
        /// <returns>Number of target results at q ≤ 0.01</returns>
        public int EstimateFdr(
            List<DiaSearchResult> results,
            DiaFeatureVector[] featureVectors,
            DiaClassifierType classifierType = DiaClassifierType.GradientBoostedTree,
            int maxIterations = 10)
        {
            if (results == null || results.Count == 0)
                return 0;
            if (featureVectors == null || featureVectors.Length != results.Count)
                throw new ArgumentException("featureVectors must be parallel to results");

            int featureCount = DiaFeatureVector.ClassifierFeatureCount;
            int n = results.Count;

            // === Build flat feature matrix (row-major) ===
            var featureMatrix = ArrayPool<float>.Shared.Rent(n * featureCount);
            var labels = ArrayPool<float>.Shared.Rent(n);
            var scores = ArrayPool<float>.Shared.Rent(n);
            var weights = ArrayPool<float>.Shared.Rent(n);

            try
            {
                // Pack feature vectors into flat array
                for (int i = 0; i < n; i++)
                {
                    featureVectors[i].WriteTo(
                        featureMatrix.AsSpan(i * featureCount, featureCount));
                }

                // Count targets and decoys
                int targetCount = 0, decoyCount = 0;
                for (int i = 0; i < n; i++)
                {
                    if (results[i].IsDecoy) decoyCount++;
                    else targetCount++;
                }

                if (targetCount == 0 || decoyCount == 0)
                    return 0;

                // === Initial labeling ===
                // Targets with ApexScore > 0.5 (or any reasonable threshold) → 1, decoys → 0
                // Targets below threshold also get 0 initially (conservative start)
                float initialThreshold = 0.5f;
                for (int i = 0; i < n; i++)
                {
                    if (results[i].IsDecoy)
                        labels[i] = 0f;
                    else
                        labels[i] = results[i].ApexScore > initialThreshold ? 1f : 0f;
                }

                // Check we have enough positive labels to start
                int posCount = 0;
                for (int i = 0; i < n; i++)
                    if (labels[i] > 0.5f) posCount++;

                if (posCount < 100)
                {
                    // Fallback: use all targets as positive (less reliable but avoids empty training)
                    for (int i = 0; i < n; i++)
                        labels[i] = results[i].IsDecoy ? 0f : 1f;
                    posCount = targetCount;
                }

                // === Iterative classifier training ===
                var summary = new DiaFdrRunSummary
                {
                    ClassifierType = classifierType,
                    TotalResults = n,
                    TargetCount = targetCount,
                    DecoyCount = decoyCount,
                    IterationDetails = new List<IterationDetail>()
                };

                int prevIdsAt1Pct = -1;
                IDiaClassifier classifier = CreateClassifier(classifierType, featureCount);

                for (int iter = 0; iter < maxIterations; iter++)
                {
                    // Compute class-balancing weights
                    int iterPosCount = 0, iterNegCount = 0;
                    for (int i = 0; i < n; i++)
                    {
                        if (labels[i] > 0.5f) iterPosCount++;
                        else iterNegCount++;
                    }

                    if (iterPosCount == 0 || iterNegCount == 0)
                        break;

                    // Balance: weight minority class up
                    float posWeight = (float)n / (2f * iterPosCount);
                    float negWeight = (float)n / (2f * iterNegCount);
                    for (int i = 0; i < n; i++)
                        weights[i] = labels[i] > 0.5f ? posWeight : negWeight;

                    // Train classifier
                    classifier = CreateClassifier(classifierType, featureCount);
                    classifier.Train(
                        featureMatrix.AsSpan(0, n * featureCount),
                        labels.AsSpan(0, n),
                        n,
                        weights.AsSpan(0, n));

                    // Predict scores
                    classifier.PredictBatch(
                        featureMatrix.AsSpan(0, n * featureCount),
                        n,
                        scores.AsSpan(0, n));

                    // Store scores on results
                    for (int i = 0; i < n; i++)
                        results[i].ClassifierScore = scores[i];

                    // Compute AUC for monitoring
                    float auc = classifier.ComputeAuc(
                        featureMatrix.AsSpan(0, n * featureCount),
                        labels.AsSpan(0, n),
                        n);

                    // Run target-decoy FDR
                    int idsAt1Pct = ComputeTargetDecoyFdr(results);

                    var detail = new IterationDetail
                    {
                        Iteration = iter + 1,
                        PositiveLabels = iterPosCount,
                        NegativeLabels = iterNegCount,
                        Auc = auc,
                        IdsAtOnePercentFdr = idsAt1Pct
                    };
                    summary.IterationDetails.Add(detail);

                    // Check convergence: IDs stabilized within 0.5%
                    if (prevIdsAt1Pct > 0 && iter > 1)
                    {
                        float changePct = MathF.Abs(idsAt1Pct - prevIdsAt1Pct) / (float)prevIdsAt1Pct;
                        if (changePct < 0.005f)
                        {
                            summary.ConvergedAtIteration = iter + 1;
                            break;
                        }
                    }

                    prevIdsAt1Pct = idsAt1Pct;

                    // Re-label for next iteration: confident targets at q ≤ 0.01 → positive
                    // All decoys remain 0. Non-confident targets → 0.
                    for (int i = 0; i < n; i++)
                    {
                        if (results[i].IsDecoy)
                        {
                            labels[i] = 0f;
                        }
                        else
                        {
                            labels[i] = (results[i].FdrInfo != null && results[i].FdrInfo.QValue <= 0.01f)
                                ? 1f : 0f;
                        }
                    }
                }

                summary.FinalIdsAtOnePercentFdr = prevIdsAt1Pct;
                LastRunSummary = summary;
                return prevIdsAt1Pct;
            }
            finally
            {
                ArrayPool<float>.Shared.Return(featureMatrix);
                ArrayPool<float>.Shared.Return(labels);
                ArrayPool<float>.Shared.Return(scores);
                ArrayPool<float>.Shared.Return(weights);
            }
        }

        /// <summary>
        /// Target-decoy FDR estimation using the ClassifierScore.
        /// Sorts all results by ClassifierScore descending, counts targets/decoys,
        /// computes q-values. Assigns DiaFdrInfo to each result.
        /// </summary>
        /// <returns>Count of target results at q ≤ 0.01</returns>
        private static int ComputeTargetDecoyFdr(List<DiaSearchResult> results)
        {
            int n = results.Count;

            // Sort by ClassifierScore descending
            var indices = ArrayPool<int>.Shared.Rent(n);
            try
            {
                for (int i = 0; i < n; i++) indices[i] = i;
                Array.Sort(indices, 0, n, new ClassifierScoreComparer(results));

                // Forward pass: compute FDR at each rank
                int cumulTargets = 0;
                int cumulDecoys = 0;
                var fdrs = ArrayPool<float>.Shared.Rent(n);
                try
                {
                    for (int rank = 0; rank < n; rank++)
                    {
                        int idx = indices[rank];
                        if (results[idx].IsDecoy)
                            cumulDecoys++;
                        else
                            cumulTargets++;

                        float fdr = cumulTargets > 0
                            ? (float)cumulDecoys / (float)cumulTargets
                            : 1f;
                        fdrs[rank] = Math.Min(fdr, 1f);
                    }

                    // Backward pass: convert FDR → q-value (monotone decreasing from tail)
                    float minFdrSoFar = 1f;
                    for (int rank = n - 1; rank >= 0; rank--)
                    {
                        minFdrSoFar = MathF.Min(minFdrSoFar, fdrs[rank]);
                        fdrs[rank] = minFdrSoFar;
                    }

                    // Assign q-values
                    int idsAt1Pct = 0;
                    for (int rank = 0; rank < n; rank++)
                    {
                        int idx = indices[rank];
                        var result = results[idx];

                        if (result.FdrInfo == null)
                            result.FdrInfo = new DiaFdrInfo();

                        result.FdrInfo.QValue = fdrs[rank];
                        result.FdrInfo.CumulativeTarget = 0; // set below
                        result.FdrInfo.CumulativeDecoy = 0;

                        if (!result.IsDecoy && fdrs[rank] <= 0.01f)
                            idsAt1Pct++;
                    }

                    // Set cumulative counts
                    cumulTargets = 0;
                    cumulDecoys = 0;
                    for (int rank = 0; rank < n; rank++)
                    {
                        int idx = indices[rank];
                        if (results[idx].IsDecoy)
                            cumulDecoys++;
                        else
                            cumulTargets++;

                        if (results[idx].FdrInfo != null)
                        {
                            results[idx].FdrInfo.CumulativeTarget = cumulTargets;
                            results[idx].FdrInfo.CumulativeDecoy = cumulDecoys;
                        }
                    }

                    return idsAt1Pct;
                }
                finally
                {
                    ArrayPool<float>.Shared.Return(fdrs);
                }
            }
            finally
            {
                ArrayPool<int>.Shared.Return(indices);
            }
        }

        /// <summary>
        /// Creates a classifier instance based on the selected type.
        /// </summary>
        private static IDiaClassifier CreateClassifier(DiaClassifierType type, int featureCount)
        {
            return type switch
            {
                DiaClassifierType.GradientBoostedTree => new DiaGbtClassifierAdapter(
                    new DiaGradientBoostedClassifier(
                        featureCount: featureCount,
                        numTrees: 100,
                        maxDepth: 4,
                        learningRate: 0.1f,
                        minSamplesLeaf: 20,
                        l2Regularization: 1.0f,
                        numBins: 64,
                        subsampleFraction: 0.8f)),

                DiaClassifierType.LinearDiscriminant => new DiaLdaClassifierAdapter(
                    featureCount),

                _ => throw new ArgumentException($"Unknown classifier type: {type}")
            };
        }

        private sealed class ClassifierScoreComparer : System.Collections.Generic.IComparer<int>
        {
            private readonly List<DiaSearchResult> _results;
            public ClassifierScoreComparer(List<DiaSearchResult> r) { _results = r; }
            public int Compare(int a, int b) =>
                _results[b].ClassifierScore.CompareTo(_results[a].ClassifierScore);
        }

        #region Summary Types

        /// <summary>Summary of an FDR estimation run.</summary>
        public class DiaFdrRunSummary
        {
            public DiaClassifierType ClassifierType { get; set; }
            public int TotalResults { get; set; }
            public int TargetCount { get; set; }
            public int DecoyCount { get; set; }
            public int FinalIdsAtOnePercentFdr { get; set; }
            public int ConvergedAtIteration { get; set; }
            public List<IterationDetail> IterationDetails { get; set; }

            public override string ToString()
            {
                var sb = new System.Text.StringBuilder();
                sb.AppendLine($"DIA FDR Summary: {ClassifierType}");
                sb.AppendLine($"  Total: {TotalResults} ({TargetCount} targets, {DecoyCount} decoys)");
                sb.AppendLine($"  Final IDs at 1% FDR: {FinalIdsAtOnePercentFdr}");
                if (ConvergedAtIteration > 0)
                    sb.AppendLine($"  Converged at iteration {ConvergedAtIteration}");
                if (IterationDetails != null)
                {
                    foreach (var d in IterationDetails)
                        sb.AppendLine($"  Iter {d.Iteration}: AUC={d.Auc:F4}, " +
                                      $"IDs@1%={d.IdsAtOnePercentFdr}, " +
                                      $"pos={d.PositiveLabels}, neg={d.NegativeLabels}");
                }
                return sb.ToString();
            }
        }

        /// <summary>Details of a single FDR iteration.</summary>
        public class IterationDetail
        {
            public int Iteration { get; set; }
            public int PositiveLabels { get; set; }
            public int NegativeLabels { get; set; }
            public float Auc { get; set; }
            public int IdsAtOnePercentFdr { get; set; }
        }

        #endregion
    }

    #region Classifier Adapters

    /// <summary>
    /// Adapter wrapping DiaGradientBoostedClassifier to IDiaClassifier interface.
    /// </summary>
    internal sealed class DiaGbtClassifierAdapter : IDiaClassifier
    {
        private readonly DiaGradientBoostedClassifier _gbt;

        public DiaGbtClassifierAdapter(DiaGradientBoostedClassifier gbt)
        {
            _gbt = gbt ?? throw new ArgumentNullException(nameof(gbt));
        }

        public int FeatureCount => _gbt.FeatureCount;
        public bool IsTrained => _gbt.IsTrained;

        public void Train(ReadOnlySpan<float> features, ReadOnlySpan<float> labels,
            int numSamples, ReadOnlySpan<float> sampleWeights = default)
            => _gbt.Train(features, labels, numSamples, sampleWeights);

        public void PredictBatch(ReadOnlySpan<float> features, int numSamples, Span<float> scores)
            => _gbt.PredictBatch(features, numSamples, scores);

        public float ComputeAuc(ReadOnlySpan<float> features, ReadOnlySpan<float> labels, int numSamples)
            => _gbt.ComputeAuc(features, labels, numSamples);
    }

    /// <summary>
    /// Adapter wrapping the existing DiaLinearDiscriminant to IDiaClassifier interface.
    /// This preserves backward compatibility — existing LDA can still be used.
    /// 
    /// NOTE: DiaLinearDiscriminant must be updated to implement Train/PredictBatch
    /// with the Span-based signatures, or this adapter bridges the gap.
    /// </summary>
    internal sealed class DiaLdaClassifierAdapter : IDiaClassifier
    {
        private readonly int _featureCount;
        private float[] _weights;
        private float _bias;
        private bool _isTrained;

        public DiaLdaClassifierAdapter(int featureCount)
        {
            _featureCount = featureCount;
        }

        public int FeatureCount => _featureCount;
        public bool IsTrained => _isTrained;

        /// <summary>
        /// Simple Fisher LDA: projects features onto the direction that maximizes
        /// between-class / within-class variance ratio.
        /// </summary>
        public void Train(ReadOnlySpan<float> features, ReadOnlySpan<float> labels,
            int numSamples, ReadOnlySpan<float> sampleWeights = default)
        {
            // Compute class means
            var mean0 = new float[_featureCount];
            var mean1 = new float[_featureCount];
            float count0 = 0f, count1 = 0f;

            for (int i = 0; i < numSamples; i++)
            {
                int offset = i * _featureCount;
                float w = sampleWeights.IsEmpty ? 1f : sampleWeights[i];

                if (labels[i] > 0.5f)
                {
                    count1 += w;
                    for (int f = 0; f < _featureCount; f++)
                    {
                        float v = features[offset + f];
                        if (float.IsNaN(v)) v = 0f;
                        mean1[f] += v * w;
                    }
                }
                else
                {
                    count0 += w;
                    for (int f = 0; f < _featureCount; f++)
                    {
                        float v = features[offset + f];
                        if (float.IsNaN(v)) v = 0f;
                        mean0[f] += v * w;
                    }
                }
            }

            if (count0 == 0f || count1 == 0f)
            {
                _weights = new float[_featureCount];
                _bias = 0f;
                _isTrained = true;
                return;
            }

            for (int f = 0; f < _featureCount; f++)
            {
                mean0[f] /= count0;
                mean1[f] /= count1;
            }

            // Compute pooled within-class variance per feature (diagonal approximation)
            var variance = new float[_featureCount];
            for (int i = 0; i < numSamples; i++)
            {
                int offset = i * _featureCount;
                float w = sampleWeights.IsEmpty ? 1f : sampleWeights[i];
                var mean = labels[i] > 0.5f ? mean1 : mean0;
                for (int f = 0; f < _featureCount; f++)
                {
                    float v = features[offset + f];
                    if (float.IsNaN(v)) v = 0f;
                    float diff = v - mean[f];
                    variance[f] += w * diff * diff;
                }
            }

            float totalWeight = count0 + count1;
            _weights = new float[_featureCount];
            _bias = 0f;
            for (int f = 0; f < _featureCount; f++)
            {
                variance[f] /= totalWeight;
                variance[f] = MathF.Max(variance[f], 1e-8f); // avoid division by zero
                _weights[f] = (mean1[f] - mean0[f]) / variance[f];
                _bias -= 0.5f * _weights[f] * (mean1[f] + mean0[f]);
            }

            _isTrained = true;
        }

        public void PredictBatch(ReadOnlySpan<float> features, int numSamples, Span<float> scores)
        {
            for (int i = 0; i < numSamples; i++)
            {
                float logit = _bias;
                int offset = i * _featureCount;
                for (int f = 0; f < _featureCount; f++)
                {
                    float v = features[offset + f];
                    if (float.IsNaN(v)) v = 0f;
                    logit += _weights[f] * v;
                }
                // Sigmoid to get probability
                scores[i] = Sigmoid(logit);
            }
        }

        public float ComputeAuc(ReadOnlySpan<float> features, ReadOnlySpan<float> labels, int numSamples)
        {
            var scores = ArrayPool<float>.Shared.Rent(numSamples);
            try
            {
                PredictBatch(features, numSamples, scores);
                return DiaGradientBoostedClassifier.CalculateAuc(
                    scores.AsSpan(0, numSamples), labels, numSamples);
            }
            finally
            {
                ArrayPool<float>.Shared.Return(scores);
            }
        }

        private static float Sigmoid(float x)
        {
            if (x > 20f) return 1f;
            if (x < -20f) return 0f;
            return 1f / (1f + MathF.Exp(-x));
        }
    }

    #endregion
}