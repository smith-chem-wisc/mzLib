// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;
using System.Runtime.CompilerServices;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Gradient-boosted decision tree classifier for DIA target/decoy discrimination.
    /// 
    /// Replaces the linear LDA classifier (DiaLinearDiscriminant) with a non-linear model
    /// that can capture complex feature interactions in the 28-dimensional feature space.
    /// This is the single highest-impact change for improving FDR-controlled identification counts.
    /// 
    /// Design choices:
    ///   - Self-contained: no ML.NET, no external dependencies
    ///   - Gradient boosting with binary cross-entropy loss (logistic regression trees)
    ///   - Shallow trees (depth 4-5) to avoid overfitting with ~50K samples
    ///   - L2 regularization (shrinkage) + min samples per leaf
    ///   - Histogram-based split finding for O(n × features × bins) training
    ///   - Prediction is O(trees × depth) per sample — extremely fast
    ///   - Uses ArrayPool for all temporary allocations during training
    ///   - Float precision throughout (matches DiaFeatureVector)
    /// 
    /// The classifier outputs a probability estimate P(target | features), which is stored
    /// as ClassifierScore on DiaSearchResult. DiaFdrEngine uses this for target-decoy FDR.
    /// 
    /// In the iterative semi-supervised FDR workflow:
    ///   1. Initial: all targets with score > threshold are "positive", all decoys "negative"
    ///   2. Train GBT on (features, label) pairs
    ///   3. Predict ClassifierScore for all results
    ///   4. Run target-decoy FDR → q-values
    ///   5. Retrain using confident targets (q ≤ 0.01) as positive
    ///   6. Repeat until convergence
    /// </summary>
    public sealed class DiaGradientBoostedClassifier
    {
        // === Model parameters (set at construction, fixed during training) ===
        private readonly int _maxDepth;
        private readonly int _numTrees;
        private readonly float _learningRate;
        private readonly int _minSamplesLeaf;
        private readonly float _l2Regularization;
        private readonly int _numBins;
        private readonly float _subsampleFraction;
        private readonly int _featureCount;

        // === Trained model state ===
        private DecisionTree[] _trees;
        private float _initialPrediction; // log-odds bias (intercept)
        private bool _isTrained;

        // === Histogram binning state (computed once from training data) ===
        private float[][] _binEdges; // [feature][bin] — thresholds for histogram bins

        /// <summary>Number of features expected per sample</summary>
        public int FeatureCount => _featureCount;

        /// <summary>Whether the model has been trained</summary>
        public bool IsTrained => _isTrained;

        /// <summary>Number of trees in the ensemble</summary>
        public int TreeCount => _numTrees;

        /// <summary>
        /// Feature importance scores (sum of gain across all splits for each feature).
        /// Available after training. Higher = more important.
        /// </summary>
        public float[] FeatureImportances { get; private set; }

        /// <summary>
        /// Creates a gradient-boosted classifier with the specified hyperparameters.
        /// </summary>
        /// <param name="featureCount">Number of input features (must match DiaFeatureVector.ClassifierFeatureCount)</param>
        /// <param name="numTrees">Number of boosting rounds (trees). More = better fit but slower. Default: 100</param>
        /// <param name="maxDepth">Maximum tree depth. 4-5 captures non-linear interactions without overfitting. Default: 4</param>
        /// <param name="learningRate">Shrinkage factor per tree. Lower = more robust but needs more trees. Default: 0.1</param>
        /// <param name="minSamplesLeaf">Minimum samples in a leaf node. Regularization against overfitting. Default: 20</param>
        /// <param name="l2Regularization">L2 penalty on leaf values (lambda in XGBoost). Default: 1.0</param>
        /// <param name="numBins">Number of histogram bins for split finding. More = finer splits. Default: 64</param>
        /// <param name="subsampleFraction">Fraction of samples used per tree (stochastic gradient boosting). Default: 0.8</param>
        public DiaGradientBoostedClassifier(
            int featureCount,
            int numTrees = 100,
            int maxDepth = 4,
            float learningRate = 0.1f,
            int minSamplesLeaf = 20,
            float l2Regularization = 1.0f,
            int numBins = 64,
            float subsampleFraction = 0.8f)
        {
            if (featureCount <= 0) throw new ArgumentOutOfRangeException(nameof(featureCount));
            if (numTrees <= 0) throw new ArgumentOutOfRangeException(nameof(numTrees));
            if (maxDepth <= 0 || maxDepth > 12) throw new ArgumentOutOfRangeException(nameof(maxDepth), "maxDepth must be 1-12");
            if (learningRate <= 0f || learningRate > 1f) throw new ArgumentOutOfRangeException(nameof(learningRate));
            if (minSamplesLeaf <= 0) throw new ArgumentOutOfRangeException(nameof(minSamplesLeaf));
            if (l2Regularization < 0f) throw new ArgumentOutOfRangeException(nameof(l2Regularization));
            if (numBins < 4 || numBins > 256) throw new ArgumentOutOfRangeException(nameof(numBins), "numBins must be 4-256");
            if (subsampleFraction <= 0f || subsampleFraction > 1f) throw new ArgumentOutOfRangeException(nameof(subsampleFraction));

            _featureCount = featureCount;
            _numTrees = numTrees;
            _maxDepth = maxDepth;
            _learningRate = learningRate;
            _minSamplesLeaf = minSamplesLeaf;
            _l2Regularization = l2Regularization;
            _numBins = numBins;
            _subsampleFraction = subsampleFraction;
        }

        /// <summary>
        /// Trains the classifier on feature matrix + binary labels.
        /// 
        /// Features are stored in row-major flat array: features[sample * featureCount + feature].
        /// Labels: 1.0f = target (positive), 0.0f = decoy (negative).
        /// 
        /// Training uses histogram-based gradient boosting with logistic loss.
        /// Time complexity: O(numTrees × numSamples × featureCount × log(numBins))
        /// </summary>
        /// <param name="features">Row-major feature matrix [numSamples × featureCount]</param>
        /// <param name="labels">Binary labels, length = numSamples</param>
        /// <param name="numSamples">Number of training samples</param>
        /// <param name="sampleWeights">Optional per-sample weights. Null = uniform. Length must equal numSamples.</param>
        public void Train(ReadOnlySpan<float> features, ReadOnlySpan<float> labels, int numSamples,
            ReadOnlySpan<float> sampleWeights = default)
        {
            if (features.Length != numSamples * _featureCount)
                throw new ArgumentException(
                    $"features length ({features.Length}) must equal numSamples ({numSamples}) × featureCount ({_featureCount})");
            if (labels.Length != numSamples)
                throw new ArgumentException($"labels length ({labels.Length}) must equal numSamples ({numSamples})");
            if (numSamples < 2 * _minSamplesLeaf)
                throw new ArgumentException($"Need at least {2 * _minSamplesLeaf} samples (have {numSamples})");
            if (!sampleWeights.IsEmpty && sampleWeights.Length != numSamples)
                throw new ArgumentException($"sampleWeights length must equal numSamples");

            // === Step 1: Build histogram bin edges from data ===
            _binEdges = BuildBinEdges(features, numSamples);

            // === Step 2: Bin the features into integer bin indices ===
            var binnedFeatures = ArrayPool<byte>.Shared.Rent(numSamples * _featureCount);
            try
            {
                BinFeatures(features, numSamples, binnedFeatures);

                // === Step 3: Initialize predictions with log-odds of positive rate ===
                float posCount = 0f;
                for (int i = 0; i < numSamples; i++)
                    posCount += labels[i];
                float negCount = numSamples - posCount;
                _initialPrediction = MathF.Log((posCount + 1f) / (negCount + 1f));

                // === Step 4: Allocate working arrays ===
                var predictions = ArrayPool<float>.Shared.Rent(numSamples);
                var gradients = ArrayPool<float>.Shared.Rent(numSamples);
                var hessians = ArrayPool<float>.Shared.Rent(numSamples);
                var sampleMask = ArrayPool<bool>.Shared.Rent(numSamples);
                var importances = new float[_featureCount];

                try
                {
                    // Initialize predictions
                    for (int i = 0; i < numSamples; i++)
                        predictions[i] = _initialPrediction;

                    _trees = new DecisionTree[_numTrees];
                    var rng = new Random(42);

                    // === Step 5: Boosting iterations ===
                    for (int t = 0; t < _numTrees; t++)
                    {
                        ComputeGradientsAndHessians(predictions, labels, numSamples,
                            sampleWeights, gradients, hessians);

                        // Stochastic subsampling
                        int subsampleCount;
                        if (_subsampleFraction < 1.0f)
                        {
                            subsampleCount = 0;
                            for (int i = 0; i < numSamples; i++)
                            {
                                sampleMask[i] = rng.NextDouble() < _subsampleFraction;
                                if (sampleMask[i]) subsampleCount++;
                            }
                        }
                        else
                        {
                            for (int i = 0; i < numSamples; i++)
                                sampleMask[i] = true;
                            subsampleCount = numSamples;
                        }

                        if (subsampleCount < 2 * _minSamplesLeaf)
                        {
                            for (int i = 0; i < numSamples; i++)
                                sampleMask[i] = true;
                            subsampleCount = numSamples;
                        }

                        _trees[t] = BuildTree(binnedFeatures, gradients, hessians,
                            sampleMask, numSamples, subsampleCount, importances);

                        // Update predictions
                        for (int i = 0; i < numSamples; i++)
                        {
                            float leafValue = _trees[t].Predict(binnedFeatures, i, _featureCount);
                            predictions[i] += _learningRate * leafValue;
                        }
                    }

                    // Normalize importances
                    float maxImp = 0f;
                    for (int f = 0; f < _featureCount; f++)
                        if (importances[f] > maxImp) maxImp = importances[f];
                    if (maxImp > 0f)
                        for (int f = 0; f < _featureCount; f++)
                            importances[f] /= maxImp;

                    FeatureImportances = importances;
                    _isTrained = true;
                }
                finally
                {
                    ArrayPool<float>.Shared.Return(predictions);
                    ArrayPool<float>.Shared.Return(gradients);
                    ArrayPool<float>.Shared.Return(hessians);
                    ArrayPool<bool>.Shared.Return(sampleMask);
                }
            }
            finally
            {
                ArrayPool<byte>.Shared.Return(binnedFeatures);
            }
        }

        /// <summary>
        /// Predicts probability P(target | features) for a single sample.
        /// Returns value in [0, 1]. Higher = more likely to be a true target.
        /// Stored as DiaSearchResult.ClassifierScore.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public float PredictProbability(ReadOnlySpan<float> features)
        {
            if (!_isTrained) throw new InvalidOperationException("Model has not been trained");
            if (features.Length != _featureCount)
                throw new ArgumentException($"Expected {_featureCount} features, got {features.Length}");

            float logit = _initialPrediction;

            // Bin the features
            Span<byte> binned = stackalloc byte[_featureCount];
            for (int f = 0; f < _featureCount; f++)
                binned[f] = FindBin(_binEdges[f], features[f]);

            for (int t = 0; t < _trees.Length; t++)
                logit += _learningRate * _trees[t].PredictFromBinned(binned);

            return Sigmoid(logit);
        }

        /// <summary>
        /// Batch prediction for multiple samples. More efficient than individual calls.
        /// Features in row-major layout: features[sample * featureCount + feature].
        /// </summary>
        public void PredictBatch(ReadOnlySpan<float> features, int numSamples, Span<float> scores)
        {
            if (!_isTrained) throw new InvalidOperationException("Model has not been trained");
            if (features.Length != numSamples * _featureCount)
                throw new ArgumentException("features size mismatch");
            if (scores.Length < numSamples)
                throw new ArgumentException("scores buffer too small");

            var binnedAll = ArrayPool<byte>.Shared.Rent(numSamples * _featureCount);
            try
            {
                for (int i = 0; i < numSamples; i++)
                {
                    int rowOffset = i * _featureCount;
                    for (int f = 0; f < _featureCount; f++)
                        binnedAll[rowOffset + f] = FindBin(_binEdges[f], features[rowOffset + f]);
                }

                for (int i = 0; i < numSamples; i++)
                {
                    float logit = _initialPrediction;
                    for (int t = 0; t < _trees.Length; t++)
                        logit += _learningRate * _trees[t].Predict(binnedAll, i, _featureCount);
                    scores[i] = Sigmoid(logit);
                }
            }
            finally
            {
                ArrayPool<byte>.Shared.Return(binnedAll);
            }
        }

        /// <summary>
        /// Computes training AUC on the given data (for monitoring convergence).
        /// </summary>
        public float ComputeAuc(ReadOnlySpan<float> features, ReadOnlySpan<float> labels, int numSamples)
        {
            if (!_isTrained) throw new InvalidOperationException("Model has not been trained");

            var scores = ArrayPool<float>.Shared.Rent(numSamples);
            try
            {
                PredictBatch(features, numSamples, scores);
                return CalculateAuc(scores.AsSpan(0, numSamples), labels, numSamples);
            }
            finally
            {
                ArrayPool<float>.Shared.Return(scores);
            }
        }

        /// <summary>
        /// Returns a human-readable summary of feature importances.
        /// </summary>
        public string GetImportanceSummary(string[] featureNames)
        {
            if (FeatureImportances == null) return "Model not trained";
            if (featureNames.Length != _featureCount)
                throw new ArgumentException("featureNames length mismatch");

            var indices = new int[_featureCount];
            for (int i = 0; i < _featureCount; i++) indices[i] = i;
            Array.Sort(indices, (a, b) => FeatureImportances[b].CompareTo(FeatureImportances[a]));

            var sb = new System.Text.StringBuilder();
            sb.AppendLine("Feature Importances (GBT):");
            for (int i = 0; i < _featureCount; i++)
            {
                int idx = indices[i];
                sb.AppendLine($"  {FeatureImportances[idx]:F4}  {featureNames[idx]}");
            }
            return sb.ToString();
        }

        #region Private: Training Internals

        /// <summary>
        /// Builds quantile-based bin edges for each feature using a sorted subsample.
        /// </summary>
        private float[][] BuildBinEdges(ReadOnlySpan<float> features, int numSamples)
        {
            var edges = new float[_featureCount][];
            int subsample = Math.Min(numSamples, 10000);
            var temp = ArrayPool<float>.Shared.Rent(subsample);

            try
            {
                var rng = new Random(123);
                for (int f = 0; f < _featureCount; f++)
                {
                    if (subsample == numSamples)
                    {
                        for (int i = 0; i < numSamples; i++)
                            temp[i] = features[i * _featureCount + f];
                    }
                    else
                    {
                        for (int i = 0; i < subsample; i++)
                        {
                            int idx = rng.Next(numSamples);
                            temp[i] = features[idx * _featureCount + f];
                        }
                    }

                    // Replace NaN with 0 for binning purposes
                    for (int i = 0; i < subsample; i++)
                        if (float.IsNaN(temp[i])) temp[i] = 0f;

                    Array.Sort(temp, 0, subsample);

                    int numEdges = _numBins - 1;
                    edges[f] = new float[numEdges];
                    for (int b = 0; b < numEdges; b++)
                    {
                        float quantile = (b + 1.0f) / _numBins;
                        int idx = (int)(quantile * (subsample - 1));
                        edges[f][b] = temp[idx];
                    }

                    DeduplicateEdges(edges[f]);
                }
            }
            finally
            {
                ArrayPool<float>.Shared.Return(temp);
            }

            return edges;
        }

        private static void DeduplicateEdges(float[] edges)
        {
            for (int i = 1; i < edges.Length; i++)
            {
                if (edges[i] <= edges[i - 1])
                    edges[i] = float.PositiveInfinity;
            }
        }

        private void BinFeatures(ReadOnlySpan<float> features, int numSamples, byte[] binnedOut)
        {
            for (int i = 0; i < numSamples; i++)
            {
                int rowOffset = i * _featureCount;
                for (int f = 0; f < _featureCount; f++)
                    binnedOut[rowOffset + f] = FindBin(_binEdges[f], features[rowOffset + f]);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static byte FindBin(float[] edges, float value)
        {
            if (float.IsNaN(value)) return 0;

            int lo = 0, hi = edges.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) >> 1;
                if (value > edges[mid])
                    lo = mid + 1;
                else
                    hi = mid;
            }
            return (byte)lo;
        }

        private static void ComputeGradientsAndHessians(
            float[] predictions, ReadOnlySpan<float> labels, int numSamples,
            ReadOnlySpan<float> weights, float[] gradients, float[] hessians)
        {
            bool hasWeights = !weights.IsEmpty;
            for (int i = 0; i < numSamples; i++)
            {
                float p = Sigmoid(predictions[i]);
                float g = labels[i] - p;
                float h = p * (1f - p);
                h = MathF.Max(h, 1e-7f);

                if (hasWeights)
                {
                    g *= weights[i];
                    h *= weights[i];
                }

                gradients[i] = g;
                hessians[i] = h;
            }
        }

        /// <summary>
        /// Builds a single regression tree on the current gradients/hessians.
        /// Histogram-based split finding: O(samples × features × bins).
        /// Stack-based (non-recursive) tree construction.
        /// </summary>
        private DecisionTree BuildTree(
            byte[] binnedFeatures, float[] gradients, float[] hessians,
            bool[] sampleMask, int totalSamples, int activeSamples,
            float[] importanceAccumulator)
        {
            int maxNodes = (1 << (_maxDepth + 1)) - 1;
            var nodes = new TreeNode[maxNodes];
            int nodeCount = 0;

            var sampleIndices = ArrayPool<int>.Shared.Rent(activeSamples);
            int sIdx = 0;
            for (int i = 0; i < totalSamples; i++)
            {
                if (sampleMask[i])
                    sampleIndices[sIdx++] = i;
            }

            var gradHist = new float[_featureCount * _numBins];
            var hessHist = new float[_featureCount * _numBins];

            // Stack: (nodeIndex, sampleStart, sampleCount, depth)
            var buildStack = new (int nodeIdx, int start, int count, int depth)[maxNodes];
            int stackSize = 0;

            nodes[0] = new TreeNode();
            nodeCount = 1;
            buildStack[stackSize++] = (0, 0, activeSamples, 0);

            while (stackSize > 0)
            {
                var (nodeIdx, start, count, depth) = buildStack[--stackSize];

                float sumGrad = 0f, sumHess = 0f;
                for (int s = start; s < start + count; s++)
                {
                    int idx = sampleIndices[s];
                    sumGrad += gradients[idx];
                    sumHess += hessians[idx];
                }
                nodes[nodeIdx].LeafValue = sumGrad / (sumHess + _l2Regularization);
                nodes[nodeIdx].IsLeaf = true;

                if (depth >= _maxDepth || count < 2 * _minSamplesLeaf)
                    continue;

                var bestSplit = FindBestSplit(binnedFeatures, gradients, hessians,
                    sampleIndices, start, count, sumGrad, sumHess,
                    gradHist, hessHist);

                if (!bestSplit.IsValid)
                    continue;

                importanceAccumulator[bestSplit.FeatureIndex] += bestSplit.Gain;

                int leftCount = PartitionSamples(binnedFeatures, sampleIndices,
                    start, count, bestSplit.FeatureIndex, bestSplit.SplitBin);

                if (leftCount < _minSamplesLeaf || (count - leftCount) < _minSamplesLeaf)
                    continue;

                int leftIdx = nodeCount++;
                int rightIdx = nodeCount++;
                if (rightIdx >= maxNodes)
                    continue;

                nodes[leftIdx] = new TreeNode();
                nodes[rightIdx] = new TreeNode();

                nodes[nodeIdx].IsLeaf = false;
                nodes[nodeIdx].SplitFeature = bestSplit.FeatureIndex;
                nodes[nodeIdx].SplitBin = bestSplit.SplitBin;
                nodes[nodeIdx].LeftChild = leftIdx;
                nodes[nodeIdx].RightChild = rightIdx;

                buildStack[stackSize++] = (rightIdx, start + leftCount, count - leftCount, depth + 1);
                buildStack[stackSize++] = (leftIdx, start, leftCount, depth + 1);
            }

            ArrayPool<int>.Shared.Return(sampleIndices);

            var compactNodes = new TreeNode[nodeCount];
            Array.Copy(nodes, compactNodes, nodeCount);
            return new DecisionTree(compactNodes);
        }

        private SplitInfo FindBestSplit(
            byte[] binnedFeatures, float[] gradients, float[] hessians,
            int[] sampleIndices, int start, int count,
            float totalGrad, float totalHess,
            float[] gradHist, float[] hessHist)
        {
            var best = new SplitInfo { Gain = 0f };
            float parentScore = (totalGrad * totalGrad) / (totalHess + _l2Regularization);

            for (int f = 0; f < _featureCount; f++)
            {
                int histOffset = f * _numBins;
                for (int b = 0; b < _numBins; b++)
                {
                    gradHist[histOffset + b] = 0f;
                    hessHist[histOffset + b] = 0f;
                }

                for (int s = start; s < start + count; s++)
                {
                    int idx = sampleIndices[s];
                    byte bin = binnedFeatures[idx * _featureCount + f];
                    gradHist[histOffset + bin] += gradients[idx];
                    hessHist[histOffset + bin] += hessians[idx];
                }

                float leftGrad = 0f, leftHess = 0f;
                for (int b = 0; b < _numBins - 1; b++)
                {
                    leftGrad += gradHist[histOffset + b];
                    leftHess += hessHist[histOffset + b];
                    float rightGrad = totalGrad - leftGrad;
                    float rightHess = totalHess - leftHess;

                    if (leftHess < _minSamplesLeaf * 0.01f || rightHess < _minSamplesLeaf * 0.01f)
                        continue;

                    float leftScore = (leftGrad * leftGrad) / (leftHess + _l2Regularization);
                    float rightScore = (rightGrad * rightGrad) / (rightHess + _l2Regularization);
                    float gain = leftScore + rightScore - parentScore;

                    if (gain > best.Gain)
                    {
                        best = new SplitInfo
                        {
                            FeatureIndex = f,
                            SplitBin = (byte)b,
                            Gain = gain,
                            IsValid = true
                        };
                    }
                }
            }

            return best;
        }

        private int PartitionSamples(byte[] binnedFeatures, int[] sampleIndices,
            int start, int count, int featureIndex, byte splitBin)
        {
            int left = start;
            int right = start + count - 1;

            while (left <= right)
            {
                int idx = sampleIndices[left];
                byte bin = binnedFeatures[idx * _featureCount + featureIndex];
                if (bin <= splitBin)
                {
                    left++;
                }
                else
                {
                    int tmp = sampleIndices[left];
                    sampleIndices[left] = sampleIndices[right];
                    sampleIndices[right] = tmp;
                    right--;
                }
            }

            return left - start;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float Sigmoid(float x)
        {
            if (x > 20f) return 1f;
            if (x < -20f) return 0f;
            return 1f / (1f + MathF.Exp(-x));
        }

        /// <summary>
        /// AUC via Wilcoxon-Mann-Whitney statistic (sorted rank method).
        /// </summary>
        public static float CalculateAuc(Span<float> scores, ReadOnlySpan<float> labels, int n)
        {
            var indices = ArrayPool<int>.Shared.Rent(n);
            try
            {
                for (int i = 0; i < n; i++) indices[i] = i;
                var scoresArr = scores.ToArray();
                Array.Sort(indices, 0, n, new ScoreComparer(scoresArr));

                long posCount = 0, negCount = 0;
                for (int i = 0; i < n; i++)
                {
                    if (labels[indices[i]] > 0.5f) posCount++;
                    else negCount++;
                }

                if (posCount == 0 || negCount == 0) return 0.5f;

                long tp = 0, areaAccum = 0;
                for (int i = 0; i < n; i++)
                {
                    if (labels[indices[i]] > 0.5f)
                        tp++;
                    else
                        areaAccum += tp;
                }

                return (float)areaAccum / (float)(posCount * negCount);
            }
            finally
            {
                ArrayPool<int>.Shared.Return(indices);
            }
        }

        private sealed class ScoreComparer : System.Collections.Generic.IComparer<int>
        {
            private readonly float[] _scores;
            public ScoreComparer(float[] scores) { _scores = scores; }
            public int Compare(int a, int b) => _scores[b].CompareTo(_scores[a]);
        }

        #endregion

        #region Internal Types

        private struct SplitInfo
        {
            public int FeatureIndex;
            public byte SplitBin;
            public float Gain;
            public bool IsValid;
        }

        private struct TreeNode
        {
            public bool IsLeaf;
            public int SplitFeature;
            public byte SplitBin;
            public int LeftChild;
            public int RightChild;
            public float LeafValue;
        }

        /// <summary>
        /// A complete decision tree. Array-of-struct nodes for cache-friendly traversal.
        /// </summary>
        private sealed class DecisionTree
        {
            private readonly TreeNode[] _nodes;

            public DecisionTree(TreeNode[] nodes) { _nodes = nodes; }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public float Predict(byte[] binnedFeatures, int sampleIndex, int featureCount)
            {
                int nodeIdx = 0;
                int rowOffset = sampleIndex * featureCount;

                while (!_nodes[nodeIdx].IsLeaf)
                {
                    ref readonly var node = ref _nodes[nodeIdx];
                    byte bin = binnedFeatures[rowOffset + node.SplitFeature];
                    nodeIdx = bin <= node.SplitBin ? node.LeftChild : node.RightChild;
                }

                return _nodes[nodeIdx].LeafValue;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public float PredictFromBinned(Span<byte> binned)
            {
                int nodeIdx = 0;

                while (!_nodes[nodeIdx].IsLeaf)
                {
                    ref readonly var node = ref _nodes[nodeIdx];
                    nodeIdx = binned[node.SplitFeature] <= node.SplitBin
                        ? node.LeftChild : node.RightChild;
                }

                return _nodes[nodeIdx].LeafValue;
            }
        }

        #endregion
    }
}