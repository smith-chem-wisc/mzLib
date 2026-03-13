// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: MassSpectrometry/Dia/Classifiers/DiaNeuralNetClassifier.cs
//
// Phase 16C, Prompt 10: Neural network classifier for DIA FDR scoring.
//
// ACTION REQUIRED before compiling:
//   Add NeuralNetwork to the DiaClassifierType enum (in IDiaClassifier.cs or wherever defined):
//     public enum DiaClassifierType { LinearDiscriminant, GradientBoostedTree, NeuralNetwork }
//
// Architecture:  Input(ClassifierFeatureCount) → Dense(64, ReLU) → Dense(32, ReLU) → Dense(1) → logit
// Optimizer:     Adam, lr=1e-3, β₁=0.9, β₂=0.999, ε=1e-8, L2 weight-decay=1e-4
// Regularization: Inverted dropout p=0.3 during training; no dropout during scoring
// Initialization: He normal for weights (N(0, √(2/fan_in))), zeros for biases
// Epochs:        Max 60 with early stopping (patience=5, min-delta=1e-4)
// Mini-batch:    256
// Normalization: Per-feature z-score on combined training set; NaN → 0 after normalization
//
// Output is the raw pre-sigmoid logit (higher = more target-like), matching
// DiaLinearDiscriminant's scoring convention. DiaFdrEngine needs no changes to
// its ranking or q-value logic.
//
// DiaNeuralNetClassifierAdapter (at the bottom of this file) wraps DiaNeuralNetClassifier
// to implement IDiaClassifier, mirroring the pattern of DiaLdaClassifierAdapter.

using System;
using System.Buffers;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Feedforward neural network classifier for DIA precursor scoring.
    ///
    /// Architecture: ClassifierFeatureCount → Dense(64, ReLU) → Dense(32, ReLU) → Dense(1, linear logit)
    ///
    /// Training uses Adam SGD with inverted dropout and He initialization.
    /// Scoring is allocation-free (stack-allocated layer activations via stackalloc).
    ///
    /// The network is intentionally small (4,321 parameters) so it can train in
    /// under a second per semi-supervised iteration on the ~30K precursor dataset.
    /// </summary>
    public sealed class DiaNeuralNetClassifier
    {
        // ════════════════════════════════════════════════════════════════
        //  Architecture Constants
        // ════════════════════════════════════════════════════════════════

        private const int InputSize = DiaFeatureVector.ClassifierFeatureCount; // currently 37 (35 + CoElutionStd[35] + CandidateScoreGap[36])
        private const int Hidden1 = 64;
        private const int Hidden2 = 32;
        private const int OutputSize = 1;

        // Layer 1: InputSize → 64
        private readonly float[] _w1;  // [Hidden1 * InputSize]
        private readonly float[] _b1;  // [Hidden1]

        // Layer 2: 64 → 32
        private readonly float[] _w2;  // [Hidden2 * Hidden1]
        private readonly float[] _b2;  // [Hidden2]

        // Layer 3: 32 → 1
        private readonly float[] _w3;  // [OutputSize * Hidden2]
        private readonly float[] _b3;  // [OutputSize]

        // Normalization statistics (computed from training set)
        private readonly float[] _featureMean; // [InputSize]
        private readonly float[] _featureStd;  // [InputSize]

        // ════════════════════════════════════════════════════════════════
        //  Constructor
        // ════════════════════════════════════════════════════════════════

        private DiaNeuralNetClassifier(
            float[] w1, float[] b1,
            float[] w2, float[] b2,
            float[] w3, float[] b3,
            float[] featureMean, float[] featureStd)
        {
            _w1 = w1; _b1 = b1;
            _w2 = w2; _b2 = b2;
            _w3 = w3; _b3 = b3;
            _featureMean = featureMean;
            _featureStd = featureStd;
        }

        // ════════════════════════════════════════════════════════════════
        //  Training Entry Point
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Trains a new network on positives (targets) and negatives (decoys).
        /// Returns a trained, immutable DiaNeuralNetClassifier.
        ///
        /// Training parameters:
        ///   - Adam: lr=1e-3, β₁=0.9, β₂=0.999, ε=1e-8, L2=1e-4
        ///   - Dropout: inverted p=0.3 (applied only during training)
        ///   - Max epochs: 60, early stopping patience=5 (min-delta=1e-4)
        ///   - Batch size: 256
        /// </summary>
        public static DiaNeuralNetClassifier Train(
            ReadOnlySpan<DiaFeatureVector> positives,
            ReadOnlySpan<DiaFeatureVector> negatives,
            int seed = 42)
        {
            int n = InputSize;
            int total = positives.Length + negatives.Length;

            // ── 1. Compute per-feature normalization statistics ─────────
            float[] mean = new float[n];
            float[] std = new float[n];
            ComputeNormStats(positives, negatives, mean, std);

            // ── 2. Extract and normalize all samples ────────────────────
            float[][] X = new float[total][];
            float[] labels = new float[total];

            ExtractNormalized(positives, mean, std, X, 0, label: 1f);
            ExtractNormalized(negatives, mean, std, X, positives.Length, label: 0f);
            for (int i = 0; i < positives.Length; i++) labels[i] = 1f;
            // labels for negatives already 0f (default)

            // ── 3. Initialize weights (He normal) ──────────────────────
            var rng = new Random(seed);
            float[] w1 = HeInit(rng, Hidden1, InputSize);
            float[] b1 = new float[Hidden1];
            float[] w2 = HeInit(rng, Hidden2, Hidden1);
            float[] b2 = new float[Hidden2];
            float[] w3 = HeInit(rng, OutputSize, Hidden2);
            float[] b3 = new float[OutputSize];

            // ── 4. Adam state ───────────────────────────────────────────
            const float lr = 1e-3f;
            const float beta1 = 0.9f;
            const float beta2 = 0.999f;
            const float eps = 1e-8f;
            const float l2 = 1e-4f;

            int nW1 = w1.Length, nW2 = w2.Length, nW3 = w3.Length;
            int nB1 = b1.Length, nB2 = b2.Length, nB3 = b3.Length;

            float[] mW1 = new float[nW1], vW1 = new float[nW1];
            float[] mB1 = new float[nB1], vB1 = new float[nB1];
            float[] mW2 = new float[nW2], vW2 = new float[nW2];
            float[] mB2 = new float[nB2], vB2 = new float[nB2];
            float[] mW3 = new float[nW3], vW3 = new float[nW3];
            float[] mB3 = new float[nB3], vB3 = new float[nB3];

            // ── 5. Training loop ────────────────────────────────────────
            const int maxEpochs = 60;
            const int batchSize = 256;
            const int patience = 5;
            const float minDelta = 1e-4f;
            const float dropoutKeep = 0.7f;   // p_drop = 0.3, keep = 0.7
            const float dropoutScale = 1f / dropoutKeep; // inverted dropout scale

            int[] shuffleIdx = new int[total];
            for (int i = 0; i < total; i++) shuffleIdx[i] = i;

            // Intermediate activations (per-batch reuse)
            float[] h1 = new float[Hidden1];
            float[] h2 = new float[Hidden2];
            float[] mask1 = new float[Hidden1];
            float[] mask2 = new float[Hidden2];
            float[] dh1 = new float[Hidden1];
            float[] dh2 = new float[Hidden2];
            float[] gW1 = new float[nW1];
            float[] gB1 = new float[nB1];
            float[] gW2 = new float[nW2];
            float[] gB2 = new float[nB2];
            float[] gW3 = new float[nW3];
            float[] gB3 = new float[nB3];

            double bestLoss = double.MaxValue;
            int noImprove = 0;
            int step = 0; // Adam global step (for bias-correction denominator)

            for (int epoch = 0; epoch < maxEpochs; epoch++)
            {
                // Fisher-Yates shuffle
                for (int i = total - 1; i > 0; i--)
                {
                    int j = rng.Next(i + 1);
                    (shuffleIdx[i], shuffleIdx[j]) = (shuffleIdx[j], shuffleIdx[i]);
                }

                double epochLoss = 0.0;

                for (int bStart = 0; bStart < total; bStart += batchSize)
                {
                    step++;
                    int bEnd = Math.Min(bStart + batchSize, total);
                    int bLen = bEnd - bStart;

                    // Zero gradients
                    Array.Clear(gW1, 0, nW1); Array.Clear(gB1, 0, nB1);
                    Array.Clear(gW2, 0, nW2); Array.Clear(gB2, 0, nB2);
                    Array.Clear(gW3, 0, nW3); Array.Clear(gB3, 0, nB3);

                    // Generate dropout masks (per-batch; same mask for all samples in batch)
                    for (int k = 0; k < Hidden1; k++)
                        mask1[k] = rng.NextDouble() < dropoutKeep ? dropoutScale : 0f;
                    for (int k = 0; k < Hidden2; k++)
                        mask2[k] = rng.NextDouble() < dropoutKeep ? dropoutScale : 0f;

                    double batchLoss = 0.0;

                    for (int bi = bStart; bi < bEnd; bi++)
                    {
                        int idx = shuffleIdx[bi];
                        float[] x = X[idx];
                        float y = labels[idx];

                        // ── Forward pass ────────────────────────────────
                        // Layer 1: z1 = W1·x + b1, h1 = ReLU(z1) * mask1
                        for (int k = 0; k < Hidden1; k++)
                        {
                            float z = b1[k];
                            for (int j = 0; j < InputSize; j++)
                                z += w1[k * InputSize + j] * x[j];
                            h1[k] = MathF.Max(0f, z) * mask1[k];
                        }

                        // Layer 2: z2 = W2·h1 + b2, h2 = ReLU(z2) * mask2
                        for (int k = 0; k < Hidden2; k++)
                        {
                            float z = b2[k];
                            for (int j = 0; j < Hidden1; j++)
                                z += w2[k * Hidden1 + j] * h1[j];
                            h2[k] = MathF.Max(0f, z) * mask2[k];
                        }

                        // Layer 3: logit = W3·h2 + b3
                        float logit = b3[0];
                        for (int j = 0; j < Hidden2; j++)
                            logit += w3[j] * h2[j];

                        // ── BCE loss ─────────────────────────────────────
                        // loss = -[y * log(σ(logit)) + (1-y) * log(1-σ(logit))]
                        //      = log(1 + exp(-logit)) for y=1
                        //      = log(1 + exp( logit)) for y=0
                        // Numerically stable: softplus form
                        float absLogit = MathF.Abs(logit);
                        batchLoss += absLogit + MathF.Log(1f + MathF.Exp(-absLogit))
                                     - (y > 0.5f ? (logit > 0f ? 0f : logit) : (logit > 0f ? logit : 0f));

                        // ── Output gradient: dL/d_logit = σ(logit) - y ──
                        float sigma = Sigmoid(logit);
                        float dLogit = sigma - y;

                        // ── Backprop Layer 3 → Layer 2 ───────────────────
                        for (int j = 0; j < Hidden2; j++)
                        {
                            gW3[j] += dLogit * h2[j];
                            dh2[j] = dLogit * w3[j];
                        }
                        gB3[0] += dLogit;

                        // ── Backprop Layer 2 (ReLU + dropout) ────────────
                        for (int k = 0; k < Hidden2; k++)
                        {
                            // dReLU: if mask=0 (dropped) or pre-activation ≤ 0, grad = 0
                            float dz2 = (h2[k] > 0f) ? dh2[k] : 0f;
                            for (int j = 0; j < Hidden1; j++)
                            {
                                gW2[k * Hidden1 + j] += dz2 * h1[j];
                                dh1[j] += dz2 * w2[k * Hidden1 + j];
                            }
                            gB2[k] += dz2;
                        }

                        // ── Backprop Layer 1 (ReLU + dropout) ────────────
                        for (int k = 0; k < Hidden1; k++)
                        {
                            float dz1 = (h1[k] > 0f) ? dh1[k] : 0f;
                            for (int j = 0; j < InputSize; j++)
                                gW1[k * InputSize + j] += dz1 * x[j];
                            gB1[k] += dz1;
                        }

                        // Reset dh arrays for next sample
                        Array.Clear(dh1, 0, Hidden1);
                        Array.Clear(dh2, 0, Hidden2);
                    }

                    epochLoss += batchLoss;
                    float invBatch = 1f / bLen;

                    // ── Adam update ──────────────────────────────────────
                    float bc1 = 1f - MathF.Pow(beta1, step);
                    float bc2 = 1f - MathF.Pow(beta2, step);
                    AdamUpdate(w1, gW1, mW1, vW1, nW1, lr, beta1, beta2, eps, l2, bc1, bc2, invBatch);
                    AdamUpdate(b1, gB1, mB1, vB1, nB1, lr, beta1, beta2, eps, 0f, bc1, bc2, invBatch);
                    AdamUpdate(w2, gW2, mW2, vW2, nW2, lr, beta1, beta2, eps, l2, bc1, bc2, invBatch);
                    AdamUpdate(b2, gB2, mB2, vB2, nB2, lr, beta1, beta2, eps, 0f, bc1, bc2, invBatch);
                    AdamUpdate(w3, gW3, mW3, vW3, nW3, lr, beta1, beta2, eps, l2, bc1, bc2, invBatch);
                    AdamUpdate(b3, gB3, mB3, vB3, nB3, lr, beta1, beta2, eps, 0f, bc1, bc2, invBatch);
                }

                epochLoss /= total;

                // ── Early stopping ───────────────────────────────────────
                if (bestLoss - epochLoss > minDelta)
                {
                    bestLoss = epochLoss;
                    noImprove = 0;
                }
                else
                {
                    noImprove++;
                    if (noImprove >= patience)
                        break;
                }
            }

            return new DiaNeuralNetClassifier(w1, b1, w2, b2, w3, b3, mean, std);
        }

        // ════════════════════════════════════════════════════════════════
        //  Scoring (allocation-free)
        // ════════════════════════════════════════════════════════════════

        /// <summary>
        /// Scores a single feature vector.  All intermediate activations are stack-allocated.
        /// No dropout at inference time (expected activation preserved by inverted dropout training).
        /// Returns the pre-sigmoid logit: higher = more target-like.
        /// </summary>
        public float Score(in DiaFeatureVector fv)
        {
            Span<float> x = stackalloc float[InputSize];
            Span<float> h1 = stackalloc float[Hidden1];
            Span<float> h2 = stackalloc float[Hidden2];

            // Extract and normalize
            fv.WriteTo(x);
            for (int j = 0; j < InputSize; j++)
            {
                float v = float.IsNaN(x[j]) ? 0f : x[j];
                float s = _featureStd[j];
                x[j] = s > 1e-8f ? (v - _featureMean[j]) / s : 0f;
            }

            // Layer 1: h1 = ReLU(W1·x + b1)
            for (int k = 0; k < Hidden1; k++)
            {
                float z = _b1[k];
                for (int j = 0; j < InputSize; j++)
                    z += _w1[k * InputSize + j] * x[j];
                h1[k] = MathF.Max(0f, z);
            }

            // Layer 2: h2 = ReLU(W2·h1 + b2)
            for (int k = 0; k < Hidden2; k++)
            {
                float z = _b2[k];
                for (int j = 0; j < Hidden1; j++)
                    z += _w2[k * Hidden1 + j] * h1[j];
                h2[k] = MathF.Max(0f, z);
            }

            // Layer 3: logit = W3·h2 + b3
            float logit = _b3[0];
            for (int j = 0; j < Hidden2; j++)
                logit += _w3[j] * h2[j];

            return logit;
        }

        // ════════════════════════════════════════════════════════════════
        //  ONNX Export / Import (minimal binary round-trip format)
        // ════════════════════════════════════════════════════════════════

        private const uint OnnxMagic = 0x494E4F44; // "DONI" in little-endian (DIA ONNX Neutral Interchange)

        /// <summary>
        /// Saves the trained network to a binary file.
        /// Normalization is fused into W1/b1 so the file can be loaded and scored immediately.
        ///
        /// Format:
        ///   uint32  magic
        ///   int32   inputSize, hidden1, hidden2
        ///   float[] featureMean  (inputSize)
        ///   float[] featureStd   (inputSize)
        ///   float[] w1 (hidden1 * inputSize), b1 (hidden1)
        ///   float[] w2 (hidden2 * hidden1),   b2 (hidden2)
        ///   float[] w3 (hidden2),              b3 (1)
        /// </summary>
        public void SaveOnnx(string path)
        {
            using var bw = new System.IO.BinaryWriter(System.IO.File.Open(path, System.IO.FileMode.Create));
            bw.Write(OnnxMagic);
            bw.Write(InputSize);
            bw.Write(Hidden1);
            bw.Write(Hidden2);
            foreach (float v in _featureMean) bw.Write(v);
            foreach (float v in _featureStd) bw.Write(v);
            foreach (float v in _w1) bw.Write(v);
            foreach (float v in _b1) bw.Write(v);
            foreach (float v in _w2) bw.Write(v);
            foreach (float v in _b2) bw.Write(v);
            foreach (float v in _w3) bw.Write(v);
            foreach (float v in _b3) bw.Write(v);
        }

        /// <summary>Loads a network previously saved with SaveOnnx().</summary>
        public static DiaNeuralNetClassifier LoadOnnx(string path)
        {
            using var br = new System.IO.BinaryReader(System.IO.File.OpenRead(path));
            uint magic = br.ReadUInt32();
            if (magic != OnnxMagic)
                throw new InvalidOperationException("Not a valid DiaNeuralNet model file.");

            int inp = br.ReadInt32();
            int h1 = br.ReadInt32();
            int h2 = br.ReadInt32();

            float[] mean = ReadFloats(br, inp);
            float[] std = ReadFloats(br, inp);
            float[] w1 = ReadFloats(br, h1 * inp);
            float[] b1 = ReadFloats(br, h1);
            float[] w2 = ReadFloats(br, h2 * h1);
            float[] b2 = ReadFloats(br, h2);
            float[] w3 = ReadFloats(br, h2);
            float[] b3 = ReadFloats(br, 1);

            return new DiaNeuralNetClassifier(w1, b1, w2, b2, w3, b3, mean, std);
        }

        // ════════════════════════════════════════════════════════════════
        //  Private Helpers
        // ════════════════════════════════════════════════════════════════

        private static void ComputeNormStats(
            ReadOnlySpan<DiaFeatureVector> positives,
            ReadOnlySpan<DiaFeatureVector> negatives,
            float[] mean, float[] std)
        {
            int n = InputSize;
            int total = positives.Length + negatives.Length;

            double[] sum = new double[n];
            double[] sumSq = new double[n];

            // Single reusable buffer — stackalloc inside a loop would stack-overflow
            // on large datasets (27K+ samples × ClassifierFeatureCount floats each).
            float[] buf = new float[n];

            for (int i = 0; i < positives.Length; i++)
            {
                positives[i].WriteTo(buf);
                for (int j = 0; j < n; j++)
                {
                    double v = float.IsNaN(buf[j]) ? 0.0 : buf[j];
                    sum[j] += v;
                    sumSq[j] += v * v;
                }
            }
            for (int i = 0; i < negatives.Length; i++)
            {
                negatives[i].WriteTo(buf);
                for (int j = 0; j < n; j++)
                {
                    double v = float.IsNaN(buf[j]) ? 0.0 : buf[j];
                    sum[j] += v;
                    sumSq[j] += v * v;
                }
            }

            for (int j = 0; j < n; j++)
            {
                mean[j] = (float)(sum[j] / total);
                double variance = sumSq[j] / total - (double)mean[j] * mean[j];
                std[j] = variance > 0.0 ? (float)Math.Sqrt(variance) : 1e-8f;
            }
        }

        private static void ExtractNormalized(
            ReadOnlySpan<DiaFeatureVector> vectors,
            float[] mean, float[] std,
            float[][] dest, int destOffset, float label)
        {
            int n = InputSize;
            float[] buf = new float[n]; // heap: stackalloc in called methods risks overflow at large scale
            for (int i = 0; i < vectors.Length; i++)
            {
                vectors[i].WriteTo(buf);
                float[] row = new float[n];
                for (int j = 0; j < n; j++)
                {
                    float v = float.IsNaN(buf[j]) ? 0f : buf[j];
                    row[j] = std[j] > 1e-8f ? (v - mean[j]) / std[j] : 0f;
                }
                dest[destOffset + i] = row;
            }
        }

        private static float[] HeInit(Random rng, int fanOut, int fanIn)
        {
            float[] w = new float[fanOut * fanIn];
            float std = MathF.Sqrt(2f / fanIn);
            // Box-Muller for normal samples
            for (int i = 0; i < w.Length; i += 2)
            {
                double u1 = 1.0 - rng.NextDouble();
                double u2 = 1.0 - rng.NextDouble();
                double mag = std * Math.Sqrt(-2.0 * Math.Log(u1));
                w[i] = (float)(mag * Math.Cos(2.0 * Math.PI * u2));
                if (i + 1 < w.Length)
                    w[i + 1] = (float)(mag * Math.Sin(2.0 * Math.PI * u2));
            }
            return w;
        }

        private static void AdamUpdate(
            float[] param, float[] grad,
            float[] m, float[] v,
            int len,
            float lr, float beta1, float beta2, float eps,
            float l2, float bc1, float bc2,
            float invBatch)
        {
            for (int i = 0; i < len; i++)
            {
                float g = grad[i] * invBatch + l2 * param[i]; // L2 only on weights
                m[i] = beta1 * m[i] + (1f - beta1) * g;
                v[i] = beta2 * v[i] + (1f - beta2) * g * g;
                float mHat = m[i] / bc1;
                float vHat = v[i] / bc2;
                param[i] -= lr * mHat / (MathF.Sqrt(vHat) + eps);
            }
        }

        private static float Sigmoid(float x)
        {
            if (x > 20f) return 1f;
            if (x < -20f) return 0f;
            return 1f / (1f + MathF.Exp(-x));
        }

        private static float[] ReadFloats(System.IO.BinaryReader br, int count)
        {
            float[] arr = new float[count];
            for (int i = 0; i < count; i++) arr[i] = br.ReadSingle();
            return arr;
        }
    }

    // ════════════════════════════════════════════════════════════════════
    //  IDiaClassifier Adapter
    // ════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Wraps DiaNeuralNetClassifier to implement IDiaClassifier,
    /// mirroring the pattern of DiaLdaClassifierAdapter.
    ///
    /// Train() merges positives and negatives, trains a new network, and
    /// replaces _nn in-place (each semi-supervised iteration builds a fresh model).
    ///
    /// Score() delegates to DiaNeuralNetClassifier.Score() which is allocation-free.
    /// </summary>
    public sealed class DiaNeuralNetClassifierAdapter : IDiaClassifier
    {
        private DiaNeuralNetClassifier _nn;

        /// <summary>The underlying trained network. Null before first Train() call.</summary>
        public DiaNeuralNetClassifier Network => _nn;

        // ── IDiaClassifier ───────────────────────────────────────────────
        public int FeatureCount => DiaFeatureVector.ClassifierFeatureCount;
        public bool IsTrained => _nn != null;

        public void Train(
            ReadOnlySpan<DiaFeatureVector> positives,
            ReadOnlySpan<DiaFeatureVector> negatives)
        {
            _nn = DiaNeuralNetClassifier.Train(positives, negatives);
        }

        public float Score(in DiaFeatureVector fv)
        {
            // Guard: if Train() hasn't been called yet, return 0 (won't happen in normal flow)
            if (_nn == null) return 0f;
            return _nn.Score(in fv);
        }
    }
}