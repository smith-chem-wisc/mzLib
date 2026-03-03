// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: MassSpectrometry/Dia/FDR/DiaGbtClassifierAdapter.cs

using System;
using System.Buffers;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Wraps <see cref="DiaGradientBoostedClassifier"/> behind the <see cref="IDiaClassifier"/>
    /// interface so it can be used interchangeably with DiaLinearDiscriminant in DiaFdrEngine.
    /// 
    /// Handles the conversion between DiaFeatureVector (struct with named fields) and
    /// the flat float[] layout that the GBT expects.
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

        /// <summary>
        /// Trains the GBT on positive (target) and negative (decoy) feature vectors.
        /// Converts DiaFeatureVector spans → flat float arrays for the GBT Train() method.
        /// </summary>
        public void Train(ReadOnlySpan<DiaFeatureVector> positives, ReadOnlySpan<DiaFeatureVector> negatives)
        {
            int nPos = positives.Length;
            int nNeg = negatives.Length;
            int n = nPos + nNeg;
            int fc = DiaFeatureVector.ClassifierFeatureCount;

            var featureMatrix = ArrayPool<float>.Shared.Rent(n * fc);
            var labels = ArrayPool<float>.Shared.Rent(n);
            try
            {
                // Pack positives (label = 1)
                for (int i = 0; i < nPos; i++)
                {
                    positives[i].WriteTo(featureMatrix.AsSpan(i * fc, fc));
                    labels[i] = 1f;
                }
                // Pack negatives (label = 0)
                for (int i = 0; i < nNeg; i++)
                {
                    negatives[i].WriteTo(featureMatrix.AsSpan((nPos + i) * fc, fc));
                    labels[nPos + i] = 0f;
                }

                _gbt.Train(
                    featureMatrix.AsSpan(0, n * fc),
                    labels.AsSpan(0, n),
                    n);
            }
            finally
            {
                ArrayPool<float>.Shared.Return(featureMatrix);
                ArrayPool<float>.Shared.Return(labels);
            }
        }

        /// <summary>
        /// Scores a single feature vector by converting to flat float span,
        /// then calling GBT PredictProbability.
        /// </summary>
        public float Score(in DiaFeatureVector features)
        {
            Span<float> flat = stackalloc float[DiaFeatureVector.ClassifierFeatureCount];
            features.WriteTo(flat);
            return _gbt.PredictProbability(flat);
        }
    }

    /// <summary>
    /// Wraps the existing <see cref="DiaLinearDiscriminant"/> behind <see cref="IDiaClassifier"/>.
    /// This is a thin pass-through since DiaLinearDiscriminant already uses DiaFeatureVector directly.
    /// </summary>
    internal sealed class DiaLdaClassifierAdapter : IDiaClassifier
    {
        private DiaLinearDiscriminant _lda;
        private readonly float _learningRate;
        private readonly float _l2Lambda;
        private readonly int _maxEpochs;
        private readonly int _batchSize;

        public DiaLdaClassifierAdapter(
            float learningRate = 0.05f,
            float l2Lambda = 5e-3f,
            int maxEpochs = 300,
            int batchSize = 256)
        {
            _learningRate = learningRate;
            _l2Lambda = l2Lambda;
            _maxEpochs = maxEpochs;
            _batchSize = batchSize;
        }

        public int FeatureCount => DiaFeatureVector.ClassifierFeatureCount;
        public bool IsTrained => _lda != null;

        public void Train(ReadOnlySpan<DiaFeatureVector> positives, ReadOnlySpan<DiaFeatureVector> negatives)
        {
            _lda = DiaLinearDiscriminant.TrainLogisticRegression(
                positives, negatives,
                learningRate: _learningRate,
                l2Lambda: _l2Lambda,
                maxEpochs: _maxEpochs,
                batchSize: _batchSize,
                useInteraction: false);
        }

        public float Score(in DiaFeatureVector features) => _lda.Score(in features);

        /// <summary>Exposes the underlying LDA for diagnostics (weights, bias, etc.)</summary>
        public DiaLinearDiscriminant Lda => _lda;
    }
}