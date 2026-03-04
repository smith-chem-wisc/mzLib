// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: MassSpectrometry/Dia/FDR/IDiaClassifier.cs

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Abstraction for DIA classifiers used in the iterative FDR workflow.
    /// Both DiaLinearDiscriminant (LDA) and DiaGradientBoostedClassifier (GBT)
    /// implement this via adapter classes.
    /// 
    /// DiaFdrEngine depends only on this interface, making classifier selection
    /// a runtime decision with no code duplication in the FDR loop.
    /// </summary>
    public interface IDiaClassifier
    {
        /// <summary>Number of features expected per sample</summary>
        int FeatureCount { get; }

        /// <summary>Whether the model has been trained</summary>
        bool IsTrained { get; }

        /// <summary>
        /// Trains on paired positive/negative feature vectors.
        /// This signature matches the existing DiaLinearDiscriminant.TrainLogisticRegression pattern.
        /// </summary>
        void Train(ReadOnlySpan<DiaFeatureVector> positives, ReadOnlySpan<DiaFeatureVector> negatives);

        /// <summary>
        /// Scores a single feature vector. Higher = more target-like.
        /// For GBT this is P(target | features). For LDA this is the logistic score.
        /// </summary>
        float Score(in DiaFeatureVector features);
    }

    /// <summary>
    /// Specifies which classifier to use in the DIA FDR workflow.
    /// </summary>
    public enum DiaClassifierType
    {
        /// <summary>Linear Discriminant Analysis (Phase 10 original — fast, linear boundaries only)</summary>
        LinearDiscriminant,

        /// <summary>Gradient-Boosted Trees (Phase 14 — non-linear, captures feature interactions)</summary>
        GradientBoostedTree,
        
        /// <summary>
        /// Represents a neural network model used for machine learning tasks.
        /// </summary>
        /// <remarks>PHase 16c, Prompt 10. This class provides methods for training the neural network, making predictions, and
        /// evaluating performance. It supports various architectures and can be configured with different
        /// hyperparameters to optimize learning.</remarks>
        NeuralNetwork
    }
}