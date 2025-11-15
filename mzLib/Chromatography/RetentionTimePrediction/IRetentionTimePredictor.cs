using Chromatography.RetentionTimePrediction.Util;

namespace Chromatography.RetentionTimePrediction
{
    /// <summary>
    /// Represents a retention time predictor that can predict RT for peptides.
    /// Designed for high-performance batch prediction with optional modification checking.
    /// </summary>
    public interface IRetentionTimePredictor
    {
        /// <summary>
        /// Name/identifier for this predictor (e.g., "SSRCalc3", "Chronologer")
        /// </summary>
        string PredictorName { get; }

        /// <summary>
        /// Indicates whether this predictor requires modification compatibility checking.
        /// If true, a IModificationCompatibilityChecker should be provided for full functionality.
        /// </summary>
        bool RequiresModificationChecking { get; }

        /// <summary>
        /// Predicts retention time for a given peptide.
        /// Returns null if prediction cannot be made.
        /// </summary>
        /// <param name="peptide">The peptide for which to predict RT</param>
        /// <returns>Predicted retention time in predictor-specific units, or null if prediction not possible</returns>
        double? PredictRetentionTime(IRetentionPredictable peptide);

        /// <summary>
        /// Checks if this predictor can handle the given peptide based on basic constraints
        /// (length, amino acid content, etc.) - does NOT check modifications.
        /// </summary>
        /// <param name="peptide">The peptide to validate</param>
        /// <param name="failureReason">Brief reason if invalid</param>
        /// <returns>True if prediction is possible based on basic constraints</returns>
        bool CanPredict(IRetentionPredictable peptide, out string? failureReason);

        /// <summary>
        /// Sets the modification compatibility checker for this predictor.
        /// Only needed if RequiresModificationChecking is true.
        /// </summary>
        /// <param name="checker">The compatibility checker to use</param>
        void SetModificationChecker(IModificationCompatibilityChecker checker);
    }
}