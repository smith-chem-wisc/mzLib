using Chromatography.RetentionTimePrediction.Util;
namespace Chromatography.RetentionTimePrediction;

/// <summary>
/// Represents a retention time predictor that can predict RT for peptides.
/// </summary>
public interface IRetentionTimePredictor
{
    /// <summary>
    /// Name/identifier for this predictor (e.g., "SSRCalc3", "Chronologer")
    /// </summary>
    string PredictorName { get; }

    /// <summary>
    /// Gets the separation type this predictor is designed for
    /// </summary>
    SeparationType SeparationType { get; }

    /// <summary>
    /// Predicts retention time for a given peptide.
    /// Returns null if prediction cannot be made.
    /// </summary>
    /// <returns>Predicted retention time in predictor-specific units, or null if prediction not possible</returns>
    double? PredictRetentionTime(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason);

    public string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason);
}