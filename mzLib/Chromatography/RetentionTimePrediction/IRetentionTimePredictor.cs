using Chromatography.RetentionTimePrediction.CZE;
using Chromatography.RetentionTimePrediction.Chronologer;
using Chromatography.RetentionTimePrediction.SSRCalc;

namespace Chromatography.RetentionTimePrediction;

/// <summary>
/// Represents a retention time predictor that can predict RT for peptides.
/// </summary>
public interface IRetentionTimePredictor : IDisposable
{
    /// <summary>
    /// The available predictor types that can be created via <see cref="Create"/>.
    /// </summary>
    public enum PredictorType
    {
        SSRCalc3,
        CZE,
        Chronologer
    }

    /// <summary>
    /// Creates a new predictor of the specified type with default parameters.
    /// The caller is responsible for disposing the returned predictor.
    /// </summary>
    public static IRetentionTimePredictor Create(PredictorType type) => type switch
    {
        PredictorType.SSRCalc3 => new SSRCalc3RetentionTimePredictor(),
        PredictorType.CZE => new CZERetentionTimePredictor(),
        PredictorType.Chronologer => new ChronologerRetentionTimePredictor(),
        _ => throw new ArgumentOutOfRangeException(nameof(type), type, null)
    };

    /// <summary>
    /// Name/identifier for this predictor (e.g., "SSRCalc3", "Chronologer")
    /// </summary>
    string PredictorName { get; }

    /// <summary>
    /// Gets the separation type this predictor is designed for
    /// </summary>
    SeparationType SeparationType { get; }

    /// <summary>
    /// Predicts a retention time equivalent for a given peptide.
    /// The value may represent a time, iRT, or hydrophobicity depending on the predictor.
    /// Returns null if prediction cannot be made.
    /// </summary>
    /// <returns>Predicted value in predictor-specific units, or null if prediction not possible</returns>
    double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason);

    /// <summary>
    /// Predicts retention time equivalents for a batch of peptides.
    /// </summary>
    IReadOnlyList<RetentionTimeEquivalentResult> PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides);

    public string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason);
}