using System.Collections.Generic;

namespace Chromatography.RetentionTimePrediction;

/// <summary>
/// Unified contract for all retention time predictors — local (Chronologer, SSRCalc3)
/// and remote (Koina/Prosit). MetaMorpheus consumes this interface directly without
/// an adapter layer.
/// </summary>
public interface IRetentionTimePredictor
{
    /// <summary>Human-readable name, e.g. "Chronologer", "Prosit2019iRT"</summary>
    string PredictorName { get; }

    /// <summary>Chromatographic or electrophoretic mode this predictor targets.</summary>
    SeparationType SeparationType { get; }

    /// <summary>
    /// Predicts retention time for a single peptide.
    /// Returns null when prediction is not possible (invalid sequence, unsupported
    /// modifications, model error, etc.).
    /// </summary>
    double? PredictRetentionTime(IRetentionPredictable peptide,
                                 out RetentionTimeFailureReason? failureReason);

    /// <summary>
    /// Predicts retention times for a batch of peptides.
    ///
    /// The default implementation calls PredictRetentionTime per peptide.
    /// Koina-backed implementations should override this to issue a single HTTP
    /// call for the entire batch.
    ///
    /// Dictionary key is peptide.FullSequence.
    /// Null values indicate prediction was not possible for that peptide.
    /// </summary>
    Dictionary<string, double?> PredictRetentionTimes(
        IEnumerable<IRetentionPredictable> peptides)
    {
        var results = new Dictionary<string, double?>();
        foreach (var peptide in peptides)
            results[peptide.FullSequence] = PredictRetentionTime(peptide, out _);
        return results;
    }

    /// <summary>
    /// Returns the predictor-specific formatted sequence string for a peptide,
    /// or null if the peptide cannot be formatted for this predictor.
    /// Useful for diagnostics and for callers that cache formatted sequences.
    /// </summary>
    string? GetFormattedSequence(IRetentionPredictable peptide,
                                  out RetentionTimeFailureReason? failureReason);
}
