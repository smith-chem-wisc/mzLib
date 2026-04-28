using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;

namespace Chromatography.RetentionTimePrediction;

/// <summary>
/// Unified contract for all retention time predictors — local (Chronologer, SSRCalc3)
/// and remote (Koina/Prosit). MetaMorpheus consumes this interface directly without
/// an adapter layer.
/// </summary>
public interface IRetentionTimePredictor
{
    /// <summary>Human-readable name, e.g. "Chronologer", "Prosit2019iRT".</summary>
    string PredictorName { get; }

    /// <summary>Chromatographic or electrophoretic mode this predictor targets.</summary>
    SeparationType SeparationType { get; }

    /// <summary>
    /// Predicts retention time for a single peptide.
    /// Returns null when prediction is not possible (invalid sequence, unsupported
    /// modifications, model error, etc.). The <paramref name="failureReason"/> out
    /// parameter is set when null is returned.
    /// </summary>
    double? PredictRetentionTime(IRetentionPredictable peptide,
                                 out RetentionTimeFailureReason? failureReason);

    /// <summary>
    /// Predicts retention times for a batch of peptides.
    ///
    /// The default implementation calls <see cref="PredictRetentionTime"/> per peptide.
    /// Koina-backed implementations override this to issue a single HTTP call for the
    /// entire batch.
    ///
    /// The returned dictionary is read-only. Key is <c>peptide.FullSequence</c>.
    /// A null value indicates prediction was not possible for that peptide; the specific
    /// failure reason is not preserved in the batch result. Callers requiring per-peptide
    /// failure diagnostics should call <see cref="PredictRetentionTime"/> individually.
    ///
    /// Callers are responsible for deduplicating input before calling this method.
    /// Duplicate <c>FullSequence</c> entries result in redundant prediction work in the
    /// default implementation (the last result silently overwrites earlier ones in the
    /// output dictionary).
    ///
    /// Null elements in <paramref name="peptides"/> are skipped.
    /// </summary>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="peptides"/> is null.
    /// </exception>
    IReadOnlyDictionary<string, double?> PredictRetentionTimes(
        IEnumerable<IRetentionPredictable> peptides)
    {
        if (peptides is null)
            throw new ArgumentNullException(nameof(peptides));

        var results = new Dictionary<string, double?>();
        foreach (var peptide in peptides)
        {
            if (peptide is null)
                continue;
            results[peptide.FullSequence] = PredictRetentionTime(peptide, out _);
        }

        return new ReadOnlyDictionary<string, double?>(results);
    }

    /// <summary>
    /// Returns the predictor-specific formatted sequence string for a peptide,
    /// or null (with <paramref name="failureReason"/> set) if the peptide cannot
    /// be formatted for this predictor.
    /// Useful for diagnostics and for callers that cache formatted sequences.
    /// </summary>
    string? GetFormattedSequence(IRetentionPredictable peptide,
                                  out RetentionTimeFailureReason? failureReason);
}
