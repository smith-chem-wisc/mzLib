using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;

namespace Chromatography.RetentionTimePrediction;

/// <summary>
/// Contract for a retention time predictor. Implementations produce a predicted
/// retention-time-equivalent value (time, iRT, or hydrophobicity, depending on the
/// predictor) for a given peptide, individually or in batches.
/// </summary>
/// <remarks>
/// <para>
/// <b>Disposal:</b> the interface extends <see cref="IDisposable"/> so callers holding
/// an <see cref="IRetentionTimePredictor"/>-typed variable can write
/// <c>using var p = RetentionTimePredictorFactory.Create(...);</c> directly. Implementors
/// that hold unmanaged resources (e.g. the TorchSharp model in
/// <see cref="Chronologer.ChronologerRetentionTimePredictor"/>) release them in
/// <see cref="IDisposable.Dispose"/>; lightweight predictors rely on the base-class no-op.
/// </para>
/// <para>
/// <b>Construction is deliberately not part of this interface.</b> Use
/// <see cref="RetentionTimePredictorFactory.Create"/> with a
/// <see cref="PredictorType"/> value to build concrete predictors. This keeps the
/// abstraction free of references to concrete types — most importantly the
/// TorchSharp-backed <see cref="Chronologer.ChronologerRetentionTimePredictor"/>,
/// which would otherwise force every consumer of this interface to transitively
/// depend on TorchSharp and ship its large native binaries.
/// </para>
/// <para>
/// <b>Migration:</b> the previously-nested <c>Create</c> method and
/// <c>PredictorType</c> enum on this interface have been moved to the top-level
/// <see cref="RetentionTimePredictorFactory"/> and <see cref="PredictorType"/>.
/// </para>
/// </remarks>
public interface IRetentionTimePredictor : IDisposable
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
    /// Predicts a retention time equivalent for a given peptide.
    /// The value may represent a time, iRT, or hydrophobicity depending on the predictor.
    /// Returns null if prediction cannot be made.
    /// </summary>
    /// <returns>Predicted value in predictor-specific units, or null if prediction not possible</returns>
    double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason);

    /// <summary>
    /// Predicts retention time equivalents for a batch of peptides.
    /// Results are materialized and safe to enumerate multiple times.
    /// </summary>
    /// <remarks>
    /// <b>Order is not guaranteed.</b> Pair predictions to inputs via the
    /// <c>Peptide</c> element of each tuple, not by index — implementations may run
    /// in parallel and emit completed items in any order.
    /// </remarks>
    IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)> PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1);

    public string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason);
}