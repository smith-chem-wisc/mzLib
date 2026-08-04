using System;
using System.Collections.Generic;
using Omics.RetentionTimePrediction;

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
/// <c>ChronologerRetentionTimePredictor</c>) release them in
/// <see cref="IDisposable.Dispose"/>; lightweight predictors rely on the base-class no-op.
/// </para>
/// <para>
/// <b>Construction is deliberately not part of this interface.</b> Use
/// <c>RetentionTimePredictorFactory.Create</c> with a <c>PredictorType</c> value to build
/// concrete predictors. This keeps the abstraction free of references to concrete types —
/// most importantly the TorchSharp-backed <c>ChronologerRetentionTimePredictor</c>, which
/// would otherwise force every consumer of this interface to transitively depend on
/// TorchSharp and ship its large native binaries.
/// </para>
/// <para>
/// <b>Assembly layout:</b> the factory, the <c>PredictorType</c> enum and the Chronologer
/// predictor all live in the separate <c>Chromatography.Chronologer</c> assembly, which is
/// why they are referenced here as plain code spans rather than <c>see cref</c> links. This
/// assembly deliberately does not reference that one; the dependency runs the other way.
/// </para>
/// </remarks>
public interface IRetentionTimePredictor : IDisposable
{
    /// <summary>Human-readable name, e.g. "Chronologer", "Prosit2019iRT".</summary>
    string PredictorName { get; }

    /// <summary>
    /// Gets the separation type this predictor is designed for
    /// </summary>
    SeparationType SeparationType { get; }

    /// <summary>
    /// Returns the predictor-specific formatted sequence string for a peptide,
    /// or null (with <paramref name="failureReason"/> set) if the peptide cannot
    /// be formatted for this predictor.
    /// Useful for diagnostics and for callers that cache formatted sequences.
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
