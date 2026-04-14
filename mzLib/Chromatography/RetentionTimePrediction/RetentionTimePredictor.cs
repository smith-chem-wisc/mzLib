using Omics.SequenceConversion;

namespace Chromatography.RetentionTimePrediction;

public abstract class RetentionTimePredictor : IRetentionTimePredictor
{
    protected static readonly char[] CanonicalAminoAcids = "ACDEFGHIKLMNPQRSTVWY".ToCharArray();

    /// <summary>
    /// Name/identifier for this predictor (e.g., "SSRCalc3", "Chronologer")
    /// </summary>
    public abstract string PredictorName { get; }

    /// <summary>
    /// Gets the separation type this predictor is designed for
    /// </summary>
    public abstract SeparationType SeparationType { get; }

    /// <summary>
    /// Determine what to do if we encounter incompatible modifications in an RT predictor where mods are important. 
    /// </summary>
    public virtual SequenceConversionHandlingMode SequenceHandlingMode { get; }
    protected virtual int MinSequenceLength => 7;
    protected virtual int MaxSequenceLength => int.MaxValue;

    protected RetentionTimePredictor(SequenceConversionHandlingMode sequenceHandlingMode = SequenceConversionHandlingMode.UsePrimarySequence)
    {
        SequenceHandlingMode = sequenceHandlingMode;
    }

    public double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        if (peptide == null)
            throw new ArgumentNullException(nameof(peptide));

        if (!ValidateBasicConstraints(peptide, out failureReason))
            return null;

        string? sequenceToPredict = GetFormattedSequence(peptide, out failureReason);
        if (sequenceToPredict == null)
            return null;
        try
        {

            return PredictCore(peptide, sequenceToPredict);
        }
        catch (Exception)
        {
            failureReason = RetentionTimeFailureReason.PredictionError;
            return null;
        }
    }

    /// <summary>
    /// Core prediction logic - called when peptide passes all validation
    /// </summary>
    protected abstract double? PredictCore(IRetentionPredictable peptide, string? formattedSequence = null);

    /// <summary>
    /// Default batch implementation — loops over singles. Predictors that support true
    /// batched inference (e.g. Chronologer) should override this method.
    /// </summary>
    // TODO: Chronologer should override this with a true batched tensor call
    public virtual IReadOnlyList<RetentionTimeEquivalentResult> PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides)
    {
        return peptides
            .Select(p => new RetentionTimeEquivalentResult(p, PredictRetentionTimeEquivalent(p, out var reason), reason))
            .ToList();
    }

    /// <summary>
    /// No-op dispose for predictors that hold no unmanaged resources.
    /// Chronologer overrides this to release its TorchSharp model.
    /// </summary>
    public virtual void Dispose() { }

    /// <summary>
    /// Format the peptide sequence appropriately for this predictor
    /// </summary>
    public abstract string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason);

    /// <summary>
    /// Validate basic constraints (length, amino acid content, etc.) and NOT modifications
    /// </summary>
    protected virtual bool ValidateBasicConstraints(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        failureReason = null;

        if (string.IsNullOrEmpty(peptide.BaseSequence))
        {
            failureReason = RetentionTimeFailureReason.EmptySequence;
            return false;
        }

        if (peptide.BaseSequence.Length < MinSequenceLength)
        {
            failureReason = RetentionTimeFailureReason.SequenceTooShort;
            return false;
        }

        if (peptide.BaseSequence.Length > MaxSequenceLength)
        {
            failureReason = RetentionTimeFailureReason.SequenceTooLong;
            return false;
        }

        return true;
    }
}
/// <summary>
/// The result of a retention time equivalent prediction for a single peptide.
/// The predicted value may represent a time, iRT, or hydrophobicity depending on the predictor used.
/// </summary>
public readonly record struct RetentionTimeEquivalentResult(
    IRetentionPredictable Peptide,
    double? PredictedValue,
    RetentionTimeFailureReason? FailureReason)
{
    /// <summary>
    /// True if a value was successfully predicted.
    /// </summary>
    public bool Success => PredictedValue.HasValue;
}