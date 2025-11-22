using Chromatography.RetentionTimePrediction.Util;

namespace Chromatography.RetentionTimePrediction;

public abstract class RetentionTimePredictor : IRetentionTimePredictor
{
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
    public virtual IncompatibleModHandlingMode ModHandlingMode { get; }
    protected virtual int MinSequenceLength => 7;
    protected virtual int MaxSequenceLength => int.MaxValue;

    protected RetentionTimePredictor(IncompatibleModHandlingMode modHandlingMode = IncompatibleModHandlingMode.UsePrimarySequence)
    {
        ModHandlingMode = modHandlingMode;
    }

    public double? PredictRetentionTime(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        if (!ValidateBasicConstraints(peptide, out failureReason))
            return null;

        string? sequenceToPredict = GetFormattedSequence(peptide, out failureReason);
        if (sequenceToPredict == null)
            return null;

        return PredictCore(peptide, sequenceToPredict);
    }

    /// <summary>
    /// Core prediction logic - called when peptide passes all validation
    /// </summary>
    protected abstract double? PredictCore(IRetentionPredictable peptide, string? formattedSequence = null);

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