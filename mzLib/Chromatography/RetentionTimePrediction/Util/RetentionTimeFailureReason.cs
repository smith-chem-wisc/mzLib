namespace Chromatography.RetentionTimePrediction;

public enum RetentionTimeFailureReason
{
    EmptySequence,
    SequenceTooShort,
    SequenceTooLong,
    InvalidAminoAcid, // Most Commonly selenocysteine (U)
    InvalidMass, // CZE-Specific
    IncompatibleModifications, // Mod-related predictors only (e.g., Chronologer)
    PredictionError, // If the predictor encounters an error during prediction
}