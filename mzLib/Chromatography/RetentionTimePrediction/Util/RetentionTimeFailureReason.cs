namespace Chromatography.RetentionTimePrediction.Util;

public enum RetentionTimeFailureReason
{
    EmptySequence,
    SequenceTooShort,
    SequenceTooLong,
    InvalidAminoAcid,
    InvalidMass, // CZE-Specific
    IncompatibleModifications, // Mod-related predictors only (e.g., Chronologer)
}