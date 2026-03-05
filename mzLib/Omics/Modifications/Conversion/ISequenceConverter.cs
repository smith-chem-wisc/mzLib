using System.Collections.Generic;

namespace Omics.Modifications.Conversion;

/// <summary>
/// Defines conversion operations for sequences and modifications.
/// </summary>
public interface ISequenceConverter
{
    /// <summary>
    /// Attempts to convert a bio-polymer sequence to the target convention.
    /// </summary>
    bool TryConvertFullSequence(
        IBioPolymerWithSetMods bioPolymer,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out string? converted,
        out SequenceConversionFailureReason? failureReason);

    /// <summary>
    /// Converts a bio-polymer sequence to the target convention.
    /// </summary>
    string ConvertFullSequence(
        IBioPolymerWithSetMods bioPolymer,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException);

    /// <summary>
    /// Attempts to convert a full sequence string between conventions.
    /// </summary>
    bool TryConvertFullSequence(
        string fullSequence,
        ModificationNamingConvention sourceConvention,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out string? converted,
        out SequenceConversionFailureReason? failureReason);

    /// <summary>
    /// Converts a full sequence string between conventions.
    /// </summary>
    string ConvertFullSequence(
        string fullSequence,
        ModificationNamingConvention sourceConvention,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException);

    /// <summary>
    /// Converts a modification dictionary to the target convention.
    /// </summary>
    Dictionary<int, Modification> ConvertModifications(
        IReadOnlyDictionary<int, Modification> source,
        string baseSequence,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException);

    /// <summary>
    /// Attempts to convert a modification dictionary to the target convention.
    /// </summary>
    bool TryConvertModifications(
        IReadOnlyDictionary<int, Modification> source,
        string baseSequence,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out Dictionary<int, Modification>? converted,
        out SequenceConversionFailureReason? failureReason);

    /// <summary>
    /// Attempts to convert a single modification at a specific position.
    /// </summary>
    bool TryConvertModification(
        Modification source,
        string baseSequence,
        int modificationKey,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out Modification? converted,
        out SequenceConversionFailureReason? reason);

    /// <summary>
    /// Converts a single modification at a specific position.
    /// </summary>
    Modification ConvertModification(
        Modification source,
        string baseSequence,
        int modificationKey,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException);

    /// <summary>
    /// Attempts to convert a modification definition to the target convention.
    /// </summary>
    bool TryConvertModificationDefinition(
        Modification source,
        ModificationNamingConvention targetConvention,
        out Modification? converted);

    /// <summary>
    /// Converts a modification definition to the target convention.
    /// </summary>
    Modification ConvertModificationDefinition(
        Modification source,
        ModificationNamingConvention targetConvention);

    /// <summary>
    /// Attempts to build a Chronologer-compatible sequence.
    /// </summary>
    bool TryGetChronologerSequence(
        IBioPolymerWithSetMods bioPolymer,
        SequenceConversionHandlingMode handlingMode,
        out string? massShiftSequence,
        out SequenceConversionFailureReason? reason);
}
