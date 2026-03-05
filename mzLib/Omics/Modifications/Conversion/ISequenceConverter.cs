using System.Collections.Generic;

namespace Omics.Modifications.Conversion;

public interface ISequenceConverter
{
    bool TryConvertFullSequence(
        IBioPolymerWithSetMods bioPolymer,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out string? converted,
        out SequenceConversionFailureReason? failureReason);

    string ConvertFullSequence(
        IBioPolymerWithSetMods bioPolymer,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException);

    bool TryConvertFullSequence(
        string fullSequence,
        ModificationNamingConvention sourceConvention,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out string? converted,
        out SequenceConversionFailureReason? failureReason);

    string ConvertFullSequence(
        string fullSequence,
        ModificationNamingConvention sourceConvention,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException);

    Dictionary<int, Modification> ConvertModifications(
        IReadOnlyDictionary<int, Modification> source,
        string baseSequence,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException);

    bool TryConvertModifications(
        IReadOnlyDictionary<int, Modification> source,
        string baseSequence,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out Dictionary<int, Modification>? converted,
        out SequenceConversionFailureReason? failureReason);

    bool TryConvertModification(
        Modification source,
        string baseSequence,
        int modificationKey,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode,
        out Modification? converted,
        out SequenceConversionFailureReason? reason);

    Modification ConvertModification(
        Modification source,
        string baseSequence,
        int modificationKey,
        ModificationNamingConvention targetConvention,
        SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException);

    bool TryConvertModificationDefinition(
        Modification source,
        ModificationNamingConvention targetConvention,
        out Modification? converted);

    Modification ConvertModificationDefinition(
        Modification source,
        ModificationNamingConvention targetConvention);

    bool TryGetChronologerSequence(
        IBioPolymerWithSetMods bioPolymer,
        SequenceConversionHandlingMode handlingMode,
        out string? massShiftSequence,
        out SequenceConversionFailureReason? reason);
}
