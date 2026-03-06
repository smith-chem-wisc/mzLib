namespace Omics.SequenceConversion;

/// <summary>
/// Provides simple string replacement behaviour based on a mapping dictionary.
/// </summary>
public abstract class DictionaryBasedConverter : SequenceConverterBase
{
    protected abstract IReadOnlyDictionary<string, string> ModificationMappings { get; }

    protected override string ConvertSequenceCore(string fullSequence, SequenceConversionHandlingMode handlingMode,
        IList<string> warnings)
    {
        string converted = fullSequence;
        foreach (var kvp in ModificationMappings)
        {
            converted = converted.Replace(kvp.Key, kvp.Value);
        }

        return converted;
    }
}
