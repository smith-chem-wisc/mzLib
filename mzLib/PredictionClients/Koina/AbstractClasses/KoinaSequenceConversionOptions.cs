using Omics.Modifications;
using Omics.Modifications.Conversion;

namespace PredictionClients.Koina.AbstractClasses;

public sealed class KoinaSequenceConversionOptions
{
    public bool Enabled { get; init; }
    public ModificationNamingConvention TargetConvention { get; init; } = ModificationNamingConvention.Unimod;
    public SequenceConversionHandlingMode HandlingMode { get; init; } = SequenceConversionHandlingMode.RemoveIncompatibleMods;
}
