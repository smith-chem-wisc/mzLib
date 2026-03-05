using Omics.Modifications;
using Omics.Modifications.Conversion;

namespace UsefulProteomicsDatabases;

/// <summary>
/// Configuration options for enabling SequenceConverter-based modification normalization
/// when exporting protein databases.
/// </summary>
public sealed class ProteinDbWriterConversionOptions
{
    /// <summary>
    /// Enables conversion when set to <c>true</c>.
    /// </summary>
    public bool Enabled { get; init; }

    /// <summary>
    /// Target modification naming convention. Defaults to UniProt.
    /// </summary>
    public ModificationNamingConvention TargetConvention { get; init; } = ModificationNamingConvention.UniProt;

    /// <summary>
    /// Handling mode for incompatible modifications. Defaults to removing them.
    /// </summary>
    public SequenceConversionHandlingMode HandlingMode { get; init; } = SequenceConversionHandlingMode.RemoveIncompatibleMods;

    /// <summary>
    /// Optional converter override. Defaults to <see cref="SequenceConverter.Default"/> when null.
    /// </summary>
    public ISequenceConverter? Converter { get; init; }
}
