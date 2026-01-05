namespace Omics.Modifications;

/// <summary>
/// Represents the naming convention to use for modification identifiers
/// </summary>
public enum ModificationNamingConvention
{
    /// <summary>
    /// Any of the available naming conventions
    /// </summary>
    Mixed,

    /// <summary>
    /// MetaMorpheus-style naming: includes modification type prefix
    /// Example: "Common Biological:Phosphorylation on S"
    /// </summary>
    MetaMorpheus,
    
    /// <summary>
    /// UniProt-style naming: uses "UniProt:" prefix
    /// Example: "UniProt:Phosphoserine on S"
    /// </summary>
    UniProt
}
