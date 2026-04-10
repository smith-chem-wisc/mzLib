using System.Globalization;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers;

/// <summary>
/// Represents a single protein-level entry in FDRBench format for protein-level FDP estimation.
/// </summary>
/// <remarks>
/// FDRBench requires the following column format for protein-level FDP calculation:
/// <list type="table">
///     <listheader>
///         <term>Column</term>
///         <description>Description</description>
///     </listheader>
///     <item>
///         <term>protein</term>
///         <description>
///             Protein identifier (multiple IDs separated by semicolons). Each row represents a protein group.
///             For target proteins, use the protein accession as-is.
///             Example: P62857 or P12345;P67890
///         </description>
///     </item>
///     <item>
///         <term>q_value</term>
///         <description>FDR (q-value) reported by the search engine</description>
///     </item>
///     <item>
///         <term>score</term>
///         <description>
///             Protein score used for ranking. Higher scores indicate better matches when using
///             default ranking (score:1 in FDRBench). The score direction can be specified
///             as score:0 (lower is better) or score:1 (higher is better).
///         </description>
///     </item>
/// </list>
/// 
/// Note: For protein-level analysis, only the three columns above are required.
/// Decoy hits should NOT be included in the input file for FDP calculation.
/// </remarks>
public class FdrBenchProtein
{
    public static CsvConfiguration CsvConfiguration { get; } = new CsvConfiguration(CultureInfo.InvariantCulture)
    {
        Delimiter = "\t",
        HasHeaderRecord = true,
        TrimOptions = TrimOptions.Trim,
        BadDataFound = null,
    };

    /// <summary>Protein identifier(s), multiple IDs separated by semicolons</summary>
    [Name("protein")]
    public string Protein { get; set; } = string.Empty;

    /// <summary>FDR (q-value) reported by the search engine</summary>
    [Name("q_value")]
    public double QValue { get; set; }

    /// <summary>Protein score used for ranking (higher = better by default)</summary>
    [Name("score")]
    public double Score { get; set; }
}
