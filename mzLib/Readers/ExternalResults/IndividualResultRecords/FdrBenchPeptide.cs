using System.Globalization;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers;

/// <summary>
/// Represents a single peptide/precursor-level entry in FDRBench format for FDP estimation.
/// </summary>
/// <remarks>
/// FDRBench requires the following column format:
/// <list type="table">
///     <listheader>
///         <term>Column</term>
///         <description>Description</description>
///     </listheader>
///     <item>
///         <term>run</term>
///         <description>File or run identifier (e.g., raw file name)</description>
///     </item>
///     <item>
///         <term>peptide</term>
///         <description>Base peptide sequence without modifications</description>
///     </item>
///     <item>
///         <term>mod_peptide</term>
///         <description>
///             Modified peptide sequence using UniMod notation. Modifications are encoded as
///             (UniMod:ID) where ID is the UniMod accession number.
///             Example: AAAPAPEEEMDEC(UniMod:4)EQALAAEPK for acetylation (UniMod:4) on Cysteine.
///         </description>
///     </item>
///     <item>
///         <term>charge</term>
///         <description>Precursor charge state</description>
///     </item>
///     <item>
///         <term>q_value</term>
///         <description>FDR (q-value) reported by the search engine</description>
///     </item>
///     <item>
///         <term>PEP</term>
///         <description>Posterior Error Probability (optional)</description>
///     </item>
///     <item>
///         <term>protein</term>
///         <description>
///             Protein identifier(s). For target peptides, append _target to the accession.
///             Multiple proteins should be separated by semicolons.
///             Example: P12345_target or P12345;P67890_target
///         </description>
///     </item>
///     <item>
///         <term>score</term>
///         <description>
///             Score used for ranking. Higher scores indicate better matches when using
///             default ranking (score:1 in FDRBench). The score direction can be specified
///             as score:0 (lower is better) or score:1 (higher is better).
///         </description>
///     </item>
/// </list>
/// 
/// Important: Decoy hits should NOT be included in the input file for FDP calculation.
/// </remarks>
public class FdrBenchPeptide
{
    public static CsvConfiguration CsvConfiguration { get; } = new CsvConfiguration(CultureInfo.InvariantCulture)
    {
        Delimiter = "\t",
        HasHeaderRecord = true,
        TrimOptions = TrimOptions.Trim,
        BadDataFound = null,
    };

    /// <summary>File or run identifier (e.g., raw file name)</summary>
    [Name("run")]
    public string Run { get; set; } = string.Empty;

    /// <summary>Base peptide sequence without modifications</summary>
    [Name("peptide")]
    public string Peptide { get; set; } = string.Empty;

    /// <summary>Modified peptide sequence using UniMod notation (e.g., AAAPAPEEEMDEC(UniMod:4)EQALAAEPK)</summary>
    [Name("mod_peptide")]
    public string ModifiedPeptide { get; set; } = string.Empty;

    /// <summary>Precursor charge state</summary>
    [Name("charge")]
    public int Charge { get; set; }

    /// <summary>FDR (q-value) reported by the search engine</summary>
    [Name("q_value")]
    public double QValue { get; set; }

    /// <summary>Posterior Error Probability (optional)</summary>
    [Name("PEP")]
    public double? Pep { get; set; }

    /// <summary>Protein identifier(s) with _target suffix for target peptides</summary>
    [Name("protein")]
    public string Protein { get; set; } = string.Empty;

    /// <summary>Score used for ranking (higher = better by default)</summary>
    [Name("score")]
    public double Score { get; set; }
}
