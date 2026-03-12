

namespace Omics.SequenceConversion;

/// <summary>
/// Schema definition for the mzLib sequence format.
/// 
/// Format examples:
/// - Simple: "PEPTIDE"
/// - With residue modification: "PEP[Oxidation on M]TIDE" or "PEPTK[Common Fixed:Carbamidomethyl on C]IDE"
/// - With N-terminal modification: "[Common Biological:Acetylation on X]PEPTIDE" (note: no separator after N-term mod)
/// - With C-terminal modification: "PEPTIDE-[Common Artifact:Amidation on E]" (note: hyphen before bracket at end)
/// 
/// mzLib format characteristics:
/// - Uses square brackets [ ] for modification annotations
/// - N-terminal modifications directly precede sequence with no separator
/// - Uses hyphen "-" as C-terminal modification separator
/// - Modification identifiers use "IdWithMotif" format (e.g., "Oxidation on M")
/// - Can include modification type prefix (e.g., "Common Fixed:Carbamidomethyl on C")
/// 
/// Note: Mass shift notation (e.g., "[+15.995]") is handled by MassShiftSequenceFormatSchema, not this schema.
/// </summary>
public class MzLibSequenceFormatSchema : SequenceFormatSchema
{
    /// <summary>
    /// Singleton instance of the mzLib schema.
    /// </summary>
    public static readonly MzLibSequenceFormatSchema Instance = new();

    /// <summary>
    /// Private constructor to enforce singleton pattern.
    /// </summary>
    protected MzLibSequenceFormatSchema(char modOpen = '[', char modClosed = ']', string? nTermSeparator = "", string? cTermSeparator = "-")
        : base(modOpen, modClosed, nTermSeparator, cTermSeparator) { }

    /// <inheritdoc />
    public override string FormatName => "mzLib";

}
