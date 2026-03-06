namespace Omics.SequenceConversion;

/// <summary>
/// Schema definition for the mass shift sequence format.
/// 
/// Format examples:
/// - Simple: "PEPTIDE"
/// - With residue modification by mass: "PEP[+15.995]TIDE" (oxidation)
/// - With N-terminal modification: "[+42.011]PEPTIDE" (acetylation)
/// - With C-terminal modification: "PEPTIDE-[+0.984]" (amidation)
/// - Negative mass shifts: "PEP[-18.011]TIDE" (water loss)
/// 
/// Mass shift format characteristics:
/// - Uses square brackets [ ] for modification annotations
/// - Modifications are specified as mass shifts with sign (+ or -)
/// - N-terminal modifications directly precede sequence with no separator
/// - Uses hyphen "-" as C-terminal modification separator
/// - Format: [±mass] where mass is a decimal number
/// 
/// This format is useful when:
/// - Exact modification identities are unknown
/// - Working with raw mass spectrometry data
/// - Converting between systems with different modification databases
/// - Performing modifications not in standard databases
/// </summary>
public sealed class MassShiftSequenceFormatSchema : SequenceFormatSchema
{
    /// <summary>
    /// Singleton instance of the mass shift schema.
    /// </summary>
    public static readonly MassShiftSequenceFormatSchema Instance = new();

    /// <summary>
    /// Private constructor to enforce singleton pattern.
    /// </summary>
    private MassShiftSequenceFormatSchema() { }

    /// <inheritdoc />
    public override string FormatName => "MassShift";

    /// <inheritdoc />
    public override char ModOpenBracket => '[';

    /// <inheritdoc />
    public override char ModCloseBracket => ']';

    /// <inheritdoc />
    /// <remarks>
    /// Empty string indicates N-terminal modifications directly precede the sequence with no separator.
    /// Example: "[+42.011]PEPTIDE" not "[+42.011]-PEPTIDE"
    /// </remarks>
    public override string? NTermSeparator => "";

    /// <inheritdoc />
    public override string? CTermSeparator => "-";
}
