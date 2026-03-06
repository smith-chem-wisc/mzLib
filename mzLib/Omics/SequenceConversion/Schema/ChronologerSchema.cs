namespace Omics.SequenceConversion;

/// <summary>
/// Schema definition for the Chronologer sequence format.
/// 
/// Chronologer is a deep learning model for retention time prediction that uses
/// a specialized single-character encoding for amino acids and their modifications.
/// 
/// Format characteristics:
/// - Uses single-character codes for modified residues (lowercase letters)
/// - Uses special characters for N/C-terminus states
/// - Maximum sequence length: 50 amino acids (+ 2 for termini)
/// - Does NOT use bracket-based modification annotations
/// 
/// Alphabet (52 + 2 positions):
/// - Positions 1-20: Canonical amino acids (ACDEFGHIKLMNPQRSTVWY)
/// - Positions 21-37: Modified residues (cmdestyabunopqrxz)
/// - Positions 38-44: N/C terminus states (-^()&*_)
/// - Positions 45-54: User-defined slots (0123456789)
/// 
/// Modification encoding:
/// - m: Oxidation on M [+15.995]
/// - c: Carbamidomethyl on C [+57.021]
/// - d: Alternative C mod [+39.99]
/// - e: PyroGlu from E/Q [-18.01/-17.02]
/// - s: Phosphorylation on S [+79.966]
/// - t: Phosphorylation on T [+79.966]
/// - y: Phosphorylation on Y [+79.966]
/// - a: Acetylation on K [+42.011]
/// - b: Succinylation on K [+100.0]
/// - u: Ubiquitination on K [+114.0]
/// - n: Methylation on K [+14.016]
/// - o: Dimethylation on K [+28.031]
/// - p: Trimethylation on K [+42.047]
/// - q: Methylation on R [+14.016]
/// - r: Dimethylation on R [+28.031]
/// - z: GlyGly on K [+224.1]
/// - x: Heavy GlyGly on K [+229.1]
/// 
/// N-terminus states:
/// - '-': Free N-terminus
/// - '^': N-terminal acetylation [+42.01]
/// - ')': PyroGlu at N-terminus (from E)
/// - '(': Cyclized CAM-Cys at N-terminus
/// - '&': N-terminal GlyGly [+224.1]
/// - '*': N-terminal heavy GlyGly [+229.1]
/// 
/// C-terminus state:
/// - '_': C-terminus
/// </summary>
public sealed class ChronologerSequenceFormatSchema : SequenceFormatSchema
{
    /// <summary>
    /// Singleton instance of the Chronologer schema.
    /// </summary>
    public static readonly ChronologerSequenceFormatSchema Instance = new();

    /// <summary>
    /// Private constructor to enforce singleton pattern.
    /// </summary>
    private ChronologerSequenceFormatSchema() { }

    /// <inheritdoc />
    public override string FormatName => "Chronologer";

    /// <inheritdoc />
    /// <remarks>
    /// Chronologer doesn't use bracket-based modification annotations.
    /// This is set to '[' for schema compatibility only.
    /// </remarks>
    public override char ModOpenBracket => '[';

    /// <inheritdoc />
    /// <remarks>
    /// Chronologer doesn't use bracket-based modification annotations.
    /// This is set to ']' for schema compatibility only.
    /// </remarks>
    public override char ModCloseBracket => ']';

    /// <inheritdoc />
    /// <remarks>
    /// Chronologer encodes N-terminal modifications as special characters,
    /// not as separate bracketed annotations, so no separator is used.
    /// </remarks>
    public override string? NTermSeparator => null;

    /// <inheritdoc />
    /// <remarks>
    /// Chronologer encodes C-terminal state as a special character,
    /// not as a separate bracketed annotation, so no separator is used.
    /// </remarks>
    public override string? CTermSeparator => null;

    /// <summary>
    /// Maximum sequence length supported by Chronologer (excluding termini tokens).
    /// </summary>
    public const int MaxSequenceLength = 50;

    /// <summary>
    /// Total encoded length including N and C terminus tokens.
    /// </summary>
    public const int EncodedLength = MaxSequenceLength + 2;

    /// <summary>
    /// The canonical amino acids recognized by Chronologer.
    /// </summary>
    public const string CanonicalAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>
    /// Free N-terminus token.
    /// </summary>
    public const char FreeNTerminus = '-';

    /// <summary>
    /// C-terminus token.
    /// </summary>
    public const char CTerminus = '_';

    /// <summary>
    /// N-terminal acetylation token.
    /// </summary>
    public const char NTermAcetyl = '^';

    /// <summary>
    /// N-terminal PyroGlu token.
    /// </summary>
    public const char NTermPyroGlu = ')';

    /// <summary>
    /// N-terminal cyclized CAM-Cys token.
    /// </summary>
    public const char NTermCyclizedCamCys = '(';

    /// <summary>
    /// N-terminal GlyGly token.
    /// </summary>
    public const char NTermGlyGly = '&';

    /// <summary>
    /// N-terminal heavy GlyGly token.
    /// </summary>
    public const char NTermHeavyGlyGly = '*';
}
