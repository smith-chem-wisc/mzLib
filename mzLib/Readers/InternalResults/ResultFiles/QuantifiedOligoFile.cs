namespace Readers;

/// <summary>
/// Reads a FlashLFQ quantified-oligo table. Its columns and row format are shared with the
/// quantified-peptide table; this type exposes the RNA-specific supported file type.
/// </summary>
public class QuantifiedOligoFile : QuantifiedPeptideFile
{
    public override SupportedFileType FileType => SupportedFileType.FlashLFQQuantifiedOligo;

    public QuantifiedOligoFile(string filePath) : base(filePath) { }

    public QuantifiedOligoFile() : base() { }
}
