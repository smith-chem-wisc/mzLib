namespace Readers;

/// <summary>
/// Reads a MetaMorpheus transcript-group table. Its columns and row format are shared with the
/// protein-group table; this type exposes the RNA-specific supported file type.
/// </summary>
public class TranscriptGroupFromTsvFile : ProteinGroupFromTsvFile
{
    public override SupportedFileType FileType => SupportedFileType.MetaMorpheusQuantifiedTranscriptGroups;

    public TranscriptGroupFromTsvFile(string filePath) : base(filePath) { }

    public TranscriptGroupFromTsvFile() : base() { }
}
