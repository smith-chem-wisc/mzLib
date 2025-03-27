namespace Readers;

/// <summary>
/// Class Representing a directory of Chimerys results. Each file type has a path, and the file is loaded when requested.
/// Each file is optional in their output and can be null.
/// </summary>
public class ChimerysResultDirectory
{
    public readonly string DirectoryPath;
    public readonly string? PsmPath;
    public readonly string? PeptidePath;
    public readonly string? ModifiedPeptidePath;
    public readonly string? ProteinGroupPath;
    public readonly string? PrecursorPath;

    private ChimerysPsmFile? _psmFile;
    private ChimerysPeptideFile? _peptideFile;
    private ChimerysModifiedPeptideFile? _modifiedPeptideFile;
    private ChimerysProteinGroupFile? _proteinGroupFile;
    private ChimerysPrecursorFile? _precursorFile;

    public ChimerysResultDirectory(string directoryPath)
    {
        DirectoryPath = directoryPath;
        var allFiles = Directory.GetFiles(directoryPath);

        foreach (var filePath in allFiles)
        {
            switch (filePath.ParseFileType())
            {
                case SupportedFileType.ChimerysPsm:
                    PsmPath = filePath;
                    break;
                case SupportedFileType.ChimerysPeptide:
                    PeptidePath = filePath;
                    break;
                case SupportedFileType.ChimerysModifiedPeptide:
                    ModifiedPeptidePath = filePath;
                    break;
                case SupportedFileType.ChimerysProteinGroup:
                    ProteinGroupPath = filePath;
                    break;
                case SupportedFileType.ChimerysPrecursor:
                    PrecursorPath = filePath;
                    break;
                default:
                    continue;
            }
        }
    }

    public ChimerysPsmFile? PsmFile
    {
        get
        {
            if (_psmFile is not null) return _psmFile;
            if (PsmPath is null)
                return null;
            _psmFile = new ChimerysPsmFile(PsmPath);
            return _psmFile;
        }
    }

    public ChimerysPeptideFile? PeptideFile
    {
        get
        {
            if (_peptideFile is not null) return _peptideFile;
            if (PeptidePath is null)
                return null;
            _peptideFile = new ChimerysPeptideFile(PeptidePath);
            return _peptideFile;
        }
    }

    public ChimerysModifiedPeptideFile? ModifiedPeptideFile
    {
        get
        {
            if (_modifiedPeptideFile is not null) return _modifiedPeptideFile;
            if (ModifiedPeptidePath is null)
                return null;
            _modifiedPeptideFile = new ChimerysModifiedPeptideFile(ModifiedPeptidePath);
            return _modifiedPeptideFile;
        }
    }

    public ChimerysProteinGroupFile? ProteinGroupFile
    {
        get
        {
            if (_proteinGroupFile is not null) return _proteinGroupFile;
            if (ProteinGroupPath is null)
                return null;
            _proteinGroupFile = new ChimerysProteinGroupFile(ProteinGroupPath);
            return _proteinGroupFile;
        }
    }

    public ChimerysPrecursorFile? PrecursorFile
    {
        get
        {
            if (_precursorFile is not null) return _precursorFile;
            if (PrecursorPath is null)
                return null;
            _precursorFile = new ChimerysPrecursorFile(PrecursorPath);
            return _precursorFile;
        }
    }
}