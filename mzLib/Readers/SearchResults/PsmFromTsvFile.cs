namespace Readers;

public class PsmFromTsvFile : ResultFile<PsmFromTsv>, IResultFile, IQuantifiableResultFile
{
    public override SupportedFileType FileType => SupportedFileType.psmtsv;
    public override Software Software { get; set; }
    public IEnumerable<IQuantifiableRecord> GetQuantifiableResults() => Results;

    /// <summary>
    /// Links the file name associated with the protein to the raw file path of MassSpec data
    /// </summary>
    /// <param name="fullFilePaths"> list of file paths associated with each distinct record </param>
    /// <returns> Dictionary of file names and their associted full paths </returns>
    public Dictionary<string, string> FileNameToFilePath(List<string> fullFilePaths)
    {
        var fileNames = Results.Select(r => r.FileName).Distinct();
        Dictionary<string, string> fileNameFilePathDictionary = new Dictionary<string, string>();
        foreach(var fileName in fileNames)
        {
            var fullPath = fullFilePaths.First(p => p.Contains(fileName));
            fileNameFilePathDictionary.Add(fileName, fullPath);
        }
        return fileNameFilePathDictionary;
    }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public PsmFromTsvFile() : base() { }

    public PsmFromTsvFile(string filePath) : base(filePath, Software.MetaMorpheus) { }

    public override void LoadResults()
    {
        Results = SpectrumMatchTsvReader.ReadPsmTsv(FilePath, out List<string> warnings);
    }

    public override void WriteResults(string outputPath) => throw new NotImplementedException();
}