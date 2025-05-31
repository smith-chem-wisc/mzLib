using Omics;
using Omics.Fragmentation;

namespace Readers.Puf;

internal class PufDirectoryResultFile : ResultFile<PufMsMsExperiment>
{
    public List<PufResultFile> PufFiles { get; } = new();
    public override SupportedFileType FileType => SupportedFileType.PufDirectory;
    public override Software Software { get; set; } = Software.ProsightPC;

    public PufDirectoryResultFile() : base("") { }
    public PufDirectoryResultFile(string directoryPath) : base(directoryPath) { }

    public override void LoadResults()
    {
        var localResults = new List<PufMsMsExperiment>();
        foreach (var file in Directory.GetFiles(FilePath, "*.puf"))
        {
            var pufFile = new PufResultFile(file);
            pufFile.LoadResults();
            PufFiles.Add(pufFile);
            localResults.AddRange(pufFile.Results); // or aggregate identifications as needed
        }
        Results = localResults;
    }

    public override void WriteResults(string outputPath)
    {
        Directory.CreateDirectory(outputPath);
        foreach (var pufFile in PufFiles)
        {
            var specificFileName = pufFile.FilePath.Replace(Path.GetDirectoryName(pufFile.FilePath) ?? "", "").TrimStart(Path.DirectorySeparatorChar);
            var outputFilePath = Path.Combine(outputPath, specificFileName);
            pufFile.WriteResults(outputFilePath);
        }
    }
}