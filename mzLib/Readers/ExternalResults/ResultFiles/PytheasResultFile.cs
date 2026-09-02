using Readers.ExternalResults.IndividualResultRecords;
using System.IO;

namespace Readers.ExternalResults.ResultFiles;

/// <summary>
/// A Pytheas match-output file: a '#theoretical_digest'-headed text file whose match lines are
/// grouped under PRECURSOR_ION headers. See <see cref="PytheasResult"/> for the line layout.
/// </summary>
public class PytheasResultFile : ResultFile<PytheasResult>
{
    public override SupportedFileType FileType => SupportedFileType.PytheasResult;
    public override Software Software { get; set; }

    /// <summary>The '#' metadata lines at the top of the file, kept so writing round-trips them.</summary>
    public List<string> HeaderLines { get; } = new();

    public PytheasResultFile(string filePath) : base(filePath, Software.Pytheas) { }

    /// <summary>Constructor used to initialize from the factory method.</summary>
    public PytheasResultFile() : base() { }

    public override void LoadResults()
    {
        var results = new List<PytheasResult>();
        string precursorIon = null;

        foreach (var line in File.ReadLines(FilePath))
        {
            if (line.Length == 0)
                continue;

            if (line.StartsWith("#", StringComparison.Ordinal))
            {
                HeaderLines.Add(line);
                continue;
            }

            if (line.StartsWith("PRECURSOR_ION=", StringComparison.Ordinal))
            {
                precursorIon = line.Substring("PRECURSOR_ION=".Length);
                continue;
            }

            results.Add(PytheasResult.Parse(line, precursorIon));
        }

        Results = results;
    }

    public override void WriteResults(string outputPath)
    {
        if (!CanRead(outputPath))
            outputPath += FileType.GetFileExtension();

        // Reading Results forces the lazy load, which populates HeaderLines.
        var results = Results;

        using var writer = new StreamWriter(outputPath);

        foreach (var header in HeaderLines)
            writer.WriteLine(header);

        string lastPrecursorIon = null;
        foreach (var result in results)
        {
            if (result.PrecursorIon != lastPrecursorIon)
            {
                writer.WriteLine();
                writer.WriteLine("PRECURSOR_ION=" + result.PrecursorIon);
                lastPrecursorIon = result.PrecursorIon;
            }
            writer.WriteLine(result.RawLine);
        }
    }
}