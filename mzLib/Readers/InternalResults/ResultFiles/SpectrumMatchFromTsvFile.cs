using System.Collections.Generic;
using System.Linq;

namespace Readers;

public abstract class SpectrumMatchFromTsvFile<T> : ResultFile<T>, IQuantifiableResultFile where T: SpectrumMatchFromTsv
{
    protected SpectrumMatchParsingParameters? ParsingParams { get; }
    public override SupportedFileType FileType => FilePath.ParseFileType();
    public override Software Software { get; set; }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    protected SpectrumMatchFromTsvFile(SpectrumMatchParsingParameters? parsingParams = null) : base()
    {
        ParsingParams = parsingParams;
    }

    protected SpectrumMatchFromTsvFile(string filePath, SpectrumMatchParsingParameters? parsingParams = null) : base(filePath, Software.MetaMorpheus)
    {
        ParsingParams = parsingParams;
    }

    public override void LoadResults()
    {
        Results = SpectrumMatchTsvReader.ReadTsv<T>(FilePath, out List<string> warnings, ParsingParams);
    }
    public override void WriteResults(string outputPath) => throw new NotImplementedException();
    public IEnumerable<IQuantifiableRecord> GetQuantifiableResults() => Results;
    public Dictionary<string, string> FileNameToFilePath(List<string> fullFilePath)
    {
        var fileNames = Results.Select(r => r.FileName).Distinct().ToList();
        var fileNameToPath = new Dictionary<string, string>();
        foreach (string fullPath in fullFilePath)
        {
            var fileName = Path.GetFileNameWithoutExtension(fullPath);
            if (fileNames.Contains(fileName))
            {
                fileNameToPath[fileName] = fullPath;
            }
        }
        return fileNameToPath;
    }

    protected internal IReadOnlyList<FdrBenchPeptide> ToFdrBenchPeptideRecords(FdrBenchWriterSettings? settings = null)
    {
        settings ??= FdrBenchWriterSettings.Default;
        var records = new List<FdrBenchPeptide>();

        foreach (var spectrumMatch in Results)
        {
            if (!settings.IncludeDecoys && spectrumMatch.IsDecoy)
            {
                continue;
            }

            var formattedSequence = FdrBenchSequenceFormatter.FormatSequence(
                spectrumMatch.FullSequence,
                spectrumMatch.BaseSequence,
                settings);

            records.Add(new FdrBenchPeptide
            {
                Run = spectrumMatch.FileNameWithoutExtension,
                Peptide = spectrumMatch.BaseSequence,
                ModifiedPeptide = formattedSequence,
                Charge = spectrumMatch.PrecursorCharge,
                QValue = spectrumMatch.QValue,
                Pep = double.IsNaN(spectrumMatch.PEP) ? null : spectrumMatch.PEP,
                Protein = FormatProteinAccessionForFdrBench(spectrumMatch.Accession),
                Score = spectrumMatch.Score
            });
        }

        return records;
    }

    protected internal IReadOnlyList<FdrBenchProtein> ToFdrBenchProteinRecords(FdrBenchWriterSettings? settings = null)
    {
        settings ??= FdrBenchWriterSettings.Default;
        var records = new List<FdrBenchProtein>();

        foreach (var spectrumMatch in Results)
        {
            if (!settings.IncludeDecoys && spectrumMatch.IsDecoy)
            {
                continue;
            }

            records.Add(new FdrBenchProtein
            {
                Protein = FormatProteinAccessionForFdrBench(spectrumMatch.Accession),
                QValue = spectrumMatch.QValue,
                Score = spectrumMatch.Score
            });
        }

        return records;
    }

    private static string FormatProteinAccessionForFdrBench(string accession)
    {
        if (string.IsNullOrWhiteSpace(accession))
        {
            return string.Empty;
        }

        var tokens = accession
            .Split('|', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);

        if (tokens.Length == 0)
        {
            tokens = new[] { accession };
        }

        return string.Join(';', tokens.Select(token =>
            token.EndsWith("_target", StringComparison.OrdinalIgnoreCase) ? token : $"{token}_target"));
    }
}

public class SpectrumMatchFromTsvFile : SpectrumMatchFromTsvFile<SpectrumMatchFromTsv> 
{

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public SpectrumMatchFromTsvFile() : base() { }

    public SpectrumMatchFromTsvFile(string filePath) : base(filePath) { }
}
