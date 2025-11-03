using CsvHelper;
using CsvHelper.Configuration;
using Omics;
using Readers.ExternalResults.IndividualResultRecords;
using System.Globalization;

namespace Readers.ExternalResults.ResultFiles;

public class CasanovoMzTabFile : ResultFile<CasanovoMzTabRecord>, IResultFile
{
    private static readonly CsvConfiguration CsvConfig = new CsvConfiguration(CultureInfo.InvariantCulture)
    {
        Delimiter = "\t",
        HasHeaderRecord = true,
        IgnoreBlankLines = true,
        TrimOptions = TrimOptions.Trim,
        BadDataFound = null,
        MissingFieldFound = null,
        HeaderValidated = null,
    };

    public bool LoadModificationInformation { get; init; }
    public override SupportedFileType FileType => SupportedFileType.CasanovoMzTab;
    public override Software Software { get; set; } = Software.Casanovo;

    public CasanovoMzTabFile(string filePath, bool loadModInformation = true) : base(filePath, Software.Casanovo)
    {
        LoadModificationInformation = loadModInformation;
        using var reader = new StreamReader(FilePath);
        string? line;

        // 1. Read header and ms_run info
        while ((line = reader.ReadLine()) != null)
        {
            if (line.StartsWith("MTD"))
            {
                HeaderLines.Add(line);

                // Parse ms_run location
                var split = line.Split('\t');
                if (split.Length >= 3 && split[1].StartsWith("ms_run["))
                {
                    var msrunKey = split[1].Split('-')[0]; // e.g., ms_run[1]
                    if (split[1].EndsWith("location"))
                    {
                        MsRunDictionary[msrunKey] = split[2].Replace("file:///", "").Replace('/', Path.DirectorySeparatorChar);
                    }
                }
            }
            else if (line.StartsWith("PSH"))
            {
                break;
            }
        }
    }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public CasanovoMzTabFile(bool loadModInformation = true) : base()
    {
        LoadModificationInformation = loadModInformation;
    }

    /// <summary>
    /// Stores all header lines (MTD) as read from the file.
    /// </summary>
    public List<string> HeaderLines { get; private set; } = new();

    /// <summary>
    /// Maps ms_run indices (e.g., ms_run[1]) to file paths.
    /// </summary>
    public Dictionary<string, string> MsRunDictionary { get; private set; } = new();

    public override void LoadResults()
    {
        string? line;
        var psmLines = new List<string>();

        // Read PSM table lines
        using var reader = new StreamReader(FilePath);
        while ((line = reader.ReadLine()) != null)
        {
            if (line.StartsWith("MTD"))
                continue;
            psmLines.Add(line);
        }

        // Parse PSM table with CsvHelper
        if (psmLines.Count > 1)
        {
            using var psmReader = new StringReader(string.Join(Environment.NewLine, psmLines));
            using var csv = new CsvReader(psmReader, CsvConfig);
            var parsedRecords = csv.GetRecords<CasanovoMzTabRecord>().ToList();

            // 4. Assign file name and scan number from spectra_ref
            foreach (var rec in parsedRecords)
            {
                if (!string.IsNullOrWhiteSpace(rec.SpectraRef))
                {
                    // Example: ms_run[1]:index=0
                    var parts = rec.SpectraRef.Split(':');
                    if (parts.Length == 2)
                    {
                        var msrunKey = parts[0];
                        var indexPart = parts[1];
                        if (MsRunDictionary.TryGetValue(msrunKey, out var filePath))
                        {
                            rec.FileNameWithoutExtension = Path.GetFileNameWithoutExtension(filePath);
                        }
                        if (indexPart.StartsWith("index=") && int.TryParse(indexPart.Substring(6), out int scanNum))
                        {
                            // mzTab indices are zero-based, but you may want to store as one-based
                            rec.OneBasedScanNumber = scanNum + 1;
                        }
                    }
                }

                rec.BaseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(rec.ModifiedSequence);
                if (LoadModificationInformation)
                {
                    rec.AllModsOneIsNterminus = ModificationConverter.GetModificationDictionaryFromFullSequence(rec.ModifiedSequence);
                    rec.FullSequence = IBioPolymerWithSetMods.DetermineFullSequence(rec.BaseSequence, rec.AllModsOneIsNterminus);
                }
            }

            Results = parsedRecords;
        }
        else
        {
            Results = new List<CasanovoMzTabRecord>();
        }
    }

    public override void WriteResults(string outputPath)
    {
        if (!CanRead(outputPath))
            outputPath += FileType.GetFileExtension();

        using var writer = new StreamWriter(outputPath);

        int headerCount = 0;

        // 1. Write header lines
        if (HeaderLines.Any())
        {
            foreach (var header in HeaderLines)
            {
                writer.WriteLine(header);
                headerCount++;
            }
        }

        // 2. Write PSM table
        if (Results.Any())
        {
            using var buffer = new StringWriter();
            using (var csv = new CsvWriter(buffer, CsvConfig))
            {
                csv.WriteHeader<CasanovoMzTabRecord>();
                csv.NextRecord();
                foreach (var rec in Results)
                {
                    csv.WriteRecord(rec);
                    csv.NextRecord();
                }
            }

            var lines = buffer.ToString().Replace("\r\n", "\n").Split('\n');
            for (int i = 0; i < lines.Length; i++)
            {
                var l = lines[i];
                if (string.IsNullOrWhiteSpace(l))
                    continue;
                if (i == 0)
                    writer.WriteLine("PSH\t" + l);
                else
                    writer.WriteLine("PSM\t" + l);
            }
        }
    }
}