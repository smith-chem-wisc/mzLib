using CsvHelper;

namespace Readers
{
    /// <summary>
    /// MetaMorpheus's DIA precursor table, AllDiaPrecursors.tsv: one row per precursor per run. MetaMorpheus writes it
    /// through this class and readers read it through this class, so the columns cannot drift apart (dataRepo's
    /// DATAREPO-65). The first line states whether a contaminant database was searched, so that "no contaminants" is never
    /// true by default.
    /// </summary>
    public class DiaPrecursorFile : ResultFile<DiaPrecursorFromTsv>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.DiaPrecursorTsv;
        public override Software Software { get; set; }

        /// <summary>Whether a contaminant database was searched; when false, no row can be labelled C.</summary>
        public bool ContaminantsAssessed { get; set; }

        public DiaPrecursorFile(string filePath) : base(filePath, Software.MetaMorpheus) { }

        /// <summary>A table to write: the rows, and whether a contaminant database was searched.</summary>
        public DiaPrecursorFile(string filePath, List<DiaPrecursorFromTsv> results, bool contaminantsAssessed) : base(filePath, Software.MetaMorpheus)
        {
            Results = results;
            ContaminantsAssessed = contaminantsAssessed;
        }

        /// <summary>Constructor used to initialize from the factory method</summary>
        public DiaPrecursorFile() : base() { }

        public override void LoadResults()
        {
            using var reader = new StreamReader(FilePath);
            ContaminantsAssessed = DiaContaminantsLine.Read(reader, FilePath, "DIA precursor table");
            using var csv = new CsvReader(reader, DiaPrecursorFromTsv.CsvConfiguration);
            Results = csv.GetRecords<DiaPrecursorFromTsv>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var writer = new StreamWriter(File.Create(outputPath));
            DiaContaminantsLine.Write(writer, ContaminantsAssessed);
            using var csv = new CsvWriter(writer, DiaPrecursorFromTsv.CsvConfiguration);
            csv.WriteHeader<DiaPrecursorFromTsv>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}
