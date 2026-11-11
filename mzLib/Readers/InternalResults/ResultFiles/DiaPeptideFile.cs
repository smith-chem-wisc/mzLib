using CsvHelper;

namespace Readers
{
    /// <summary>
    /// MetaMorpheus's DIA peptide table, AllDiaPeptides.tsv: one row per peptide across the runs searched together.
    /// MetaMorpheus writes it through this class and readers read it through this class, so the columns cannot drift
    /// apart (dataRepo's DATAREPO-65). Like <see cref="DiaPrecursorFile"/>, its first line states whether a contaminant
    /// database was searched. The DIA protein-group table needs no class of its own: it is MetaMorpheus's protein-group
    /// format, read by <see cref="ProteinGroupFromTsvFile"/>.
    /// </summary>
    public class DiaPeptideFile : ResultFile<DiaPeptideFromTsv>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.DiaPeptideTsv;
        public override Software Software { get; set; }

        /// <summary>Whether a contaminant database was searched; when false, no row can be labelled C.</summary>
        public bool ContaminantsAssessed { get; set; }

        public DiaPeptideFile(string filePath) : base(filePath, Software.MetaMorpheus) { }

        /// <summary>A table to write: the rows, and whether a contaminant database was searched.</summary>
        public DiaPeptideFile(string filePath, List<DiaPeptideFromTsv> results, bool contaminantsAssessed) : base(filePath, Software.MetaMorpheus)
        {
            Results = results;
            ContaminantsAssessed = contaminantsAssessed;
        }

        /// <summary>Constructor used to initialize from the factory method</summary>
        public DiaPeptideFile() : base() { }

        public override void LoadResults()
        {
            using var reader = new StreamReader(FilePath);
            ContaminantsAssessed = DiaContaminantsLine.Read(reader, FilePath, "DIA peptide table");
            using var csv = new CsvReader(reader, DiaPeptideFromTsv.CsvConfiguration);
            Results = csv.GetRecords<DiaPeptideFromTsv>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var writer = new StreamWriter(File.Create(outputPath));
            DiaContaminantsLine.Write(writer, ContaminantsAssessed);
            using var csv = new CsvWriter(writer, DiaPeptideFromTsv.CsvConfiguration);
            csv.WriteHeader<DiaPeptideFromTsv>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}
