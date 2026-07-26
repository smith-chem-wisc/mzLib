using CsvHelper;

namespace Readers
{
    /// <summary>
    /// DIA-NN's main report -- the long-format report.tsv written by DIA-NN 1.8/1.9, with one row
    /// per precursor per run. DIA-NN gives this file no distinguishing name or extension, so
    /// <see cref="SupportedFileTypeExtensions.ParseFileType"/> recognizes it by its header columns.
    /// </summary>
    public class DiaNnReportFile : ResultFile<DiaNnPrecursor>, IQuantifiableResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.DiaNnReport;
        public override Software Software { get; set; }

        public DiaNnReportFile(string filePath) : base(filePath, Software.DiaNn) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public DiaNnReportFile() : base() { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), DiaNnPrecursor.CsvConfiguration);
            Results = csv.GetRecords<DiaNnPrecursor>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            // Not CanRead: this type's "extension" is the whole conventional name "report.tsv", so
            // CanRead is false for any validly-named tsv that is not called *report.tsv, and
            // appending would turn "myResults.tsv" into "myResults.tsvreport.tsv".
            if (!outputPath.EndsWith(".tsv", StringComparison.InvariantCultureIgnoreCase))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), DiaNnPrecursor.CsvConfiguration);

            csv.WriteHeader<DiaNnPrecursor>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public IEnumerable<IQuantifiableRecord> GetQuantifiableResults() => Results;

        /// <summary>
        /// Creates a dictionary linking each DIA-NN run name to the full path of its spectra file.
        /// </summary>
        /// <param name="fullFilePath"> list of all full file paths associated with a given result </param>
        /// <returns> dictionary with key run name and value fullFilePath </returns>
        public Dictionary<string, string> FileNameToFilePath(List<string> fullFilePath)
        {
            List<string> runNames = Results.Select(precursor => precursor.FileName).Distinct().ToList();
            fullFilePath = fullFilePath.Distinct().ToList();
            Dictionary<string, string> allFiles = new Dictionary<string, string>();

            foreach (var runName in runNames)
            {
                string match = null;
                int matchLength = -1;
                bool ambiguous = false;

                foreach (var file in fullFilePath)
                {
                    string fileNameWithoutExtension = DiaNnPrecursor.StripPathAndExtension(file);

                    // An empty name is a prefix of everything and would match every run
                    if (fileNameWithoutExtension.Length == 0)
                        continue;

                    // An exact match is unambiguous, so take it and stop looking
                    if (fileNameWithoutExtension.Equals(runName, StringComparison.InvariantCultureIgnoreCase))
                    {
                        match = file;
                        ambiguous = false;
                        break;
                    }

                    if (!IsSuffixExtensionOf(runName, fileNameWithoutExtension))
                        continue;

                    // Prefer the longest candidate: given files "sample" and "sample_rep1" for run
                    // "sample_rep1_centroid", the longer one shares more of the name and is the
                    // likelier pairing. Equal-length candidates give no reason to prefer either.
                    if (fileNameWithoutExtension.Length > matchLength)
                    {
                        match = file;
                        matchLength = fileNameWithoutExtension.Length;
                        ambiguous = false;
                    }
                    else if (fileNameWithoutExtension.Length == matchLength)
                    {
                        ambiguous = true;
                    }
                }

                // Leaving an ambiguous run unmatched makes MakeIdentifications throw by name, which
                // is far easier to diagnose than quantifying against the wrong file
                if (match is not null && !ambiguous)
                    allFiles.Add(runName, match);
            }

            // One spectra file cannot be the source of two runs. If it were claimed twice, FlashLFQ
            // would merge both runs' identifications into a single SpectraFileInfo and quantify them
            // as one file, so drop the whole contested set rather than pick a winner.
            foreach (var contested in allFiles.GroupBy(pair => pair.Value).Where(group => group.Count() > 1).ToList())
                foreach (var pair in contested)
                    allFiles.Remove(pair.Key);

            return allFiles;
        }

        /// <summary>
        /// True when <paramref name="runName"/> is <paramref name="fileName"/> followed by a
        /// delimited suffix, which is the shape DIA-NN produces when the search ran on a converted
        /// copy of the caller's file ("sample" -> run "sample_centroid").
        /// <para>
        /// The suffix has to start at a delimiter. Plain containment would make the file "A1" match
        /// the run "A10", and 96-well plate naming produces exactly that collision.
        /// </para>
        /// </summary>
        private static bool IsSuffixExtensionOf(string runName, string fileName)
        {
            // Either side may carry the extra suffix: the caller may hold the converted file instead
            string longer = runName.Length > fileName.Length ? runName : fileName;
            string shorter = runName.Length > fileName.Length ? fileName : runName;

            if (!longer.StartsWith(shorter, StringComparison.InvariantCultureIgnoreCase))
                return false;

            char boundary = longer[shorter.Length];
            return boundary is '_' or '-' or '.' or ' ';
        }
    }
}
