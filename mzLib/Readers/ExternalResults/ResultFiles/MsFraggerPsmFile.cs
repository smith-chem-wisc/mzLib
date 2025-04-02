using CsvHelper;

namespace Readers
{
    public class MsFraggerPsmFile : ResultFile<MsFraggerPsm>, IQuantifiableResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.MsFraggerPsm;
        public override Software Software { get; set; }
        public MsFraggerPsmFile(string filePath) : base(filePath, Software.MsFragger) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public MsFraggerPsmFile() : base() { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), MsFraggerPsm.CsvConfiguration);
            Results = csv.GetRecords<MsFraggerPsm>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), MsFraggerPsm.CsvConfiguration);

            csv.WriteHeader<MsFraggerPsm>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

        public IEnumerable<IQuantifiableRecord> GetQuantifiableResults() => Results;

        /// <summary>
        /// Creates a dictionary linking a shortened file name to its corresponding full file path
        /// </summary>
        /// <param name="fullFilePath"> list of all full file paths associted with a given result </param>
        /// <returns> dictionary with key fileName and value fullFilePath </returns>
        public Dictionary<string, string> FileNameToFilePath (List<string> fullFilePath)
        {
            List<string> rawFileNames = Results.Select(psm => psm.FileName).Distinct().ToList();
            fullFilePath = fullFilePath.Distinct().ToList();
            Dictionary<string, string> allFiles = new Dictionary<string, string>();

            foreach(var fileName in rawFileNames)
            {
                string shortFileName = Path.GetFileName(fileName);

                // MSFragger results append the raw file with "interact-" and replace .raw with .pep.xml
                // In order to correctly match the file names, these changes must be removed
                shortFileName = shortFileName.Replace("interact-", "").Replace(".pep.xml", "");

                foreach(var file in fullFilePath)
                {
                    if (file.Contains(shortFileName) && !allFiles.ContainsKey(fileName))
                    {
                        allFiles.Add(fileName, file);
                        break;
                    }
                }
            }

            return allFiles;
        }
    }
}