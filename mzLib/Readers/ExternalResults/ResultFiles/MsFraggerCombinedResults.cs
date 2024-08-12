using CsvHelper;
using Readers.ExternalResults.BaseClasses;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.ExternalResults.ResultFiles
{
    public class MsFraggerCombinedResults : ResultFile<MsFraggerPsm>, IResultFile, IQuantifiableResultFile
    {
        public string fullFolderPath;
        public List<MsFraggerPsmFile> allFiles;
        private List<string> allFilePaths;

        public override SupportedFileType FileType => SupportedFileType.MsFraggerPsm;
        public override Software Software { get; set; }
        public MsFraggerCombinedResults(string filePath) : base(filePath, Software.MsFragger) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public MsFraggerCombinedResults() : base() { }

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

        public Dictionary<string, string> FileNameToFilePath(List<string> fullFilePath)
        {
            List<string> rawFileNames = Results.Select(psm => psm.FileName).Distinct().ToList();
            fullFilePath = fullFilePath.Distinct().ToList();
            Dictionary<string, string> allFiles = new Dictionary<string, string>();

            foreach (var fileName in rawFileNames)
            {
                string shortFileName = Path.GetFileNameWithoutExtension(fileName);
                if (shortFileName.Contains("."))
                {
                    shortFileName = Path.GetFileNameWithoutExtension(shortFileName);
                }

                foreach (var file in fullFilePath)
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

        public Dictionary<string, string> FileNameToFilePath()
        {
            return null;
        }

        /// <summary>
        /// Adds all the MSFragger files of each sample to AllFiles
        /// </summary>
        public void FindAllFiles()
        {

        }

        /// <summary>
        /// Adds the path to each MSFragger file to AllFilePaths
        /// </summary>
        private void FindAllFilesPaths()
        {

        }


    }
}
