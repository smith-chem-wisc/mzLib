using CsvHelper;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public class MaxQuantEvidenceFile : ResultFile<MaxQuantEvidence>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.MaxQuantEvidence;
        public override Software Software { get; set; }
        public MaxQuantEvidenceFile(string filePath) : base(filePath, Software.MaxQuant) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public MaxQuantEvidenceFile() : base() { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), MaxQuantEvidence.CsvConfiguration);
            Results = csv.GetRecords<MaxQuantEvidence>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), MaxQuantEvidence.CsvConfiguration);

            csv.WriteHeader<MaxQuantEvidence>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}
