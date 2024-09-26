using CsvHelper;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public class MsFraggerPeptideFile : ResultFile<MsFraggerPeptide>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.MsFraggerPeptide;
        public override Software Software { get; set; }

        public MsFraggerPeptideFile(string filePath) : base(filePath, Software.Unspecified) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public MsFraggerPeptideFile() : base() { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), MsFraggerPeptide.CsvConfiguration);
            Results = csv.GetRecords<MsFraggerPeptide>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)),
                MsFraggerPeptide.CsvConfiguration);
            csv.WriteHeader<MsFraggerPeptide>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

    }
}
