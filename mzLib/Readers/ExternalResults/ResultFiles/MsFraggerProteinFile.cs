using CsvHelper;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public class MsFraggerProteinFile : ResultFile<MsFraggerProtein>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.MsFraggerProtein;
        public override Software Software { get; set; }

        public MsFraggerProteinFile(string filePath) : base(filePath, Software.MsFragger) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public MsFraggerProteinFile() : base() { }

        public override void LoadResults()
        {
            var csv = new CsvReader(new StreamReader(FilePath), MsFraggerProtein.CsvConfiguration);
            Results = csv.GetRecords<MsFraggerProtein>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), MsFraggerProtein.CsvConfiguration);

            csv.WriteHeader<MsFraggerProtein>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}
