using CsvHelper;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using Readers.ExternalResults.IndividualResultRecords;

namespace Readers.ExternalResults.ResultFiles
{
    internal class MsFraggerSpeclibFile : ResultFile<MsFraggerSpeclib>, IResultFile
    {
        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }

        public MsFraggerSpeclibFile(string filePath) : base(filePath, Software.MsFragger)
        {
            FileType = filePath.ParseFileType();
        }

        public MsFraggerSpeclibFile() : base()
        {
            FileType = FilePath.IsNullOrEmpty() ? SupportedFileType.MsFraggerSpeclib : FilePath.ParseFileType();
        }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), MsFraggerSpeclib.CsvConfiguration);
            Results = csv.GetRecords<MsFraggerSpeclib>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using (var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), MsFraggerSpeclib.CsvConfiguration))
            {
                csv.WriteHeader<MsFraggerSpeclib>();
                foreach (var result in Results)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }
        }
    }
    
}
