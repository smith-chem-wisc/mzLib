using CsvHelper;
using Easy.Common.Extensions;
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
            var headers = File.ReadLines(FilePath).First().Split('\t');
            // If it is a consensus file with additional information
            if (headers.Any(p => p.Contains("Spectral Count")))
            {
                var spectraCount = headers.Where(p => p.Contains("Spectral Count"))
                    .ToDictionary(p => p, p => headers.IndexOf(p));
                var intensityCount = headers.Where(p => p.Contains("Intensity"))
                    .ToDictionary(p => p, p => headers.IndexOf(p));
                var results = new List<MsFraggerPeptide>();
                bool readHeader = false;
                while (csv.Read())
                {
                    if (readHeader == false)
                    {
                        csv.ReadHeader();
                        readHeader = true;
                        continue;
                    }
                    var record = csv.GetRecord<MsFraggerPeptide>();
                    if (record is null)
                        continue;
                    foreach (var kvp in spectraCount)
                    {
                        var file = kvp.Key.Replace(" Spectral Count", "");
                        record.FileToPsmCount[file] = csv.GetField<int>(kvp.Value);
                    }
                    foreach (var kvp in intensityCount)
                    {
                        var file = kvp.Key.Replace(" Intensity", "");
                        record.IntensityByFile[file] = csv.GetField<double>(kvp.Value);
                    }
                    results.Add(record);
                }
                Results = results;
            }
            else
                Results = csv.GetRecords<MsFraggerPeptide>().ToList();
            //try
            //{
            //    string dirName = Path.GetDirectoryName(FilePath);
            //    Results.ForEach(p => p.FileNameWithoutExtension = Path.GetFileNameWithoutExtension(dirName));
            //}
            //catch (Exception e)
            //{
            //    Console.WriteLine(e);
            //}
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
