using CsvHelper;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    /// <summary>
    /// Concrete Product for reading and representing a experiment annotation file
    /// </summary>
    public class ExperimentAnnotationFile: ResultFile<ExperimentAnnotation>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.ExperimentAnnotation;

        public override Software Software { get; set; }

        public ExperimentAnnotationFile(string filePath) : base(filePath, Software.MsFragger) { }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public ExperimentAnnotationFile() : base() { }

        /// <summary>
        /// Load Results to the Results List from the given filepath
        /// </summary>
        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), ExperimentAnnotation.CsvConfiguration);
            Results = csv.GetRecords<ExperimentAnnotation>().ToList();
        }

        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), ExperimentAnnotation.CsvConfiguration);

            csv.WriteHeader<ExperimentAnnotation>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }
}