using CsvHelper;
using Omics.SpectrumMatch;

namespace Readers.SpectrumLibraries
{
    public class LibrarySpectrumFile : ResultFile<LibrarySpectrum>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.msp;
        public override Software Software { get; set; }

        public LibrarySpectrumFile() : base()
        {
        }

        public LibrarySpectrumFile(string filePath) : base(filePath, Software.MetaMorpheus)
        {
        }

        public override void LoadResults()
        {
            Results = MspSpectrumLibraryReader.ReadMspMsp(FilePath, out List<string> warnings);
        }
        /// <summary>
        /// Writes results to a specific output path
        /// </summary>
        /// <param name="outputPath">destination path</param>
        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var fs = new FileStream(outputPath, FileMode.Create, FileAccess.Write);


            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), FlashDeconvTsv.CsvConfiguration);

            csv.WriteHeader<FlashDeconvTsv>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }
    }

}


