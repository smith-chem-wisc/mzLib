using CsvHelper;
using Omics.SpectrumMatch;

namespace Readers.SpectrumLibraries
{
    public class MspSpectrumFile : SpectrumLibraryFile
    {
        public override SupportedFileType FileType => SupportedFileType.msp;
        public override Software Software { get; set; }

        public MspSpectrumFile() : base()
        {
        }

        public MspSpectrumFile(string filePath) : base(filePath, Software.MetaMorpheus)
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
        public override void WriteResults(string outputPath) => WriteMsp(outputPath);

    }

}


