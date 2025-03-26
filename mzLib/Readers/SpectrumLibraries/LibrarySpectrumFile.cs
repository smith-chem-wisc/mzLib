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

        public override void WriteResults(string outputPath) => throw new NotImplementedException();
    }

}


