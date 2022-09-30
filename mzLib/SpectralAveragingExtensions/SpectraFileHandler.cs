using Readers;
using MassSpectrometry;
using MzLibUtil;
using SourceFile = Readers.SourceFile;

namespace SpectralAveragingExtensions
{
    public static class SpectraFileHandler
    {
        /// <summary>
        /// Creates a List of MsDataScans from a spectra file
        /// Currently supports MzML and raw
        /// </summary>
        /// <param name="filepath"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException">thrown if file type is not supported</exception>
        public static List<MsDataScan> LoadAllScansFromFile(string filepath)
        {
            var reader = ReaderCreator.CreateReader(filepath);
            reader.LoadAllStaticData();
            return reader.GetAllScansList();
        }

        /// <summary>
        /// Gets the source file for the spectra file at designated path
        /// </summary>
        /// <param name="filepath"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException">throws if file type is not supported</exception>
        public static SourceFile GetSourceFile(string filepath)
        {
            var reader = ReaderCreator.CreateReader(filepath);
            reader.GetSourceFile(); 
            return reader.SourceFile;

        }
    }
}