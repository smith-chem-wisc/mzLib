using Omics.SpectrumMatch;

namespace Readers.SpectrumLibraries
{
    internal class MspSpectrumLibraryReader
    {
        public static List<LibrarySpectrum> ReadMsp(string filePath, out List<string> warnings)
        {
            warnings = new List<string>();
            return new List<LibrarySpectrum>();
        }
    }
}
