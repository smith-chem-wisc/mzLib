using CsvHelper;
using Omics.SpectrumMatch;

namespace Readers.SpectrumLibraries
{
    public class MsFraggerSpeclibFile : SpectrumLibraryFile
    {
        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }

        public List<MsFraggerSpeclib> OriginalRecords { get; private set; } //single line in the .speclib file. combine ultiple lines to create a single library spectrum


        public MsFraggerSpeclibFile() : base()
        {
        }

        public MsFraggerSpeclibFile(string filePath) : base(filePath, Software.MsFragger)
        {
        }


        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), MsFraggerSpeclib.CsvConfiguration);
            OriginalRecords = csv.GetRecords<MsFraggerSpeclib>().ToList();
            List<LibrarySpectrum> librarySpectra = new List<LibrarySpectrum>();
            foreach (var spectrumFragmentGroup in OriginalRecords.GroupBy(p=>p))
            {
                //TODO populate the library spectrum from the spectrum fragment group
                LibrarySpectrum ls = null;
                librarySpectra.Add(ls);
            }
            Results = librarySpectra;
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using (var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), MsFraggerSpeclib.CsvConfiguration))
            {
                csv.WriteHeader<MsFraggerSpeclib>();
                foreach (var result in OriginalRecords)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }
        }
    }

}
