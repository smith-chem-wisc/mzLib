using CsvHelper;
using Omics.SpectrumMatch;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.SpectrumLibraries
{
    public abstract class SpectrumLibraryFile : ResultFile<LibrarySpectrum>, IResultFile
    {
        protected SpectrumLibraryFile(string filePath, Software software = Software.Unspecified) : base(filePath, software)
        {
            
        }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        protected internal SpectrumLibraryFile() : base()
        {

        }
        public void WriteMsp(string FilePath)
        {
            if (!CanRead(FilePath))
                FilePath += FileType.GetFileExtension();

            using var fs = new FileStream(FilePath, FileMode.Create, FileAccess.Write);

            //TODO: Implement MspSpectrumLibraryWriter
            //using var csv = new CsvWriter(new StreamWriter(File.Create(FilePath)), FlashDeconvTsv.CsvConfiguration);

            //csv.WriteHeader<FlashDeconvTsv>();
            //foreach (var result in Results)
            //{
            //    csv.NextRecord();
            //    csv.WriteRecord(result);
            //}
        }
        public void WriteMsFraggerSpecLibFile(string FilePath)
        {
            if (!CanRead(FilePath))
                FilePath += FileType.GetFileExtension();
            if (this is MsFraggerSpeclibFile)
            {
                //just write the file in the original fragger format with one fragment per line and no grouped spectra
                this.WriteResults(FilePath);
            }
            else
            {
                var sw = new StreamWriter(File.Create(FilePath));
                using (var csv = new CsvWriter(sw, MsFraggerSpeclib.CsvConfiguration))
                {
                    csv.WriteHeader<MsFraggerSpeclib>();
                    foreach (LibrarySpectrum result in Results)
                    {
                        //TODO: fix this
                        sw.WriteLine(result.ToFraggerLibraryString("", ""));
                    }

                }

            }
        }
    }
}
