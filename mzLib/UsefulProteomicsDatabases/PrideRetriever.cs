using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace UsefulProteomicsDatabases
{
    public class PrideRetriever
    {
        private string b = "https://www.ebi.ac.uk/pride/ws/archive/v2/status/PXD048177";


        public static string RetrieveMassSpecProject(string PrideProjectAccession, string outputFullFilePath)
        {
            string htmlQueryString = "https://www.ebi.ac.uk/pride/ws/archive/v2/status/" + PrideProjectAccession;

            if (htmlQueryString.Length > 0)
            {
                Loaders.DownloadContent(htmlQueryString, outputFullFilePath);
                return outputFullFilePath;
            }

            return "";
        }
        public static string RetrieveProjectFileByFilename(string PrideProjectAccession, string filename, string outputFullFilePath)
        {
            string htmlQueryString =
                "https://www.ebi.ac.uk/pride/ws/archive/v2/files/fileByName?fileName=d_atg1_d_atg11_proteome_data_analysis.7z&projectAccession=PXD048176";

            if (htmlQueryString.Length > 0)
            {
                Loaders.DownloadContent(htmlQueryString, outputFullFilePath);
                return outputFullFilePath;
            }

            return "";
        }
    }
}
