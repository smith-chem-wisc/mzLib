using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
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
        public static string RetrieveProjectFileByFilename(string prideProjectAccession, string prideFilename, string outputFileDirectory)
        {
            string htmlQueryString =
                "https://www.ebi.ac.uk/pride/ws/archive/v2/files/fileByName?fileName=" + prideFilename + "&projectAccession=" + prideProjectAccession;

            if (htmlQueryString.Length > 0)
            {
                Loaders.DownloadContent(htmlQueryString, Path.Join(outputFileDirectory,prideFilename));
                return outputFileDirectory;
            }

            return "";
        }

        public static void PrideFtp()
        {
            // Get the object used to communicate with the server.
            FtpWebRequest request = (FtpWebRequest)WebRequest.Create("ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2023/12/PXD048176/d_atg1_d_atg11_proteome_data_analysis.7z");
            request.Method = WebRequestMethods.Ftp.DownloadFile;

            // This example assumes the FTP site uses anonymous logon.
            request.Credentials = new NetworkCredential("anonymous", "anonymous");

            FtpWebResponse response = (FtpWebResponse)request.GetResponse();

            Stream responseStream = response.GetResponseStream();
            StreamReader reader = new StreamReader(responseStream);

            var j = reader.ReadToEnd();
            var r = response.StatusDescription;

            //Console.WriteLine(reader.ReadToEnd());

            //Console.WriteLine($"Download Complete, status {response.StatusDescription}");

            reader.Close();
            response.Close();
        }
    }
}
