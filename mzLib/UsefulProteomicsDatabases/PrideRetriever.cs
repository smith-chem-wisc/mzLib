using Omics.Fragmentation;
using Omics.SpectrumMatch;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using System.Reflection.PortableExecutable;
using System.Text;
using System.Threading.Tasks;
using static System.Net.WebRequestMethods;
using Newtonsoft.Json;
using System.Net.Http;
using Newtonsoft.Json.Linq;
using System.ComponentModel;
using NUnit.Framework;

namespace UsefulProteomicsDatabases
{
    public class PrideRetriever
    {
        public static List<PRIDEEntry> PrideEntries;
        public static List<string> PrideFileNames;
        public static PrideProject PrideProject;
        private const string ApiBaseUrl = "https://www.ebi.ac.uk/pride/ws/archive/v2/";

        public PrideRetriever (string PrideProjectAccession, string outputFullFileFolder)
        {
            PrideEntries = RetrieveMassSpecProject( PrideProjectAccession, outputFullFileFolder + "/prideEntries.json");
            PrideProject = RetrievePrideProjectExperimentalInformation(PrideProjectAccession, outputFullFileFolder + "/PrideProject.json");
            PrideFileNames = PrideEntries.Select(p => p.FileName).ToList();
        }

        public List<PRIDEEntry> GetPrideEntries()
        {
            return PrideEntries;
        }

        public List<string> GetPrideFileNames()
        {
            return PrideFileNames;
        }

        public PrideProject GetPrideProjectInformation()
        {
            return PrideProject;
        }

        public  List<PRIDEEntry> RetrieveMassSpecProject(string PrideProjectAccession, string outputFullFilePath)
        {
            string htmlQueryString = ApiBaseUrl+ "files/byProject?accession=" + PrideProjectAccession;

            if (htmlQueryString.Length > 0)
            {
                Loaders.DownloadContent(htmlQueryString, outputFullFilePath);
                using (StreamReader r = new StreamReader(outputFullFilePath))
                {
                    string json = r.ReadToEnd();
                    PrideEntries = JsonConvert.DeserializeObject<List<PRIDEEntry>>(json);
                }
                return PrideEntries;
            }

            return null;
        }

        public static PrideProject RetrievePrideProjectExperimentalInformation(string PrideProjectAccession, string outputFullFilePath)
        {
            string htmlQueryString = ApiBaseUrl+"projects/" + PrideProjectAccession;

            if (htmlQueryString.Length > 0)
            {
                Loaders.DownloadContent(htmlQueryString, outputFullFilePath);
                using (StreamReader r = new StreamReader(outputFullFilePath))
                {
                    string json = r.ReadToEnd();
                    PrideProject = JsonConvert.DeserializeObject<PrideProject>(json);
                }
                return PrideProject;
            }

            return null;
        }

        public string GetFtpLinkByFileNmae(string prideFilename)
        {
            PRIDEEntry targetEntry = GetPrideEntries().Where(p=>p.FileName == prideFilename).ToList().First();
            List<PublicFileLocation> publicFileLocations = targetEntry.PublicFileLocations;
            String FTPlink = publicFileLocations.Where(p=>p.Name.Contains("FTP", StringComparison.InvariantCultureIgnoreCase)).ToList().First().Value;
            return FTPlink;
        }

        public string RetrieveProjectFileByFilename( string prideFilename, string outputFileDirectory)
        {
            string FTPlink = GetFtpLinkByFileNmae(prideFilename);
            if (FTPlink.Length > 0)
            {
                string fullFilePath = outputFileDirectory +"\\" +prideFilename;
                PrideRetriever.SimpleWebClientDownload(FTPlink, fullFilePath);
                return fullFilePath;
            }
            return "";
        }

        public List<string> GetAllFilesByAccession(string PrideProjectAccession, string outputFileDirectory)
        {
            List<string> fullFilePathList = new List<string>();
            foreach ( var prideFilename in PrideFileNames)
            {
                string FTPlink = GetFtpLinkByFileNmae(prideFilename);
                if (FTPlink.Length > 0)
                {
                    string fullFilePath = outputFileDirectory + "\\" + prideFilename;
                    PrideRetriever.SimpleWebClientDownload(FTPlink, fullFilePath);
                    fullFilePathList.Add(fullFilePath);
                }
            }
            return fullFilePathList;
        }

        public static List<string> SearchByKeywordsAndFilters(string keyword, string queryFilter, int pageSize, int page, string dateGap, string sortDirection, string sortFields)
        {
            List<string> accessionList = new List<string>();
            var requestUrl = $"{ApiBaseUrl}search/projects?keyword={keyword}&pageSize={pageSize}&page={page}&sortDirection={sortDirection}&sortFields={sortFields}";
            if (!string.IsNullOrEmpty(queryFilter))
            {
                requestUrl += $"&filter={queryFilter}";
            }

            if (!string.IsNullOrEmpty(dateGap))
            {
                requestUrl += $"&dateGap={dateGap}";
            }

            var outputFullFileFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData");
            var outputDir = Path.Combine(outputFullFileFolder, @"PrideRetrieveSearchBykeys");
            Directory.CreateDirectory(outputDir);

            string outputFullFilePath = Path.Combine(outputDir, @"searchBykeys.json");

            Loaders.DownloadContent(requestUrl, outputFullFilePath);

            using (StreamReader r = new StreamReader(outputFullFilePath))
            {
                string json = r.ReadToEnd();
                string jsonContent = System.IO.File.ReadAllText(outputFullFilePath);
                var jsonData = JObject.Parse(jsonContent);
                accessionList = jsonData["_embedded"]["compactprojects"].Select(p=>p["accession"].ToString()).ToList();
            }
            Directory.Delete(outputDir, true);
            return accessionList;
        }

        public static void PrideFtp()
        {
            // Get the object used to communicate with the server.
            FtpWebRequest request = (FtpWebRequest)WebRequest.Create("ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2023/12/PXD048176/d_atg1_d_atg11_proteome_data_analysis.7z");
            request.Method = WebRequestMethods.Ftp.DownloadFile;

            // This example assumes the FTP site uses anonymous logon.
            request.Credentials = new NetworkCredential("anonymous", "anonymous");

            request.KeepAlive = true;
            request.UsePassive = true; 
            request.UseBinary = true;

            //// Read the file from the server & write to destination                
            //using (FtpWebResponse response = (FtpWebResponse)request.GetResponse()) // Error here
            //using (Stream responseStream = response.GetResponseStream())
            //using (StreamReader reader = new StreamReader(responseStream))
            //using (StreamWriter destination = new StreamWriter(destinationFile))
            //{
            //    destination.Write(reader.ReadToEnd());
            //    destination.Flush();
            //}

        }

        public static void SimpleWebClientDownload(string url, string fullFilePath)
        {
            WebClient client = new WebClient();
            client.Credentials = new NetworkCredential("anonymous", "anonymous");
            client.DownloadFile(
                url, fullFilePath);
        }
    }
}
