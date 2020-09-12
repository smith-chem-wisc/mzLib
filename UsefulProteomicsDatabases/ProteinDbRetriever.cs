using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Net;

namespace UsefulProteomicsDatabases
{
    public static class ProteinDbRetriever
    {
        /// <summary>
        /// Constructor forms a UniProt search query capable of downloading information from UniProt for storage in a file
        /// It's primary function is to download a proteome, but it could have many other purposes.
        /// A poorly formed query will not return anything, but it won't crash.
        /// Successful search returns full path.
        /// Failed search returns null.
        /// </summary>
        /// <param name="proteomeID">valid proteome ID corresponding to a specific organism (e.g. UP000005640)</param>
        /// <param name="formt">format of retrieved proteome (e.g. fasta or xml)</param>
        /// <param name="reviewed">if yes file contains only reviewd proteins</param>
        /// <param name="compress">if yes file is saved as .gz</param>
        /// <param name="absolutePathToStorageDirectory"></param>
        public static string RetrieveProteome(string proteomeID, string absolutePathToStorageDirectory, ProteomeFormat format, Reviewed reviewed, Compress compress, IncludeIsoforms include)
        {
            if (Directory.Exists(absolutePathToStorageDirectory))
            {
                string htmlQueryString = "";
                string filename = "\\" + proteomeID;
                if (format == ProteomeFormat.fasta)
                {
                    if (reviewed == Reviewed.yes)
                    {
                        filename += "_reviewed";
                    }
                    else
                    {
                        filename += "_unreviewed";
                    }
                    //only fasta style proteome allows retrieval of extra isoforms
                    if (include == IncludeIsoforms.yes)
                    {
                        filename += "_isoform";
                    }
                    filename += ".fasta";
                    if (compress == Compress.yes)
                    {
                        filename += ".gz";
                    }
                    htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:" + reviewed + "&compress=" + compress + "&format=" + format + "&include:" + include;
                }
                else if (format == ProteomeFormat.xml)
                {
                    if (reviewed == Reviewed.yes)
                    {
                        filename += "_reviewed";
                    }
                    else
                    {
                        filename += "_unreviewed";
                    }
                    filename += ".xml";
                    if (compress == Compress.yes)
                    {
                        filename += ".gz";
                    }
                    htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:" + reviewed + "&compress=" + compress + "&format=" + format;
                }
                if (htmlQueryString.Length > 0)
                {
                    using (WebClient Client = new WebClient())
                    {
                        Client.DownloadFile(htmlQueryString, absolutePathToStorageDirectory + filename);
                        return absolutePathToStorageDirectory + filename;
                    }
                }
                //we don't support other file types yet.
                return null;
            }
            //a successful search will return the full path
            //a failed search will return null
            return null;
        }

        /// <summary>
        /// downloades and then returns the filepath to a compressed (.gz), tab-delimited text file of the available proteomes. Line one is the header.
        /// </summary>
        /// <param name="filepath">filepath to the downloaded filefilepath</param>
        /// <returns></returns>
        public static string DownloadAvailableUniProtProteomes(string filepath)
        {
            if (Directory.Exists(filepath))
            {
                string htmlQueryString = "https://www.uniprot.org/proteomes/?query=*&format=tab&compress=yes&columns=id,name,organism-id,proteincount,busco,cpd,assembly%20representation";
                string filename = "\\availableUniProtProteomes.txt.gz";

                filepath += filename;
                using (WebClient Client = new WebClient())
                {
                    Client.DownloadFile(htmlQueryString, filepath);
                }

                if (File.Exists(filepath))
                {
                    return filepath;
                }

                //the download was not successful
                else
                {
                    return null;
                }
            }

            //the filepath did not exist
            return null;
        }

        /// <summary>
        /// The list of available UniProt Proteomes was downloaded 09/03/2020
        /// This is returned as a dictionary
        /// key = Proteome ID
        /// value = Organism
        /// </summary>
        /// <returns></returns>
        public static Dictionary<string, string> UniprotProteomesList(string completePathToAvailableUniProtProteomes)
        {
            if (File.Exists(completePathToAvailableUniProtProteomes))
            {
                Dictionary<string, string> dictionaryOfAvailableProteomes = new Dictionary<string, string>();
                List<string> idNameList = new List<string>();
                if (Path.GetExtension(completePathToAvailableUniProtProteomes) == ".gz")
                {
                    idNameList = ReadAllZippedLines(completePathToAvailableUniProtProteomes).ToList();
                    foreach (string item in idNameList)
                    {
                        var lineValuesArray = item.Split("\t");
                        dictionaryOfAvailableProteomes.Add(lineValuesArray[0], lineValuesArray[1]);
                    }
                    return dictionaryOfAvailableProteomes;
                }
                return null;
            }
            //file does not exist
            return null;
        }

        /// <summary>
        /// One can retrieve a table of information about specific proteins.
        /// This method returns the names of columns that could be included in the table
        /// </summary>
        /// <returns></returns>
        public static Dictionary<string, string> UniprotColumnsList()
        {
            string currentDirectory = Directory.GetCurrentDirectory();
            string filePath = Path.Combine(currentDirectory, "UniProtKB_columnNamesForProgrammaticAccess.txt");
            Dictionary<string, string> d = new Dictionary<string, string>();
            string[] idNameList = File.ReadAllLines(filePath);
            foreach (string item in idNameList)
            {
                if (!item.StartsWith('#'))
                {
                    var j = item.Split("\t");
                    d.Add(j[0], j[1]);
                }
            }
            return d;
        }

        /// <summary>
        /// This is an enumerated list of a file types available for download at UniProt
        /// </summary>
        public enum ProteomeFormat
        {
            html,
            tab,
            xls,
            fasta,
            gff,
            txt,
            xml,
            rdf,
            list,
            rss
        }

        /// <summary>
        /// This designates whether or not unreviewed proteins are included in the downloaded proteome.
        /// </summary>
        public enum Reviewed
        {
            yes,
            no
        }

        /// <summary>
        /// Compress will cause the table to be downloaded in .gz format
        /// </summary>
        public enum Compress
        {
            yes,
            no
        }

        /// <summary>
        /// Include isoform sequences when the format parameter is set to fasta.
        /// Include description of referenced data when the format parameter is set to rdf.
        /// </summary>
        public enum IncludeIsoforms
        {
            yes,
            no
        }

        /// <summary>
        /// Columns to select for retrieving results in tab or xls format.
        /// https://www.uniprot.org/help/uniprotkb_column_names
        /// </summary>
        public enum Columns
        {
        }

        public static IEnumerable<string> ReadAllZippedLines(string filename)
        {
            using (var fileStream = File.OpenRead(filename))
            {
                using (var gzipStream = new GZipStream(fileStream, CompressionMode.Decompress))
                {
                    using (var reader = new StreamReader(gzipStream))
                    {
                        string currentLine;
                        while ((currentLine = reader.ReadLine()) != null)
                        {
                            yield return currentLine;
                        }
                    }
                }
            }
        }
    }
}