using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;

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
        public static string RetrieveProteome(string proteomeID, string absolutePathToStorageDirectory, ProteomeFormat format, 
            Reviewed reviewed, Compress compress, IncludeIsoforms include)
        {
            if (Directory.Exists(absolutePathToStorageDirectory))
            {
                string htmlQueryString = "";
                string filename = "\\" + proteomeID;
                bool compressBool = false; 
                bool isoformBool = false; 
                bool reviewedBool = false; 
                if (format == ProteomeFormat.fasta)
                {
                    if (reviewed == Reviewed.yes)
                    {
                        filename += "_reviewed";
                        reviewedBool = true; 
                    }
                    else
                    {
                        filename += "_unreviewed";
                    }
                    //only fasta style proteome allows retrieval of extra isoforms
                    if (include == IncludeIsoforms.yes)
                    {
                        filename += "_isoform";
                        isoformBool = true;
                    }
                    filename += ".fasta";
                    if (compress == Compress.yes)
                    {
                        filename += ".gz";
                        compressBool = true;
                    }

                    htmlQueryString = "https://rest.uniprot.org/uniprot/search?query=" + proteomeID + "+AND+" + "reviewed:" + reviewedBool.ToString().ToLower() + 
                        "&compressed=" + compressBool.ToString().ToLower() + "&format=" + format + "&includeIsoforms:" + isoformBool.ToString().ToLower();

                }
                else if (format == ProteomeFormat.xml)
                {
                    if (reviewed == Reviewed.yes)
                    {
                        filename += "_reviewed";
                        reviewedBool = true; 
                    }
                    else
                    {
                        filename += "_unreviewed";
                    }
                    filename += ".xml";
                    if (compress == Compress.yes)
                    {
                        filename += ".gz";
                        compressBool = true; 
                    }
                    htmlQueryString = "https://rest.uniprot.org/proteome/search?query=" + proteomeID + "+AND+reviewed:" + reviewedBool.ToString().ToLower()
                        + "&compressed=" + compressBool.ToString().ToLower() + "&format=" + format;

                }
                if (htmlQueryString.Length > 0)
                {
                    Loaders.DownloadContent(htmlQueryString, absolutePathToStorageDirectory + filename);
                    return absolutePathToStorageDirectory + filename;
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
        /// <param name="destinationFolder">filepath to the downloaded filefilepath</param>
        /// <returns></returns>
        public static string DownloadAvailableUniProtProteomes(string destinationFolder)
        {
            if (Directory.Exists(destinationFolder))
            {   
                string htmlQueryString = "https://rest.uniprot.org/proteomes/search?query=*&format=tsv&compressed=true";

                string filename = "availableUniProtProteomes.txt.gz";

                string filepath = Path.Combine(destinationFolder, filename);
                Loaders.DownloadContent(htmlQueryString, filepath);

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
                Dictionary<string, string> dictionaryOfAvailableProteomes = new();
                string fileExtension = Path.GetExtension(completePathToAvailableUniProtProteomes);
                switch (fileExtension)
                {
                    case ".gz":
                        foreach (string item in ReadAllGZippedLines(completePathToAvailableUniProtProteomes).ToList())
                        {
                            var lineValuesArray = item.Split("\t");
                            dictionaryOfAvailableProteomes.Add(lineValuesArray[0], lineValuesArray[1]);
                        }
                        return dictionaryOfAvailableProteomes;

                    case ".zip":
                        foreach (string item in ReadAllZippedLines(completePathToAvailableUniProtProteomes).ToList())
                        {
                            var lineValuesArray = item.Split("\t");
                            dictionaryOfAvailableProteomes.Add(lineValuesArray[0], lineValuesArray[1]);
                        }
                        return dictionaryOfAvailableProteomes;

                    case ".txt":
                        foreach (string item in File.ReadAllLines(completePathToAvailableUniProtProteomes).ToList())
                        {
                            var lineValuesArray = item.Split("\t");
                            dictionaryOfAvailableProteomes.Add(lineValuesArray[0], lineValuesArray[1]);
                        }
                        return dictionaryOfAvailableProteomes;

                    default:
                        return null; //no file with the appropriate extension is present
                };
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
            Dictionary<string, string> d = new();
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
        /// https://www.uniprot.org/help/return_fields
        /// </summary>
        public enum Columns
        {
        }

        public static IEnumerable<string> ReadAllGZippedLines(string filename)
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

        public static IEnumerable<string> ReadAllZippedLines(string filename)
        {
            using (ZipArchive archive = ZipFile.OpenRead(filename))
            {
                foreach (ZipArchiveEntry entry in archive.Entries)
                {
                    using (var reader = new StreamReader(entry.Open()))
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