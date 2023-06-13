// Copyright 2016 Stefan Solntsev
//
// This file (Loaders.cs) is part of UsefulProteomicsDatabases.
//
// UsefulProteomicsDatabases is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// UsefulProteomicsDatabases is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with UsefulProteomicsDatabases. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Security.Cryptography;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using System.Xml.Serialization;
using UsefulProteomicsDatabases.Generated;
using TopDownProteomics.IO.Obo;

namespace UsefulProteomicsDatabases
{
    public static class Loaders
    {
        static Loaders()
        {
            ServicePointManager.SecurityProtocol = SecurityProtocolType.Tls12 | SecurityProtocolType.Tls11 | SecurityProtocolType.Tls;
            SetCultureAsInvariantCulture();
        }

        public static void SetCultureAsInvariantCulture()
        {
            CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;
            CultureInfo.DefaultThreadCurrentUICulture = CultureInfo.InvariantCulture;
            Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
            Thread.CurrentThread.CurrentUICulture = CultureInfo.InvariantCulture;
        }

        public static void UpdateUniprot(string uniprotLocation)
        {
            DownloadUniprot(uniprotLocation);
            if (!File.Exists(uniprotLocation))
            {
                Console.WriteLine("Uniprot database did not exist, writing to disk");
                File.Move(uniprotLocation + ".temp", uniprotLocation);
                return;
            }
            bool ye = FilesAreEqual_Hash(uniprotLocation + ".temp", uniprotLocation);
            if (ye)
            {
                Console.WriteLine("Uniprot database is up to date, doing nothing");
                File.Delete(uniprotLocation + ".temp");
            }
            else
            {
                Console.WriteLine("Uniprot database updated, saving old version as backup");
                File.Move(uniprotLocation, uniprotLocation + DateTime.Now.ToString("dd-MMM-yyyy-HH-mm-ss"));
                File.Move(uniprotLocation + ".temp", uniprotLocation);
            }
        }

        public static void UpdateUnimod(string unimodLocation)
        {
            DownloadUnimod(unimodLocation);
            if (!File.Exists(unimodLocation))
            {
                Console.WriteLine("Unimod database did not exist, writing to disk");
                File.Move(unimodLocation + ".temp", unimodLocation);
                return;
            }
            bool ye = FilesAreEqual_Hash(unimodLocation + ".temp", unimodLocation);
            if (ye)
            {
                Console.WriteLine("Unimod database is up to date, doing nothing");
                File.Delete(unimodLocation + ".temp");
            }
            else
            {
                Console.WriteLine("Unimod database updated, saving old version as backup");
                File.Move(unimodLocation, unimodLocation + DateTime.Now.ToString("dd-MMM-yyyy-HH-mm-ss"));
                File.Move(unimodLocation + ".temp", unimodLocation);
            }
        }

        public static void UpdatePsiMod(string psimodLocation)
        {
            DownloadPsiMod(psimodLocation);
            if (!File.Exists(psimodLocation))
            {
                Console.WriteLine("PSI-MOD database did not exist, writing to disk");
                File.Move(psimodLocation + ".temp", psimodLocation);
                return;
            }
            if (FilesAreEqual_Hash(psimodLocation + ".temp", psimodLocation))
            {
                Console.WriteLine("PSI-MOD database is up to date, doing nothing");
                File.Delete(psimodLocation + ".temp");
            }
            else
            {
                Console.WriteLine("PSI-MOD database updated, saving old version as backup");
                File.Move(psimodLocation, psimodLocation + DateTime.Now.ToString("dd-MMM-yyyy-HH-mm-ss"));
                File.Move(psimodLocation + ".temp", psimodLocation);
            }
        }

        public static void UpdatePsiModObo(string psiModOboLocation)
        {
            DownloadPsiModObo(psiModOboLocation);
            if (!File.Exists(psiModOboLocation))
            {
                Console.WriteLine("psi-mod.obo database did not exist, writing to disk");
                File.Move(psiModOboLocation + ".temp", psiModOboLocation);
                return;
            }
            if (FilesAreEqual_Hash(psiModOboLocation + ".temp", psiModOboLocation))
            {
                Console.WriteLine("psi-mod.obo database is up to date, doing nothing");
                File.Delete(psiModOboLocation + ".temp");
            }
            else
            {
                Console.WriteLine("psi-mod.obo database updated, saving old version as backup");
                File.Move(psiModOboLocation, psiModOboLocation + DateTime.Now.ToString("dd-MMM-yyyy-HH-mm-ss"));
                File.Move(psiModOboLocation + ".temp", psiModOboLocation);
            }
        }

        public static IEnumerable<OboTerm> ReadPsiModFile(string psiModOboLocation)
        {
            OboParser oboParser = new();
            return oboParser.Parse(psiModOboLocation); 
        }

        public static Dictionary<string, int> GetFormalChargesDictionary(obo psiModDeserialized)
        {
            var modsWithFormalCharges = psiModDeserialized.Items.OfType<UsefulProteomicsDatabases.Generated.oboTerm>()
                .Where(b => b.xref_analog != null && b.xref_analog.Any(c => c.dbname.Equals("FormalCharge")));
            Regex digitsOnly = new(@"[^\d]");
            return modsWithFormalCharges.ToDictionary(b => "PSI-MOD; " + b.id, b => int.Parse(digitsOnly.Replace(b.xref_analog.First(c => c.dbname.Equals("FormalCharge")).name, "")));
        }

        public static string GetFormalChargeString(this OboTagValuePair tvPair)
        {
            return GetStringFromOboTagValuePairValue(tvPair, pattern: @"[\d](?:\+|\-)"); 
        }

        private static string GetStringFromOboTagValuePairValue(OboTagValuePair oboTerms, 
            string pattern)
        {
            var matchGroup = Regex.Match(oboTerms.Value, pattern); 
            return matchGroup.Groups[0].Value; 
        }

        public static Dictionary<string, int> GetFormalChargesDictionary(IEnumerable<OboTerm> terms)
        {
            var formalCharges =
                from i in terms
                from j in i.ValuePairs
                where j.Value.Contains("FormalCharge")
                select (i.Id, j.Value);
            Regex digitsOnly = new(@"[^\d]");
            var modifiedResults = formalCharges
                .Select
                (
                    id => ("PSI-MOD; " + id.Id,
                        int.Parse(GetChargeString(id.Value) + digitsOnly.Replace(id.Value, "")))
                );
            return modifiedResults.ToDictionary(i => i.Item1, i => i.Item2);
        }
        // This method is used to fix an issue where the polarity of the formal charge is not read correctly. 
        private static string GetChargeString(string entry)
        {
            var chargeMatch = Regex.Match(entry, @"(\+|\-)");
            return chargeMatch.Groups[1].Value;
        }

        public static void LoadElements()
        {
            // has the periodic table already been loaded?
            if (PeriodicTable.GetElement(1) != null)
            {
                return;
            }

            // periodic table has not been loaded yet - load it
            PeriodicTableLoader.Load();
        }

        public static IEnumerable<Modification> LoadUnimod(string unimodLocation)
        {
            if (!File.Exists(unimodLocation))
            {
                UpdateUnimod(unimodLocation);
            }
            return UnimodLoader.ReadMods(unimodLocation);
        }

        public static Generated.obo LoadPsiMod(string psimodLocation)
        {
            var psimodSerializer = new XmlSerializer(typeof(Generated.obo));

            if (!File.Exists(psimodLocation))
            {
                UpdatePsiMod(psimodLocation);
            }
            using (FileStream stream = new FileStream(psimodLocation, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                return psimodSerializer.Deserialize(stream) as Generated.obo;
            }
        }

        public static IEnumerable<Modification> LoadUniprot(string uniprotLocation, Dictionary<string, int> formalChargesDictionary)
        {
            if (!File.Exists(uniprotLocation))
            {
                UpdateUniprot(uniprotLocation);
            }
            return PtmListLoader.ReadModsFromFile(uniprotLocation, formalChargesDictionary, out var _).OfType<Modification>();
        }

        /// <summary>
        /// Retrieves data using async/await
        /// </summary>
        /// <param name="url">path to retrieve data from</param>
        /// <returns></returns>
        public static async Task<HttpResponseMessage> AwaitAsync_GetSomeData(string url)
        {
            var client = new HttpClient();
            var response = await client.GetAsync(url).ConfigureAwait(false);
            return response;
        }

        /// <summary>
        /// Downloads content from the web and saves it as a new file
        /// </summary>
        /// <param name="url">path to retrieve data from</param>
        /// <param name="outputFile">path to write data to</param>
        public static void DownloadContent(string url, string outputFile)
        {
            var httpResponseMessage = AwaitAsync_GetSomeData(url).Result;

            using (FileStream stream = new(outputFile, FileMode.CreateNew))
            {
                Task.Run(() => httpResponseMessage.Content.CopyToAsync(stream)).Wait();
            }
        }

        private static bool FilesAreEqual_Hash(string first, string second)
        {
            using (FileStream a = File.Open(first, FileMode.Open, FileAccess.Read))
            using (FileStream b = File.Open(second, FileMode.Open, FileAccess.Read))
            {
                byte[] firstHash = MD5.Create().ComputeHash(a);
                byte[] secondHash = MD5.Create().ComputeHash(b);
                for (int i = 0; i < firstHash.Length; i++)
                {
                    if (firstHash[i] != secondHash[i])
                    {
                        return false;
                    }
                }
                return true;
            }
        }

        private static void DownloadPsiMod(string psimodLocation)
        {
            DownloadContent(@"https://github.com/smith-chem-wisc/psi-mod-CV/blob/master/PSI-MOD.obo.xml?raw=true", psimodLocation + ".temp");
        }
        private static void DownloadPsiModObo(string psiModOboLocation)
        {
            DownloadContent(@"https://github.com/HUPO-PSI/psi-mod-CV/blob/master/PSI-MOD.obo?raw=true", psiModOboLocation + ".temp");
        }
        private static void DownloadUnimod(string unimodLocation)
        {
            DownloadContent(@"http://www.unimod.org/xml/unimod.xml", unimodLocation + ".temp");
        }

        private static void DownloadElements(string elementLocation)
        {
            DownloadContent(@"http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some", elementLocation + ".temp");
        }

        private static void DownloadUniprot(string uniprotLocation)
        {
            DownloadContent(@"http://uniprot.org/docs/ptmlist.txt", uniprotLocation + ".temp");
        }
    }
}