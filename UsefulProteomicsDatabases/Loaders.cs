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

using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using System.Security.Cryptography;
using System.Text.RegularExpressions;
using System.Xml.Serialization;
using UsefulProteomicsDatabases.Generated;

namespace UsefulProteomicsDatabases
{
    public static class Loaders
    {
        static Loaders()
        {
            ServicePointManager.SecurityProtocol = SecurityProtocolType.Tls12 | SecurityProtocolType.Tls11 | SecurityProtocolType.Tls;
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

        public static void UpdateElements(string elementLocation)
        {
            DownloadElements(elementLocation);
            if (!File.Exists(elementLocation))
            {
                Console.WriteLine("Element database did not exist, writing to disk");
                File.Move(elementLocation + ".temp", elementLocation);
                return;
            }
            if (FilesAreEqual_Hash(elementLocation + ".temp", elementLocation))
            {
                Console.WriteLine("Element database is up to date, doing nothing");
                File.Delete(elementLocation + ".temp");
            }
            else
            {
                Console.WriteLine("Element database updated, saving old version as backup");
                File.Move(elementLocation, elementLocation + DateTime.Now.ToString("dd-MMM-yyyy-HH-mm-ss"));
                File.Move(elementLocation + ".temp", elementLocation);
            }
        }

        public static Dictionary<string, int> GetFormalChargesDictionary(obo psiModDeserialized)
        {
            var modsWithFormalCharges = psiModDeserialized.Items.OfType<UsefulProteomicsDatabases.Generated.oboTerm>().Where(b => b.xref_analog != null && b.xref_analog.Any(c => c.dbname.Equals("FormalCharge")));
            Regex digitsOnly = new Regex(@"[^\d]");
            return modsWithFormalCharges.ToDictionary(b => "PSI-MOD; " + b.id, b => int.Parse(digitsOnly.Replace(b.xref_analog.First(c => c.dbname.Equals("FormalCharge")).name, "")));
        }

        public static void LoadElements(string elementLocation)
        {
            if (!File.Exists(elementLocation))
            {
                UpdateElements(elementLocation);
            }
            PeriodicTableLoader.Load(elementLocation);
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
            return psimodSerializer.Deserialize(new FileStream(psimodLocation, FileMode.Open, FileAccess.Read, FileShare.Read)) as Generated.obo;
        }

        public static IEnumerable<Modification> LoadUniprot(string uniprotLocation, Dictionary<string, int> formalChargesDictionary)
        {
            if (!File.Exists(uniprotLocation))
            {
                UpdateUniprot(uniprotLocation);
            }
            return PtmListLoader.ReadModsFromFile(uniprotLocation, formalChargesDictionary, out var errors).OfType<Modification>();
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
            using (WebClient Client = new WebClient())
            {
                Client.DownloadFile(@"https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/master/PSI-MOD.obo.xml", psimodLocation + ".temp");
            }
        }

        private static void DownloadUnimod(string unimodLocation)
        {
            using (WebClient Client = new WebClient())
            {
                Client.DownloadFile(@"http://www.unimod.org/xml/unimod.xml", unimodLocation + ".temp");
            }
        }

        private static void DownloadElements(string elementLocation)
        {
            using (WebClient Client = new WebClient())
            {
                Client.DownloadFile(@"http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some", elementLocation + ".temp");
            }
        }

        private static void DownloadUniprot(string uniprotLocation)
        {
            using (WebClient Client = new WebClient())
            {
                Client.DownloadFile(@"http://www.uniprot.org/docs/ptmlist.txt", uniprotLocation + ".temp");
            }
        }
    }
}