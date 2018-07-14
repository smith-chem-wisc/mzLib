using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Proteomics.ProteolyticDigestion
{
    public static class ProteaseDictionary
    {
        static ProteaseDictionary()
        {
            var pathToProgramFiles = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
            string dataDirectory = !String.IsNullOrWhiteSpace(pathToProgramFiles) && AppDomain.CurrentDomain.BaseDirectory.Contains(pathToProgramFiles)
                    && !AppDomain.CurrentDomain.BaseDirectory.Contains("Jenkins") ?
                Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "MetaMorpheus") :
                AppDomain.CurrentDomain.BaseDirectory;
            string p = Path.Combine(dataDirectory, "ProteolyticDigestion", "proteases.tsv");
            Dictionary = LoadProteaseDictionary(Path.Combine(dataDirectory, "ProteolyticDigestion", "proteases.tsv"));
        }

        public static Dictionary<string, Protease> Dictionary { get; private set; }

        public static Dictionary<string, Protease> LoadProteaseDictionary(string proteasesLocation)
        {
            Dictionary<string, Protease> dict = new Dictionary<string, Protease>();
            using (StreamReader proteases = new StreamReader(proteasesLocation))
            {
                proteases.ReadLine();

                while (proteases.Peek() != -1)
                {
                    string line = proteases.ReadLine();
                    string[][] fields = line.Split('\t').Select(x => x.Split('|')).ToArray();
                    string name = fields[0][0];
                    string[] preventing;
                    List<Tuple<string, TerminusType>> sequencesInducingCleavage = new List<Tuple<string, TerminusType>>();
                    List<Tuple<string, TerminusType>> sequencePreventingCleavage = new List<Tuple<string, TerminusType>>();
                    for (int i = 0; i < fields[1].Length; i++)
                    {
                        if (!fields[1][i].Equals(""))
                        {
                            sequencesInducingCleavage.Add(new Tuple<string, TerminusType>(fields[1][i], ((TerminusType)Enum.Parse(typeof(TerminusType), fields[3][i], true))));
                            if (!fields[2].Contains(""))
                            {
                                preventing = (fields[2][i].Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries));
                                for (int j = 0; j < preventing.Length; j++)
                                {
                                    sequencePreventingCleavage.Add(new Tuple<string, TerminusType>(preventing[j], (TerminusType)Enum.Parse(typeof(TerminusType), fields[3][i], true)));
                                }
                            }
                        }
                    }
                    var cleavageSpecificity = ((CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), fields[4][0], true));
                    string psiMsAccessionNumber = fields[5][0];
                    string psiMsName = fields[6][0];
                    string siteRegexp = fields[7][0];
                    var protease = new Protease(name, sequencesInducingCleavage, sequencePreventingCleavage, cleavageSpecificity, psiMsAccessionNumber, psiMsName, siteRegexp);
                    dict.Add(protease.Name, protease);
                }
            }
            return dict;
        }
    }
}