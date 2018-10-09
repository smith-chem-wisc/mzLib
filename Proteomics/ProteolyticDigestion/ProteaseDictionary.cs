using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Proteomics.Fragmentation;

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

            string path = Path.Combine(dataDirectory, "ProteolyticDigestion", "proteases.tsv");          
            Dictionary = LoadProteaseDictionary(path);
        
        }

        public static Dictionary<string, Protease> Dictionary { get; private set; }

        public static Dictionary<string, Protease> LoadProteaseDictionary(string path)
        {

            Dictionary<string, Protease> dict = new Dictionary<string, Protease>();

            string[] myLines = File.ReadAllLines(path);
            myLines = myLines.Skip(1).ToArray();

            foreach (string line in myLines)
            {
                if (line.Trim() != string.Empty) // skip empty lines
                {
                    string[] fields = line.Split('\t');
                    List<DigestionMotif> motifList = DigestionMotif.ParseProteaseFromString(fields[1]);

                    string name = fields[0];
                    var cleavageSpecificity = ((CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), fields[4], true));
                    string psiMsAccessionNumber = fields[5];
                    string psiMsName = fields[6];
                    var protease = new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList);
                    dict.Add(protease.Name, protease);
                }
            }

            return dict;

        }
    }
}