using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using MzLibUtil;
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

        public static Dictionary<string, Protease> Dictionary { get; set; }

        public static Dictionary<string, Protease> LoadProteaseDictionary(string path, List<Modification> proteaseMods = null)
        {

            Dictionary<string, Protease> dict = new Dictionary<string, Protease>();

            string[] myLines = File.ReadAllLines(path);
            myLines = myLines.Skip(1).ToArray();

            foreach (string line in myLines)
            {
                if (line.Trim() != string.Empty) // skip empty lines
                {
                    string[] fields = line.Split('\t');
                    List<DigestionMotif> motifList = DigestionMotif.ParseDigestionMotifsFromString(fields[1]);
                    string name = fields[0];
                    var cleavageSpecificity = ((CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), fields[4], true));
                    string psiMsAccessionNumber = fields[5];
                    string psiMsName = fields[6];
                    //name of the modification that is associated with proteolytic cleavage
                    string proteaseModDetails = fields[8];  
                    //if this protease has an associated modification, look it up in the list of mods loaded fro the protease mods file
                    if (proteaseModDetails != "" && proteaseMods != null)
                    {
                        if (proteaseMods.Select(p => p.IdWithMotif).ToList().Contains(proteaseModDetails))
                        {
                            Modification proteaseModification = proteaseMods.Where(p => p.IdWithMotif == proteaseModDetails).First();
                            var protease = new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList, proteaseModification);
                            if (!dict.ContainsKey(protease.Name))
                            {
                                dict.Add(protease.Name, protease);
                            }
                            else 
                            {
                                throw new MzLibException("More than one protease named "+ protease.Name +" exists");
                            }
                            
                        }
                        else 
                        {                            
                            var protease = new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList);
                            if (!dict.ContainsKey(protease.Name))
                            {
                                dict.Add(protease.Name, protease);
                            }
                            else
                            {
                                throw new MzLibException("More than one protease named " + protease.Name + " exists");
                            }
                            throw new MzLibException(proteaseModDetails + " is not a valid modification");
                        }
                        
                    }
                    else
                    {
                        var protease = new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList);
                        if (!dict.ContainsKey(protease.Name))
                        {
                            dict.Add(protease.Name, protease);
                        }
                        else
                        {
                            throw new MzLibException("More than one protease named " + protease.Name + " exists");
                        }
                    }
                    
                }
            }

            return dict;

        }

        
    }
}