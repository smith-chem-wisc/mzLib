using System;
using System.Collections.Generic;
using System.Globalization;
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
                    List<DigestionMotif> motifList = DigestionMotif.ParseDigestionMotifsFromString(fields[1]);

                    string name = fields[0];
                    var cleavageSpecificity = ((CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), fields[4], true));
                    string psiMsAccessionNumber = fields[5];
                    string psiMsName = fields[6];
                    string massShiftString = fields[8];
                    Dictionary<string, double> massShifts = new Dictionary<string, double>();
                    var numberFormat = new NumberFormatInfo();
                    numberFormat.NegativeSign = "-";
                    numberFormat.NumberDecimalSeparator = ".";
                    string proteaseModDetails = fields[8];
                    bool addH_Nterminal = true;
                    bool addOH_Cterminal = true;
                    if (proteaseModDetails != "" && proteaseModDetails.Contains(':'))
                    {
                        var tempArray = proteaseModDetails.Split(':');
                        var modDetails = tempArray[0];
                        var waterDetails = tempArray[1];
                        Modification mod = Modification.ParseProteaseModificationsFromString(modDetails);
                        var H2Obools = waterDetails.Split(',');
                        if (H2Obools[0].ToLower().Contains("true"))
                        {
                            addH_Nterminal = true;
                        }
                        if (H2Obools[1].ToLower().Contains("true"))
                        {
                            addOH_Cterminal = true;
                        }
                        if (H2Obools[0].ToLower().Contains("false"))
                        {
                            addH_Nterminal = false;
                        }
                        if (H2Obools[1].ToLower().Contains("false"))
                        {
                            addOH_Cterminal = false;
                        }                        
                        Tuple<Modification, bool, bool> proteaseModificationDetails = new Tuple<Modification, bool, bool>(mod, addH_Nterminal, addOH_Cterminal);
                        var protease = new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList, proteaseModificationDetails);
                        dict.Add(protease.Name, protease);
                    }
                    else
                    {
                        var protease = new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList);
                        dict.Add(protease.Name, protease);
                    }
                    
                }
            }

            return dict;

        }

        
    }
}