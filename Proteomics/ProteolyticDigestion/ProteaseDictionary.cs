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
                    if (massShiftString != "")
                    {
                        if (massShiftString.Contains(','))
                        {
                            var shiftStrings = massShiftString.Split(',');
                            foreach (var shift in shiftStrings)
                            {
                                if (shift.Contains(':'))
                                {
                                    var content = shift.Split(':');
                                    var residue = content[0];
                                    var massShift = double.Parse(content[1], numberFormat);
                                    if (!massShifts.ContainsKey(residue))
                                    {
                                        massShifts.Add(residue, massShift);
                                    }
                                }
                            }
                        }
                        else
                        {
                            if (massShiftString.Contains(':'))
                            {
                                var content = massShiftString.Split(':');
                                var residue = content[0];
                                var shift = double.Parse(content[1], numberFormat);
                                if (!massShifts.ContainsKey(residue))
                                {
                                    massShifts.Add(residue, shift);
                                }
                            }
                        }
                    }
                    var protease = new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList, massShifts);
                    dict.Add(protease.Name, protease);
                }
            }

            return dict;

        }
    }
}