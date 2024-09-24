using Easy.Common.Extensions;
using MzLibUtil;
using Omics.Digestion;

namespace Transcriptomics.Digestion
{
    public static class RnaseDictionary
    {
        static RnaseDictionary()
        {
            // TODO: Load dictionary automatically on first call to RnaseDictioanry


            var pathToProgramFiles = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
            string dataDirectory = !String.IsNullOrWhiteSpace(pathToProgramFiles) &&
                                   AppDomain.CurrentDomain.BaseDirectory.Contains(pathToProgramFiles)
                                   && !AppDomain.CurrentDomain.BaseDirectory.Contains("Jenkins")
                ? Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
                    "MetaMorpheus")
                : AppDomain.CurrentDomain.BaseDirectory;
            string path = Path.Combine(dataDirectory, "Digestion", "rnases.tsv");
            Dictionary = LoadRnaseDictionary(path);
        }

        public static Dictionary<string, Rnase> Dictionary { get; set; }


        public static Dictionary<string, Rnase> LoadRnaseDictionary(string path)
        {
            Dictionary<string, Rnase> dict = new();
            string[] myLines = File.ReadAllLines(path).Skip(1).ToArray();

            foreach (var line in File.ReadAllLines(path).Skip(1))
            {
                if (line.Trim() != string.Empty)
                {
                    string[] fields = line.Split('\t');
                    List<DigestionMotif> motifList = DigestionMotif.ParseDigestionMotifsFromString(fields[1]);
                    string name = fields[0];
                    CleavageSpecificity cleavage = Enum.Parse<CleavageSpecificity>(fields[4], true);

                    var rnase = new Rnase(name, cleavage, motifList);
                    if (!dict.ContainsKey(rnase.Name))
                    {
                        dict.Add(rnase.Name, rnase);
                    }
                    else
                    {
                        throw new MzLibException("More than one Rnase named {rnase.Name} exits");
                    }
                }
            }

            if (!Dictionary.IsNotNullOrEmpty())
                Dictionary = dict;

            return dict;
        }
    }
}