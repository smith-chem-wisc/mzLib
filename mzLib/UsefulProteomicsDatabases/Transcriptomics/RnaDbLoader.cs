#nullable enable
using System;
using System.Collections.Generic;
using System.IO.Compression;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Transcriptomics;
using Easy.Common;
using System.Text.RegularExpressions;
using Chemistry;
using Easy.Common.Extensions;

namespace UsefulProteomicsDatabases.Transcriptomics
{
    public enum RnaFastaHeaderType
    {
        Modomics,
        Unknown,
    }

    public static class RnaDbLoader
    {

        public static readonly Dictionary<string, FastaHeaderFieldRegex> ModomicsFieldRegexes =
            new Dictionary<string, FastaHeaderFieldRegex>()
            {
                { "Id", new FastaHeaderFieldRegex("Id", @"id:(?<id>.+?)\|", 0, 1) },
                { "Name", new FastaHeaderFieldRegex("Name", @"Name:(?<Name>.+?)\|", 0, 1) },
                { "SOterm", new FastaHeaderFieldRegex("SOterm", @"SOterm:(?<SOterm>.+?)\|", 0, 1) },
                { "Type", new FastaHeaderFieldRegex("Type", @"Type:(?<Type>.+?)\|", 0, 1) },
                { "Subtype", new FastaHeaderFieldRegex("Subtype", @"Subtype:(?<Subtype>.+?)\|", 0, 1) },
                { "Feature", new FastaHeaderFieldRegex("Feature", @"Feature:(?<Feature>.+?)\|", 0, 1) },
                { "Organism", new FastaHeaderFieldRegex("Organism", @"Species:(?<Species>.+?)$", 0, 1) },
                { "Cellular Localization", new FastaHeaderFieldRegex("CellularLocalization", @"Cellular_Localization:(?<Cellular_Localization>.+?)\|", 0, 1) },
            };


        public static List<RNA> LoadRnaFasta(string rnaDbLocation, bool generateTargets, DecoyType decoyType,
            bool isContaminant, out List<string> errors, IHasChemicalFormula? fivePrimeTerm = null, IHasChemicalFormula? threePrimeTerm = null)
        {
            RnaFastaHeaderType? headerType = null;
            Regex substituteWhitespace = new Regex(@"\s+");
            errors = new List<string>();
            List<RNA> targets = new List<RNA>();
            string identifierHeader = null;

            string name = null;
            string organism = null;
            string identifier = null;

            string newDbLocation = rnaDbLocation;

            //we had trouble decompressing and streaming on the fly so we decompress completely first, then stream the file, then delete the decompressed file
            if (rnaDbLocation.EndsWith(".gz"))
            {
                newDbLocation = Path.Combine(Path.GetDirectoryName(rnaDbLocation), "temp.fasta");
                using var stream = new FileStream(rnaDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read);
                using FileStream outputFileStream = File.Create(newDbLocation);
                using var decompressor = new GZipStream(stream, CompressionMode.Decompress);
                decompressor.CopyTo(outputFileStream);
            }

            using (var fastaFileStream = new FileStream(newDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                StringBuilder sb = null;
                StreamReader fasta = new StreamReader(fastaFileStream);
                Dictionary<string, string> regexResults = new();
                Dictionary<string, FastaHeaderFieldRegex> regexes = null;

                while (true)
                {
                    string line = "";
                    line = fasta.ReadLine();
                    if (line == null) { break; }

                    if (line.StartsWith(">"))
                    {
                        if (headerType is null)
                        {
                            headerType = DetectFastaHeaderType(line);

                            switch (headerType)
                            {
                                case RnaFastaHeaderType.Modomics:
                                    regexes = ModomicsFieldRegexes;
                                    identifierHeader = "SOterm";
                                    break;

                                case RnaFastaHeaderType.Unknown:
                                case null:
                                default:
                                    throw new MzLibUtil.MzLibException("Unknown fasta header format: " + line);
                            }
                        }


                        regexResults = ParseRegexFields(line, regexes);
                        name = regexResults["Name"];
                        regexResults.Remove("Name");
                        organism = regexResults["Organism"];
                        regexResults.Remove("Organism");
                        identifier = regexResults[identifierHeader];
                        regexResults.Remove(identifierHeader);

                        sb = new StringBuilder();
                    }
                    else if (sb is not null)
                    {
                        sb.Append(line.Trim());
                    }

                    if ((fasta.Peek() == '>' || fasta.Peek() == -1) /*&& accession != null*/ && sb != null)
                    {
                        string sequence = substituteWhitespace.Replace(sb.ToString(), "");
                        Dictionary<string, string> additonalDatabaseFields =
                            regexResults.ToDictionary(x => x.Key, x => x.Value);

                        // Do we need to sanitize the sequence? 

                        RNA rna = new RNA(sequence, name, identifier, organism, rnaDbLocation,
                            fivePrimeTerm, threePrimeTerm, null,
                            isContaminant, false, additonalDatabaseFields );
                        if (rna.Length == 0)
                            errors.Add("Line" + line + ", Rna length of 0: " + rna.Name + "was skipped from database: " + rnaDbLocation);
                        else
                            targets.Add(rna);

                        name = null;
                        organism = null;
                        identifier = null;
                        regexResults.Clear();
                    }

                    // no input left
                    if (fasta.Peek() == -1)
                    {
                        break;
                    }
                }
            }

            if (newDbLocation != rnaDbLocation)
                File.Delete(newDbLocation);

            if (!targets.Any())
                errors.Add("No targets were loaded from database: " + rnaDbLocation);
            
            List<RNA> decoys = RnaDecoyGenerator.GenerateDecoys(targets, decoyType);
            return generateTargets ? targets.Concat(decoys).ToList() : decoys;
        }

        private static RnaFastaHeaderType DetectFastaHeaderType(string line)
        {
            if (!line.StartsWith(">"))
                return RnaFastaHeaderType.Unknown;

            // modomics -> >id:1|Name:tdbR00000010|SOterm:SO:0000254

            return RnaFastaHeaderType.Modomics;
        }

        private static Dictionary<string, string> ParseRegexFields(string line,
            Dictionary<string, FastaHeaderFieldRegex> regexes)
        {
            Dictionary<string, string> fields = new Dictionary<string, string>();

            foreach (var regex in regexes)
            {
                string match = ProteinDbLoader.ApplyRegex(regex.Value, line);
                fields.Add(regex.Key, match);
            }

            return fields;
        }
    }
}
