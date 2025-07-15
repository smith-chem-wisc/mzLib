using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.IO.Compression;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Xml;
using Chemistry;
using Transcriptomics;
using Omics.BioPolymer;

namespace UsefulProteomicsDatabases.Transcriptomics
{
    public enum RnaFastaHeaderType
    {
        Modomics,
        Ensembl,
        NcbiRefSeq,
        NcbiAssembly,
        Unknown,
    }

    public enum SequenceTransformationOnRead
    {
        None,
        ConvertAllTtoU
    }

    public static class RnaDbLoader
    {

        #region Header Detection and Property Regexes

        public static string[] RnaFastaExtensions = [".fasta", ".fna", ".fa"];

        /// <summary>
        /// Checks for XX_YYYY where x is a capital letter and y is a number. 
        /// </summary>
        private static readonly Regex _ncbiRefSeqHeaderRegex = new Regex(@"\b[A-Z]{2}_[0-9]{4}\b");
        private static readonly Regex SubstituteWhitespace = new Regex(@"\s+");
        private static readonly Regex _ncbiAssemblyHeaderRegex = new Regex(@"^>NM_\d+\.\d+ ", RegexOptions.Compiled);
        private static readonly Regex _ncbiRefSeqGeneHeaderRegex = new Regex(@"^>NC_\d+\.\d+:", RegexOptions.Compiled);


        public static RnaFastaHeaderType DetectRnaFastaHeaderType(string line)
        {
            if (line.StartsWith(">id"))
                return RnaFastaHeaderType.Modomics;
            if (line.StartsWith(">ENST"))
                return RnaFastaHeaderType.Ensembl;
            if (_ncbiAssemblyHeaderRegex.IsMatch(line))
                return RnaFastaHeaderType.NcbiAssembly;
            if (_ncbiRefSeqGeneHeaderRegex.IsMatch(line))
                return RnaFastaHeaderType.NcbiRefSeq;

            return RnaFastaHeaderType.Unknown;
        }

        /// <summary>
        /// Dictionary that extract accession number, species, name, and additional dataField of modomics
        /// </summary>
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

        public static readonly Dictionary<string, FastaHeaderFieldRegex> EnsemblFieldRegexes =
            new()
            {
                { "Accession", new FastaHeaderFieldRegex("Accession", @"^>(ENST\d+\.\d+)", 0, 1) },
                { "Name", new FastaHeaderFieldRegex("Name", @"scaffold:([^\s]+)", 0, 1) },
                { "Organism", new FastaHeaderFieldRegex("Organism", @"scaffold:([^\s:]+):([^\s:]+):", 0, 1) }, 
                { "Gene", new FastaHeaderFieldRegex("Gene", @"gene:(ENSG\d+\.\d+)", 0, 1) },
                { "GeneBiotype", new FastaHeaderFieldRegex("GeneBiotype", @"gene_biotype:([^\s]+)", 0, 1) },
                { "TranscriptBiotype", new FastaHeaderFieldRegex("TranscriptBiotype", @"transcript_biotype:([^\s]+)", 0, 1) },
                { "GeneSymbol", new FastaHeaderFieldRegex("GeneSymbol", @"gene_symbol:([^\s]+)", 0, 1) },
                { "Description", new FastaHeaderFieldRegex("Description", @"description:([^\[]+)", 0, 1) },
            };

        public static readonly Dictionary<string, FastaHeaderFieldRegex> NcbiAssemblyFieldRegexes =
            new()
            {
                { "Accession", new FastaHeaderFieldRegex("Accession", @"^>(NM_\d+\.\d+)", 0, 1) },
                { "Name", new FastaHeaderFieldRegex("Name", @"^>NM_\d+\.\d+ ([^(]+)", 0, 1) },
                // Gene: everything after "Homo sapiens" and before the first comma
                { "Gene", new FastaHeaderFieldRegex("Gene", @"^>NM_\d+\.\d+ ([A-Za-z]+ [A-Za-z]+) ([^,]+)", 0, 2) },
                {
                    "Organism", new FastaHeaderFieldRegex("Organism", @"^>NM_\d+\.\d+ (([A-Za-z]+ ){1}[A-Za-z]+)", 0, 1)
                },
            };

        public static readonly Dictionary<string, FastaHeaderFieldRegex> NcbiRefSeqGeneFieldRegexes =
            new()
            {
                { "Accession", new FastaHeaderFieldRegex("Accession", @"^>(NC_\d+\.\d+)", 0, 1) },
                { "Name", new FastaHeaderFieldRegex("Name", @"^\S+:\d+-\d+ ([^\[]+)", 0, 1) },
                { "Organism", new FastaHeaderFieldRegex("Organism", @"\[organism=([^\]]+)\]", 0, 1) },
                { "Gene", new FastaHeaderFieldRegex("Gene", @"\[GeneID=(\d+)\]", 0, 1) },
                { "Chromosome", new FastaHeaderFieldRegex("Chromosome", @"\[chromosome=([^\]]+)\]", 0, 1) },
            };

        #endregion

        /// <summary>
        /// Loads an RNA file from the specified location, optionally generating decoys and adding error tracking
        /// </summary>
        /// <param name="rnaDbLocation">The file path to the RNA FASTA database</param>
        /// <param name="generateTargets">Flag indicating whether to generate targets or not</param>
        /// <param name="decoyType">The type of decoy generation to apply</param>
        /// <param name="isContaminant">Indicates if the RNA sequence is a contaminant</param>
        /// <param name="errors">Outputs any errors encountered during the process</param>
        /// <param name="fivePrimeTerm">An optional 5' prime chemical modification term</param>
        /// <param name="threePrimeTerm">An optional 3' prime chemical modification term</param>
        /// <returns>A list of RNA sequences loaded from the FASTA database</returns>
        /// <exception cref="MzLibUtil.MzLibException">Thrown if the FASTA header format is unknown or other issues occur during loading.</exception>

        public static List<RNA> LoadRnaFasta(string rnaDbLocation, bool generateTargets, DecoyType decoyType,
            bool isContaminant, out List<string> errors, IHasChemicalFormula? fivePrimeTerm = null, IHasChemicalFormula? threePrimeTerm = null, 
            int maxThreads = 1, string decoyIdentifier = "DECOY")
        {
            RnaFastaHeaderType? headerType = null;
            SequenceTransformationOnRead sequenceTransformation = SequenceTransformationOnRead.None;
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
                            headerType = DetectRnaFastaHeaderType(line);

                            switch (headerType)
                            {
                                case RnaFastaHeaderType.Modomics:
                                    regexes = ModomicsFieldRegexes;
                                    identifierHeader = "SOterm";
                                    break;
                                case RnaFastaHeaderType.Ensembl:
                                    regexes = EnsemblFieldRegexes;
                                    identifierHeader = "Accession";
                                    sequenceTransformation = SequenceTransformationOnRead.ConvertAllTtoU;
                                    break;
                                case RnaFastaHeaderType.NcbiAssembly:
                                    regexes = NcbiAssemblyFieldRegexes;
                                    identifierHeader = "Accession";
                                    sequenceTransformation = SequenceTransformationOnRead.ConvertAllTtoU;
                                    break;
                                case RnaFastaHeaderType.NcbiRefSeq:
                                    regexes = NcbiRefSeqGeneFieldRegexes;
                                    identifierHeader = "Accession";
                                    sequenceTransformation = SequenceTransformationOnRead.ConvertAllTtoU;
                                    break;
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
                    else
                    {
                        sb?.Append(line.Trim());
                    }

                    if ((fasta.Peek() == '>' || fasta.Peek() == -1) /*&& accession != null*/ && sb != null)
                    {
                        Dictionary<string, string> additonalDatabaseFields =
                            regexResults.ToDictionary(x => x.Key, x => x.Value);

                        List<Tuple<string, string>> geneNames = null!;
                        if (regexResults.ContainsKey("Gene"))
                        {
                            string geneName = regexResults["Gene"];
                            regexResults.Remove("Gene");
                            if (!string.IsNullOrEmpty(geneName))
                            {
                                geneNames = new List<Tuple<string, string>> { new Tuple<string, string>(geneName, geneName) };
                            }
                        }

                        var sequence = SanitizeAndTransform(sb.ToString(), sequenceTransformation);

                        RNA rna = new RNA(sequence, identifier,
                            null, fivePrimeTerminus: fivePrimeTerm, threePrimeTerminus: threePrimeTerm, name: name, organism: organism, databaseFilePath: rnaDbLocation, isContaminant: isContaminant, isDecoy: false, geneNames: geneNames, databaseAdditionalFields: additonalDatabaseFields);
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

            List<RNA> decoys = RnaDecoyGenerator.GenerateDecoys(targets, decoyType, maxThreads, decoyIdentifier);
            return generateTargets ? targets.Concat(decoys).ToList() : decoys;
        }


        private static Dictionary<string, string> ParseRegexFields(string line,
            Dictionary<string, FastaHeaderFieldRegex> regexes)
        {
            Dictionary<string, string> fields = new Dictionary<string, string>();

            foreach (var regex in regexes)
            {
                string match = regex.Value.ApplyRegex(line);
                fields.Add(regex.Key, match);
            }

            return fields;
        }

        public static Dictionary<string, IList<Modification>> IdToPossibleMods = new Dictionary<string, IList<Modification>>();
        public static Dictionary<string, Modification> IdWithMotifToMod = new Dictionary<string, Modification>();

        public static List<RNA> LoadRnaXML(string rnaDbLocation, bool generateTargets, DecoyType decoyType,
            bool isContaminant, IEnumerable<Modification> allKnownModifications,
            IEnumerable<string> modTypesToExclude, out Dictionary<string, Modification> unknownModifications,
            int maxHeterozygousVariants = 4, int minAlleleDepth = 1,
            int maxThreads = 1, IHasChemicalFormula? fivePrimeTerm = null, IHasChemicalFormula? threePrimeTerm = null,
            string decoyIdentifier = "DECOY")
        {
            var prespecified = ProteinDbLoader.GetPtmListFromProteinXml(rnaDbLocation);
            allKnownModifications = allKnownModifications ?? new List<Modification>();
            modTypesToExclude = modTypesToExclude ?? new List<string>();

            if (prespecified.Count > 0 || allKnownModifications.Count() > 0)
            {
                //modsDictionary = GetModificationDict(new HashSet<Modification>(prespecified.Concat(allKnownModifications)));
                IdToPossibleMods = ProteinDbLoader.GetModificationDict(new HashSet<Modification>(prespecified.Concat(allKnownModifications)));
                IdWithMotifToMod = ProteinDbLoader.GetModificationDictWithMotifs(new HashSet<Modification>(prespecified.Concat(allKnownModifications)));
            }
            List<RNA> targets = new List<RNA>();
            unknownModifications = new Dictionary<string, Modification>();

            string newProteinDbLocation = rnaDbLocation;

            //we had trouble decompressing and streaming on the fly so we decompress completely first, then stream the file, then delete the decompressed file
            if (rnaDbLocation.EndsWith(".gz"))
            {
                newProteinDbLocation = Path.Combine(Path.GetDirectoryName(rnaDbLocation), "temp.xml");
                using var stream = new FileStream(rnaDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read);
                using FileStream outputFileStream = File.Create(newProteinDbLocation);
                using var decompressor = new GZipStream(stream, CompressionMode.Decompress);
                decompressor.CopyTo(outputFileStream);
            }

            using (var uniprotXmlFileStream = new FileStream(newProteinDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                Regex substituteWhitespace = new Regex(@"\s+");

                ProteinXmlEntry block = new ProteinXmlEntry();

                using (XmlReader xml = XmlReader.Create(uniprotXmlFileStream))
                {
                    while (xml.Read())
                    {
                        if (xml.NodeType == XmlNodeType.Element)
                        {
                            block.ParseElement(xml.Name, xml);
                        }
                        if (xml.NodeType == XmlNodeType.EndElement || xml.IsEmptyElement)
                        {
                            RNA newProtein = block.ParseRnaEndElement(xml, modTypesToExclude, unknownModifications, isContaminant, rnaDbLocation);
                            if (newProtein != null)
                            {
                                targets.Add(newProtein);
                            }
                        }
                    }
                }
            }
            if (newProteinDbLocation != rnaDbLocation)
            {
                File.Delete(newProteinDbLocation);
            }

            List<RNA> decoys = RnaDecoyGenerator.GenerateDecoys(targets, decoyType, maxThreads, decoyIdentifier);
            IEnumerable<RNA> proteinsToExpand = generateTargets ? targets.Concat(decoys) : decoys;
            return proteinsToExpand.SelectMany(p => p.GetVariantBioPolymers(maxHeterozygousVariants, minAlleleDepth)).ToList();
        }


        // TODO: Some oligo databases may have the reverse strand, this is currently not handled yet and this code assumes we are always reading in the strand to search against. 
        public static string SanitizeAndTransform(string rawSequence, SequenceTransformationOnRead sequenceTransformation)
        {
            var cleanedSequence = SubstituteWhitespace.Replace(rawSequence, "");

            switch (sequenceTransformation)
            {
                case SequenceTransformationOnRead.ConvertAllTtoU:
                    return cleanedSequence.Replace('T', 'U');
                case SequenceTransformationOnRead.None:
                    return cleanedSequence;
                default:
                    throw new ArgumentOutOfRangeException(nameof(sequenceTransformation), sequenceTransformation, null);
            }
        }
    }
}
