using Easy.Common.Extensions;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml;

namespace UsefulProteomicsDatabases
{
    public enum FastaHeaderType
    { UniProt, Ensembl, Gencode, Unknown }

    public static class ProteinDbLoader
    {
        public static readonly FastaHeaderFieldRegex UniprotAccessionRegex = new FastaHeaderFieldRegex("accession", @"[|](.+)[|]", 0, 1);
        public static readonly FastaHeaderFieldRegex UniprotFullNameRegex = new FastaHeaderFieldRegex("fullName", @"\s(.*?)\s(OS=|GN=|PE=|SV=|OX=)", 0, 1);
        public static readonly FastaHeaderFieldRegex UniprotNameRegex = new FastaHeaderFieldRegex("name", @"\|(?:.+)\|(.*?)(\s|$)", 0, 1);
        public static readonly FastaHeaderFieldRegex UniprotGeneNameRegex = new FastaHeaderFieldRegex("geneName", @"GN=(.*?)(\s|$)", 0, 1);
        public static readonly FastaHeaderFieldRegex UniprotOrganismRegex = new FastaHeaderFieldRegex("organism", @"OS=(.*?)\s(GN=|PE=|SV=|OX=)", 0, 1);

        public static readonly FastaHeaderFieldRegex EnsemblAccessionRegex = new FastaHeaderFieldRegex("accession", @"([A-Z0-9_.]+)", 0, 1);
        public static readonly FastaHeaderFieldRegex EnsemblFullNameRegex = new FastaHeaderFieldRegex("fullName", @"(pep:.*)", 0, 1);
        public static readonly FastaHeaderFieldRegex EnsemblGeneNameRegex = new FastaHeaderFieldRegex("geneName", @"gene:([^ ]+)", 0, 1);

        public static readonly FastaHeaderFieldRegex GencodeAccessionRegex = new FastaHeaderFieldRegex("accession", @">(.*?)\|", 0, 1);
        public static readonly FastaHeaderFieldRegex GencodeFullNameRegex = new FastaHeaderFieldRegex("fullName", @">(?:.*?)\|(?:.*?)\|(?:.*?)\|(?:.*?)\|(?:.*?)\|(.*?)\|", 0, 1);
        public static readonly FastaHeaderFieldRegex GencodeGeneNameRegex = new FastaHeaderFieldRegex("geneName", @">(?:.*?)\|(?:.*?)\|(?:.*?)\|(?:.*?)\|(?:.*?)\|(?:.*?)\|(.*?)\|", 0, 1);

        public static Dictionary<string, IList<Modification>> IdToPossibleMods = new Dictionary<string, IList<Modification>>();
        public static Dictionary<string, Modification> IdWithMotifToMod = new Dictionary<string, Modification>();

        private static Regex UnicodeRegex = new Regex(@"[^\u0000-\u007F]+");

        /// <summary>
        /// Stores the last database file path.
        /// </summary>
        private static string last_database_location;

        /// <summary>
        /// Stores the modification list read during LoadProteinXML
        /// </summary>
        private static List<Modification> protein_xml_modlist_general = new List<Modification>();

        /// <summary>
        /// Load a mzLibProteinDb or UniProt XML file. Protein modifications may be specified before the protein entries (mzLibProteinDb format).
        /// If so, this modification list can be acquired with GetPtmListFromProteinXml after using this method.
        /// They may also be read in separately from a ptmlist text file, and then input as allKnownModifications.
        /// If protein modifications are specified both in the mzLibProteinDb XML file and in allKnownModifications, they are collapsed into a HashSet of Modifications before generating Protein entries.
        /// </summary>
        [SuppressMessage("Microsoft.Usage", "CA2202:Do not dispose objects multiple times")]
        public static List<Protein> LoadProteinXML(string proteinDbLocation, bool generateTargets, DecoyType decoyType, IEnumerable<Modification> allKnownModifications,
            bool isContaminant, IEnumerable<string> modTypesToExclude, out Dictionary<string, Modification> unknownModifications, int maxThreads = -1,
            int maxHeterozygousVariants = 4, int minAlleleDepth = 1, bool addTruncations = false)
        {
            List<Modification> prespecified = GetPtmListFromProteinXml(proteinDbLocation);
            allKnownModifications = allKnownModifications ?? new List<Modification>();
            modTypesToExclude = modTypesToExclude ?? new List<string>();

            //Dictionary<string, IList<Modification>> modsDictionary = new Dictionary<string, IList<Modification>>();
            if (prespecified.Count > 0 || allKnownModifications.Count() > 0)
            {
                //modsDictionary = GetModificationDict(new HashSet<Modification>(prespecified.Concat(allKnownModifications)));
                IdToPossibleMods = GetModificationDict(new HashSet<Modification>(prespecified.Concat(allKnownModifications)));
                IdWithMotifToMod = GetModificationDictWithMotifs(new HashSet<Modification>(prespecified.Concat(allKnownModifications)));
            }
            List<Protein> targets = new List<Protein>();
            unknownModifications = new Dictionary<string, Modification>();

            string newProteinDbLocation = proteinDbLocation;

            //we had trouble decompressing and streaming on the fly so we decompress completely first, then stream the file, then delete the decompressed file
            if (proteinDbLocation.EndsWith(".gz"))
            {
                newProteinDbLocation = Path.Combine(Path.GetDirectoryName(proteinDbLocation),"temp.xml");
                using var stream = new FileStream(proteinDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read);
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
                            Protein newProtein = block.ParseEndElement(xml, modTypesToExclude, unknownModifications, isContaminant, proteinDbLocation);
                            if (newProtein != null)
                            {
                                if (addTruncations)
                                {
                                    newProtein.AddTruncations();
                                }
                                targets.Add(newProtein);
                            }
                        }
                    }
                    
                }

            }

            if (newProteinDbLocation != proteinDbLocation)
            {
                File.Delete(newProteinDbLocation);
            }

            List<Protein> decoys = DecoyProteinGenerator.GenerateDecoys(targets, decoyType, maxThreads);
            IEnumerable<Protein> proteinsToExpand = generateTargets ? targets.Concat(decoys) : decoys;
            return proteinsToExpand.SelectMany(p => p.GetVariantProteins(maxHeterozygousVariants, minAlleleDepth)).ToList();
        }

        /// <summary>
        /// Get the modification entries specified in a mzLibProteinDb XML file (.xml or .xml.gz).
        /// </summary>
        [SuppressMessage("Microsoft.Usage", "CA2202:Do not dispose objects multiple times")]
        public static List<Modification> GetPtmListFromProteinXml(string proteinDbLocation)
        {
            if (proteinDbLocation.Equals(last_database_location))
            {
                return protein_xml_modlist_general;
            }
            else
            {
                last_database_location = proteinDbLocation;
                StringBuilder storedKnownModificationsBuilder = new StringBuilder();
                using (var stream = new FileStream(proteinDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    Regex startingWhitespace = new Regex(@"/^\s+/gm");
                    Stream uniprotXmlFileStream = proteinDbLocation.EndsWith(".gz") ?
                        (Stream)(new GZipStream(stream, CompressionMode.Decompress)) :
                        stream;

                    using (XmlReader xml = XmlReader.Create(uniprotXmlFileStream))
                    {
                        while (xml.Read())
                        {
                            if (xml.NodeType == XmlNodeType.Element)
                            {
                                if (xml.Name == "modification")
                                {
                                    string modification = startingWhitespace.Replace(xml.ReadElementString(), "");
                                    storedKnownModificationsBuilder.AppendLine(modification);
                                }
                                else if (xml.Name == "entry")
                                {
                                    //if we are up to entry fields in the protein database, then there no more prespecified modifications to read and we
                                    //can begin processing all the lines we have read.
                                    //This block of code does not process information in any of the entries.
                                    protein_xml_modlist_general = storedKnownModificationsBuilder.Length <= 0 ?
                                        new List<Modification>() :
                                        PtmListLoader.ReadModsFromString(storedKnownModificationsBuilder.ToString(), out var errors).ToList();
                                    break;
                                }
                            }
                        }
                    }
                    return protein_xml_modlist_general;
                }
            }
        }

        /// <summary>
        /// Load a protein fasta database, using regular expressions to get various aspects of the headers. The first regex capture group is used as each field.
        /// </summary>
        public static List<Protein> LoadProteinFasta(string proteinDbLocation, bool generateTargets, DecoyType decoyType, bool isContaminant, out List<string> errors,
            FastaHeaderFieldRegex accessionRegex = null, FastaHeaderFieldRegex fullNameRegex = null, FastaHeaderFieldRegex nameRegex = null,
            FastaHeaderFieldRegex geneNameRegex = null, FastaHeaderFieldRegex organismRegex = null, int maxThreads = -1, bool addTruncations = false)
        {
            FastaHeaderType? HeaderType = null;
            HashSet<string> unique_accessions = new HashSet<string>();
            int unique_identifier = 2; //for isoforms. the first will be "accession", the next will be "accession_2"
            string accession = null;
            string name = null;
            string fullName = null;
            string organism = null;
            List<Tuple<string, string>> geneName = new List<Tuple<string, string>>();
            errors = new List<string>();
            Regex substituteWhitespace = new Regex(@"\s+");

            List<Protein> targets = new List<Protein>();

            string newProteinDbLocation = proteinDbLocation;

            //we had trouble decompressing and streaming on the fly so we decompress completely first, then stream the file, then delete the decompressed file
            if (proteinDbLocation.EndsWith(".gz"))
            {
                newProteinDbLocation = Path.Combine(Path.GetDirectoryName(proteinDbLocation), "temp.fasta");
                using var stream = new FileStream(proteinDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read);
                using FileStream outputFileStream = File.Create(newProteinDbLocation);
                using var decompressor = new GZipStream(stream, CompressionMode.Decompress);
                decompressor.CopyTo(outputFileStream);
            }

            using (var fastaFileStream = new FileStream(newProteinDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                StringBuilder sb = null;
                StreamReader fasta = new StreamReader(fastaFileStream);

                while (true)
                {
                    string line = "";
                    line = fasta.ReadLine();
                    if (line == null) { break; }

                    if (line.StartsWith(">"))
                    {
                        // determine the header type if necessary (this figures out how to parse the different fields in the fasta headers)
                        if (HeaderType == null && accessionRegex == null)
                        {
                            HeaderType = DetectFastaHeaderFormat(line);

                            if (HeaderType == FastaHeaderType.Unknown)
                            {
                                throw new MzLibUtil.MzLibException("Unknown fasta header format: " + line);
                            }

                            switch (HeaderType)
                            {
                                case FastaHeaderType.UniProt:
                                    accessionRegex = UniprotAccessionRegex;
                                    fullNameRegex = UniprotFullNameRegex;
                                    nameRegex = UniprotNameRegex;
                                    organismRegex = UniprotOrganismRegex;
                                    geneNameRegex = UniprotGeneNameRegex;
                                    break;

                                case FastaHeaderType.Ensembl:
                                    accessionRegex = EnsemblAccessionRegex;
                                    fullNameRegex = EnsemblFullNameRegex;
                                    geneNameRegex = EnsemblGeneNameRegex;
                                    break;

                                case FastaHeaderType.Gencode:
                                    accessionRegex = GencodeAccessionRegex;
                                    fullNameRegex = GencodeFullNameRegex;
                                    geneNameRegex = GencodeGeneNameRegex;
                                    break;
                            }
                        }

                        accession = ApplyRegex(accessionRegex, line);
                        fullName = ApplyRegex(fullNameRegex, line);
                        name = ApplyRegex(nameRegex, line);
                        organism = ApplyRegex(organismRegex, line);
                        string geneNameString = ApplyRegex(geneNameRegex, line);
                        if (geneNameString != null)
                        {
                            geneName.Add(new Tuple<string, string>("primary", geneNameString));
                        }

                        if (accession == null || accession == "")
                        {
                            accession = line.Substring(1).TrimEnd();
                        }

                        sb = new StringBuilder();
                    }
                    else if (sb != null)
                    {
                        sb.Append(line.Trim());
                    }

                    if ((fasta.Peek() == '>' || fasta.Peek() == -1) && accession != null && sb != null)
                    {
                        string sequence = substituteWhitespace.Replace(sb.ToString(), "");

                        // sanitize the sequence to replace unexpected characters with X (unknown amino acid)
                        // sometimes strange characters get added by RNA sequencing software, etc.
                        sequence = SanitizeAminoAcidSequence(sequence, 'X');

                        if (unique_accessions.Contains(accession)) //this will happen for isoforms
                        {
                            string originalAccession = accession; //save the original
                            accession += "_" + unique_identifier.ToString(); //add a number onto it
                            while (unique_accessions.Contains(accession)) //if that number was already added
                            {
                                unique_identifier++; //keep increasing it
                                accession = originalAccession + "_" + unique_identifier.ToString(); //try the new number
                            }
                            unique_identifier = 2; //reset
                        }
                        unique_accessions.Add(accession);
                        Protein protein = new Protein(sequence, accession, organism, geneName, name: name, fullName: fullName,
                            isContaminant: isContaminant, databaseFilePath: proteinDbLocation, addTruncations: addTruncations);
                        if (protein.Length == 0)
                        {
                            errors.Add("Line" + line + ", Protein Length of 0: " + protein.Name + " was skipped from database: " + proteinDbLocation);
                        }
                        else
                        {
                            targets.Add(protein);
                        }

                        accession = null;
                        name = null;
                        fullName = null;
                        organism = null;
                        geneName = new List<Tuple<string, string>>();
                    }

                    // no input left
                    if (fasta.Peek() == -1)
                    {
                        break;
                    }
                }
            }

            if (newProteinDbLocation != proteinDbLocation)
            {
                File.Delete(newProteinDbLocation);
            }

            if (!targets.Any())
            {
                errors.Add("Error: No proteins could be read from the database: " + proteinDbLocation);
            }
            List<Protein> decoys = DecoyProteinGenerator.GenerateDecoys(targets, decoyType, maxThreads);
            return generateTargets ? targets.Concat(decoys).ToList() : decoys;
        }

        /// <summary>
        /// Merge proteins that have the same accession, sequence, and contaminant designation.
        /// </summary>
        public static IEnumerable<Protein> MergeProteins(IEnumerable<Protein> mergeThese)
        {
            Dictionary<Tuple<string, string, bool, bool>, List<Protein>> proteinsByAccessionSequenceContaminant = new Dictionary<Tuple<string, string, bool, bool>, List<Protein>>();
            foreach (Protein p in mergeThese)
            {
                Tuple<string, string, bool, bool> key = new Tuple<string, string, bool, bool>(p.Accession, p.BaseSequence, p.IsContaminant, p.IsDecoy);
                if (!proteinsByAccessionSequenceContaminant.TryGetValue(key, out List<Protein> bundled))
                {
                    proteinsByAccessionSequenceContaminant.Add(key, new List<Protein> { p });
                }
                else
                {
                    bundled.Add(p);
                }
            }

            foreach (KeyValuePair<Tuple<string, string, bool, bool>, List<Protein>> proteins in proteinsByAccessionSequenceContaminant)
            {
                HashSet<string> names = new HashSet<string>(proteins.Value.Select(p => p.Name));
                HashSet<string> fullnames = new HashSet<string>(proteins.Value.Select(p => p.FullName));
                HashSet<string> descriptions = new HashSet<string>(proteins.Value.Select(p => p.FullDescription));
                HashSet<Tuple<string, string>> genenames = new HashSet<Tuple<string, string>>(proteins.Value.SelectMany(p => p.GeneNames));
                HashSet<ProteolysisProduct> proteolysis = new HashSet<ProteolysisProduct>(proteins.Value.SelectMany(p => p.ProteolysisProducts));
                HashSet<SequenceVariation> variants = new HashSet<SequenceVariation>(proteins.Value.SelectMany(p => p.SequenceVariations));
                HashSet<DatabaseReference> references = new HashSet<DatabaseReference>(proteins.Value.SelectMany(p => p.DatabaseReferences));
                HashSet<DisulfideBond> bonds = new HashSet<DisulfideBond>(proteins.Value.SelectMany(p => p.DisulfideBonds));
                HashSet<SpliceSite> splices = new HashSet<SpliceSite>(proteins.Value.SelectMany(p => p.SpliceSites));

                Dictionary<int, HashSet<Modification>> mod_dict = new Dictionary<int, HashSet<Modification>>();
                foreach (KeyValuePair<int, List<Modification>> nice in proteins.Value.SelectMany(p => p.OneBasedPossibleLocalizedModifications).ToList())
                {
                    if (!mod_dict.TryGetValue(nice.Key, out HashSet<Modification> val))
                    {
                        mod_dict.Add(nice.Key, new HashSet<Modification>(nice.Value));
                    }
                    else
                    {
                        foreach (Modification mod in nice.Value)
                        {
                            val.Add(mod);
                        }
                    }
                }
                Dictionary<int, List<Modification>> mod_dict2 = mod_dict.ToDictionary(kv => kv.Key, kv => kv.Value.ToList());

                yield return new Protein(
                    proteins.Key.Item2,
                    proteins.Key.Item1,
                    isContaminant: proteins.Key.Item3,
                    isDecoy: proteins.Key.Item4,
                    geneNames: genenames.ToList(),
                    oneBasedModifications: mod_dict2,
                    proteolysisProducts: proteolysis.ToList(),
                    name: names.FirstOrDefault(),
                    fullName: fullnames.FirstOrDefault(),
                    databaseReferences: references.ToList(),
                    disulfideBonds: bonds.ToList(),
                    sequenceVariations: variants.ToList(),
                    spliceSites: splices.ToList()
                    );
            }
        }

        private static string ApplyRegex(FastaHeaderFieldRegex regex, string line)
        {
            string result = null;
            if (regex != null)
            {
                var matches = regex.Regex.Matches(line);
                if (matches.Count > regex.Match && matches[regex.Match].Groups.Count > regex.Group)
                {
                    result = matches[regex.Match].Groups[regex.Group].Value;
                }
            }
            return result;
        }

        private static Dictionary<string, IList<Modification>> GetModificationDict(IEnumerable<Modification> mods)
        {
            var mod_dict = new Dictionary<string, IList<Modification>>();

            foreach (Modification mod in mods.Where(m => m.ValidModification))
            {
                string modIdWithoutMotif = mod.OriginalId.Split(new[] { " on " }, StringSplitOptions.None).First();
                if (mod_dict.TryGetValue(modIdWithoutMotif, out IList<Modification> val))
                {
                    val.Add(mod);
                }
                else
                {
                    mod_dict.Add(modIdWithoutMotif, new List<Modification> { mod });
                }
            }

            return mod_dict;
        }

        private static Dictionary<string, Modification> GetModificationDictWithMotifs(IEnumerable<Modification> mods)
        {
            var mod_dict = new Dictionary<string, Modification>();

            foreach (Modification mod in mods.Where(m => m.ValidModification))
            {
                if (!mod_dict.ContainsKey(mod.IdWithMotif))
                {
                    mod_dict.Add(mod.IdWithMotif, mod);
                }
            }

            return mod_dict;
        }

        public static string SanitizeAminoAcidSequence(string originalSequence, char replacementCharacter)
        {
            string cleaned = UnicodeRegex.Replace(originalSequence, replacementCharacter.ToString());

            for (int r = 0; r < cleaned.Length; r++)
            {
                if (!Residue.TryGetResidue(cleaned[r], out Residue res))
                {
                    cleaned = cleaned.Replace(cleaned[r], replacementCharacter);
                }
            }

            return cleaned;
        }

        public static FastaHeaderType DetectFastaHeaderFormat(string line)
        {
            if (!line.StartsWith(">"))
            {
                return FastaHeaderType.Unknown;
            }

            // check for uniprot header format
            // uniprot will look like ">sp|ACCESSION ....."
            // i.e., ">", followed by a two-letter lowercase code, then "|"
            if (line.Length >= 4 && char.IsLower(line[1]) && char.IsLower(line[2]) && line[3] == '|')
            {
                return FastaHeaderType.UniProt;
            }

            // check for ensembl header format
            // ensembl will have fields specified by a field name and then a colon, then a field value
            if (line.Contains("gene_symbol:", StringComparison.InvariantCultureIgnoreCase)
                || line.Contains("gene:", StringComparison.InvariantCultureIgnoreCase)
                || line.Contains("Acc:", StringComparison.InvariantCultureIgnoreCase)
                || line.Contains("description:", StringComparison.InvariantCultureIgnoreCase))
            {
                return FastaHeaderType.Ensembl;
            }

            // check for gencode header format
            // gencode will have 8 fields, each separated by "|".
            // missing fields are still present but have a "-" character as a placeholder.
            if (line.Count(p => p == '|') == 7)
            {
                return FastaHeaderType.Gencode;
            }

            return FastaHeaderType.Unknown;
        }
    }
}