using Proteomics;
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
    public static class ProteinDbLoader
    {
        #region Public Fields

        public static readonly FastaHeaderFieldRegex UniprotAccessionRegex = new FastaHeaderFieldRegex("accession", @"([A-Z0-9_.]+)", 0, 1);

        public static readonly FastaHeaderFieldRegex UniprotFullNameRegex = new FastaHeaderFieldRegex("fullName", @"\s(.*?)\sOS=", 0, 1);

        public static readonly FastaHeaderFieldRegex UniprotNameRegex = new FastaHeaderFieldRegex("name", @"\|([^\|][A-Z0-9_]+)", 1, 1);

        public static readonly FastaHeaderFieldRegex UniprotGeneNameRegex = new FastaHeaderFieldRegex("geneName", @"GN=([^ ]+)", 0, 1);

        public static readonly FastaHeaderFieldRegex UniprotOrganismRegex = new FastaHeaderFieldRegex("organism", @"OS=(.*?)\sGN=", 0, 1);

        public static readonly FastaHeaderFieldRegex EnsemblAccessionRegex = new FastaHeaderFieldRegex("accession", @"([A-Z0-9_.]+)", 0, 1);

        public static readonly FastaHeaderFieldRegex EnsemblFullNameRegex = new FastaHeaderFieldRegex("fullName", @"(pep:.*)", 0, 1);

        public static readonly FastaHeaderFieldRegex EnsemblGeneNameRegex = new FastaHeaderFieldRegex("geneName", @"gene:([^ ]+)", 0, 1);

        #endregion Public Fields

        #region Private Fields

        /// <summary>
        /// Stores the last database file path.
        /// </summary>
        private static string last_database_location;

        /// <summary>
        /// Stores the modification list read during LoadProteinXML
        /// </summary>
        private static List<Modification> protein_xml_modlist;

        #endregion Private Fields

        #region Public Methods

        /// <summary>
        /// Load a mzLibProteinDb or UniProt XML file. Protein modifications may be specified before the protein entries (mzLibProteinDb format).
        /// If so, this modification list can be acquired with GetPtmListFromProteinXml after using this method.
        /// They may also be read in separately from a ptmlist text file, and then input as allKnownModifications.
        /// If protein modifications are specified both in the mzLibProteinDb XML file and in allKnownModifications, they are collapsed into a HashSet of Modifications before generating Protein entries.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="proteinDbLocation"></param>
        /// <param name="generateTargetProteins"></param>
        /// <param name="generateDecoyProteins"></param>
        /// <param name="allKnownModifications"></param>
        /// <param name="isContaminant"></param>
        /// <param name="dbRefTypesToKeep"></param>
        /// <param name="modTypesToExclude"></param>
        /// <param name="unknownModifications"></param>
        /// <returns></returns>
        [SuppressMessage("Microsoft.Usage", "CA2202:Do not dispose objects multiple times")]
        public static List<Protein> LoadProteinXML(string proteinDbLocation, bool generateTargetProteins, DecoyType decoyType, IEnumerable<Modification> allKnownModifications,
            bool isContaminant, IEnumerable<string> modTypesToExclude, out Dictionary<string, Modification> unknownModifications)
        {
            List<Modification> prespecified = GetPtmListFromProteinXml(proteinDbLocation);

            Dictionary<string, IList<Modification>> mod_dict = new Dictionary<string, IList<Modification>>();
            if (prespecified.Count > 0 || allKnownModifications.Count() > 0)
            {
                mod_dict = GetModificationDict(new HashSet<Modification>(prespecified.Concat(allKnownModifications)));
            }

            List<Protein> result = new List<Protein>();
            unknownModifications = new Dictionary<string, Modification>();
            using (var stream = new FileStream(proteinDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                Regex substituteWhitespace = new Regex(@"\s+");

                Stream uniprotXmlFileStream = proteinDbLocation.EndsWith(".gz") ?
                    (Stream)(new GZipStream(stream, CompressionMode.Decompress)) :
                    stream;

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
                            var newProteinEntries = block.ParseEndElement(xml, mod_dict, modTypesToExclude, unknownModifications,
                                generateTargetProteins, decoyType, isContaminant, proteinDbLocation);
                            result.AddRange(newProteinEntries);
                        }
                    }
                }
            }
            return result;
        }

        /// <summary>
        /// Get the modification entries specified in a mzLibProteinDb XML file (.xml or .xml.gz).
        /// </summary>
        /// <param name="proteinDbLocation"></param>
        /// <returns></returns>
        [SuppressMessage("Microsoft.Usage", "CA2202:Do not dispose objects multiple times")]
        public static List<Modification> GetPtmListFromProteinXml(string proteinDbLocation)
        {
            if (proteinDbLocation.Equals(last_database_location))
            {
                return protein_xml_modlist;
            }
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
                                protein_xml_modlist = storedKnownModificationsBuilder.Length <= 0 ?
                                    new List<Modification>() :
                                    PtmListLoader.ReadModsFromString(storedKnownModificationsBuilder.ToString()).ToList();
                                return protein_xml_modlist;
                            }
                        }
                    }
                }
            }
            protein_xml_modlist = new List<Modification>();
            return protein_xml_modlist;
        }

        /// <summary>
        /// Load a protein fasta database, using regular expressions to get various aspects of the headers. The first regex capture group is used as each field.
        /// </summary>
        /// <param name="proteinDbLocation"></param>
        /// <param name="originalTarget"></param>
        /// <param name="onTheFlyDecoys"></param>
        /// <param name="isContaminant"></param>
        /// <param name="accessionRegex"></param>
        /// <param name="fullNameRegex"></param>
        /// <param name="nameRegex"></param>
        /// <param name="geneNameRegex"></param>
        /// <param name="organismRegex"></param>
        /// <param name="errors"></param>
        /// <returns></returns>
        public static List<Protein> LoadProteinFasta(string proteinDbLocation, bool originalTarget, DecoyType onTheFlyDecoys, bool isContaminant,
            FastaHeaderFieldRegex accessionRegex, FastaHeaderFieldRegex fullNameRegex, FastaHeaderFieldRegex nameRegex,
            FastaHeaderFieldRegex geneNameRegex, FastaHeaderFieldRegex organismRegex, out List<string> errors)
        {
            HashSet<string> unique_accessions = new HashSet<string>();
            int unique_identifier = 1;
            string accession = null;
            string name = null;
            string fullName = null;
            string organism = null;
            List<Tuple<string, string>> geneName = new List<Tuple<string, string>>();
            errors = new List<string>();
            Regex substituteWhitespace = new Regex(@"\s+");

            List<Protein> result = new List<Protein>();

            using (var stream = new FileStream(proteinDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                Stream fastaFileStream = proteinDbLocation.EndsWith(".gz") ?
                    (Stream)(new GZipStream(stream, CompressionMode.Decompress)) :
                    stream;

                StringBuilder sb = null;
                StreamReader fasta = new StreamReader(fastaFileStream);

                while (true)
                {
                    string line = "";
                    line = fasta.ReadLine();
                    if (line == null) { break; }

                    if (line.StartsWith(">"))
                    {
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
                        while (unique_accessions.Contains(accession))
                        {
                            accession += "_" + unique_identifier.ToString();
                            unique_identifier++;
                        }
                        unique_accessions.Add(accession);
                        if (originalTarget)
                        {
                            Protein protein = new Protein(sequence, accession, organism, geneName, name: name, full_name: fullName,
                                isContaminant: isContaminant, databaseFilePath: proteinDbLocation);
                            if (protein.Length == 0)
                            {
                                errors.Add("Line" + line + ", Protein Length of 0: " + protein.Name + " was skipped from database: " + proteinDbLocation);
                            }
                            else
                            {
                                result.Add(protein);
                            }
                        }

                        if (onTheFlyDecoys == DecoyType.Reverse)
                        {
                            char[] sequence_array = sequence.ToCharArray();
                            int starts_with_met = sequence.StartsWith("M", StringComparison.Ordinal) ? 1 : 0;
                            Array.Reverse(sequence_array, starts_with_met, sequence.Length - starts_with_met); // Do not include the initiator methionine in reversal!!!
                            var reversed_sequence = new string(sequence_array);
                            Protein decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, organism, geneName, name: name, full_name: fullName,
                                isDecoy: true, isContaminant: isContaminant, databaseFilePath: proteinDbLocation);
                            result.Add(decoy_protein);
                        }
                        else if (onTheFlyDecoys == DecoyType.Slide)
                        {
                            int numSlides = 20;
                            char[] sequence_array_unslide = sequence.ToCharArray();
                            char[] sequence_array_slide = sequence.ToCharArray();
                            bool starts_with_met_slide = sequence.StartsWith("M", StringComparison.Ordinal);
                            for (int i = starts_with_met_slide ? 1 : 0; i < sequence.Length; i++)
                            {
                                sequence_array_slide[i] = sequence_array_unslide[DecoyProteinGenerator.GetOldShuffleIndex(i, numSlides, sequence.Length, starts_with_met_slide)];
                            }
                            string slide_sequence = new string(sequence_array_slide);
                            Protein decoy_protein_slide = new Protein(slide_sequence, "DECOY_" + accession, organism, geneName, name: name, full_name: fullName,
                                isDecoy: true, isContaminant: isContaminant, databaseFilePath: proteinDbLocation);
                            result.Add(decoy_protein_slide);
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
            if (!result.Any())
            {
                errors.Add("Error: No proteins could be read from the database: " + proteinDbLocation);
            }
            return result;
        }

        /// <summary>
        /// Merge proteins that have the same accession, sequence, and contaminant designation.
        /// </summary>
        /// <param name="merge_these"></param>
        /// <returns></returns>
        public static IEnumerable<Protein> Merge_proteins(IEnumerable<Protein> merge_these)
        {
            Dictionary<Tuple<string, string, bool, bool>, List<Protein>> proteinsByAccessionSequenceContaminant = new Dictionary<Tuple<string, string, bool, bool>, List<Protein>>();
            foreach (Protein p in merge_these)
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
                    gene_names: genenames.ToList(),
                    oneBasedModifications: mod_dict2,
                    proteolysisProducts: proteolysis.ToList(),
                    name: names.FirstOrDefault(),
                    full_name: fullnames.FirstOrDefault(),
                    databaseReferences: references.ToList(),
                    disulfideBonds: bonds.ToList(),
                    sequenceVariations: variants.ToList()
                    );
            }
        }

        #endregion Public Methods

        #region Private Methods

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
            foreach (Modification nice in mods)
            {
                if (mod_dict.TryGetValue(nice.id, out IList<Modification> val))
                {
                    val.Add(nice);
                }
                else
                {
                    mod_dict.Add(nice.id, new List<Modification> { nice });
                }
            }
            return mod_dict;
        }

        #endregion Private Methods
    }
}