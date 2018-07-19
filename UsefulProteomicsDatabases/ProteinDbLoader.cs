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
        public static readonly FastaHeaderFieldRegex UniprotAccessionRegex = new FastaHeaderFieldRegex("accession", @"([A-Z0-9_.]+)", 0, 1);
        public static readonly FastaHeaderFieldRegex UniprotFullNameRegex = new FastaHeaderFieldRegex("fullName", @"\s(.*?)\sOS=", 0, 1);
        public static readonly FastaHeaderFieldRegex UniprotNameRegex = new FastaHeaderFieldRegex("name", @"\|([^\|][A-Z0-9_]+)", 1, 1);
        public static readonly FastaHeaderFieldRegex UniprotGeneNameRegex = new FastaHeaderFieldRegex("geneName", @"GN=([^ ]+)", 0, 1);
        public static readonly FastaHeaderFieldRegex UniprotOrganismRegex = new FastaHeaderFieldRegex("organism", @"OS=(.*?)\sGN=", 0, 1);

        public static readonly FastaHeaderFieldRegex EnsemblAccessionRegex = new FastaHeaderFieldRegex("accession", @"([A-Z0-9_.]+)", 0, 1);
        public static readonly FastaHeaderFieldRegex EnsemblFullNameRegex = new FastaHeaderFieldRegex("fullName", @"(pep:.*)", 0, 1);
        public static readonly FastaHeaderFieldRegex EnsemblGeneNameRegex = new FastaHeaderFieldRegex("geneName", @"gene:([^ ]+)", 0, 1);

        /// <summary>
        /// Stores the last database file path.
        /// </summary>
        private static string last_database_location;

        /// <summary>
        /// Stores the modification list read during LoadProteinXML
        /// </summary>
        private static List<ModificationGeneral> protein_xml_modlist_general;

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
        public static List<Protein> LoadProteinXML(string proteinDbLocation, bool generateTargets, DecoyType decoyType, IEnumerable<ModificationGeneral> allKnownModifications,
            bool isContaminant, IEnumerable<string> modTypesToExclude, out Dictionary<string, ModificationGeneral> unknownModifications, int maxThreads = -1)
        {
            List<ModificationGeneral> prespecified = GetPtmListFromProteinXml(proteinDbLocation);
            allKnownModifications = allKnownModifications ?? new List<ModificationGeneral>();
            modTypesToExclude = modTypesToExclude ?? new List<string>();

            Dictionary<string, IList<ModificationGeneral>> mod_dict = new Dictionary<string, IList<ModificationGeneral>>();
            if (prespecified.Count > 0 || allKnownModifications.Count() > 0)
            {
                mod_dict = GetModificationDict(new HashSet<ModificationGeneral>(prespecified.Concat(allKnownModifications)));
            }

            List<Protein> targets = new List<Protein>();
            unknownModifications = new Dictionary<string, ModificationGeneral>();
            using (var stream = new FileStream(proteinDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                Regex substituteWhitespace = new Regex(@"\s+");

                Stream uniprotXmlFileStream = proteinDbLocation.EndsWith("gz") ? // allow for .bgz and .tgz, which are (rarely) used
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
                                isContaminant, proteinDbLocation);
                            targets.AddRange(newProteinEntries);
                        }
                    }
                }
            }
            List<Protein> decoys = DecoyProteinGenerator.GenerateDecoys(targets, decoyType, maxThreads);
            return (generateTargets ? targets : new List<Protein>()).Concat(decoyType != DecoyType.None ? decoys : new List<Protein>()).ToList();
        }

        /// <summary>
        /// Get the modification entries specified in a mzLibProteinDb XML file (.xml or .xml.gz).
        /// </summary>
        /// <param name="proteinDbLocation"></param>
        /// <returns></returns>
        [SuppressMessage("Microsoft.Usage", "CA2202:Do not dispose objects multiple times")]
        public static List<ModificationGeneral> GetPtmListFromProteinXml(string proteinDbLocation)
        {
            if (proteinDbLocation.Equals(last_database_location))
            {
                return protein_xml_modlist_general;
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
                                protein_xml_modlist_general = storedKnownModificationsBuilder.Length <= 0 ?
                                    new List<ModificationGeneral>() :
                                    PtmListLoaderGeneral.ReadModsFromString(storedKnownModificationsBuilder.ToString()).ToList();
                                return protein_xml_modlist_general;
                            }
                        }
                    }
                }
            }
            protein_xml_modlist_general = new List<ModificationGeneral>();
            return protein_xml_modlist_general;
        }

        /// <summary>
        /// Load a protein fasta database, using regular expressions to get various aspects of the headers. The first regex capture group is used as each field.
        /// </summary>
        /// <param name="proteinDbLocation"></param>
        /// <param name="generateTargets"></param>
        /// <param name="decoyType"></param>
        /// <param name="isContaminant"></param>
        /// <param name="accessionRegex"></param>
        /// <param name="fullNameRegex"></param>
        /// <param name="nameRegex"></param>
        /// <param name="geneNameRegex"></param>
        /// <param name="organismRegex"></param>
        /// <param name="errors"></param>
        /// <returns></returns>
        public static List<Protein> LoadProteinFasta(string proteinDbLocation, bool generateTargets, DecoyType decoyType, bool isContaminant,
            FastaHeaderFieldRegex accessionRegex, FastaHeaderFieldRegex fullNameRegex, FastaHeaderFieldRegex nameRegex,
            FastaHeaderFieldRegex geneNameRegex, FastaHeaderFieldRegex organismRegex, out List<string> errors, int maxThreads = -1)
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

            List<Protein> targets = new List<Protein>();

            using (var stream = new FileStream(proteinDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                Stream fastaFileStream = proteinDbLocation.EndsWith("gz") ? // allow for .bgz and .tgz, which are (rarely) used
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
                        Protein protein = new Protein(sequence, accession, organism, geneName, name: name, full_name: fullName,
                            isContaminant: isContaminant, databaseFilePath: proteinDbLocation);
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
            if (!targets.Any())
            {
                errors.Add("Error: No proteins could be read from the database: " + proteinDbLocation);
            }
            List<Protein> decoys = DecoyProteinGenerator.GenerateDecoys(targets, decoyType, maxThreads);
            return (generateTargets ? targets : new List<Protein>()).Concat(decoyType != DecoyType.None ? decoys : new List<Protein>()).ToList();
        }

        /// <summary>
        /// Merge proteins that have the same accession, sequence, and contaminant designation.
        /// </summary>
        /// <param name="mergeThese"></param>
        /// <returns></returns>
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

                Dictionary<int, HashSet<ModificationGeneral>> mod_dict = new Dictionary<int, HashSet<ModificationGeneral>>();
                foreach (KeyValuePair<int, List<ModificationGeneral>> nice in proteins.Value.SelectMany(p => p.OneBasedPossibleLocalizedModifications).ToList())
                {
                    if (!mod_dict.TryGetValue(nice.Key, out HashSet<ModificationGeneral> val))
                    {
                        mod_dict.Add(nice.Key, new HashSet<ModificationGeneral>(nice.Value));
                    }
                    else
                    {
                        foreach (ModificationGeneral mod in nice.Value)
                        {
                            val.Add(mod);
                        }
                    }
                }
                Dictionary<int, List<ModificationGeneral>> mod_dict2 = mod_dict.ToDictionary(kv => kv.Key, kv => kv.Value.ToList());

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

        private static Dictionary<string, IList<ModificationGeneral>> GetModificationDict(IEnumerable<ModificationGeneral> mods)
        {
            var mod_dict = new Dictionary<string, IList<ModificationGeneral>>();
            foreach (ModificationGeneral nice in mods)
            {
                if (mod_dict.TryGetValue(nice.Id, out IList<ModificationGeneral> val))
                {
                    val.Add(nice);
                }
                else
                {
                    mod_dict.Add(nice.Id, new List<ModificationGeneral> { nice });
                }
            }
            return mod_dict;
        }
    }
}