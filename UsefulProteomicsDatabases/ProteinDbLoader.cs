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

        public static Regex uniprot_accession_expression = new Regex(@"([A-Z0-9_]+)");

        public static Regex uniprot_fullName_expression = new Regex(@"\|([^\|]+)\sOS=");

        public static Regex uniprot_gene_expression = new Regex(@"GN=([^ ]+)");

        public static Regex ensembl_accession_expression = new Regex(@"([A-Z0-9_]+)");

        public static Regex ensembl_fullName_expression = new Regex(@"(pep:.*)");

        public static Regex ensembl_gene_expression = new Regex(@"gene:([^ ]+)");

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
        /// <param name="IsContaminant"></param>
        /// <param name="dbRefTypesToKeep"></param>
        /// <param name="modTypesToExclude"></param>
        /// <param name="unknownModifications"></param>
        /// <returns></returns>
        [SuppressMessage("Microsoft.Usage", "CA2202:Do not dispose objects multiple times")]
        public static List<Protein> LoadProteinXML(string proteinDbLocation, bool generateTargetProteins, bool generateDecoyProteins, IEnumerable<Modification> allKnownModifications, bool IsContaminant, IEnumerable<string> modTypesToExclude, out Dictionary<string, Modification> unknownModifications)
        {
            List<Modification> prespecified = GetPtmListFromProteinXml(proteinDbLocation);

            Dictionary<string, IList<Modification>> mod_dict = new Dictionary<string, IList<Modification>>();
            if (prespecified.Count > 0 || allKnownModifications.Count() > 0)
                mod_dict = GetModificationDict(new HashSet<Modification>(prespecified.Concat(allKnownModifications)));

            List<Protein> result = new List<Protein>();
            unknownModifications = new Dictionary<string, Modification>();
            using (var stream = new FileStream(proteinDbLocation, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                Regex substituteWhitespace = new Regex(@"\s+");

                Stream uniprotXmlFileStream = proteinDbLocation.EndsWith(".gz") ?
                    (Stream)(new GZipStream(stream, CompressionMode.Decompress)) :
                    stream;

                string[] nodes = new string[6];

                string accession = null;
                string name = null;
                string full_name = null;
                string sequence = null;
                string feature_type = null;
                string feature_description = null;
                string original_value = ""; // if no content is found, assume it is empty, not null (e.g. <original>A</original><variation/> for a deletion event)
                string variation_value = "";
                string dbReference_type = null;
                string dbReference_id = null;
                List<string> property_types = new List<string>();
                List<string> property_values = new List<string>();
                int oneBasedfeature_position = -1;
                int? oneBasedbeginPosition = null;
                int? oneBasedendPosition = null;
                List<ProteolysisProduct> proteolysisProducts = new List<ProteolysisProduct>();
                List<SequenceVariation> sequenceVariations = new List<SequenceVariation>();
                List<DisulfideBond> disulfideBonds = new List<DisulfideBond>();
                var oneBasedModifications = new Dictionary<int, List<Modification>>();
                List<Tuple<string, string>> gene_names = new List<Tuple<string, string>>();
                bool reading_gene = false;
                List<DatabaseReference> databaseReferences = new List<DatabaseReference>();

                using (XmlReader xml = XmlReader.Create(uniprotXmlFileStream))
                {
                    while (xml.Read())
                    {
                        switch (xml.NodeType)
                        {
                            case XmlNodeType.Element:
                                nodes[xml.Depth] = xml.Name;
                                int outValue;
                                switch (xml.Name)
                                {
                                    case "accession":
                                        if (accession == null)
                                        {
                                            accession = xml.ReadElementString();
                                        }
                                        break;

                                    case "name":
                                        if (xml.Depth == 2)
                                        {
                                            name = xml.ReadElementString();
                                        }
                                        if (reading_gene) gene_names.Add(new Tuple<string, string>(xml.GetAttribute("type"), xml.ReadElementString()));
                                        break;

                                    case "gene":
                                        reading_gene = true;
                                        break;

                                    case "fullName":
                                        if (full_name == null)
                                        {
                                            full_name = xml.ReadElementString();
                                        }
                                        break;

                                    case "feature":
                                        feature_type = xml.GetAttribute("type");
                                        feature_description = xml.GetAttribute("description");
                                        break;

                                    case "original":
                                        original_value = xml.ReadElementString();
                                        break;

                                    case "variation":
                                        variation_value = xml.ReadElementString();
                                        break;

                                    case "dbReference":
                                        property_types.Clear();
                                        property_values.Clear();
                                        dbReference_type = xml.GetAttribute("type");
                                        dbReference_id = xml.GetAttribute("id");
                                        break;

                                    case "property":
                                        property_types.Add(xml.GetAttribute("type"));
                                        property_values.Add(xml.GetAttribute("value"));
                                        break;

                                    case "position":
                                        oneBasedfeature_position = int.Parse(xml.GetAttribute("position"));
                                        break;

                                    case "begin":
                                        oneBasedbeginPosition = int.TryParse(xml.GetAttribute("position"), out outValue) ? (int?)outValue : null;
                                        break;

                                    case "end":
                                        oneBasedendPosition = int.TryParse(xml.GetAttribute("position"), out outValue) ? (int?)outValue : null;
                                        break;

                                    case "sequence":
                                        sequence = substituteWhitespace.Replace(xml.ReadElementString(), "");
                                        break;
                                }
                                break;

                            case XmlNodeType.EndElement:
                                switch (xml.Name)
                                {
                                    case "feature":
                                        if (feature_type == "modified residue")
                                        {
                                            feature_description = feature_description.Split(';')[0];

                                            // Create new entry for this residue, if needed
                                            if (!oneBasedModifications.TryGetValue(oneBasedfeature_position, out List<Modification> residue_modifications))
                                            {
                                                residue_modifications = new List<Modification>();
                                                oneBasedModifications.Add(oneBasedfeature_position, residue_modifications);
                                            }
                                            if (mod_dict.ContainsKey(feature_description))
                                            {
                                                // Known and not of a type in the exclusion list
                                                List<Modification> mods = mod_dict[feature_description].Where(m => !modTypesToExclude.Contains(m.modificationType)).ToList();
                                                if (mods.Count == 0 && oneBasedModifications[oneBasedfeature_position].Count == 0)
                                                {
                                                    oneBasedModifications.Remove(oneBasedfeature_position);
                                                }
                                                else
                                                {
                                                    oneBasedModifications[oneBasedfeature_position].AddRange(mods);
                                                }
                                            }
                                            else if (unknownModifications.ContainsKey(feature_description))
                                            {
                                                // Not known but seen
                                                residue_modifications.Add(unknownModifications[feature_description]);
                                            }
                                            else
                                            {
                                                // Not known and not seen
                                                unknownModifications[feature_description] = new Modification(feature_description, "unknown");
                                                residue_modifications.Add(unknownModifications[feature_description]);
                                            }
                                        }
                                        else if (feature_type == "peptide" || feature_type == "propeptide" || feature_type == "chain" || feature_type == "signal peptide")
                                        {
                                            proteolysisProducts.Add(new ProteolysisProduct(oneBasedbeginPosition, oneBasedendPosition, feature_type));
                                        }
                                        else if (feature_type == "sequence variant" // Only keep if there is variant sequence information and position information
                                            && variation_value != null
                                            && variation_value != "")
                                        {
                                            if (oneBasedbeginPosition != null && oneBasedendPosition != null)
                                                sequenceVariations.Add(new SequenceVariation((int)oneBasedbeginPosition, (int)oneBasedendPosition, original_value, variation_value, feature_description));
                                            else if (oneBasedfeature_position >= 1)
                                                sequenceVariations.Add(new SequenceVariation(oneBasedfeature_position, original_value, variation_value, feature_description));
                                        }
                                        else if (feature_type == "disulfide bond")
                                        {
                                            if (oneBasedbeginPosition != null && oneBasedendPosition != null)
                                                disulfideBonds.Add(new DisulfideBond((int)oneBasedbeginPosition, (int)oneBasedendPosition, feature_description));
                                            else if (oneBasedfeature_position >= 1)
                                                disulfideBonds.Add(new DisulfideBond(oneBasedfeature_position, feature_description));
                                        }
                                        oneBasedbeginPosition = null;
                                        oneBasedendPosition = null;
                                        oneBasedfeature_position = -1;
                                        original_value = "";
                                        variation_value = "";
                                        break;

                                    case "dbReference":
                                        databaseReferences.Add(new DatabaseReference(dbReference_type, dbReference_id, Enumerable.Range(0, property_types.Count).Select(i => new Tuple<string, string>(property_types[i], property_values[i])).ToList()));
                                        property_types = new List<string>();
                                        property_values = new List<string>();
                                        dbReference_type = null;
                                        dbReference_id = null;
                                        break;

                                    case "gene":
                                        reading_gene = false;
                                        break;

                                    case "entry":
                                        if (accession != null && sequence != null)
                                        {
                                            if (generateTargetProteins)
                                            {
                                                var protein = new Protein(sequence, accession, gene_names, oneBasedModifications, proteolysisProducts, name, full_name, false, IsContaminant, databaseReferences, sequenceVariations, disulfideBonds, proteinDbLocation);
                                                result.Add(protein);
                                            }

                                            if (generateDecoyProteins)
                                            {
                                                char[] sequence_array = sequence.ToCharArray();
                                                Dictionary<int, List<Modification>> decoy_modifications = null;
                                                List<DisulfideBond> decoy_disulfides = new List<DisulfideBond>();
                                                if (sequence.StartsWith("M", StringComparison.InvariantCulture))
                                                {
                                                    // Do not include the initiator methionine in reversal!!!
                                                    Array.Reverse(sequence_array, 1, sequence.Length - 1);
                                                    decoy_modifications = new Dictionary<int, List<Modification>>(oneBasedModifications.Count);
                                                    foreach (var kvp in oneBasedModifications)
                                                    {
                                                        if (kvp.Key == 1)
                                                        {
                                                            decoy_modifications.Add(1, kvp.Value);
                                                        }
                                                        else if (kvp.Key > 1)
                                                        {
                                                            decoy_modifications.Add(sequence.Length - kvp.Key + 2, kvp.Value);
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    Array.Reverse(sequence_array);
                                                    decoy_modifications = new Dictionary<int, List<Modification>>(oneBasedModifications.Count);
                                                    foreach (var kvp in oneBasedModifications)
                                                    {
                                                        decoy_modifications.Add(sequence.Length - kvp.Key + 1, kvp.Value);
                                                    }
                                                }
                                                var reversed_sequence = new string(sequence_array);

                                                List<ProteolysisProduct> decoyPP = new List<ProteolysisProduct>();
                                                foreach (ProteolysisProduct pp in proteolysisProducts)
                                                {
                                                    decoyPP.Add(new ProteolysisProduct(sequence.Length - pp.OneBasedEndPosition + 1, sequence.Length - pp.OneBasedBeginPosition, pp.Type));
                                                }
                                                foreach (DisulfideBond disulfideBond in disulfideBonds)
                                                {
                                                    decoy_disulfides.Add(new DisulfideBond(reversed_sequence.Length - disulfideBond.OneBasedBeginPosition + 1, reversed_sequence.Length - disulfideBond.OneBasedEndPosition + 1, "DECOY DISULFIDE BOND: " + disulfideBond.Description));
                                                }

                                                List<SequenceVariation> decoy_variations = new List<SequenceVariation>();
                                                foreach (SequenceVariation sv in sequenceVariations)
                                                {
                                                    char[] original_array = sv.OriginalSequence.ToArray();
                                                    char[] variation_array = sv.VariantSequence.ToArray();
                                                    if (sv.OneBasedBeginPosition == 1)
                                                    {
                                                        bool orig_init_m = sv.OriginalSequence.StartsWith("M", StringComparison.InvariantCulture);
                                                        bool var_init_m = sv.VariantSequence.StartsWith("M", StringComparison.InvariantCulture);
                                                        if (orig_init_m && !var_init_m)
                                                            decoy_variations.Add(new SequenceVariation(1, "M", "", "DECOY VARIANT: Initiator Methionine Change in " + sv.Description));
                                                        original_array = sv.OriginalSequence.Substring(Convert.ToInt32(orig_init_m)).ToArray();
                                                        variation_array = sv.VariantSequence.Substring(Convert.ToInt32(var_init_m)).ToArray();
                                                    }
                                                    int decoy_end = sequence.Length - sv.OneBasedBeginPosition + 2 + Convert.ToInt32(sv.OneBasedEndPosition == reversed_sequence.Length) - Convert.ToInt32(sv.OneBasedBeginPosition == 1);
                                                    int decoy_begin = decoy_end - original_array.Length + 1;
                                                    Array.Reverse(original_array);
                                                    Array.Reverse(variation_array);
                                                    decoy_variations.Add(new SequenceVariation(decoy_begin, decoy_end, new string(original_array), new string(variation_array), "DECOY VARIANT: " + sv.Description));
                                                }
                                                var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, gene_names, decoy_modifications, decoyPP, name, full_name, true, IsContaminant, null, decoy_variations, decoy_disulfides, proteinDbLocation);

                                                result.Add(decoy_protein);
                                            }
                                        }
                                        accession = null;
                                        name = null;
                                        full_name = null;
                                        sequence = null;
                                        feature_type = null;
                                        feature_description = null;
                                        original_value = "";
                                        variation_value = "";
                                        dbReference_type = null;
                                        dbReference_id = null;
                                        property_types = new List<string>();
                                        property_values = new List<string>();
                                        oneBasedfeature_position = -1;
                                        oneBasedModifications = new Dictionary<int, List<Modification>>();
                                        proteolysisProducts = new List<ProteolysisProduct>();
                                        sequenceVariations = new List<SequenceVariation>();
                                        disulfideBonds = new List<DisulfideBond>();
                                        databaseReferences = new List<DatabaseReference>();
                                        gene_names = new List<Tuple<string, string>>();
                                        reading_gene = false;
                                        break;
                                }
                                break;
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
                return protein_xml_modlist;
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
                        switch (xml.NodeType)
                        {
                            case XmlNodeType.Element:
                                switch (xml.Name)
                                {
                                    case "modification":
                                        string modification = startingWhitespace.Replace(xml.ReadElementString(), "");
                                        storedKnownModificationsBuilder.AppendLine(modification);
                                        break;

                                    case "entry":
                                        if (storedKnownModificationsBuilder.Length <= 0)
                                            protein_xml_modlist = new List<Modification>();
                                        else
                                            protein_xml_modlist = PtmListLoader.ReadModsFromString(storedKnownModificationsBuilder.ToString()).ToList<Modification>();
                                        return protein_xml_modlist;
                                }
                                break;
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
        /// <param name="onTheFlyDecoys"></param>
        /// <param name="IsContaminant"></param>
        /// <param name="accession_expression"></param>
        /// <param name="full_name_expression"></param>
        /// <param name="name_expression"></param>
        /// <param name="gene_expression"></param>
        /// <returns></returns>
        public static List<Protein> LoadProteinFasta(string proteinDbLocation, bool originalTarget, bool onTheFlyDecoys, bool IsContaminant, Regex accession_expression, Regex full_name_expression, Regex name_expression, Regex gene_expression)
        {
            HashSet<string> unique_accessions = new HashSet<string>();
            int unique_identifier = 1;
            string accession = null;
            string name = null;
            string full_name = null;
            List<Tuple<string, string>> gene_name = new List<Tuple<string, string>>();

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
                    string line = fasta.ReadLine();

                    if (line.StartsWith(">"))
                    {
                        var accession_match = accession_expression.Match(line);
                        var full_name_match = full_name_expression.Match(line);
                        var name_match = name_expression.Match(line);
                        var gene_name_match = gene_expression.Match(line);

                        if (accession_match.Groups.Count > 1) accession = accession_expression.Match(line).Groups[1].Value;
                        if (full_name_match.Groups.Count > 1) full_name = full_name_expression.Match(line).Groups[1].Value;
                        if (name_match.Groups.Count > 1) name = name_expression.Match(line).Groups[1].Value;
                        if (gene_name_match.Groups.Count > 1) gene_name.Add(new Tuple<string, string>("primary", gene_expression.Match(line).Groups[1].Value));

                        if (accession == null || accession == "")
                            accession = line.Substring(1).TrimEnd();

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
                            Protein protein = new Protein(sequence, accession, gene_name, name: name, full_name: full_name, isContaminant: IsContaminant, databaseFilePath: proteinDbLocation);
                            result.Add(protein);
                        }

                        if (onTheFlyDecoys)
                        {
                            char[] sequence_array = sequence.ToCharArray();
                            int starts_with_met = Convert.ToInt32(sequence.StartsWith("M", StringComparison.InvariantCulture));
                            Array.Reverse(sequence_array, starts_with_met, sequence.Length - starts_with_met); // Do not include the initiator methionine in reversal!!!
                            var reversed_sequence = new string(sequence_array);
                            Protein decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, gene_name, name: name, full_name: full_name, isDecoy: true, isContaminant: IsContaminant, databaseFilePath: proteinDbLocation);
                            result.Add(decoy_protein);
                        }

                        accession = null;
                        name = null;
                        full_name = null;
                        gene_name = new List<Tuple<string, string>>();
                    }

                    // no input left
                    if (fasta.Peek() == -1)
                    {
                        break;
                    }
                }
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
                    proteinsByAccessionSequenceContaminant.Add(key, new List<Protein> { p });
                else
                    bundled.Add(p);
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
                        mod_dict.Add(nice.Key, new HashSet<Modification>(nice.Value));
                    else
                        foreach (Modification mod in nice.Value)
                        {
                            val.Add(mod);
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

        private static Dictionary<string, IList<Modification>> GetModificationDict(IEnumerable<Modification> mods)
        {
            var mod_dict = new Dictionary<string, IList<Modification>>();
            foreach (Modification nice in mods)
            {
                if (mod_dict.TryGetValue(nice.id, out IList<Modification> val))
                    val.Add(nice);
                else
                    mod_dict.Add(nice.id, new List<Modification> { nice });
            }
            return mod_dict;
        }

        #endregion Private Methods
    }
}