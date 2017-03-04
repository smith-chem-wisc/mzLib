using Ionic.Zlib;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml;
using System.Linq;

namespace UsefulProteomicsDatabases
{
    public static class ProteinDbLoader
    {

        #region Public Methods

        public static List<Protein> LoadProteinXML<T>(string proteinDbLocation, bool onTheFlyDecoys, IEnumerable<T> allKnownModifications, bool IsContaminant, IEnumerable<string> dbRefTypesToKeep, out Dictionary<string, Modification> unknownModifications)
            where T : Modification
        {
            var mod_dict = new Dictionary<string, IList<Modification>>();
            foreach (var nice in allKnownModifications)
            {
                IList<Modification> val;
                if (mod_dict.TryGetValue(nice.id, out val))
                    val.Add(nice);
                else
                    mod_dict.Add(nice.id, new List<Modification> { nice });
            }

            List<Protein> result = new List<Protein>();
            unknownModifications = new Dictionary<string, Modification>();
            using (var stream = new FileStream(proteinDbLocation, FileMode.Open))
            {
                string accession = null;
                string name = null;
                string full_name = null;

                Regex substituteWhitespace = new Regex(@"\s+");

                var oneBasedBeginPositions = new List<int?>();
                var oneBasedEndPositions = new List<int?>();
                var peptideTypes = new List<string>();
                var oneBasedModifications = new Dictionary<int, List<Modification>>();

                Stream uniprotXmlFileStream = proteinDbLocation.EndsWith(".gz") ?
                    (Stream)(new GZipStream(stream, CompressionMode.Decompress)) :
                    stream;

                string[] nodes = new string[6];

                string sequence = null;
                string feature_type = null;
                string feature_description = null;
                string dbReference_type = null;
                string dbReference_id = null;
                List<string> property_types = new List<string>();
                List<string> property_values = new List<string>();
                int oneBasedfeature_position = -1;
                int? oneBasedbeginPosition = null;
                int? oneBasedendPosition = null;
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
                                            List<Modification> residue_modifications;
                                            // Create new entry for this residue, if needed
                                            if (!oneBasedModifications.TryGetValue(oneBasedfeature_position, out residue_modifications))
                                            {
                                                residue_modifications = new List<Modification>();
                                                oneBasedModifications.Add(oneBasedfeature_position, residue_modifications);
                                            }
                                            if (mod_dict.ContainsKey(feature_description))
                                            {
                                                // Known
                                                residue_modifications.AddRange(mod_dict[feature_description]);
                                            }
                                            else if (unknownModifications.ContainsKey(feature_description))
                                            {
                                                // Not known but seen
                                                residue_modifications.Add(unknownModifications[feature_description]);
                                            }
                                            else
                                            {
                                                // Not known and not seen
                                                unknownModifications[feature_description] = new Modification(feature_description);
                                                residue_modifications.Add(unknownModifications[feature_description]);
                                            }
                                        }
                                        else if (feature_type == "peptide" || feature_type == "propeptide" || feature_type == "chain" || feature_type == "signal peptide")
                                        {
                                            oneBasedBeginPositions.Add(oneBasedbeginPosition);
                                            oneBasedEndPositions.Add(oneBasedendPosition);
                                            peptideTypes.Add(feature_type);
                                        }
                                        oneBasedbeginPosition = null;
                                        oneBasedendPosition = null;
                                        oneBasedfeature_position = -1;
                                        break;

                                    case "dbReference":
                                        if (dbRefTypesToKeep != null && dbRefTypesToKeep.Contains(dbReference_type))
                                        {
                                            databaseReferences.Add(new DatabaseReference(dbReference_type, dbReference_id, Enumerable.Range(0, property_types.Count).Select(i => new Tuple<string, string>(property_types[i], property_values[i])).ToList()));
                                        }
                                        property_types = new List<string>();
                                        property_values = new List<string>();
                                        dbReference_type = null;
                                        dbReference_id = null;
                                        break;

                                    case "entry":
                                        if (accession != null && sequence != null)
                                        {
                                            var protein = new Protein(sequence, accession, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, false, IsContaminant, databaseReferences);

                                            result.Add(protein);

                                            if (onTheFlyDecoys)
                                            {
                                                char[] sequence_array = sequence.ToCharArray();
                                                Dictionary<int, List<Modification>> decoy_modifications = null;
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
                                                int?[] decoybeginPositions = new int?[oneBasedBeginPositions.Count];
                                                int?[] decoyendPositions = new int?[oneBasedEndPositions.Count];
                                                string[] decoyBigPeptideTypes = new string[oneBasedEndPositions.Count];
                                                for (int i = 0; i < decoybeginPositions.Length; i++)
                                                {
                                                    decoybeginPositions[oneBasedBeginPositions.Count - i - 1] = sequence.Length - oneBasedEndPositions[i] + 1;
                                                    decoyendPositions[oneBasedBeginPositions.Count - i - 1] = sequence.Length - oneBasedBeginPositions[i] + 1;
                                                    decoyBigPeptideTypes[oneBasedBeginPositions.Count - i - 1] = peptideTypes[i];
                                                }
                                                var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, decoy_modifications, decoybeginPositions, decoyendPositions, decoyBigPeptideTypes, name, full_name, true, IsContaminant, null);

                                                result.Add(decoy_protein);
                                            }
                                        }
                                        accession = null;
                                        name = null;
                                        full_name = null;
                                        sequence = null;
                                        feature_type = null;
                                        feature_description = null;
                                        dbReference_type = null;
                                        dbReference_id = null;
                                        property_types = new List<string>();
                                        property_values = new List<string>();
                                        oneBasedfeature_position = -1;
                                        oneBasedModifications = new Dictionary<int, List<Modification>>();
                                        oneBasedBeginPositions = new List<int?>();
                                        oneBasedEndPositions = new List<int?>();
                                        peptideTypes = new List<string>();
                                        databaseReferences = new List<DatabaseReference>();
                                        break;
                                }
                                break;
                        }
                    }
                }
            }
            return result;
        }

        public static Regex uniprot_accession_expression = new Regex(@"([A-Z0-9_]+)\w");

        public static Regex uniprot_fullName_expression = new Regex(@"([a-zA-Z0-9_ ]+)(OS=)\w");

        public static Regex ensembl_accession_expression = new Regex(@"([A-Z0-9_]+)\w");

        public static Regex ensembl_fullName_expression = new Regex(@"(pep:.*)");

        public static List<Protein> LoadProteinFasta(string proteinDbLocation, bool onTheFlyDecoys, bool IsContaminant, Regex accession_expression, Regex full_name_expression, Regex name_expression)
        {
            string accession = null;
            string name = null;
            string full_name = null;

            Regex substituteWhitespace = new Regex(@"\s+");

            int?[] oneBasedBeginPositions = null;
            int?[] oneBasedEndPositions = null;
            string[] productTypes = null;
            Dictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>();
            List<DatabaseReference> databaseReferences = new List<DatabaseReference>();
            List<Protein> result = new List<Protein>();

            using (var stream = new FileStream(proteinDbLocation, FileMode.Open))
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
                        accession = accession_expression.Match(line).Value;
                        if (accession == null || accession == "") accession = line.Substring(1);
                        full_name = full_name_expression.Match(line).Value;
                        name = name_expression.Match(line).Value;
                        sb = new StringBuilder();
                    }

                    else if (sb != null)
                    {
                        sb.Append(line.Trim());
                    }

                    if ((fasta.Peek() == '>' || fasta.Peek() == -1) && accession != null && sb != null)
                    {
                        string sequence = substituteWhitespace.Replace(sb.ToString(), "");
                        Protein protein = new Protein(sequence, accession, oneBasedModifications, oneBasedBeginPositions, oneBasedEndPositions, productTypes, name, full_name, false, IsContaminant, databaseReferences);
                        result.Add(protein);

                        if (onTheFlyDecoys)
                        {
                            char[] sequence_array = sequence.ToCharArray();
                            int starts_with_met = Convert.ToInt32(sequence.StartsWith("M", StringComparison.InvariantCulture));
                            Array.Reverse(sequence_array, starts_with_met, sequence.Length - starts_with_met); // Do not include the initiator methionine in reversal!!!
                            var reversed_sequence = new string(sequence_array);
                            Protein decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, oneBasedModifications, oneBasedBeginPositions, oneBasedEndPositions, productTypes, name, full_name, true, IsContaminant, null);
                            result.Add(decoy_protein);
                        }
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

        #endregion Public Methods

    }
}