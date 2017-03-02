using Ionic.Zlib;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml;

namespace UsefulProteomicsDatabases
{
    public static class ProteinDbLoader
    {

        #region Public Methods

        public static List<Protein> LoadProteinDb<T>(string proteinDbLocation, bool onTheFlyDecoys, IEnumerable<T> allKnownModifications, bool IsContaminant, out Dictionary<string, Modification> unknownModifications)
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

                // xml db
                if (!proteinDbLocation.EndsWith(".fasta"))
                {
                    var oneBasedBeginPositions = new List<int?>();
                    var oneBasedEndPositions = new List<int?>();
                    var peptideTypes = new List<string>();
                    var oneBasedModifications = new Dictionary<int, List<Modification>>();

                    Stream uniprotXmlFileStream = stream;
                    if (proteinDbLocation.EndsWith(".gz"))
                        uniprotXmlFileStream = new GZipStream(stream, CompressionMode.Decompress);

                    string[] nodes = new string[6];

                    string sequence = null;
                    string feature_type = null;
                    string feature_description = null;
                    string dbReference_type = null;
                    string dbReference_id = null;
                    string property_type = null;
                    string property_value = null;
                    int oneBasedfeature_position = -1;
                    int? oneBasedbeginPosition = null;
                    int? oneBasedendPosition = null;
                    List<GoTerm> goTerms = new List<GoTerm>();

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
                                            dbReference_type = xml.GetAttribute("type");
                                            dbReference_id = xml.GetAttribute("id");
                                            break;

                                        case "property":
                                            property_type = xml.GetAttribute("type");
                                            property_value = xml.GetAttribute("value");
                                            if (dbReference_type == "GO" && property_type == "term")
                                            {
                                                switch (property_value.Split(':')[0].ToString())
                                                {
                                                    case "C":
                                                        goTerms.Add(new GoTerm(dbReference_id.Split(':')[1].ToString(),
                                                                               property_value.Split(':')[1].ToString(),
                                                                               Aspect.cellularComponent));
                                                        break;

                                                    case "F":
                                                        goTerms.Add(new GoTerm(dbReference_id.Split(':')[1].ToString(),
                                                                               property_value.Split(':')[1].ToString(),
                                                                               Aspect.molecularFunction));
                                                        break;

                                                    case "P":
                                                        goTerms.Add(new GoTerm(dbReference_id.Split(':')[1].ToString(),
                                                                               property_value.Split(':')[1].ToString(),
                                                                               Aspect.biologicalProcess));
                                                        break;
                                                }
                                            }
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
                                            dbReference_type = null;
                                            dbReference_id = null;
                                            break;

                                        case "entry":
                                            if (accession != null && sequence != null)
                                            {
                                                var protein = new Protein(sequence, accession, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, false, IsContaminant, goTerms);

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
                                            property_type = null;
                                            property_value = null;
                                            oneBasedfeature_position = -1;
                                            oneBasedModifications = new Dictionary<int, List<Modification>>();
                                            oneBasedBeginPositions = new List<int?>();
                                            oneBasedEndPositions = new List<int?>();
                                            peptideTypes = new List<string>();
                                            goTerms = new List<GoTerm>();
                                            break;
                                    }
                                    break;
                            }
                        }
                    }
                }

                // fasta db
                else
                {
                    StreamReader fasta = new StreamReader(stream);

                    StringBuilder sb = null;
                    while (true)
                    {
                        string line = fasta.ReadLine();

                        if (line.StartsWith(">"))
                        {
                            // fasta protein only has accession, fullname, sequence (no mods)
                            string[] delimiters = { ">", "|", " OS=" };
                            string[] output = line.Split(delimiters, StringSplitOptions.None);
                            if (output.Length > 4)
                            {
                                accession = output[2];
                                name = accession;
                                full_name = output[3];
                            }
                            else
                            {
                                // can't read protein description
                                full_name = line.Substring(1);
                                accession = line.Substring(1);
                            }

                            // new protein
                            sb = new StringBuilder();
                        }
                        else if (sb != null)
                        {
                            sb.Append(line.Trim());
                        }

                        if ((fasta.Peek() == '>' || fasta.Peek() == -1) && accession != null && sb != null)
                        {
                            var sequence = substituteWhitespace.Replace(sb.ToString(), "");
                            var protein = new Protein(sequence, accession, name, full_name, false, IsContaminant);

                            result.Add(protein);

                            if (onTheFlyDecoys)
                            {
                                char[] sequence_array = sequence.ToCharArray();
                                int starts_with_met = Convert.ToInt32(sequence.StartsWith("M", StringComparison.InvariantCulture));
                                Array.Reverse(sequence_array, starts_with_met, sequence.Length - starts_with_met); // Do not include the initiator methionine in reversal!!!
                                var reversed_sequence = new string(sequence_array);
                                var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, name, full_name, true, IsContaminant);

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
            }
            return result;
        }

        #endregion Public Methods

    }
}