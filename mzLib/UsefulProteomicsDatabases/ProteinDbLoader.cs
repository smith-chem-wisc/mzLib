using Ionic.Zlib;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Xml;

namespace UsefulProteomicsDatabases
{
    public static class ProteinDbLoader
    {
        #region Public Methods

        public static IEnumerable<Protein> LoadProteinDb(string proteinDbLocation, bool onTheFlyDecoys, IDictionary<string, HashSet<BaseModification>> allKnownModifications, bool IsContaminant)
        {
            using (var stream = new FileStream(proteinDbLocation, FileMode.Open))
            {
                string accession = null;
                string name = null;
                string full_name = null;
                int offset = 0;
                string sequence = null;

                // xml db
                if (!proteinDbLocation.EndsWith(".fasta"))
                {
                    Stream uniprotXmlFileStream = stream;
                    if (proteinDbLocation.EndsWith(".gz"))
                        uniprotXmlFileStream = new GZipStream(stream, CompressionMode.Decompress);

                    string[] nodes = new string[6];

                    string feature_type = null;
                    string feature_description = null;
                    int oneBasedfeature_position = -1;
                    int oneBasedbeginPosition = -1;
                    int oneBasedendPosition = -1;
                    var oneBasedBeginPositions = new List<int>();
                    var oneBasedEndPositions = new List<int>();
                    var peptideTypes = new List<string>();
                    var oneBasedModifications = new Dictionary<int, HashSet<BaseModification>>();
                    string db = null;

                    using (XmlReader xml = XmlReader.Create(uniprotXmlFileStream))
                    {
                        while (xml.Read())
                        {
                            switch (xml.NodeType)
                            {
                                case XmlNodeType.Element:
                                    nodes[xml.Depth] = xml.Name;
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

                                        case "position":
                                            oneBasedfeature_position = int.Parse(xml.GetAttribute("position"));
                                            break;

                                        case "begin":
                                            try
                                            {
                                                oneBasedbeginPosition = int.Parse(xml.GetAttribute("position"));
                                            }
                                            catch (ArgumentNullException)
                                            {
                                            }
                                            break;

                                        case "end":
                                            try
                                            {
                                                oneBasedendPosition = int.Parse(xml.GetAttribute("position"));
                                            }
                                            catch (ArgumentNullException)
                                            {
                                            }
                                            break;

                                        case "sequence":
                                            sequence = xml.ReadElementString().Replace("\n", null);
                                            break;

                                        case "db":
                                            db = xml.ReadElementString();
                                            break;
                                    }
                                    break;

                                case XmlNodeType.EndElement:
                                    switch (xml.Name)
                                    {
                                        case "feature":
                                            if (feature_type == "modified residue")
                                            {
                                                // Create new entry for this residue, if needed
                                                HashSet<BaseModification> residue_modifications;
                                                if (!oneBasedModifications.TryGetValue(oneBasedfeature_position, out residue_modifications))
                                                {
                                                    residue_modifications = new HashSet<BaseModification>();
                                                    oneBasedModifications.Add(oneBasedfeature_position, residue_modifications);
                                                }
                                                // Known modification
                                                if (!allKnownModifications.ContainsKey(feature_description))
                                                    allKnownModifications.Add(feature_description, new HashSet<BaseModification> { new BaseModification(feature_description) });
                                                residue_modifications.UnionWith(allKnownModifications[feature_description]);
                                            }
                                            else if ((feature_type == "peptide" || feature_type == "propeptide" || feature_type == "chain") && oneBasedbeginPosition >= 0 && oneBasedendPosition >= 0)
                                            {
                                                oneBasedBeginPositions.Add(oneBasedbeginPosition);
                                                oneBasedEndPositions.Add(oneBasedendPosition);
                                                peptideTypes.Add(feature_type);
                                            }
                                            oneBasedbeginPosition = -1;
                                            oneBasedendPosition = -1;
                                            db = null;

                                            break;

                                        case "entry":
                                            if (accession != null && sequence != null)
                                            {
                                                var protein = new Protein(sequence, accession, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, offset, false, IsContaminant);

                                                yield return protein;

                                                offset += protein.Length;
                                                if (onTheFlyDecoys)
                                                {
                                                    char[] sequence_array = sequence.ToCharArray();
                                                    Dictionary<int, HashSet<BaseModification>> decoy_modifications = null;
                                                    if (sequence.StartsWith("M", StringComparison.InvariantCulture))
                                                    {
                                                        // Do not include the initiator methionine in reversal!!!
                                                        Array.Reverse(sequence_array, 1, sequence.Length - 1);
                                                        if (oneBasedModifications != null)
                                                        {
                                                            decoy_modifications = new Dictionary<int, HashSet<BaseModification>>(oneBasedModifications.Count);
                                                            foreach (KeyValuePair<int, HashSet<BaseModification>> kvp in oneBasedModifications)
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
                                                    }
                                                    else
                                                    {
                                                        Array.Reverse(sequence_array);
                                                        if (oneBasedModifications != null)
                                                        {
                                                            decoy_modifications = new Dictionary<int, HashSet<BaseModification>>(oneBasedModifications.Count);
                                                            foreach (KeyValuePair<int, HashSet<BaseModification>> kvp in oneBasedModifications)
                                                            {
                                                                decoy_modifications.Add(sequence.Length - kvp.Key + 1, kvp.Value);
                                                            }
                                                        }
                                                    }
                                                    var reversed_sequence = new string(sequence_array);
                                                    int[] decoybeginPositions = new int[oneBasedBeginPositions.Count];
                                                    int[] decoyendPositions = new int[oneBasedEndPositions.Count];
                                                    string[] decoyBigPeptideTypes = new string[oneBasedEndPositions.Count];
                                                    for (int i = 0; i < decoybeginPositions.Length; i++)
                                                    {
                                                        decoybeginPositions[oneBasedBeginPositions.Count - i - 1] = sequence.Length - oneBasedEndPositions[i] + 1;
                                                        decoyendPositions[oneBasedBeginPositions.Count - i - 1] = sequence.Length - oneBasedBeginPositions[i] + 1;
                                                        decoyBigPeptideTypes[oneBasedBeginPositions.Count - i - 1] = peptideTypes[i];
                                                    }
                                                    var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, decoy_modifications, decoybeginPositions, decoyendPositions, decoyBigPeptideTypes, name, full_name, offset, true, IsContaminant);
                                                    yield return decoy_protein;
                                                    offset += protein.Length;
                                                }
                                            }
                                            accession = null;
                                            name = null;
                                            full_name = null;
                                            sequence = null;
                                            feature_type = null;
                                            feature_description = null;
                                            oneBasedfeature_position = -1;
                                            oneBasedModifications = new Dictionary<int, HashSet<BaseModification>>();

                                            oneBasedBeginPositions = new List<int>();
                                            oneBasedEndPositions = new List<int>();
                                            peptideTypes = new List<string>();
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
                                accession = "";
                            }

                            // new protein
                            sequence = "";
                        }
                        else
                        {
                            sequence += line.Trim();
                        }

                        if (fasta.Peek() == '>' || fasta.Peek() == -1)
                        {
                            if (accession != null && sequence != null)
                            {
                                var protein = new Protein(sequence, accession, null, null, null, null, name, full_name, offset, false, IsContaminant);
                                yield return protein;

                                if (onTheFlyDecoys)
                                {
                                    char[] sequence_array = sequence.ToCharArray();
                                    if (sequence.StartsWith("M", StringComparison.InvariantCulture))
                                    {
                                        // Do not include the initiator methionine in reversal!!!
                                        Array.Reverse(sequence_array, 1, sequence.Length - 1);
                                    }
                                    else
                                    {
                                        Array.Reverse(sequence_array);
                                    }
                                    var reversed_sequence = new string(sequence_array);
                                    var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, null, null, null, null, name, full_name, offset, true, IsContaminant);
                                    yield return decoy_protein;
                                    offset += protein.Length;
                                }
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
        }

        #endregion Public Methods
    }
}