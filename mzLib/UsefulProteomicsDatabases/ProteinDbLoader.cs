using Ionic.Zlib;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Xml;
using System.Xml.Linq;
using System.Linq;

namespace UsefulProteomicsDatabases
{
    public static class ProteinDbLoader
    {

        #region Public Methods

        public static List<Protein> LoadProteinDb(string proteinDbLocation, bool onTheFlyDecoys, IDictionary<string, IList<Modification>> allKnownModifications, bool IsContaminant, out Dictionary<string, Modification> unknownModifications)
        {
            List<Protein> result = new List<Protein>();
            unknownModifications = new Dictionary<string, Modification>();
            using (var stream = new FileStream(proteinDbLocation, FileMode.Open))
            {
<<<<<<< HEAD
                string accession = null;
                string name = null;
                string full_name = null;
                int offset = 0;

                var oneBasedBeginPositions = new List<int?>();
                var oneBasedEndPositions = new List<int?>();
                var peptideTypes = new List<string>();
                var oneBasedModifications = new Dictionary<int, List<Modification>>();

=======
>>>>>>> ee91ea2... Added GoTerms and changed LoadDb over to using XElement objects
                // xml db
                if (!proteinDbLocation.EndsWith(".fasta"))
                {
                    Stream uniprotXmlFileStream = stream;
                    if (proteinDbLocation.EndsWith(".gz"))
                        uniprotXmlFileStream = new GZipStream(stream, CompressionMode.Decompress);

                    string sequence = null;
<<<<<<< HEAD
                    string feature_type = null;
                    string feature_description = null;
                    int oneBasedfeature_position = -1;
                    int? oneBasedbeginPosition = null;
                    int? oneBasedendPosition = null;
=======
                    var oneBasedBeginPositions = new List<int?>();
                    var oneBasedEndPositions = new List<int?>();
                    var peptideTypes = new List<string>();
                    int offset = 0;
>>>>>>> ee91ea2... Added GoTerms and changed LoadDb over to using XElement objects

                    List<XElement> entries = new List<XElement>();
                    using (XmlReader uniprotXmlReader = XmlReader.Create(uniprotXmlFileStream))
                    {
                        uniprotXmlReader.MoveToContent();
                        while (uniprotXmlReader.Read())
                        {
                            if (uniprotXmlReader.NodeType == XmlNodeType.Element && uniprotXmlReader.Name == "entry")
                                entries.Add(XElement.ReadFrom(uniprotXmlReader) as XElement);
                        }
                    }

                    foreach (XElement entry in entries)
                    {
                        //Used fields
                        string accession = GetChild(entry, "accession").Value;
                        string full_name = GetDescendant(entry, "fullName").Value;
                        List<XElement> features = entry.Elements().Where(node => node.Name.LocalName == "feature").ToList();
                        List<XElement> dbReferences = entry.Elements().Where(node => node.Name.LocalName == "dbReference").ToList();
                        List<GoTerm> goTerms = new List<GoTerm>();
                        XElement sequence_elem = GetChild(entry, "sequence");
                        sequence = sequence_elem.Value.Replace("\r", null).Replace("\n", null);
                        Dictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>();

                        //Other fields
                        //Full name follows desired order: recommendedName > submittedName > alternativeName because these appear in this order in UniProt-XMLs
                        string name = GetChild(entry, "name").Value;
                        string organism = GetChild(GetChild(entry, "organism"), "name").Value;
                        string gene_name = GetChild(GetChild(entry, "gene"), "name").Value;

                        //Unused fields
                        //string dataset = GetAttribute(entry, "dataset");
                        //string fragment = GetAttribute(sequence_elem, "fragment");
                        //string protein_existence = GetAttribute(GetChild(entry, "proteinExistence"), "type");
                        //string sequence_version = GetAttribute(sequence_elem, "version");

                        //Process dbReferences to retrieve Gene Ontology Terms
                        foreach (XElement dbReference in dbReferences)
                        {
                            string dbReference_type = GetAttribute(dbReference, "type");
                            if (dbReference_type == "GO")
                            {
                                GoTerm go = new GoTerm();
                                string ID = GetAttribute(dbReference, "id");
                                go.id = ID.Split(':')[1].ToString();
                                IEnumerable<XElement> dbProperties = from XElement in dbReference.Elements() where XElement.Name.LocalName == "property" select XElement;

                                foreach (XElement property in dbProperties)
                                {
                                    string type = GetAttribute(property, "type");
                                    if (type == "term")
                                    {
                                        string description = GetAttribute(property, "value");
                                        switch (description.Split(':')[0].ToString())
                                        {
                                            case "C":
                                                go.aspect = Aspect.cellularComponent;
                                                go.description = description.Split(':')[1].ToString();
                                                break;
                                            case "F":
                                                go.aspect = Aspect.molecularFunction;
                                                go.description = description.Split(':')[1].ToString();
                                                break;
                                            case "P":
                                                go.aspect = Aspect.biologicalProcess;
                                                go.description = description.Split(':')[1].ToString();
                                                break;
                                        }
                                        goTerms.Add(go);
                                    }
                                }
                            }
                        }

                        //Process the modified residues
                        foreach (XElement feature in features)
                        {
                            string feature_type = GetAttribute(feature, "type");
                            switch (feature_type)
                            {
                                case "modified residue":
                                    XElement feature_position_elem = GetDescendant(feature, "position");
                                    int feature_position = ConvertPositionElem(feature_position_elem);
                                    string feature_description = GetAttribute(feature, "description").Split(';')[0];
                                    if (feature_position >= 0 && allKnownModifications.ContainsKey(feature_description))
                                    {
                                        List<Modification> modListAtPos;
                                        if (oneBasedModifications.TryGetValue(feature_position, out modListAtPos))
                                            modListAtPos.AddRange(allKnownModifications[feature_description]);
                                        else
                                        {
                                            modListAtPos = new List<Modification>(allKnownModifications[feature_description]);
                                            oneBasedModifications.Add(feature_position, modListAtPos);
                                        }
                                    }
                                    break;
                                case "chain":
                                case "signal peptide":
                                case "propeptide":
                                case "peptide":
                                    XElement feature_begin_elem = GetDescendant(feature, "begin");
                                    XElement feature_end_elem = GetDescendant(feature, "end");
                                    int feature_begin = ConvertPositionElem(feature_begin_elem);
                                    int feature_end = ConvertPositionElem(feature_end_elem);
                                    peptideTypes.Add(feature_type);
                                    oneBasedBeginPositions.Add(feature_begin);
                                    oneBasedEndPositions.Add(feature_end);
                                    break;
                                case "splice variant":
                                case "sequence variant":
                                    break;
                            }
                        }

                        //Add the full length protein, and then add the fragments with segments of the above modification dictionary
                        var protein = new Protein(sequence, accession, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, offset, false, IsContaminant, goTerms);
                        result.Add(protein);
                        offset += protein.Length;

                        if (onTheFlyDecoys)
                        {
                            char[] sequence_array = sequence.ToCharArray();
                            Dictionary<int, List<Modification>> decoy_modifications = null;
                            int starts_with_met = Convert.ToInt32(sequence.StartsWith("M", StringComparison.InvariantCulture));
                            Array.Reverse(sequence_array, starts_with_met, sequence_array.Length - starts_with_met); // Do not include the initiator methionine in reversal!!!
                            if (oneBasedModifications != null)
                            {
                                decoy_modifications = new Dictionary<int, List<Modification>>(oneBasedModifications.Count);
                                foreach (var kvp in oneBasedModifications)
                                {
                                    if (starts_with_met == 1 && kvp.Key == 1)
                                    {
                                        decoy_modifications.Add(1, kvp.Value);
                                    }
                                    else if (starts_with_met == 0 || kvp.Key > 1)
                                    {
                                        decoy_modifications.Add(sequence.Length - kvp.Key + 1 + starts_with_met, kvp.Value);
                                    }
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
                            var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, decoy_modifications, decoybeginPositions, decoyendPositions, decoyBigPeptideTypes, name, full_name, offset, true, IsContaminant, null);

                            result.Add(decoy_protein);
                            offset += decoy_protein.Length;
                        }
                    }
                }

                // fasta db
                else
                {
                    string accession = null;
                    string name = null;
                    string full_name = null;
                    int offset = 0;

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
                            var sequence = sb.ToString();
<<<<<<< HEAD
                            var protein = new Protein(sequence, accession, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, offset, false, IsContaminant);
=======
                            var protein = new Protein(sequence, accession, null, null, null, null, name, full_name, offset, false, IsContaminant, null);
>>>>>>> ee91ea2... Added GoTerms and changed LoadDb over to using XElement objects
                            offset += protein.Length;

                            result.Add(protein);

                            if (onTheFlyDecoys)
                            {
                                char[] sequence_array = sequence.ToCharArray();
                                int starts_with_met = Convert.ToInt32(sequence.StartsWith("M", StringComparison.InvariantCulture));
                                Array.Reverse(sequence_array, starts_with_met, sequence.Length - starts_with_met); // Do not include the initiator methionine in reversal!!!
                                var reversed_sequence = new string(sequence_array);
<<<<<<< HEAD
                                var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, offset, true, IsContaminant);
=======
                                var decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, null, null, null, null, name, full_name, offset, true, IsContaminant, null);
>>>>>>> ee91ea2... Added GoTerms and changed LoadDb over to using XElement objects

                                result.Add(decoy_protein);
                                offset += protein.Length;
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

        #region Private Methods

        private static string GetAttribute(XElement element, string attribute_name)
        {
            XAttribute attribute = element.Attributes().FirstOrDefault(a => a.Name.LocalName == attribute_name);
            return attribute == null ? "" : attribute.Value;
        }

        private static XElement GetChild(XElement element, string name)
        {
            XElement e = element.Elements().FirstOrDefault(elem => elem.Name.LocalName == name);
            if (e != null) return e;
            else return new XElement("dummy_node");
        }

        private static XElement GetDescendant(XElement element, string name)
        {
            XElement e = element.Descendants().FirstOrDefault(elem => elem.Name.LocalName == name);
            if (e != null) return e;
            else return new XElement("dummy_node");
        }

        private static int ConvertPositionElem(XElement position_elem)
        {
            string feature_position = GetAttribute(position_elem, "position");
            string feature_position_status = GetAttribute(position_elem, "status"); //positionType elements have default 'status' of certain
            if (feature_position != "" && feature_position_status == "")
                return Convert.ToInt32(feature_position) - 1;
            else
                return -1;
        }

        #endregion Private Methods
    }
}