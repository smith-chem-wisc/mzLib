using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Xml;

namespace UsefulProteomicsDatabases
{
    public class ProteinDbWriter
    {

        #region Public Methods

        /// <summary>
        /// Writes a protein database in mzLibProteinDb format, with additional modifications from the Mods list.
        /// </summary>
        /// <param name="Mods"></param>
        /// <param name="proteinList"></param>
        /// <param name="outputFileName"></param>
        /// <returns>The new "modified residue" entries that are added due to being in the Mods dictionary</returns>
        public static Dictionary<string, int> WriteXmlDatabase(Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>> Mods, List<Protein> proteinList, string outputFileName)
        {
            var xmlWriterSettings = new XmlWriterSettings
            {
                Indent = true,
                IndentChars = "  "
            };

            Dictionary<string, int> newModResEntries = new Dictionary<string, int>();

            using (XmlWriter writer = XmlWriter.Create(outputFileName, xmlWriterSettings))
            {
                writer.WriteStartDocument();
                writer.WriteStartElement("mzLibProteinDb");

                HashSet<Modification> all_relevant_modifications = new HashSet<Modification>(
                    Mods.Where(kv => proteinList.Select(p => p.Accession).Contains(kv.Key))
                    .SelectMany(kv => kv.Value.Select(v => v.Item2))
                    .Concat(proteinList.SelectMany(p => p.OneBasedPossibleLocalizedModifications.Values.SelectMany(list => list))));

                foreach (Modification mod in all_relevant_modifications.OrderBy(m => m.id))
                {
                    writer.WriteStartElement("modification");
                    writer.WriteString(mod.ToString() + Environment.NewLine + "//");
                    writer.WriteEndElement();
                }
                foreach (Protein protein in proteinList)
                {
                    writer.WriteStartElement("entry");
                    writer.WriteStartElement("accession");
                    writer.WriteString(protein.Accession);
                    writer.WriteEndElement();
                    writer.WriteStartElement("name");
                    writer.WriteString(protein.Name);
                    writer.WriteEndElement();

                    writer.WriteStartElement("protein");
                    writer.WriteStartElement("recommendedName");
                    writer.WriteStartElement("fullName");
                    writer.WriteString(protein.FullName);
                    writer.WriteEndElement();
                    writer.WriteEndElement();
                    writer.WriteEndElement();

                    writer.WriteStartElement("gene");
                    foreach (var gene_name in protein.GeneNames)
                    {
                        writer.WriteStartElement("name");
                        writer.WriteAttributeString("type", gene_name.Item1);
                        writer.WriteString(gene_name.Item2);
                        writer.WriteEndElement();
                    }
                    writer.WriteEndElement();

                    foreach (var dbRef in protein.DatabaseReferences)
                    {
                        writer.WriteStartElement("dbReference");
                        writer.WriteAttributeString("type", dbRef.Type);
                        writer.WriteAttributeString("id", dbRef.Id);
                        foreach (Tuple<string, string> property in dbRef.Properties)
                        {
                            writer.WriteStartElement("property");
                            writer.WriteAttributeString("type", property.Item1);
                            writer.WriteAttributeString("value", property.Item2);
                            writer.WriteEndElement();
                        }
                        writer.WriteEndElement();
                    }
                    foreach (var proteolysisProduct in protein.ProteolysisProducts)
                    {
                        writer.WriteStartElement("feature");
                        writer.WriteAttributeString("type", proteolysisProduct.Type);
                        writer.WriteStartElement("location");
                        writer.WriteStartElement("begin");
                        writer.WriteAttributeString("position", proteolysisProduct.OneBasedBeginPosition.ToString());
                        writer.WriteEndElement();
                        writer.WriteStartElement("end");
                        writer.WriteAttributeString("position", proteolysisProduct.OneBasedEndPosition.ToString());
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }

                    Dictionary<int, HashSet<string>> modsToWrite = new Dictionary<int, HashSet<string>>();

                    foreach (var ye in protein.OneBasedPossibleLocalizedModifications)
                    {
                        foreach (var nice in ye.Value)
                        {
                            if (modsToWrite.TryGetValue(ye.Key, out HashSet<string> val))
                                val.Add(nice.id);
                            else
                                modsToWrite.Add(ye.Key, new HashSet<string> { nice.id });
                        }
                    }

                    if (Mods.ContainsKey(protein.Accession))
                        foreach (var ye in Mods[protein.Accession])
                        {
                            bool modAdded = false;
                            if (modsToWrite.TryGetValue(ye.Item1, out HashSet<string> val))
                                modAdded = val.Add(ye.Item2.id);
                            else
                            {
                                modsToWrite.Add(ye.Item1, new HashSet<string> { ye.Item2.id });
                                modAdded = true;
                            }
                            if (modAdded)
                            {
                                if (newModResEntries.ContainsKey(ye.Item2.id))
                                    newModResEntries[ye.Item2.id]++;
                                else
                                    newModResEntries.Add(ye.Item2.id, 1);
                            }
                        }

                    foreach (var hm in modsToWrite.OrderBy(b => b.Key))
                    {
                        foreach (var modId in hm.Value)
                        {
                            writer.WriteStartElement("feature");
                            writer.WriteAttributeString("type", "modified residue");
                            writer.WriteAttributeString("description", modId);
                            writer.WriteStartElement("location");
                            writer.WriteStartElement("position");
                            writer.WriteAttributeString("position", hm.Key.ToString(CultureInfo.InvariantCulture));
                            writer.WriteEndElement();
                            writer.WriteEndElement();
                            writer.WriteEndElement();
                        }
                    }

                    writer.WriteStartElement("sequence");
                    writer.WriteAttributeString("length", protein.Length.ToString(CultureInfo.InvariantCulture));
                    writer.WriteString(protein.BaseSequence);
                    writer.WriteEndElement();

                    writer.WriteEndElement();
                }

                writer.WriteEndElement();
                writer.WriteEndDocument();
            }
            return newModResEntries;
        }

        public static void WriteFastaDatabase(List<Protein> proteinList, string outputFileName, string delimeter)
        {
            using (StreamWriter writer = new StreamWriter(outputFileName))
            {
                foreach (Protein protein in proteinList)
                {
                    string header = protein.FullName != protein.Accession ?
                        protein.Accession + delimeter + protein.FullName :
                        header = protein.Accession; // we read in full name and accession to be the same string if the format isn't recognized
                    writer.WriteLine(">" + header);
                    writer.WriteLine(protein.BaseSequence);
                }
            }
        }

        #endregion Public Methods

    }
}