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
        /// Writes a protein database in mzLibProteinDb format, with additional modifications from the AdditionalModsToAddToProteins list.
        /// </summary>
        /// <param name="AdditionalModsToAddToProteins"></param>
        /// <param name="proteinList"></param>
        /// <param name="outputFileName"></param>
        /// <returns>The new "modified residue" entries that are added due to being in the Mods dictionary</returns>
        public static Dictionary<string, int> WriteXmlDatabase(Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>> AdditionalModsToAddToProteins, List<Protein> proteinList, string outputFileName)
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
                    AdditionalModsToAddToProteins.Where(kv => proteinList.Select(p => p.Accession).Contains(kv.Key))
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

                    Dictionary<int, HashSet<string>> modsToWriteForThisSpecificProtein = new Dictionary<int, HashSet<string>>();

                    foreach (var ye in protein.OneBasedPossibleLocalizedModifications)
                    {
                        foreach (var nice in ye.Value)
                        {
                            if (modsToWriteForThisSpecificProtein.TryGetValue(ye.Key, out HashSet<string> val))
                                val.Add(nice.id);
                            else
                                modsToWriteForThisSpecificProtein.Add(ye.Key, new HashSet<string> { nice.id });
                        }
                    }

                    if (AdditionalModsToAddToProteins.ContainsKey(protein.Accession))
                        foreach (var ye in AdditionalModsToAddToProteins[protein.Accession])
                        {
                            int additionalModResidueIndex = ye.Item1;
                            string additionalModId = ye.Item2.id;
                            bool modAdded = false;

                            // If we already have modifications that need to be written to the specific residue, get the hash set of those mods
                            if (modsToWriteForThisSpecificProtein.TryGetValue(additionalModResidueIndex, out HashSet<string> val))
                                // Try to add the new mod to that hash set. If it's not there, modAdded=true, and it is added.
                                modAdded = val.Add(additionalModId);

                            // Otherwise, no modifications currently need to be written to the residue at residueIndex, so need to create new hash set for that residue
                            else
                            {
                                modsToWriteForThisSpecificProtein.Add(additionalModResidueIndex, new HashSet<string> { additionalModId });
                                modAdded = true;
                            }

                            // Finally, if a new modification has in fact been deemed worthy of being added to the database, mark that in the output dictionary
                            if (modAdded)
                            {
                                if (newModResEntries.ContainsKey(additionalModId))
                                    newModResEntries[additionalModId]++;
                                else
                                    newModResEntries.Add(additionalModId, 1);
                            }
                        }

                    foreach (var hm in modsToWriteForThisSpecificProtein.OrderBy(b => b.Key))
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

                    foreach (var hm in protein.SequenceVariations)
                    {
                        //Don't write if there is no position or sequence information
                        if (hm.OneBasedPosition < 1 && hm.OneBasedBeginPosition == null && hm.OneBasedEndPosition == null 
                            || (hm.OriginalSequence == null || hm.OriginalSequence == "") && (hm.VariantSequence == null || hm.VariantSequence == ""))
                            continue;

                        writer.WriteStartElement("feature");
                        writer.WriteAttributeString("type", "sequence variant");
                        writer.WriteAttributeString("description", hm.Description);
                        writer.WriteStartElement("original");
                        writer.WriteString(hm.OriginalSequence);
                        writer.WriteEndElement(); // original
                        writer.WriteStartElement("variation");
                        writer.WriteString(hm.VariantSequence);
                        writer.WriteEndElement(); // variation
                        writer.WriteStartElement("location");
                        if (hm.OneBasedPosition >= 1)
                        {
                            writer.WriteStartElement("position");
                            writer.WriteAttributeString("position", hm.OneBasedPosition.ToString());
                            writer.WriteEndElement();
                        }
                        else
                        {
                            writer.WriteStartElement("begin");
                            writer.WriteAttributeString("position", hm.OneBasedBeginPosition.ToString());
                            writer.WriteEndElement();
                            writer.WriteStartElement("end");
                            writer.WriteAttributeString("position", hm.OneBasedEndPosition.ToString());
                            writer.WriteEndElement();
                        }
                        writer.WriteEndElement(); // location
                        writer.WriteEndElement(); // feature
                    }

                    writer.WriteStartElement("sequence");
                    writer.WriteAttributeString("length", protein.Length.ToString(CultureInfo.InvariantCulture));
                    writer.WriteString(protein.BaseSequence);
                    writer.WriteEndElement(); // sequence
                    writer.WriteEndElement(); // entry
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