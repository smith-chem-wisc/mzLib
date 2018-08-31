using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System.Xml;

namespace UsefulProteomicsDatabases
{
    public class ProteinXmlEntry
    {
        private Regex SubstituteWhitespace = new Regex(@"\s+");

        public string Accession { get; private set; }
        public string Name { get; private set; }
        public string FullName { get; private set; }
        public string Organism { get; private set; }
        public string Sequence { get; private set; }
        public string FeatureType { get; private set; }
        public string FeatureDescription { get; private set; }
        public string OriginalValue { get; private set; } = ""; // if no content is found, assume it is empty, not null (e.g. <original>A</original><variation/> for a deletion event)
        public string VariationValue { get; private set; } = "";
        public string DBReferenceType { get; private set; }
        public string DBReferenceId { get; private set; }
        public List<string> PropertyTypes { get; private set; } = new List<string>();
        public List<string> PropertyValues { get; private set; } = new List<string>();
        public int OneBasedFeaturePosition { get; private set; } = -1;
        public int? OneBasedBeginPosition { get; private set; }
        public int? OneBasedEndPosition { get; private set; }
        public List<ProteolysisProduct> ProteolysisProducts { get; private set; } = new List<ProteolysisProduct>();
        public List<SequenceVariation> SequenceVariations { get; private set; } = new List<SequenceVariation>();
        public List<DisulfideBond> DisulfideBonds { get; private set; } = new List<DisulfideBond>();
        public Dictionary<int, List<Modification>> OneBasedModifications { get; private set; } = new Dictionary<int, List<Modification>>();
        public List<Tuple<string, string>> GeneNames { get; private set; } = new List<Tuple<string, string>>();
        public List<DatabaseReference> DatabaseReferences { get; private set; } = new List<DatabaseReference>();
        public bool ReadingGene { get; set; }
        public bool ReadingOrganism { get; set; }

        private List<(int, string)> AnnotatedMods { get; set; } = new List<(int position, string originalModificationID)>();

        /// <summary>
        /// Start parsing a protein XML element
        /// </summary>
        public void ParseElement(string elementName, XmlReader xml)
        {
            int outValue;
            switch (elementName)
            {
                case "accession":
                    if (Accession == null)
                    {
                        Accession = xml.ReadElementString();
                    }
                    break;

                case "name":
                    if (xml.Depth == 2 && !ReadingGene && !ReadingOrganism)
                    {
                        Name = xml.ReadElementString();
                    }
                    if (ReadingGene && !ReadingOrganism)
                    {
                        GeneNames.Add(new Tuple<string, string>(xml.GetAttribute("type"), xml.ReadElementString()));
                    }
                    if (ReadingOrganism)
                    {
                        if (xml.GetAttribute("type").Equals("scientific"))
                        {
                            Organism = xml.ReadElementString();
                        }
                    }
                    break;

                case "gene":
                    ReadingGene = true;
                    break;

                case "organism":
                    if (Organism == null)
                    {
                        ReadingOrganism = true;
                    }
                    break;

                case "fullName":
                    if (FullName == null)
                    {
                        FullName = xml.ReadElementString();
                    }
                    break;

                case "feature":
                    FeatureType = xml.GetAttribute("type");
                    FeatureDescription = xml.GetAttribute("description");
                    break;

                case "original":
                    OriginalValue = xml.ReadElementString();
                    break;

                case "variation":
                    VariationValue = xml.ReadElementString();
                    break;

                case "dbReference":
                    PropertyTypes.Clear();
                    PropertyValues.Clear();
                    DBReferenceType = xml.GetAttribute("type");
                    DBReferenceId = xml.GetAttribute("id");
                    break;

                case "property":
                    PropertyTypes.Add(xml.GetAttribute("type"));
                    PropertyValues.Add(xml.GetAttribute("value"));
                    break;

                case "position":
                    OneBasedFeaturePosition = int.Parse(xml.GetAttribute("position"));
                    break;

                case "begin":
                    OneBasedBeginPosition = int.TryParse(xml.GetAttribute("position"), out outValue) ? (int?)outValue : null;
                    break;

                case "end":
                    OneBasedEndPosition = int.TryParse(xml.GetAttribute("position"), out outValue) ? (int?)outValue : null;
                    break;

                case "sequence":
                    Sequence = SubstituteWhitespace.Replace(xml.ReadElementString(), "");
                    break;
            }
        }

        /// <summary>
        /// Finish parsing at the end of an element
        /// </summary>
        public List<Protein> ParseEndElement(XmlReader xml, IEnumerable<string> modTypesToExclude,
            Dictionary<string, Modification> unknownModifications, bool isContaminant, string proteinDbLocation)
        {
            List<Protein> result = new List<Protein>();
            if (xml.Name == "feature")
            {
                ParseFeatureEndElement(xml, modTypesToExclude, unknownModifications);
            }
            else if (xml.Name == "dbReference")
            {
                ParseDatabaseReferenceEndElement(xml);
            }
            else if (xml.Name == "gene")
            {
                ReadingGene = false;
            }
            else if (xml.Name == "organism")
            {
                ReadingOrganism = false;
            }
            else if (xml.Name == "entry")
            {
                result.AddRange(ParseEntryEndElement(xml, isContaminant, proteinDbLocation, modTypesToExclude, unknownModifications));
            }
            return result;
        }

        /// <summary>
        /// Finish parsing an entry
        /// </summary>
        public List<Protein> ParseEntryEndElement(XmlReader xml, bool isContaminant, string proteinDbLocation, IEnumerable<string> modTypesToExclude, Dictionary<string, Modification> unknownModifications)
        {
            List<Protein> result = new List<Protein>();
            if (Accession != null && Sequence != null)
            {
                ParseAnnotatedMods(modTypesToExclude, unknownModifications);
                var protein = new Protein(Sequence, Accession, Organism, GeneNames, OneBasedModifications, ProteolysisProducts, Name, FullName,
                    false, isContaminant, DatabaseReferences, SequenceVariations, DisulfideBonds, proteinDbLocation);
                result.Add(protein);
            }
            Clear();
            return result;
        }

        /// <summary>
        /// Finish parsing a feature element
        /// </summary>
        public void ParseFeatureEndElement(XmlReader xml, IEnumerable<string> modTypesToExclude,
            Dictionary<string, Modification> unknownModifications)
        {
            if (FeatureType == "modified residue")
            {
                FeatureDescription = FeatureDescription.Split(';')[0];
                AnnotatedMods.Add((OneBasedFeaturePosition, FeatureDescription));
            }
            else if (FeatureType == "peptide" || FeatureType == "propeptide" || FeatureType == "chain" || FeatureType == "signal peptide")
            {
                ProteolysisProducts.Add(new ProteolysisProduct(OneBasedBeginPosition, OneBasedEndPosition, FeatureType));
            }
            else if (FeatureType == "sequence variant" && VariationValue != null && VariationValue != "") // Only keep if there is variant sequence information and position information
            {
                if (OneBasedBeginPosition != null && OneBasedEndPosition != null)
                {
                    SequenceVariations.Add(new SequenceVariation((int)OneBasedBeginPosition, (int)OneBasedEndPosition, OriginalValue, VariationValue, FeatureDescription));
                }
                else if (OneBasedFeaturePosition >= 1)
                {
                    SequenceVariations.Add(new SequenceVariation(OneBasedFeaturePosition, OriginalValue, VariationValue, FeatureDescription));
                }
            }
            else if (FeatureType == "disulfide bond")
            {
                if (OneBasedBeginPosition != null && OneBasedEndPosition != null)
                {
                    DisulfideBonds.Add(new DisulfideBond((int)OneBasedBeginPosition, (int)OneBasedEndPosition, FeatureDescription));
                }
                else if (OneBasedFeaturePosition >= 1)
                {
                    DisulfideBonds.Add(new DisulfideBond(OneBasedFeaturePosition, FeatureDescription));
                }
            }
            OneBasedBeginPosition = null;
            OneBasedEndPosition = null;
            OneBasedFeaturePosition = -1;
            OriginalValue = "";
            VariationValue = "";
        }

        private void ParseAnnotatedMods(IEnumerable<string> modTypesToExclude, Dictionary<string, Modification> unknownModifications)
        {
            foreach ((int, string) annotatedMod in AnnotatedMods)
            {
                string annotatedId = annotatedMod.Item2;
                int annotatedModLocation = annotatedMod.Item1;

                if (ProteinDbLoader.IdWithMotifToMod.TryGetValue(annotatedId, out Modification foundMod))
                {
                    // if the list of known mods contains this IdWithMotif
                    if (ModificationLocalization.ModFits(foundMod, Sequence, 0, 0, annotatedModLocation)
                        && !modTypesToExclude.Contains(foundMod.ModificationType))
                    {
                        if (OneBasedModifications.TryGetValue(annotatedModLocation, out var listOfModsAtThisLocation))
                        {
                            listOfModsAtThisLocation.Add(foundMod);
                        }
                        else
                        {
                            OneBasedModifications.Add(annotatedModLocation, new List<Modification> { foundMod });
                        }
                    }
                    // else - the mod ID was found but the motif didn't fit the annotated location
                }

                // no known mod - try looking it up in the dictionary of mods without motif appended
                else if (ProteinDbLoader.IdToPossibleMods.TryGetValue(annotatedId, out IList<Modification> mods))
                {
                    foreach (Modification mod in mods)
                    {
                        if (ModificationLocalization.ModFits(mod, Sequence, 0, 0, annotatedModLocation)
                            && !modTypesToExclude.Contains(mod.ModificationType))
                        {
                            if (OneBasedModifications.TryGetValue(annotatedModLocation, out var listOfModsAtThisLocation))
                            {
                                listOfModsAtThisLocation.Add(mod);
                            }
                            else
                            {
                                OneBasedModifications.Add(annotatedModLocation, new List<Modification> { mod });
                            }
                            break;
                        }
                    }
                }
                else
                {
                    // could not find the annotated mod's ID in our list of known mods - it's an unknown mod
                    // I don't think this really does anything...
                    if (!unknownModifications.ContainsKey(annotatedId))
                    {
                        unknownModifications.Add(annotatedId, new Modification(annotatedId));
                    }
                }

                //}
                //foreach ((int, string) modifiedResidueEntryFromXml in temp)
                //{
                //    ModificationMotif motif = GetMotif(Sequence, modifiedResidueEntryFromXml.Item1);
                //    if (modifiedResidueEntryFromXml.Item2.Contains(" on "))
                //    {
                //        //ID probably contains the motif
                //        if (modifiedResidueEntryFromXml.Item2.Contains(" on " + motif.ToString()) || modifiedResidueEntryFromXml.Item2.Contains(" on X"))
                //        {
                //            // Create new entry for this residue, if needed
                //            if (!OneBasedModifications.TryGetValue(modifiedResidueEntryFromXml.Item1, out List<Modification> residue_modifications))
                //            {
                //                residue_modifications = new List<Modification>();
                //                OneBasedModifications.Add(modifiedResidueEntryFromXml.Item1, residue_modifications);
                //            }

                //            //TODO deal with " on X"
                //            if (mod_dict.ContainsKey(modifiedResidueEntryFromXml.Item2))
                //            {
                //                // Known and not of a type in the exclusion list
                //                List<Modification> mods = mod_dict[modifiedResidueEntryFromXml.Item2].Where(m => !modTypesToExclude.Contains(m.ModificationType)).ToList();
                //                if (mods.Count == 0 && OneBasedModifications[modifiedResidueEntryFromXml.Item1].Count == 0)
                //                {
                //                    OneBasedModifications.Remove(modifiedResidueEntryFromXml.Item1);
                //                }
                //                else
                //                {
                //                    OneBasedModifications[modifiedResidueEntryFromXml.Item1].AddRange(mods);
                //                }
                //            }
                //            else if (unknownModifications.ContainsKey(modifiedResidueEntryFromXml.Item2))
                //            {
                //                // Not known but seen
                //                residue_modifications.Add(unknownModifications[modifiedResidueEntryFromXml.Item2]);
                //            }
                //            else
                //            {
                //                // Not known and not seen
                //                //we don't know what this is so we're not doing anything with the mod at this point.
                //                continue;

                //                //unknownModifications[FeatureDescription] = new Modification(_originalId: FeatureDescription, _modificationType: "unknown", _target: motif);
                //                //residue_modifications.Add(unknownModifications[modifiedResidueEntryFromXml.Item2]);
                //            }
                //        }
                //    }
                //    else
                //    {
                //        //modification doesn't contain the motif so let's add it to the name and see if we can find it.
                //        string featureName = modifiedResidueEntryFromXml.Item2 + " on " + motif.ToString();
                //        string featureNameX = modifiedResidueEntryFromXml.Item2 + " on X";

                //        // Create new entry for this residue, if needed
                //        if (!OneBasedModifications.TryGetValue(modifiedResidueEntryFromXml.Item1, out List<Modification> residue_modifications))
                //        {
                //            residue_modifications = new List<Modification>();
                //            OneBasedModifications.Add(modifiedResidueEntryFromXml.Item1, residue_modifications);
                //        }
                //        if (mod_dict.ContainsKey(featureName) || mod_dict.ContainsKey(featureNameX))
                //        {
                //            List<Modification> mods = new List<Modification>();
                //            // Known and not of a type in the exclusion list

                //            //We look for the modification on the specified amino acid motif first
                //            if (mod_dict.ContainsKey(featureName))
                //            {
                //                mods = mod_dict[featureName].Where(m => !modTypesToExclude.Contains(m.ModificationType)).ToList();
                //            }

                //            //then if we don't find it there, we go up a level and see if that mod is described as " on X"
                //            else
                //            {
                //                mods = mod_dict[featureNameX].Where(m => !modTypesToExclude.Contains(m.ModificationType)).ToList();
                //            }
                //            if (mods.Count == 0 && OneBasedModifications[modifiedResidueEntryFromXml.Item1].Count == 0)
                //            {
                //                OneBasedModifications.Remove(modifiedResidueEntryFromXml.Item1);
                //            }
                //            else
                //            {
                //                OneBasedModifications[modifiedResidueEntryFromXml.Item1].AddRange(mods);
                //            }
                //        }
                //        else if (unknownModifications.ContainsKey(featureName))
                //        {
                //            //not sure what this is about. but since we're not adding mods we don't know anything about...
                //            continue;
                //            //// Not known but seen
                //            //residue_modifications.Add(unknownModifications[featureName]);
                //        }
                //        else
                //        {
                //            //we're not adding anything if we can't figure it out.
                //            continue;
                //            // Not known and not seen
                //            //unknownModifications[FeatureDescription] = new Modification(_originalId: FeatureDescription, _modificationType: "unknown", _target: motif);
                //            //residue_modifications.Add(unknownModifications[featureName]);
                //        }
                //    }
            }
        }

        private static ModificationMotif GetMotif(string proteinSequence, int position)
        {
            string aminoAcid = proteinSequence.Substring(position - 1, 1);
            if (ModificationMotif.TryGetMotif(aminoAcid, out ModificationMotif motif))
            {
                return motif;
            }
            else
            {
                return null;
            }
        }

        /// <summary>
        /// Finish parsing a database reference element
        /// </summary>
        /// <param name="xml"></param>
        public void ParseDatabaseReferenceEndElement(XmlReader xml)
        {
            DatabaseReferences.Add(
                new DatabaseReference(DBReferenceType, DBReferenceId,
                    Enumerable.Range(0, PropertyTypes.Count).Select(i => new Tuple<string, string>(PropertyTypes[i], PropertyValues[i])).ToList()));
            PropertyTypes = new List<string>();
            PropertyValues = new List<string>();
            DBReferenceType = null;
            DBReferenceId = null;
        }

        /// <summary>
        /// Clear this object's properties
        /// </summary>
        private void Clear()
        {
            Accession = null;
            Name = null;
            FullName = null;
            Sequence = null;
            Organism = null;
            FeatureType = null;
            FeatureDescription = null;
            OriginalValue = "";
            VariationValue = "";
            DBReferenceType = null;
            DBReferenceId = null;
            PropertyTypes = new List<string>();
            PropertyValues = new List<string>();
            OneBasedFeaturePosition = -1;
            AnnotatedMods = new List<(int, string)>();
            OneBasedModifications = new Dictionary<int, List<Modification>>();
            ProteolysisProducts = new List<ProteolysisProduct>();
            SequenceVariations = new List<SequenceVariation>();
            DisulfideBonds = new List<DisulfideBond>();
            DatabaseReferences = new List<DatabaseReference>();
            GeneNames = new List<Tuple<string, string>>();
            ReadingGene = false;
            ReadingOrganism = false;
        }
    }
}