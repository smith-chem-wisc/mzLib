using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System.Xml;
using Omics.BioPolymer;
using Omics.Modifications;
using Transcriptomics;
using UsefulProteomicsDatabases.Transcriptomics;
using System.Data;
using Proteomics.ProteolyticDigestion;
using System.Diagnostics;

namespace UsefulProteomicsDatabases
{
    public class ProteinXmlEntry
    {
        private static readonly Regex SubstituteWhitespace = new Regex(@"\s+");
        public string DatasetEntryTag { get; private set; }
        public string DatabaseCreatedEntryTag { get; private set; }
        public string DatabaseModifiedEntryTag { get; private set; }
        public string DatabaseVersionEntryTag { get; private set; }
        public string XmlnsEntryTag { get; private set; }
        public string Accession { get; private set; }
        public string Name { get; private set; }
        public string FullName { get; private set; }
        public string Organism { get; private set; }
        public string Sequence { get; private set; }
        public string FeatureType { get; private set; }
        public string FeatureDescription { get; private set; }
        public string SubFeatureType { get; private set; }
        public string SubFeatureDescription { get; private set; }
        public string OriginalValue { get; private set; } = "";
        public string VariationValue { get; private set; } = "";
        public string DBReferenceType { get; private set; }
        public string DBReferenceId { get; private set; }
        public List<string> PropertyTypes { get; private set; } = new List<string>();
        public List<string> PropertyValues { get; private set; } = new List<string>();
        public int OneBasedFeaturePosition { get; private set; } = -1;
        public int OneBasedFeatureSubPosition { get; private set; } = -1;
        public int? OneBasedBeginPosition { get; private set; }
        public int? OneBasedEndPosition { get; private set; }
        public List<TruncationProduct> ProteolysisProducts { get; private set; } = new List<TruncationProduct>();
        public List<SequenceVariation> SequenceVariations { get; private set; } = new List<SequenceVariation>();
        public List<DisulfideBond> DisulfideBonds { get; private set; } = new List<DisulfideBond>();
        public List<SpliceSite> SpliceSites { get; private set; } = new List<SpliceSite>();
        public Dictionary<int, List<Modification>> OneBasedModifications { get; private set; } = new Dictionary<int, List<Modification>>();
        public Dictionary<int, List<Modification>> OneBasedVariantModifications { get; private set; } = new Dictionary<int, List<Modification>>();
        public List<Tuple<string, string>> GeneNames { get; private set; } = new List<Tuple<string, string>>();
        public List<DatabaseReference> DatabaseReferences { get; private set; } = new List<DatabaseReference>();
        public bool ReadingGene { get; set; }
        public bool ReadingOrganism { get; set; }
        public UniProtSequenceAttributes SequenceAttributes { get; set; } = null;
        private List<(int, string)> AnnotatedMods = new List<(int position, string originalModificationID)>();
        private List<(int, string)> AnnotatedVariantMods = new List<(int position, string originalModificationID)>();

        // Captured isoform/sequence identifier from <location sequence="...">
        private string LocationSequenceId;

        public void ParseElement(string elementName, XmlReader xml)
        {
            int outValue;
            switch (elementName)
            {
                case "entry":
                    ParseEntryAttributes(xml);
                    break;
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
                case "subfeature":
                    SubFeatureType = xml.GetAttribute("type");
                    SubFeatureDescription = xml.GetAttribute("description");
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
                case "location":
                    LocationSequenceId = xml.GetAttribute("sequence");
                    break;
                case "position":
                    OneBasedFeaturePosition = int.Parse(xml.GetAttribute("position"));
                    break;
                case "subposition":
                    OneBasedFeatureSubPosition = int.Parse(xml.GetAttribute("subposition"));
                    break;
                case "begin":
                    OneBasedBeginPosition = int.TryParse(xml.GetAttribute("position"), out outValue) ? (int?)outValue : null;
                    break;
                case "end":
                    OneBasedEndPosition = int.TryParse(xml.GetAttribute("position"), out outValue) ? (int?)outValue : null;
                    break;
                case "sequence":
                    ParseSequenceAttributes(xml);
                    break;
            }
        }

        private void ParseEntryAttributes(XmlReader xml)
        {
            DatasetEntryTag = xml.GetAttribute("dataset");
            DatabaseCreatedEntryTag = xml.GetAttribute("created");
            DatabaseModifiedEntryTag = xml.GetAttribute("modified");
            DatabaseVersionEntryTag = xml.GetAttribute("version");
            XmlnsEntryTag = xml.GetAttribute("xmlns");
        }

        private void ParseSequenceAttributes(XmlReader xml)
        {
            string checksumAttr = xml.GetAttribute("checksum");
            string modifiedAttr = xml.GetAttribute("modified");
            string sequenceVersionAttribute = xml.GetAttribute("version");
            string precursorAttr = xml.GetAttribute("precursor");
            string fragmentAttrString = xml.GetAttribute("fragment");

            string checksum = string.IsNullOrEmpty(checksumAttr) ? "" : checksumAttr;
            DateTime entryModified = ParseModifiedDate(modifiedAttr);
            int sequenceVersion = ParseSequenceVersion(sequenceVersionAttribute);
            bool isPrecursor = ParseIsPrecursor(precursorAttr);
            UniProtSequenceAttributes.FragmentType fragment = ParseFragmentType(fragmentAttrString);

            Sequence = SubstituteWhitespace.Replace(xml.ReadElementString(), "");
            int length = Sequence.Length;
            int mass = ComputeSequenceMass(Sequence);

            SequenceAttributes = new UniProtSequenceAttributes(length, mass, checksum, entryModified, sequenceVersion, isPrecursor, fragment);
        }

        private static DateTime ParseModifiedDate(string modifiedAttr)
        {
            if (!string.IsNullOrEmpty(modifiedAttr))
            {
                try
                {
                    return DateTime.ParseExact(modifiedAttr, "yyyy-MM-dd", System.Globalization.CultureInfo.InvariantCulture);
                }
                catch
                {
                    Trace.TraceWarning($"Warning: Failed to parse modified date '{modifiedAttr}'. Using DateTime.Now.");
                }
            }
            return DateTime.Now;
        }

        private static int ParseSequenceVersion(string versionAttr)
        {
            if (int.TryParse(versionAttr, out int version))
            {
                return version;
            }
            return -1;
        }

        private static bool ParseIsPrecursor(string precursorAttr)
        {
            return !string.IsNullOrEmpty(precursorAttr) && precursorAttr.Equals("true", StringComparison.OrdinalIgnoreCase);
        }

        private static UniProtSequenceAttributes.FragmentType ParseFragmentType(string fragmentAttr)
        {
            if (!string.IsNullOrEmpty(fragmentAttr) &&
                Enum.TryParse(fragmentAttr, true, out UniProtSequenceAttributes.FragmentType fragment))
            {
                return fragment;
            }
            return UniProtSequenceAttributes.FragmentType.unspecified;
        }

        private static int ComputeSequenceMass(string sequence)
        {
            if (string.IsNullOrEmpty(sequence))
                return 0;
            return (int)Math.Round(new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>()).MonoisotopicMass);
        }

        public Protein ParseEndElement(XmlReader xml, IEnumerable<string> modTypesToExclude, Dictionary<string, Modification> unknownModifications,
            bool isContaminant, string proteinDbLocation, string decoyIdentifier = "DECOY")
        {
            Protein protein = null;
            if (xml.Name == "feature")
            {
                ParseFeatureEndElement(xml, modTypesToExclude, unknownModifications);
            }
            if (xml.Name == "subfeature")
            {
                ParseSubFeatureEndElement(xml, modTypesToExclude, unknownModifications);
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
                protein = ParseEntryEndElement(xml, isContaminant, proteinDbLocation, modTypesToExclude, unknownModifications, decoyIdentifier);
            }
            return protein;
        }

        internal RNA ParseRnaEndElement(XmlReader xml, IEnumerable<string> modTypesToExclude,
            Dictionary<string, Modification> unknownModifications,
            bool isContaminant, string rnaDbLocation,string decoyIdentifier = "DECOY")
        {
            RNA result = null;
            if (xml.Name == "feature")
            {
                ParseFeatureEndElement(xml, modTypesToExclude, unknownModifications);
            }
            if (xml.Name == "subfeature")
            {
                ParseSubFeatureEndElement(xml, modTypesToExclude, unknownModifications);
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
                result = ParseRnaEntryEndElement(xml, isContaminant, rnaDbLocation, modTypesToExclude, unknownModifications, decoyIdentifier);
            }
            return result;
        }

        public Protein ParseEntryEndElement(XmlReader xml, bool isContaminant, string proteinDbLocation,
            IEnumerable<string> modTypesToExclude, Dictionary<string, Modification> unknownModifications, string decoyIdentifier = "DECOY")
        {
            Protein result = null;
            bool isDecoy = false;
            if (Accession != null && Sequence != null)
            {
                Sequence = ProteinDbLoader.SanitizeAminoAcidSequence(Sequence, 'X');

                // NEW: prune any sequence variants whose coordinates exceed the now-known sequence length
                PruneOutOfRangeSequenceVariants();

                ParseAnnotatedMods(OneBasedModifications, modTypesToExclude, unknownModifications, AnnotatedMods);
                if (Accession.StartsWith(decoyIdentifier))
                {
                    isDecoy = true;
                }
                result = new Protein(Sequence, Accession, Organism, GeneNames, OneBasedModifications, ProteolysisProducts, Name, FullName,
                    isDecoy, isContaminant, DatabaseReferences, SequenceVariations, null, null, DisulfideBonds, SpliceSites, proteinDbLocation,
                    false, DatasetEntryTag, DatabaseCreatedEntryTag, DatabaseModifiedEntryTag, DatabaseVersionEntryTag, XmlnsEntryTag, SequenceAttributes);
            }
            Clear();
            return result;
        }

        internal RNA ParseRnaEntryEndElement(XmlReader xml, bool isContaminant, string rnaDbLocation,
            IEnumerable<string> modTypesToExclude, Dictionary<string, Modification> unknownModifications, string decoyIdentifier = "DECOY")
        {
            RNA result = null;
            bool isDecoy = false;
            if (Accession != null && Sequence != null)
            {
                Sequence = ProteinDbLoader.SanitizeAminoAcidSequence(Sequence, 'X');
                if (Accession.StartsWith(decoyIdentifier))
                {
                    isDecoy = true;
                }

                // Prune for RNA as well (shared logic)
                PruneOutOfRangeSequenceVariants();

                ParseAnnotatedMods(OneBasedModifications, modTypesToExclude, unknownModifications, AnnotatedMods);
                if (Accession.StartsWith(decoyIdentifier))
                {
                    isDecoy = true;
                }
                result = new RNA(Sequence, Accession, OneBasedModifications, null, null, Name, Organism, rnaDbLocation,
                    isContaminant, isDecoy, GeneNames, [], ProteolysisProducts, SequenceVariations, null, null, FullName);
            }
            Clear();
            return result;
        }

        public void ParseSubFeatureEndElement(XmlReader xml, IEnumerable<string> modTypesToExclude, Dictionary<string, Modification> unknownModifications)
        {
            if (SubFeatureType == "modified residue")
            {
                SubFeatureDescription = SubFeatureDescription.Split(';')[0];
                AnnotatedVariantMods.Add((OneBasedFeatureSubPosition, SubFeatureDescription));
            }
        }

        public void ParseFeatureEndElement(XmlReader xml, IEnumerable<string> modTypesToExclude, Dictionary<string, Modification> unknownModifications)
        {
            if (FeatureType == "modified residue")
            {
                FeatureDescription = FeatureDescription.Split(';')[0];
                AnnotatedMods.Add((OneBasedFeaturePosition, FeatureDescription));
            }
            else if (FeatureType == "lipid moiety-binding region")
            {
                FeatureDescription = FeatureDescription.Split(';')[0];
                AnnotatedMods.Add((OneBasedFeaturePosition, FeatureDescription));
            }
            else if (FeatureType == "peptide" || FeatureType == "propeptide" || FeatureType == "chain" || FeatureType == "signal peptide")
            {
                string type = FeatureType;
                if (OneBasedBeginPosition.HasValue)
                {
                    type = type + "(" + (int)OneBasedBeginPosition.Value;
                    if (OneBasedEndPosition.HasValue)
                    {
                        type = type + "-" + (int)OneBasedEndPosition.Value + ")";
                    }
                    else
                    {
                        type += "-null)";
                    }
                }
                else
                {
                    if (OneBasedEndPosition.HasValue)
                    {
                        type = type + "(null-" + (int)OneBasedEndPosition.Value + ")";
                    }
                    else
                    {
                        type += "null-null";
                    }
                }
                ProteolysisProducts.Add(new TruncationProduct(OneBasedBeginPosition, OneBasedEndPosition, type));
            }
            else if (FeatureType == "sequence variant" && VariationValue != null && VariationValue != "")
            {
                bool appliesToThisSequence = true;
                if (!string.IsNullOrEmpty(LocationSequenceId))
                {
                    string acc = Accession ?? "";
                    appliesToThisSequence =
                        LocationSequenceId.Equals(acc, StringComparison.OrdinalIgnoreCase)
                        || (!string.IsNullOrEmpty(acc) && LocationSequenceId.Equals($"{acc}-1", StringComparison.OrdinalIgnoreCase));
                }

                if (appliesToThisSequence)
                {
                    ParseAnnotatedMods(OneBasedVariantModifications, modTypesToExclude, unknownModifications, AnnotatedVariantMods);

                    // NOTE: We can NOT validate coordinate vs sequence length here because sequence is usually parsed later.
                    // Validation is deferred to PruneOutOfRangeSequenceVariants() during ParseEntryEndElement.

                    if (OneBasedBeginPosition != null && OneBasedEndPosition != null)
                    {
                        SequenceVariations.Add(
                            new SequenceVariation(
                                (int)OneBasedBeginPosition,
                                (int)OneBasedEndPosition,
                                OriginalValue,
                                VariationValue,
                                FeatureDescription,
                                variantCallFormatDataString: null,
                                oneBasedModifications: OneBasedVariantModifications));
                    }
                    else if (OneBasedFeaturePosition >= 1)
                    {
                        SequenceVariations.Add(
                            new SequenceVariation(
                                OneBasedFeaturePosition,
                                OriginalValue,
                                VariationValue,
                                FeatureDescription,
                                variantCallFormatDataString: null,
                                oneBasedModifications: OneBasedVariantModifications));
                    }

                    AnnotatedVariantMods = new List<(int, string)>();
                    OneBasedVariantModifications = new Dictionary<int, List<Modification>>();
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
            else if (FeatureType == "splice site")
            {
                if (OneBasedBeginPosition != null && OneBasedEndPosition != null)
                {
                    SpliceSites.Add(new SpliceSite((int)OneBasedBeginPosition, (int)OneBasedEndPosition, FeatureDescription));
                }
                else if (OneBasedFeaturePosition >= 1)
                {
                    SpliceSites.Add(new SpliceSite(OneBasedFeaturePosition, FeatureDescription));
                }
            }

            OneBasedBeginPosition = null;
            OneBasedEndPosition = null;
            OneBasedFeaturePosition = -1;
            OriginalValue = "";
            VariationValue = "";
            LocationSequenceId = null;
        }

        private void ParseDatabaseReferenceEndElement(XmlReader xml)
        {
            DatabaseReferences.Add(
                new DatabaseReference(DBReferenceType, DBReferenceId,
                    Enumerable.Range(0, PropertyTypes.Count).Select(i => new Tuple<string, string>(PropertyTypes[i], PropertyValues[i])).ToList()));
            PropertyTypes = new List<string>();
            PropertyValues = new List<string>();
            DBReferenceType = null;
            DBReferenceId = null;
        }

        private void Clear()
        {
            DatasetEntryTag = null;
            DatabaseCreatedEntryTag = null;
            DatabaseModifiedEntryTag = null;
            DatabaseVersionEntryTag = null;
            XmlnsEntryTag = null;
            Accession = null;
            Name = null;
            FullName = null;
            Sequence = null;
            SequenceAttributes = null;
            Organism = null;
            FeatureType = null;
            FeatureDescription = null;
            SubFeatureType = null;
            SubFeatureDescription = null;
            OriginalValue = "";
            VariationValue = "";
            DBReferenceType = null;
            DBReferenceId = null;
            PropertyTypes = new List<string>();
            PropertyValues = new List<string>();
            OneBasedFeaturePosition = -1;
            OneBasedFeatureSubPosition = -1;
            AnnotatedMods = new List<(int, string)>();
            OneBasedModifications = new Dictionary<int, List<Modification>>();
            ProteolysisProducts = new List<TruncationProduct>();
            SequenceVariations = new List<SequenceVariation>();
            DisulfideBonds = new List<DisulfideBond>();
            SpliceSites = new List<SpliceSite>();
            DatabaseReferences = new List<DatabaseReference>();
            GeneNames = new List<Tuple<string, string>>();
            ReadingGene = false;
            ReadingOrganism = false;
            LocationSequenceId = null;
            AnnotatedVariantMods = new List<(int, string)>();
            OneBasedVariantModifications = new Dictionary<int, List<Modification>>();
        }

        private void PruneOutOfRangeSequenceVariants()
        {
            if (string.IsNullOrEmpty(Sequence) || SequenceVariations.Count == 0)
                return;

            int len = Sequence.Length;
            int removed = SequenceVariations.RemoveAll(v =>
                v.OneBasedBeginPosition > len || v.OneBasedEndPosition > len);

            if (removed > 0)
            {
                Trace.TraceWarning($"Pruned {removed} out-of-range sequence variant(s) for accession {Accession} (protein length {len}).");
            }
        }

        private static void ParseAnnotatedMods(
            Dictionary<int, List<Modification>> destination,
            IEnumerable<string> modTypesToExclude,
            Dictionary<string, Modification> unknownModifications,
            List<(int position, string originalModificationID)> annotatedMods)
        {
            foreach (var (annotatedModLocation, annotatedId) in annotatedMods)
            {
                if (ProteinDbLoader.IdWithMotifToMod.TryGetValue(annotatedId, out Modification foundMod)
                    || RnaDbLoader.IdWithMotifToMod.TryGetValue(annotatedId, out foundMod))
                {
                    if (!modTypesToExclude.Contains(foundMod.ModificationType))
                    {
                        if (destination.TryGetValue(annotatedModLocation, out var listOfModsAtThisLocation))
                        {
                            listOfModsAtThisLocation.Add(foundMod);
                        }
                        else
                        {
                            destination.Add(annotatedModLocation, new List<Modification> { foundMod });
                        }
                    }
                }
                else if (ProteinDbLoader.IdToPossibleMods.TryGetValue(annotatedId, out IList<Modification> mods)
                         || RnaDbLoader.IdToPossibleMods.TryGetValue(annotatedId, out mods))
                {
                    foreach (Modification mod in mods)
                    {
                        if (!modTypesToExclude.Contains(mod.ModificationType))
                        {
                            if (destination.TryGetValue(annotatedModLocation, out var listOfModsAtThisLocation))
                            {
                                listOfModsAtThisLocation.Add(mod);
                            }
                            else
                            {
                                destination.Add(annotatedModLocation, new List<Modification> { mod });
                            }
                            break;
                        }
                    }
                }
                else
                {
                    if (!unknownModifications.ContainsKey(annotatedId))
                    {
                        unknownModifications.Add(annotatedId, new Modification(annotatedId));
                    }
                }
            }
        }
    }
}