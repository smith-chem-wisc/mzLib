using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text.RegularExpressions;
using System.Xml;
using Transcriptomics;
using UsefulProteomicsDatabases.Transcriptomics;

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
        public string OriginalValue { get; private set; } = ""; // if no content is found, assume it is empty, not null (e.g. <original>A</original><variation/> for a deletion event)
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
        public UniProtSequenceAttributes SequenceAttributes { get; set; } = null; // this is used to store the sequence attributes from the <sequence> element, if present
        private List<(int, string)> AnnotatedMods = new List<(int position, string originalModificationID)>();
        private List<(int, string)> AnnotatedVariantMods = new List<(int position, string originalModificationID)>();
        // Captured isoform/sequence identifier from <location sequence="...">
        private string LocationSequenceId;

        /// <summary>
        /// Finalizes the parsing of a protein XML entry and constructs a <see cref="Protein"/> object.
        /// This method is called when the end of an &lt;entry&gt; element is reached during XML parsing.
        /// It sanitizes the sequence, prunes out-of-range sequence variants, resolves and attaches modifications,
        /// and aggregates all parsed data (such as gene names, proteolysis products, sequence variations, disulfide bonds, and splice sites)
        /// into a new <see cref="Protein"/> instance.
        /// After construction, the internal state is cleared to prepare for the next entry.
        /// </summary>
        /// <param name="xml">The <see cref="XmlReader"/> positioned at the end of the &lt;entry&gt; element.</param>
        /// <param name="isContaminant">Indicates whether the protein is a contaminant.</param>
        /// <param name="proteinDbLocation">The file path or identifier of the protein database source.</param>
        /// <param name="modTypesToExclude">A collection of modification types to exclude from the protein.</param>
        /// <param name="unknownModifications">A dictionary to collect modifications that could not be resolved.</param>
        /// <param name="decoyIdentifier">A string used to identify decoy proteins (default: "DECOY").</param>
        /// <returns>
        /// A constructed <see cref="Protein"/> object containing all parsed and resolved information,
        /// or <c>null</c> if the entry is incomplete.
        /// </returns>
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

        /// <summary>
        /// Parses and stores key metadata attributes from the current &lt;entry&gt; element in the XML.
        /// This includes dataset, creation date, modification date, version, and XML namespace information.
        /// The extracted values are assigned to the corresponding properties of the <see cref="ProteinXmlEntry"/> instance.
        /// This method is typically called when the parser encounters the start of a protein entry in a UniProt or similar XML file.
        /// </summary>
        /// <param name="xml">The <see cref="XmlReader"/> positioned at the &lt;entry&gt; element whose attributes are to be read.</param>
        private void ParseEntryAttributes(XmlReader xml)
        {
            DatasetEntryTag = xml.GetAttribute("dataset");
            DatabaseCreatedEntryTag = xml.GetAttribute("created");
            DatabaseModifiedEntryTag = xml.GetAttribute("modified");
            DatabaseVersionEntryTag = xml.GetAttribute("version");
            XmlnsEntryTag = xml.GetAttribute("xmlns");
        }
        /// <summary>
        /// Parses some attributes of a &lt;sequence&gt; XML element and assigns their values to the corresponding properties of the ProteinXmlEntry.
        /// Note: the Length and Mass of the sequence are computed based on the sequence string after parsing it.
        /// 
        /// Attribute definitions:
        /// - length: (string) The length of the protein sequence.
        /// - mass: (string) The mass of the protein sequence.
        /// - checksum: (string) The checksum value for the sequence.
        /// - modified: (string) The date the sequence was last modified.
        /// - version: (string) The version of the sequence.
        /// - precursor: (string) Indicates if the sequence is a precursor.
        /// - fragment: (FragmentType) Indicates the type of fragment (unspecified, single, multiple).
        /// </summary>
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

            // Read sequence and compute length/mass
            Sequence = SubstituteWhitespace.Replace(xml.ReadElementString(), "");
            int length = Sequence.Length;
            int mass = ComputeSequenceMass(Sequence);

            SequenceAttributes = new UniProtSequenceAttributes(length, mass, checksum, entryModified, sequenceVersion, isPrecursor, fragment);

        }
        // Helper method to parse the modified date attribute, with fallback to DateTime.Now if parsing fails.
        /// <summary>
        /// Parses the modified date attribute from the sequence element.
        /// Returns DateTime.Now if parsing fails or the attribute is missing.
        /// </summary>
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
                    // Parsing failed; falling back to current date.
                    System.Diagnostics.Trace.TraceWarning($"Warning: Failed to parse modified date '{modifiedAttr}'. Using DateTime.Now.");
                }
            }
            return DateTime.Now;
        }

        // Helper method to parse the sequence version attribute.
        /// <summary>
        /// Parses the version attribute from the sequence element.
        /// Returns -1 if parsing fails or the attribute is missing.
        /// </summary>
        private static int ParseSequenceVersion(string versionAttr)
        {
            if (int.TryParse(versionAttr, out int version))
            {
                return version;
            }
            return -1;
        }

        // Helper method to parse the precursor attribute.
        /// <summary>
        /// Parses the precursor attribute from the sequence element.
        /// Returns false if the attribute is missing or not "true".
        /// </summary>
        private static bool ParseIsPrecursor(string precursorAttr)
        {
            return !string.IsNullOrEmpty(precursorAttr) && precursorAttr.Equals("true", StringComparison.OrdinalIgnoreCase);
        }

        // Helper method to parse the fragment type attribute.
        /// <summary>
        /// Parses the fragment attribute from the sequence element.
        /// Returns FragmentType.unspecified if parsing fails or the attribute is missing.
        /// </summary>
        private static UniProtSequenceAttributes.FragmentType ParseFragmentType(string fragmentAttr)
        {
            if (!string.IsNullOrEmpty(fragmentAttr) &&
                Enum.TryParse(fragmentAttr, true, out UniProtSequenceAttributes.FragmentType fragment))
            {
                return fragment;
            }
            return UniProtSequenceAttributes.FragmentType.unspecified;
        }

        /// <summary>
        /// Computes the monoisotopic mass of a protein or nucleic acid sequence without modifications.
        /// If the input sequence is null or empty, returns 0.
        /// Internally, constructs a <see cref="PeptideWithSetModifications"/> using the provided sequence and an empty modification dictionary,
        /// then returns the rounded monoisotopic mass as an integer.
        /// This method is used to populate sequence attributes such as mass during XML parsing.
        /// </summary>
        private static int ComputeSequenceMass(string sequence)
        {
            if (string.IsNullOrEmpty(sequence))
                return 0;
            return (int)Math.Round(new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>()).MonoisotopicMass);
        }
        /// <summary>
        /// Handles the end of an XML element during protein database parsing, updating the internal state or finalizing objects as needed.
        /// Depending on the element name, this method processes and stores feature, subfeature, database reference, gene, and organism information,
        /// or, if the end of an &lt;entry&gt; element is reached, constructs and returns a fully populated <see cref="Protein"/> object.
        /// For &lt;feature&gt; and &lt;subfeature&gt; elements, it attaches modifications or proteolytic products.
        /// For &lt;dbReference&gt;, it records database cross-references.
        /// For &lt;gene&gt; and &lt;organism&gt;, it updates parsing state flags.
        /// For &lt;entry&gt;, it aggregates all parsed data, resolves modifications, and returns a new <see cref="Protein"/> instance,
        /// clearing the internal state for the next entry.
        /// </summary>
        /// <param name="xml">The <see cref="XmlReader"/> positioned at the end of the current XML element.</param>
        /// <param name="modTypesToExclude">A collection of modification types to exclude from the protein.</param>
        /// <param name="unknownModifications">A dictionary to collect modifications that could not be resolved.</param>
        /// <param name="isContaminant">Indicates whether the protein is a contaminant.</param>
        /// <param name="proteinDbLocation">The file path or identifier of the protein database source.</param>
        /// <param name="decoyIdentifier">A string used to identify decoy proteins (default: "DECOY").</param>
        /// <returns>
        /// A constructed <see cref="Protein"/> object if the end of an &lt;entry&gt; element is reached and all required data is present;
        /// otherwise, <c>null</c>.
        /// </returns>
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
        /// <summary>
        /// Handles the end of an XML element during RNA database parsing, updating the internal state or finalizing objects as needed.
        /// Depending on the element name, this method processes and stores feature, subfeature, and database reference information,
        /// or, if the end of an &lt;entry&gt; element is reached, constructs and returns a fully populated <see cref="RNA"/> object.
        /// For &lt;feature&gt; and &lt;subfeature&gt; elements, it attaches modifications or truncation products.
        /// For &lt;dbReference&gt;, it records database cross-references.
        /// For &lt;gene&gt; and &lt;organism&gt;, it updates parsing state flags.
        /// For &lt;entry&gt;, it aggregates all parsed data, resolves modifications, and returns a new <see cref="RNA"/> instance,
        /// clearing the internal state for the next entry.
        /// </summary>
        /// <param name="xml">The <see cref="XmlReader"/> positioned at the end of the current XML element.</param>
        /// <param name="modTypesToExclude">A collection of modification types to exclude from the RNA.</param>
        /// <param name="unknownModifications">A dictionary to collect modifications that could not be resolved.</param>
        /// <param name="isContaminant">Indicates whether the RNA is a contaminant.</param>
        /// <param name="rnaDbLocation">The file path or identifier of the RNA database source.</param>
        /// <param name="decoyIdentifier">A string used to identify decoy RNAs (default: "DECOY").</param>
        /// <returns>
        /// A constructed <see cref="RNA"/> object if the end of an &lt;entry&gt; element is reached and all required data is present;
        /// otherwise, <c>null</c>.
        /// </returns>
        internal RNA ParseRnaEndElement(XmlReader xml, IEnumerable<string> modTypesToExclude,
            Dictionary<string, Modification> unknownModifications,
            bool isContaminant, string rnaDbLocation, string decoyIdentifier = "DECOY")
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

        /// <summary>
        /// Finalizes the parsing of a protein XML entry and constructs a <see cref="Protein"/> object from the accumulated data.
        /// This method is called when the end of an &lt;entry&gt; element is reached during XML parsing.
        /// It performs several key tasks:
        /// <list type="bullet">
        ///   <item>Sanitizes the parsed sequence (e.g., replacing invalid amino acids with 'X').</item>
        ///   <item>Prunes any sequence variants whose coordinates exceed the sequence length.</item>
        ///   <item>Resolves and attaches all annotated modifications, excluding those of specified types or unknowns.</item>
        ///   <item>Determines if the protein is a decoy based on the accession and decoy identifier.</item>
        ///   <item>Aggregates all parsed data (gene names, proteolysis products, sequence variations, disulfide bonds, splice sites, database references, and sequence attributes) into a new <see cref="Protein"/> instance.</item>
        ///   <item>Clears the internal state of the <see cref="ProteinXmlEntry"/> to prepare for parsing the next entry.</item>
        /// </list>
        /// If either the accession or sequence is missing, returns <c>null</c>.
        /// </summary>
        /// <param name="xml">The <see cref="XmlReader"/> positioned at the end of the &lt;entry&gt; element.</param>
        /// <param name="isContaminant">Indicates whether the protein is a contaminant.</param>
        /// <param name="proteinDbLocation">The file path or identifier of the protein database source.</param>
        /// <param name="modTypesToExclude">A collection of modification types to exclude from the protein.</param>
        /// <param name="unknownModifications">A dictionary to collect modifications that could not be resolved.</param>
        /// <param name="decoyIdentifier">A string used to identify decoy proteins (default: "DECOY").</param>
        /// <returns>
        /// A constructed <see cref="Protein"/> object containing all parsed and resolved information,
        /// or <c>null</c> if the entry is incomplete.
        /// </returns>

        public Protein ParseEntryEndElement(XmlReader xml, bool isContaminant, string proteinDbLocation, IEnumerable<string> modTypesToExclude, Dictionary<string, Modification> unknownModifications, string decoyIdentifier = "DECOY")
        {
            Protein result = null;
            bool isDecoy = false;
            if (Accession != null && Sequence != null)
            {
                // sanitize the sequence to replace unexpected characters with X (unknown amino acid)
                // sometimes strange characters get added by RNA sequencing software, etc.
                Sequence = ProteinDbLoader.SanitizeAminoAcidSequence(Sequence, 'X');

                ParseAnnotatedMods(OneBasedModifications, modTypesToExclude, unknownModifications, AnnotatedMods);
                //prune any sequence variants whose coordinates exceed the known sequence length
                PruneOutOfRangeSequenceVariants();
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
        /// <summary>
        /// Finalizes the parsing of an RNA XML entry and constructs an <see cref="RNA"/> object from the accumulated data.
        /// This method is called when the end of an &lt;entry&gt; element is reached during XML parsing for RNA records.
        /// It performs several key tasks:
        /// <list type="bullet">
        ///   <item>Sanitizes the parsed sequence (e.g., replacing invalid characters with 'X').</item>
        ///   <item>Prunes any sequence variants whose coordinates exceed the sequence length.</item>
        ///   <item>Resolves and attaches all annotated modifications, excluding those of specified types or unknowns.</item>
        ///   <item>Determines if the RNA is a decoy based on the accession and decoy identifier.</item>
        ///   <item>Aggregates all parsed data (gene names, proteolysis products, sequence variations, and other metadata) into a new <see cref="RNA"/> instance.</item>
        ///   <item>Clears the internal state of the <see cref="ProteinXmlEntry"/> to prepare for parsing the next entry.</item>
        /// </list>
        /// If either the accession or sequence is missing, returns <c>null</c>.
        /// </summary>
        /// <param name="xml">The <see cref="XmlReader"/> positioned at the end of the &lt;entry&gt; element.</param>
        /// <param name="isContaminant">Indicates whether the RNA is a contaminant.</param>
        /// <param name="rnaDbLocation">The file path or identifier of the RNA database source.</param>
        /// <param name="modTypesToExclude">A collection of modification types to exclude from the RNA.</param>
        /// <param name="unknownModifications">A dictionary to collect modifications that could not be resolved.</param>
        /// <param name="decoyIdentifier">A string used to identify decoy RNAs (default: "DECOY").</param>
        /// <returns>
        /// A constructed <see cref="RNA"/> object containing all parsed and resolved information,
        /// or <c>null</c> if the entry is incomplete.
        /// </returns>
        internal RNA ParseRnaEntryEndElement(XmlReader xml, bool isContaminant, string rnaDbLocation,
            IEnumerable<string> modTypesToExclude, Dictionary<string, Modification> unknownModifications, string decoyIdentifier = "DECOY")
        {
            RNA result = null;
            bool isDecoy = false;
            if (Accession != null && Sequence != null)
            {
                // sanitize the sequence to replace unexpected characters with X (unknown amino acid)
                // sometimes strange characters get added by RNA sequencing software, etc.
                Sequence = ProteinDbLoader.SanitizeAminoAcidSequence(Sequence, 'X');
                //prune any sequence variants whose coordinates exceed the known sequence length
                PruneOutOfRangeSequenceVariants();
                if (Accession.StartsWith(decoyIdentifier))
                {
                    isDecoy = true;
                }

                ParseAnnotatedMods(OneBasedModifications, modTypesToExclude, unknownModifications, AnnotatedMods);
                result = new RNA(Sequence, Accession, OneBasedModifications, null, null, Name, Organism, rnaDbLocation,
                    isContaminant, isDecoy, GeneNames, [], ProteolysisProducts, SequenceVariations, null, null, FullName);
            }
            Clear();
            return result;
        }

        /// <summary>
        /// Finish parsing a subfeature element
        /// </summary>
        public void ParseSubFeatureEndElement(XmlReader xml, IEnumerable<string> modTypesToExclude, Dictionary<string, Modification> unknownModifications)
        {
            if (SubFeatureType == "modified residue")
            {
                SubFeatureDescription = SubFeatureDescription.Split(';')[0];
                AnnotatedVariantMods.Add((OneBasedFeatureSubPosition, SubFeatureDescription));
            }
        }

        /// <summary>
        /// Processes the end of a &lt;feature&gt; element during XML parsing and updates the internal state with the parsed feature information.
        /// Depending on the feature type, this method:
        /// <list type="bullet">
        ///   <item>Adds modification annotations for "modified residue" and "lipid moiety-binding region" features.</item>
        ///   <item>Creates and adds <see cref="TruncationProduct"/> objects for proteolytic features such as "peptide", "propeptide", "chain", and "signal peptide".</item>
        ///   <item>Handles "sequence variant" features by creating <see cref="SequenceVariation"/> objects, including variant-specific modifications, and ensures they apply to the correct sequence or isoform.</item>
        ///   <item>Creates and adds <see cref="DisulfideBond"/> or <see cref="SpliceSite"/> objects for their respective feature types, using available position information.</item>
        /// </list>
        /// After processing, resets feature-related state variables to prepare for the next feature.
        /// </summary>
        /// <param name="xml">The <see cref="XmlReader"/> positioned at the end of the &lt;feature&gt; element.</param>
        /// <param name="modTypesToExclude">A collection of modification types to exclude from the protein.</param>
        /// <param name="unknownModifications">A dictionary to collect modifications that could not be resolved.</param>
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
                //next we are going to add test descrbing the begin and end positions (if any) of the feature. This results in increased information in the output about feature location in the protein
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
                                //variantCallFormatDataString: null,
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
                                //variantCallFormatDataString: null,
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

        /// <summary>
        /// Finish parsing a database reference element
        /// </summary>
        /// <param name="xml"></param>
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

        /// <summary>
        /// Clear this object's properties
        /// </summary>
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
        /// <summary>
        /// Resolves and attaches annotated modifications to the specified destination dictionary based on parsed feature or variant annotations.
        /// For each annotated modification, attempts to look up the modification by its identifier (with motif) in both protein and RNA modification dictionaries.
        /// If found and not excluded by <paramref name="modTypesToExclude"/>, the modification is added to the destination at the specified position.
        /// If not found by identifier, attempts to resolve the modification by possible matches (without motif) and adds the first non-excluded match.
        /// If no match is found, records the modification as unknown in <paramref name="unknownModifications"/> to avoid repeated warnings.
        /// This method is used to populate the protein or variant modification dictionaries during XML parsing.
        /// </summary>
        /// <param name="destination">Dictionary mapping one-based positions to lists of modifications to be populated.</param>
        /// <param name="modTypesToExclude">A collection of modification types to exclude from assignment.</param>
        /// <param name="unknownModifications">A dictionary to collect modifications that could not be resolved by identifier or type.</param>
        /// <param name="annotatedMods">List of (position, modification identifier) tuples parsed from XML features or subfeatures.</param>
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
    }
}