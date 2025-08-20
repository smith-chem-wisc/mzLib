using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using MzLibUtil;
using Omics.BioPolymer;
using System.Data;
using Chemistry;
using Easy.Common.Extensions;


namespace Proteomics
{
    public class Protein : IBioPolymer, IEquatable<Protein>, IComparable<Protein>
    {
        private List<TruncationProduct> _proteolysisProducts;

        /// <summary>
        /// Protein. Filters out modifications that do not match their amino acid target site.
        /// </summary>
        /// <param name="sequence">Base sequence of the protein.</param>
        /// <param name="accession">Unique accession for the protein.</param>
        /// <param name="organism">Organism with this protein.</param>
        /// <param name="geneNames">List of gene names as tuple of (nameType, name), e.g. (primary, HLA-A)</param>
        /// <param name="oneBasedModifications">Modifications at positions along the sequence.</param>
        /// <param name="proteolysisProducts"></param>
        /// <param name="name"></param>
        /// <param name="fullName"></param>
        /// <param name="isDecoy"></param>
        /// <param name="isContaminant"></param>
        /// <param name="databaseReferences"></param>
        /// <param name="sequenceVariations"></param>
        /// <param name="disulfideBonds"></param>
        /// <param name="spliceSites"></param>
        /// <param name="databaseFilePath"></param>
        public Protein(string sequence, string accession, string organism = null, List<Tuple<string, string>> geneNames = null,
            IDictionary<int, List<Modification>> oneBasedModifications = null, List<TruncationProduct> proteolysisProducts = null,
            string name = null, string fullName = null, bool isDecoy = false, bool isContaminant = false, List<DatabaseReference> databaseReferences = null,
            List<SequenceVariation> sequenceVariations = null, List<SequenceVariation> appliedSequenceVariations = null, string sampleNameForVariants = null,
            List<DisulfideBond> disulfideBonds = null, List<SpliceSite> spliceSites = null, string databaseFilePath = null, bool addTruncations = false, 
            string dataset = "unknown", string created = "unknown", string modified = "unknown", string version = "unknown", string xmlns = "http://uniprot.org/uniprot",
            UniProtSequenceAttributes uniProtSequenceAttributes = null)
        {
            BaseSequence = sequence;
            NonVariantProtein = this;
            Accession = accession;

            Name = name;
            Organism = organism;
            FullName = fullName;
            IsDecoy = isDecoy;
            IsContaminant = isContaminant;
            DatabaseFilePath = databaseFilePath;
            SampleNameForVariants = sampleNameForVariants;

            GeneNames = geneNames ?? new List<Tuple<string, string>>();
            _proteolysisProducts = proteolysisProducts ?? new List<TruncationProduct>();
            SequenceVariations = sequenceVariations ?? new List<SequenceVariation>();
            AppliedSequenceVariations = appliedSequenceVariations ?? new List<SequenceVariation>();
            OriginalNonVariantModifications = oneBasedModifications ?? new Dictionary<int, List<Modification>>();
            if (oneBasedModifications != null)
            {
                OneBasedPossibleLocalizedModifications = ((IBioPolymer)this).SelectValidOneBaseMods(oneBasedModifications);
            }
            else
            {
                OneBasedPossibleLocalizedModifications = new Dictionary<int, List<Modification>>();
            }
            DatabaseReferences = databaseReferences ?? new List<DatabaseReference>();
            DisulfideBonds = disulfideBonds ?? new List<DisulfideBond>();
            SpliceSites = spliceSites ?? new List<SpliceSite>();

            if (addTruncations)
            {
                this.AddTruncations();
            }
            DatasetEntryTag = dataset;
            CreatedEntryTag = created;
            ModifiedEntryTag = modified;
            VersionEntryTag = version;
            XmlnsEntryTag = xmlns;
            UniProtSequenceAttributes = uniProtSequenceAttributes ?? new UniProtSequenceAttributes(Length, (int)Math.Round(new PeptideWithSetModifications(BaseSequence, new Dictionary<string,Modification>()).MonoisotopicMass), "unknown", DateTime.Now, -1);
        }

        /// <summary>
        /// Protein construction that clones a protein but assigns a different base sequence
        /// For use in SILAC experiments and in decoy construction
        /// </summary>
        /// <param name="originalProtein"></param>
        /// <param name="newBaseSequence"></param>
        /// <param name="silacAccession"></param>
        public Protein(Protein originalProtein, string newBaseSequence)
        {
            BaseSequence = newBaseSequence;
            Accession = originalProtein.Accession;
            NonVariantProtein = originalProtein.ConsensusVariant as Protein;
            Name = originalProtein.Name;
            Organism = originalProtein.Organism;
            FullName = originalProtein.FullName;
            IsDecoy = originalProtein.IsDecoy;
            IsContaminant = originalProtein.IsContaminant;
            DatabaseFilePath = originalProtein.DatabaseFilePath;
            SampleNameForVariants = originalProtein.SampleNameForVariants;
            GeneNames = originalProtein.GeneNames;
            _proteolysisProducts = originalProtein._proteolysisProducts;
            SequenceVariations = originalProtein.SequenceVariations;
            AppliedSequenceVariations = originalProtein.AppliedSequenceVariations;
            OriginalNonVariantModifications = originalProtein.OriginalNonVariantModifications;
            OneBasedPossibleLocalizedModifications = originalProtein.OneBasedPossibleLocalizedModifications;
            DatabaseReferences = originalProtein.DatabaseReferences;
            DisulfideBonds = originalProtein.DisulfideBonds;
            SpliceSites = originalProtein.SpliceSites;
            DatabaseFilePath = originalProtein.DatabaseFilePath;
            DatasetEntryTag = originalProtein.DatasetEntryTag;
            CreatedEntryTag = originalProtein.CreatedEntryTag;
            ModifiedEntryTag = originalProtein.ModifiedEntryTag;
            VersionEntryTag = originalProtein.VersionEntryTag;
            XmlnsEntryTag = originalProtein.XmlnsEntryTag;
            UniProtSequenceAttributes = originalProtein.UniProtSequenceAttributes;
        }

        /// <summary>
        /// Protein construction with applied variations
        /// </summary>
        /// <param name="variantBaseSequence"></param>
        /// <param name="protein"></param>
        /// <param name="appliedSequenceVariations"></param>
        /// <param name="applicableProteolysisProducts"></param>
        /// <param name="oneBasedModifications"></param>
        /// <param name="sampleNameForVariants"></param>
        public Protein(string variantBaseSequence, Protein protein, IEnumerable<SequenceVariation> appliedSequenceVariations,
            IEnumerable<TruncationProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
            : this(
                  variantBaseSequence,
                  VariantApplication.GetAccession(protein, appliedSequenceVariations),
                  organism: protein.Organism,
                  geneNames: new List<Tuple<string, string>>(protein.GeneNames),
                  oneBasedModifications: oneBasedModifications != null ? oneBasedModifications.ToDictionary(x => x.Key, x => x.Value) : new Dictionary<int, List<Modification>>(),
                  proteolysisProducts: new List<TruncationProduct>(applicableProteolysisProducts ?? new List<TruncationProduct>()),
                  name: VariantApplication.GetVariantName(protein.Name, appliedSequenceVariations),
                  fullName: VariantApplication.GetVariantName(protein.FullName, appliedSequenceVariations), 
                  isDecoy: protein.IsDecoy,
                  isContaminant: protein.IsContaminant,
                  databaseReferences: new List<DatabaseReference>(protein.DatabaseReferences),
                  sequenceVariations: new List<SequenceVariation>(protein.SequenceVariations),
                  disulfideBonds: new List<DisulfideBond>(protein.DisulfideBonds),
                  spliceSites: new List<SpliceSite>(protein.SpliceSites),
                  databaseFilePath: protein.DatabaseFilePath,
                  dataset: protein.DatasetEntryTag,
                  created: protein.CreatedEntryTag,
                  modified: protein.ModifiedEntryTag,
                  version: protein.VersionEntryTag,
                  xmlns: protein.XmlnsEntryTag,
                  uniProtSequenceAttributes: protein.UniProtSequenceAttributes)
        {
            NonVariantProtein = protein.ConsensusVariant as Protein;
            OriginalNonVariantModifications = ConsensusVariant.OriginalNonVariantModifications;
            AppliedSequenceVariations = (appliedSequenceVariations ?? new List<SequenceVariation>()).ToList();
            SampleNameForVariants = sampleNameForVariants;
        }

        /// <summary>
        /// Modifications (values) located at one-based protein positions (keys)
        /// </summary>
        public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; private set; }

        /// <summary>
        /// The list of gene names consists of tuples, where Item1 is the type of gene name, and Item2 is the name. There may be many genes and names of a certain type produced when reading an XML protein database.
        /// </summary>
        public List<Tuple<string, string>> GeneNames { get; }

        /// <summary>
        /// Unique accession for this protein.
        /// </summary>
        public string Accession { get; }

        /// <summary>
        /// Base sequence, which may contain applied sequence variations.
        /// </summary>
        public string BaseSequence { get; private set; }

        public string Organism { get; }
        public bool IsDecoy { get; }
        public int Length => BaseSequence.Length;
        public string FullDescription => Accession + "|" + Name + "|" + FullName;
        public string Name { get; }
        public string FullName { get; }
        public bool IsContaminant { get; }
        public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

        #region Database Handling and XML Parsed Fields

        /// <summary>
        /// Sequence Variants as defined in the parsed XML database
        /// </summary>
        public List<SequenceVariation> SequenceVariations { get; }

        /// <summary>
        /// Disulfide Bonds as defined in the parsed XML database
        /// </summary>
        public List<DisulfideBond> DisulfideBonds { get; }

        /// <summary>
        /// Splice Sites as defined in the parsed XML Database
        /// </summary>
        public List<SpliceSite> SpliceSites { get; }

        //TODO: Generate all the proteolytic products as distinct proteins during XML reading and delete the TruncationProducts parameter
        /// <summary>
        /// Truncation products as defined in the parsed XML Database
        /// </summary>
        public List<TruncationProduct> TruncationProducts => _proteolysisProducts;

        /// <summary>
        /// The references for a protein in the parsed XML Database
        /// </summary>
        public List<DatabaseReference> DatabaseReferences { get; }

        public string DatabaseFilePath { get; }
        public string DatasetEntryTag { get; private set; }
        public string CreatedEntryTag { get; private set; }
        public string ModifiedEntryTag { get; private set; }
        public string VersionEntryTag { get; private set; }
        public string XmlnsEntryTag { get; private set; }
        public UniProtSequenceAttributes UniProtSequenceAttributes { get; private set; }
        /// <summary>
        /// Formats a string for a UniProt fasta header. See https://www.uniprot.org/help/fasta-headers.
        /// Note that the db field isn't very applicable here, so mz is placed in to denote written by mzLib.
        /// </summary>
        public string GetUniProtFastaHeader()
        {
            var n = GeneNames.FirstOrDefault();
            string geneName = n == null ? "" : n.Item2;

            return string.Format("mz|{0}|{1} {2} OS={3} GN={4}", Accession, Name, FullName, Organism, geneName);
        }

        /// <summary>
        /// Formats a string for an ensembl header
        /// </summary>
        public string GetEnsemblFastaHeader()
        {
            return string.Format("{0} {1}", Accession, FullName);
        }

        #endregion

        /// <summary>
        /// Gets peptides for digestion of a protein
        /// Legacy
        /// </summary>
        public IEnumerable<PeptideWithSetModifications> Digest(DigestionParams digestionParams,
            List<Modification> allKnownFixedModifications, List<Modification> variableModifications,
            List<SilacLabel> silacLabels = null, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null,
            bool topDownTruncationSearch = false) =>
            Digest((IDigestionParams)digestionParams, allKnownFixedModifications, variableModifications, silacLabels, turnoverLabels, topDownTruncationSearch)
                .Cast<PeptideWithSetModifications>();

        /// <summary>
        /// Gets peptides for digestion of a protein
        /// Implemented with interfaces to allow for use of both Proteomics and Omics classes
        /// </summary>
        public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams, List<Modification> allKnownFixedModifications,
            List<Modification> variableModifications, List<SilacLabel> silacLabels = null, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null, bool topDownTruncationSearch = false)
        {

            if (digestionParams is not DigestionParams digestionParameters)
                throw new ArgumentException(
                    "DigestionParameters must be of type DigestionParams for protein digestion");


            //can't be null
            allKnownFixedModifications = allKnownFixedModifications ?? new List<Modification>();
            // add in any modifications that are caused by protease digestion
            if (digestionParameters.Protease.CleavageMod != null && !allKnownFixedModifications.Contains(digestionParameters.Protease.CleavageMod))
            {
                allKnownFixedModifications.Add(digestionParameters.Protease.CleavageMod);
            }
            variableModifications = variableModifications ?? new List<Modification>();
            CleavageSpecificity searchModeType = digestionParameters.SearchModeType;

            ProteinDigestion digestion = new(digestionParameters, allKnownFixedModifications, variableModifications);
            IEnumerable<ProteolyticPeptide> unmodifiedPeptides =
                searchModeType == CleavageSpecificity.Semi ?
                digestion.SpeedySemiSpecificDigestion(this) :
                    digestion.Digestion(this, topDownTruncationSearch);

            if (digestionParameters.KeepNGlycopeptide || digestionParameters.KeepOGlycopeptide)
            {
                unmodifiedPeptides = GetGlycoPeptides(unmodifiedPeptides, digestionParameters.KeepNGlycopeptide, digestionParameters.KeepOGlycopeptide);
            }

            IEnumerable<PeptideWithSetModifications> modifiedPeptides = unmodifiedPeptides.SelectMany(peptide => 
                peptide.GetModifiedPeptides(allKnownFixedModifications, digestionParameters, variableModifications));

            //Remove terminal modifications (if needed)
            if (searchModeType == CleavageSpecificity.SingleN ||
                searchModeType == CleavageSpecificity.SingleC ||
                (searchModeType == CleavageSpecificity.None && (digestionParams.FragmentationTerminus == FragmentationTerminus.N || digestionParams.FragmentationTerminus == FragmentationTerminus.C)))
            {
                modifiedPeptides = RemoveTerminalModifications(modifiedPeptides, digestionParams.FragmentationTerminus, allKnownFixedModifications);
            }

            //add silac labels (if needed)
            if (silacLabels != null)
            {
                return GetSilacPeptides(modifiedPeptides, silacLabels, digestionParameters.GeneratehUnlabeledProteinsForSilac, turnoverLabels);
            }

            return modifiedPeptides;
        }

        /// <summary>
        /// Remove terminal modifications from the C-terminus of SingleN peptides and the N-terminus of SingleC peptides/
        /// These terminal modifications create redundant entries and increase search time
        /// </summary>
        internal static IEnumerable<PeptideWithSetModifications> RemoveTerminalModifications(IEnumerable<PeptideWithSetModifications> modifiedPeptides, FragmentationTerminus fragmentationTerminus, IEnumerable<Modification> allFixedMods)
        {
            string terminalStringToLookFor = fragmentationTerminus == FragmentationTerminus.N ? "C-terminal" : "N-terminal";
            List<Modification> fixedTerminalMods = allFixedMods.Where(x => x.LocationRestriction.Contains(terminalStringToLookFor)).ToList();
            foreach (PeptideWithSetModifications pwsm in modifiedPeptides)
            {
                if (!pwsm.AllModsOneIsNterminus.Values.Any(x => x.LocationRestriction.Contains(terminalStringToLookFor) && !fixedTerminalMods.Contains(x)))
                {
                    yield return pwsm;
                }
            }
        }

        /// <summary>
        /// Add additional peptides with SILAC amino acids
        /// </summary>
        internal IEnumerable<PeptideWithSetModifications> GetSilacPeptides(IEnumerable<PeptideWithSetModifications> originalPeptides, List<SilacLabel> silacLabels, bool generateUnlabeledProteins, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels)
        {
            //if this is a multiplex experiment (pooling multiple samples, not a turnover), then only create the fully unlabeled/labeled peptides
            if (turnoverLabels == null)
            {
                //unlabeled peptides
                if (generateUnlabeledProteins)
                {
                    foreach (PeptideWithSetModifications pwsm in originalPeptides)
                    {
                        yield return pwsm;
                    }
                }

                //fully labeled peptides
                foreach (SilacLabel label in silacLabels)
                {
                    Protein silacProtein = GenerateFullyLabeledSilacProtein(label);
                    foreach (PeptideWithSetModifications pwsm in originalPeptides)
                    {
                        //duplicate the peptides with the updated protein sequence that contains only silac labels
                        yield return new PeptideWithSetModifications(silacProtein, pwsm.DigestionParams, pwsm.OneBasedStartResidueInProtein, pwsm.OneBasedEndResidueInProtein, pwsm.CleavageSpecificityForFdrCategory, pwsm.PeptideDescription, pwsm.MissedCleavages, pwsm.AllModsOneIsNterminus, pwsm.NumFixedMods);
                    }
                }
            }
            else //if this is a turnover experiment, we want to be able to look for peptides containing mixtures of heavy and light amino acids (typically occurs for missed cleavages)
            {
                (SilacLabel startLabel, SilacLabel endLabel) turnoverLabelsValue = turnoverLabels.Value;
                SilacLabel startLabel = turnoverLabelsValue.startLabel;
                SilacLabel endLabel = turnoverLabelsValue.endLabel;

                //This allows you to move from one label to another (rather than unlabeled->labeled or labeled->unlabeled). Useful for when your lab is swimming in cash and you have stock in a SILAC company
                if (startLabel != null && endLabel != null) //if neither the start nor end conditions are unlabeled, then generate fully labeled proteins using the "startLabel" (otherwise maintain the unlabeled)
                {
                    Protein silacStartProtein = GenerateFullyLabeledSilacProtein(startLabel);
                    PeptideWithSetModifications[] originalPeptideArray = originalPeptides.ToArray();
                    for (int i = 0; i < originalPeptideArray.Length; i++)
                    {
                        PeptideWithSetModifications pwsm = originalPeptideArray[i];
                        //duplicate the peptides with the updated protein sequence that contains only silac labels
                        originalPeptideArray[i] = new PeptideWithSetModifications(silacStartProtein, pwsm.DigestionParams, pwsm.OneBasedStartResidueInProtein, pwsm.OneBasedEndResidueInProtein, pwsm.CleavageSpecificityForFdrCategory, pwsm.PeptideDescription, pwsm.MissedCleavages, pwsm.AllModsOneIsNterminus, pwsm.NumFixedMods);
                    }
                    originalPeptides = originalPeptideArray;

                    //modify the end label amino acids to recognize the new "original" amino acid
                    //get the residues that were changed
                    List<SilacLabel> originalLabels = new List<SilacLabel> { startLabel };
                    if (startLabel.AdditionalLabels != null)
                    {
                        originalLabels.AddRange(startLabel.AdditionalLabels);
                    }
                    SilacLabel startLabelWithSharedOriginalAminoAcid = originalLabels.Where(x => x.OriginalAminoAcid == endLabel.OriginalAminoAcid).FirstOrDefault();
                    SilacLabel updatedEndLabel = startLabelWithSharedOriginalAminoAcid == null ?
                        endLabel :
                        new SilacLabel(startLabelWithSharedOriginalAminoAcid.AminoAcidLabel, endLabel.AminoAcidLabel, endLabel.LabelChemicalFormula, endLabel.ConvertMassDifferenceToDouble());
                    if (endLabel.AdditionalLabels != null)
                    {
                        foreach (SilacLabel additionalLabel in endLabel.AdditionalLabels)
                        {
                            startLabelWithSharedOriginalAminoAcid = originalLabels.Where(x => x.OriginalAminoAcid == additionalLabel.OriginalAminoAcid).FirstOrDefault();
                            updatedEndLabel.AddAdditionalSilacLabel(
                                startLabelWithSharedOriginalAminoAcid == null ?
                                additionalLabel :
                                new SilacLabel(startLabelWithSharedOriginalAminoAcid.AminoAcidLabel, additionalLabel.AminoAcidLabel, additionalLabel.LabelChemicalFormula, additionalLabel.ConvertMassDifferenceToDouble()));
                        }
                    }

                    //double check that all labeled amino acids can become unlabeled/relabeled
                    if (startLabel.AdditionalLabels != null)
                    {
                        foreach (SilacLabel originalLabel in originalLabels)
                        {
                            if (updatedEndLabel.OriginalAminoAcid != originalLabel.AminoAcidLabel &&
                                (updatedEndLabel.AdditionalLabels == null || !updatedEndLabel.AdditionalLabels.Any(x => x.OriginalAminoAcid == originalLabel.AminoAcidLabel)))
                            {
                                updatedEndLabel.AddAdditionalSilacLabel(new SilacLabel(originalLabel.AminoAcidLabel, originalLabel.OriginalAminoAcid, originalLabel.LabelChemicalFormula, originalLabel.ConvertMassDifferenceToDouble()));
                            }
                        }
                    }
                    endLabel = updatedEndLabel;
                }

                //add all unlabeled (or if no unlabeled, then the startLabeled) peptides
                foreach (PeptideWithSetModifications pwsm in originalPeptides)
                {
                    yield return pwsm;
                }

                //the order (below) matters when neither labels are null, because the fully labeled "start" has already been created above, so we want to use the end label here if it's not unlabeled (null)
                SilacLabel label = endLabel ?? startLabel; //pick the labeled (not the unlabeled). If no unlabeled, take the endLabel

                Protein silacEndProtein = GenerateFullyLabeledSilacProtein(label);

                //add all peptides containing any label (may also contain unlabeled)
                if (label.AdditionalLabels == null) //if there's only one (which is common)
                {
                    //get the residues to change
                    char originalResidue = label.OriginalAminoAcid;
                    char labeledResidue = label.AminoAcidLabel;

                    //label peptides
                    foreach (PeptideWithSetModifications pwsm in originalPeptides)
                    {
                        //find the indexes in the base sequence for labeling
                        char[] baseSequenceArray = pwsm.BaseSequence.ToArray();
                        List<int> indexesOfResiduesToBeLabeled = new List<int>();
                        for (int c = 0; c < baseSequenceArray.Length; c++)
                        {
                            if (baseSequenceArray[c] == originalResidue)
                            {
                                indexesOfResiduesToBeLabeled.Add(c);
                            }
                        }
                        //if there's something to label
                        if (indexesOfResiduesToBeLabeled.Count != 0)
                        {
                            List<PeptideWithSetModifications> pwsmsForCombinatorics = new List<PeptideWithSetModifications> { pwsm };
                            for (int a = 0; a < indexesOfResiduesToBeLabeled.Count; a++)
                            {
                                List<PeptideWithSetModifications> localPwsmsForCombinatorics = new List<PeptideWithSetModifications>();
                                foreach (PeptideWithSetModifications pwsmCombination in pwsmsForCombinatorics)
                                {
                                    char[] combinatoricBaseSequenceArray = pwsmCombination.BaseSequence.ToArray();
                                    combinatoricBaseSequenceArray[indexesOfResiduesToBeLabeled[a]] = labeledResidue;
                                    string updatedBaseSequence = string.Concat(combinatoricBaseSequenceArray);

                                    PeptideWithSetModifications labeledPwsm = new PeptideWithSetModifications(silacEndProtein, pwsm.DigestionParams,
                                        pwsm.OneBasedStartResidueInProtein, pwsm.OneBasedEndResidueInProtein, pwsm.CleavageSpecificityForFdrCategory,
                                        pwsm.PeptideDescription, pwsm.MissedCleavages, pwsm.AllModsOneIsNterminus, pwsm.NumFixedMods, updatedBaseSequence);
                                    yield return labeledPwsm; //return
                                    localPwsmsForCombinatorics.Add(labeledPwsm); //add so it can be used again
                                }
                                pwsmsForCombinatorics.AddRange(localPwsmsForCombinatorics);
                            }
                        }
                    }
                }
                else //if there are more than one (i.e. K and R are labeled)
                {
                    //get the residues to change
                    char[] originalResidues = new char[label.AdditionalLabels.Count + 1];
                    char[] labeledResidues = new char[label.AdditionalLabels.Count + 1];
                    originalResidues[0] = label.OriginalAminoAcid;
                    labeledResidues[0] = label.AminoAcidLabel;
                    for (int i = 0; i < label.AdditionalLabels.Count; i++)
                    {
                        originalResidues[i + 1] = label.AdditionalLabels[i].OriginalAminoAcid;
                        labeledResidues[i + 1] = label.AdditionalLabels[i].AminoAcidLabel;
                    }

                    //label peptides
                    foreach (PeptideWithSetModifications pwsm in originalPeptides)
                    {
                        //find the indexes in the base sequence for labeling
                        char[] baseSequenceArray = pwsm.BaseSequence.ToArray();
                        Dictionary<int, char> indexesOfResiduesToBeLabeled = new Dictionary<int, char>();
                        for (int peptideResidueIndex = 0; peptideResidueIndex < baseSequenceArray.Length; peptideResidueIndex++)
                        {
                            for (int silacResidue = 0; silacResidue < originalResidues.Length; silacResidue++)
                            {
                                if (baseSequenceArray[peptideResidueIndex] == originalResidues[silacResidue])
                                {
                                    indexesOfResiduesToBeLabeled.Add(peptideResidueIndex, labeledResidues[silacResidue]);
                                }
                            }
                        }
                        //if there's something to label
                        if (indexesOfResiduesToBeLabeled.Count != 0)
                        {
                            List<PeptideWithSetModifications> pwsmsForCombinatorics = new List<PeptideWithSetModifications> { pwsm };
                            foreach (KeyValuePair<int, char> kvp in indexesOfResiduesToBeLabeled)
                            {
                                List<PeptideWithSetModifications> localPwsmsForCombinatorics = new List<PeptideWithSetModifications>();
                                foreach (PeptideWithSetModifications pwsmCombination in pwsmsForCombinatorics)
                                {
                                    char[] combinatoricBaseSequenceArray = pwsmCombination.BaseSequence.ToArray();
                                    combinatoricBaseSequenceArray[kvp.Key] = kvp.Value;
                                    string updatedBaseSequence = string.Concat(combinatoricBaseSequenceArray);

                                    PeptideWithSetModifications labeledPwsm = new PeptideWithSetModifications(silacEndProtein, pwsm.DigestionParams,
                                        pwsm.OneBasedStartResidueInProtein, pwsm.OneBasedEndResidueInProtein, pwsm.CleavageSpecificityForFdrCategory,
                                        pwsm.PeptideDescription, pwsm.MissedCleavages, pwsm.AllModsOneIsNterminus, pwsm.NumFixedMods, updatedBaseSequence);
                                    yield return labeledPwsm; //return
                                    localPwsmsForCombinatorics.Add(labeledPwsm); //add so it can be used again
                                }
                                pwsmsForCombinatorics.AddRange(localPwsmsForCombinatorics);
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Only keep glycopeptides by filtering the NGlycopeptide motif 'NxS || NxT' or OGlycopeptide motif 'S || T'
        /// </summary>
        internal IEnumerable<ProteolyticPeptide> GetGlycoPeptides(IEnumerable<ProteolyticPeptide> originalPeptides, bool keepNGlycopeptide, bool keepOGlycopeptide)
        {
            Regex rgx = new Regex("N[A-Z][ST]");
            foreach (ProteolyticPeptide pwsm in originalPeptides)
            {
                bool yielded = false;
                if (keepNGlycopeptide)
                {
                    if (rgx.IsMatch(pwsm.BaseSequence))
                    {
                        yielded = true;
                        yield return pwsm;
                    }
                }

                if (keepOGlycopeptide && !yielded)
                {
                    if (pwsm.BaseSequence.Contains('S') || pwsm.BaseSequence.Contains('T'))
                    {
                        yield return pwsm;
                    }
                }
            }
        }

        /// <summary>
        /// Generates a protein that is fully labeled with the specified silac label
        /// </summary>
        private Protein GenerateFullyLabeledSilacProtein(SilacLabel label)
        {
            string updatedBaseSequence = BaseSequence.Replace(label.OriginalAminoAcid, label.AminoAcidLabel);
            if (label.AdditionalLabels != null) //if there is more than one label per replicate (i.e both R and K were labeled in a sample before pooling)
            {
                foreach (SilacLabel additionalLabel in label.AdditionalLabels)
                {
                    updatedBaseSequence = updatedBaseSequence.Replace(additionalLabel.OriginalAminoAcid, additionalLabel.AminoAcidLabel);
                }
            }
            return CloneWithNewSequenceAndMods(updatedBaseSequence, null) as Protein;
        }

        #region Sequence Variants

        public IBioPolymer ConsensusVariant => NonVariantProtein;

        /// <summary>
        /// Protein before applying variations.
        /// </summary>
        public Protein NonVariantProtein { get; }

        /// <summary>
        /// Sequence variations that have been applied to the base sequence.
        /// </summary>
        public List<SequenceVariation> AppliedSequenceVariations { get; }

        /// <summary>
        /// Sample name from which applied variants came, e.g. tumor or normal.
        /// </summary>
        public string SampleNameForVariants { get; }

        /// <summary>
        /// Original modifications as defined in the Parsed XML database
        /// </summary>
        public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }

        public TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original, IEnumerable<SequenceVariation> appliedSequenceVariants,
            IEnumerable<TruncationProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
            where TBioPolymerType : IHasSequenceVariants
        {
            if (original is not Protein originalProtein)
                throw new ArgumentException("The original BioPolymer must be Protein to create a protein variant");

            var variantProtein =  new Protein(variantBaseSequence, originalProtein, appliedSequenceVariants, 
                applicableProteolysisProducts, oneBasedModifications, sampleNameForVariants);
            return (TBioPolymerType)(IHasSequenceVariants)variantProtein;
        }
        #endregion

        #region Truncation Products

        /// <summary>
        /// Protein XML files contain annotated proteolysis products for many proteins (e.g. signal peptides, chain peptides).
        /// This method adds N- and C-terminal truncations to these products.
        /// </summary>

        public void AddTruncationsToExistingProteolysisProducts(int fullProteinOneBasedBegin, int fullProteinOneBasedEnd, bool addNterminalDigestionTruncations, bool addCterminalDigestionTruncations, int minProductBaseSequenceLength, int lengthOfProteolysis, string proteolyisisProductName)
        {
            bool sequenceContainsNterminus = (fullProteinOneBasedBegin == 1);

            if (sequenceContainsNterminus)
            {
                //Digest N-terminus
                if (addNterminalDigestionTruncations)
                {
                    if (BaseSequence.Substring(0, 1) == "M")
                    {
                        AddNterminalTruncations(lengthOfProteolysis + 1, fullProteinOneBasedBegin, fullProteinOneBasedEnd, minProductBaseSequenceLength, proteolyisisProductName);
                    }
                    else
                    {
                        AddNterminalTruncations(lengthOfProteolysis, fullProteinOneBasedBegin, fullProteinOneBasedEnd, minProductBaseSequenceLength, proteolyisisProductName);
                    }
                }
                //Digest C-terminus -- not effected by variable N-terminus behavior
                if (addCterminalDigestionTruncations)
                {
                    // if first residue is M, then we have to add c-terminal markers for both with and without the M
                    if (BaseSequence.Substring(0, 1) == "M")
                    {
                        //add sequences WITHOUT methionine
                        AddCterminalTruncations(lengthOfProteolysis, fullProteinOneBasedEnd, fullProteinOneBasedBegin + 1, minProductBaseSequenceLength, proteolyisisProductName);
                    }
                    //add sequences with methionine
                    AddCterminalTruncations(lengthOfProteolysis, fullProteinOneBasedEnd, fullProteinOneBasedBegin, minProductBaseSequenceLength, proteolyisisProductName);
                }
            }
            else // sequence does not contain N-terminus
            {
                //Digest C-terminus
                if (addCterminalDigestionTruncations)
                {
                    AddCterminalTruncations(lengthOfProteolysis, fullProteinOneBasedEnd, fullProteinOneBasedBegin, minProductBaseSequenceLength, proteolyisisProductName);
                }

                //Digest N-terminus
                if (addNterminalDigestionTruncations)
                {
                    AddNterminalTruncations(lengthOfProteolysis, fullProteinOneBasedBegin, fullProteinOneBasedEnd, minProductBaseSequenceLength, proteolyisisProductName);
                }
            }
        }
        /// <summary>
        /// Returns of list of proteoforms with the specified number of C-terminal amino acid truncations subject to minimum length criteria
        /// </summary>
        private void AddCterminalTruncations(int lengthOfProteolysis, int fullProteinOneBasedEnd, int fullProteinOneBasedBegin, int minProductBaseSequenceLength, string proteolyisisProductName)
        {
            for (int i = 1; i <= lengthOfProteolysis; i++)
            {
                int newEnd = fullProteinOneBasedEnd - i;
                int length = newEnd - fullProteinOneBasedBegin + 1;
                if (length >= minProductBaseSequenceLength)
                {
                    _proteolysisProducts.Add(new TruncationProduct(fullProteinOneBasedBegin, newEnd, proteolyisisProductName));
                }
            }
        }
        /// <summary>
        /// Returns of list of proteoforms with the specified number of N-terminal amino acid truncations subject to minimum length criteria
        /// </summary>

        private void AddNterminalTruncations(int lengthOfProteolysis, int fullProteinOneBasedBegin, int fullProteinOneBasedEnd, int minProductBaseSequenceLength, string proteolyisisProductName)
        {
            for (int i = 1; i <= lengthOfProteolysis; i++)
            {
                int newBegin = fullProteinOneBasedBegin + i;
                int length = fullProteinOneBasedEnd - newBegin + 1;
                if (length >= minProductBaseSequenceLength)
                {
                    _proteolysisProducts.Add(new TruncationProduct(newBegin, fullProteinOneBasedEnd, proteolyisisProductName));
                }
            }
        }

        /// <summary>
        /// This the main entry point for adding sequences in a top-down truncation search.
        /// The way this is designed is such at all base sequences to be searched end up in the list Protein.TruncationProducts
        /// This includes the intact protein. IT DOES NOT INCLUDE ANY DOUBLY (BOTH ENDS) DIGESTED PRODUCTS.
        /// The original proteolysis products (if any) are already in that list. These are annotated in protein.xml files.
        /// The options to keep in mind are present in the following variables
        /// </summary>
        /// <param name="addFullProtein"> This needs to be added to the proteolysisProducts list to be searched </param>
        /// <param name="addForEachOrigninalProteolysisProduct"> the original products are there but those resulting from N- or C-terminal degradation still need to be added</param>
        /// <param name="addNterminalDigestionTruncations"></param>
        /// <param name="addCterminalDigestionTruncations"></param>
        /// <param name="minProductBaseSequenceLength"> the same as the min detectable peptide</param>
        /// <param name="lengthOfProteolysis"> the number of amino acids that can be removed from either end.</param>
        public void AddTruncations(bool addFullProtein = true, bool addForEachOrigninalProteolysisProduct = true, bool addNterminalDigestionTruncations = true, bool addCterminalDigestionTruncations = true, int minProductBaseSequenceLength = 7, int lengthOfProteolysis = 5)
        {
            if (addFullProtein) //this loop adds the intact protoeoform and its proteolysis products to the proteolysis products list
            {
                AddIntactProteoformToTruncationsProducts(minProductBaseSequenceLength);
                if (addNterminalDigestionTruncations)
                {
                    AddTruncationsToExistingProteolysisProducts(1, BaseSequence.Length, true, false, minProductBaseSequenceLength, lengthOfProteolysis, "full-length proteoform N-terminal digestion truncation");
                }
                if (addCterminalDigestionTruncations)
                {
                    AddTruncationsToExistingProteolysisProducts(1, BaseSequence.Length, false, true, minProductBaseSequenceLength, lengthOfProteolysis, "full-length proteoform C-terminal digestion truncation");
                }
            }

            if (addForEachOrigninalProteolysisProduct) // this does not include the original intact proteoform
            {
                List<TruncationProduct> existingProducts = TruncationProducts.Where(p => !p.Type.Contains("truncation") && !p.Type.Contains("full-length proteoform")).ToList();
                foreach (TruncationProduct product in existingProducts)
                {
                    if (product.OneBasedBeginPosition.HasValue && product.OneBasedEndPosition.HasValue)
                    {
                        string proteolyisisProductName = "truncation";

                        if (!String.IsNullOrEmpty(product.Type))
                        {
                            proteolyisisProductName = product.Type + " " + proteolyisisProductName;
                        }
                        //the original proteolysis product is already on the list so we don't need to duplicate
                        if (addNterminalDigestionTruncations)
                        {
                            AddTruncationsToExistingProteolysisProducts(product.OneBasedBeginPosition.Value, product.OneBasedEndPosition.Value, true, false, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName);
                        }
                        if (addCterminalDigestionTruncations)
                        {
                            AddTruncationsToExistingProteolysisProducts(product.OneBasedBeginPosition.Value, product.OneBasedEndPosition.Value, false, true, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName);
                        }
                    }
                }
            }
            CleaveOnceBetweenProteolysisProducts();
        }
        /// <summary>
        /// This method adds proteoforms with N- and C-terminal amino acid loss to the list of species included in top-down search
        /// </summary>
        public void AddIntactProteoformToTruncationsProducts(int minProductBaseSequenceLength)
        {
            if (BaseSequence.Length >= minProductBaseSequenceLength)
            {
                _proteolysisProducts.Add(new TruncationProduct(1, BaseSequence.Length, "full-length proteoform"));
            }
        }

        /// <summary>
        /// proteins with multiple proteolysis products are not always full cleaved. we observed proteolysis products w/ missed cleavages.
        /// This method allows for one missed cleavage between proteolysis products.
        /// </summary>

        public void CleaveOnceBetweenProteolysisProducts(int minimumProductLength = 7)
        {
            List<int> cleavagePostions = new();
            List<TruncationProduct> localProducts = _proteolysisProducts.Where(p => !p.Type.Contains("truncation") && !p.Type.Contains("full-length proteoform")).ToList();
            List<int> proteolysisProductEndPositions = localProducts.Where(p => p.OneBasedEndPosition.HasValue).Select(p => p.OneBasedEndPosition.Value).ToList();
            if (proteolysisProductEndPositions.Count > 0)
            {
                foreach (int proteolysisProductEndPosition in proteolysisProductEndPositions)
                {
                    if (localProducts.Any(p => p.OneBasedBeginPosition == (proteolysisProductEndPosition + 1)))
                    {
                        cleavagePostions.Add(proteolysisProductEndPosition);
                    }
                }
            }

            foreach (int position in cleavagePostions)
            {
                if (position - 1 >= minimumProductLength)
                {
                    string leftType = $"N-terminal Portion of Singly Cleaved Protein(1-{position})";
                    TruncationProduct leftProduct = new(1, position, leftType);

                    //here we're making sure a product with these begin/end positions isn't already present
                    if (!_proteolysisProducts.Any(p => p.OneBasedBeginPosition == leftProduct.OneBasedBeginPosition && p.OneBasedEndPosition == leftProduct.OneBasedEndPosition))
                    {
                        _proteolysisProducts.Add(leftProduct);
                    }
                }

                if (BaseSequence.Length - position - 1 >= minimumProductLength)
                {
                    string rightType = $"C-terminal Portion of Singly Cleaved Protein({position + 1}-{BaseSequence.Length})";
                    TruncationProduct rightProduct = new(position + 1, BaseSequence.Length, rightType);

                    //here we're making sure a product with these begin/end positions isn't already present
                    if (!_proteolysisProducts.Any(p => p.OneBasedBeginPosition == rightProduct.OneBasedBeginPosition && p.OneBasedEndPosition == rightProduct.OneBasedEndPosition))
                    {
                        _proteolysisProducts.Add(rightProduct);
                    }
                }
            }
        }

        #endregion

        public IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence, IDictionary<int, List<Modification>>? newMods = null)
        {
            // Create a new protein with the new base sequence and modifications
            Protein newProtein = new Protein(this, newBaseSequence);
            if (newMods.IsNullOrEmpty()) 
                return newProtein;

            // If new modifications are provided, use them
            newProtein.OriginalNonVariantModifications = this.OriginalNonVariantModifications;
            newProtein.OneBasedPossibleLocalizedModifications = ((IBioPolymer)newProtein).SelectValidOneBaseMods(newMods!);

            return newProtein;
        }

        public int CompareTo(Protein other)
        {
            //permits sorting of proteins
            return this.Accession.CompareTo(other.Accession);
        }

        //not sure if we require any additional fields for equality
        public override bool Equals(object obj)
        {
            if (obj is Protein bioPol)
            {
                return Equals(bioPol);
            }
            return false;
        }

        public bool Equals(Protein other)
        {
            return (this as IBioPolymer).Equals(other);   
        }
        
        /// <summary>
        /// The protein object uses the default hash code method for speed,
        /// but note that two protein objects with the same information will give two different hash codes.
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode()
        {
            return this.BaseSequence.GetHashCode();
        }

        public override string ToString()
        {
            return this.Accession.ToString();
        }
    }
}