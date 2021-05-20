using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public class Protein
    {
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
            IDictionary<int, List<Modification>> oneBasedModifications = null, List<ProteolysisProduct> proteolysisProducts = null,
            string name = null, string fullName = null, bool isDecoy = false, bool isContaminant = false, List<DatabaseReference> databaseReferences = null,
            List<SequenceVariation> sequenceVariations = null, List<SequenceVariation> appliedSequenceVariations = null, string sampleNameForVariants = null,
            List<DisulfideBond> disulfideBonds = null, List<SpliceSite> spliceSites = null, string databaseFilePath = null)
        {
            // Mandatory
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
            ProteolysisProducts = proteolysisProducts ?? new List<ProteolysisProduct>();
            SequenceVariations = sequenceVariations ?? new List<SequenceVariation>();
            AppliedSequenceVariations = appliedSequenceVariations ?? new List<SequenceVariation>();
            OriginalNonVariantModifications = oneBasedModifications ?? new Dictionary<int, List<Modification>>();
            if (oneBasedModifications != null)
            {
                OneBasedPossibleLocalizedModifications = SelectValidOneBaseMods(oneBasedModifications);
            }
            else
            {
                OneBasedPossibleLocalizedModifications = new Dictionary<int, List<Modification>>();
            }
            DatabaseReferences = databaseReferences ?? new List<DatabaseReference>();
            DisulfideBonds = disulfideBonds ?? new List<DisulfideBond>();
            SpliceSites = spliceSites ?? new List<SpliceSite>();
        }

        /// <summary>
        /// Protein construction that clones a protein but assigns a different base sequence
        /// For use in SILAC experiments
        /// </summary>
        /// <param name="originalProtein"></param>
        /// <param name="silacSequence"></param>
        /// <param name="silacAccession"></param>
        public Protein(Protein originalProtein, string silacSequence)
        {
            BaseSequence = silacSequence;
            Accession = originalProtein.Accession;
            NonVariantProtein = originalProtein.NonVariantProtein;
            Name = originalProtein.Name;
            Organism = originalProtein.Organism;
            FullName = originalProtein.FullName;
            IsDecoy = originalProtein.IsDecoy;
            IsContaminant = originalProtein.IsContaminant;
            DatabaseFilePath = originalProtein.DatabaseFilePath;
            SampleNameForVariants = originalProtein.SampleNameForVariants;
            GeneNames = originalProtein.GeneNames;
            ProteolysisProducts = originalProtein.ProteolysisProducts;
            SequenceVariations = originalProtein.SequenceVariations;
            AppliedSequenceVariations = originalProtein.AppliedSequenceVariations;
            OriginalNonVariantModifications = originalProtein.OriginalNonVariantModifications;
            OneBasedPossibleLocalizedModifications = originalProtein.OneBasedPossibleLocalizedModifications;
            DatabaseReferences = originalProtein.DatabaseReferences;
            DisulfideBonds = originalProtein.DisulfideBonds;
            SpliceSites = originalProtein.SpliceSites;
            DatabaseFilePath = originalProtein.DatabaseFilePath;
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
            IEnumerable<ProteolysisProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
            : this(variantBaseSequence,
                  VariantApplication.GetAccession(protein, appliedSequenceVariations),
                  organism: protein.Organism,
                  geneNames: new List<Tuple<string, string>>(protein.GeneNames),
                  oneBasedModifications: oneBasedModifications != null ? oneBasedModifications.ToDictionary(x => x.Key, x => x.Value) : new Dictionary<int, List<Modification>>(),
                  proteolysisProducts: new List<ProteolysisProduct>(applicableProteolysisProducts ?? new List<ProteolysisProduct>()),
                  name: GetName(appliedSequenceVariations, protein.Name),
                  fullName: GetName(appliedSequenceVariations, protein.FullName),
                  isDecoy: protein.IsDecoy,
                  isContaminant: protein.IsContaminant,
                  databaseReferences: new List<DatabaseReference>(protein.DatabaseReferences),
                  sequenceVariations: new List<SequenceVariation>(protein.SequenceVariations),
                  disulfideBonds: new List<DisulfideBond>(protein.DisulfideBonds),
                  spliceSites: new List<SpliceSite>(protein.SpliceSites),
                  databaseFilePath: protein.DatabaseFilePath)
        {
            NonVariantProtein = protein.NonVariantProtein;
            OriginalNonVariantModifications = NonVariantProtein.OriginalNonVariantModifications;
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
        public IEnumerable<Tuple<string, string>> GeneNames { get; }

        /// <summary>
        /// Unique accession for this protein.
        /// </summary>
        public string Accession { get; }

        /// <summary>
        /// Base sequence, which may contain applied sequence variations.
        /// </summary>
        public string BaseSequence { get; }

        public string Organism { get; }
        public bool IsDecoy { get; }
        public IEnumerable<SequenceVariation> SequenceVariations { get; }
        public IEnumerable<DisulfideBond> DisulfideBonds { get; }
        public IEnumerable<SpliceSite> SpliceSites { get; }
        //TODO: Generate all the proteolytic products as distinct proteins during XML reading and delete the ProteolysisProducts parameter
        public IEnumerable<ProteolysisProduct> ProteolysisProducts { get; }
        public IEnumerable<DatabaseReference> DatabaseReferences { get; }
        public string DatabaseFilePath { get; }

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

        public double Probability { get; set; } // for protein pep project

        public int Length
        {
            get
            {
                return BaseSequence.Length;
            }
        }

        public string FullDescription
        {
            get
            {
                return Accession + "|" + Name + "|" + FullName;
            }
        }

        public string Name { get; }
        public string FullName { get; }
        public bool IsContaminant { get; }
        internal IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }

        public char this[int zeroBasedIndex]
        {
            get
            {
                return BaseSequence[zeroBasedIndex];
            }
        }

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

        /// <summary>
        /// Gets peptides for digestion of a protein
        /// </summary>
        public IEnumerable<PeptideWithSetModifications> Digest(DigestionParams digestionParams, List<Modification> allKnownFixedModifications,
            List<Modification> variableModifications, List<SilacLabel> silacLabels = null, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null)
        {            
            //can't be null
            allKnownFixedModifications = allKnownFixedModifications ?? new List<Modification>(); 
            // add in any modifications that are caused by protease digestion
            if (digestionParams.Protease.CleavageMod!= null && !allKnownFixedModifications.Contains(digestionParams.Protease.CleavageMod))
            {
                allKnownFixedModifications.Add(digestionParams.Protease.CleavageMod);                
            }                      
            variableModifications = variableModifications ?? new List<Modification>();
            CleavageSpecificity searchModeType = digestionParams.SearchModeType;

            ProteinDigestion digestion = new ProteinDigestion(digestionParams, allKnownFixedModifications, variableModifications);
            IEnumerable<ProteolyticPeptide> unmodifiedPeptides =
                searchModeType == CleavageSpecificity.Semi ?
                digestion.SpeedySemiSpecificDigestion(this) :
                digestion.Digestion(this);

            IEnumerable<PeptideWithSetModifications> modifiedPeptides = unmodifiedPeptides.SelectMany(peptide => peptide.GetModifiedPeptides(allKnownFixedModifications, digestionParams, variableModifications));

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
               return GetSilacPeptides(modifiedPeptides, silacLabels, digestionParams.GeneratehUnlabeledProteinsForSilac, turnoverLabels);
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
                    if(startLabel.AdditionalLabels!=null)
                    {
                        foreach(SilacLabel originalLabel in originalLabels)
                        {
                            if(updatedEndLabel.OriginalAminoAcid!= originalLabel.AminoAcidLabel && 
                                (updatedEndLabel.AdditionalLabels==null || !updatedEndLabel.AdditionalLabels.Any(x=>x.OriginalAminoAcid == originalLabel.AminoAcidLabel)))
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
            return new Protein(this, updatedBaseSequence);
        }

        /// <summary>
        /// Gets proteins with applied variants from this protein
        /// </summary>
        public List<Protein> GetVariantProteins(int maxAllowedVariantsForCombinitorics = 4, int minAlleleDepth = 1)
        {
            return VariantApplication.ApplyVariants(this, SequenceVariations, maxAllowedVariantsForCombinitorics, minAlleleDepth);
        }

        /// <summary>
        /// Restore all modifications that were read in, including those that did not match their target amino acid.
        /// </summary>
        public void RestoreUnfilteredModifications()
        {
            OneBasedPossibleLocalizedModifications = OriginalNonVariantModifications;
        }

        /// <summary>
        /// Filters modifications that do not match their target amino acid.
        /// </summary>
        /// <param name="dict"></param>
        /// <returns></returns>
        private IDictionary<int, List<Modification>> SelectValidOneBaseMods(IDictionary<int, List<Modification>> dict)
        {
            Dictionary<int, List<Modification>> validModDictionary = new Dictionary<int, List<Modification>>();
            foreach (KeyValuePair<int, List<Modification>> entry in dict)
            {
                List<Modification> validMods = new List<Modification>();
                foreach (Modification m in entry.Value)
                {
                    //mod must be valid mod and the motif of the mod must be present in the protein at the specified location
                    if (m.ValidModification && ModificationLocalization.ModFits(m, BaseSequence, 0, BaseSequence.Length, entry.Key))
                    {
                        validMods.Add(m);
                    }
                }

                if (validMods.Any())
                {
                    if (validModDictionary.Keys.Contains(entry.Key))
                    {
                        validModDictionary[entry.Key].AddRange(validMods);
                    }
                    else
                    {
                        validModDictionary.Add(entry.Key, validMods);
                    }
                }
            }
            return validModDictionary;
        }

        private static string GetName(IEnumerable<SequenceVariation> appliedVariations, string name)
        {
            bool emptyVars = appliedVariations == null || appliedVariations.Count() == 0;
            if (name == null && emptyVars)
            {
                return null;
            }
            else
            {
                string variantTag = emptyVars ? "" : $" variant:{VariantApplication.CombineDescriptions(appliedVariations)}";
                return name + variantTag;
            }
        }

        public int CompareTo(Protein other)
        {
            //permits sorting of proteins
            return this.Accession.CompareTo(other.Accession);
        }

        //not sure if we require any additional fields for equality
        public override bool Equals(object obj)
        {
            Protein otherProtein = (Protein)obj;
            return otherProtein != null && otherProtein.Accession.Equals(Accession) && otherProtein.BaseSequence.Equals(BaseSequence);
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