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
        public Protein(Protein originalProtein, string silacSequence, string silacAccession)
        {
            BaseSequence = silacSequence;
            Accession = silacAccession;
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
        public IEnumerable<PeptideWithSetModifications> Digest(DigestionParams digestionParams, IEnumerable<Modification> allKnownFixedModifications,
            List<Modification> variableModifications, List<SilacLabel> silacLabels = null)
        {
            //can't be null
            allKnownFixedModifications = allKnownFixedModifications ?? new List<Modification>();
            variableModifications = variableModifications ?? new List<Modification>();
            CleavageSpecificity searchModeType = digestionParams.SearchModeType;

            ProteinDigestion digestion = new ProteinDigestion(digestionParams, allKnownFixedModifications, variableModifications);
            IEnumerable<ProteolyticPeptide> unmodifiedPeptides =
                searchModeType == CleavageSpecificity.Semi ?
                digestion.SpeedySemiSpecificDigestion(this) :
                digestion.Digestion(this);

            IEnumerable<PeptideWithSetModifications> modifiedPeptides = unmodifiedPeptides.SelectMany(peptide => peptide.GetModifiedPeptides(allKnownFixedModifications, digestionParams, variableModifications));

            //Remove terminal modifications (if needed)
            if (searchModeType == CleavageSpecificity.SingleN || searchModeType == CleavageSpecificity.SingleC)
            {
                modifiedPeptides = RemoveTerminalModifications(modifiedPeptides, searchModeType, allKnownFixedModifications);
            }

            //add silac labels (if needed)
            if (silacLabels != null)
            {
                return GetSilacPeptides(modifiedPeptides, silacLabels, digestionParams.GeneratehUnlabeledProteinsForSilac);
            }

            return modifiedPeptides;
        }

        /// <summary>
        /// Remove terminal modifications from the C-terminus of SingleN peptides and the N-terminus of SingleC peptides/
        /// These terminal modifications create redundant entries and increase search time
        /// </summary>
        internal static IEnumerable<PeptideWithSetModifications> RemoveTerminalModifications(IEnumerable<PeptideWithSetModifications> modifiedPeptides, CleavageSpecificity searchModeType, IEnumerable<Modification> allFixedMods)
        {
            string terminalStringToLookFor = searchModeType == CleavageSpecificity.SingleN ? "C-terminal" : "N-terminal";
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
        internal IEnumerable<PeptideWithSetModifications> GetSilacPeptides(IEnumerable<PeptideWithSetModifications> originalPeptides, List<SilacLabel> silacLabels, bool generateUnlabeledProteins)
        {
            //unlabeled peptides
            if (generateUnlabeledProteins)
            {
                foreach (PeptideWithSetModifications pwsm in originalPeptides)
                {
                    yield return pwsm;
                }
            }

            //labeled peptides
            foreach (SilacLabel label in silacLabels)
            {
                string updatedBaseSequence = BaseSequence.Replace(label.OriginalAminoAcid, label.AminoAcidLabel);
                string updatedAccession = Accession + label.MassDifference;
                if (label.AdditionalLabels != null) //if there is more than one label per replicate (i.e both R and K were labeled in a sample before pooling)
                {
                    foreach (SilacLabel additionalLabel in label.AdditionalLabels)
                    {
                        updatedBaseSequence = updatedBaseSequence.Replace(additionalLabel.OriginalAminoAcid, additionalLabel.AminoAcidLabel);
                        updatedAccession += additionalLabel.MassDifference;
                    }
                }
                Protein silacProtein = new Protein(this, updatedBaseSequence, updatedAccession);

                foreach (PeptideWithSetModifications pwsm in originalPeptides)
                {
                    //duplicate the peptides with the updated protein sequence that contains only silac labels
                    yield return new PeptideWithSetModifications(silacProtein, pwsm.DigestionParams, pwsm.OneBasedStartResidueInProtein, pwsm.OneBasedEndResidueInProtein, pwsm.CleavageSpecificityForFdrCategory, pwsm.PeptideDescription, pwsm.MissedCleavages, pwsm.AllModsOneIsNterminus, pwsm.NumFixedMods);
                }
            }
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

            if (otherProtein == null)
            {
                return false;
            }

            return otherProtein != null && otherProtein.Accession == this.Accession && otherProtein.BaseSequence == this.BaseSequence;
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