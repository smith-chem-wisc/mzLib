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
            List<SequenceVariation> sequenceVariations = null, List<DisulfideBond> disulfideBonds = null, List<SpliceSite> spliceSites = null, string databaseFilePath = null)
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

            GeneNames = geneNames ?? new List<Tuple<string, string>>();
            ProteolysisProducts = proteolysisProducts ?? new List<ProteolysisProduct>();
            SequenceVariations = sequenceVariations ?? new List<SequenceVariation>();
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
        public List<SequenceVariation> AppliedSequenceVariations { get; } = new List<SequenceVariation>();

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

        //public override bool Equals(object obj)
        //{
        //    Protein p = obj as Protein;

        //    return p != null
        //        && p.BaseSequence == BaseSequence
        //        && p.Name == Name
        //        && p.Accession == Accession
        //        && p.FullName == FullName
        //        && p.FullDescription == FullDescription
        //        && p.IsContaminant == IsContaminant
        //        && p.IsDecoy == IsDecoy
        //        && p.Organism == Organism
        //        && p.SampleNameForVariants == SampleNameForVariants
        //        && p.GeneNames.OrderBy(x => $"{x.Item1}{x.Item2}").SequenceEqual(GeneNames.OrderBy(x => $"{x.Item1}{x.Item2}"))
        //        && p.SequenceVariations.OrderBy(x => $"{x.SimpleString()}{x.OneBasedEndPosition}{x.Description.ToString()}").SequenceEqual(SequenceVariations.OrderBy(x => $"{x.SimpleString()}{x.OneBasedEndPosition}{x.Description.ToString()}"))
        //        && p.SpliceSites.OrderBy(x => $"{x.OneBasedBeginPosition}-{x.OneBasedEndPosition}:{x.Description}").SequenceEqual(SpliceSites.OrderBy(x => $"{x.OneBasedBeginPosition}-{x.OneBasedEndPosition}:{x.Description}"))
        //        && p.AppliedSequenceVariations.OrderBy(x => $"{x.SimpleString()}{x.OneBasedEndPosition}{x.Description.ToString()}").SequenceEqual(AppliedSequenceVariations.OrderBy(x => $"{x.SimpleString()}{x.OneBasedEndPosition}{x.Description.ToString()}"))
        //        && p.DatabaseReferences.OrderBy(x => x).SequenceEqual(DatabaseReferences.OrderBy(x => x))
        //        && p.DisulfideBonds.OrderBy(x => x).SequenceEqual(DisulfideBonds.OrderBy(x => x))
        //        && p.ProteolysisProducts.OrderBy(x => x).SequenceEqual(ProteolysisProducts.OrderBy(x => x))
        //        && p.OneBasedPossibleLocalizedModifications.OrderBy(x => x.Key).SelectMany(x => $"{x.Key.ToString()}{string.Join("", x.Value.Select(mod => mod.ToString()).OrderBy(mod => mod))}")
        //            .SequenceEqual(OneBasedPossibleLocalizedModifications.OrderBy(x => x.Key).SelectMany(x => $"{x.Key.ToString()}{string.Join("", x.Value.Select(mod => mod.ToString()).OrderBy(mod => mod))}"));
        //}

        public override bool Equals(object obj)
        {
            Protein p = obj as Protein;

            return p != null
                && p.BaseSequence == BaseSequence
                && p.Name == Name
                && p.Accession == Accession
                && p.FullName == FullName
                && p.FullDescription == FullDescription
                && p.IsContaminant == IsContaminant
                && p.IsDecoy == IsDecoy
                && p.Organism == Organism
                && p.SampleNameForVariants == SampleNameForVariants
                && p.GeneNames.SequenceEqual(GeneNames)
                && p.SequenceVariations.SequenceEqual(SequenceVariations)
                && p.SpliceSites.SequenceEqual(SpliceSites)
                && p.AppliedSequenceVariations.SequenceEqual(AppliedSequenceVariations)
                && p.DatabaseReferences.SequenceEqual(DatabaseReferences)
                && p.DisulfideBonds.SequenceEqual(DisulfideBonds)
                && p.ProteolysisProducts.SequenceEqual(ProteolysisProducts)
                && p.OneBasedPossibleLocalizedModifications.OrderBy(x => x.Key).SelectMany(x => $"{x.Key.ToString()}{string.Join("", x.Value.Select(mod => mod.ToString()).OrderBy(mod => mod))}")
                    .SequenceEqual(OneBasedPossibleLocalizedModifications.OrderBy(x => x.Key).SelectMany(x => $"{x.Key.ToString()}{string.Join("", x.Value.Select(mod => mod.ToString()).OrderBy(mod => mod))}"));
        }

        public override int GetHashCode()
        {
            int hash = (BaseSequence ?? "").GetHashCode();
            hash ^= (Name ?? "").GetHashCode();
            hash ^= (Accession ?? "").GetHashCode();
            hash ^= (FullName ?? "").GetHashCode();
            hash ^= (FullDescription ?? "").GetHashCode();
            hash ^= IsContaminant.GetHashCode();
            hash ^= IsDecoy.GetHashCode();
            hash ^= (Organism ?? "").GetHashCode();
            hash ^= (SampleNameForVariants ?? "").GetHashCode();

            foreach (Tuple<string, string> gn in GeneNames)
            {
                hash ^= gn.GetHashCode();
            }
            foreach (SequenceVariation sv in SequenceVariations)
            {
                hash ^= sv.GetHashCode();
            }
            foreach (SequenceVariation sv in AppliedSequenceVariations)
            {
                hash ^= sv.GetHashCode();
            }
            foreach (DatabaseReference dr in DatabaseReferences)
            {
                hash ^= dr.GetHashCode();
            }
            foreach (DisulfideBond db in DisulfideBonds)
            {
                hash ^= db.GetHashCode();
            }
            foreach (ProteolysisProduct pp in ProteolysisProducts)
            {
                hash ^= pp.GetHashCode();
            }
            foreach (var kv in OneBasedPossibleLocalizedModifications)
            {
                foreach (Modification mod in kv.Value)
                {
                    hash ^= kv.Key.GetHashCode() ^ mod.GetHashCode();
                }
            }

            return hash;
        }

        /// <summary>
        /// Gets peptides for digestion of a protein
        /// </summary>
        public IEnumerable<PeptideWithSetModifications> Digest(DigestionParams digestionParams, IEnumerable<Modification> allKnownFixedModifications,
            List<Modification> variableModifications)
        {
            ProteinDigestion digestion = new ProteinDigestion(digestionParams, allKnownFixedModifications, variableModifications);
            return digestionParams.SearchModeType == CleavageSpecificity.Semi ? digestion.SpeedySemiSpecificDigestion(this) : digestion.Digestion(this);
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
    }
}