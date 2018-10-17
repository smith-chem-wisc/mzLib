using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public class ProteinWithAppliedVariants
        : Protein
    {
        public Protein Protein { get; }
        public List<SequenceVariation> AppliedSequenceVariations { get; }
        public string Individual { get; }

        public ProteinWithAppliedVariants(string variantBaseSequence, Protein protein, IEnumerable<SequenceVariation> appliedSequenceVariations,
            IEnumerable<ProteolysisProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications, string individual)
            : base(variantBaseSequence,
                  GetAccession(protein, appliedSequenceVariations),
                  organism: protein.Organism,
                  geneNames: new List<Tuple<string, string>>(protein.GeneNames),
                  oneBasedModifications: oneBasedModifications.ToDictionary(x => x.Key, x => x.Value),
                  proteolysisProducts: new List<ProteolysisProduct>(applicableProteolysisProducts),
                  name: protein.Name + (appliedSequenceVariations == null ? "" : " variant:" + CombineDescriptions(appliedSequenceVariations)),
                  fullName: protein.FullName + (appliedSequenceVariations == null ? "" : " variant:" + CombineDescriptions(appliedSequenceVariations)),
                  isDecoy: protein.IsDecoy,
                  isContaminant: protein.IsContaminant,
                  databaseReferences: new List<DatabaseReference>(protein.DatabaseReferences),
                  sequenceVariations: new List<SequenceVariation>(protein.SequenceVariations),
                  disulfideBonds: new List<DisulfideBond>(protein.DisulfideBonds),
                  databaseFilePath: protein.DatabaseFilePath)
        {
            Protein = protein;
            AppliedSequenceVariations = appliedSequenceVariations != null ? appliedSequenceVariations.ToList() : new List<SequenceVariation>();
            Individual = individual;
        }

        public override bool Equals(object obj)
        {
            ProteinWithAppliedVariants p = obj as ProteinWithAppliedVariants;
            return p != null &&
                p.Individual.Equals(Individual) &&
                p.SequenceVariations.OrderBy(x => x).SequenceEqual(SequenceVariations.OrderBy(x => x)) &&
                p.Protein.Equals(Protein);
        }

        public override int GetHashCode()
        {
            int hash = Individual.GetHashCode() ^ Protein.GetHashCode();
            foreach (SequenceVariation sv in AppliedSequenceVariations)
            {
                hash ^= sv.GetHashCode();
            }
            return hash;
        }

        /// <summary>
        /// Applies variant changes to protein sequence
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="uniqueEffectsToApply"></param>
        /// <returns></returns>
        internal List<ProteinWithAppliedVariants> ApplyVariants(IEnumerable<SequenceVariation> sequenceVariations, int maxAllowedVariantsForCombinitorics)
        {
            List<SequenceVariation> uniqueEffectsToApply = sequenceVariations
                .GroupBy(v => v.SimpleString())
                .Select(x => x.First())
                .Where(v => v.Description.Genotypes.Count > 0) // this is a VCF line
                .OrderByDescending(v => v.OneBasedBeginPosition) // apply variants at the end of the protein sequence first
                .ToList();

            ProteinWithAppliedVariants proteinCopy = new ProteinWithAppliedVariants(BaseSequence, Protein, AppliedSequenceVariations, ProteolysisProducts, OneBasedPossibleLocalizedModifications, Individual);
            
            // If there aren't any variants to apply, just return the base protein
            if (uniqueEffectsToApply.Count == 0)
            {
                return new List<ProteinWithAppliedVariants> { proteinCopy };
            }

            HashSet<string> individuals = new HashSet<string>(uniqueEffectsToApply.SelectMany(v => v.Description.Genotypes.Keys));
            List<ProteinWithAppliedVariants> variantProteins = new List<ProteinWithAppliedVariants>();

            // loop through genotypes for each sample/individual (e.g. tumor and normal)
            foreach (string individual in individuals)
            {
                bool tooManyHeterozygousVariants = uniqueEffectsToApply.Count(v => v.Description.Heterozygous[individual]) > maxAllowedVariantsForCombinitorics;
                List<ProteinWithAppliedVariants> newVariantProteins = new List<ProteinWithAppliedVariants> { proteinCopy };
                foreach (var variant in uniqueEffectsToApply)
                {
                    // homozygous alternate
                    if (variant.Description.Homozygous[individual] && variant.Description.Genotypes[individual].All(x => x != "0"))
                    {
                        newVariantProteins = newVariantProteins.Select(p => ApplyVariant(variant, p, individual)).ToList();
                    }

                    // heterozygous basic
                    // first protein with variants contains all homozygous variation, second contains all variations
                    else if (variant.Description.Heterozygous[individual] && tooManyHeterozygousVariants)
                    {
                        if (newVariantProteins.Count == 1)
                        {
                            ProteinWithAppliedVariants variantProtein = ApplyVariant(variant, newVariantProteins[0], individual);
                            newVariantProteins.Add(variantProtein);
                        }
                        else
                        {
                            newVariantProteins[1] = ApplyVariant(variant, newVariantProteins[1], individual);
                        }
                    }
                    
                    // heterozygous combinitorics
                    else if (variant.Description.Heterozygous[individual] && !tooManyHeterozygousVariants)
                    {
                        List<ProteinWithAppliedVariants> combinitoricProteins = new List<ProteinWithAppliedVariants>();

                        foreach (var protein in newVariantProteins)
                        {
                            // keep reference allele
                            if (variant.Description.Genotypes[individual].Contains("0"))
                            {
                                combinitoricProteins.Add(protein);
                            }

                            // alternate allele (replace all, since in heterozygous with two alternates, both alternates are included)
                            combinitoricProteins.Add(ApplyVariant(variant, proteinCopy, individual));
                        }
                        newVariantProteins = combinitoricProteins;
                    }
                }
                variantProteins.AddRange(newVariantProteins);
            }

            return variantProteins.GroupBy(x => x.BaseSequence).Select(x => x.First()).ToList();
        }
       
        /// <summary>
        /// Applies a variant to a protein sequence
        /// </summary>
        /// <param name="variant"></param>
        /// <returns></returns>
        internal static ProteinWithAppliedVariants ApplyVariant(SequenceVariation variant, ProteinWithAppliedVariants protein, string individual)
        {
            string seqBefore = protein.BaseSequence.Substring(0, variant.OneBasedBeginPosition - 1);
            string seqVariant = variant.VariantSequence;
            List<ProteolysisProduct> adjustedProteolysisProducts = protein.AdjustProteolysisProductIndices(variant, protein.ProteolysisProducts);
            Dictionary<int, List<Modification>> adjustedModifications = protein.AdjustModificationIndices(variant, protein);
            List<SequenceVariation> adjustedAppliedVariations = protein.AdjustSequenceVariations(variant, protein.AppliedSequenceVariations);
            int afterIdx = variant.OneBasedBeginPosition + variant.OriginalSequence.Length - 1;

            // check to see if there is incomplete indel overlap, which would lead to weird variant sequences
            // complete overlap is okay, since it will be overwritten; this can happen if there are two alternate alleles,
            //    e.g. reference sequence is wrong at that point
            bool intersectsAppliedRegionIncompletely = protein.AppliedSequenceVariations.Any(x => variant.Intersects(x) && !variant.Includes(x));
            if (intersectsAppliedRegionIncompletely)
            {
                // use original protein sequence for the remaining sequence
                string seqAfter = protein.Protein.BaseSequence.Length - afterIdx <= 0 ? "" : protein.Protein.BaseSequence.Substring(afterIdx);
                return new ProteinWithAppliedVariants(seqBefore + seqVariant + seqAfter, protein.Protein, new[] { variant }, adjustedProteolysisProducts, adjustedModifications, individual);
            }
            else
            {
                List<SequenceVariation> variations = protein.AppliedSequenceVariations
                    .Where(x => !variant.Includes(x))
                    .Concat(new[] { variant })
                    .ToList();
                // use this variant protein sequence for the remaining sequence
                string seqAfter = protein.BaseSequence.Length - afterIdx <= 0 ? "" : protein.BaseSequence.Substring(afterIdx);
                return new ProteinWithAppliedVariants(seqBefore + seqVariant + seqAfter, protein.Protein, variations, adjustedProteolysisProducts, adjustedModifications, individual);
            }
        }

        internal List<SequenceVariation> AdjustSequenceVariations(SequenceVariation variantGettingApplied, IEnumerable<SequenceVariation> alreadyAppliedVariations)
        {
            List<SequenceVariation> variations = new List<SequenceVariation>();
            if (alreadyAppliedVariations == null) { return variations; }
            foreach (SequenceVariation v in alreadyAppliedVariations)
            {
                int addedIdx = alreadyAppliedVariations
                    .Where(applied => applied.OneBasedEndPosition < v.OneBasedBeginPosition)
                    .Sum(applied => applied.VariantSequence.Length - applied.OriginalSequence.Length);
                
                // variant was entirely before the one being applied (shouldn't happen because of order of applying variants)
                if (v.OneBasedEndPosition - addedIdx < variantGettingApplied.OneBasedBeginPosition)
                {
                    variations.Add(v);
                }

                // adjust indices based on new included sequencek, minding possible overlaps to be filtered later
                else
                {
                    int overlap = variantGettingApplied.OneBasedEndPosition - v.OneBasedBeginPosition + 1;
                    int sequenceLengthChange = variantGettingApplied.VariantSequence.Length - variantGettingApplied.OriginalSequence.Length;
                    variations.Add(new SequenceVariation(
                        v.OneBasedBeginPosition + sequenceLengthChange - overlap,
                        v.OneBasedEndPosition + sequenceLengthChange - overlap,
                        v.VariantSequence,
                        v.OriginalSequence,
                        v.Description.Description,
                        v.OneBasedModifications.ToDictionary(kv => kv.Key, kv => kv.Value)));
                }
            }
            return variations;
        }

        /// <summary>
        /// Eliminates proteolysis products that overlap sequence variations.
        /// Since frameshift indels are written across the remaining sequence,
        /// this eliminates proteolysis products that conflict with large deletions and other structural variations.
        /// </summary>
        /// <param name="variants"></param>
        /// <param name="proteolysisProducts"></param>
        /// <returns></returns>
        internal List<ProteolysisProduct> AdjustProteolysisProductIndices(SequenceVariation variant, IEnumerable<ProteolysisProduct> proteolysisProducts)
        {
            List<ProteolysisProduct> products = new List<ProteolysisProduct>();
            if (proteolysisProducts == null) { return products; }
            int sequenceLengthChange = variant.VariantSequence.Length - variant.OriginalSequence.Length;
            foreach (ProteolysisProduct p in proteolysisProducts.Where(p => p.OneBasedEndPosition.HasValue && p.OneBasedBeginPosition.HasValue))
            {
                // proteolysis product is entirely before the variant
                if (variant.OneBasedBeginPosition > p.OneBasedEndPosition)
                {
                    products.Add(p);
                }
                // proteolysis product straddles the variant, but the cleavage site(s) are still intact; the ends aren't considered cleavage sites
                else if ((p.OneBasedBeginPosition < variant.OneBasedBeginPosition || p.OneBasedBeginPosition == 1 || p.OneBasedBeginPosition == 2)
                    && (p.OneBasedEndPosition > variant.OneBasedEndPosition || p.OneBasedEndPosition == Protein.BaseSequence.Length)
                    && p.OneBasedEndPosition + sequenceLengthChange <= BaseSequence.Length)
                {
                    products.Add(new ProteolysisProduct(p.OneBasedBeginPosition, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
                }
                // proteolysis product is after the variant
                else if (p.OneBasedBeginPosition > variant.OneBasedEndPosition
                    && p.OneBasedBeginPosition + sequenceLengthChange <= BaseSequence.Length && p.OneBasedEndPosition + sequenceLengthChange <= BaseSequence.Length)
                {
                    products.Add(new ProteolysisProduct(p.OneBasedBeginPosition + sequenceLengthChange, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
                }
                else // sequence variant conflicts with proteolysis cleavage site (cleavage site was lost)
                {
                    continue;
                }
            }
            return products;
        }

        /// <summary>
        /// Adjusts modification indices.
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="modificationDictionary"></param>
        /// <returns></returns>
        internal Dictionary<int, List<Modification>> AdjustModificationIndices(SequenceVariation variant, ProteinWithAppliedVariants protein)
        {
            IDictionary<int, List<Modification>> modificationDictionary = protein.Protein.OneBasedPossibleLocalizedModifications;
            IDictionary<int, List<Modification>> variantModificationDictionary = variant.OneBasedModifications;
            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            int sequenceLengthChange = variant.VariantSequence.Length - variant.OriginalSequence.Length;

            // change modification indices for variant sequence
            if (modificationDictionary != null)
            {
                foreach (KeyValuePair<int, List<Modification>> kv in modificationDictionary)
                {
                    if (variant.OneBasedBeginPosition > kv.Key)
                    {
                        mods.Add(kv.Key, kv.Value);
                    }
                    else if (variant.OneBasedEndPosition < kv.Key && kv.Key + sequenceLengthChange <= BaseSequence.Length)
                    {
                        mods.Add(kv.Key + sequenceLengthChange, kv.Value);
                    }
                    else // sequence variant conflicts with modification site (modification site substitution)
                    {
                        continue;
                    }
                }
            }

            // sequence variant modifications are indexed to the variant sequence 
            //    NOTE: this code assumes variants are added from end to beginning of protein, so that previously added variant mods are adjusted above
            if (variantModificationDictionary != null)
            {
                foreach (var kv in variantModificationDictionary)
                {
                    if (mods.TryGetValue(kv.Key, out var modsAtPos))
                    {
                        modsAtPos.AddRange(kv.Value);
                    }
                    else
                    {
                        mods.Add(kv.Key, kv.Value);
                    }
                }
            }

            return mods;
        }

        /// <summary>
        /// Gets the accession for a protein with applied variations
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="sequenceVariation"></param>
        public static string GetAccession(Protein protein, IEnumerable<SequenceVariation> appliedSequenceVariations)
        {
            return protein.Accession + (appliedSequenceVariations == null ? "" : "_" + CombineSimpleStrings(appliedSequenceVariations));
        }

        /// <summary>
        /// Format string to append to accession
        /// </summary>
        /// <param name="variations"></param>
        /// <returns></returns>
        internal static string CombineSimpleStrings(IEnumerable<SequenceVariation> variations)
        {
            return variations == null ? "" : string.Join("_", variations.Select(v => v.SimpleString()));
        }

        /// <summary>
        /// Format string to append to protein names
        /// </summary>
        /// <param name="variations"></param>
        /// <returns></returns>
        internal static string CombineDescriptions(IEnumerable<SequenceVariation> variations)
        {
            return variations == null ? "" : string.Join(", variant:", variations.Select(d => d.Description));
        }
    }
}