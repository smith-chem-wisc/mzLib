using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public static class VariantApplication
    {
        /// <summary>
        /// Applies multiple variant changes to a protein sequence
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="uniqueEffectsToApply"></param>
        /// <returns></returns>
        internal static List<Protein> ApplyVariants(Protein protein, IEnumerable<SequenceVariation> sequenceVariations, int maxAllowedVariantsForCombinitorics, int minAlleleDepth)
        {
            List<SequenceVariation> uniqueEffectsToApply = sequenceVariations
                .GroupBy(v => v.SimpleString())
                .Select(x => x.First())
                .Where(v => v.Description.Genotypes.Count > 0) // this is a VCF line
                .OrderByDescending(v => v.OneBasedBeginPosition) // apply variants at the end of the protein sequence first
                .ToList();

            Protein proteinCopy = new Protein(protein.BaseSequence, protein, null, protein.ProteolysisProducts, protein.OneBasedPossibleLocalizedModifications, null);

            // If there aren't any variants to apply, just return the base protein
            if (uniqueEffectsToApply.Count == 0)
            {
                return new List<Protein> { proteinCopy };
            }

            HashSet<string> individuals = new HashSet<string>(uniqueEffectsToApply.SelectMany(v => v.Description.Genotypes.Keys));
            List<Protein> variantProteins = new List<Protein>();

            // loop through genotypes for each sample/individual (e.g. tumor and normal)
            foreach (string individual in individuals)
            {
                bool tooManyHeterozygousVariants = uniqueEffectsToApply.Count(v => v.Description.Heterozygous[individual]) > maxAllowedVariantsForCombinitorics;
                List<Protein> newVariantProteins = new List<Protein> { proteinCopy };
                foreach (var variant in uniqueEffectsToApply)
                {
                    bool isHomozygousAlternate = variant.Description.Homozygous[individual] && variant.Description.Genotypes[individual].All(d => d != "0");
                    bool isDeepAlternateAllele = variant.Description.GenotypeAlleleDepthMap[individual].All(d => d.Item1 != "0" || int.TryParse(d.Item2, out int depth) && depth >= minAlleleDepth);

                    // homozygous alternate
                    if (isHomozygousAlternate && isDeepAlternateAllele)
                    {
                        newVariantProteins = newVariantProteins.Select(p => ApplySingleVariant(variant, p, individual)).ToList();
                    }

                    // heterozygous basic
                    // first protein with variants contains all homozygous variation, second contains all variations
                    else if (variant.Description.Heterozygous[individual] && isDeepAlternateAllele && tooManyHeterozygousVariants)
                    {
                        if (newVariantProteins.Count == 1)
                        {
                            Protein variantProtein = ApplySingleVariant(variant, newVariantProteins[0], individual);
                            newVariantProteins.Add(variantProtein);
                        }
                        else
                        {
                            newVariantProteins[1] = ApplySingleVariant(variant, newVariantProteins[1], individual);
                        }
                    }

                    // heterozygous combinitorics
                    else if (variant.Description.Heterozygous[individual] && isDeepAlternateAllele && !tooManyHeterozygousVariants)
                    {
                        List<Protein> combinitoricProteins = new List<Protein>();

                        foreach (Protein ppp in newVariantProteins)
                        {
                            // keep reference allele
                            if (variant.Description.Genotypes[individual].Contains("0"))
                            {
                                combinitoricProteins.Add(ppp);
                            }

                            // alternate allele (replace all, since in heterozygous with two alternates, both alternates are included)
                            combinitoricProteins.Add(ApplySingleVariant(variant, ppp, individual));
                        }
                        newVariantProteins = combinitoricProteins;
                    }
                }
                variantProteins.AddRange(newVariantProteins);
            }

            return variantProteins.GroupBy(x => x.BaseSequence).Select(x => x.First()).ToList();
        }

        /// <summary>
        /// Applies a single variant to a protein sequence
        /// </summary>
        /// <param name="variant"></param>
        /// <returns></returns>
        internal static Protein ApplySingleVariant(SequenceVariation variant, Protein protein, string individual)
        {
            string seqBefore = protein.BaseSequence.Substring(0, variant.OneBasedBeginPosition - 1);
            string seqVariant = variant.VariantSequence;
            List<ProteolysisProduct> adjustedProteolysisProducts = AdjustProteolysisProductIndices(variant, protein, protein.ProteolysisProducts);
            Dictionary<int, List<Modification>> adjustedModifications = AdjustModificationIndices(variant, protein);
            List<SequenceVariation> adjustedAppliedVariations = AdjustSequenceVariations(variant, protein.AppliedSequenceVariations);
            int afterIdx = variant.OneBasedBeginPosition + variant.OriginalSequence.Length - 1;

            // check to see if there is incomplete indel overlap, which would lead to weird variant sequences
            // complete overlap is okay, since it will be overwritten; this can happen if there are two alternate alleles,
            //    e.g. reference sequence is wrong at that point
            bool intersectsAppliedRegionIncompletely = protein.AppliedSequenceVariations.Any(x => variant.Intersects(x) && !variant.Includes(x));
            if (intersectsAppliedRegionIncompletely)
            {
                // use original protein sequence for the remaining sequence
                string seqAfter = protein.BaseSequence.Length - afterIdx <= 0 ? "" : protein.NonVariantBaseSequence.Substring(afterIdx);
                return new Protein(seqBefore + seqVariant + seqAfter, protein, new[] { variant }, adjustedProteolysisProducts, adjustedModifications, individual);
            }
            else
            {
                List<SequenceVariation> variations = protein.AppliedSequenceVariations
                    .Where(x => !variant.Includes(x))
                    .Concat(new[] { variant })
                    .ToList();
                // use this variant protein sequence for the remaining sequence
                string seqAfter = protein.BaseSequence.Length - afterIdx <= 0 ? "" : protein.BaseSequence.Substring(afterIdx);
                return new Protein(seqBefore + seqVariant + seqAfter, protein, variations, adjustedProteolysisProducts, adjustedModifications, individual);
            }
        }

        /// <summary>
        /// Adjusts the indices of sequence variations due to applying a single additional variant
        /// </summary>
        /// <param name="variantGettingApplied"></param>
        /// <param name="alreadyAppliedVariations"></param>
        /// <returns></returns>
        internal static List<SequenceVariation> AdjustSequenceVariations(SequenceVariation variantGettingApplied, IEnumerable<SequenceVariation> alreadyAppliedVariations)
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
        internal static List<ProteolysisProduct> AdjustProteolysisProductIndices(SequenceVariation variant, Protein protein, IEnumerable<ProteolysisProduct> proteolysisProducts)
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
                    && (p.OneBasedEndPosition > variant.OneBasedEndPosition || p.OneBasedEndPosition == protein.NonVariantBaseSequence.Length)
                    && p.OneBasedEndPosition + sequenceLengthChange <= protein.BaseSequence.Length)
                {
                    products.Add(new ProteolysisProduct(p.OneBasedBeginPosition, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
                }
                // proteolysis product is after the variant
                else if (p.OneBasedBeginPosition > variant.OneBasedEndPosition
                    && p.OneBasedBeginPosition + sequenceLengthChange <= protein.BaseSequence.Length
                    && p.OneBasedEndPosition + sequenceLengthChange <= protein.BaseSequence.Length)
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
        internal static Dictionary<int, List<Modification>> AdjustModificationIndices(SequenceVariation variant, Protein protein)
        {
            IDictionary<int, List<Modification>> modificationDictionary = protein.OneBasedPossibleLocalizedModifications;
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
                    else if (variant.OneBasedEndPosition < kv.Key && kv.Key + sequenceLengthChange <= protein.BaseSequence.Length)
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