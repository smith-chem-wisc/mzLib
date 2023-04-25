using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public static class VariantApplication
    {
        /// <summary>
        /// Gets the accession for a protein with applied variations
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="sequenceVariation"></param>
        public static string GetAccession(Protein protein,IEnumerable<SequenceVariation> appliedSequenceVariations)
        {
            return protein.NonVariantProtein.Accession + 
                (appliedSequenceVariations == null || appliedSequenceVariations.Count() == 0 ? "" : $"_{CombineSimpleStrings(appliedSequenceVariations)}");
        }

        /// <summary>
        /// Determines if the modification falls on a variant amino acid
        /// </summary>
        /// <param name="protein"></param>
        /// <param name=""></param>
        /// <returns>true if a modification index on the protein falls within the applied variant</returns>
        /// <remarks>
        /// A. Cesnik - 4/25/23 
        /// Variants annotated in protein entries can be applied to a sequence, i.e. a change is made to the sequence. 
        /// One of the things Spritz can do that no other tool can do is enable finding modifications on these sites of variation, 
        /// since I amended the sequence variant XML entries to have modifications.
        /// </remarks>
        public static bool IsSequenceVariantModification(SequenceVariation appliedVariant, int variantProteinIndex)
        {
            return appliedVariant != null && appliedVariant.Includes(variantProteinIndex);
        }

        /// <summary>
        /// Restores modification index on a variant protein to the index on the nonvariant protein,
        /// or if it falls on a variant, this restores the position on the protein with only that variant
        /// </summary>
        /// <param name="protein">Protein containing applied sequence variations</param>
        /// <param name="variantProteinModificationIndex">The one-based index of the amino acid residue bearing a modification</param>
        /// <returns></returns>
        /// <remarks>
        /// A. Cesnik - 4/25/23 
        /// Useful for comparing modification indices on variant proteins to the original protein.
        /// Variations can introduce length changes and other changes to the sequence, 
        /// so the indices of the modifications aren’t directly comparable, but this method makes that possible.
        /// </remarks>
        public static int RestoreModificationIndex(Protein protein, int variantProteinModificationIndex)
        {
            return variantProteinModificationIndex - protein.AppliedSequenceVariations
                .Where(v => v.OneBasedEndPosition < variantProteinModificationIndex)
                .Sum(v => v.VariantSequence.Length - v.OriginalSequence.Length);
        }

        /// <summary>
        /// Format string to append to accession
        /// </summary>
        /// <param name="variations"></param>
        /// <returns></returns>
        internal static string CombineSimpleStrings(IEnumerable<SequenceVariation> variations)
        {            
            return variations == null || variations.Count() == 0? "" : string.Join("_", variations.Select(v => v.SimpleString()));
        }

        /// <summary>
        /// Format string to append to protein names
        /// </summary>
        /// <param name="variations"></param>
        /// <returns></returns>
        internal static string CombineDescriptions(IEnumerable<SequenceVariation> variations)
        {
            return variations == null || variations.Count() == 0 ? "" : string.Join(", variant:", variations.Select(d => d.Description));
        }

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
                    bool variantAlleleIsInTheGenotype = variant.Description.Genotypes[individual].Contains(variant.Description.AlleleIndex.ToString()); // should catch the case where it's -1 if the INFO isn't from SnpEff
                    if (!variantAlleleIsInTheGenotype)
                    {
                        continue;
                    }
                    bool isHomozygousAlternate = variant.Description.Homozygous[individual] && variant.Description.Genotypes[individual].All(d => d == variant.Description.AlleleIndex.ToString()); // note this isn't a great test for homozygosity, since the genotype could be 1/2 and this would still return true. But currently, alleles 1 and 2 will be included as separate variants, so this is fine for now.
                    bool isDeepReferenceAllele = int.TryParse(variant.Description.AlleleDepths[individual][0], out int depthRef) && depthRef >= minAlleleDepth;
                    bool isDeepAlternateAllele = int.TryParse(variant.Description.AlleleDepths[individual][variant.Description.AlleleIndex], out int depthAlt) && depthAlt >= minAlleleDepth;

                    // homozygous alternate
                    if (isHomozygousAlternate && isDeepAlternateAllele)
                    {
                        newVariantProteins = newVariantProteins.Select(p => ApplySingleVariant(variant, p, individual)).ToList();
                    }

                    // heterozygous basic
                    // first protein with variants contains all homozygous variation, second contains all variations
                    else if (variant.Description.Heterozygous[individual] && tooManyHeterozygousVariants)
                    {
                        if (isDeepAlternateAllele && isDeepReferenceAllele)
                        {
                            if (newVariantProteins.Count == 1 && maxAllowedVariantsForCombinitorics > 0)
                            {
                                Protein variantProtein = ApplySingleVariant(variant, newVariantProteins[0], individual);
                                newVariantProteins.Add(variantProtein);
                            }
                            else if (maxAllowedVariantsForCombinitorics > 0)
                            {
                                newVariantProteins[1] = ApplySingleVariant(variant, newVariantProteins[1], individual);
                            }
                            else
                            {
                                // no heterozygous variants
                            }
                        }
                        else if (isDeepAlternateAllele && maxAllowedVariantsForCombinitorics > 0)
                        {
                            newVariantProteins = newVariantProteins.Select(p => ApplySingleVariant(variant, p, individual)).ToList();
                        }
                        else
                        {
                            // keep reference only
                        }
                    }

                    // heterozygous combinitorics
                    else if (variant.Description.Heterozygous[individual] && isDeepAlternateAllele && !tooManyHeterozygousVariants)
                    {
                        List<Protein> combinitoricProteins = new List<Protein>();

                        foreach (Protein ppp in newVariantProteins)
                        {
                            if (isDeepAlternateAllele && maxAllowedVariantsForCombinitorics > 0 && isDeepReferenceAllele)
                            {
                                // keep reference allele
                                if (variant.Description.Genotypes[individual].Contains("0"))
                                {
                                    combinitoricProteins.Add(ppp);
                                }

                                // alternate allele (replace all, since in heterozygous with two alternates, both alternates are included)
                                combinitoricProteins.Add(ApplySingleVariant(variant, ppp, individual));
                            }
                            else if (isDeepAlternateAllele && maxAllowedVariantsForCombinitorics > 0)
                            {
                                combinitoricProteins.Add(ApplySingleVariant(variant, ppp, individual));
                            }
                            else if (variant.Description.Genotypes[individual].Contains("0"))
                            {
                                combinitoricProteins.Add(ppp);
                            }
                            else
                            {
                                // must be two alternate alleles with not enough depth
                            }
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
        /// <param name="variantGettingApplied"></param>
        /// <returns></returns>
        internal static Protein ApplySingleVariant(SequenceVariation variantGettingApplied, Protein protein, string individual)
        {
            string seqBefore = protein.BaseSequence.Substring(0, variantGettingApplied.OneBasedBeginPosition - 1);
            string seqVariant = variantGettingApplied.VariantSequence;
            int afterIdx = variantGettingApplied.OneBasedBeginPosition + variantGettingApplied.OriginalSequence.Length - 1;

            SequenceVariation variantAfterApplication = new SequenceVariation(
                variantGettingApplied.OneBasedBeginPosition, 
                variantGettingApplied.OneBasedBeginPosition + variantGettingApplied.VariantSequence.Length - 1, 
                variantGettingApplied.OriginalSequence, 
                variantGettingApplied.VariantSequence, 
                variantGettingApplied.Description.Description, 
                variantGettingApplied.OneBasedModifications.ToDictionary(kv => kv.Key, kv => kv.Value));

            // check to see if there is incomplete indel overlap, which would lead to weird variant sequences
            // complete overlap is okay, since it will be overwritten; this can happen if there are two alternate alleles,
            //    e.g. reference sequence is wrong at that point
            bool intersectsAppliedRegionIncompletely = protein.AppliedSequenceVariations.Any(x => variantGettingApplied.Intersects(x) && !variantGettingApplied.Includes(x));
            IEnumerable<SequenceVariation> appliedVariations = new[] { variantAfterApplication };
            string seqAfter = null;
            if (intersectsAppliedRegionIncompletely)
            {
                // use original protein sequence for the remaining sequence
                seqAfter = protein.BaseSequence.Length - afterIdx <= 0 ? "" : protein.NonVariantProtein.BaseSequence.Substring(afterIdx);
            }
            else
            {
                // use this variant protein sequence for the remaining sequence
                seqAfter = protein.BaseSequence.Length - afterIdx <= 0 ? "" : protein.BaseSequence.Substring(afterIdx);
                appliedVariations = appliedVariations
                    .Concat(protein.AppliedSequenceVariations.Where(x => !variantGettingApplied.Includes(x)))
                    .ToList();
            }
            string variantSequence = (seqBefore + seqVariant + seqAfter).Split('*')[0]; // there may be a stop gained

            // adjust indices
            List<ProteolysisProduct> adjustedProteolysisProducts = AdjustProteolysisProductIndices(variantGettingApplied, variantSequence, protein, protein.ProteolysisProducts);
            Dictionary<int, List<Modification>> adjustedModifications = AdjustModificationIndices(variantGettingApplied, variantSequence, protein);
            List<SequenceVariation> adjustedAppliedVariations = AdjustSequenceVariationIndices(variantGettingApplied, variantSequence, appliedVariations);

            return new Protein(variantSequence, protein, adjustedAppliedVariations, adjustedProteolysisProducts, adjustedModifications, individual);
        }

        /// <summary>
        /// Adjusts the indices of sequence variations due to applying a single additional variant
        /// </summary>
        /// <param name="variantGettingApplied"></param>
        /// <param name="alreadyAppliedVariations"></param>
        /// <returns></returns>
        internal static List<SequenceVariation> AdjustSequenceVariationIndices(SequenceVariation variantGettingApplied, string variantAppliedProteinSequence, IEnumerable<SequenceVariation> alreadyAppliedVariations)
        {
            List<SequenceVariation> variations = new List<SequenceVariation>();
            if (alreadyAppliedVariations == null) { return variations; }
            foreach (SequenceVariation v in alreadyAppliedVariations)
            {
                int addedIdx = alreadyAppliedVariations
                    .Where(applied => applied.OneBasedEndPosition < v.OneBasedBeginPosition)
                    .Sum(applied => applied.VariantSequence.Length - applied.OriginalSequence.Length);

                // variant was entirely before the one being applied (shouldn't happen because of order of applying variants)
                // or it's the current variation
                if (v.Description.Equals(variantGettingApplied.Description) || v.OneBasedEndPosition - addedIdx < variantGettingApplied.OneBasedBeginPosition)
                {
                    variations.Add(v);
                }

                // adjust indices based on new included sequence, minding possible overlaps to be filtered later
                else
                {
                    int intersectOneBasedStart = Math.Max(variantGettingApplied.OneBasedBeginPosition, v.OneBasedBeginPosition);
                    int intersectOneBasedEnd = Math.Min(variantGettingApplied.OneBasedEndPosition, v.OneBasedEndPosition);
                    int overlap = intersectOneBasedEnd < intersectOneBasedStart ? 0 : // no overlap
                        intersectOneBasedEnd - intersectOneBasedStart + 1; // there's some overlap
                    int sequenceLengthChange = variantGettingApplied.VariantSequence.Length - variantGettingApplied.OriginalSequence.Length;
                    int begin = v.OneBasedBeginPosition + sequenceLengthChange - overlap;
                    if (begin > variantAppliedProteinSequence.Length)
                    {
                        continue; // cut out by a stop gain
                    }
                    int end = v.OneBasedEndPosition + sequenceLengthChange - overlap;
                    if (end > variantAppliedProteinSequence.Length)
                    {
                        end = variantAppliedProteinSequence.Length; // end shortened by a stop gain
                    }
                    variations.Add(new SequenceVariation(
                        begin,
                        end,
                        v.OriginalSequence,
                        v.VariantSequence,
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
        internal static List<ProteolysisProduct> AdjustProteolysisProductIndices(SequenceVariation variant, string variantAppliedProteinSequence, Protein protein, IEnumerable<ProteolysisProduct> proteolysisProducts)
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
                    && (p.OneBasedEndPosition > variant.OneBasedEndPosition || p.OneBasedEndPosition == protein.NonVariantProtein.BaseSequence.Length))
                {
                    if (variant.VariantSequence.EndsWith("*"))
                    {
                        products.Add(new ProteolysisProduct(p.OneBasedBeginPosition, variantAppliedProteinSequence.Length, p.Type));
                    }
                    else if (p.OneBasedEndPosition + sequenceLengthChange <= variantAppliedProteinSequence.Length)
                    {
                        products.Add(new ProteolysisProduct(p.OneBasedBeginPosition, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
                    }
                    else
                    {
                        // cleavage site is not intact
                    }
                }
                // proteolysis product is after the variant and there is no stop gain
                else if (p.OneBasedBeginPosition > variant.OneBasedEndPosition
                    && p.OneBasedBeginPosition + sequenceLengthChange <= variantAppliedProteinSequence.Length
                    && p.OneBasedEndPosition + sequenceLengthChange <= variantAppliedProteinSequence.Length
                    && !variant.VariantSequence.EndsWith("*"))
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
        internal static Dictionary<int, List<Modification>> AdjustModificationIndices(SequenceVariation variant, string variantAppliedProteinSequence, Protein protein)
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
                    if (kv.Key > variantAppliedProteinSequence.Length)
                    {
                        continue; // it was cut out by a stop gain
                    }
                    // mod is before the variant
                    else if (kv.Key < variant.OneBasedBeginPosition)
                    {
                        mods.Add(kv.Key, kv.Value);
                    }
                    // mod is after the variant and not affected by a stop gain
                    else if (variant.OneBasedEndPosition < kv.Key && kv.Key + sequenceLengthChange <= variantAppliedProteinSequence.Length)
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
    }
}