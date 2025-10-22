using MzLibUtil;
using Omics.BioPolymer;
using Omics.Modifications;

namespace Omics.BioPolymer
{
    /// <summary>
    /// Provides methods for applying sequence variations to proteins and handling modifications on variant sequences.
    /// </summary>
    /// <remarks>
    /// Originally by A. Cesnik on 11/2/18, updated on 4/25/23. NB moved it and generalized for use in Transcriptomics on 3/25/25.
    /// </remarks>
    public static class VariantApplication
    {
        /// <summary>
        /// Creates a list of IBioPolymers of the same type as the original protein, each with applied variants from this protein.
        /// </summary>
        /// <typeparam name="TBioPolymerType">Type of BioPolymer to create variants of</typeparam>
        /// <param name="protein">original to generate variants of</param>
        /// <param name="maxAllowedVariantsForCombinatorics"></param>
        /// <param name="minAlleleDepth"></param>
        /// <remarks>This replaces a method call that was previously an instance method in Protein</remarks>
        public static List<TBioPolymerType> GetVariantBioPolymers<TBioPolymerType>(this TBioPolymerType protein, int maxAllowedVariantsForCombinatorics = 4, int minAlleleDepth = 1)
            where TBioPolymerType : IHasSequenceVariants
        {
            protein.ConsensusVariant.ConvertNucleotideSubstitutionModificationsToSequenceVariants();
            protein.ConvertNucleotideSubstitutionModificationsToSequenceVariants();
            if (protein.SequenceVariations.All(v => v.AreValid()) && protein.SequenceVariations.Any(v => v.Description == null || v.Description.Genotypes.Count == 0))
            {
                // this is a protein with either no VCF lines or a mix of VCF and non-VCF lines
                return ApplyAllVariantCombinations(protein, protein.SequenceVariations, maxAllowedVariantsForCombinatorics).ToList();
            }
            // this is a protein with only VCF lines
            return ApplyVariants(protein, protein.SequenceVariations, maxAllowedVariantsForCombinatorics, minAlleleDepth);
        }

        /// <summary>
        /// Gets the name of a protein with applied variations
        /// </summary>
        public static string? GetVariantName(string? name, IEnumerable<SequenceVariation>? appliedVariations)
        {
            bool emptyVars = appliedVariations.IsNullOrEmpty();
            if (name == null && emptyVars)
                return null;

            string variantTag = emptyVars ? "" : $" variant:{CombineDescriptions(appliedVariations)}";
            return name + variantTag;
        }

        /// <summary>
        /// Gets the accession for a protein with applied variations
        /// </summary>
        public static string GetAccession(IHasSequenceVariants protein, IEnumerable<SequenceVariation>? appliedSequenceVariations)
        {
            return protein.ConsensusVariant.Accession +
                   (appliedSequenceVariations.IsNullOrEmpty() ? "" : $"_{CombineSimpleStrings(appliedSequenceVariations)}");
        }

        /// <summary>
        /// Determines if the modification falls on a variant amino acid
        /// </summary>
        /// <returns>true if a modification index on the protein falls within the applied variant</returns>
        /// <remarks>
        /// A. Cesnik - 4/25/23 
        /// Variants annotated in protein entries can be applied to a sequence, i.e. a change is made to the sequence. 
        /// One of the things Spritz can do that no other tool can do is enable finding modifications on these sites of variation, 
        /// since I amended the sequence variant XML entries to have modifications.
        /// </remarks>
        public static bool IsSequenceVariantModification(SequenceVariation? appliedVariant, int variantProteinIndex)
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
        public static int RestoreModificationIndex(IHasSequenceVariants protein, int variantProteinModificationIndex)
        {
            return variantProteinModificationIndex - protein.AppliedSequenceVariations
                .Where(v => v.OneBasedEndPosition < variantProteinModificationIndex)
                .Sum(v => v.VariantSequence.Length - v.OriginalSequence.Length);
        }

        /// <summary>
        /// Applies multiple variant changes to a protein sequence
        /// </summary>
        public static List<TBioPolymerType> ApplyVariants<TBioPolymerType>(TBioPolymerType protein, IEnumerable<SequenceVariation> sequenceVariations, int maxAllowedVariantsForCombinitorics, int minAlleleDepth)
            where TBioPolymerType : IHasSequenceVariants
        {
            List<SequenceVariation> uniqueEffectsToApply = sequenceVariations
                .GroupBy(v => v.SimpleString())
                .Select(x => x.First())
                .Where(v => v.Description.Genotypes.Count > 0) // this is a VCF line
                .OrderByDescending(v => v.OneBasedBeginPosition) // apply variants at the end of the protein sequence first
                .ToList();

            TBioPolymerType proteinCopy = protein.CreateVariant(protein.BaseSequence, protein, null, protein.TruncationProducts, protein.OneBasedPossibleLocalizedModifications, null);

            // If there aren't any variants to apply, just return the base protein
            if (uniqueEffectsToApply.Count == 0)
            {
                return new List<TBioPolymerType> { proteinCopy };
            }

            HashSet<string> individuals = new HashSet<string>(uniqueEffectsToApply.SelectMany(v => v.Description.Genotypes.Keys));
            List<TBioPolymerType> variantProteins = new();
            List<TBioPolymerType> newVariantProteins = new();
            // loop through genotypes for each sample/individual (e.g. tumor and normal)
            foreach (string individual in individuals)
            {
                newVariantProteins.Clear();
                newVariantProteins.Add(proteinCopy);

                bool tooManyHeterozygousVariants = uniqueEffectsToApply.Count(v => v.Description.Heterozygous[individual]) > maxAllowedVariantsForCombinitorics;
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
                                TBioPolymerType variantProtein = ApplySingleVariant(variant, newVariantProteins[0], individual);
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
                        List<TBioPolymerType> combinitoricProteins = new();

                        foreach (var ppp in newVariantProteins)
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
        private static TBioPolymerType ApplySingleVariant<TBioPolymerType>(SequenceVariation variantGettingApplied, TBioPolymerType protein, string individual)
            where TBioPolymerType : IHasSequenceVariants
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
                seqAfter = protein.BaseSequence.Length - afterIdx <= 0 ? "" : protein.ConsensusVariant.BaseSequence.Substring(afterIdx);
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
            List<TruncationProduct> adjustedProteolysisProducts = AdjustTruncationProductIndices(variantGettingApplied, variantSequence, protein, protein.TruncationProducts);
            Dictionary<int, List<Modification>> adjustedModifications = AdjustModificationIndices(variantGettingApplied, variantSequence, protein);
            List<SequenceVariation> adjustedAppliedVariations = AdjustSequenceVariationIndices(variantGettingApplied, variantSequence, appliedVariations);

            return protein.CreateVariant(variantSequence, protein, adjustedAppliedVariations, adjustedProteolysisProducts, adjustedModifications, individual);
        }

        /// <summary>
        /// Adjusts the indices of sequence variations due to applying a single additional variant
        /// </summary>
        private static List<SequenceVariation> AdjustSequenceVariationIndices(SequenceVariation variantGettingApplied, string variantAppliedProteinSequence, IEnumerable<SequenceVariation> alreadyAppliedVariations)
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
        private static List<TruncationProduct> AdjustTruncationProductIndices(SequenceVariation variant, string variantAppliedProteinSequence, IHasSequenceVariants protein, IEnumerable<TruncationProduct> proteolysisProducts)
        {
            List<TruncationProduct> products = new List<TruncationProduct>();
            if (proteolysisProducts == null) { return products; }
            int sequenceLengthChange = variant.VariantSequence.Length - variant.OriginalSequence.Length;
            foreach (TruncationProduct p in proteolysisProducts.Where(p => p.OneBasedEndPosition.HasValue && p.OneBasedBeginPosition.HasValue))
            {
                // proteolysis product is entirely before the variant
                if (variant.OneBasedBeginPosition > p.OneBasedEndPosition)
                {
                    products.Add(p);
                }
                // proteolysis product straddles the variant, but the cleavage site(s) are still intact; the ends aren't considered cleavage sites
                else if ((p.OneBasedBeginPosition < variant.OneBasedBeginPosition || p.OneBasedBeginPosition == 1 || p.OneBasedBeginPosition == 2)
                         && (p.OneBasedEndPosition > variant.OneBasedEndPosition || p.OneBasedEndPosition == protein.ConsensusVariant.BaseSequence.Length))
                {
                    if (variant.VariantSequence.EndsWith("*"))
                    {
                        products.Add(new TruncationProduct(p.OneBasedBeginPosition, variantAppliedProteinSequence.Length, p.Type));
                    }
                    else if (p.OneBasedEndPosition + sequenceLengthChange <= variantAppliedProteinSequence.Length)
                    {
                        products.Add(new TruncationProduct(p.OneBasedBeginPosition, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
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
                    products.Add(new TruncationProduct(p.OneBasedBeginPosition + sequenceLengthChange, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
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
        private static Dictionary<int, List<Modification>> AdjustModificationIndices(SequenceVariation variant, string variantAppliedProteinSequence, IHasSequenceVariants protein)
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

        /// <summary>
        /// Format string to append to accession
        /// </summary>
        private static string CombineSimpleStrings(IEnumerable<SequenceVariation>? variations)
        {
            return variations.IsNullOrEmpty() ? "" : string.Join("_", variations.Select(v => v.SimpleString()));
        }

        /// <summary>
        /// Format string to append to protein names
        /// </summary>
        public static string CombineDescriptions(IEnumerable<SequenceVariation>? variations)
        {
            return variations.IsNullOrEmpty() ? "" : string.Join(", variant:", variations.Select(d => d.Description));
        }
        /// <summary>
        /// Applies all possible combinations of the provided SequenceVariation list to the base TBioPolymerType object,
        /// starting with the fewest single variations and up to the specified maximum number of combinations.
        /// </summary>
        /// <typeparam name="TBioPolymerType">The type of the biopolymer object.</typeparam>
        /// <param name="baseBioPolymer">The base biopolymer object to apply variations to.</param>
        /// <param name="variations">List of SequenceVariation objects to combine and apply. Assumed not null or empty.</param>
        /// <param name="maxCombinations">Maximum number of combinations to return.</param>
        /// <returns>
        /// An IEnumerable of TBioPolymerType objects, each with a unique combination of variations applied.
        /// </returns>
        public static IEnumerable<TBioPolymerType> ApplyAllVariantCombinations<TBioPolymerType>(
            TBioPolymerType baseBioPolymer,
            List<SequenceVariation> variations,
            int maxCombinations)
            where TBioPolymerType : IHasSequenceVariants
        {
            int count = 0;

            // Always yield the base biopolymer first
            yield return baseBioPolymer;
            count++;
            if (count >= maxCombinations)
                yield break;

            int n = variations.Count;
            for (int size = 1; size <= n; size++)
            {
                foreach (var combo in GetCombinations(variations, size))
                {
                    var result = baseBioPolymer;
                    foreach (var variant in combo)
                    {
                        result = ApplySingleVariant(variant, result, string.Empty);
                    }
                    if (result != null)
                    {
                        yield return result;
                        count++;
                        if (count >= maxCombinations)
                            yield break;
                    }
                }
            }
        }

        /// <summary>
        /// Generates all possible combinations of the specified size from the input list.
        /// </summary>
        /// <param name="variations">List of SequenceVariation objects to combine. Assumed not null or empty.</param>
        /// <param name="size">The size of each combination.</param>
        /// <returns>
        /// An IEnumerable of IList&lt;SequenceVariation&gt; representing each combination.
        /// </returns>
        private static IEnumerable<IList<SequenceVariation>> GetCombinations(List<SequenceVariation> variations, int size)
        {
            int n = variations.Count;
            var indices = new int[size];
            for (int i = 0; i < size; i++) indices[i] = i;

            while (true)
            {
                var combo = new List<SequenceVariation>(size);
                for (int i = 0; i < size; i++)
                    combo.Add(variations[indices[i]]);
                yield return combo;

                int pos = size - 1;
                while (pos >= 0 && indices[pos] == n - size + pos)
                    pos--;
                if (pos < 0) break;
                indices[pos]++;
                for (int i = pos + 1; i < size; i++)
                    indices[i] = indices[i - 1] + 1;
            }
        }
        public static void ConvertNucleotideSubstitutionModificationsToSequenceVariants<TBioPolymerType>(this TBioPolymerType protein)
            where TBioPolymerType : IHasSequenceVariants
        {
            List<KeyValuePair<int, Modification>> modificationsToRemove = new();
            //convert substitution modifications to sequence variations
            foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
            {
                foreach (Modification mod in kvp.Value)
                {
                    if (mod.ModificationType.Contains("nucleotide substitution") && mod.OriginalId.Contains("->"))
                    {
                        string[] originalAndSubstitutedAminoAcids = mod.OriginalId.Split(new[] { "->" }, StringSplitOptions.RemoveEmptyEntries);
                        SequenceVariation sequenceVariation = new SequenceVariation(kvp.Key, kvp.Key, originalAndSubstitutedAminoAcids[0], originalAndSubstitutedAminoAcids[1], "Putative GPTMD Substitution");
                        if (!protein.SequenceVariations.Contains(sequenceVariation))
                        {
                            protein.SequenceVariations.Add(sequenceVariation);
                        }
                        KeyValuePair<int, Modification> pair = new(kvp.Key, mod);
                        modificationsToRemove.Add(pair);
                    }
                }
            }
            //remove the modifications that were converted to sequence variations
            foreach (KeyValuePair<int, Modification> pair in modificationsToRemove)
            {
                if (protein.OneBasedPossibleLocalizedModifications.ContainsKey(pair.Key))
                {
                    List<Modification> modList = protein.OneBasedPossibleLocalizedModifications[pair.Key];
                    var modToRemove = modList.FirstOrDefault(m =>
                        m.ModificationType == pair.Value.ModificationType &&
                        m.OriginalId == pair.Value.OriginalId);
                    if (modToRemove != null)
                    {
                        modList.Remove(modToRemove);
                        if (modList.Count == 0)
                        {
                            protein.OneBasedPossibleLocalizedModifications.Remove(pair.Key);
                            protein.ConsensusVariant.OneBasedPossibleLocalizedModifications.Remove(pair.Key);
                        }
                    }
                }
                if (protein.OriginalNonVariantModifications.ContainsKey(pair.Key))
                {
                    List<Modification> modList = protein.OriginalNonVariantModifications[pair.Key];
                    var modToRemove = modList.FirstOrDefault(m =>
                        m.ModificationType == pair.Value.ModificationType &&
                        m.OriginalId == pair.Value.OriginalId);
                    if (modToRemove != null)
                    {
                        modList.Remove(modToRemove);
                        if (modList.Count == 0)
                        {
                            protein.OriginalNonVariantModifications.Remove(pair.Key);
                            protein.ConsensusVariant.OriginalNonVariantModifications.Remove(pair.Key);
                        }
                    }
                }
            }
        }
    }
}