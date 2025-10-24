using MzLibUtil;
using Omics.BioPolymer;
using Omics.Modifications;

namespace Omics.BioPolymer
{
    public static class VariantApplication
    {
        public static List<TBioPolymerType> GetVariantBioPolymers<TBioPolymerType>(this TBioPolymerType protein, int maxAllowedVariantsForCombinatorics = 4, int minAlleleDepth = 1)
            where TBioPolymerType : IHasSequenceVariants
        {
            protein.ConsensusVariant.ConvertNucleotideSubstitutionModificationsToSequenceVariants();
            protein.ConvertNucleotideSubstitutionModificationsToSequenceVariants();
            if (protein.SequenceVariations.All(v => v.AreValid()) && protein.SequenceVariations.Any(v => v.VariantCallFormatDataString == null || v.VariantCallFormatDataString.Genotypes.Count == 0))
            {
                return ApplyAllVariantCombinations(protein, protein.SequenceVariations, maxAllowedVariantsForCombinatorics).ToList();
            }
            return ApplyVariants(protein, protein.SequenceVariations, maxAllowedVariantsForCombinitorics: maxAllowedVariantsForCombinatorics, minAlleleDepth);
        }

        public static string? GetVariantName(string? name, IEnumerable<SequenceVariation>? appliedVariations)
        {
            bool emptyVars = appliedVariations.IsNullOrEmpty();
            if (name == null && emptyVars)
                return null;

            string variantTag = emptyVars ? "" : $" variant:{CombineDescriptions(appliedVariations)}";
            return name + variantTag;
        }

        public static string GetAccession(IHasSequenceVariants protein, IEnumerable<SequenceVariation>? appliedSequenceVariations)
        {
            return protein.ConsensusVariant.Accession +
                   (appliedSequenceVariations.IsNullOrEmpty() ? "" : $"_{CombineSimpleStrings(appliedSequenceVariations)}");
        }

        public static bool IsSequenceVariantModification(SequenceVariation? appliedVariant, int variantProteinIndex)
        {
            return appliedVariant != null && appliedVariant.Includes(variantProteinIndex);
        }

        public static int RestoreModificationIndex(IHasSequenceVariants protein, int variantProteinModificationIndex)
        {
            return variantProteinModificationIndex - protein.AppliedSequenceVariations
                .Where(v => v.OneBasedEndPosition < variantProteinModificationIndex)
                .Sum(v => v.VariantSequence.Length - v.OriginalSequence.Length);
        }

        public static List<TBioPolymerType> ApplyVariants<TBioPolymerType>(TBioPolymerType protein, IEnumerable<SequenceVariation> sequenceVariations, int maxAllowedVariantsForCombinitorics, int minAlleleDepth)
            where TBioPolymerType : IHasSequenceVariants
        {
            List<SequenceVariation> uniqueEffectsToApply = sequenceVariations
                .GroupBy(v => v.SimpleString())
                .Select(x => x.First())
                .Where(v => v.VariantCallFormatDataString.Genotypes.Count > 0)
                .OrderByDescending(v => v.OneBasedBeginPosition)
                .ToList();

            TBioPolymerType proteinCopy = protein.CreateVariant(protein.BaseSequence, protein, null, protein.TruncationProducts, protein.OneBasedPossibleLocalizedModifications, null);

            if (uniqueEffectsToApply.Count == 0)
            {
                return new List<TBioPolymerType> { proteinCopy };
            }

            HashSet<string> individuals = new HashSet<string>(uniqueEffectsToApply.SelectMany(v => v.VariantCallFormatDataString.Genotypes.Keys));
            List<TBioPolymerType> variantProteins = new();
            List<TBioPolymerType> newVariantProteins = new();
            foreach (string individual in individuals)
            {
                newVariantProteins.Clear();
                newVariantProteins.Add(proteinCopy);

                bool tooManyHeterozygousVariants = uniqueEffectsToApply.Count(v => v.VariantCallFormatDataString.Heterozygous[individual]) > maxAllowedVariantsForCombinitorics;
                foreach (var variant in uniqueEffectsToApply)
                {
                    bool variantAlleleIsInTheGenotype = variant.VariantCallFormatDataString.Genotypes[individual].Contains(variant.VariantCallFormatDataString.AlleleIndex.ToString());
                    if (!variantAlleleIsInTheGenotype)
                    {
                        continue;
                    }
                    bool isHomozygousAlternate = variant.VariantCallFormatDataString.Homozygous[individual] && variant.VariantCallFormatDataString.Genotypes[individual].All(d => d == variant.VariantCallFormatDataString.AlleleIndex.ToString());
                    bool isDeepReferenceAllele = int.TryParse(variant.VariantCallFormatDataString.AlleleDepths[individual][0], out int depthRef) && depthRef >= minAlleleDepth;
                    bool isDeepAlternateAllele = int.TryParse(variant.VariantCallFormatDataString.AlleleDepths[individual][variant.VariantCallFormatDataString.AlleleIndex], out int depthAlt) && depthAlt >= minAlleleDepth;

                    if (isHomozygousAlternate && isDeepAlternateAllele)
                    {
                        newVariantProteins = newVariantProteins.Select(p => ApplySingleVariant(variant, p, individual)).ToList();
                    }
                    else if (variant.VariantCallFormatDataString.Heterozygous[individual] && tooManyHeterozygousVariants)
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
                        }
                        else if (isDeepAlternateAllele && maxAllowedVariantsForCombinitorics > 0)
                        {
                            newVariantProteins = newVariantProteins.Select(p => ApplySingleVariant(variant, p, individual)).ToList();
                        }
                    }
                    else if (variant.VariantCallFormatDataString.Heterozygous[individual] && isDeepAlternateAllele && !tooManyHeterozygousVariants)
                    {
                        List<TBioPolymerType> combinitoricProteins = new();

                        foreach (var ppp in newVariantProteins)
                        {
                            if (isDeepAlternateAllele && maxAllowedVariantsForCombinitorics > 0 && isDeepReferenceAllele)
                            {
                                if (variant.VariantCallFormatDataString.Genotypes[individual].Contains("0"))
                                {
                                    combinitoricProteins.Add(ppp);
                                }
                                combinitoricProteins.Add(ApplySingleVariant(variant, ppp, individual));
                            }
                            else if (isDeepAlternateAllele && maxAllowedVariantsForCombinitorics > 0)
                            {
                                combinitoricProteins.Add(ApplySingleVariant(variant, ppp, individual));
                            }
                            else if (variant.VariantCallFormatDataString.Genotypes[individual].Contains("0"))
                            {
                                combinitoricProteins.Add(ppp);
                            }
                        }
                        newVariantProteins = combinitoricProteins;
                    }
                }
                variantProteins.AddRange(newVariantProteins);
            }

            return variantProteins.GroupBy(x => x.BaseSequence).Select(x => x.First()).ToList();
        }

        private static TBioPolymerType ApplySingleVariant<TBioPolymerType>(SequenceVariation variantGettingApplied, TBioPolymerType protein, string individual)
            where TBioPolymerType : IHasSequenceVariants
        {
            string seqBefore = protein.BaseSequence.Substring(0, variantGettingApplied.OneBasedBeginPosition - 1);
            string seqVariant = variantGettingApplied.VariantSequence;
            int afterIdx = variantGettingApplied.OneBasedBeginPosition + variantGettingApplied.OriginalSequence.Length - 1;

            SequenceVariation variantAfterApplication;
            var vcf = variantGettingApplied.VariantCallFormatDataString;
            if (vcf != null)
            {
                variantAfterApplication = new SequenceVariation(
                    variantGettingApplied.OneBasedBeginPosition,
                    variantGettingApplied.OneBasedBeginPosition + variantGettingApplied.VariantSequence.Length - 1,
                    variantGettingApplied.OriginalSequence,
                    variantGettingApplied.VariantSequence,
                    vcf.Description,
                    vcf,
                    variantGettingApplied.OneBasedModifications.ToDictionary(kv => kv.Key, kv => kv.Value));
            }
            else
            {
                variantAfterApplication = new SequenceVariation(
                    variantGettingApplied.OneBasedBeginPosition,
                    variantGettingApplied.OneBasedBeginPosition + variantGettingApplied.VariantSequence.Length - 1,
                    variantGettingApplied.OriginalSequence,
                    variantGettingApplied.VariantSequence,
                    variantGettingApplied.Description,
                    variantGettingApplied.OneBasedModifications.ToDictionary(kv => kv.Key, kv => kv.Value));
            }

            // NULL-SAFE: AppliedSequenceVariations can be null for base proteins
            bool intersectsAppliedRegionIncompletely =
                protein.AppliedSequenceVariations != null
                && protein.AppliedSequenceVariations.Any(x => variantGettingApplied.Intersects(x) && !variantGettingApplied.Includes(x));

            IEnumerable<SequenceVariation> appliedVariations = new[] { variantAfterApplication };
            string seqAfter;
            if (intersectsAppliedRegionIncompletely)
            {
                seqAfter = protein.BaseSequence.Length - afterIdx <= 0 ? "" : protein.ConsensusVariant.BaseSequence.Substring(afterIdx);
            }
            else
            {
                seqAfter = protein.BaseSequence.Length - afterIdx <= 0 ? "" : protein.BaseSequence.Substring(afterIdx);
                appliedVariations = appliedVariations
                    .Concat((protein.AppliedSequenceVariations ?? Enumerable.Empty<SequenceVariation>()).Where(x => !variantGettingApplied.Includes(x)))
                    .ToList();
            }
            string variantSequence = (seqBefore + seqVariant + seqAfter).Split('*')[0];

            List<TruncationProduct> adjustedProteolysisProducts = AdjustTruncationProductIndices(variantGettingApplied, variantSequence, protein, protein.TruncationProducts);
            Dictionary<int, List<Modification>> adjustedModifications = AdjustModificationIndices(variantGettingApplied, variantSequence, protein);
            List<SequenceVariation> adjustedAppliedVariations = AdjustSequenceVariationIndices(variantGettingApplied, variantSequence, appliedVariations);

            return protein.CreateVariant(variantSequence, protein, adjustedAppliedVariations, adjustedProteolysisProducts, adjustedModifications, individual);
        }

        private static List<SequenceVariation> AdjustSequenceVariationIndices(SequenceVariation variantGettingApplied, string variantAppliedProteinSequence, IEnumerable<SequenceVariation> alreadyAppliedVariations)
        {
            List<SequenceVariation> variations = new List<SequenceVariation>();
            if (alreadyAppliedVariations == null) { return variations; }
            foreach (SequenceVariation v in alreadyAppliedVariations)
            {
                int addedIdx = alreadyAppliedVariations
                    .Where(applied => applied.OneBasedEndPosition < v.OneBasedBeginPosition)
                    .Sum(applied => applied.VariantSequence.Length - applied.OriginalSequence.Length);

                // NULL-SAFE compare; or it is the current variation
                if ((v.VariantCallFormatDataString != null && v.VariantCallFormatDataString.Equals(variantGettingApplied.VariantCallFormatDataString))
                    || v.OneBasedEndPosition - addedIdx < variantGettingApplied.OneBasedBeginPosition)
                {
                    variations.Add(v);
                }
                else
                {
                    int intersectOneBasedStart = Math.Max(variantGettingApplied.OneBasedBeginPosition, v.OneBasedBeginPosition);
                    int intersectOneBasedEnd = Math.Min(variantGettingApplied.OneBasedEndPosition, v.OneBasedEndPosition);
                    int overlap = intersectOneBasedEnd < intersectOneBasedStart ? 0 :
                        intersectOneBasedEnd - intersectOneBasedStart + 1;
                    int sequenceLengthChange = variantGettingApplied.VariantSequence.Length - variantGettingApplied.OriginalSequence.Length;
                    int begin = v.OneBasedBeginPosition + sequenceLengthChange - overlap;
                    if (begin > variantAppliedProteinSequence.Length)
                    {
                        continue;
                    }
                    int end = v.OneBasedEndPosition + sequenceLengthChange - overlap;
                    if (end > variantAppliedProteinSequence.Length)
                    {
                        end = variantAppliedProteinSequence.Length;
                    }

                    var vcf = v.VariantCallFormatDataString;
                    if (vcf != null)
                    {
                        variations.Add(new SequenceVariation(
                            begin,
                            end,
                            v.OriginalSequence,
                            v.VariantSequence,
                            vcf.Description,
                            vcf,
                            v.OneBasedModifications.ToDictionary(kv => kv.Key, kv => kv.Value)));
                    }
                    else
                    {
                        variations.Add(new SequenceVariation(
                            begin,
                            end,
                            v.OriginalSequence,
                            v.VariantSequence,
                            v.Description,
                            v.OneBasedModifications.ToDictionary(kv => kv.Key, kv => kv.Value)));
                    }
                }
            }
            return variations;
        }

        private static List<TruncationProduct> AdjustTruncationProductIndices(SequenceVariation variant, string variantAppliedProteinSequence, IHasSequenceVariants protein, IEnumerable<TruncationProduct> proteolysisProducts)
        {
            List<TruncationProduct> products = new List<TruncationProduct>();
            if (proteolysisProducts == null) { return products; }
            int sequenceLengthChange = variant.VariantSequence.Length - variant.OriginalSequence.Length;
            foreach (TruncationProduct p in proteolysisProducts.Where(p => p.OneBasedEndPosition.HasValue && p.OneBasedBeginPosition.HasValue))
            {
                if (variant.OneBasedBeginPosition > p.OneBasedEndPosition)
                {
                    products.Add(p);
                }
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
                }
                else if (p.OneBasedBeginPosition > variant.OneBasedEndPosition
                         && p.OneBasedBeginPosition + sequenceLengthChange <= variantAppliedProteinSequence.Length
                         && p.OneBasedEndPosition + sequenceLengthChange <= variantAppliedProteinSequence.Length
                         && !variant.VariantSequence.EndsWith("*"))
                {
                    products.Add(new TruncationProduct(p.OneBasedBeginPosition + sequenceLengthChange, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
                }
            }
            return products;
        }

        private static Dictionary<int, List<Modification>> AdjustModificationIndices(SequenceVariation variant, string variantAppliedProteinSequence, IHasSequenceVariants protein)
        {
            IDictionary<int, List<Modification>> modificationDictionary = protein.OneBasedPossibleLocalizedModifications;
            IDictionary<int, List<Modification>> variantModificationDictionary = variant.OneBasedModifications;
            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            int sequenceLengthChange = variant.VariantSequence.Length - variant.OriginalSequence.Length;

            if (modificationDictionary != null)
            {
                foreach (KeyValuePair<int, List<Modification>> kv in modificationDictionary)
                {
                    if (kv.Key > variantAppliedProteinSequence.Length)
                    {
                        continue;
                    }
                    else if (kv.Key < variant.OneBasedBeginPosition)
                    {
                        mods.Add(kv.Key, kv.Value);
                    }
                    else if (variant.OneBasedEndPosition < kv.Key && kv.Key + sequenceLengthChange <= variantAppliedProteinSequence.Length)
                    {
                        mods.Add(kv.Key + sequenceLengthChange, kv.Value);
                    }
                }
            }

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

        private static string CombineSimpleStrings(IEnumerable<SequenceVariation>? variations)
        {
            return variations.IsNullOrEmpty() ? "" : string.Join("_", variations.Select(v => v.SimpleString()));
        }

        public static string CombineDescriptions(IEnumerable<SequenceVariation>? variations)
        {
            return variations.IsNullOrEmpty() ? "" : string.Join(", variant:", variations.Select(d => d.VariantCallFormatDataString?.Description ?? d.Description));
        }

        public static IEnumerable<TBioPolymerType> ApplyAllVariantCombinations<TBioPolymerType>(
            TBioPolymerType baseBioPolymer,
            List<SequenceVariation> variations,
            int maxCombinations)
            where TBioPolymerType : IHasSequenceVariants
        {
            int count = 0;

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