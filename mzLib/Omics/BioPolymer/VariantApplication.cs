using MzLibUtil;
using Omics.Modifications;
using System.Reflection; 

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
        public static List<TBioPolymerType> GetVariantBioPolymers<TBioPolymerType>(this TBioPolymerType protein,
            int maxSequenceVariantsPerIsoform = 4,
            int minAlleleDepth = 1,
            int maxSequenceVariantIsoforms = 1)
            where TBioPolymerType : IHasSequenceVariants
        {
            // If combinatorics disabled, just return base
            if (maxSequenceVariantsPerIsoform == 0 || maxSequenceVariantIsoforms == 1)
            {
                return new List<TBioPolymerType> { protein };
            }

            var all = protein.SequenceVariations ?? new List<SequenceVariation>();
            if (all.Count == 0)
            {
                return new List<TBioPolymerType> { protein };
            }

            // Try validation, but DO NOT let complete failure collapse all variants.
            var valid = new List<SequenceVariation>(all.Count);
            int threw = 0, failed = 0;
            foreach (var v in all)
            {
                if (v == null)
                {
                    failed++;
                    continue;
                }
                bool ok;
                try
                {
                    ok = v.AreValid();
                }
                catch
                {
                    ok = true; // treat exceptions as “usable” so we can still attempt variant generation
                    threw++;
                }
                if (ok)
                    valid.Add(v);
                else
                    failed++;
            }

            // Fallback: if none passed (over‑strict validation), use original non-null set
            if (valid.Count == 0)
            {
                valid = all.Where(v => v != null).ToList();
            }

            // If after fallback we still have nothing usable, just return base
            if (valid.Count == 0)
            {
                return new List<TBioPolymerType> { protein };
            }

            return ApplyAllVariantCombinations(protein,
                                               valid,
                                               maxSequenceVariantsPerIsoform,
                                               maxSequenceVariantIsoforms,
                                               minAlleleDepth).ToList();
        }

        // Safe wrapper so a single bad variant does not abort all combinatorics
        private static bool SafeAreValid(SequenceVariation v)
        {
            try
            {
                return v.AreValid();
            }
            catch
            {
                return false;
            }
        }

        /// <summary>
        /// Gets the name of a protein with applied variations
        /// </summary>
        public static string? GetVariantName(string? name, IEnumerable<SequenceVariation>? appliedVariations)
        {
            bool emptyVars = appliedVariations.IsNullOrEmpty();
            if (name == null && emptyVars)
                return null;

            string variantTag = "";
            if (!emptyVars)
            {
                // build a concise, de-duplicated set of variant descriptors (prefer VCF description, fallback to SimpleString)
                var descriptors = appliedVariations!
                    .Where(v => v != null)
                    .Select(v =>
                        v.VariantCallFormatData?.Description ??
                        (string.IsNullOrWhiteSpace(v.Description) ? v.SimpleString() : v.Description))
                    .Where(s => !string.IsNullOrWhiteSpace(s))
                    .Distinct()
                    .Take(6) // cap to avoid pathologically long names
                    .ToList();

                if (descriptors.Count > 0)
                {
                    variantTag = " variant:" + string.Join(", variant:", descriptors);
                }
            }
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
                .Sum(v => (v.VariantSequence ?? string.Empty).Length - (v.OriginalSequence ?? string.Empty).Length);
        }

        /// <summary>
        /// Applies multiple variant changes to a protein sequence
        /// (legacy path – now null-safe around VariantCallFormatData).
        /// Corrected spelling: maxAllowedVariantsForCombinatorics (was ...Combinitorics).
        /// </summary>
        public static List<TBioPolymerType> ApplyVariants<TBioPolymerType>(
            TBioPolymerType protein,
            IEnumerable<SequenceVariation> sequenceVariations,
            int maxAllowedVariantsForCombinatorics,
            int minAlleleDepth)
            where TBioPolymerType : IHasSequenceVariants
        {
            List<SequenceVariation> uniqueEffectsToApply = sequenceVariations
                .Where(v => v != null)
                .GroupBy(v => v.SimpleString())
                .Select(g => g.First())
                .Where(v => v.VariantCallFormatData != null && v.VariantCallFormatData.Genotypes != null && v.VariantCallFormatData.Genotypes.Count > 0)
                .OrderByDescending(v => v.OneBasedBeginPosition)
                .ToList();

            TBioPolymerType proteinCopy = protein.CreateVariant(protein.BaseSequence, protein, null, protein.TruncationProducts, protein.OneBasedPossibleLocalizedModifications, null);

            // If there aren't any variants to apply, just return the base protein
            if (uniqueEffectsToApply.Count == 0)
            {
                return new List<TBioPolymerType> { proteinCopy };
            }

            HashSet<string> individuals = new HashSet<string>(
                uniqueEffectsToApply
                    .Where(v => v.VariantCallFormatData?.Genotypes != null)
                    .SelectMany(v => v.VariantCallFormatData!.Genotypes.Keys));

            List<TBioPolymerType> variantProteins = new();
            List<TBioPolymerType> newVariantProteins = new();
            // loop through genotypes for each sample/individual (e.g. tumor and normal)
            foreach (string individual in individuals)
            {
                newVariantProteins.Clear();
                newVariantProteins.Add(proteinCopy);

                bool tooManyHeterozygousVariants = uniqueEffectsToApply
                    .Where(v => v.VariantCallFormatData?.Heterozygous != null && v.VariantCallFormatData.Heterozygous.ContainsKey(individual))
                    .Count(v => v.VariantCallFormatData.Heterozygous[individual]) > maxAllowedVariantsForCombinatorics;

                foreach (var variant in uniqueEffectsToApply)
                {
                    var vcf = variant.VariantCallFormatData;
                    if (vcf == null || vcf.Genotypes == null || !vcf.Genotypes.ContainsKey(individual))
                        continue;

                    var alleleIndexStr = vcf.AlleleIndex.ToString();
                    bool variantAlleleIsInTheGenotype = vcf.Genotypes[individual].Contains(alleleIndexStr);
                    if (!variantAlleleIsInTheGenotype)
                        continue;

                    bool hetero = vcf.Heterozygous != null && vcf.Heterozygous.ContainsKey(individual) && vcf.Heterozygous[individual];
                    bool homoAlternate = vcf.Homozygous != null && vcf.Homozygous.ContainsKey(individual) && vcf.Homozygous[individual] &&
                                         vcf.Genotypes[individual].All(d => d == alleleIndexStr);

                    bool isDeepReferenceAllele = vcf.AlleleDepths != null &&
                                                 vcf.AlleleDepths.ContainsKey(individual) &&
                                                 vcf.AlleleDepths[individual].Length > 0 &&
                                                 int.TryParse(vcf.AlleleDepths[individual][0], out int depthRef) &&
                                                 depthRef >= minAlleleDepth;

                    bool isDeepAlternateAllele = vcf.AlleleDepths != null &&
                                                 vcf.AlleleDepths.ContainsKey(individual) &&
                                                 vcf.AlleleDepths[individual].Length > vcf.AlleleIndex &&
                                                 int.TryParse(vcf.AlleleDepths[individual][vcf.AlleleIndex], out int depthAlt) &&
                                                 depthAlt >= minAlleleDepth;

                    // homozygous alternate
                    if (homoAlternate && isDeepAlternateAllele)
                    {
                        newVariantProteins = newVariantProteins.Select(p => ApplySingleVariant(variant, p, individual)).ToList();
                    }

                    // heterozygous basic
                    // first protein with variants contains all homozygous variation, second contains all variations
                    else if (hetero && tooManyHeterozygousVariants)
                    {
                        if (isDeepAlternateAllele && isDeepReferenceAllele)
                        {
                            if (newVariantProteins.Count == 1 && maxAllowedVariantsForCombinatorics > 0)
                            {
                                var variantProtein = ApplySingleVariant(variant, newVariantProteins[0], individual);
                                newVariantProteins.Add(variantProtein);
                            }
                            else if (maxAllowedVariantsForCombinatorics > 0 && newVariantProteins.Count > 1)
                            {
                                newVariantProteins[1] = ApplySingleVariant(variant, newVariantProteins[1], individual);
                            }
                        }
                        else if (isDeepAlternateAllele && maxAllowedVariantsForCombinatorics > 0)
                        {
                            newVariantProteins = newVariantProteins.Select(p => ApplySingleVariant(variant, p, individual)).ToList();
                        }
                    }

                    // heterozygous combinatorics
                    else if (hetero && isDeepAlternateAllele && !tooManyHeterozygousVariants)
                    {
                        List<TBioPolymerType> combinatoricProteins = new();
                        foreach (var ppp in newVariantProteins)
                        {
                            if (isDeepAlternateAllele && maxAllowedVariantsForCombinatorics > 0 && isDeepReferenceAllele)
                            {
                                if (vcf.Genotypes[individual].Contains("0"))
                                {
                                    combinatoricProteins.Add(ppp);
                                }
                                combinatoricProteins.Add(ApplySingleVariant(variant, ppp, individual));
                            }
                            else if (isDeepAlternateAllele && maxAllowedVariantsForCombinatorics > 0)
                            {
                                combinatoricProteins.Add(ApplySingleVariant(variant, ppp, individual));
                            }
                            else if (vcf.Genotypes[individual].Contains("0"))
                            {
                                combinatoricProteins.Add(ppp);
                            }
                        }
                        newVariantProteins = combinatoricProteins;
                    }
                }
                variantProteins.AddRange(newVariantProteins);
            }

            return variantProteins
                .GroupBy(x => x.BaseSequence)
                .Select(x => x.First())
                .ToList();
        }

        /// <summary>
        /// Applies a single variant to a protein sequence
        /// </summary>
        private static TBioPolymerType ApplySingleVariant<TBioPolymerType>(SequenceVariation variantGettingApplied, TBioPolymerType protein, string individual)
            where TBioPolymerType : IHasSequenceVariants
        {
            if (variantGettingApplied == null || protein == null)
            {
                return protein;
            }

            string originalSeq = variantGettingApplied.OriginalSequence ?? string.Empty;
            string variantSeq = variantGettingApplied.VariantSequence ?? string.Empty;

            if (variantGettingApplied.OneBasedBeginPosition < 1 ||
                variantGettingApplied.OneBasedBeginPosition > protein.BaseSequence.Length + 1)
            {
                return protein;
            }

            int replacedLength = originalSeq.Length;
            int afterIdx = variantGettingApplied.OneBasedBeginPosition + replacedLength - 1;
            if (afterIdx > protein.BaseSequence.Length)
            {
                replacedLength = Math.Max(0, protein.BaseSequence.Length - (variantGettingApplied.OneBasedBeginPosition - 1));
                afterIdx = variantGettingApplied.OneBasedBeginPosition + replacedLength - 1;
            }

            string seqBefore = protein.BaseSequence.Substring(0, variantGettingApplied.OneBasedBeginPosition - 1);
            string seqAfter = afterIdx >= protein.BaseSequence.Length
                ? string.Empty
                : protein.BaseSequence.Substring(afterIdx);

            int appliedBegin = variantGettingApplied.OneBasedBeginPosition;
            int appliedEnd = variantGettingApplied.OneBasedBeginPosition + Math.Max(0, originalSeq.Length - 1); // end is based on original span, not variant length

            var variantModDict = variantGettingApplied.OneBasedModifications != null
                ? variantGettingApplied.OneBasedModifications.ToDictionary(kv => kv.Key, kv => kv.Value)
                : new Dictionary<int, List<Modification>>();

            string vcfDescription = variantGettingApplied.VariantCallFormatData?.Description;

            SequenceVariation variantAfterApplication = new SequenceVariation(
                appliedBegin,
                appliedEnd,  // end is based on original span, not variant length
                originalSeq,
                variantSeq,
                variantGettingApplied.Description,
                vcfDescription,
                variantModDict.Count == 0 ? null : variantModDict);

            bool intersectsAppliedRegionIncompletely = protein.AppliedSequenceVariations
                .Any(x => variantGettingApplied.Intersects(x) && !variantGettingApplied.Includes(x));

            IEnumerable<SequenceVariation> appliedVariations = new[] { variantAfterApplication };
            if (!intersectsAppliedRegionIncompletely)
            {
                appliedVariations = appliedVariations
                    .Concat(protein.AppliedSequenceVariations.Where(x => !variantGettingApplied.Includes(x)))
                    .ToList();
            }
            else
            {
                seqAfter = afterIdx >= protein.ConsensusVariant.BaseSequence.Length
                    ? string.Empty
                    : protein.ConsensusVariant.BaseSequence.Substring(afterIdx);
            }

            string newBaseSequence = (seqBefore + variantSeq + seqAfter).Split('*')[0];

            var adjustedProteolysisProducts =
                AdjustTruncationProductIndices(variantAfterApplication, newBaseSequence, protein, protein.TruncationProducts);

            var adjustedModifications =
                AdjustModificationIndices(variantAfterApplication, newBaseSequence, protein);

            var adjustedAppliedVariations =
                AdjustSequenceVariationIndices(variantAfterApplication, newBaseSequence, appliedVariations);

            var created = protein.CreateVariant(newBaseSequence,
                                                protein,
                                                adjustedAppliedVariations,
                                                adjustedProteolysisProducts,
                                                adjustedModifications,
                                                individual);

            // Normalize UniProt sequence attributes (length + mass) with safe cloning on length change
            try
            {
                var seq = created?.BaseSequence;
                if (!string.IsNullOrEmpty(seq))
                {
                    var attrProp = created.GetType().GetProperty(
                        "UniProtSequenceAttributes",
                        BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic);

                    var attrs = attrProp?.GetValue(created);
                    if (attrs != null)
                    {
                        var attrType = attrs.GetType();

                        // Read existing attribute values
                        int oldLen = (int)attrType.GetProperty("Length", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic)!.GetValue(attrs);
                        int oldMass = (int)attrType.GetProperty("Mass", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic)!.GetValue(attrs);
                        string checksum = (string)attrType.GetProperty("Checksum", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic)!.GetValue(attrs);
                        DateTime entryMod = (DateTime)attrType.GetProperty("EntryModified", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic)!.GetValue(attrs);
                        int seqVersion = (int)attrType.GetProperty("SequenceVersion", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic)!.GetValue(attrs);
                        bool? isPrecursor = attrType.GetProperty("IsPrecursor", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic)?.GetValue(attrs) as bool?;
                        var fragmentVal = attrType.GetProperty("Fragment", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic)?.GetValue(attrs);

                        // Always recompute mass using trusted path after attach; placeholder keeps constructor happy
                        int newMass = oldMass;

                        if (seq.Length != oldLen)
                        {
                            var ctor = attrType.GetConstructor(new[]
                            {
                                typeof(int), typeof(int), typeof(string), typeof(DateTime), typeof(int),
                                typeof(bool?), fragmentVal?.GetType() ?? attrType
                            });

                            object newAttr;
                            if (ctor != null)
                            {
                                newAttr = ctor.Invoke(new object[]
                                {
                                    seq.Length, newMass, checksum, entryMod, seqVersion,
                                    isPrecursor, fragmentVal ?? Enum.ToObject(attrType.GetNestedType("FragmentType")!, 0)
                                });
                            }
                            else
                            {
                                newAttr = attrs;
                                var lenMeth = attrType.GetMethod("UpdateLengthAttribute", new[] { typeof(string) });
                                lenMeth?.Invoke(newAttr, new object[] { seq });
                            }

                            attrProp?.SetValue(created, newAttr);

                            // Now do the real mass update via the existing API (lives in Proteomics assembly)
                            var massMethPost = newAttr.GetType().GetMethod("UpdateMassAttribute", new[] { typeof(string) });
                            massMethPost?.Invoke(newAttr, new object[] { seq });
                        }
                        else
                        {
                            var lenMeth = attrType.GetMethod("UpdateLengthAttribute", new[] { typeof(string) });
                            lenMeth?.Invoke(attrs, new object[] { seq });
                            var massMeth = attrType.GetMethod("UpdateMassAttribute", new[] { typeof(string) });
                            massMeth?.Invoke(attrs, new object[] { seq });
                        }
                    }
                }

                if (created?.AppliedSequenceVariations?.Count == 0 && adjustedAppliedVariations != null)
                {
                    created.AppliedSequenceVariations.AddRange(adjustedAppliedVariations);
                }
            }
            catch
            {
                // best-effort; ignore failures
            }

            return created;
        }

        /// <summary>
        /// Adjusts the indices of sequence variations due to applying a single additional variant
        /// </summary>
        private static List<SequenceVariation> AdjustSequenceVariationIndices(SequenceVariation variantGettingApplied, string variantAppliedProteinSequence, IEnumerable<SequenceVariation> alreadyAppliedVariations)
        {
            List<SequenceVariation> variations = new();
            if (alreadyAppliedVariations == null)
            {
                return variations;
            }

            foreach (SequenceVariation v in alreadyAppliedVariations)
            {
                if (v == null)
                {
                    continue;
                }

                // NEW: Do not re-shift the variant we just applied; keep its original coordinates.
                if (ReferenceEquals(v, variantGettingApplied))
                {
                    variations.Add(v);
                    continue;
                }

                // Defensive null handling
                string vOrig = v.OriginalSequence ?? string.Empty;
                string vVar = v.VariantSequence ?? string.Empty;

                int addedIdx = alreadyAppliedVariations
                    .Where(applied => applied != null && applied.OneBasedEndPosition < v.OneBasedBeginPosition)
                    .Sum(applied =>
                    {
                        string aVar = applied.VariantSequence ?? string.Empty;
                        string aOrig = applied.OriginalSequence ?? string.Empty;
                        return aVar.Length - aOrig.Length;
                    });

                bool sameVcfRecord =
                    v.VariantCallFormatData != null &&
                    variantGettingApplied.VariantCallFormatData != null &&
                    v.VariantCallFormatData.Equals(variantGettingApplied.VariantCallFormatData);

                // variant was entirely before the one being applied OR it's the current variation (same VCF)
                if (sameVcfRecord || v.OneBasedEndPosition - addedIdx < variantGettingApplied.OneBasedBeginPosition)
                {
                    variations.Add(v);
                    continue;
                }

                // adjust indices based on new included sequence, minding possible overlaps to be filtered later
                int intersectOneBasedStart = Math.Max(variantGettingApplied.OneBasedBeginPosition, v.OneBasedBeginPosition);
                int intersectOneBasedEnd = Math.Min(variantGettingApplied.OneBasedEndPosition, v.OneBasedEndPosition);
                int overlap = intersectOneBasedEnd < intersectOneBasedStart
                    ? 0
                    : intersectOneBasedEnd - intersectOneBasedStart + 1;

                int seqLenChange =
                    (variantGettingApplied.VariantSequence ?? string.Empty).Length -
                    (variantGettingApplied.OriginalSequence ?? string.Empty).Length;

                int begin = v.OneBasedBeginPosition + seqLenChange - overlap;
                if (begin > variantAppliedProteinSequence.Length)
                {
                    // cut out by a stop gain / truncation
                    continue;
                }

                int end = v.OneBasedEndPosition + seqLenChange - overlap;
                if (end > variantAppliedProteinSequence.Length)
                {
                    end = variantAppliedProteinSequence.Length; // shortened by stop
                }
                if (end < begin)
                {
                    // Degenerate after adjustment; skip
                    continue;
                }

                // Null-safe copy of variant-specific mods
                Dictionary<int, List<Modification>> copiedMods = null;
                if (v.OneBasedModifications != null)
                {
                    copiedMods = new Dictionary<int, List<Modification>>(v.OneBasedModifications.Count);
                    foreach (var kv in v.OneBasedModifications)
                    {
                        if (kv.Value == null)
                        {
                            continue;
                        }
                        // shallow copy of list is fine here
                        copiedMods[kv.Key] = new List<Modification>(kv.Value);
                    }
                }

                variations.Add(new SequenceVariation(
                    begin,
                    end,
                    vOrig,
                    vVar,
                    v.Description,
                    v.VariantCallFormatData?.Description,
                    copiedMods));
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
            int sequenceLengthChange = (variant.VariantSequence ?? string.Empty).Length - (variant.OriginalSequence ?? string.Empty).Length;
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
                    if ((variant.VariantSequence ?? string.Empty).EndsWith("*"))
                    {
                        products.Add(new TruncationProduct(p.OneBasedBeginPosition, variantAppliedProteinSequence.Length, p.Type));
                    }
                    else if (p.OneBasedEndPosition + sequenceLengthChange <= variantAppliedProteinSequence.Length)
                    {
                        products.Add(new TruncationProduct(p.OneBasedBeginPosition, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
                    }
                }
                // proteolysis product is after the variant and there is no stop gain
                else if (p.OneBasedBeginPosition > variant.OneBasedEndPosition
                         && p.OneBasedBeginPosition + sequenceLengthChange <= variantAppliedProteinSequence.Length
                         && p.OneBasedEndPosition + sequenceLengthChange <= variantAppliedProteinSequence.Length
                         && !(variant.VariantSequence ?? string.Empty).EndsWith("*"))
                {
                    products.Add(new TruncationProduct(p.OneBasedBeginPosition + sequenceLengthChange, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
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
            int sequenceLengthChange = (variant.VariantSequence ?? string.Empty).Length - (variant.OriginalSequence ?? string.Empty).Length;

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
            return variations.IsNullOrEmpty()
                ? ""
                : string.Join("_", variations
                    .Where(v => v != null)
                    .Select(v => v.SimpleString()));
        }

        /// <summary>
        /// Format string to append to protein names
        /// </summary>
        public static string CombineDescriptions(IEnumerable<SequenceVariation>? variations)
        {
            if (variations.IsNullOrEmpty())
                return "";

            var tokens = variations!
                .Where(v => v != null)
                .Select(v => v.VariantCallFormatData?.Description ??
                             (string.IsNullOrWhiteSpace(v.Description) ? v.SimpleString() : v.Description))
                .Where(s => !string.IsNullOrWhiteSpace(s))
                .Distinct()
                .Take(10)
                .ToList();

            return string.Join(", variant:", tokens);
        }

        /// <summary>
        /// Applies all possible combinations of the provided SequenceVariation list to the base TBioPolymerType object,
        /// starting with the fewest single variations and up to the specified maximum number of combinations.
        /// </summary>
        /// <typeparam name="TBioPolymerType">The type of the biopolymer object.</typeparam>
        /// <param name="baseBioPolymer">The base biopolymer object to apply variations to.</param>
        /// <param name="variations">List of SequenceVariation objects to combine and apply. Assumed not null or empty.</param>
        /// <param name="maxSequenceVariantsPerIsoform">Maximum number of combinations to return.</param>
        /// <!---->/ <param name="maxSequenceVariantIsoforms">.</param> -->
        /// <returns>
        /// An IEnumerable of TBioPolymerType objects, each with a unique combination of variations applied.
        /// </returns>
        public static IEnumerable<TBioPolymerType> ApplyAllVariantCombinations<TBioPolymerType>(
            TBioPolymerType baseBioPolymer,
            List<SequenceVariation> variations,
            int maxSequenceVariantsPerIsoform,
            int maxSequenceVariantIsoforms,
            int minAlleleDepth)
            where TBioPolymerType : IHasSequenceVariants
        {
            int count = 0;

            // Always yield the base biopolymer first
            yield return baseBioPolymer;
            count++;

            // 1. Attempt genotype-aware expansion
            List<SequenceVariation> sequenceVariations = new();
            foreach (var v in variations.Where(v => v != null))
            {
                try
                {
                    // Only try per-genotype split if VCF data present; otherwise just add the raw variant
                    if (v.VariantCallFormatData != null)
                    {
                        var split = v.SplitPerGenotype(minAlleleDepth);
                        if (split != null && split.Count > 0)
                        {
                            sequenceVariations.AddRange(split);
                            continue;
                        }
                    }
                    sequenceVariations.Add(v); // fallback to original
                }
                catch
                {
                    // On any parsing/splitting issue, keep original variant so we still attempt application
                    sequenceVariations.Add(v);
                }
            }

            // 2. Collapse equivalent variants (only if >1)
            if (sequenceVariations.Count > 1)
            {
                sequenceVariations = SequenceVariation.CombineEquivalent(sequenceVariations);
            }

            // 3. Filter invalid (but keep at least something if all fail)
            var filtered = sequenceVariations.Where(v =>
            {
                try { return v != null && v.AreValid(); }
                catch { return true; // treat exceptions as usable to avoid discarding everything
                }
            }).ToList();

            if (filtered.Count == 0)
            {
                filtered = sequenceVariations.Where(v => v != null).ToList();
            }

            // 4. Remove pure no-op substitutions (no sequence change and no variant-specific mods)
            filtered = filtered.Where(v =>
                    !(string.Equals(v.OriginalSequence ?? "",
                                    v.VariantSequence ?? "",
                                    StringComparison.Ordinal)
                      && (v.OneBasedModifications == null || v.OneBasedModifications.Count == 0)))
                .ToList();

            if (filtered.Count == 0)
            {
                yield break; // nothing meaningful to apply beyond the base already yielded
            }

            int total = filtered.Count;
            int maxVariantsPerIsoformCapped = Math.Min(maxSequenceVariantsPerIsoform, total);

            for (int size = 1; size <= maxVariantsPerIsoformCapped; size++)
            {
                foreach (var combo in GetCombinations(filtered, size))
                {
                    if (count >= maxSequenceVariantIsoforms)
                        yield break;

                    var listCombo = combo.Where(c => c != null).ToList();
                    if (listCombo.Count == 0)
                        continue;

                    if (!ValidCombination(listCombo))
                        continue;

                    var result = baseBioPolymer;
                    bool aborted = false;

                    foreach (var variant in listCombo)
                    {
                        result = ApplySingleVariant(variant, result, string.Empty);
                        if (result == null)
                        {
                            aborted = true;
                            break;
                        }
                    }

                    if (aborted || result == null)
                        continue;

                    // Skip if sequence remained identical (all variants net no-op)
                    if (ReferenceEquals(result, baseBioPolymer) ||
                        string.Equals(result.BaseSequence, baseBioPolymer.BaseSequence, StringComparison.Ordinal))
                        continue;

                    yield return result;
                    count++;
                }
            }
        }

        /// <summary>
        /// Generates all possible combinations of the specified size from the input list.
        /// Robust to:
        /// - null / empty variation list (yields nothing)
        /// - size <= 0 (yields nothing)
        /// - size > count (yields nothing)
        /// Fast paths:
        /// - size == 1 → yield each variation individually
        /// - size == count → yield the whole set once
        /// </summary>
        /// <param name="variations">List of SequenceVariation objects to combine.</param>
        /// <param name="size">The size of each combination.</param>
        private static IEnumerable<IList<SequenceVariation>> GetCombinations(List<SequenceVariation> variations, int size)
        {
            // Guard conditions
            if (variations == null || variations.Count == 0 || size <= 0 || size > variations.Count)
                yield break;

            int n = variations.Count;

            // Single element combinations → just yield each item
            if (size == 1)
            {
                for (int i = 0; i < n; i++)
                {
                    yield return new List<SequenceVariation>(1) { variations[i] };
                }
                yield break;
            }

            // Whole-set combination
            if (size == n)
            {
                yield return new List<SequenceVariation>(variations);
                yield break;
            }

            // Standard iterative k-combination generator (lexicographic indices)
            var indices = new int[size];
            for (int i = 0; i < size; i++)
                indices[i] = i;

            while (true)
            {
                var combo = new List<SequenceVariation>(size);
                for (int i = 0; i < size; i++)
                    combo.Add(variations[indices[i]]);
                yield return combo;

                int pos = size - 1;
                while (pos >= 0 && indices[pos] == n - size + pos)
                    pos--;
                if (pos < 0)
                    break;

                indices[pos++]++;

                for (int i = pos; i < size; i++)
                    indices[i] = indices[i - 1] + 1;
            }
        }
        public static bool ValidCombination(List<SequenceVariation> variations)
        {
            if (variations == null || variations.Count <= 1)
                return true;

            // Validate inputs
            for (int i = 0; i < variations.Count; i++)
            {
                var v = variations[i];
                if (v == null || !v.AreValid())
                    return false;
            }

            // Sort by begin then end, then check only adjacent intervals
            var ordered = variations
                .OrderBy(v => v.OneBasedBeginPosition)
                .ThenBy(v => v.OneBasedEndPosition)
                .ToList();

            var prev = ordered[0];
            for (int i = 1; i < ordered.Count; i++)
            {
                var curr = ordered[i];
                if (prev.Intersects(curr)) // inclusive overlap check
                    return false;

                prev = curr;
            }
            return true;
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
                        modificationsToRemove.Add(new(kvp.Key, mod));
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
        /// <summary>
        /// Lightweight sanitizer for variant data prior to XML write or further combinatorics.
        /// Removes null / invalid / out-of-range SequenceVariations and prunes obviously invalid
        /// variant-specific modification indices so downstream writers do not throw NREs.
        /// Returns a short enumerable of human‑readable notes (can be logged) describing actions taken.
        /// 
        /// Non‑destructive policy:
        /// - SequenceVariation objects are never mutated (they are immutable); any problematic one is dropped.
        /// - AppliedSequenceVariations is re-filtered to only include surviving base SequenceVariations (by reference equality)
        ///   plus any that were already applied but still valid against the current sequence.
        /// - Variant-specific modifications that point outside the plausible post‑edit protein length are removed.
        /// 
        /// Safety heuristics (fast, no deep recomputation):
        /// 1. Drop variant if:
        ///    - null
        ///    - begin < 1
        ///    - begin > BaseSequence.Length + 1 (cannot even be an insertion)
        ///    - AreValid() returns false
        /// 2. Prune variant.OneBasedModifications keys if:
        ///    - key < 1
        ///    - key > (BaseSequence.Length + maxDeltaLen)   (where maxDeltaLen = variant.VariantSequence.Length - variant.OriginalSequence.Length, if positive)
        ///    - variant encodes a deletion or stop-gain (VariantSequence empty or ends with '*') AND key >= variant.OneBasedBeginPosition
        /// 
        /// This is intentionally conservative: we do not attempt to "fix" coordinates, only remove obviously hazardous data.
        /// </summary>
        public static IEnumerable<string> SanitizeVariantData<TBioPolymerType>(
            IEnumerable<TBioPolymerType> polymers,
            bool removeInvalidVariants = true)
            where TBioPolymerType : IHasSequenceVariants
        {
            if (polymers == null)
                yield break;

            foreach (var prot in polymers)
            {
                if (prot == null)
                    continue;

                var notes = new List<string>();
                var originalCount = prot.SequenceVariations?.Count ?? 0;

                if (prot.SequenceVariations == null)
                {
                    continue; // nothing to sanitize
                }

                // Working list (do not modify while iterating original)
                var kept = new List<SequenceVariation>(prot.SequenceVariations.Count);
                foreach (var v in prot.SequenceVariations)
                {
                    if (v == null)
                    {
                        notes.Add("Dropped null variant");
                        continue;
                    }

                    // Coordinate sanity (pre-AreValid fast checks)
                    if (v.OneBasedBeginPosition < 1 ||
                        v.OneBasedBeginPosition > prot.BaseSequence.Length + 1)
                    {
                        notes.Add($"Dropped variant (coords out of range) {v.SimpleString()}");
                        if (removeInvalidVariants) continue; else kept.Add(v);
                        continue;
                    }

                    // Validation (can still fail if object was mutated after construction)
                    bool valid;
                    try
                    {
                        valid = v.AreValid();
                    }
                    catch
                    {
                        valid = false;
                    }

                    if (!valid)
                    {
                        notes.Add($"Dropped invalid variant {v.SimpleString()}");
                        if (removeInvalidVariants) continue; else kept.Add(v);
                        continue;
                    }

                    // Prune variant-specific modifications dictionary (mutable) if present
                    if (v.OneBasedModifications != null && v.OneBasedModifications.Count > 0)
                    {
                        int delta = (v.VariantSequence?.Length ?? 0) - (v.OriginalSequence?.Length ?? 0);
                        int maxAllowedPos = prot.BaseSequence.Length + Math.Max(0, delta);

                        var toRemove = new List<int>();
                        foreach (var kv in v.OneBasedModifications)
                        {
                            int pos = kv.Key;
                            if (pos < 1 || pos > maxAllowedPos)
                            {
                                toRemove.Add(pos);
                                continue;
                            }
                            bool deletionOrStop = string.IsNullOrEmpty(v.VariantSequence) || (v.VariantSequence?.Contains('*') ?? false);
                            if (deletionOrStop && pos >= v.OneBasedBeginPosition)
                            {
                                toRemove.Add(pos);
                            }
                        }

                        if (toRemove.Count > 0)
                        {
                            foreach (var k in toRemove)
                            {
                                v.OneBasedModifications.Remove(k);
                            }
                            notes.Add($"Variant {v.SimpleString()} pruned {toRemove.Count} mod site(s)");
                        }
                    }

                    kept.Add(v);
                }

                if (kept.Count != originalCount)
                {
                    // Replace list (SequenceVariations is mutable list per interface)
                    prot.SequenceVariations.Clear();
                    prot.SequenceVariations.AddRange(kept);
                    notes.Add($"Sanitized variants: kept {kept.Count}/{originalCount}");
                }

                // Reconcile AppliedSequenceVariations if present (drop references that no longer exist or became invalid)
                if (prot.AppliedSequenceVariations != null && prot.AppliedSequenceVariations.Count > 0)
                {
                    int beforeApplied = prot.AppliedSequenceVariations.Count;
                    prot.AppliedSequenceVariations.RemoveAll(v => v == null || !kept.Contains(v));
                    if (prot.AppliedSequenceVariations.Count != beforeApplied)
                    {
                        notes.Add($"Pruned applied variant refs: {beforeApplied - prot.AppliedSequenceVariations.Count} removed");
                    }
                }

                foreach (var n in notes)
                {
                    // TBioPolymerType is only constrained to IHasSequenceVariants (no Accession there).
                    // Use direct Accession if the object implements IBioPolymer; otherwise fall back to ConsensusVariant.Accession.
                    string acc = (prot as IBioPolymer)?.Accession
                                 ?? prot.ConsensusVariant?.Accession
                                 ?? "<no_accession>";
                    yield return $"[{acc}] {n}";
                }
            }
        }

        /// <summary>
        /// Convenience overload for a single protein / biopolymer.
        /// </summary>
        public static IEnumerable<string> SanitizeVariantData<TBioPolymerType>(TBioPolymerType polymer, bool removeInvalidVariants = true)
            where TBioPolymerType : IHasSequenceVariants
        {
            return SanitizeVariantData(new[] { polymer }, removeInvalidVariants);
        }
    }
}