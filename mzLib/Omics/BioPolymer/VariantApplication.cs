using MzLibUtil;
using Omics.BioPolymer;
using Omics.Modifications;

namespace Omics.BioPolymer
{
    /// <summary>
    /// Utilities for constructing and applying sequence variants to biopolymers (proteins/RNA) in a 1-based, inclusive coordinate system.
    /// This includes:
    /// - Expanding a base biopolymer into concrete variant instances (genotype-aware or combinatorial).
    /// - Renaming/accession tagging based on applied variants.
    /// - Re-mapping indices for truncation products and localized modifications after edits (insertions/deletions/substitutions).
    /// - Converting certain annotation-style modifications (e.g., nucleotide substitution markers) into true sequence variants.
    /// </summary>
    public static class VariantApplication
    {
        /// <summary>
        /// Builds concrete variant biopolymers from a base entity.
        /// If any known variation lacks genotype information (no VCF or empty genotypes), a safe combinatorial expansion is used (bounded by <paramref name="maxAllowedVariantsForCombinatorics"/>).
        /// Otherwise, variants are applied in a genotype/allele-depth-aware manner via <see cref="ApplyVariants{TBioPolymerType}(TBioPolymerType, IEnumerable{SequenceVariation}, int, int)"/>.
        /// </summary>
        /// <typeparam name="TBioPolymerType">A biopolymer type that supports sequence variants.</typeparam>
        /// <param name="protein">The base biopolymer to expand into variant instances.</param>
        /// <param name="maxAllowedVariantsForCombinatorics">
        /// Maximum number of variants to consider when creating combinations for genotype-less scenarios.
        /// This caps explosion in <see cref="ApplyAllVariantCombinations{TBioPolymerType}(TBioPolymerType, System.Collections.Generic.List{SequenceVariation}, int)"/>.
        /// </param>
        /// <param name="minAlleleDepth">
        /// Minimum AD (Allele Depth) per sample required for an allele (reference or alternate) to be considered "deep" enough to participate in genotype-aware application.
        /// </param>
        /// <returns>A list of concrete variants derived from <paramref name="protein"/>.</returns>
        public static List<TBioPolymerType> GetVariantBioPolymers<TBioPolymerType>(this TBioPolymerType protein, int maxAllowedVariantsForCombinatorics = 4, int minAlleleDepth = 1)
            where TBioPolymerType : IHasSequenceVariants
        {
            // Materialize any substitution-like annotations into concrete sequence variants on both the consensus and the instance
            protein.ConsensusVariant.ConvertNucleotideSubstitutionModificationsToSequenceVariants();
            protein.ConvertNucleotideSubstitutionModificationsToSequenceVariants();

            // If all variants are positionally valid and any is missing VCF/genotype data, fall back to bounded combinatorial application
            if (protein.SequenceVariations.All(v => v.AreValid()) && protein.SequenceVariations.Any(v => v.VariantCallFormatDataString == null || v.VariantCallFormatDataString.Genotypes.Count == 0))
            {
                // this is a protein with either no VCF lines or a mix of VCF and non-VCF lines
                return ApplyAllVariantCombinations(protein, protein.SequenceVariations, maxAllowedVariantsForCombinatorics).ToList(); 
            }

            // Otherwise, do genotype/allele-depth-aware application with combinatorics limited for heterozygous sites
            return ApplyVariants(protein, protein.SequenceVariations, maxAllowedVariantsForCombinitorics: maxAllowedVariantsForCombinatorics, minAlleleDepth);
        }

        /// <summary>
        /// Produces a name with an appended variant tag built from applied variation descriptions.
        /// If both <paramref name="name"/> and <paramref name="appliedVariations"/> are effectively empty, returns null.
        /// </summary>
        /// <param name="name">Base name (e.g., protein name). May be null.</param>
        /// <param name="appliedVariations">Variations applied to the instance. If null/empty, no variant tag is appended.</param>
        /// <returns>
        /// The base <paramref name="name"/> plus " variant:{...}" when variations exist; null if both inputs are empty.
        /// </returns>
        public static string? GetVariantName(string? name, IEnumerable<SequenceVariation>? appliedVariations)
        {
            bool emptyVars = appliedVariations.IsNullOrEmpty();
            if (name == null && emptyVars)
                return null;

            string variantTag = emptyVars ? "" : $" variant:{CombineDescriptions(appliedVariations)}";
            return name + variantTag;
        }

        /// <summary>
        /// Constructs a variant-aware accession by appending a compact variant string to the consensus accession.
        /// </summary>
        /// <param name="protein">The variant-capable biopolymer.</param>
        /// <param name="appliedSequenceVariations">Applied variations to encode into the suffix; may be null/empty.</param>
        /// <returns>
        /// The consensus accession or consensus accession with "_{SimpleString}_..." suffix if variations are present.
        /// </returns>
        public static string GetAccession(IHasSequenceVariants protein, IEnumerable<SequenceVariation>? appliedSequenceVariations)
        {
            return protein.ConsensusVariant.Accession +
                   (appliedSequenceVariations.IsNullOrEmpty() ? "" : $"_{CombineSimpleStrings(appliedSequenceVariations)}");
        }

        /// <summary>
        /// Determines if a specific 1-based position in the variant biopolymer lies within a particular variation's range.
        /// </summary>
        /// <param name="appliedVariant">The variation of interest; may be null.</param>
        /// <param name="variantProteinIndex">1-based position in the current (possibly already-edited) variant sequence.</param>
        /// <returns>True if the position is included by the variation; otherwise false.</returns>
        public static bool IsSequenceVariantModification(SequenceVariation? appliedVariant, int variantProteinIndex)
        {
            return appliedVariant != null && appliedVariant.Includes(variantProteinIndex);
        }

        /// <summary>
        /// Maps a modification index from the edited (variant) sequence back to the original consensus index by subtracting the
        /// net length changes of all applied variations that ended before the queried position.
        /// </summary>
        /// <param name="protein">The variant-capable biopolymer containing applied variations.</param>
        /// <param name="variantProteinModificationIndex">1-based index in the variant sequence.</param>
        /// <returns>The corresponding 1-based index in the original consensus sequence.</returns>
        public static int RestoreModificationIndex(IHasSequenceVariants protein, int variantProteinModificationIndex)
        {
            return variantProteinModificationIndex - protein.AppliedSequenceVariations
                .Where(v => v.OneBasedEndPosition < variantProteinModificationIndex)
                .Sum(v => v.VariantSequence.Length - v.OriginalSequence.Length);
        }

        /// <summary>
        /// Applies a set of sequence variations in a genotype- and allele-depth-aware fashion for all individuals found in the VCF payloads.
        /// Heterozygous sites can produce combinatorial branches up to <paramref name="maxAllowedVariantsForCombinitorics"/> per individual.
        /// Results are deduplicated by final base sequence.
        /// </summary>
        /// <typeparam name="TBioPolymerType">A biopolymer type that supports sequence variants.</typeparam>
        /// <param name="protein">The base biopolymer to which variations will be applied.</param>
        /// <param name="sequenceVariations">Candidate variations. Duplicates by effect are collapsed using <see cref="SequenceVariation.SimpleString"/>.</param>
        /// <param name="maxAllowedVariantsForCombinitorics">
        /// Upper cap for heterozygous combinatorial branching. If an individual has more heterozygous variants than this number,
        /// the algorithm limits branching to control explosion.
        /// </param>
        /// <param name="minAlleleDepth">Minimum AD (Allele Depth) per sample for an allele to be considered in application.</param>
        /// <returns>A list of concrete variant biopolymers across all individuals encoded in the VCF payloads.</returns>
        public static List<TBioPolymerType> ApplyVariants<TBioPolymerType>(TBioPolymerType protein, IEnumerable<SequenceVariation> sequenceVariations, int maxAllowedVariantsForCombinitorics, int minAlleleDepth)
            where TBioPolymerType : IHasSequenceVariants
        {
            // Remove duplicate effects (by SimpleString), require variants with genotype data, apply from higher to lower positions
            List<SequenceVariation> uniqueEffectsToApply = sequenceVariations
                .GroupBy(v => v.SimpleString())
                .Select(x => x.First())
                .Where(v => v.VariantCallFormatDataString.Genotypes.Count > 0)
                .OrderByDescending(v => v.OneBasedBeginPosition)
                .ToList();

            // A shallow "base" variant to branch from (no applied variants yet)
            TBioPolymerType proteinCopy = protein.CreateVariant(protein.BaseSequence, protein, null, protein.TruncationProducts, protein.OneBasedPossibleLocalizedModifications, null);

            if (uniqueEffectsToApply.Count == 0)
            {
                return new List<TBioPolymerType> { proteinCopy };
            }

            // All per-sample identifiers present in the VCF objects
            HashSet<string> individuals = new HashSet<string>(uniqueEffectsToApply.SelectMany(v => v.VariantCallFormatDataString.Genotypes.Keys));
            List<TBioPolymerType> variantProteins = new();
            List<TBioPolymerType> newVariantProteins = new();

            foreach (string individual in individuals)
            {
                newVariantProteins.Clear();
                newVariantProteins.Add(proteinCopy);

                // Whether to limit combinatorial branching for this individual
                bool tooManyHeterozygousVariants = uniqueEffectsToApply.Count(v => v.VariantCallFormatDataString.Heterozygous[individual]) > maxAllowedVariantsForCombinitorics;

                foreach (var variant in uniqueEffectsToApply)
                {
                    // Only proceed if the individual's genotype references this variant's alternate allele index
                    bool variantAlleleIsInTheGenotype = variant.VariantCallFormatDataString.Genotypes[individual].Contains(variant.VariantCallFormatDataString.AlleleIndex.ToString());
                    if (!variantAlleleIsInTheGenotype)
                    {
                        continue;
                    }

                    // Zygosity and depth checks for this individual
                    bool isHomozygousAlternate = variant.VariantCallFormatDataString.Homozygous[individual] && variant.VariantCallFormatDataString.Genotypes[individual].All(d => d == variant.VariantCallFormatDataString.AlleleIndex.ToString());
                    bool isDeepReferenceAllele = int.TryParse(variant.VariantCallFormatDataString.AlleleDepths[individual][0], out int depthRef) && depthRef >= minAlleleDepth;
                    bool isDeepAlternateAllele = int.TryParse(variant.VariantCallFormatDataString.AlleleDepths[individual][variant.VariantCallFormatDataString.AlleleIndex], out int depthAlt) && depthAlt >= minAlleleDepth;

                    if (isHomozygousAlternate && isDeepAlternateAllele)
                    {
                        // Deterministic application: all branches take the alt allele
                        newVariantProteins = newVariantProteins.Select(p => ApplySingleVariant(variant, p, individual)).ToList();
                    }
                    else if (variant.VariantCallFormatDataString.Heterozygous[individual] && tooManyHeterozygousVariants)
                    {
                        // Limit branching: either keep ref, take alt, or update second branch if already present
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
                        // Full combinatorics (bounded) for heterozygous case: keep ref and/or take alt depending on depths
                        List<TBioPolymerType> combinitoricProteins = new();

                        foreach (var ppp in newVariantProteins)
                        {
                            if (isDeepAlternateAllele && maxAllowedVariantsForCombinitorics > 0 && isDeepReferenceAllele)
                            {
                                if (variant.VariantCallFormatDataString.Genotypes[individual].Contains("0"))
                                {
                                    combinitoricProteins.Add(ppp); // keep reference branch
                                }
                                combinitoricProteins.Add(ApplySingleVariant(variant, ppp, individual)); // alternate branch
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

            // De-duplicate by final sequence to avoid identical variants from different branches
            return variantProteins.GroupBy(x => x.BaseSequence).Select(x => x.First()).ToList();
        }

        /// <summary>
        /// Applies a single <see cref="SequenceVariation"/> to a biopolymer and returns the updated instance with
        /// sequence, truncation products, localized modifications, and applied-variation indices all re-mapped.
        /// </summary>
        /// <typeparam name="TBioPolymerType">A biopolymer type that supports sequence variants.</typeparam>
        /// <param name="variantGettingApplied">Variation to apply (coordinates refer to the current protein's sequence).</param>
        /// <param name="protein">Source biopolymer to mutate.</param>
        /// <param name="individual">Sample identifier used to annotate the created variant (may be empty).</param>
        /// <returns>A new variant instance created via <see cref="IHasSequenceVariants.CreateVariant{TBioPolymerType}"/>.</returns>
        private static TBioPolymerType ApplySingleVariant<TBioPolymerType>(SequenceVariation variantGettingApplied, TBioPolymerType protein, string individual)
            where TBioPolymerType : IHasSequenceVariants
        {
            // Sequence prefix before the edit region
            string seqBefore = protein.BaseSequence.Substring(0, variantGettingApplied.OneBasedBeginPosition - 1);
            // The replacement (alternate) sequence
            string seqVariant = variantGettingApplied.VariantSequence;
            // First index in the original sequence after the edited region
            int afterIdx = variantGettingApplied.OneBasedBeginPosition + variantGettingApplied.OriginalSequence.Length - 1;

            // Reify a "post-application" variation object pinned to the inserted length
            var vcf = variantGettingApplied.VariantCallFormatDataString;
            SequenceVariation variantAfterApplication = new SequenceVariation(
                variantGettingApplied.OneBasedBeginPosition,
                variantGettingApplied.OneBasedBeginPosition + variantGettingApplied.VariantSequence.Length - 1,
                variantGettingApplied.OriginalSequence,
                variantGettingApplied.VariantSequence,
                vcf != null ? vcf.Description : variantGettingApplied.Description,
                vcf,
                variantGettingApplied.OneBasedModifications.ToDictionary(kv => kv.Key, kv => kv.Value));

            // If an already-applied variation partially overlaps the current edit, use the consensus tail to avoid index corruption
            bool intersectsAppliedRegionIncompletely =
                protein.AppliedSequenceVariations != null
                && protein.AppliedSequenceVariations.Any(x => variantGettingApplied.Intersects(x) && !variantGettingApplied.Includes(x));

            IEnumerable<SequenceVariation> appliedVariations = new[] { variantAfterApplication };
            string seqAfter;
            if (intersectsAppliedRegionIncompletely)
            {
                // Tail from consensus (not the possibly already-mutated BaseSequence)
                seqAfter = protein.BaseSequence.Length - afterIdx <= 0 ? "" : protein.ConsensusVariant.BaseSequence.Substring(afterIdx);
            }
            else
            {
                // Tail from the current BaseSequence; keep any previously applied variations that are not fully contained by the current edit
                seqAfter = protein.BaseSequence.Length - afterIdx <= 0 ? "" : protein.BaseSequence.Substring(afterIdx);
                appliedVariations = appliedVariations
                    .Concat((protein.AppliedSequenceVariations ?? Enumerable.Empty<SequenceVariation>()).Where(x => !variantGettingApplied.Includes(x)))
                    .ToList();
            }

            // Clip at stop (*) if any
            string variantSequence = (seqBefore + seqVariant + seqAfter).Split('*')[0];

            // Re-map dependent structures after the sequence length change
            List<TruncationProduct> adjustedProteolysisProducts = AdjustTruncationProductIndices(variantGettingApplied, variantSequence, protein, protein.TruncationProducts);
            Dictionary<int, List<Modification>> adjustedModifications = AdjustModificationIndices(variantGettingApplied, variantSequence, protein);
            List<SequenceVariation> adjustedAppliedVariations = AdjustSequenceVariationIndices(variantGettingApplied, variantSequence, appliedVariations);

            return protein.CreateVariant(variantSequence, protein, adjustedAppliedVariations, adjustedProteolysisProducts, adjustedModifications, individual);
        }

        /// <summary>
        /// Adjusts (re-bases) the coordinates of already-applied variations after applying a new variation.
        /// Handles overlaps by trimming and shifts by adding the net length change.
        /// </summary>
        /// <param name="variantGettingApplied">The newly applied variation causing coordinate shifts.</param>
        /// <param name="variantAppliedProteinSequence">The updated sequence after the application (used to clamp ends).</param>
        /// <param name="alreadyAppliedVariations">Variations that were previously applied (may be null).</param>
        /// <returns>A new list of variations with updated coordinates, filtered for validity.</returns>
        private static List<SequenceVariation> AdjustSequenceVariationIndices(SequenceVariation variantGettingApplied, string variantAppliedProteinSequence, IEnumerable<SequenceVariation> alreadyAppliedVariations)
        {
            List<SequenceVariation> variations = new List<SequenceVariation>();
            if (alreadyAppliedVariations == null) { return variations; }

            foreach (SequenceVariation v in alreadyAppliedVariations)
            {
                // Net length already introduced before the start of v
                int addedIdx = alreadyAppliedVariations
                    .Where(applied => applied.OneBasedEndPosition < v.OneBasedBeginPosition)
                    .Sum(applied => applied.VariantSequence.Length - applied.OriginalSequence.Length);

                // Either same VCF payload (null-safe) or fully before the new application region after compensating for prior shifts
                if ((v.VariantCallFormatDataString != null && v.VariantCallFormatDataString.Equals(variantGettingApplied.VariantCallFormatDataString))
                    || v.OneBasedEndPosition - addedIdx < variantGettingApplied.OneBasedBeginPosition)
                {
                    variations.Add(v);
                }
                else
                {
                    // Compute overlap with the newly applied edit and shift by the net length change
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

        /// <summary>
        /// Re-bases proteolysis/truncation product ranges after applying a sequence variation.
        /// Shifts segments to the right if the edit occurs before them and expands/contracts segments that span the edit.
        /// If a stop (*) is introduced by the variant, downstream products are clamped to the new sequence length.
        /// </summary>
        /// <param name="variant">The applied variation (source of positional shifts and length change).</param>
        /// <param name="variantAppliedProteinSequence">The updated sequence after application.</param>
        /// <param name="protein">The source biopolymer (used for consensus length in certain boundary checks).</param>
        /// <param name="proteolysisProducts">Existing truncation products to be re-based; may be null.</param>
        /// <returns>Updated list of truncation products within the new coordinate system.</returns>
        private static List<TruncationProduct> AdjustTruncationProductIndices(SequenceVariation variant, string variantAppliedProteinSequence, IHasSequenceVariants protein, IEnumerable<TruncationProduct> proteolysisProducts)
        {
            List<TruncationProduct> products = new List<TruncationProduct>();
            if (proteolysisProducts == null) { return products; }

            int sequenceLengthChange = variant.VariantSequence.Length - variant.OriginalSequence.Length;

            foreach (TruncationProduct p in proteolysisProducts.Where(p => p.OneBasedEndPosition.HasValue && p.OneBasedBeginPosition.HasValue))
            {
                // Entirely before the edit: unchanged
                if (variant.OneBasedBeginPosition > p.OneBasedEndPosition)
                {
                    products.Add(p);
                }
                // Segment spans the edit or is clamped at boundaries: extend/contract or clamp to stop
                else if ((p.OneBasedBeginPosition < variant.OneBasedBeginPosition || p.OneBasedBeginPosition == 1 || p.OneBasedBeginPosition == 2)
                         && (p.OneBasedEndPosition > variant.OneBasedEndPosition || p.OneBasedEndPosition == protein.ConsensusVariant.BaseSequence.Length))
                {
                    if (variant.VariantSequence.EndsWith("*"))
                    {
                        // Introduced stop codon/terminator: clamp to the new sequence length
                        products.Add(new TruncationProduct(p.OneBasedBeginPosition, variantAppliedProteinSequence.Length, p.Type));
                    }
                    else if (p.OneBasedEndPosition + sequenceLengthChange <= variantAppliedProteinSequence.Length)
                    {
                        products.Add(new TruncationProduct(p.OneBasedBeginPosition, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
                    }
                }
                // Entirely after the edit: shift right/left by the net length change, if still within bounds and not terminated
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

        /// <summary>
        /// Produces a new localized modification dictionary re-based to the post-variation coordinate system and merged with
        /// any one-based modifications carried by the applied variant itself.
        /// </summary>
        /// <param name="variant">The variation causing index shifts.</param>
        /// <param name="variantAppliedProteinSequence">The updated sequence after application (for bounds checks).</param>
        /// <param name="protein">The source biopolymer providing the original modifications map.</param>
        /// <returns>A new dictionary of one-based modification lists keyed by position in the updated sequence.</returns>
        private static Dictionary<int, List<Modification>> AdjustModificationIndices(SequenceVariation variant, string variantAppliedProteinSequence, IHasSequenceVariants protein)
        {
            // Original per-position modifications (pre-application)
            IDictionary<int, List<Modification>> modificationDictionary = protein.OneBasedPossibleLocalizedModifications;
            // Modifications contributed by the variant itself (coordinated in the same one-based space)
            IDictionary<int, List<Modification>> variantModificationDictionary = variant.OneBasedModifications;

            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            int sequenceLengthChange = variant.VariantSequence.Length - variant.OriginalSequence.Length;

            // Re-base original modifications
            if (modificationDictionary != null)
            {
                foreach (KeyValuePair<int, List<Modification>> kv in modificationDictionary)
                {
                    if (kv.Key > variantAppliedProteinSequence.Length)
                    {
                        continue; // drop if beyond new end
                    }
                    else if (kv.Key < variant.OneBasedBeginPosition)
                    {
                        mods.Add(kv.Key, kv.Value); // unaffected positions
                    }
                    else if (variant.OneBasedEndPosition < kv.Key && kv.Key + sequenceLengthChange <= variantAppliedProteinSequence.Length)
                    {
                        mods.Add(kv.Key + sequenceLengthChange, kv.Value); // shift after the edit
                    }
                }
            }

            // Merge-in variant-borne modifications (may share positions)
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
        /// Concatenates the compact representations (<see cref="SequenceVariation.SimpleString"/>) of the provided variations with underscores.
        /// </summary>
        /// <param name="variations">Variations to stringify; may be null/empty.</param>
        /// <returns>An underscore-joined string or empty string if none.</returns>
        private static string CombineSimpleStrings(IEnumerable<SequenceVariation>? variations)
        {
            return variations.IsNullOrEmpty() ? "" : string.Join("_", variations.Select(v => v.SimpleString()));
        }

        /// <summary>
        /// Concatenates human-readable descriptions for the provided variations.
        /// Prefers VCF description when available, otherwise falls back to the variation's own description.
        /// </summary>
        /// <param name="variations">Variations to describe; may be null/empty.</param>
        /// <returns>A comma-delimited description string or empty string if none.</returns>
        public static string CombineDescriptions(IEnumerable<SequenceVariation>? variations)
        {
            return variations.IsNullOrEmpty() ? "" : string.Join(", variant:", variations.Select(d => d.VariantCallFormatDataString?.Description ?? d.Description));
        }

        /// <summary>
        /// Applies all combinations of the provided variations to the base biopolymer up to a maximum number of yielded results.
        /// The base (no-variant) biopolymer is yielded first, followed by combinations in increasing size.
        /// </summary>
        /// <typeparam name="TBioPolymerType">A biopolymer type that supports sequence variants.</typeparam>
        /// <param name="baseBioPolymer">The starting biopolymer (no variations applied).</param>
        /// <param name="variations">Candidate variations to combine.</param>
        /// <param name="maxCombinations">Maximum number of variants (including the base) to yield to bound combinatorial growth.</param>
        /// <returns>An enumerable of applied-variant biopolymers.</returns>
        public static IEnumerable<TBioPolymerType> ApplyAllVariantCombinations<TBioPolymerType>(
            TBioPolymerType baseBioPolymer,
            List<SequenceVariation> variations,
            int maxCombinations)
            where TBioPolymerType : IHasSequenceVariants
        {
            int count = 0; // number of variants yielded so far

            yield return baseBioPolymer;
            count++;
            if (count >= maxCombinations)
                yield break;

            int n = variations.Count; // total variation count
            for (int size = 1; size <= n; size++)
            {
                foreach (var combo in GetCombinations(variations, size))
                {
                    var result = baseBioPolymer; // start from base and apply this combination in order
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
        /// Generates all k-combinations (order-independent, no repetition) of the given list.
        /// This is a standard lexicographic index-based combinator that yields increasing index tuples.
        /// </summary>
        /// <param name="variations">Source list to combine.</param>
        /// <param name="size">Combination size k (0 &lt; k &lt;= n).</param>
        /// <returns>An enumerable of read-only lists containing the selected variations.</returns>
        private static IEnumerable<IList<SequenceVariation>> GetCombinations(List<SequenceVariation> variations, int size)
        {
            int n = variations.Count;
            var indices = new int[size];
            for (int i = 0; i < size; i++) indices[i] = i; // initial 0..k-1

            while (true)
            {
                // Materialize current combination
                var combo = new List<SequenceVariation>(size);
                for (int i = 0; i < size; i++)
                    combo.Add(variations[indices[i]]);
                yield return combo;

                // Advance to next lexicographic combination
                int pos = size - 1;
                while (pos >= 0 && indices[pos] == n - size + pos)
                    pos--;
                if (pos < 0) break;
                indices[pos]++;
                for (int i = pos + 1; i < size; i++)
                    indices[i] = indices[i - 1] + 1;
            }
        }

        /// <summary>
        /// Scans localized modifications for "nucleotide substitution" annotations of the form "X-&gt;Y" and converts them
        /// into concrete <see cref="SequenceVariation"/> objects at the associated positions.
        /// The originating annotation-style modifications are removed from both the possible-localized and original-non-variant collections
        /// (and from consensus mirrors) when they are fully consumed by the conversion.
        /// </summary>
        /// <typeparam name="TBioPolymerType">A biopolymer type that supports sequence variants.</typeparam>
        /// <param name="protein">The biopolymer whose modifications/variants are to be updated in-place.</param>
        public static void ConvertNucleotideSubstitutionModificationsToSequenceVariants<TBioPolymerType>(this TBioPolymerType protein)
            where TBioPolymerType : IHasSequenceVariants
        {
            // Collect mods to remove after converting them to sequence variants
            List<KeyValuePair<int, Modification>> modificationsToRemove = new();

            foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
            {
                foreach (Modification mod in kvp.Value)
                {
                    // Look for annotation-style nucleotide substitutions (e.g., "A->G")
                    if (mod.ModificationType.Contains("nucleotide substitution") && mod.OriginalId.Contains("->"))
                    {
                        string[] originalAndSubstitutedAminoAcids = mod.OriginalId.Split(new[] { "->" }, StringSplitOptions.RemoveEmptyEntries);
                        SequenceVariation sequenceVariation = new SequenceVariation(kvp.Key, kvp.Key, originalAndSubstitutedAminoAcids[0], originalAndSubstitutedAminoAcids[1], "Putative GPTMD Substitution");
                        if (!protein.SequenceVariations.Contains(sequenceVariation))
                        {
                            protein.SequenceVariations.Add(sequenceVariation);
                        }
                        // Defer removal to avoid mutating collection during enumeration
                        KeyValuePair<int, Modification> pair = new(kvp.Key, mod);
                        modificationsToRemove.Add(pair);
                    }
                }
            }

            // Remove the consumed annotation-style modifications from both live and consensus dictionaries
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