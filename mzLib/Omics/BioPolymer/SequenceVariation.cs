using System.Collections.Generic;
using System.Linq;
using Omics.Modifications;

namespace Omics.BioPolymer
{
    /// <summary>
    /// Represents a contiguous amino-acid sequence change (substitution, insertion, deletion, truncation, etc.).
    /// Coordinates are 1-based and inclusive. For point substitutions, begin == end.
    /// <para>
    /// Optional <see cref="VariantCallFormatData"/> (multi‑sample VCF line) can describe the genomic origin,
    /// allelic depth, genotypes, etc. Variant-specific PTMs can be attached via <see cref="OneBasedModifications"/>.
    /// </para>
    /// Validation ensures coordinates are logical and that any supplied variant‑specific modifications
    /// still fall within the valid residue span after the variation is applied (e.g. a premature stop “*”
    /// or a deletion invalidates modifications at and after the replaced region).
    /// </summary>
    public class SequenceVariation
    {
        /// <summary>
        /// Create a sequence variation replacing the span [oneBasedBeginPosition, oneBasedEndPosition]
        /// with <paramref name="variantSequence"/>. The <paramref name="originalSequence"/> is optional
        /// (empty string treated as unknown). A VCF line string may be supplied to initialize
        /// <see cref="VariantCallFormatData"/>. Variant-specific modifications can be provided keyed by
        /// 1-based residue position (post-variation coordinates).
        /// </summary>
        public SequenceVariation(int oneBasedBeginPosition,
            int oneBasedEndPosition,
            string originalSequence,
            string variantSequence,
            string description,
            string? variantCallFormatDataString = null,
            Dictionary<int, List<Modification>>? oneBasedModifications = null)
        {
            OneBasedBeginPosition = oneBasedBeginPosition;
            OneBasedEndPosition = oneBasedEndPosition;
            OriginalSequence = originalSequence ?? "";
            VariantSequence = variantSequence ?? "";
            Description = description;
            VariantCallFormatData = variantCallFormatDataString is null ? null : new VariantCallFormat(variantCallFormatDataString);
            OneBasedModifications = oneBasedModifications ?? new Dictionary<int, List<Modification>>();

            var invalid = GetInvalidModificationPositions().ToList();
            if (invalid.Count > 0)
            {
                throw new ArgumentException($"SequenceVariation contains modification positions that are invalid after applying the variation: {string.Join(", ", invalid)}");
            }

            if (!AreValid())
            {
                throw new ArgumentException("SequenceVariation coordinates are invalid.");
            }
        }

        /// <summary>
        /// For variations with only position information (not begin and end).
        /// Sets the end to the end of the original protein sequence to which this variation applies.
        /// </summary>
        /// <param name="oneBasedPosition"></param>
        /// <param name="originalSequence"></param>
        /// <param name="variantSequence"></param>
        /// <param name="description"></param>
        /// <param name="oneBasedModifications"></param>
        public SequenceVariation(int oneBasedBeginPosition,
            int oneBasedEndPosition,
            string originalSequence,
            string variantSequence,
            string description,
            VariantCallFormat vcf,
            Dictionary<int, List<Modification>>? oneBasedModifications = null)
        {
            OneBasedBeginPosition = oneBasedBeginPosition;
            OneBasedEndPosition = oneBasedEndPosition;
            OriginalSequence = originalSequence ?? "";
            VariantSequence = variantSequence ?? "";
            Description = description;
            VariantCallFormatData = vcf;
            OneBasedModifications = oneBasedModifications ?? new Dictionary<int, List<Modification>>();

            var invalid = GetInvalidModificationPositions().ToList();
            if (invalid.Count > 0)
            {
                throw new ArgumentException($"SequenceVariation contains modification positions that are invalid after applying the variation: {string.Join(", ", invalid)}");
            }

            if (!AreValid())
            {
                throw new ArgumentException("SequenceVariation coordinates are invalid.");
            }
        }
        /// <summary>
        /// Convenience constructor when only a single position is provided (point change or insertion).
        /// If <paramref name="originalSequence"/> is null the end position equals the start; otherwise
        /// it spans the length of <paramref name="originalSequence"/>.
        /// </summary>
        public SequenceVariation(int oneBasedPosition,
            string? originalSequence,
            string variantSequence,
            string description,
            string? variantCallFormatDataString = null,
            Dictionary<int, List<Modification>>? oneBasedModifications = null)
            : this(oneBasedPosition,
                originalSequence == null ? oneBasedPosition : oneBasedPosition + originalSequence.Length - 1,
                originalSequence,
                variantSequence,
                description,
                variantCallFormatDataString,
                oneBasedModifications)
        { }
        /// <summary>
        /// Beginning position of original sequence to be replaced
        /// </summary>
        public int OneBasedBeginPosition { get; }

        /// <summary>
        /// End position of original sequence to be replaced
        /// </summary>
        public int OneBasedEndPosition { get; }

        /// <summary>
        /// Original sequence information (optional)
        /// </summary>
        public string OriginalSequence { get; }

        /// <summary>
        /// Variant sequence information (required)
        /// </summary>
        public string VariantSequence { get; }
        /// <summary>Free-form description (may aggregate provenance / sample info).</summary>
        public string Description { get; }
        /// <summary>
        /// Description of this variation (optional)
        /// </summary>
        public VariantCallFormat? VariantCallFormatData { get; }

        public Dictionary<int, List<Modification>> OneBasedModifications { get; }

        /// <summary>
        /// Unified annotation text for free-form searching/classification.
        /// Prefers the raw VCF line if available, otherwise the free-form Description.
        /// </summary>
        public string SearchableAnnotation => VariantCallFormatData?.Description ?? Description ?? string.Empty;

        /// <summary>
        /// Reference allele (REF) convenience passthrough (null if no VCF).
        /// </summary>
        public string? ReferenceAllele => VariantCallFormatData?.ReferenceAlleleString;

        /// <summary>
        /// First (primary) alternate allele convenience passthrough if available.
        /// Returns null if no VCF or ALT not parsable. (Implement inside VariantCallFormat if not already present.)
        /// </summary>
        public string? AlternateAllele => VariantCallFormatData?.AlternateAlleleString; // ensure VariantCallFormat exposes this; if not, remove.
        public override bool Equals(object obj)
        {
            if (obj is not SequenceVariation s)
                return false;

            if (OneBasedBeginPosition != s.OneBasedBeginPosition
                || OneBasedEndPosition != s.OneBasedEndPosition
                || !string.Equals(OriginalSequence, s.OriginalSequence, StringComparison.Ordinal)
                || !string.Equals(VariantSequence, s.VariantSequence, StringComparison.Ordinal))
            {
                return false;
            }

            // VCF metadata
            if (!((VariantCallFormatData?.Equals(s.VariantCallFormatData)) ?? s.VariantCallFormatData == null))
            {
                return false;
            }

            // Modifications (both constructors ensure dictionary is non-null)
            return ModificationDictionariesEqual(OneBasedModifications, s.OneBasedModifications);
        }

        public override int GetHashCode()
        {
            var hash = new HashCode();
            hash.Add(OneBasedBeginPosition);
            hash.Add(OneBasedEndPosition);
            hash.Add(OriginalSequence);
            hash.Add(VariantSequence);
            hash.Add(VariantCallFormatData?.GetHashCode() ?? 0);

            if (OneBasedModifications != null && OneBasedModifications.Count > 0)
            {
                // Stable ordering
                foreach (var site in OneBasedModifications.OrderBy(k => k.Key))
                {
                    var siteHash = new HashCode();
                    siteHash.Add(site.Key);

                    if (site.Value != null && site.Value.Count > 0)
                    {
                        foreach (var key in site.Value
                                     .Select(m => m.IdWithMotif ?? m.OriginalId ?? m.ToString())
                                     .OrderBy(k => k, StringComparer.Ordinal))
                        {
                            siteHash.Add(key);
                        }
                    }

                    hash.Add(siteHash.ToHashCode());
                }
            }

            return hash.ToHashCode();
        }

        /// <summary>
        /// Returns a simple string represantation of this amino acid change
        /// </summary>
        /// <returns></returns>
        public string SimpleString()
        {
            return OriginalSequence + OneBasedBeginPosition.ToString() + VariantSequence;
        }

        /// <summary>
        /// Determines whether this interval overlaps the queried interval
        /// </summary>
        /// <param name="segment"></param>
        /// <returns></returns>
        internal bool Intersects(SequenceVariation segment)
        {
            return segment.OneBasedEndPosition >= OneBasedBeginPosition && segment.OneBasedBeginPosition <= OneBasedEndPosition;
        }

        /// <summary>
        /// Determines whether this interval overlaps the queried interval
        /// </summary>
        /// <param name="segment"></param>
        /// <returns></returns>
        internal bool Intersects(TruncationProduct segment)
        {
            return segment.OneBasedEndPosition >= OneBasedBeginPosition && segment.OneBasedBeginPosition <= OneBasedEndPosition;
        }

        /// <summary>
        /// Determines whether this interval overlaps the queried position
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        internal bool Intersects(int pos)
        {
            return OneBasedBeginPosition <= pos && pos <= OneBasedEndPosition;
        }

        /// <summary>
        /// Determines whether this interval includes the queried interval
        /// </summary>
        /// <param name="segment"></param>
        /// <returns></returns>
        internal bool Includes(SequenceVariation segment)
        {
            return OneBasedBeginPosition <= segment.OneBasedBeginPosition && OneBasedEndPosition >= segment.OneBasedEndPosition;
        }
        // Commented out by AVC on 4/5/23. Unused and untested in current code base,
        // but can't rule out that it could be useful in the future.  
        /// <summary>
        /// Determines whether this interval includes the queried interval
        /// </summary>
        /// <param name="segment"></param>
        /// <returns></returns>
        // internal bool Includes(TruncationProduct segment)
        // {
        //     return OneBasedBeginPosition <= segment.OneBasedBeginPosition && OneBasedEndPosition >= segment.OneBasedEndPosition;
        // }

        /// <summary>
        /// Determines whether this interval overlaps the queried position
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        internal bool Includes(int pos)
        {
            return OneBasedBeginPosition <= pos && pos <= OneBasedEndPosition;
        }
        public bool AreValid()
        {
            return OneBasedBeginPosition > 0 && OneBasedEndPosition >= OneBasedBeginPosition;
        }
        /// <summary>
        /// Order-insensitive multiset comparison of modification dictionaries.
        /// </summary>
        private static bool ModificationDictionariesEqual(
            Dictionary<int, List<Modification>> a,
            Dictionary<int, List<Modification>> b)
        {
            if (ReferenceEquals(a, b))
                return true;
            if (a is null || b is null)
                return false;
            if (a.Count != b.Count)
                return false;

            // Compare position sets
            if (!a.Keys.OrderBy(k => k).SequenceEqual(b.Keys.OrderBy(k => k)))
                return false;

            foreach (var pos in a.Keys)
            {
                var listA = a[pos];
                var listB = b[pos];

                if (listA is null && listB is null)
                    continue;
                if (listA is null || listB is null)
                    return false;
                if (listA.Count != listB.Count)
                    return false;

                // Build frequency maps for multiset compare
                var freqA = listA
                    .GroupBy(m => m.IdWithMotif ?? m.OriginalId ?? m.ToString())
                    .ToDictionary(g => g.Key, g => g.Count(), StringComparer.Ordinal);
                var freqB = listB
                    .GroupBy(m => m.IdWithMotif ?? m.OriginalId ?? m.ToString())
                    .ToDictionary(g => g.Key, g => g.Count(), StringComparer.Ordinal);

                if (freqA.Count != freqB.Count)
                    return false;

                foreach (var kv in freqA)
                {
                    if (!freqB.TryGetValue(kv.Key, out int countB) || countB != kv.Value)
                        return false;
                }
            }

            return true;
        }
        /// <summary>
        /// Yields modification positions deemed invalid under the current edit semantics.
        /// </summary>
        private IEnumerable<int> GetInvalidModificationPositions()
        {
            if (OneBasedModifications == null || OneBasedModifications.Count == 0)
            {
                yield break;
            }

            bool isTermination = VariantSequence == "*" || VariantSequence.Length == 0;

            if (isTermination)
            {
                foreach (var kvp in OneBasedModifications)
                {
                    if (kvp.Key >= OneBasedBeginPosition)
                    {
                        yield return kvp.Key;
                    }
                }
                yield break;
            }

            int newSpanEnd = OneBasedBeginPosition + VariantSequence.Length - 1;

            foreach (var kvp in OneBasedModifications)
            {
                int pos = kvp.Key;
                if (pos <= 0)
                {
                    yield return pos;
                    continue;
                }

                if (pos >= OneBasedBeginPosition
                    && pos <= OneBasedEndPosition
                    && pos > newSpanEnd)
                {
                    yield return pos;
                }
            }
        }
        /// <summary>
        /// Split multi-sample VCF metadata into per-sample <see cref="SequenceVariation"/> objects.
        /// Produces genotype-aware variants (e.g. optionally yields “no-op” for homozygous reference or
        /// both ref+alt for heterozygous). See XML remarks in implementation for decision matrix.
        /// </summary>
        public List<SequenceVariation> SplitPerGenotype(
            int minDepth = 0,
            bool includeReferenceForHeterozygous = false,
            bool emitReferenceForHomozygousRef = false,
            bool skipIfAltIndexMismatch = true)
        {
            var result = new List<SequenceVariation>();

            if (VariantCallFormatData == null ||
                VariantCallFormatData.Genotypes == null ||
                VariantCallFormatData.Genotypes.Count == 0)
            {
                return result;
            }

            string originalVcfLine = VariantCallFormatData.Description;
            string[] vcfFields = originalVcfLine.Split('\t');
            if (vcfFields.Length < 10)
            {
                return result;
            }

            var fixedCols = vcfFields.Take(9).ToArray();
            string format = fixedCols[8];
            string[] formatTokens = format.Split(':');
            int dpIndex = Array.IndexOf(formatTokens, "DP");
            int sampleCount = vcfFields.Length - 9;
            int storedAltIndex = VariantCallFormatData.AlleleIndex; // 1..N alt, 0 ref, -1 unknown

            for (int sampleIdx = 0; sampleIdx < sampleCount; sampleIdx++)
            {
                string sampleKey = sampleIdx.ToString();
                if (!VariantCallFormatData.Genotypes.TryGetValue(sampleKey, out var gtTokens) || gtTokens.Length == 0)
                {
                    continue;
                }

                int depth = 0;
                if (VariantCallFormatData.AlleleDepths != null &&
                    VariantCallFormatData.AlleleDepths.TryGetValue(sampleKey, out var adTokens) &&
                    adTokens != null && adTokens.Length > 0)
                {
                    foreach (var tok in adTokens)
                    {
                        if (tok == "." || string.IsNullOrWhiteSpace(tok)) continue;
                        if (int.TryParse(tok, out int val) && val >= 0) depth += val;
                    }
                }
                else if (dpIndex >= 0)
                {
                    string sampleColumnRaw = vcfFields[9 + sampleIdx];
                    var parts = sampleColumnRaw.Split(':');
                    if (parts.Length == formatTokens.Length &&
                        int.TryParse(parts[dpIndex], out int dpVal) && dpVal >= 0)
                    {
                        depth = dpVal;
                    }
                }
                if (depth < minDepth)
                {
                    continue;
                }

                VariantCallFormat.Zygosity zyg;
                if (!VariantCallFormatData.ZygosityBySample.TryGetValue(sampleKey, out zyg))
                {
                    var called = gtTokens.Where(a => a != ".").Distinct().ToArray();
                    zyg = called.Length == 0 ? VariantCallFormat.Zygosity.Unknown :
                          called.Length == 1 ? VariantCallFormat.Zygosity.Homozygous :
                          VariantCallFormat.Zygosity.Heterozygous;
                }

                var numericAlleles = new List<int>();
                bool parseError = false;
                foreach (var a in gtTokens)
                {
                    if (a == ".") continue;
                    if (int.TryParse(a, out int ai)) numericAlleles.Add(ai); else { parseError = true; break; }
                }
                if (parseError || numericAlleles.Count == 0)
                {
                    continue;
                }

                bool allRef = numericAlleles.All(a => a == 0);
                bool allStoredAlt = storedAltIndex > 0 && numericAlleles.All(a => a == storedAltIndex);
                bool containsDifferentAlt = storedAltIndex > 0 && numericAlleles.Any(a => a > 0 && a != storedAltIndex);
                if (containsDifferentAlt && skipIfAltIndexMismatch)
                {
                    continue;
                }

                string sampleColumn = vcfFields[9 + sampleIdx];
                string singleSampleLine = string.Join("\t", fixedCols) + "\t" + sampleColumn;

                Dictionary<int, List<Modification>> CloneMods()
                {
                    if (OneBasedModifications == null || OneBasedModifications.Count == 0) return null;
                    var clone = new Dictionary<int, List<Modification>>(OneBasedModifications.Count);
                    foreach (var kv in OneBasedModifications)
                        clone[kv.Key] = new List<Modification>(kv.Value);
                    return clone;
                }

                void TryAdd(int begin, int end, string refSeq, string altSeq, string descTag)
                {
                    string annotatedDesc = $"{Description} | Sample={sampleIdx} Zygosity={zyg} Depth={depth} Mode={descTag}";
                    try
                    {
                        var sv = new SequenceVariation(
                            begin,
                            end,
                            refSeq,
                            altSeq,
                            annotatedDesc,
                            singleSampleLine,
                            CloneMods());
                        if (sv.AreValid())
                        {
                            result.Add(sv);
                        }
                    }
                    catch
                    {
                        // ignore invalid candidate
                    }
                }

                if (allRef)
                {
                    if (emitReferenceForHomozygousRef)
                    {
                        TryAdd(OneBasedBeginPosition, OneBasedEndPosition, OriginalSequence, OriginalSequence, "HomozygousRef");
                    }
                }
                else if (allStoredAlt)
                {
                    TryAdd(OneBasedBeginPosition, OneBasedEndPosition, OriginalSequence, VariantSequence, "HomozygousAlt");
                }
                else
                {
                    if (containsDifferentAlt && storedAltIndex > 0 && !skipIfAltIndexMismatch)
                    {
                        TryAdd(OneBasedBeginPosition, OneBasedEndPosition, OriginalSequence, VariantSequence, "MixedAltIndex(StoredAltOnly)");
                    }
                    else
                    {
                        if (includeReferenceForHeterozygous)
                        {
                            TryAdd(OneBasedBeginPosition, OneBasedEndPosition, OriginalSequence, OriginalSequence, "HeterozygousRef");
                        }
                        TryAdd(OneBasedBeginPosition, OneBasedEndPosition, OriginalSequence, VariantSequence, "HeterozygousAlt");
                    }
                }
            }
            return result;
        }
        /// <summary>
        /// Collapse equivalent variations (same coordinates, original sequence, and variant sequence)
        /// into a single representative per unique key.
        /// <para>
        /// Merging rules:
        /// <list type="bullet">
        /// <item><description><b>Keying:</b> (Begin, End, OriginalSequence, VariantSequence).</description></item>
        /// <item><description><b>Modifications:</b> dictionaries are merged; for each position, modification lists are de-duplicated (using <see cref="Modification.Equals(object)"/>).</description></item>
        /// <item><description><b>VariantCallFormatData:</b> one representative (first non-null) is retained. If multiple distinct non-null instances exist, the first is chosen silently.</description></item>
        /// <item><description><b>Description:</b> If a single source → kept verbatim; if multiple sources → a concise aggregate:
        ///   <c>Combined(n): desc1 | desc2 | desc3 (+k more)</c> (showing at most 3 unique descriptions).</description></item>
        /// <item><description><b>Validation:</b> Each merged candidate is constructed and only returned if <see cref="AreValid"/> passes.</description></item>
        /// </list>
        /// Output is deterministically ordered by Begin, End, OriginalSequence, VariantSequence.
        /// </para>
        /// </summary>
        public static List<SequenceVariation> CombineEquivalent(IEnumerable<SequenceVariation> variations)
        {
            var result = new List<SequenceVariation>();
            if (variations == null)
            {
                return result;
            }

            var groups = variations.GroupBy(v => new
            {
                v.OneBasedBeginPosition,
                v.OneBasedEndPosition,
                Orig = v.OriginalSequence ?? "",
                Var = v.VariantSequence ?? ""
            });

            foreach (var g in groups)
            {
                var members = g.ToList();

                var uniqueDescs = members
                    .Select(v => v.Description)
                    .Where(d => !string.IsNullOrWhiteSpace(d))
                    .Distinct()
                    .ToList();

                string description;
                if (uniqueDescs.Count <= 1)
                {
                    description = uniqueDescs.FirstOrDefault() ?? "";
                }
                else
                {
                    const int maxShow = 3;
                    if (uniqueDescs.Count <= maxShow)
                    {
                        description = $"Combined({uniqueDescs.Count}): " + string.Join(" | ", uniqueDescs);
                    }
                    else
                    {
                        int remain = uniqueDescs.Count - maxShow;
                        description = $"Combined({uniqueDescs.Count}): " +
                                      string.Join(" | ", uniqueDescs.Take(maxShow)) +
                                      $" (+{remain} more)";
                    }
                }

                VariantCallFormat? representativeVcf = members
                    .Select(m => m.VariantCallFormatData)
                    .FirstOrDefault(v => v != null);

                Dictionary<int, List<Modification>>? mergedMods = null;
                foreach (var mv in members)
                {
                    if (mv.OneBasedModifications == null || mv.OneBasedModifications.Count == 0)
                    {
                        continue;
                    }

                    mergedMods ??= new Dictionary<int, List<Modification>>();

                    foreach (var kvp in mv.OneBasedModifications)
                    {
                        if (!mergedMods.TryGetValue(kvp.Key, out var existingList))
                        {
                            mergedMods[kvp.Key] = kvp.Value == null
                                ? new List<Modification>()
                                : kvp.Value.Distinct().ToList();
                        }
                        else
                        {
                            if (kvp.Value != null && kvp.Value.Count > 0)
                            {
                                existingList.AddRange(kvp.Value);
                                mergedMods[kvp.Key] = existingList.Distinct().ToList();
                            }
                        }
                    }
                }

                try
                {
                    var combined = representativeVcf == null
                        ? new SequenceVariation(
                            g.Key.OneBasedBeginPosition,
                            g.Key.OneBasedEndPosition,
                            g.Key.Orig,
                            g.Key.Var,
                            description,
                            (string?)null,
                            mergedMods)
                        : new SequenceVariation(
                            g.Key.OneBasedBeginPosition,
                            g.Key.OneBasedEndPosition,
                            g.Key.Orig,
                            g.Key.Var,
                            description,
                            representativeVcf,
                            mergedMods);

                    if (combined.AreValid())
                    {
                        result.Add(combined);
                    }
                }
                catch
                {
                    // skip invalid merged candidate
                }
            }

            return result
                .OrderBy(v => v.OneBasedBeginPosition)
                .ThenBy(v => v.OneBasedEndPosition)
                .ThenBy(v => v.OriginalSequence)
                .ThenBy(v => v.VariantSequence)
                .ToList();
        }
    }
}