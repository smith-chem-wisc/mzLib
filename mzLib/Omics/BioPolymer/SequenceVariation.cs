using System;
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
        #region Constructors

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
        /// Overload accepting an already parsed <see cref="VariantCallFormat"/> instance.
        /// </summary>
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

        #endregion

        #region Public Properties

        /// <summary>1-based inclusive begin coordinate.</summary>
        public int OneBasedBeginPosition { get; }

        /// <summary>1-based inclusive end coordinate.</summary>
        public int OneBasedEndPosition { get; }

        /// <summary>Original (replaced) amino acid sequence segment (may be empty for insertions).</summary>
        public string OriginalSequence { get; }

        /// <summary>New amino acid sequence inserted in place of <see cref="OriginalSequence"/> (empty for deletions).</summary>
        public string VariantSequence { get; }

        /// <summary>Free-form description (may aggregate provenance / sample info).</summary>
        public string Description { get; }

        /// <summary>Optional multi-sample VCF record describing the variant (can be null or collapsed).</summary>
        public VariantCallFormat? VariantCallFormatData { get; }

        /// <summary>
        /// Variant-specific modifications keyed by 1-based residue positions in the sequence AFTER variation application.
        /// Positions are validated in <see cref="AreValid"/> against the altered span (<see cref="VariantSequence"/>).
        /// </summary>
        public Dictionary<int, List<Modification>> OneBasedModifications { get; }

        #endregion

        #region Equality / Hash

        public override bool Equals(object obj)
        {
            SequenceVariation s = obj as SequenceVariation;
            return s != null
                && OneBasedBeginPosition == s.OneBasedBeginPosition
                && OneBasedEndPosition == s.OneBasedEndPosition
                && (s.OriginalSequence == null && OriginalSequence == null || OriginalSequence.Equals(s.OriginalSequence))
                && (s.VariantSequence == null && VariantSequence == null || VariantSequence.Equals(s.VariantSequence))
                && ((VariantCallFormatData?.Equals(s.VariantCallFormatData)) ?? s.VariantCallFormatData == null)
                && (s.OneBasedModifications == null && OneBasedModifications == null ||
                    s.OneBasedModifications.Keys.ToList().SequenceEqual(OneBasedModifications.Keys.ToList())
                    && s.OneBasedModifications.Values.SelectMany(m => m).ToList().SequenceEqual(OneBasedModifications.Values.SelectMany(m => m).ToList()));
        }

        public override int GetHashCode()
        {
            return OneBasedBeginPosition.GetHashCode()
                ^ OneBasedEndPosition.GetHashCode()
                ^ OriginalSequence.GetHashCode()
                ^ VariantSequence.GetHashCode()
                ^ (VariantCallFormatData?.GetHashCode() ?? 0);
        }

        #endregion

        #region Convenience / Interval Logic

        /// <summary>Simple concatenated representation (Original + Begin + Variant).</summary>
        public string SimpleString() => OriginalSequence + OneBasedBeginPosition + VariantSequence;

        internal bool Intersects(SequenceVariation segment) =>
            segment.OneBasedEndPosition >= OneBasedBeginPosition && segment.OneBasedBeginPosition <= OneBasedEndPosition;

        internal bool Intersects(TruncationProduct segment) =>
            segment.OneBasedEndPosition >= OneBasedBeginPosition && segment.OneBasedBeginPosition <= OneBasedEndPosition;

        internal bool Intersects(int pos) => OneBasedBeginPosition <= pos && pos <= OneBasedEndPosition;

        internal bool Includes(SequenceVariation segment) =>
            OneBasedBeginPosition <= segment.OneBasedBeginPosition && OneBasedEndPosition >= segment.OneBasedEndPosition;

        internal bool Includes(int pos) => OneBasedBeginPosition <= pos && pos <= OneBasedEndPosition;

        #endregion

        #region Validation

        /// <summary>
        /// Validates coordinate ordering (begin &gt;= 1 and end &gt;= begin) and ensures
        /// that any variant-specific modifications remain addressable after the edit:
        /// <list type="number">
        /// <item>Deletion (VariantSequence length == 0) or termination (“*”): disallow modifications at/after begin.</item>
        /// <item>Otherwise: modifications inside the replaced span must fall within the new substituted span.</item>
        /// </list>
        /// </summary>
        public bool AreValid()
        {
            if (OneBasedBeginPosition <= 0 || OneBasedEndPosition < OneBasedBeginPosition)
            {
                return false;
            }

            if (OneBasedModifications == null || OneBasedModifications.Count == 0)
            {
                return true;
            }

            return !GetInvalidModificationPositions().Any();
        }

        #endregion

        #region Genotype Splitting

        /// <summary>
        /// Split multi-sample VCF metadata into per-sample <see cref="SequenceVariation"/> objects.
        /// Produces genotype-aware variants (e.g. optionally yields “no-op” for homozygous reference or
        /// both ref+alt for heterozygous). See XML remarks in source for decision matrix.
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

                // Depth
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

                // Zygosity
                VariantCallFormat.Zygosity zyg;
                if (!VariantCallFormatData.ZygosityBySample.TryGetValue(sampleKey, out zyg))
                {
                    var called = gtTokens.Where(a => a != ".").Distinct().ToArray();
                    zyg = called.Length == 0 ? VariantCallFormat.Zygosity.Unknown :
                          called.Length == 1 ? VariantCallFormat.Zygosity.Homozygous :
                          VariantCallFormat.Zygosity.Heterozygous;
                }

                // Alleles
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
                        // ignore
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

        #endregion

        #region Combination / Collapsing

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
        /// <param name="variations">Input collection (may be null or empty).</param>
        /// <returns>Collapsed list of <see cref="SequenceVariation"/> objects.</returns>
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

                // Collect distinct descriptions (ignore null/whitespace)
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

                // Choose representative VCF (first non-null)
                VariantCallFormat? representativeVcf = members
                    .Select(m => m.VariantCallFormatData)
                    .FirstOrDefault(v => v != null);

                // Merge modifications
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

                // Construct new merged variation
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
                    // Skip invalid merged candidate
                }
            }

            return result
                .OrderBy(v => v.OneBasedBeginPosition)
                .ThenBy(v => v.OneBasedEndPosition)
                .ThenBy(v => v.OriginalSequence)
                .ThenBy(v => v.VariantSequence)
                .ToList();
        }

        #endregion

        #region Modification Management

        /// <summary>
        /// Attempt to add a single variant-specific modification at the supplied 1-based position
        /// (post-variation coordinate system). Applies the same validity rules enforced during
        /// construction and by <see cref="AreValid"/> / internal <c>GetInvalidModificationPositions</c>.
        /// </summary>
        /// <param name="oneBasedPosition">1-based residue position AFTER applying this variation.</param>
        /// <param name="modification">Modification to add (must be non-null).</param>
        /// <param name="error">
        /// Populated with a short reason when the addition fails; null when successful.
        /// </param>
        /// <returns>true if the modification was added (or was already present at that position); false otherwise.</returns>
        public bool TryAddModification(int oneBasedPosition, Modification modification, out string? error)
        {
            error = null;

            if (modification is null)
            {
                error = "Modification is null.";
                return false;
            }

            if (oneBasedPosition <= 0)
            {
                error = "Position must be > 0.";
                return false;
            }

            bool isTermination = VariantSequence == "*" || VariantSequence.Length == 0;

            if (isTermination)
            {
                // No modifications allowed at or after the variation begin for termination/deletion
                if (oneBasedPosition >= OneBasedBeginPosition)
                {
                    error = "Position invalid for a termination or deletion at/after the begin coordinate.";
                    return false;
                }
            }
            else
            {
                // NEW LOGIC:
                // Only enforce the "beyond new variant span" restriction for coordinates that were actually
                // inside the ORIGINAL replaced span (i.e. <= original end). This allows adding modifications
                // immediately after an insertion expansion, which was previously (incorrectly) rejected.
                // Original replaced span = [OneBasedBeginPosition, OneBasedEndPosition]
                // New variant span        = [OneBasedBeginPosition, OneBasedBeginPosition + VariantSequence.Length - 1]
                int newSpanEnd = OneBasedBeginPosition + VariantSequence.Length - 1;

                if (oneBasedPosition >= OneBasedBeginPosition
                    && oneBasedPosition <= OneBasedEndPosition   // ensure it was in the original replaced region
                    && oneBasedPosition > newSpanEnd)            // but lies past the substituted span
                {
                    error = "Position lies beyond the new variant span inside the edited region.";
                    return false;
                }
            }

            if (!OneBasedModifications.TryGetValue(oneBasedPosition, out var list))
            {
                list = new List<Modification>();
                OneBasedModifications[oneBasedPosition] = list;
            }

            if (!list.Contains(modification))
            {
                list.Add(modification);
            }

            return true;
        }

        /// <summary>
        /// Bulk-add multiple modifications. Each entry is validated with <see cref="TryAddModification"/>.
        /// </summary>
        /// <param name="modifications">
        /// Sequence of (position, modification) pairs (positions are 1-based post-variation).
        /// </param>
        /// <param name="throwOnFirstInvalid">
        /// If true, throws on the first invalid modification (nothing is rolled back).
        /// If false, silently skips invalid entries and records them in <paramref name="skipped"/>.
        /// </param>
        /// <param name="skipped">
        /// Returns a list of (position, reason) pairs for invalid entries when not throwing.
        /// Null when all succeeded or when <paramref name="throwOnFirstInvalid"/> is true and no invalid encountered.
        /// </param>
        /// <returns>The number of successfully added (new or deduplicated) modification positions affected.</returns>
        public int AddModifications(
            IEnumerable<(int position, Modification modification)> modifications,
            bool throwOnFirstInvalid,
            out List<(int position, string reason)>? skipped)
        {
            skipped = null;
            if (modifications == null)
            {
                return 0;
            }

            int affectedPositions = 0;

            foreach (var (pos, mod) in modifications)
            {
                if (TryAddModification(pos, mod, out var reason))
                {
                    affectedPositions++;
                }
                else
                {
                    if (throwOnFirstInvalid)
                    {
                        throw new ArgumentException($"Invalid modification at position {pos}: {reason}");
                    }

                    skipped ??= new List<(int, string)>();
                    skipped.Add((pos, reason ?? "Unknown reason"));
                }
            }

            return affectedPositions;
        }

        #endregion

        #region Internal Helpers

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

                // Updated to match TryAddModification logic: only invalidate when the position was inside
                // the ORIGINAL replaced span but past the substituted (shorter) variant span.
                if (pos >= OneBasedBeginPosition
                    && pos <= OneBasedEndPosition
                    && pos > newSpanEnd)
                {
                    yield return pos;
                }
            }
        }

        #endregion
    }
}