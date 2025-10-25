using System.Collections.Generic;
using System.Linq;
using Omics.Modifications;

namespace Omics.BioPolymer
{
    /// <summary>
    /// Represents a contiguous sequence variation on a 1-based, inclusive coordinate system.
    /// A variation spans <see cref="OneBasedBeginPosition"/>.. <see cref="OneBasedEndPosition"/> in the parent sequence,
    /// replaces the <see cref="OriginalSequence"/> (which may be empty for an insertion) with <see cref="VariantSequence"/>
    /// (which may be empty for a deletion), and can optionally carry site-specific <see cref="OneBasedModifications"/>.
    /// When available, a parsed <see cref="VariantCallFormat"/> is attached via <see cref="VariantCallFormatDataString"/>.
    /// 
    /// Typical interpretations:
    /// - Substitution: non-empty <see cref="OriginalSequence"/> and non-empty <see cref="VariantSequence"/> of equal length.
    /// - Insertion: empty <see cref="OriginalSequence"/> and non-empty <see cref="VariantSequence"/>.
    /// - Deletion: non-empty <see cref="OriginalSequence"/> and empty <see cref="VariantSequence"/>.
    /// </summary>
    public class SequenceVariation
    {
        /// <summary>
        /// Create a variation with an explicit VCF object.
        /// </summary>
        /// <param name="oneBasedBeginPosition">
        /// 1-based, inclusive start position in the parent sequence where the variation begins.
        /// Must be &gt;= 1. See <see cref="AreValid"/> for validity conditions.
        /// </param>
        /// <param name="oneBasedEndPosition">
        /// 1-based, inclusive end position in the parent sequence where the variation ends.
        /// Must satisfy <c>oneBasedEndPosition &gt;= oneBasedBeginPosition</c>.
        /// </param>
        /// <param name="originalSequence">
        /// Reference subsequence being replaced. Null is coerced to an empty string.
        /// Empty string typically indicates an insertion at <paramref name="oneBasedBeginPosition"/>.
        /// </param>
        /// <param name="variantSequence">
        /// Alternate subsequence to insert in place of <paramref name="originalSequence"/>.
        /// Null is coerced to an empty string. Empty string typically indicates a deletion.
        /// </param>
        /// <param name="description">
        /// Free-form description of the variation. Often the original VCF line or human-readable note.
        /// </param>
        /// <param name="variantCallFormat">
        /// Parsed VCF wrapper for the originating record. Used for downstream analysis of genotype/allele metadata.
        /// </param>
        /// <param name="oneBasedModifications">
        /// Optional mapping from absolute 1-based residue positions to one or more <see cref="Modification"/> objects
        /// applied at that position (in the same coordinate system as <paramref name="oneBasedBeginPosition"/>/<paramref name="oneBasedEndPosition"/>).
        /// If null, an empty dictionary is created.
        /// </param>
        public SequenceVariation(int oneBasedBeginPosition, int oneBasedEndPosition, string originalSequence, string variantSequence, string description, VariantCallFormat variantCallFormat, Dictionary<int, List<Modification>>? oneBasedModifications = null)
        {
            OneBasedBeginPosition = oneBasedBeginPosition;
            OneBasedEndPosition = oneBasedEndPosition;
            OriginalSequence = originalSequence ?? "";
            VariantSequence = variantSequence ?? "";
            Description = description;
            VariantCallFormatDataString = variantCallFormat;
            OneBasedModifications = oneBasedModifications ?? new Dictionary<int, List<Modification>>();
        }

        /// <summary>
        /// Create a variation by providing a raw VCF line (string representation) which will be parsed into a <see cref="VariantCallFormat"/>.
        /// </summary>
        /// <param name="oneBasedBeginPosition">
        /// 1-based, inclusive start position in the parent sequence where the variation begins.
        /// Must be &gt;= 1. See <see cref="AreValid"/> for validity conditions.
        /// </param>
        /// <param name="oneBasedEndPosition">
        /// 1-based, inclusive end position in the parent sequence where the variation ends.
        /// Must satisfy <c>oneBasedEndPosition &gt;= oneBasedBeginPosition</c>.
        /// </param>
        /// <param name="originalSequence">
        /// Reference subsequence being replaced. Null is coerced to an empty string.
        /// Empty string typically indicates an insertion at <paramref name="oneBasedBeginPosition"/>.
        /// </param>
        /// <param name="variantSequence">
        /// Alternate subsequence to insert in place of <paramref name="originalSequence"/>.
        /// Null is coerced to an empty string. Empty string typically indicates a deletion.
        /// </param>
        /// <param name="description">
        /// Free-form description of the variation. Often the original VCF line or human-readable note.
        /// </param>
        /// <param name="variantCallFormatStringRepresentation">
        /// Raw VCF record (a single, tab-delimited line). It is parsed into <see cref="VariantCallFormatDataString"/>.
        /// </param>
        /// <param name="oneBasedModifications">
        /// Optional mapping from absolute 1-based residue positions to one or more <see cref="Modification"/> objects
        /// applied at that position (in the same coordinate system as <paramref name="oneBasedBeginPosition"/>/<paramref name="oneBasedEndPosition"/>).
        /// If null, an empty dictionary is created.
        /// </param>
        public SequenceVariation(int oneBasedBeginPosition, int oneBasedEndPosition, string originalSequence, string variantSequence, string description, string variantCallFormatStringRepresentation, Dictionary<int, List<Modification>>? oneBasedModifications = null)
        {
            OneBasedBeginPosition = oneBasedBeginPosition;
            OneBasedEndPosition = oneBasedEndPosition;
            OriginalSequence = originalSequence ?? "";
            VariantSequence = variantSequence ?? "";
            Description = description;
            VariantCallFormatDataString = new VariantCallFormat(variantCallFormatStringRepresentation);
            OneBasedModifications = oneBasedModifications ?? new Dictionary<int, List<Modification>>();
        }

        /// <summary>
        /// Create a variation without a separate VCF string; a <see cref="VariantCallFormat"/> is still constructed
        /// from <paramref name="description"/> to maintain a non-null object for tests and downstream consumers.
        /// </summary>
        /// <param name="oneBasedBeginPosition">
        /// 1-based, inclusive start position in the parent sequence where the variation begins.
        /// Must be &gt;= 1. See <see cref="AreValid"/> for validity conditions.
        /// </param>
        /// <param name="oneBasedEndPosition">
        /// 1-based, inclusive end position in the parent sequence where the variation ends.
        /// Must satisfy <c>oneBasedEndPosition &gt;= oneBasedBeginPosition</c>.
        /// </param>
        /// <param name="originalSequence">
        /// Reference subsequence being replaced. Null is coerced to an empty string.
        /// Empty string typically indicates an insertion at <paramref name="oneBasedBeginPosition"/>.
        /// </param>
        /// <param name="variantSequence">
        /// Alternate subsequence to insert in place of <paramref name="originalSequence"/>.
        /// Null is coerced to an empty string. Empty string typically indicates a deletion.
        /// </param>
        /// <param name="description">
        /// Free-form description of the variation. Also used to initialize <see cref="VariantCallFormatDataString"/>.
        /// </param>
        /// <param name="oneBasedModifications">
        /// Optional mapping from absolute 1-based residue positions to one or more <see cref="Modification"/> objects
        /// applied at that position (in the same coordinate system as <paramref name="oneBasedBeginPosition"/>/<paramref name="oneBasedEndPosition"/>).
        /// If null, an empty dictionary is created.
        /// </param>
        public SequenceVariation(int oneBasedBeginPosition, int oneBasedEndPosition, string originalSequence, string variantSequence, string description, Dictionary<int, List<Modification>>? oneBasedModifications = null)
        {
            OneBasedBeginPosition = oneBasedBeginPosition;
            OneBasedEndPosition = oneBasedEndPosition;
            OriginalSequence = originalSequence ?? "";
            VariantSequence = variantSequence ?? "";
            Description = description;
            // Always construct a VariantCallFormat so tests relying on non-null VCF objects pass.
            VariantCallFormatDataString = new VariantCallFormat(description);
            OneBasedModifications = oneBasedModifications ?? new Dictionary<int, List<Modification>>();
        }

        /// <summary>
        /// Convenience constructor for single-position variations. The end position is inferred from
        /// <paramref name="oneBasedPosition"/> and the length of <paramref name="originalSequence"/>:
        /// <c>end = position + length(originalSequence) - 1</c> (or <c>end = position</c> if <paramref name="originalSequence"/> is null).
        /// </summary>
        /// <param name="oneBasedPosition">
        /// 1-based, inclusive position of the variation start. Used to infer <see cref="OneBasedEndPosition"/>.
        /// </param>
        /// <param name="originalSequence">
        /// Reference subsequence being replaced. If null, treated as empty for end-position inference.
        /// </param>
        /// <param name="variantSequence">
        /// Alternate subsequence to insert in place of <paramref name="originalSequence"/>. May be empty for deletions.
        /// </param>
        /// <param name="description">
        /// Free-form description of the variation. Also used to initialize <see cref="VariantCallFormatDataString"/>.
        /// </param>
        /// <param name="oneBasedModifications">
        /// Optional mapping from absolute 1-based residue positions to one or more <see cref="Modification"/> objects
        /// applied at that position. If null, an empty dictionary is created.
        /// </param>
        public SequenceVariation(int oneBasedPosition, string originalSequence, string variantSequence, string description, Dictionary<int, List<Modification>>? oneBasedModifications = null)
            : this(oneBasedPosition, originalSequence == null ? oneBasedPosition : oneBasedPosition + originalSequence.Length - 1, originalSequence, variantSequence, description, oneBasedModifications)
        { }

        /// <summary>
        /// 1-based, inclusive start position of this variation within the parent sequence.
        /// </summary>
        public int OneBasedBeginPosition { get; }

        /// <summary>
        /// 1-based, inclusive end position of this variation within the parent sequence.
        /// For single-site variations with non-null <see cref="OriginalSequence"/>, this is
        /// <c>OneBasedBeginPosition + OriginalSequence.Length - 1</c>.
        /// </summary>
        public int OneBasedEndPosition { get; }

        /// <summary>
        /// The reference subsequence replaced by this variation. Empty string implies an insertion.
        /// </summary>
        public string OriginalSequence { get; }

        /// <summary>
        /// The alternate subsequence inserted by this variation. Empty string implies a deletion.
        /// </summary>
        public string VariantSequence { get; }

        /// <summary>
        /// Free-form description of the variation. Often the raw VCF line or a human-readable summary.
        /// </summary>
        public string Description { get; }

        /// <summary>
        /// Optional parsed VCF wrapper providing structured access to the originating VCF record.
        /// May be null in some construction paths; in the provided constructors it is initialized.
        /// </summary>
        public VariantCallFormat? VariantCallFormatDataString { get; }

        /// <summary>
        /// Mapping from absolute 1-based residue positions to a list of <see cref="Modification"/> objects
        /// to apply at each position. Never null; defaults to an empty dictionary.
        /// </summary>
        public Dictionary<int, List<Modification>> OneBasedModifications { get; }

        /// <summary>
        /// Determines value equality with another object.
        /// Two <see cref="SequenceVariation"/> objects are equal when:
        /// - Begin and end positions are equal
        /// - Original and variant sequences are equal (nulls treated as equal only if both are null)
        /// - <see cref="VariantCallFormatDataString"/> are both null or equal
        /// - <see cref="OneBasedModifications"/> have identical key sets and identical flattened modification lists (sequence-equal)
        /// </summary>
        /// <param name="obj">Object to compare against.</param>
        /// <returns>True if equal by the criteria above; otherwise false.</returns>
        public override bool Equals(object obj)
        {
            SequenceVariation s = obj as SequenceVariation;
            return s != null
                && OneBasedBeginPosition == s.OneBasedBeginPosition
                && OneBasedEndPosition == s.OneBasedEndPosition
                && (s.OriginalSequence == null && OriginalSequence == null || OriginalSequence.Equals(s.OriginalSequence))
                && (s.VariantSequence == null && VariantSequence == null || VariantSequence.Equals(s.VariantSequence))
                && ((s.VariantCallFormatDataString == null && VariantCallFormatDataString == null)
                    || (VariantCallFormatDataString != null && VariantCallFormatDataString.Equals(s.VariantCallFormatDataString)))
                && (s.OneBasedModifications == null && OneBasedModifications == null ||
                    s.OneBasedModifications.Keys.ToList().SequenceEqual(OneBasedModifications.Keys.ToList())
                    && s.OneBasedModifications.Values.SelectMany(m => m).ToList().SequenceEqual(OneBasedModifications.Values.SelectMany(m => m).ToList()));
        }

        /// <summary>
        /// Computes a hash code from begin/end positions, sequences, and the VCF wrapper (if present).
        /// </summary>
        /// <returns>A hash code suitable for hash-based collections.</returns>
        public override int GetHashCode()
        {
            return OneBasedBeginPosition.GetHashCode()
                ^ OneBasedEndPosition.GetHashCode()
                ^ OriginalSequence.GetHashCode()
                ^ VariantSequence.GetHashCode()
                ^ (VariantCallFormatDataString?.GetHashCode() ?? 0);
        }

        /// <summary>
        /// Produces a compact, human-readable representation: <c>{OriginalSequence}{OneBasedBeginPosition}{VariantSequence}</c>.
        /// Example: substitution A-&gt;T at position 12 yields <c>"A12T"</c>.
        /// </summary>
        /// <returns>The compact representation string.</returns>
        public string SimpleString()
        {
            return OriginalSequence + OneBasedBeginPosition.ToString() + VariantSequence;
        }

        /// <summary>
        /// Tests whether the current variation intersects (overlaps) another variation in coordinate space.
        /// Intersection is inclusive: any shared position in the 1-based, inclusive ranges is considered overlap.
        /// </summary>
        /// <param name="segment">The other <see cref="SequenceVariation"/> to test.</param>
        /// <returns>True if the ranges overlap; otherwise false.</returns>
        internal bool Intersects(SequenceVariation segment)
        {
            return segment.OneBasedEndPosition >= OneBasedBeginPosition && segment.OneBasedBeginPosition <= OneBasedEndPosition;
        }

        /// <summary>
        /// Tests whether the current variation intersects (overlaps) a truncation product range.
        /// Intersection is inclusive: any shared position in the 1-based, inclusive ranges is considered overlap.
        /// </summary>
        /// <param name="segment">The <see cref="TruncationProduct"/> segment to test.</param>
        /// <returns>True if the ranges overlap; otherwise false.</returns>
        internal bool Intersects(TruncationProduct segment)
        {
            return segment.OneBasedEndPosition >= OneBasedBeginPosition && segment.OneBasedBeginPosition <= OneBasedEndPosition;
        }

        /// <summary>
        /// Tests whether the current variation intersects a single 1-based position.
        /// </summary>
        /// <param name="pos">A 1-based, inclusive position in the parent sequence.</param>
        /// <returns>True if <paramref name="pos"/> lies within the variation’s range; otherwise false.</returns>
        internal bool Intersects(int pos)
        {
            return OneBasedBeginPosition <= pos && pos <= OneBasedEndPosition;
        }

        /// <summary>
        /// Tests whether the current variation fully includes another variation’s range.
        /// Inclusion is inclusive on both ends.
        /// </summary>
        /// <param name="segment">The other <see cref="SequenceVariation"/> to test.</param>
        /// <returns>True if the current range fully contains <paramref name="segment"/>; otherwise false.</returns>
        internal bool Includes(SequenceVariation segment)
        {
            return OneBasedBeginPosition <= segment.OneBasedBeginPosition && OneBasedEndPosition >= segment.OneBasedEndPosition;
        }

        /// <summary>
        /// Tests whether the current variation includes a single 1-based position.
        /// Inclusion is inclusive on both ends.
        /// </summary>
        /// <param name="pos">A 1-based, inclusive position in the parent sequence.</param>
        /// <returns>True if <paramref name="pos"/> lies within the variation’s range; otherwise false.</returns>
        internal bool Includes(int pos)
        {
            return OneBasedBeginPosition <= pos && pos <= OneBasedEndPosition;
        }

        /// <summary>
        /// Validates this variation.
        /// Rules:
        /// 1. Coordinates must be sensible (begin >= 1 and end >= begin).
        /// 2. Variation must represent a meaningful change:
        ///    - Either the sequence actually changes (insertion, deletion, substitution, stop, frameshift),
        ///    - OR there are variant-specific modifications.
        ///    A “no-op” (OriginalSequence == VariantSequence with no variant-specific mods) is invalid.
        /// 3. If variant-specific modifications exist, they must not violate positional constraints
        ///    (see <see cref="GetInvalidModificationPositions"/>).
        /// </summary>
        public bool AreValid()
        {
            if (OneBasedBeginPosition <= 0 || OneBasedEndPosition < OneBasedBeginPosition)
            {
                return false;
            }

            bool noSequenceChange = string.Equals(OriginalSequence ?? string.Empty,
                                                  VariantSequence ?? string.Empty,
                                                  StringComparison.Ordinal);

            bool hasMods = OneBasedModifications != null && OneBasedModifications.Count > 0;

            if (noSequenceChange && !hasMods)
            {
                return false;
            }

            if (!hasMods)
            {
                return true;
            }

            return !GetInvalidModificationPositions().Any();
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
    }
}