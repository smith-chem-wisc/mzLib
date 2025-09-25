using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Modifications;

namespace Omics.BioPolymer
{
    public class SequenceVariation
    {
        /// <summary>
        /// For longer sequence variations, where a range of sequence is replaced. Point mutations should be specified with the same begin and end positions.
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
        /// Overload that takes an already-parsed VariantCallFormat.
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
        /// For variations with only position information (not begin and end).
        /// Sets the end to the end of the original sequence span this variation replaces.
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

        public int OneBasedBeginPosition { get; }
        public int OneBasedEndPosition { get; }
        public string OriginalSequence { get; }
        public string VariantSequence { get; }
        public string Description { get; }
        public VariantCallFormat? VariantCallFormatData { get; }
        public Dictionary<int, List<Modification>> OneBasedModifications { get; }

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

        public string SimpleString()
        {
            return OriginalSequence + OneBasedBeginPosition + VariantSequence;
        }

        internal bool Intersects(SequenceVariation segment)
        {
            return segment.OneBasedEndPosition >= OneBasedBeginPosition && segment.OneBasedBeginPosition <= OneBasedEndPosition;
        }

        internal bool Intersects(TruncationProduct segment)
        {
            return segment.OneBasedEndPosition >= OneBasedBeginPosition && segment.OneBasedBeginPosition <= OneBasedEndPosition;
        }

        internal bool Intersects(int pos)
        {
            return OneBasedBeginPosition <= pos && pos <= OneBasedEndPosition;
        }

        internal bool Includes(SequenceVariation segment)
        {
            return OneBasedBeginPosition <= segment.OneBasedBeginPosition && OneBasedEndPosition >= segment.OneBasedEndPosition;
        }

        internal bool Includes(int pos)
        {
            return OneBasedBeginPosition <= pos && pos <= OneBasedEndPosition;
        }

        /// <summary>
        /// Validates coordinate logic AND that all modification positions remain valid after applying the variation.
        /// Rules / assumptions:
        /// 1. Coordinates must be positive and ordered.
        /// 2. The region [Begin, End] of the original sequence is replaced by VariantSequence.
        /// 3. If VariantSequence == "*" (termination) OR VariantSequence length == 0 (deletion) then
        ///    no modification at or beyond OneBasedBeginPosition is allowed (the sequence terminates or is removed there).
        /// 4. Otherwise, modifications inside the replaced span must fall within the new span:
        ///       Allowed internal range: [Begin, Begin + VariantSequence.Length - 1]
        ///    Modifications before Begin are always allowed (unchanged prefix).
        ///    (We do not attempt to remap downstream positions here because
        ///     keys are assumed to represent positions in the post-variation sequence.)
        /// </summary>
        public bool AreValid()
        {
            if (OneBasedBeginPosition <= 0 || OneBasedEndPosition < OneBasedBeginPosition)
            {
                return false;
            }

            // If no modifications, coordinate validation above is enough
            if (OneBasedModifications == null || OneBasedModifications.Count == 0)
            {
                return true;
            }

            return !GetInvalidModificationPositions().Any();
        }

        /// <summary>
        /// Returns modification positions that are invalid under the current variation assumptions (see AreValid()).
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
                // Any modification at or after the begin position becomes invalid
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
                // negative or zero always invalid
                if (pos <= 0)
                {
                    yield return pos;
                    continue;
                }

                // Inside replaced region AFTER applying variation must lie in the new span
                if (pos >= OneBasedBeginPosition && pos > newSpanEnd)
                {
                    yield return pos;
                }
            }
        }
    }
}