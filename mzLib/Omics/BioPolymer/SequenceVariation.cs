using System.Collections.Generic;
using System.Linq;
using Omics.Modifications;

namespace Omics.BioPolymer
{
    public class SequenceVariation
    {
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

        public SequenceVariation(int oneBasedPosition, string originalSequence, string variantSequence, string description, Dictionary<int, List<Modification>>? oneBasedModifications = null)
            : this(oneBasedPosition, originalSequence == null ? oneBasedPosition : oneBasedPosition + originalSequence.Length - 1, originalSequence, variantSequence, description, oneBasedModifications)
        { }

        public int OneBasedBeginPosition { get; }
        public int OneBasedEndPosition { get; }
        public string OriginalSequence { get; }
        public string VariantSequence { get; }
        public string Description { get; }
        public VariantCallFormat? VariantCallFormatDataString { get; }
        public Dictionary<int, List<Modification>> OneBasedModifications { get; }

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

        public override int GetHashCode()
        {
            return OneBasedBeginPosition.GetHashCode()
                ^ OneBasedEndPosition.GetHashCode()
                ^ OriginalSequence.GetHashCode()
                ^ VariantSequence.GetHashCode()
                ^ (VariantCallFormatDataString?.GetHashCode() ?? 0);
        }

        public string SimpleString()
        {
            return OriginalSequence + OneBasedBeginPosition.ToString() + VariantSequence;
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

        public bool AreValid()
        {
            return OneBasedBeginPosition > 0 && OneBasedEndPosition >= OneBasedBeginPosition;
        }
    }
}