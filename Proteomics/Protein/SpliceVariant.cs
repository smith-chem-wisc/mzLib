using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace Proteomics
{
    public class SpliceVariant
    {
        /// <summary>
        /// For variants with sequence elements.
        /// </summary>
        /// <param name="oneBasedBegin"></param>
        /// <param name="oneBasedEnd"></param>
        /// <param name="originalSequence"></param>
        /// <param name="variantSequence"></param>
        /// <param name="description"></param>
        /// <param name="oneBasedModifications"></param>
        public SpliceVariant(int oneBasedBegin, int oneBasedEnd, string originalSequence, string variantSequence, string description, Dictionary<int, List<Modification>> oneBasedModifications = null)
        {
            OneBasedBeginPosition = oneBasedBegin;
            OneBasedEndPosition = oneBasedEnd;
            OriginalSequence = originalSequence ?? "";
            VariantSequence = variantSequence ?? "";
            Description = description;
            OneBasedModifications = oneBasedModifications ?? new Dictionary<int, List<Modification>>();
        }

        /// <summary>
        /// For variants with only position information (not begin and end).
        /// Sets the end to the end of the original protein sequence to which this variation applies.
        /// </summary>
        /// <param name="oneBasedPosition"></param>
        /// <param name="originalSequence"></param>
        /// <param name="variantSequence"></param>
        /// <param name="description"></param>
        /// <param name="oneBasedModifications"></param>
        public SpliceVariant(int oneBasedPosition, string originalSequence, string variantSequence, string description, Dictionary<int, List<Modification>> oneBasedModifications = null)
            : this(oneBasedPosition, originalSequence == null ? oneBasedPosition : oneBasedPosition + originalSequence.Length - 1, originalSequence, variantSequence, description, oneBasedModifications)
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
        /// Original sequence information (optional, represents a missing sequence if absent)
        /// </summary>
        public string OriginalSequence { get; }

        /// <summary>
        /// Variant sequence information (optional, represents a missing sequence if absent)
        /// </summary>
        public string VariantSequence { get; }

        /// <summary>
        /// Description of this variation (optional)
        /// </summary>
        public string Description { get; }

        /// <summary>
        /// Modifications specifically for this variant
        /// </summary>
        public Dictionary<int, List<Modification>> OneBasedModifications { get; }

        public override bool Equals(object obj)
        {
            SpliceVariant s = obj as SpliceVariant;
            return s != null 
                && OneBasedBeginPosition == s.OneBasedBeginPosition
                && OneBasedEndPosition == s.OneBasedEndPosition
                && (OriginalSequence == null && s.OriginalSequence == null || OriginalSequence.Equals(s.OriginalSequence))
                && (VariantSequence == null && s.VariantSequence == null || VariantSequence.Equals(s.VariantSequence))
                && (Description == null && s.Description == null || Description.Equals(s.Description))
                && (OneBasedModifications == null && s.OneBasedModifications == null ||
                    OneBasedModifications.Keys.ToList().SequenceEqual(s.OneBasedModifications.Keys.ToList())
                    && OneBasedModifications.Values.SelectMany(m => m).ToList().SequenceEqual(s.OneBasedModifications.Values.SelectMany(m => m).ToList()));
        }

        public override int GetHashCode()
        {
            return OneBasedBeginPosition.GetHashCode()
                ^ OneBasedEndPosition.GetHashCode()
                ^ OriginalSequence.GetHashCode() // null handled in constructor
                ^ VariantSequence.GetHashCode() // null handled in constructor
                ^ Description.GetHashCode(); // always constructed in constructor
        }
    }
}
