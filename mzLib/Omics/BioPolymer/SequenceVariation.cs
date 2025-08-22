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
        /// <param name="oneBasedBeginPosition"></param>
        /// <param name="oneBasedEndPosition"></param>
        /// <param name="originalSequence"></param>
        /// <param name="variantSequence"></param>
        /// <param name="oneBasedModifications"></param>
        public SequenceVariation(int oneBasedBeginPosition, int oneBasedEndPosition, string originalSequence, string variantSequence, string description, Dictionary<int, List<Modification>>? oneBasedModifications = null)
        {
            OneBasedBeginPosition = oneBasedBeginPosition;
            OneBasedEndPosition = oneBasedEndPosition;
            OriginalSequence = originalSequence ?? "";
            VariantSequence = variantSequence ?? "";
            Description = new SequenceVariantDescription(description);
            OneBasedModifications = oneBasedModifications ?? [];
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
        public SequenceVariation(int oneBasedPosition, string originalSequence, string variantSequence, string description, Dictionary<int, List<Modification>>? oneBasedModifications = null)
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
        /// Original sequence information (optional)
        /// </summary>
        public string OriginalSequence { get; }

        /// <summary>
        /// Variant sequence information (required)
        /// </summary>
        public string VariantSequence { get; }

        /// <summary>
        /// Description of this variation (optional)
        /// </summary>
        public SequenceVariantDescription Description { get; }

        /// <summary>
        /// Modifications specifically for this variant
        /// </summary>
        public Dictionary<int, List<Modification>> OneBasedModifications { get; }

        public override bool Equals(object obj)
        {
            SequenceVariation s = obj as SequenceVariation;
            return s != null
                && OneBasedBeginPosition == s.OneBasedBeginPosition
                && OneBasedEndPosition == s.OneBasedEndPosition
                && (s.OriginalSequence == null && OriginalSequence == null || OriginalSequence.Equals(s.OriginalSequence))
                && (s.VariantSequence == null && VariantSequence == null || VariantSequence.Equals(s.VariantSequence))
                && (s.Description == null && Description == null || Description.Equals(s.Description))
                && (s.OneBasedModifications == null && OneBasedModifications == null ||
                    s.OneBasedModifications.Keys.ToList().SequenceEqual(OneBasedModifications.Keys.ToList())
                    && s.OneBasedModifications.Values.SelectMany(m => m).ToList().SequenceEqual(OneBasedModifications.Values.SelectMany(m => m).ToList()));
        }

        public override int GetHashCode()
        {
            return OneBasedBeginPosition.GetHashCode()
                ^ OneBasedEndPosition.GetHashCode()
                ^ OriginalSequence.GetHashCode() // null handled in constructor
                ^ VariantSequence.GetHashCode() // null handled in constructor
                ^ Description.GetHashCode(); // always constructed in constructor
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
    }
}