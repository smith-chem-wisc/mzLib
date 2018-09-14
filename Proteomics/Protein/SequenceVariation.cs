namespace Proteomics
{
    public class SequenceVariation
    {
        /// <summary>
        /// For longer sequence variations, where a range of sequence is replaced. Point mutations should be specified with the same begin and end positions.
        /// </summary>
        /// <param name="OneBasedBeginPosition"></param>
        /// <param name="OneBasedEndPosition"></param>
        /// <param name="OriginalSequence"></param>
        /// <param name="VariantSequence"></param>
        public SequenceVariation(int OneBasedBeginPosition, int OneBasedEndPosition, string OriginalSequence, string VariantSequence, string Description)
        {
            this.OneBasedBeginPosition = OneBasedBeginPosition;
            this.OneBasedEndPosition = OneBasedEndPosition;
            this.OriginalSequence = OriginalSequence ?? "";
            this.VariantSequence = VariantSequence ?? "";
            this.Description = Description ?? "";
        }

        /// <summary>
        /// For variations with only position information (not begin and end).
        /// Sets the end to the end of the original protein sequence to which this variation applies.
        /// </summary>
        /// <param name="OneBasedPosition"></param>
        /// <param name="OriginalSequence"></param>
        /// <param name="VariantSequence"></param>
        /// <param name="Description"></param>
        public SequenceVariation(int OneBasedPosition, string OriginalSequence, string VariantSequence, string Description)
            : this(OneBasedPosition, OneBasedPosition + OriginalSequence.Length - 1, OriginalSequence, VariantSequence, Description)
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
        public string Description { get; }

        public override bool Equals(object obj)
        {
            SequenceVariation s = obj as SequenceVariation;
            return s != null
                && OneBasedBeginPosition == s.OneBasedBeginPosition
                && OneBasedEndPosition == s.OneBasedEndPosition
                && OriginalSequence == s.OriginalSequence
                && VariantSequence == s.VariantSequence
                && Description == s.Description;
        }

        public override int GetHashCode()
        {
            return OneBasedBeginPosition.GetHashCode()
                ^ OneBasedEndPosition.GetHashCode()
                ^ OriginalSequence.GetHashCode()
                ^ VariantSequence.GetHashCode()
                ^ Description.GetHashCode();
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
        internal bool Intersects(ProteolysisProduct segment)
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

        /// <summary>
        /// Determines whether this interval includes the queried interval
        /// </summary>
        /// <param name="segment"></param>
        /// <returns></returns>
        internal bool Includes(ProteolysisProduct segment)
        {
            return OneBasedBeginPosition <= segment.OneBasedBeginPosition && OneBasedEndPosition >= segment.OneBasedEndPosition;
        }

        /// <summary>
        /// Determines whether this interval overlaps the queried position
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        internal bool Includes(int pos)
        {
            return OneBasedBeginPosition <= pos && pos <= OneBasedEndPosition;
        }
    }
}