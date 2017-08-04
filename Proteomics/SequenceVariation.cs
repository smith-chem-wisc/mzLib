namespace Proteomics
{
    public class SequenceVariation
    {

        #region Public Properties

        /// <summary>
        /// Beginning position of original sequence to be replaced
        /// </summary>
        public int OneBasedBeginPosition { get; set; }

        /// <summary>
        /// End position of original sequence to be replaced
        /// </summary>
        public int OneBasedEndPosition { get; set; }

        /// <summary>
        /// Original sequence information (optional)
        /// </summary>
        public string OriginalSequence { get; set; }

        /// <summary>
        /// Variant sequence information (required)
        /// </summary>
        public string VariantSequence { get; set; }

        /// <summary>
        /// Description of this variation (optional)
        /// </summary>
        public string Description { get; set; }

        #endregion Public Properties

        #region Public Constructor

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
            this.OriginalSequence = OriginalSequence;
            this.VariantSequence = VariantSequence;
            this.Description = Description;
        }

        /// <summary>
        /// For point mutations. Sets begin and end to the same position.
        /// </summary>
        /// <param name="OneBasedPosition"></param>
        /// <param name="OriginalSequence"></param>
        /// <param name="VariantSequence"></param>
        /// <param name="Description"></param>
        public SequenceVariation(int OneBasedPosition, string OriginalSequence, string VariantSequence, string Description)
            : this(OneBasedPosition, OneBasedPosition, OriginalSequence, VariantSequence, Description)
        { }

        #endregion Public Constructor

        #region Public Methods

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

        #endregion Public Methods

    }
}
