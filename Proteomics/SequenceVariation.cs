namespace Proteomics
{
    public class SequenceVariation
    {

        #region Public Properties

        public int? OneBasedBeginPosition { get; set; }

        public int? OneBasedEndPosition { get; set; }

        public int OneBasedPosition { get; set; }

        public string OriginalSequence { get; set; }

        public string VariantSequence { get; set; }

        public string Description { get; set; }

        #endregion Public Properties

        #region Public Constructor

        /// <summary>
        /// For point mutations, position will be used, and begin and end will be null. For longer sequence variations, position will be -1, and the begin and end positions will be used.
        /// </summary>
        /// <param name="OneBasedBeginPosition"></param>
        /// <param name="OneBasedEndPosition"></param>
        /// <param name="OneBasedPosition"></param>
        /// <param name="OriginalSequence"></param>
        /// <param name="VariantSequence"></param>
        public SequenceVariation(int? OneBasedBeginPosition, int? OneBasedEndPosition, int OneBasedPosition, string OriginalSequence, string VariantSequence, string Description)
        {
            this.OneBasedBeginPosition = OneBasedBeginPosition;
            this.OneBasedEndPosition = OneBasedEndPosition;
            this.OneBasedPosition = OneBasedPosition;
            this.OriginalSequence = OriginalSequence;
            this.VariantSequence = VariantSequence;
            this.Description = Description;
        }

        #endregion Public Constructor

    }
}
