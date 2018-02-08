using System.Collections.Generic;

namespace FlashLFQ
{
    public class Identification
    {
        #region Public Fields

        public readonly string BaseSequence;
        public readonly string ModifiedSequence;
        public readonly double ms2RetentionTimeInMinutes;
        public readonly double monoisotopicMass;
        public int precursorChargeState;
        public readonly List<string> proteinGroupNames;
        public HashSet<ProteinGroup> proteinGroups;
        public double massToLookFor;
        public readonly RawFileInfo fileInfo;

        #endregion Public Fields

        #region Public Constructors

        public Identification(RawFileInfo fileInfo, string BaseSequence, string ModifiedSequence, double monoisotopicMass, double ms2RetentionTimeInMinutes, int chargeState, List<string> proteinGroupNames)
        {
            this.fileInfo = fileInfo;
            this.BaseSequence = BaseSequence;
            this.ModifiedSequence = ModifiedSequence;
            this.monoisotopicMass = monoisotopicMass;
            this.ms2RetentionTimeInMinutes = ms2RetentionTimeInMinutes;
            this.precursorChargeState = chargeState;
            this.proteinGroupNames = proteinGroupNames;
            proteinGroups = new HashSet<ProteinGroup>();
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return ModifiedSequence;
        }

        #endregion Public Methods
    }
}