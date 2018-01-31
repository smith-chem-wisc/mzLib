﻿using System.Collections.Generic;

namespace FlashLFQ
{
    public class Identification
    {
        #region Public Fields

        public readonly string BaseSequence;
        public readonly string FullSequence;
        public readonly double ms2RetentionTime;
        public readonly double monoisotopicMass;
        public int chargeState;
        public List<ProteinGroup> proteinGroups;
        public double massToLookFor;
        public string fileName = "";

        #endregion Public Fields

        #region Public Constructors

        public Identification(string fileName, string BaseSequence, string FullSequence, double monoisotopicMass, double ms2RetentionTime, int chargeState)
        {
            this.fileName = fileName;
            this.BaseSequence = BaseSequence;
            this.FullSequence = FullSequence;
            this.monoisotopicMass = monoisotopicMass;
            this.ms2RetentionTime = ms2RetentionTime;
            this.chargeState = chargeState;
            this.proteinGroups = new List<ProteinGroup>();
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return FullSequence;
        }

        #endregion Public Methods
    }
}