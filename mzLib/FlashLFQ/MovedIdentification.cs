#nullable enable 
using System;
using System.Collections.Generic;

namespace FlashLFQ
{
    /// <summary>
    /// A peuso-identification object that stores the information of a peptide identification, use for peak tracking
    /// The only difference is that it has a time shift from alignment
    /// </summary>
    internal class MovedIdentification : Identification, IEquatable<Identification>
    {
        /// <summary>
        /// The time shift from alignment
        /// ex. if the time shift is 0.5, the identification is 0.5 minutes later than the original identification
        /// </summary>
        public readonly double TimeShift;
        public MovedIdentification(SpectraFileInfo fileInfo, string BaseSequence, string ModifiedSequence,
            double monoisotopicMass,
            double ms2RetentionTimeInMinutes, int chargeState, List<ProteinGroup> proteinGroups, double timeShift) : base(fileInfo, BaseSequence, ModifiedSequence,
            monoisotopicMass,ms2RetentionTimeInMinutes + timeShift, chargeState, proteinGroups ) 
        {
            
            this.TimeShift = timeShift;
        }

        /// <summary>
        /// Check moved identification is equal to the original identification
        /// Equalities are based on the base sequence, modified sequence, and use for protein quant
        /// </summary>
        /// <param name="id"></param>
        /// <returns></returns>
        public bool Equals(Identification? id)
        {
            if (id is null) return false;
            if (ReferenceEquals(this, id)) return true;
            if (BaseSequence != id.BaseSequence || 
                ModifiedSequence != id.ModifiedSequence || 
                UseForProteinQuant != id.UseForProteinQuant)
            {
                return false;
            }
            return true;
        }
    }
}
