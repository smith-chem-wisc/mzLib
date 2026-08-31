using System;
using System.Collections.Generic;

namespace FlashLFQ.IsoTracker
{
    /// <summary>
    /// Represents a group of peptides that share the same base sequence and monoisotopic mass. Only used for isobaric quantification.
    /// </summary>
    public class IsobaricPeptideGroup
    {
        public string BaseSequence { get; set; }
        /// <summary>
        /// The monoisotopic mass is got by dividing the monoMass with 0.01 in GroupPeptideForIsoTracker.
        /// The number is used for peptide grouping.
        /// </summary>
        public double MonoisotopicMass { get; set; }
        public List<Identification> Identifications { get; set; }
        public IsobaricPeptideGroup(string baseSequence, double monoisotopicMass, List<Identification> identifications)
        {
            BaseSequence = baseSequence;
            MonoisotopicMass = monoisotopicMass;
            Identifications = identifications ?? new List<Identification>();
        }
        public override bool Equals(object obj)
        {
            if (obj is IsobaricPeptideGroup other)
            {
                // Compare BaseSequence and MonoisotopicMass (rounded to 4 decimal places for floating point safety)
                return BaseSequence == other.BaseSequence && MonoisotopicMass == other.MonoisotopicMass;
            }

            return false;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(BaseSequence, MonoisotopicMass);
        }
    }
}
