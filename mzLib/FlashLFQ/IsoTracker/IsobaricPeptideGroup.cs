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
                return BaseSequence == other.BaseSequence &&
                       Math.Round(Convert.ToDouble(MonoisotopicMass), 4) == Math.Round(Convert.ToDouble(other.MonoisotopicMass), 4);
            }

            return false;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(BaseSequence, MonoisotopicMass);
        }
    }
}
