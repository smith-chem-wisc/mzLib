using MzLibUtil;
using System.Collections.Generic;


namespace FlashLFQ.IsoTracker
{
    /// <summary>
    /// This class is used to group isobaric peptides together based on their mass.
    /// </summary>
    public class PeptideMassBin
    {
        public double MaxMass { get; private set; } // The mass of the isobaric peptide
        public List<Identification> Ids { get; set; } // The identification of the isobaric peptide

        public PeptideMassBin(Identification id, Tolerance tolerance)
        {
            Ids = new List<Identification>() { id };
            MaxMass = tolerance.GetMaximumValue(id.PeakfindingMass);
        }

        public override bool Equals(object obj)
        {
            return base.Equals(obj);
        }
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
    }
}
