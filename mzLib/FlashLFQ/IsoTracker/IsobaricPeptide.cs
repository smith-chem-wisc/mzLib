using MzLibUtil;
using System.Collections.Generic;


namespace FlashLFQ.IsoTracker
{
    public class IsobaricPeptide
    {
        public double MaxMass { get; private set; } // The mass of the isobaric peptide
        public List<Identification> Ids { get; set; } // The identification of the isobaric peptide

        public IsobaricPeptide(Identification id, string ppm)
        {
            Ids = new List<Identification>() {id};
            Tolerance tolerance = Tolerance.ParseToleranceString(ppm);
            MaxMass = tolerance.GetMaximumValue(id.PeakfindingMass);
        }
    }
}
