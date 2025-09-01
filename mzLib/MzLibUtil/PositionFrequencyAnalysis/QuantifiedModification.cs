namespace MzLibUtil.PositionFrequencyAnalysis
{
    public class QuantifiedModification
    {
        public string IdWithMotif { get; set; }
        public string ModificationLocalization { get; set; } // e.g. "N-terminus", "C-terminus", or amino acid name
        public int PeptidePositionZeroIsNTerminus { get; set; }
        public int ProteinPositionZeroIsNTerminus { get; set; }
        public double Intensity { get; set; }

        public QuantifiedModification(string idWithMotif, int positionInPeptide, int? positionInProtein = null, string modLocalization = null, double intensity = 0)
        {
            IdWithMotif = idWithMotif;
            PeptidePositionZeroIsNTerminus = positionInPeptide;
            ProteinPositionZeroIsNTerminus = positionInProtein ?? -1; // -1 means that the position in the protein is unknown
            ModificationLocalization = modLocalization ?? "Unknown"; 
            Intensity = intensity;
        }
    }
}
