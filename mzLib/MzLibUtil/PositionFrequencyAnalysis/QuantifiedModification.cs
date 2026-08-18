namespace MzLibUtil.PositionFrequencyAnalysis
{
    /// <summary>
    /// A class to store information about a quantified modification.
    /// </summary>
    public class QuantifiedModification
    {
        public string Name { get; set; }
        public int PeptidePositionZeroIsNTerminus { get; set; }
        public int ProteinPositionZeroIsNTerminus { get; set; }
        public double Intensity { get; set; }

        /// <summary>
        /// Constructor for a QuantifiedModification object.
        /// </summary>
        /// <param name="name">Full name of the modification, in the format "MODTYPE: MODID on MOTIF" </param>
        /// <param name="positionInPeptide">Zero-based position in the peptide.</param>
        /// <param name="positionInProtein">Zero-based position in the peptide's parent protein.</param>
        /// <param name="intensity"></param>
        public QuantifiedModification(string name, int positionInPeptide, int? positionInProtein = null, double intensity = 0)
        {
            Name = name;
            PeptidePositionZeroIsNTerminus = positionInPeptide;
            ProteinPositionZeroIsNTerminus = positionInProtein ?? -1; // -1 means that the position in the protein is unknown
            Intensity = intensity;
        }
    }
}
