using System.Collections.Generic;

namespace MzLibUtil.PositionFrequencyAnalysis
{
    /// <summary>
    /// A lightweight record of a quantified peptide, storing its full sequence (with modifications),
    /// base sequence (without modifications), the protein groups it maps to, and its observed intensity.
    /// The base sequence is derived automatically from the full sequence.
    /// </summary>
    public class QuantifiedPeptideRecord
    {
        public string FullSequence { get; set; }
        public string BaseSequence { get; private set; } 
        public HashSet<string> ProteinGroups { get; set; }
        public double Intensity { get; set; }
        /// <summary>
        /// Initializes a new <see cref="QuantifiedPeptideRecord"/>.
        /// </summary>
        /// <param name="fullSequence">Full peptide sequence with embedded modification notation.</param>
        /// <param name="proteinGroups">Protein groups this peptide maps to.</param>
        /// <param name="intensity">Observed quantification intensity.</param>
        public QuantifiedPeptideRecord(string fullSequence, HashSet<string> proteinGroups, double intensity)
        {
            FullSequence = fullSequence;
            ProteinGroups = proteinGroups;
            Intensity = intensity;
            BaseSequence = fullSequence.GetBaseSequenceFromFullSequence();
        }
    }
}
