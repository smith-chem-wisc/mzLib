using System.Collections.Generic;
using System.Text.RegularExpressions;

namespace MzLibUtil.PositionFrequencyAnalysis
{
    public class QuantifiedPeptideRecord
    {
        public string FullSequence { get; set; }
        public string BaseSequence { get; set; }
        public HashSet<string> ProteinGroups { get; set; }
        public double Intensity { get; set; }
        /// <summary>
        /// A record of a quantified peptide, storing its full sequence (with modifications), base sequence (without modifications),
        /// protein groups it maps to, and intensity. The base sequence is derived from the full sequence and is not passed 
        /// as initialization parameter.
        /// </summary>
        /// <param name="fullSequence"></param>
        /// <param name="proteinGroups"></param>
        /// <param name="intensity"></param>
        public QuantifiedPeptideRecord(string fullSequence, HashSet<string> proteinGroups, double intensity)
        {
            FullSequence = fullSequence;
            ProteinGroups = proteinGroups;
            Intensity = intensity;
            BaseSequence = fullSequence.GetBaseSequenceFromFullSequence();
        }
    }
}
