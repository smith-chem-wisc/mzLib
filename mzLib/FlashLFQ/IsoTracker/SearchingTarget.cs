using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ.IsoTracker
{
    public class SearchingTarget
    {
        public List<char> TargetMotifs; // motif string list

        /// <summary>
        /// Directly use the motif string to create a SearchingTarget
        /// </summary>
        /// <param name="motifs"></param>
        /// <param name="option2"></param>
        public SearchingTarget(List<char> motifs)
        {
            TargetMotifs = motifs;
        }

        /// <summary>
        /// Filter the peptide sequence based on the modification motifs
        /// </summary>
        /// <param name="peptideSequence"></param>
        /// <returns></returns>
        public bool MotifFilter(string peptideSequence)
        {
            if (peptideSequence == null)
                return false;

            List<char> motifList = new List<char>(); // the motif list in the peptide sequence
            for (int i = 0; i < peptideSequence.Length; i++)
            {
                if (peptideSequence[i] == '[' && i != 0)
                {
                    motifList.Add(peptideSequence[i - 1]);
                }
            }
            // if any motif in the peptide sequence is in the targetList, return true
            foreach (var motif in motifList)
            {
                if (TargetMotifs.Any(p => p == motif))
                {
                    return true;
                }
            }
            return false;
        }

    }
}
