using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.Metadata.Ecma335;
using System.Runtime.InteropServices.Marshalling;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using Omics.Modifications;

namespace FlashLFQ.IsoTracker
{
    public class SearchingTarget
    {
        public HashSet<ModificationMotif> ModificationMotifs { get; set; }
        public List<string> ModList;
        // Define a default list with some example values


        /// <summary>
        /// Directly use the motif string to create a SearchingTarget
        /// </summary>
        /// <param name="motifs"></param>
        /// <param name="option2"></param>
        public SearchingTarget(List<string> motifs)
        {
            ModificationMotifs = new HashSet<ModificationMotif>();
            ModList = motifs;
            if (motifs == null)
                return;

            foreach (var motif in motifs)
            {
                if (ModificationMotif.TryGetMotif(motif, out var modificationMotif))
                {
                    ModificationMotifs.Add(modificationMotif);
                }
            }
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

            return ModificationMotifs.Any(p=> p.ModficationMotifPattern.IsMatch(peptideSequence));
        }
    }
}
