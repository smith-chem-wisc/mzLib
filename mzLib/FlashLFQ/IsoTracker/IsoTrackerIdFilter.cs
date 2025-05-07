using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using Omics.Modifications;

namespace FlashLFQ.IsoTracker
{
    public class IsoTrackerIdFilter
    {
        public List<ModificationMotif> TargetMotifs; // motif string list
        public List<Regex> TargetMotifPattern; // motif regex list for sequence checking

        /// <summary>
        /// Directly use the motif string to create a SearchingTarget
        /// </summary>
        /// <param name="motifs"></param>
        /// <param name="option2"></param>
        public IsoTrackerIdFilter(List<char> motifs)
        {
            TargetMotifs = new List<ModificationMotif>();
            if (motifs == null) return; // If the user does not set the motif list, then no filter will be applied and no regex will be created.
            foreach (var motif in motifs)
            {
                ModificationMotif.TryGetMotif(motif.ToString(), out var motifObj);
                if (motifObj != null)
                {
                    TargetMotifs.Add(motifObj);
                }
            }
            if (TargetMotifs.Count == 0) return; // Try to handle any cases where the user input is not valid.
            PatternBuilding();
        }

        /// <summary>
        /// The parameterless constructor for toml file reader.
        /// But we still need to deal with the toml stuff for my new parameter in the future.
        /// </summary>
        public IsoTrackerIdFilter(): this(null)
        {
        }

        public void PatternBuilding()
        {
            TargetMotifPattern = new List<Regex>();
            foreach (var motif in TargetMotifs)
            {
                string pattern = Regex.Escape(motif.ToString()) + @"\[[^\]]+\]";
                TargetMotifPattern.Add(new Regex(pattern));
            }
        }

        /// <summary>
        /// Filter the peptide sequence based on the modification motifs
        /// </summary>
        /// <param name="peptideSequence"></param>
        /// <returns></returns>
        public bool ContainsAcceptableModifiedResidue(string peptideSequence)
        {
            if (peptideSequence == null)
                return false;
            if (TargetMotifs.Count == 0)
            {
                return true; // no motif filter
            }
            return TargetMotifPattern.Any(p=> p.IsMatch(peptideSequence));
        }

    }
}
