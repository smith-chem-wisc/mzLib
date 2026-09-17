using MzLibUtil;
using Omics.Modifications;

namespace Omics.Digestion
{
    public abstract class DigestionAgent
    {
        protected static readonly HashSetPool<int> HashSetPool = new HashSetPool<int>(8);

        protected DigestionAgent(string name, CleavageSpecificity cleavageSpecificity, List<DigestionMotif> motifList, Modification cleavageMod)
        {
            Name = name;
            CleavageSpecificity = cleavageSpecificity;
            DigestionMotifs = motifList ?? new List<DigestionMotif>();
            CleavageMod = cleavageMod;
        }

        public readonly string Name;
        public CleavageSpecificity CleavageSpecificity { get; init; }
        public List<DigestionMotif> DigestionMotifs { get; init; }
        public Modification CleavageMod { get; set; }

        /// <summary>
        /// This method is used to determine cleavage specificity if the cleavage specificity is unknown
        /// This occurs in the speedy nonspecific/semispecific searches when digesting post-search
        /// </summary>
        /// <returns></returns>
        public CleavageSpecificity GetCleavageSpecificity(IBioPolymer bioPolymer, int startIndex, int endIndex, bool retainMethionine)
        {
            int cleavableMatches = 0;
            if (CleavageSpecificity != CleavageSpecificity.SingleN && CleavageSpecificity != CleavageSpecificity.SingleC) //if it's single protease, don't bother
            {
                List<int> indicesToCleave = GetDigestionSiteIndices(bioPolymer.BaseSequence);
                //if the start index is a cleavable index (-1 because one based) OR if the start index is after a cleavable methionine
                if (indicesToCleave.Contains(startIndex - 1) ||
                    (startIndex == 2 && bioPolymer.BaseSequence[0] == 'M' && !retainMethionine) ||
                    bioPolymer.TruncationProducts.Any(x => x.OneBasedBeginPosition == startIndex))
                {
                    cleavableMatches++;
                }
                //if the end index is a cleavable index
                if (indicesToCleave.Contains(endIndex) ||
                    bioPolymer.TruncationProducts.Any(x => x.OneBasedEndPosition == endIndex))
                {
                    cleavableMatches++;
                }
            }
            if (cleavableMatches == 0) //if neither were cleavable, (or it's singleN/C) then it's nonspecific
            {
                return CleavageSpecificity.None;
            }
            else if (cleavableMatches == 1) //if one index was cleavable, then it's semi specific
            {
                return CleavageSpecificity.Semi;
            }
            else //2 if both, then it's fully speific
            {
                return CleavageSpecificity.Full;
            }
        }

        public override string ToString()
        {
            return Name;
        }

        public override bool Equals(object? obj)
        {
            return obj is DigestionAgent agent && agent.Name == Name;
        }

        public override int GetHashCode()
        {
            return Name.GetHashCode();
        }

        /// <summary>
        /// Is length of given peptide okay, given minimum and maximum?
        /// </summary>
        /// <param name="length"></param>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <returns></returns>
        protected static bool ValidLength(int length, int minLength, int maxLength)
        {
            return ValidMinLength(length, minLength) && ValidMaxLength(length, maxLength);
        }

        /// <summary>
        /// Is length of given peptide okay, given minimum?
        /// </summary>
        /// <param name="length"></param>
        /// <param name="minLength"></param>
        /// <returns></returns>
        protected static bool ValidMinLength(int length, int minLength)
        {
            return length >= minLength;
        }

        /// <summary>
        /// Is length of given peptide okay, given maximum?
        /// </summary>
        /// <param name="length"></param>
        /// <param name="maxLength"></param>
        /// <returns></returns>
        protected static bool ValidMaxLength(int? length, int maxLength)
        {
            return !length.HasValue || length <= maxLength;
        }

        /// <summary>
        /// Gets the indices after which this protease will cleave a given protein sequence
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        public List<int> GetDigestionSiteIndices(string sequence)
        {
            var indices = HashSetPool.Get(); // use hash set to ensure no duplicates
            try // Try block is to ensure that, even if an error gets thrown, the hashset is returned to the pool
            {
                indices.Add(0); // The start of the protein is treated as a cleavage site to retain the n-terminal peptide

                for (int r = 0; r < sequence.Length; r++)
                {
                    var cutSiteIndex = -1;
                    bool cleavagePrevented = false;

                    foreach (DigestionMotif motif in DigestionMotifs)
                    {
                        var motifResults = motif.Fits(sequence, r);
                        bool motifFits = motifResults.Item1;
                        bool motifPreventsCleavage = motifResults.Item2;

                        if (motifFits && r + motif.CutIndex < sequence.Length)
                        {
                            cutSiteIndex = Math.Max(r + motif.CutIndex, cutSiteIndex);
                        }

                        if (motifPreventsCleavage) // if any motif prevents cleave
                        {
                            cleavagePrevented = true;
                        }
                    }

                    // if no motif prevents cleave
                    if (!cleavagePrevented && cutSiteIndex != -1)
                    {
                        indices.Add(cutSiteIndex);
                    }
                }

                indices.Add(sequence.Length); // The end of the protein is treated as a cleavage site to retain the c-terminal peptide
                return indices.ToList(); // convert the hashset to a list for return. 
            }
            finally
            {
                // return hashset to pool. This clears it and gets it ready for the next time it is needed from the pool.
                HashSetPool.Return(indices);
            }
        }
    }
}
