using Omics.Modifications;

namespace Omics.Digestion
{
    public abstract class DigestionAgent
    {
        protected DigestionAgent(string name, CleavageSpecificity cleavageSpecificity, List<DigestionMotif> motifList, Modification cleavageMod)
        {
            Name = name;
            CleavageSpecificity = cleavageSpecificity;
            DigestionMotifs = motifList ?? new List<DigestionMotif>();
            CleavageMod = cleavageMod;
        }

        public string Name { get; init; }
        public CleavageSpecificity CleavageSpecificity { get; init; }
        public List<DigestionMotif> DigestionMotifs { get; init; }
        public Modification CleavageMod { get; set; }

        public override string ToString()
        {
            return Name;
        }

        /// <summary>
        /// Is length of given peptide okay, given minimum and maximum?
        /// </summary>
        /// <param name="peptideLength"></param>
        /// <param name="minPeptideLength"></param>
        /// <param name="maxPeptideLength"></param>
        /// <returns></returns>
        protected static bool OkayLength(int peptideLength, int minPeptideLength, int maxPeptideLength)
        {
            return OkayMinLength(peptideLength, minPeptideLength) && OkayMaxLength(peptideLength, maxPeptideLength);
        }

        /// <summary>
        /// Is length of given peptide okay, given minimum?
        /// </summary>
        /// <param name="peptideLength"></param>
        /// <param name="minPeptideLength"></param>
        /// <returns></returns>
        protected static bool OkayMinLength(int peptideLength, int minPeptideLength)
        {
            return peptideLength >= minPeptideLength;
        }

        /// <summary>
        /// Is length of given peptide okay, given maximum?
        /// </summary>
        /// <param name="peptideLength"></param>
        /// <param name="maxPeptideLength"></param>
        /// <returns></returns>
        protected static bool OkayMaxLength(int? peptideLength, int maxPeptideLength)
        {
            return !peptideLength.HasValue || peptideLength <= maxPeptideLength;
        }

        /// <summary>
        /// Gets the indices after which this protease will cleave a given protein sequence
        /// </summary>
        /// <param name="proteinSequence"></param>
        /// <returns></returns>
        public List<int> GetDigestionSiteIndices(string proteinSequence)
        {
            var indices = new List<int>();

            for (int r = 0; r < proteinSequence.Length; r++)
            {
                var cutSiteIndex = -1;
                bool cleavagePrevented = false;

                foreach (DigestionMotif motif in DigestionMotifs)
                {
                    var motifResults = motif.Fits(proteinSequence, r);
                    bool motifFits = motifResults.Item1;
                    bool motifPreventsCleavage = motifResults.Item2;

                    if (motifFits && r + motif.CutIndex < proteinSequence.Length)
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

            indices.Add(0); // The start of the protein is treated as a cleavage site to retain the n-terminal peptide
            indices.Add(proteinSequence.Length); // The end of the protein is treated as a cleavage site to retain the c-terminal peptide
            return indices.Distinct().OrderBy(i => i).ToList();
        }

       
    }
}
