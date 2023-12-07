using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
            var indices = new List<int>();

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

            indices.Add(0); // The start of the protein is treated as a cleavage site to retain the n-terminal peptide
            indices.Add(sequence.Length); // The end of the protein is treated as a cleavage site to retain the c-terminal peptide
            return indices.Distinct().OrderBy(i => i).ToList();
        }
    }
}
