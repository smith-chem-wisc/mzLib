using MzLibUtil;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics.ProteolyticDigestion
{
    public class DigestionMotif
    {
        private static char[] B = new char[] { 'D', 'N' };
        private static char[] J = new char[] { 'I', 'L' };
        private static char[] Z = new char[] { 'E', 'Q' };

        public readonly string InducingCleavage;
        public readonly string PreventingCleavage;
        public readonly int CutIndex;

        public DigestionMotif(string inducingCleavage, string preventingCleavage, int cutIndex)
        {
            this.InducingCleavage = inducingCleavage;
            this.PreventingCleavage = preventingCleavage;
            this.CutIndex = cutIndex;
        }

        // parsing cleavage rules syntax
        // TODO: add more error checking exceptions/informative messages
        public static List<DigestionMotif> ParseDigestionMotifsFromString(string motifsString)
        {
            string[] motifStrings = motifsString.Replace("\"", string.Empty).Replace(" ", string.Empty).Split(',');
            var motifs = new List<DigestionMotif>();

            for (int i = 0; i < motifStrings.Length; i++)
            {
                string motifString = motifStrings[i];

                // error checking for syntax: missing "|" for cleavage terminus
                if (!string.IsNullOrEmpty(motifString) && !motifString.Contains("|"))
                {
                    throw new MzLibException("Digestion motif " + motifString + " did not contain a | to specify the cleavage terminus.");
                }


                motifs.Add(ParseDigestionMotifFromString(motifString));
            }

            return motifs;
        }

        private static DigestionMotif ParseDigestionMotifFromString(string motifString)
        {
            string inducingCleavage = null;
            string preventingCleavage = null;
            int cutIndex = 0;

            // finds motif cut index
            for (int j = 0; j < motifString.Length; j++)
            {
                if (motifString[j] == '[' || motifString[j] == '|')
                {
                    cutIndex = j;
                    break;
                }
            }

            motifString = motifString.Replace("|", string.Empty);

            // identify inducing and preventing cleavage sequences 
            if (!(motifString.Contains("[")))
            {
                inducingCleavage = motifString;
                preventingCleavage = null;
            }
            else
            {
                int start = motifString.IndexOf("[") + 1;
                int end = motifString.IndexOf("]");
                inducingCleavage = motifString.Substring(0, start - 1);
                preventingCleavage = motifString.Substring(start, end - start);
            }

            return new DigestionMotif(inducingCleavage, preventingCleavage, cutIndex);
        }

        /// <summary>
        /// Returns non-distinct digestion sites.
        /// </summary>
        public IEnumerable<int> GetDigestionIndicies(string sequence)
        {
            yield return 0;
            yield return sequence.Length - 1;

            for (int r = 0; r < sequence.Length; r++)
            {
                if (MotifFits(sequence, r))
                {
                    yield return r;
                }
            }
        }

        private bool MotifFits(string sequence, int location)
        {
            bool fits = true;

            // check for inducing cleavage
            for (int r = location; r < sequence.Length && fits; r++)
            {
                for (int m = 0; m < InducingCleavage.Length && fits; m++)
                {
                    char currentResidue = sequence[r + m];

                    if (!MotifMatches(InducingCleavage[m], currentResidue))
                    {
                        fits = false;
                    }
                }
            }

            // check for preventing cleavage
            if (fits)
            {
                
            }

            return fits;
        }

        private static bool MotifMatches(char motifChar, char sequenceChar)
        {
            return motifChar.Equals('X')
                || motifChar.Equals(sequenceChar)
                || motifChar.Equals('B') && B.Contains(sequenceChar)
                || motifChar.Equals('J') && J.Contains(sequenceChar)
                || motifChar.Equals('Z') && Z.Contains(sequenceChar);
        }
    }
}
