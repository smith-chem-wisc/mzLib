using MzLibUtil;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

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
        public readonly string ExcludingWC;

        public DigestionMotif(string inducingCleavage, string preventingCleavage, int cutIndex, string excludingWC)
        {
            this.InducingCleavage = inducingCleavage;
            this.PreventingCleavage = preventingCleavage;
            this.CutIndex = cutIndex;
            this.ExcludingWC = excludingWC;
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
            string excludingWC = null;
            int cutIndex = 0;

            // find preventing cleavage
            if (motifString.Contains("["))
            {
                int start = motifString.IndexOf("[") + 1;
                int end = motifString.IndexOf("]");

                preventingCleavage = motifString.Substring(start, end - start);
                motifString = Regex.Replace(motifString, @"\[[a-zA-Z]+\]", string.Empty);
            }

            // finds wildcard exceptions
            if (motifString.Contains("{"))
            {
                int start = motifString.IndexOf("{") + 1;
                int end = motifString.IndexOf("}");

                excludingWC = motifString.Substring(start, end - start);
                motifString = Regex.Replace(motifString, @"\{[a-zA-Z]+\}", string.Empty);
            }

            // finds motif cut index
            for (int j = 0; j < motifString.Length; j++)
            {
                if (motifString[j] == '|')
                {
                    cutIndex = j;
                    break;
                }
            }

            motifString = motifString.Replace("|", string.Empty);
            inducingCleavage = motifString;

            return new DigestionMotif(inducingCleavage, preventingCleavage, cutIndex, excludingWC);
        }

        public bool Fits(string sequence, int location)
        {
            bool fits = true;

            // check for inducing cleavage
            int m;
            for (m = 0; m < InducingCleavage.Length && fits; m++) // handle patterns
            {
                char currentResidue = sequence[location + m];

                if (!MotifMatches(InducingCleavage[m], currentResidue))
                {
                    fits = false;
                }
            }

            // check for preventing cleavage
            if (fits && PreventingCleavage != null)
            {
                bool match = true;
                char currentResidue;
                for (int n = 0; n < PreventingCleavage.Length && match; n++)
                {
                    if (location + m + n > sequence.Length || location - PreventingCleavage.Length + 1 + n < 0)
                    {
                        match = false;
                    } else
                    {
                        currentResidue = CutIndex != 0 ? sequence[location + m + n] : sequence[location - PreventingCleavage.Length + 1 + n];                      
                        if (!PreventingCleavage[n].Equals(currentResidue))
                        {
                            match = false;
                        }
                    }
                }

                fits = match ? false : true;
            }

            return fits;
        }
        
        private bool MotifMatches(char motifChar, char sequenceChar)
        {
            return motifChar.Equals('X') && !sequenceChar.ToString().Equals(ExcludingWC)
                || motifChar.Equals(sequenceChar)
                || motifChar.Equals('B') && B.Contains(sequenceChar)
                || motifChar.Equals('J') && J.Contains(sequenceChar)
                || motifChar.Equals('Z') && Z.Contains(sequenceChar);
        }
    }
}
