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
        public readonly string ExcludeFromWildcard;

        public DigestionMotif(string inducingCleavage, string preventingCleavage, int cutIndex, string excludeFromWildcard)
        {
            this.InducingCleavage = inducingCleavage;
            this.PreventingCleavage = preventingCleavage;
            this.CutIndex = cutIndex;
            this.ExcludeFromWildcard = excludeFromWildcard;
        }

        // parsing cleavage rules syntax
        public static List<DigestionMotif> ParseDigestionMotifsFromString(string motifsString)
        {
            motifsString = motifsString.Replace("\"", string.Empty).Replace(" ", string.Empty);

            // throws exception if non-supported characters are used
            if (Regex.Match(motifsString, @"[^a-zA-Z0-9|,[\]{}]+").Success)
            {
                throw new MzLibException("Unrecognized protease syntax. The digestion motif can only contain letters and {}[]|");
            }
            // throws exception if user attempts separate multiple preventing cleavages using commas
            if (Regex.Match(motifsString, @"\[([\w]*,+[\w]*)*\]").Success)
            {
                throw new MzLibException("Unrecognized protease syntax. Please create a separate motif for each sequence preventing cleavage (comma separated).");
            }
            // throws exception if user attempts separate multiple wildcard exclusions
            if (Regex.Match(motifsString, @"\{([\w]*,+[\w]*)*\}").Success)
            {
                throw new MzLibException("Unrecognized protease syntax. Please create a separate motif for each wildcard exclusion (comma separated).");
            }

            string[] motifStrings = motifsString.Split(',');
            var motifs = new List<DigestionMotif>();

            for (int i = 0; i < motifStrings.Length; i++)
            {
                string motifString = motifStrings[i];
                motifs.Add(ParseDigestionMotifFromString(motifString));
            }
            return motifs;
        }

        private static DigestionMotif ParseDigestionMotifFromString(string motifString)
        {
            string inducingCleavage;
            string preventingCleavage = null;
            string excludingWC = null;
            int cutIndex = 0;

            if (motifString.Contains("{") && !motifString.Contains("}")
                || !motifString.Contains("{") && motifString.Contains("}")
                || motifString.Contains("[") && !motifString.Contains("]")
                || !motifString.Contains("[") && motifString.Contains("]"))
            {
                throw new MzLibException("Unrecognized protease syntax. Please close any brackets used.");
            }

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
                if (Regex.Matches(motifString.ToUpper(), "X").Count != excludingWC.Length)
                {
                    throw new MzLibException("Unrecognized protease syntax. Please have equal number of wildcards for multi-letter wildcard exclusions.");
                }
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

        public (bool, bool) Fits(string sequence, int location)
        {
            bool fits = true;
            char currentResidue;
            int m;

            // check for inducing cleavage
            for (m = 0; m < InducingCleavage.Length && fits; m++) // handle patterns
            {
                if (location + m >= sequence.Length)
                {
                    fits = false;
                }
                else
                {
                    currentResidue = sequence[location + m];
                    if (!MotifMatches(InducingCleavage[m], currentResidue))
                    {
                        fits = false;
                    }
                }
            }

            bool prevents = false;
            // check for preventing cleavage
            if (fits && PreventingCleavage != null)
            {
                prevents = true;
                for (int n = 0; n < PreventingCleavage.Length && prevents; n++)
                {
                    if (location + m + n >= sequence.Length || location - PreventingCleavage.Length + 1 + n < 0)
                    {
                        prevents = false;
                    }
                    else
                    {
                        currentResidue = CutIndex != 0 ? sequence[location + m + n] : sequence[location - PreventingCleavage.Length + 1 + n];
                        if (!MotifMatches(PreventingCleavage[n], currentResidue))
                        {
                            prevents = false;
                        }
                    }
                }

                fits = prevents ? false : true;
            }

            return (fits, prevents);
        }

        private bool MotifMatches(char motifChar, char sequenceChar)
        {
            return motifChar.Equals('X') && !sequenceChar.ToString().Equals(ExcludeFromWildcard)
                || motifChar.Equals(sequenceChar)
                || motifChar.Equals('B') && B.Contains(sequenceChar)
                || motifChar.Equals('J') && J.Contains(sequenceChar)
                || motifChar.Equals('Z') && Z.Contains(sequenceChar);
        }
    }
}