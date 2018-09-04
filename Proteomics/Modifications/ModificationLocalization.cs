using System.Linq;

namespace Proteomics
{
    public static class ModificationLocalization
    {
        public static bool ModFits(Modification attemptToLocalize, string proteinSequence, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex)
        {
            // First find the capital letter...
            var motif = attemptToLocalize.Target;
            var motifStartLocation = motif.ToString().IndexOf(motif.ToString().First(b => char.IsUpper(b)));

            // Look up starting at and including the capital letter
            var proteinToMotifOffset = proteinOneBasedIndex - motifStartLocation - 1;
            var indexUp = 0;
            while (indexUp < motif.ToString().Length)
            {
                if (indexUp + proteinToMotifOffset < 0 || indexUp + proteinToMotifOffset >= proteinSequence.Length
                    || !MotifMatches(motif.ToString()[indexUp], proteinSequence[indexUp + proteinToMotifOffset]))
                {
                    return false;
                }
                indexUp++;
            }
            if (attemptToLocalize.LocationRestriction == "N-terminal." && proteinOneBasedIndex > 2)
            {
                return false;
            }
            if (attemptToLocalize.LocationRestriction == "Peptide N-terminal." && peptideOneBasedIndex > 1)
            {
                return false;
            }
            if (attemptToLocalize.LocationRestriction == "C-terminal." && proteinOneBasedIndex < proteinSequence.Length)
            {
                return false;
            }
            if (attemptToLocalize.LocationRestriction == "Peptide C-terminal." && peptideOneBasedIndex < peptideLength)
            {
                return false;
            }

            // I guess Anywhere. and Unassigned. are true since how do you localize anywhere or unassigned.

            return true;
        }

        private static bool MotifMatches(char motifChar, char sequenceChar)
        {
            char upperMotifChar = char.ToUpper(motifChar);
            return upperMotifChar.Equals('X')
                || upperMotifChar.Equals(sequenceChar)
                || upperMotifChar.Equals('B') && new[] { 'D', 'N' }.Contains(sequenceChar)
                || upperMotifChar.Equals('J') && new[] { 'I', 'L' }.Contains(sequenceChar)
                || upperMotifChar.Equals('Z') && new[] { 'E', 'Q' }.Contains(sequenceChar);
        }
    }
}