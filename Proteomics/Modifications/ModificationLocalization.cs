using System.Linq;

namespace Proteomics
{
    public static class ModificationLocalization
    {
        public static bool ModFits(ModificationGeneral attemptToLocalize, Protein protein, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex)
        {
            // First find the capital letter...
            var motif = attemptToLocalize.Target;
            var hehe = motif.ToString().IndexOf(motif.ToString().First(b => char.IsUpper(b)));

            // Look up starting at and including the capital letter
            var proteinToMotifOffset = proteinOneBasedIndex - hehe - 1;
            var indexUp = 0;
            while (indexUp < motif.ToString().Length)
            {
                if (indexUp + proteinToMotifOffset < 0 || indexUp + proteinToMotifOffset >= protein.Length
                    || (!char.ToUpper(motif.ToString()[indexUp]).Equals('X')
                        && !char.ToUpper(motif.ToString()[indexUp]).Equals(protein.BaseSequence[indexUp + proteinToMotifOffset])))
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
            if (attemptToLocalize.LocationRestriction == "C-terminal." && peptideOneBasedIndex < peptideLength)
            {
                return false;
            }
            if (attemptToLocalize.LocationRestriction == "Peptide C-terminal." && proteinOneBasedIndex < protein.Length)
            {
                return false;
            }

            // I guess Anywhere. and Unassigned. are true since how do you localize anywhere or unassigned.

            return true;
        }
    }
}