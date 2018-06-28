using System.Linq;

namespace Proteomics
{
    public static class ModificationLocalization
    {
        public static bool ModFits(ModificationWithMass attemptToLocalize, Protein protein, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex)
        {
            // First find the capital letter...
            var motif = attemptToLocalize.motif;
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
            if (attemptToLocalize.terminusLocalization == TerminusLocalization.NProt && proteinOneBasedIndex > 2)
            {
                return false;
            }
            if (attemptToLocalize.terminusLocalization == TerminusLocalization.NPep && peptideOneBasedIndex > 1)
            {
                return false;
            }
            if (attemptToLocalize.terminusLocalization == TerminusLocalization.PepC && peptideOneBasedIndex < peptideLength)
            {
                return false;
            }
            if (attemptToLocalize.terminusLocalization == TerminusLocalization.ProtC && proteinOneBasedIndex < protein.Length)
            {
                return false;
            }
            return true;
        }
    }
}