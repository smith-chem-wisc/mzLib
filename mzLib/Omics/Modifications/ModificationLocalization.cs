namespace Omics.Modifications
{
    public static class ModificationLocalization
    {
        // This method is called a ton (8.8 billion times in Bottom-Up Jenkins as of 1.0.6) in MetaMorpheus. If changes are made, ensure they are efficient. 
        public static bool ModFits(Modification attemptToLocalize, string sequence, int digestionProductOneBasedIndex, int digestionProductLength, int bioPolymerOneBasedIndex)
        {
            // First find the capital letter...
            var motif = attemptToLocalize.Target.ToString();
            var motifStartLocation = -1;
            for (int i = 0; i < motif.Length; i++)
            {
                if (!char.IsUpper(motif[i])) 
                    continue;

                motifStartLocation = i;
                break;
            }

            // Look up starting at and including the capital letter
            var proteinToMotifOffset = bioPolymerOneBasedIndex - motifStartLocation - 1;
            var indexUp = 0;
            while (indexUp < motif.Length)
            {
                if (indexUp + proteinToMotifOffset < 0 || indexUp + proteinToMotifOffset >= sequence.Length
                    || !MotifMatches(motif[indexUp], sequence[indexUp + proteinToMotifOffset]))
                {
                    return false;
                }
                indexUp++;
            }
            switch (attemptToLocalize.LocationRestriction)
            {
                case "N-terminal." when bioPolymerOneBasedIndex > 2:
                case "Peptide N-terminal." when digestionProductOneBasedIndex > 1 || bioPolymerOneBasedIndex == 1:
                case "C-terminal." when bioPolymerOneBasedIndex < sequence.Length:
                case "Peptide C-terminal." when digestionProductOneBasedIndex < digestionProductLength || bioPolymerOneBasedIndex == sequence.Length:
                case "5'-terminal." when bioPolymerOneBasedIndex > 2:
                // first residue in oligo but not first in nucleic acid
                case "Oligo 5'-terminal." when digestionProductOneBasedIndex > 1
                                               || bioPolymerOneBasedIndex == 1:
                case "3'-terminal." when bioPolymerOneBasedIndex < sequence.Length:
                // not the last residue in oligo but not in nucleic acid
                case "Oligo 3'-terminal." when digestionProductOneBasedIndex < digestionProductLength
                                               || bioPolymerOneBasedIndex == sequence.Length:
                    return false;

                default:
                    // I guess Anywhere. and Unassigned. are true since how do you localize anywhere or unassigned.

                    return true;
            }
        }

        public static bool UniprotModExists(IBioPolymer bioPolymer, int i, Modification attemptToLocalize)
        {
            // uniprot mods with same mass takes precedence over variable mods
            if (bioPolymer.OneBasedPossibleLocalizedModifications.TryGetValue(i, out List<Modification> modsAtThisLocation)) {
                return modsAtThisLocation.Any(p => Math.Abs((double)(p.MonoisotopicMass - attemptToLocalize.MonoisotopicMass)) < 0.001 && p.ModificationType == "UniProt");
            }

            return false;
        }

        private static bool MotifMatches(char motifChar, char sequenceChar)
        {
            char upperMotifChar = char.ToUpper(motifChar);
            return upperMotifChar switch
            {
                'X' => true,
                'B' => sequenceChar is 'D' or 'N',
                'J' => sequenceChar is 'I' or 'L',
                'Z' => sequenceChar is 'E' or 'Q',
                _ => upperMotifChar == sequenceChar
            };
        }
    }
}