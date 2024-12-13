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
                if (char.IsUpper(motif[i]))
                {
                    motifStartLocation = i;
                    break;
                }
            }

            // Look up starting at and including the capital letter
            var proteinToMotifOffset = bioPolymerOneBasedIndex - motifStartLocation - 1;
            for (int indexUp = 0; indexUp < motif.Length; indexUp++)
            {
                int sequenceIndex = indexUp + proteinToMotifOffset;
                if (sequenceIndex < 0 || sequenceIndex >= sequence.Length || !MotifMatches(motif[indexUp], sequence[sequenceIndex]))
                {
                    return false;
                }
            }

            return attemptToLocalize.LocationRestriction switch
            {
                "N-terminal." when bioPolymerOneBasedIndex > 2 => false,
                "Peptide N-terminal." when digestionProductOneBasedIndex > 1 => false,
                "C-terminal." when bioPolymerOneBasedIndex < sequence.Length => false,
                "Peptide C-terminal." when digestionProductOneBasedIndex < digestionProductLength => false,
                "5'-terminal." when bioPolymerOneBasedIndex > 2 => false,
                // first residue in oligo but not first in nucleic acid
                "Oligo 5'-terminal." when digestionProductOneBasedIndex > 1 || bioPolymerOneBasedIndex == 1 => false,
                "3'-terminal." when bioPolymerOneBasedIndex < sequence.Length => false,
                // not the last residue in oligo but not in nucleic acid
                "Oligo 3'-terminal." when digestionProductOneBasedIndex < digestionProductLength ||
                                          bioPolymerOneBasedIndex == sequence.Length => false,
                _ => true
            };
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