using Omics.Modifications;

namespace Omics.Digestion
{
    public class LysisProduct
    {
        protected string _baseSequence;

        public LysisProduct(IBioPolymer parent, int oneBasedStartResidue, int oneBasedEndResidue, int missedCleavages, 
            CleavageSpecificity cleavageSpecificityForFdrCategory, string? description = null, string? baseSequence = null)
        {
            Parent = parent;
            OneBasedStartResidue = oneBasedStartResidue;
            OneBasedEndResidue = oneBasedEndResidue;
            MissedCleavages = missedCleavages;
            CleavageSpecificityForFdrCategory = cleavageSpecificityForFdrCategory;
            Description = description;
            _baseSequence = baseSequence;
        }

        [field: NonSerialized] public IBioPolymer Parent { get; protected set; } // BioPolymer that this lysis product is a digestion product of
        public string Description { get; protected set; } //unstructured explanation of source
        public int OneBasedStartResidue { get; }// the residue number at which the peptide begins (the first residue in a protein is 1)
        public int OneBasedEndResidue { get; }// the residue number at which the peptide ends
        public int MissedCleavages { get; } // the number of missed cleavages this peptide has with respect to the digesting protease
        public virtual char PreviousResidue => OneBasedStartResidue > 1 ? Parent[OneBasedStartResidue - 2] : '-';

        public virtual char NextResidue => OneBasedEndResidue < Parent.Length ? Parent[OneBasedEndResidue] : '-';
        public string BaseSequence =>
            _baseSequence ??= Parent.BaseSequence.Substring(OneBasedStartResidue - 1,
                OneBasedEndResidue - OneBasedStartResidue + 1);
        public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; } //structured explanation of source
        public int Length => BaseSequence.Length; //how many residues long the peptide is
        public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

        protected static IEnumerable<Dictionary<int, Modification>> GetVariableModificationPatterns(Dictionary<int, List<Modification>> possibleVariableModifications, int maxModsForPeptide, int peptideLength)
        {
            if (possibleVariableModifications.Count == 0)
            {
                yield return null;
            }
            else
            {
                var possible_variable_modifications = new Dictionary<int, List<Modification>>(possibleVariableModifications);

                int[] base_variable_modification_pattern = new int[peptideLength + 4];
                var totalAvailableMods = possible_variable_modifications.Sum(b => b.Value == null ? 0 : b.Value.Count);
                for (int variable_modifications = 0; variable_modifications <= Math.Min(totalAvailableMods, maxModsForPeptide); variable_modifications++)
                {
                    foreach (int[] variable_modification_pattern in GetVariableModificationPatterns(new List<KeyValuePair<int, List<Modification>>>(possible_variable_modifications),
                        possible_variable_modifications.Count - variable_modifications, base_variable_modification_pattern, 0))
                    {
                        yield return GetNewVariableModificationPattern(variable_modification_pattern, possible_variable_modifications);
                    }
                }
            }
        }

        protected Dictionary<int, Modification> GetFixedModsOneIsNterminusOrFivePrime(int peptideLength,
            IEnumerable<Modification> allKnownFixedModifications)
        {
            var fixedModsOneIsNterminus = new Dictionary<int, Modification>(peptideLength + 3);
            foreach (Modification mod in allKnownFixedModifications)
            {
                switch (mod.LocationRestriction)
                {
                    case "5'-terminal.":
                    case "Oligo 5'-terminal.":
                    case "N-terminal.":
                    case "Peptide N-terminal.":
                        //the modification is protease associated and is applied to the n-terminal cleaved residue, not at the beginign of the protein
                        if (mod.ModificationType == "Protease" && ModificationLocalization.ModFits(mod, Parent.BaseSequence, 1, peptideLength, OneBasedStartResidue))
                        {
                            if (OneBasedStartResidue != 1)
                            {
                                fixedModsOneIsNterminus[2] = mod;
                            }
                        }
                        //Normal N-terminal peptide modification
                        else if (ModificationLocalization.ModFits(mod, Parent.BaseSequence, 1, peptideLength, OneBasedStartResidue))
                        {
                            fixedModsOneIsNterminus[1] = mod;
                        }
                        break;

                    case "Anywhere.":
                        for (int i = 2; i <= peptideLength + 1; i++)
                        {
                            if (ModificationLocalization.ModFits(mod, Parent.BaseSequence, i - 1, peptideLength, OneBasedStartResidue + i - 2))
                            {
                                fixedModsOneIsNterminus[i] = mod;
                            }
                        }
                        break;

                    case "3'-terminal.":
                    case "Oligo 3'-terminal.":
                    case "C-terminal.":
                    case "Peptide C-terminal.":
                        //the modification is protease associated and is applied to the c-terminal cleaved residue, not if it is at the end of the protein
                        if (mod.ModificationType == "Protease" && ModificationLocalization.ModFits(mod, Parent.BaseSequence, peptideLength, peptideLength, OneBasedStartResidue + peptideLength - 1))
                        {
                            if (OneBasedEndResidue != Parent.Length)
                            {
                                fixedModsOneIsNterminus[peptideLength + 1] = mod;
                            }

                        }
                        //Normal C-terminal peptide modification 
                        else if (ModificationLocalization.ModFits(mod, Parent.BaseSequence, peptideLength, peptideLength, OneBasedStartResidue + peptideLength - 1))
                        {
                            fixedModsOneIsNterminus[peptideLength + 2] = mod;
                        }
                        break;

                    default:
                        throw new NotSupportedException("This terminus localization is not supported.");
                }
            }
            return fixedModsOneIsNterminus;
        }


        private static IEnumerable<int[]> GetVariableModificationPatterns(List<KeyValuePair<int, List<Modification>>> possibleVariableModifications,
            int unmodifiedResiduesDesired, int[] variableModificationPattern, int index)
        {
            if (index < possibleVariableModifications.Count - 1)
            {
                if (unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    foreach (int[] new_variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications,
                        unmodifiedResiduesDesired - 1, variableModificationPattern, index + 1))
                    {
                        yield return new_variable_modification_pattern;
                    }
                }
                if (unmodifiedResiduesDesired < possibleVariableModifications.Count - index)
                {
                    for (int i = 1; i <= possibleVariableModifications[index].Value.Count; i++)
                    {
                        variableModificationPattern[possibleVariableModifications[index].Key] = i;
                        foreach (int[] new_variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications,
                            unmodifiedResiduesDesired, variableModificationPattern, index + 1))
                        {
                            yield return new_variable_modification_pattern;
                        }
                    }
                }
            }
            else
            {
                if (unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    yield return variableModificationPattern;
                }
                else
                {
                    for (int i = 1; i <= possibleVariableModifications[index].Value.Count; i++)
                    {
                        variableModificationPattern[possibleVariableModifications[index].Key] = i;
                        yield return variableModificationPattern;
                    }
                }
            }
        }

        private static Dictionary<int, Modification> GetNewVariableModificationPattern(int[] variableModificationArray,
            IEnumerable<KeyValuePair<int, List<Modification>>> possibleVariableModifications)
        {
            var modification_pattern = new Dictionary<int, Modification>();

            foreach (KeyValuePair<int, List<Modification>> kvp in possibleVariableModifications)
            {
                if (variableModificationArray[kvp.Key] > 0)
                {
                    modification_pattern.Add(kvp.Key, kvp.Value[variableModificationArray[kvp.Key] - 1]);
                }
            }

            return modification_pattern;
        }


    }
}
