using MzLibUtil;
using Omics.Modifications;

namespace Omics.Digestion
{
    public abstract class DigestionProduct
    {
        protected static readonly DictionaryPool<int, List<Modification>> DictionaryPool = new();
        protected static readonly DictionaryPool<int, Modification> FixedModDictionaryPool = new(8);

        protected string _baseSequence;

        protected DigestionProduct(IBioPolymer parent, int oneBasedStartResidue, int oneBasedEndResidue, int missedCleavages, 
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

        public virtual char PreviousResidue => Parent is null ? '-' : OneBasedStartResidue > 1 ? Parent[OneBasedStartResidue - 2] : '-';

        public virtual char NextResidue => Parent is null ? '-' : OneBasedEndResidue < Parent.Length ? Parent[OneBasedEndResidue] : '-';

        public string BaseSequence =>
            _baseSequence ??= Parent.BaseSequence.Substring(OneBasedStartResidue - 1,
                OneBasedEndResidue - OneBasedStartResidue + 1);
        public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; } //structured explanation of source
        public int Length => BaseSequence.Length; //how many residues long the peptide is
        public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

        /// <summary>
        /// Generates all possible variable modification patterns for a peptide.
        /// </summary>
        /// <param name="possibleVariableModifications">A dictionary of possible variable modifications with their positions.</param>
        /// <param name="maxModsForPeptide">The maximum number of modifications allowed for the peptide.</param>
        /// <param name="peptideLength">The length of the peptide.</param>
        /// <returns>An enumerable of dictionaries representing different modification patterns.</returns>
        /// <remarks>
        /// This method generates all possible combinations of variable modifications for a given peptide. 
        /// It first calculates the total number of available modifications and the maximum number of variable modifications allowed.
        /// Then, it iterates through all possible numbers of modifications and generates the corresponding modification patterns.
        /// The returned dictionary is then appended with fixed modifications and used to construct a peptide with set mods
        /// </remarks>
        protected static IEnumerable<Dictionary<int, Modification>> GetVariableModificationPatterns(Dictionary<int, List<Modification>> possibleVariableModifications, int maxModsForPeptide, int peptideLength)
        {
            if (possibleVariableModifications.Count <= 0) 
                yield break;

            int[] baseVariableModificationPattern = new int[peptideLength + 4];
            int totalAvailableMods = possibleVariableModifications.Values.Sum(modList => modList?.Count ?? 0);
            int maxVariableMods = Math.Min(totalAvailableMods, maxModsForPeptide);

            for (int variable_modifications = 0; variable_modifications <= maxVariableMods; variable_modifications++)
            {
                foreach (int[] variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications.ToList(),
                             possibleVariableModifications.Count - variable_modifications, baseVariableModificationPattern, 0))
                {
                    // use modification pattern to construct a dictionary of modifications for the peptide
                    var modificationPattern = new Dictionary<int, Modification>(possibleVariableModifications.Count);

                    foreach (KeyValuePair<int, List<Modification>> kvp in possibleVariableModifications)
                    {
                        int modIndex = variable_modification_pattern[kvp.Key] - 1;
                        if (modIndex >= 0)
                        {
                            modificationPattern.Add(kvp.Key, kvp.Value[modIndex]);
                        }
                    }

                    yield return modificationPattern;
                }
            }
        }

        /// <summary>
        /// Sets the fixed modifications for the peptide, considering the N-terminal and C-terminal positions, by populating the <paramref name="fixedModsOneIsNterminus"/> dictionary.
        /// </summary>
        /// <param name="length">The length of the peptide.</param>
        /// <param name="allKnownFixedModifications">A collection of all known fixed modifications.</param>
        /// <param name="fixedModsOneIsNterminus">A reference to a dictionary that will hold the fixed modifications, with the key representing the position.</param>
        /// <remarks>
        /// This method iterates through all known fixed modifications and assigns them to the appropriate positions in the peptide.
        /// It considers different location restrictions such as N-terminal, C-terminal, and anywhere within the peptide.
        /// </remarks>
        protected void PopulateFixedModsOneIsNorFivePrimeTerminus(int length,
            IEnumerable<Modification> allKnownFixedModifications, ref Dictionary<int, Modification> fixedModsOneIsNterminus)
        {
            foreach (Modification mod in allKnownFixedModifications)
            {
                switch (mod.LocationRestriction)
                {
                    case "5'-terminal.":
                    case "Oligo 5'-terminal.":
                    case "N-terminal.":
                    case "Peptide N-terminal.":
                        //the modification is protease associated and is applied to the n-terminal cleaved residue, not at the beginning of the protein
                        if (ModificationLocalization.ModFits(mod, Parent.BaseSequence, 1, length, OneBasedStartResidue))
                        {
                            if (mod.ModificationType == "Protease")
                            {
                                if (OneBasedStartResidue != 1)
                                    fixedModsOneIsNterminus[2] = mod;
                            }
                            else //Normal N-terminal peptide modification
                                fixedModsOneIsNterminus[1] = mod;
                        }
                        break;

                    case "Anywhere.":
                        for (int i = 2; i <= length + 1; i++)
                        {
                            if (ModificationLocalization.ModFits(mod, Parent.BaseSequence, i - 1, length, OneBasedStartResidue + i - 2))
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
                        if (ModificationLocalization.ModFits(mod, Parent.BaseSequence, length, length, OneBasedStartResidue + length - 1))
                        {
                            if (mod.ModificationType == "Protease")
                            {
                                if (OneBasedEndResidue != Parent.Length)
                                    fixedModsOneIsNterminus[length + 1] = mod;
                            }
                            else //Normal C-terminal peptide modification 
                                fixedModsOneIsNterminus[length + 2] = mod;
                        }
                        break;

                    default:
                        throw new NotSupportedException("This terminus localization is not supported.");
                }
            }
        }

        /// <summary>
        /// Recursively generates all possible variable modification patterns for a peptide.
        /// </summary>
        /// <param name="possibleVariableModifications">A list of key-value pairs representing possible variable modifications and their positions.</param>
        /// <param name="unmodifiedResiduesDesired">The number of unmodified residues desired in the pattern.</param>
        /// <param name="variableModificationPattern">An array representing the current modification pattern.</param>
        /// <param name="index">The current index in the list of possible modifications.</param>
        /// <returns>An enumerable of arrays representing different modification patterns. The array index corresponds to the location of the modification
        /// in the peptide, while the value at that index determines which index in the <paramref name="possibleVariableModifications"/> list of modifications 
        /// to add to the final variable modification pattern </returns>
        /// <remarks>
        /// This method uses recursion to generate all possible combinations of variable modifications for a given peptide.
        /// It considers both modified and unmodified residues and generates patterns accordingly.
        /// </remarks>
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

        /// <summary>
        /// Determines if a modification can be applied to the N-terminal or 5' end of the peptide.
        /// </summary>
        /// <param name="mod">The modification to check.</param>
        /// <param name="peptideLength">The length of the peptide.</param>
        /// <returns>True if the modification can be applied to the N-terminal or 5' end; otherwise, false.</returns>
        protected bool CanBeNTerminalOrFivePrime(Modification mod, int peptideLength)
        {
            return mod.LocationRestriction is "5'-terminal." or "Oligo 5'-terminal." or "N-terminal." or "Peptide N-terminal."
                   && ModificationLocalization.ModFits(mod, Parent.BaseSequence, 1, peptideLength, OneBasedStartResidue);
        }

        /// <summary>
        /// Determines if a modification can be applied to the C-terminal or 3' end of the peptide.
        /// </summary>
        /// <param name="mod">The modification to check.</param>
        /// <param name="peptideLength">The length of the peptide.</param>
        /// <returns>True if the modification can be applied to the C-terminal or 3' end; otherwise, false.</returns>
        protected bool CanBeCTerminalOrThreePrime(Modification mod, int peptideLength)
        {
            return mod.LocationRestriction is "3'-terminal." or "Oligo 3'-terminal." or "C-terminal." or "Peptide C-terminal."
                   && ModificationLocalization.ModFits(mod, Parent.BaseSequence, peptideLength, peptideLength, OneBasedStartResidue + peptideLength - 1);
        }
    }
}
