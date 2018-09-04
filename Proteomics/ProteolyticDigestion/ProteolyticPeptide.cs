using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics.ProteolyticDigestion
{
    /// <summary>
    /// Product of digesting a protein
    /// Contains methods for modified peptide combinitorics
    /// </summary>
    [Serializable]
    public class ProteolyticPeptide
    {
        protected string _baseSequence;

        internal ProteolyticPeptide(Protein protein, int oneBasedStartResidueInProtein, int oneBasedEndResidueInProtein, int missedCleavages, string peptideDescription = null)
        {
            _protein = protein;
            OneBasedStartResidueInProtein = oneBasedStartResidueInProtein;
            OneBasedEndResidueInProtein = oneBasedEndResidueInProtein;
            MissedCleavages = missedCleavages;
            PeptideDescription = peptideDescription;
        }

        [NonSerialized] private Protein _protein; // protein that this peptide is a digestion product of
        public int OneBasedStartResidueInProtein { get; } // the residue number at which the peptide begins (the first residue in a protein is 1)
        public int OneBasedEndResidueInProtein { get; } // the residue number at which the peptide ends
        public int MissedCleavages { get; } // the number of missed cleavages this peptide has with respect to the digesting protease
        public string PeptideDescription { get; }
        public int Length { get { return BaseSequence.Length; } } //how many residues long the peptide is

        public virtual char PreviousAminoAcid
        {
            get
            {
                return OneBasedStartResidueInProtein > 1 ? Protein[OneBasedStartResidueInProtein - 2] : '-';
            }
        }

        public virtual char NextAminoAcid
        {
            get
            {
                return OneBasedEndResidueInProtein < Protein.Length ? Protein[OneBasedEndResidueInProtein] : '-';
            }
        }

        public Protein Protein
        {
            get { return _protein; }
            protected set { _protein = value; }
        }

        public string BaseSequence
        {
            get
            {
                if (_baseSequence == null)
                {
                    _baseSequence = Protein.BaseSequence.Substring(OneBasedStartResidueInProtein - 1, OneBasedEndResidueInProtein - OneBasedStartResidueInProtein + 1);
                }
                return _baseSequence;
            }
        }

        public char this[int zeroBasedIndex]
        {
            get
            {
                return BaseSequence[zeroBasedIndex];
            }
        }

        /// <summary>
        /// Gets the peptides for a specific protein interval
        /// </summary>
        /// <param name="interval"></param>
        /// <param name="allKnownFixedModifications"></param>
        /// <param name="digestionParams"></param>
        /// <param name="variableModifications"></param>
        /// <returns></returns>
        internal IEnumerable<PeptideWithSetModifications> GetModifiedPeptides(IEnumerable<Modification> allKnownFixedModifications,
            DigestionParams digestionParams, List<Modification> variableModifications)
        {
            int peptideLength = OneBasedEndResidueInProtein - OneBasedStartResidueInProtein + 1;
            int maximumVariableModificationIsoforms = digestionParams.MaxModificationIsoforms;
            int maxModsForPeptide = digestionParams.MaxModsForPeptide;
            var twoBasedPossibleVariableAndLocalizeableModifications = new Dictionary<int, List<Modification>>(peptideLength + 4);

            var pepNTermVariableMods = new List<Modification>();
            twoBasedPossibleVariableAndLocalizeableModifications.Add(1, pepNTermVariableMods);

            var pepCTermVariableMods = new List<Modification>();
            twoBasedPossibleVariableAndLocalizeableModifications.Add(peptideLength + 2, pepCTermVariableMods);

            foreach (Modification variableModification in variableModifications)
            {
                // Check if can be a n-term mod
                if (CanBeNTerminalMod(variableModification, peptideLength))
                {
                    pepNTermVariableMods.Add(variableModification);
                }

                for (int r = 0; r < peptideLength; r++)
                {
                    if (ModificationLocalization.ModFits(variableModification, Protein.BaseSequence, r + 1, peptideLength, OneBasedStartResidueInProtein + r)
                        && variableModification.LocationRestriction == "Anywhere.")
                    {
                        if (!twoBasedPossibleVariableAndLocalizeableModifications.TryGetValue(r + 2, out List<Modification> residueVariableMods))
                        {
                            residueVariableMods = new List<Modification> { variableModification };
                            twoBasedPossibleVariableAndLocalizeableModifications.Add(r + 2, residueVariableMods);
                        }
                        else
                        {
                            residueVariableMods.Add(variableModification);
                        }
                    }
                }
                // Check if can be a c-term mod
                if (CanBeCTerminalMod(variableModification, peptideLength))
                {
                    pepCTermVariableMods.Add(variableModification);
                }
            }

            // LOCALIZED MODS
            foreach (var kvp in Protein.OneBasedPossibleLocalizedModifications)
            {
                bool inBounds = kvp.Key >= OneBasedStartResidueInProtein && kvp.Key <= OneBasedEndResidueInProtein;
                if (!inBounds)
                {
                    continue;
                }

                int locInPeptide = kvp.Key - OneBasedStartResidueInProtein + 1;
                foreach (Modification modWithMass in kvp.Value)
                {
                    if (modWithMass is Modification variableModification)
                    {
                        // Check if can be a n-term mod
                        if (locInPeptide == 1 && CanBeNTerminalMod(variableModification, peptideLength) && !Protein.IsDecoy)
                        {
                            pepNTermVariableMods.Add(variableModification);
                        }

                        int r = locInPeptide - 1;
                        if (r >= 0 && r < peptideLength
                            && (Protein.IsDecoy ||
                            (ModificationLocalization.ModFits(variableModification, Protein.BaseSequence, r + 1, peptideLength, OneBasedStartResidueInProtein + r)
                                && variableModification.LocationRestriction == "Anywhere.")))
                        {
                            if (!twoBasedPossibleVariableAndLocalizeableModifications.TryGetValue(r + 2, out List<Modification> residueVariableMods))
                            {
                                residueVariableMods = new List<Modification> { variableModification };
                                twoBasedPossibleVariableAndLocalizeableModifications.Add(r + 2, residueVariableMods);
                            }
                            else
                            {
                                residueVariableMods.Add(variableModification);
                            }
                        }

                        // Check if can be a c-term mod
                        if (locInPeptide == peptideLength && CanBeCTerminalMod(variableModification, peptideLength) && !Protein.IsDecoy)
                        {
                            pepCTermVariableMods.Add(variableModification);
                        }
                    }
                }
            }

            int variable_modification_isoforms = 0;

            foreach (Dictionary<int, Modification> kvp in GetVariableModificationPatterns(twoBasedPossibleVariableAndLocalizeableModifications, maxModsForPeptide, peptideLength))
            {
                int numFixedMods = 0;
                foreach (var ok in GetFixedModsOneIsNterminus(peptideLength, allKnownFixedModifications))
                {
                    if (!kvp.ContainsKey(ok.Key))
                    {
                        numFixedMods++;
                        kvp.Add(ok.Key, ok.Value);
                    }
                }
                yield return new PeptideWithSetModifications(Protein, digestionParams, OneBasedStartResidueInProtein, OneBasedEndResidueInProtein,
                    PeptideDescription, MissedCleavages, kvp, numFixedMods);
                variable_modification_isoforms++;
                if (variable_modification_isoforms == maximumVariableModificationIsoforms)
                {
                    yield break;
                }
            }
        }

        /// <summary>
        /// Determines whether given modification can be an N-terminal modification
        /// </summary>
        /// <param name="variableModification"></param>
        /// <param name="peptideLength"></param>
        /// <returns></returns>
        private bool CanBeNTerminalMod(Modification variableModification, int peptideLength)
        {
            return ModificationLocalization.ModFits(variableModification, Protein.BaseSequence, 1, peptideLength, OneBasedStartResidueInProtein)
                && (variableModification.LocationRestriction == "N-terminal." || variableModification.LocationRestriction == "Peptide N-terminal.");
        }

        /// <summary>
        /// Determines whether given modification can be a C-terminal modification
        /// </summary>
        /// <param name="variableModification"></param>
        /// <param name="peptideLength"></param>
        /// <returns></returns>
        private bool CanBeCTerminalMod(Modification variableModification, int peptideLength)
        {
            return ModificationLocalization.ModFits(variableModification, Protein.BaseSequence, peptideLength, peptideLength, OneBasedStartResidueInProtein + peptideLength - 1)
                && (variableModification.LocationRestriction == "C-terminal." || variableModification.LocationRestriction == "Peptide C-terminal.");
        }

        private static IEnumerable<Dictionary<int, Modification>> GetVariableModificationPatterns(Dictionary<int, List<Modification>> possibleVariableModifications, int maxModsForPeptide, int peptideLength)
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

        private Dictionary<int, Modification> GetFixedModsOneIsNterminus(int peptideLength,
            IEnumerable<Modification> allKnownFixedModifications)
        {
            var fixedModsOneIsNterminus = new Dictionary<int, Modification>(peptideLength + 3);
            foreach (Modification mod in allKnownFixedModifications)
            {
                switch (mod.LocationRestriction)
                {
                    case "N-terminal.":
                    case "Peptide N-terminal.":
                        if (ModificationLocalization.ModFits(mod, Protein.BaseSequence, 1, peptideLength, OneBasedStartResidueInProtein))
                        {
                            fixedModsOneIsNterminus[1] = mod;
                        }
                        break;

                    case "Anywhere.":
                        for (int i = 2; i <= peptideLength + 1; i++)
                        {
                            if (ModificationLocalization.ModFits(mod, Protein.BaseSequence, i - 1, peptideLength, OneBasedStartResidueInProtein + i - 2))
                            {
                                fixedModsOneIsNterminus[i] = mod;
                            }
                        }
                        break;

                    case "C-terminal.":
                    case "Peptide C-terminal.":
                        if (ModificationLocalization.ModFits(mod, Protein.BaseSequence, peptideLength, peptideLength, OneBasedStartResidueInProtein + peptideLength - 1))
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
    }
}