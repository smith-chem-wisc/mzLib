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

        internal ProteolyticPeptide(Protein protein, int oneBasedStartResidueInProtein, int oneBasedEndResidueInProtein, int missedCleavages, CleavageSpecificity cleavageSpecificityForFdrCategory, string peptideDescription = null)
        {
            _protein = protein;
            OneBasedStartResidueInProtein = oneBasedStartResidueInProtein;
            OneBasedEndResidueInProtein = oneBasedEndResidueInProtein;
            MissedCleavages = missedCleavages;
            CleavageSpecificityForFdrCategory = cleavageSpecificityForFdrCategory;
            PeptideDescription = peptideDescription;
        }

        [NonSerialized] private Protein _protein; // protein that this peptide is a digestion product of
        public int OneBasedStartResidueInProtein { get; } // the residue number at which the peptide begins (the first residue in a protein is 1)
        public int OneBasedEndResidueInProtein { get; } // the residue number at which the peptide ends
        public int MissedCleavages { get; } // the number of missed cleavages this peptide has with respect to the digesting protease
        public string PeptideDescription { get; internal set; } //unstructured explanation of source
        public CleavageSpecificity CleavageSpecificityForFdrCategory { get; internal set; } //structured explanation of source
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
            var possibleVariableAndLocalizeableModifications = new Dictionary<int, List<Modification>>(); // possible mods at each index in peptide

            for (int r = -1; r <= peptideLength; r++)
            {
                possibleVariableAndLocalizeableModifications.Add(r + 1, new List<Modification>());
            }

            // Haven't changed this part much
            // FIXME: VARIABLE MODS
            foreach (Modification variableModification in variableModifications)
            {
                if (CanBeNTerminalMod(variableModification, peptideLength))
                {
                    possibleVariableAndLocalizeableModifications[0].Add(variableModification);
                }
                else if (CanBeCTerminalMod(variableModification, peptideLength))
                {
                    possibleVariableAndLocalizeableModifications[peptideLength + 1].Add(variableModification);
                }
                else
                {
                    for (int r = 0; r < peptideLength; ++r)
                    {
                        if (variableModification.LocationRestriction == "Anywhere."
                            && ModificationLocalization.ModFits(variableModification, Protein.BaseSequence, r + 1, peptideLength, OneBasedStartResidueInProtein + r))
                        {
                            possibleVariableAndLocalizeableModifications[r + 1].Add(variableModification);
                        }
                    }
                }
            }

            // FIXME: LOCALIZED MODS
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
                            possibleVariableAndLocalizeableModifications[0].Add(variableModification);
                        }
                        // check if can be a c-term mod
                        else if (locInPeptide == peptideLength && CanBeCTerminalMod(variableModification, peptideLength) && !Protein.IsDecoy)
                        {
                            possibleVariableAndLocalizeableModifications[peptideLength + 1].Add(variableModification);
                        }
                        else
                        {
                            int r = locInPeptide - 1;
                            if (r >= 0 && r < peptideLength && (Protein.IsDecoy || (variableModification.LocationRestriction == "Anywhere."
                                && ModificationLocalization.ModFits(variableModification, Protein.BaseSequence, r + 1, peptideLength, OneBasedStartResidueInProtein + r))))
                            {
                                possibleVariableAndLocalizeableModifications[r + 1].Add(variableModification);
                            }
                        }
                    }
                }
            }

            int variable_modification_isoforms = 0;
            var allFixedMods = GetFixedModsOneIsNterminus(peptideLength, allKnownFixedModifications);
            int totalAvailableMods = possibleVariableAndLocalizeableModifications.Sum(b => b.Value == null ? 0 : b.Value.Count);

            foreach (Dictionary<int, Modification> kvp in GetModificationPatterns(peptideLength, Math.Min(maxModsForPeptide, totalAvailableMods), possibleVariableAndLocalizeableModifications))
            {
                int numFixedMods = 0;
                allFixedMods.Where(ok => !kvp.ContainsKey(ok.Key)).ToList().ForEach(ok =>
                {
                    numFixedMods++;
                    kvp.Add(ok.Key, ok.Value);
                });

                yield return new PeptideWithSetModifications(Protein, digestionParams, OneBasedStartResidueInProtein, OneBasedEndResidueInProtein,
                    CleavageSpecificityForFdrCategory, PeptideDescription, MissedCleavages, kvp, numFixedMods);

                variable_modification_isoforms++;
                if (variable_modification_isoforms == maximumVariableModificationIsoforms)
                {
                    yield break;
                }
            }
        }

        // gets all modification patterns for a peptide
        private IEnumerable<Dictionary<int, Modification>> GetModificationPatterns(int peptideLength, int numMods, Dictionary<int, List<Modification>> allMods)
        {
            for (int i = 0; i <= numMods; ++i)
            {
                foreach (var modPattern in GeneratePatterns(peptideLength + 2, i, allMods))
                {
                    yield return modPattern; 
                }
            }
        }

        private List<Dictionary<int, Modification>> GeneratePatterns(int n, int k, Dictionary<int, List<Modification>> possibleMods)
        {
            var allPatterns = new List<Dictionary<int, Modification>>();

            if (n <= 0 || n < k)
            {
                return allPatterns;
            }

            var pattern = new Dictionary<int, Modification>();
            GeneratePatterns(n, k, 0, pattern, allPatterns, possibleMods);

            return allPatterns;
        }

        // a recursive method that generates all unique modification patterns for a peptide
        private void GeneratePatterns(int n, int desiredNumMods, int start, Dictionary<int, Modification> pattern,
                List<Dictionary<int, Modification>> allPatterns, Dictionary<int, List<Modification>> possibleMods)
        {
            if (pattern.Count == desiredNumMods)
            {
                allPatterns.Add(new Dictionary<int, Modification>(pattern));
                return;
            }

            for (int i = start; i < n; i++) // n is peptideLength + 2
            {
                int j = 0;
                while (j < possibleMods[i].Count)
                {
                    if (i == 0) // n-term mod
                    {
                        pattern.Add(1, possibleMods[i][j]);
                    }
                    else if (i + 1 == n) // c-term mod
                    {
                        pattern.Add(n, possibleMods[i][j]);
                    }
                    else
                    {
                        pattern.Add(i + 1, possibleMods[i][j]);
                    }

                    GeneratePatterns(n, desiredNumMods, i + 1, pattern, allPatterns, possibleMods); 
                    pattern.Remove(pattern.Keys.Last());
                    ++j;
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