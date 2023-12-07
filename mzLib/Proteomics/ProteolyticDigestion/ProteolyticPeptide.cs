using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using Omics.Digestion;
using Omics.Modifications;

namespace Proteomics.ProteolyticDigestion
{
    /// <summary>
    /// Product of digesting a protein
    /// Contains methods for modified peptide combinitorics
    /// </summary>
    [Serializable]
    public class ProteolyticPeptide : DigestionProduct
    {
        
        internal ProteolyticPeptide(Protein protein, int oneBasedStartResidueInProtein, int oneBasedEndResidueInProtein, int missedCleavages, CleavageSpecificity cleavageSpecificityForFdrCategory, string peptideDescription = null, string baseSequence = null) :
            base(protein, oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, cleavageSpecificityForFdrCategory, peptideDescription, baseSequence)
        {

        }

        
        public Protein Protein
        {
            get => Parent as Protein;
            protected set => Parent = value;
        }

        #region Properties overridden by more generic interface

        public int OneBasedEndResidueInProtein => OneBasedEndResidue;
        public int OneBasedStartResidueInProtein => OneBasedStartResidue;
        public virtual char PreviousAminoAcid => PreviousResidue;
        public virtual char NextAminoAcid => NextResidue;

        public string PeptideDescription
        {
            get => Description;
            set => Description = value;
        }

        #endregion

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
            int peptideLength = OneBasedEndResidue - OneBasedStartResidue + 1;
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
                if (CanBeNTerminalMod(variableModification, peptideLength) && !ModificationLocalization.UniprotModExists(Protein, 1, variableModification))
                {
                    pepNTermVariableMods.Add(variableModification);
                }

                for (int r = 0; r < peptideLength; r++)
                {
                    if (ModificationLocalization.ModFits(variableModification, Protein.BaseSequence, r + 1, peptideLength, OneBasedStartResidue + r)
                        && variableModification.LocationRestriction == "Anywhere." && !ModificationLocalization.UniprotModExists(Protein, r + 1, variableModification))
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
                if (CanBeCTerminalMod(variableModification, peptideLength) && !ModificationLocalization.UniprotModExists(Protein, peptideLength, variableModification))
                {
                    pepCTermVariableMods.Add(variableModification);
                }
            }

            // LOCALIZED MODS
            foreach (var kvp in Protein.OneBasedPossibleLocalizedModifications)
            {
                bool inBounds = kvp.Key >= OneBasedStartResidue && kvp.Key <= OneBasedEndResidue;
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
                foreach (var ok in GetFixedModsOneIsNorFivePrimeTerminus(peptideLength, allKnownFixedModifications))
                {
                    if (!kvp.ContainsKey(ok.Key))
                    {
                        numFixedMods++;
                        kvp.Add(ok.Key, ok.Value);
                    }
                }
                yield return new PeptideWithSetModifications(Protein, digestionParams, OneBasedStartResidue, OneBasedEndResidue,
                    CleavageSpecificityForFdrCategory, PeptideDescription, MissedCleavages, kvp, numFixedMods);
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
            return ModificationLocalization.ModFits(variableModification, Protein.BaseSequence, 1, peptideLength, OneBasedStartResidue)
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
            return ModificationLocalization.ModFits(variableModification, Protein.BaseSequence, peptideLength, peptideLength, OneBasedStartResidue + peptideLength - 1)
                && (variableModification.LocationRestriction == "C-terminal." || variableModification.LocationRestriction == "Peptide C-terminal.");
        }
    }
}