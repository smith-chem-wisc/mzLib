using System;
using System.Collections.Generic;
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
        internal IEnumerable<PeptideWithSetModifications> GetModifiedPeptides(List<Modification> allKnownFixedModifications,
            DigestionParams digestionParams, List<Modification> variableModifications)
        {
            int variable_modification_isoforms = 0;
            int peptideLength = OneBasedEndResidue - OneBasedStartResidue + 1;
            int maximumVariableModificationIsoforms = digestionParams.MaxModificationIsoforms;
            int maxModsForPeptide = digestionParams.MaxModsForPeptide;
            var twoBasedPossibleVariableAndLocalizeableModifications = DictionaryPool.Get();
            var fixedModDictionary = FixedModDictionaryPool.Get();

            try
            {
                PopulateVariableModifications(variableModifications, in twoBasedPossibleVariableAndLocalizeableModifications);
                PopulateFixedModsOneIsNorFivePrimeTerminus(peptideLength, allKnownFixedModifications, in fixedModDictionary);

                foreach (Dictionary<int, Modification> variableModPattern in GetVariableModificationPatterns(twoBasedPossibleVariableAndLocalizeableModifications, maxModsForPeptide, peptideLength))
                {
                    AppendFixedModificationsToVariable(in fixedModDictionary, in variableModPattern, out int numFixedMods);

                    yield return new PeptideWithSetModifications(Protein, digestionParams, OneBasedStartResidue, OneBasedEndResidue,
                        CleavageSpecificityForFdrCategory, PeptideDescription, MissedCleavages, variableModPattern, numFixedMods);

                    variable_modification_isoforms++;
                    if (variable_modification_isoforms == maximumVariableModificationIsoforms)
                    {
                        yield break;
                    }
                }
            }
            finally
            {
                FixedModDictionaryPool.Return(fixedModDictionary);
                DictionaryPool.Return(twoBasedPossibleVariableAndLocalizeableModifications);
            }
        }
    }
}