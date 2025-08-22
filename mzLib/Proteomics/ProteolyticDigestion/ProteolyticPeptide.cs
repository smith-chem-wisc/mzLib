using System;
using System.Collections.Generic;
using System.Linq;
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
        internal IEnumerable<PeptideWithSetModifications> GetModifiedPeptides(List<Modification> allKnownFixedModifications, //for each peptide, peptide3.getvariablemods? 
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
        /// <summary>
        /// Calculates the number of possible proteoforms for a given peptide sequence. Both fixed and variable modifications are considered.
        /// </summary>
        /// <param name="allKnownFixedModifications">Number of Fixed Modifications</param>
        /// <param name="digestionParams">Used to ascertain the user chosen max mods per peptide</param>
        /// <param name="variableModifications">Number of Variable Modifications</param>
        /// <returns></returns>
        public long GetNumberOfPossibleProteoforms(List<Modification> allKnownFixedModifications,
            DigestionParams digestionParams, List<Modification> variableModifications)
        {
            int peptideLength = OneBasedEndResidue - OneBasedStartResidue + 1;
            var twoBasedPossibleVariableAndLocalizeableModifications = DictionaryPool.Get();
            PopulateVariableModifications(variableModifications, in twoBasedPossibleVariableAndLocalizeableModifications);
            List<int> modsPerResidue = new List<int>();  //list of possible mods at residue, without attention to the order. {2,6,3,4,5}
            foreach (var x in twoBasedPossibleVariableAndLocalizeableModifications)
            {
                //iterate through each residue, add in the number of possible modifications for that residue.
                modsPerResidue.Add(x.Value.Count);
            }
            //given a maximum number of mods, the upper bound on proteoforms we can have with a base sequence.
            return DigestionProduct.PossibleProteoformRecursive(digestionParams.MaxMods, modsPerResidue); 
        }
    }
}