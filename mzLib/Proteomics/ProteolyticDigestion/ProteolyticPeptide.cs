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

                    // Modifications are placed AFTER cleavage, so nothing so far has checked whether a
                    // modification abolishes the very site this peptide was cut at. Skip the peptidoforms
                    // the protease could not have produced.
                    if (digestionParams.RespectCleavageBlockingModifications
                        && IsUnreachableThroughBlockedCleavage(variableModPattern, peptideLength, digestionParams.MaxMissedCleavages))
                    {
                        continue;
                    }

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
        /// True when this peptidoform describes a digestion the protease could not have performed,
        /// given where its cleavage-blocking modifications landed.
        ///
        /// Two ways that happens:
        /// (1) The C-terminal residue carries a cleavage-blocking modification and is an internal cut.
        ///     Trypsin cannot cleave after an acylated lysine, so this peptidoform -- typically reported
        ///     with zero missed cleavages -- describes an event that does not occur. Drop it.
        /// (2) The peptidoform only exists because of the generation slack (see
        ///     <see cref="DigestionParams.CleavageBlockingReadThroughSlack"/>) and, once blocked sites
        ///     are discounted, still has more OPEN missed cleavages than the caller allowed. A blocked
        ///     internal residue is not a cleavage site for this peptidoform, so it must not be counted
        ///     as a missed cleavage -- which is exactly what lets the read-through form of a blocked
        ///     cleavage survive at MaxMissedCleavages = 0.
        /// </summary>
        /// <remarks>
        /// Position keys follow the two-based scheme used by the modification pattern: key 1 is the
        /// peptide N-terminus, residue i (1-based) is key i + 1, and key peptideLength + 2 is the
        /// peptide C-terminus. So the C-terminal RESIDUE is key peptideLength + 1, and the residues
        /// whose C-side bond is an internal missed cleavage (positions 1 .. peptideLength - 1) are
        /// keys 2 .. peptideLength.
        ///
        /// Known approximation: a blocking modification on an internal residue is counted as covering
        /// an internal cleavage site without re-deriving the protease's site list, so a modified K or R
        /// that was never a site (trypsin's K|P rule, say) can discount a missed cleavage it did not
        /// occupy. The count is clamped so it can never go negative, and case (1) -- the correctness
        /// fix this is here for -- is exact.
        /// </remarks>
        private bool IsUnreachableThroughBlockedCleavage(Dictionary<int, Modification> variableModPattern,
            int peptideLength, int maxMissedCleavagesAllowed)
        {
            bool cTerminalResidueBlocked = false;
            int blockedInternalSites = 0;

            foreach (KeyValuePair<int, Modification> positionAndMod in variableModPattern)
            {
                if (positionAndMod.Value is null || !positionAndMod.Value.BlocksCleavage)
                {
                    continue;
                }

                if (positionAndMod.Key == peptideLength + 1)
                {
                    cTerminalResidueBlocked = true;
                }
                else if (positionAndMod.Key >= 2 && positionAndMod.Key <= peptideLength)
                {
                    blockedInternalSites++;
                }
            }

            // A peptide ending at the protein's own C-terminus was not produced by a cleavage, so no
            // modification there can invalidate it.
            bool cTerminusIsAnInternalCut = OneBasedEndResidue < Protein.Length;
            if (cTerminalResidueBlocked && cTerminusIsAnInternalCut)
            {
                return true;
            }

            int openMissedCleavages = MissedCleavages - Math.Min(blockedInternalSites, MissedCleavages);
            return openMissedCleavages > maxMissedCleavagesAllowed;
        }
    }
}