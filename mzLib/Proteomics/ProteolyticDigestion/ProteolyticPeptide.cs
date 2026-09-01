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
                    //
                    // Both gates are required, and they are not the same gate.
                    //
                    // SearchModeType == Full is what ties the drop to the generation slack, which
                    // ProteinDigestion adds under exactly that condition. The drop and the slack are two
                    // halves of one exchange: the impossible peptidoform goes out and the read-through
                    // form that replaces it comes in. Applying the drop where no slack was granted
                    // performs half the trade -- a semi search at a missed-cleavage budget too small to
                    // reach the read-through would lose the impossible form AND never gain the real one,
                    // leaving that peptide unidentifiable. So a semi or nonspecific search is left
                    // entirely alone rather than half-corrected.
                    //
                    // CleavageSpecificityForFdrCategory == Full then restricts the drop to peptides whose
                    // C-terminus is a protease cut at all: for a semi or single-terminus peptide it is a
                    // length-driven truncation, and a blocking modification there invalidates nothing.
                    //
                    // reportedMissedCleavages is what the surviving peptidoform CARRIES, and it is the
                    // blocked-discounted count, not the modification-blind one this peptide was
                    // enumerated under. A blocked residue is not a cleavage site for this peptidoform,
                    // so counting it as a missed cleavage would report a cleavage that cannot occur --
                    // and would let the generation slack leak out as peptides claiming more missed
                    // cleavages than the caller asked for.
                    int reportedMissedCleavages = MissedCleavages;
                    if (digestionParams.RespectCleavageBlockingModifications
                        && digestionParams.SearchModeType == CleavageSpecificity.Full
                        && CleavageSpecificityForFdrCategory == CleavageSpecificity.Full
                        && IsUnreachableThroughBlockedCleavage(variableModPattern, peptideLength, digestionParams.MaxMissedCleavages,
                            digestionParams.Protease, out reportedMissedCleavages))
                    {
                        continue;
                    }

                    yield return new PeptideWithSetModifications(Protein, digestionParams, OneBasedStartResidue, OneBasedEndResidue,
                        CleavageSpecificityForFdrCategory, PeptideDescription, reportedMissedCleavages, variableModPattern, numFixedMods);

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
        /// (2) The peptidoform only exists because of the generation slack that
        ///     <see cref="ProteinDigestion"/> adds when this flag is on, and -- once blocked sites are
        ///     discounted -- still has more OPEN missed cleavages than the caller allowed. A blocked
        ///     internal residue is not a cleavage site for this peptidoform, so it must not be counted
        ///     as a missed cleavage -- which is exactly what lets the read-through form of a blocked
        ///     cleavage survive at MaxMissedCleavages = 0.
        /// </summary>
        /// <param name="openMissedCleavages">
        /// The missed cleavages this peptidoform actually has once blocked residues are discounted --
        /// the count a surviving peptidoform must REPORT. Always assigned. For any peptidoform that
        /// survives, this is guaranteed to be &lt;= the caller's MaxMissedCleavages, so the generation
        /// slack stays an enumeration detail and never reaches the caller as an out-of-budget count.
        /// </param>
        /// <remarks>
        /// Position keys follow the two-based scheme used by the modification pattern: key 1 is the
        /// peptide N-terminus, residue i (1-based) is key i + 1, and key peptideLength + 2 is the
        /// peptide C-terminus. So the C-terminal RESIDUE is key peptideLength + 1, and the residues
        /// whose C-side bond is an internal missed cleavage (positions 1 .. peptideLength - 1) are
        /// keys 2 .. peptideLength.
        ///
        /// The protease is consulted, so a modification only blocks a cleavage this protease would
        /// have made: an acylated lysine counts under trypsin, Lys-C, elastase or subtilisin, and does
        /// not under Glu-C, Asp-N, Lys-N or chymotrypsin, none of which cut after a lysine. A protease
        /// that cuts N-TERMINAL to its recognition residue (Asp-N, Lys-N) is left inert rather than
        /// half-corrected: there the blocked residue invalidates a peptide's N-terminus, which is the
        /// mirror image of the drop below and is not modelled here.
        ///
        /// Known approximation, now narrowed to sequence context: cleavage sites are matched at the
        /// residue level without re-deriving the protease's site list for this protein, so a modified
        /// K or R whose site was prevented by surrounding sequence (trypsin|P's K[P] rule) can still
        /// discount a missed cleavage it did not occupy. The count is clamped so it can never go
        /// negative, and case (1) -- the correctness fix this is here for -- is exact.
        /// </remarks>
        private bool IsUnreachableThroughBlockedCleavage(Dictionary<int, Modification> variableModPattern,
            int peptideLength, int maxMissedCleavagesAllowed, DigestionAgent protease, out int openMissedCleavages)
        {
            bool cTerminalResidueBlocked = false;
            int blockedInternalSites = 0;

            foreach (KeyValuePair<int, Modification> positionAndMod in variableModPattern)
            {
                // BlocksCleavageBy, not BlocksCleavage: a modification blocks a cleavage only if this
                // protease was going to make one there. An acetylated lysine abolishes a trypsin site
                // and abolishes nothing in a Glu-C digest, where discounting it would hide a missed
                // cleavage the peptide genuinely has.
                if (!CleavageBlockingModifications.BlocksCleavageBy(positionAndMod.Value, protease))
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

            // Assigned before any early return: the caller reads this on the surviving path, and a
            // dropped peptidoform's count is simply never used.
            openMissedCleavages = MissedCleavages - Math.Min(blockedInternalSites, MissedCleavages);

            // A peptide ending at the protein's own C-terminus was not produced by a cleavage, so no
            // modification there can invalidate it.
            bool cTerminusIsAnInternalCut = OneBasedEndResidue < Protein.Length;
            if (cTerminalResidueBlocked && cTerminusIsAnInternalCut)
            {
                return true;
            }

            return openMissedCleavages > maxMissedCleavagesAllowed;
        }
    }
}