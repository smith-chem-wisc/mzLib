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
        /// <param name="allKnownFixedModifications"></param>
        /// <param name="digestionParams"></param>
        /// <param name="variableModifications"></param>
        /// <param name="cleavageBlockingPolicy">
        /// The other half of the exchange <see cref="Protein.Digest(Omics.Digestion.IDigestionParams,List{Modification},List{Modification},List{SilacLabel},System.ValueTuple{SilacLabel,SilacLabel}?,bool)"/>
        /// bought generation slack for. Passed in rather than rebuilt here, so the drop below and that
        /// slack are decided by one object and cannot disagree about the protease or the budget. The
        /// default is the inert policy, which is exactly the historical behaviour.
        /// </param>
        /// <returns></returns>
        internal IEnumerable<PeptideWithSetModifications> GetModifiedPeptides(List<Modification> allKnownFixedModifications,
            DigestionParams digestionParams, List<Modification> variableModifications,
            CleavageBlockingPolicy cleavageBlockingPolicy = default)
        {
            int variable_modification_isoforms = 0;
            int peptideLength = OneBasedEndResidue - OneBasedStartResidue + 1;
            int maximumVariableModificationIsoforms = digestionParams.MaxModificationIsoforms;
            int maxModsForPeptide = digestionParams.MaxModsForPeptide;

            // Hoisted out of the pattern loop: none of these can change between peptidoforms of the same
            // peptide, and the loop below is the ~8.8-billion-call hot path.
            //
            // There is deliberately no "can anything satisfy the requirement" gate -- see Protein.Digest
            // for why the promoting correction must NOT go inert when nothing configured can satisfy it.
            List<Modification> configuredModifications = (variableModifications ?? Enumerable.Empty<Modification>())
                .Concat(allKnownFixedModifications ?? Enumerable.Empty<Modification>())
                .ToList();
            bool respectCleavageRequirements = digestionParams.RespectCleavagePromotingModifications
                && digestionParams.SearchModeType == CleavageSpecificity.Full
                && CleavageSpecificityForFdrCategory == CleavageSpecificity.Full
                && digestionParams.DigestionAgent is not null
                && digestionParams.DigestionAgent.HasCleavageRequirement;

            // Which internal positions are SITES is a property of the sequence and the configured
            // modifications, identical for every peptidoform of this peptide, so it is found once here
            // rather than rescanned inside the pattern loop. Only occupancy varies per peptidoform.
            List<int> internalFeasibleSites = respectCleavageRequirements
                ? FindInternalFeasibleCleavageSites(digestionParams.DigestionAgent, configuredModifications)
                : null;

            // When both corrections are on, the blocking drop hands the residues it discounted to the
            // promoting drop, which applies the ONE missed-cleavage budget over both. Reused across
            // peptidoforms rather than allocated per pattern.
            List<int> blockedInternalResidues = respectCleavageRequirements && cleavageBlockingPolicy.IsActive
                ? new List<int>()
                : null;
            var twoBasedPossibleVariableAndLocalizeableModifications = DictionaryPool.Get();
            var fixedModDictionary = FixedModDictionaryPool.Get();

            try
            {
                PopulateVariableModifications(variableModifications, in twoBasedPossibleVariableAndLocalizeableModifications);
                PopulateFixedModsOneIsNorFivePrimeTerminus(peptideLength, allKnownFixedModifications, in fixedModDictionary);

                foreach (Dictionary<int, Modification> variableModPattern in GetVariableModificationPatterns(twoBasedPossibleVariableAndLocalizeableModifications, maxModsForPeptide, peptideLength))
                {
                    // Modifications are placed AFTER cleavage, so nothing so far has checked whether a
                    // modification abolishes the very site this peptide was cut at. Skip the peptidoforms
                    // the protease could not have produced, and report what the survivors really carry.
                    //
                    // Asked BEFORE the fixed modifications are merged in, so the pattern read here holds
                    // only the VARIABLE modifications. That is deliberate and is the scope the policy
                    // itself defines: a fixed blocking modification buys no slack, so it must not be able
                    // to spend any either. See CleavageBlockingPolicy.For.
                    //
                    // CleavageSpecificityForFdrCategory == Full is the per-peptide half of the gate the
                    // policy cannot know: it restricts the drop to peptides whose C-terminus is a protease
                    // cut at all. For a semi or single-terminus peptide it is a length-driven truncation,
                    // and a blocking modification there invalidates nothing.
                    //
                    // reportedMissedCleavages is what the surviving peptidoform CARRIES, and it is the
                    // blocked-discounted count, not the modification-blind one this peptide was
                    // enumerated under. A blocked residue is not a cleavage site for this peptidoform,
                    // so counting it as a missed cleavage would report a cleavage that cannot occur --
                    // and would let the generation slack leak out as peptides claiming more missed
                    // cleavages than the caller asked for.
                    //
                    // With the promoting correction also on, the budget is NOT tested here: a read-through
                    // may be over budget by blocked sites alone and back under it once the unoccupied
                    // glycoprotease sites are discounted too, so the blocked residues are passed on and the
                    // promoting drop below tests the budget once, over both.
                    int reportedMissedCleavages = MissedCleavages;
                    blockedInternalResidues?.Clear();
                    if (cleavageBlockingPolicy.IsActive
                        && CleavageSpecificityForFdrCategory == CleavageSpecificity.Full
                        && IsUnreachableThroughBlockedCleavage(variableModPattern, peptideLength, cleavageBlockingPolicy,
                            blockedInternalResidues, out reportedMissedCleavages))
                    {
                        continue;
                    }

                    AppendFixedModificationsToVariable(in fixedModDictionary, in variableModPattern, out int numFixedMods);

                    // The mirror gate: a protease whose motif REQUIRES a modification at one of its
                    // subsites cannot have made a cut where that modification is absent. Skip the
                    // peptidoforms it could not have produced.
                    //
                    // Gated exactly like the blocking drop above, and for the same reason on the second
                    // clause: CleavageSpecificityForFdrCategory == Full restricts this to peptides whose
                    // termini are protease cuts at all, since for a semi or single-terminus peptide a
                    // terminus is a length-driven truncation and no glycan can be expected to justify it.
                    //
                    // Asked AFTER the fixed modifications are merged in, unlike the blocking drop above: a
                    // fixed modification is unavoidable, and that is evidence the promoting rules use (a
                    // fixed modification can satisfy a requirement, and is the one way a forbidden
                    // condition becomes judgeable from the product that starts at the bond).
                    //
                    // This drop refines OCCUPANCY, and it is the second half of a two-stage correction.
                    // DigestionAgent.FullDigestion has already removed, from the site list itself, every site where
                    // the required modification could not be -- so the read-through across an impossible
                    // site is an ordinary peptide here and needed no slack to reach. What survives that
                    // filter is a site that COULD carry the modification; this gate removes the
                    // peptidoforms in which it does not.
                    //
                    // The peptide spanning a feasible but unoccupied site is a real read-through, reached
                    // by the generation slack DigestionAgent.FullDigestion adds. The unoccupied sites it
                    // spans are discounted from its missed cleavages together with any sites the blocking
                    // drop above discounted, and the budget is tested once over both.
                    if (respectCleavageRequirements
                        && IsUnreachableWithoutRequiredModification(variableModPattern, peptideLength,
                            digestionParams.DigestionAgent, configuredModifications, allKnownFixedModifications,
                            internalFeasibleSites, blockedInternalResidues, digestionParams.MaxMissedCleavages,
                            out reportedMissedCleavages))
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
        ///     <see cref="Protein.Digest(Omics.Digestion.IDigestionParams,List{Modification},List{Modification},List{SilacLabel},System.ValueTuple{SilacLabel,SilacLabel}?,bool)"/>
        ///     adds when this policy is active, and -- once blocked sites are
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
        /// The protease is consulted through <paramref name="policy"/>, so a modification only blocks a cleavage this protease would
        /// have made: an acylated lysine counts under trypsin, Lys-C, elastase or subtilisin, and does
        /// not under Glu-C, Asp-N, Lys-N or chymotrypsin, none of which cut after a lysine. A protease
        /// that cuts N-TERMINAL to its recognition residue (Asp-N, Lys-N) is left inert rather than
        /// half-corrected: there the blocked residue invalidates a peptide's N-terminus, which is the
        /// mirror image of the drop below and is not modelled here.
        ///
        /// <paramref name="variableModPattern"/> holds VARIABLE modifications only -- it is read before
        /// the fixed modifications are merged into it. A fixed blocking modification therefore cannot
        /// drop a peptidoform, which is the scope <see cref="CleavageBlockingPolicy.For"/> defines and
        /// explains: no slack computed from MaxMods can pay for the read-through a fixed blocking
        /// modification would require, so allowing it to drop would lose the peptide outright.
        ///
        /// Known approximation, now narrowed to sequence context: cleavage sites are matched at the
        /// residue level without re-deriving the protease's site list for this protein, so a modified
        /// K or R whose site was prevented by surrounding sequence (trypsin|P's K[P] rule) can still
        /// discount a missed cleavage it did not occupy. The count is clamped so it can never go
        /// negative, and case (1) -- the correctness fix this is here for -- is exact.
        /// </remarks>
        /// <param name="blockedInternalResidues">
        /// Null when this drop owns the budget. Otherwise it receives the one-based parent residue of
        /// every blocked internal site, the budget test is skipped, and the caller applies it once over
        /// these and the promoting correction's unoccupied sites together.
        /// </param>
        private bool IsUnreachableThroughBlockedCleavage(Dictionary<int, Modification> variableModPattern,
            int peptideLength, CleavageBlockingPolicy policy, List<int> blockedInternalResidues,
            out int openMissedCleavages)
        {
            bool cTerminalResidueBlocked = false;
            int blockedInternalSites = 0;

            foreach (KeyValuePair<int, Modification> positionAndMod in variableModPattern)
            {
                // The policy, not Modification.BlocksCleavage: a modification blocks a cleavage only if
                // this protease was going to make one there. An acetylated lysine abolishes a trypsin
                // site and abolishes nothing in a Glu-C digest, where discounting it would hide a missed
                // cleavage the peptide genuinely has.
                if (!policy.Blocks(positionAndMod.Value))
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
                    blockedInternalResidues?.Add(OneBasedStartResidue + positionAndMod.Key - 2);
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

            return blockedInternalResidues is null && openMissedCleavages > policy.MaxMissedCleavages;
        }
    }
}