using System;
using MzLibUtil;
using Omics.Modifications;

namespace Omics.Digestion
{
    public abstract class DigestionProduct
    {
        protected static readonly DictionaryPool<int, SortedSet<Modification>> DictionaryPool = new();
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

        #region Cleavage Requirement Discharge

        /// <summary>
        /// True when this peptidoform describes a digestion <paramref name="agent"/> could NOT have
        /// performed, because a modification one of its motifs REQUIRES is absent from the subsite that
        /// motif names. The mirror of a cleavage-blocking drop: that one removes a peptidoform the agent
        /// could not produce BECAUSE of a modification, this one removes a peptidoform it could not
        /// produce WITHOUT one.
        /// </summary>
        /// <param name="variableModPattern">
        /// The modifications this peptidoform carries, in the two-based key scheme the pattern generator
        /// mints: key 1 is the N-terminus, key 2 the FIRST residue, key <paramref name="productLength"/>
        /// + 1 the LAST residue, key + 2 the C-terminus. The parent residue behind key k is therefore
        /// <c>OneBasedStartResidue + k - 2</c>.
        /// </param>
        /// <param name="productLength">
        /// Passed rather than read from <see cref="Length"/>, because the concrete callers already hold
        /// the cheap arithmetic form while <see cref="Length"/> walks a substring of the parent.
        /// </param>
        /// <remarks>
        /// <para><b>Each cut is answered by exactly one of the two products that touch it.</b> A subsite
        /// on the non-prime side lies in the product ENDING at the cut; one on the prime side lies in the
        /// product STARTING at it. So StcE (P2) is discharged by the peptide whose C-terminus is the cut,
        /// and the OgpA family (P1') by the peptide whose N-terminus is. Across a whole digest every
        /// internal cut is therefore checked once, from the side that can see it.</para>
        ///
        /// <para><b>The far end is checked for feasibility, not occupancy.</b> A product still has to
        /// answer for the cut at its other terminus, whose constrained residue lies in the neighbouring
        /// product and is not in this pattern -- the pattern generator discards modifications outside the
        /// product outright. There the parent's list of possible localized modifications decides: if no
        /// modification of the required class can EVER sit at that residue, the cut is impossible and the
        /// peptidoform goes; if one can, the peptidoform is kept even though it may not itself carry one.
        /// That is deliberately one-sided. It never drops a real peptide and it leaves some impossible
        /// ones standing, which is the safe direction to be wrong in.</para>
        ///
        /// <para><b>A cut is justified if ANY motif justifies it.</b> An agent may mix motifs that carry a
        /// requirement with motifs that do not -- StcE-trypsin is exactly that -- so a tryptic cut is
        /// answered by the tryptic motif and needs no glycan, while a cut only the StcE motif explains
        /// does. Each motif is re-matched against the parent sequence at the cut, so a motif that does not
        /// fit there cannot vouch for it.</para>
        ///
        /// <para>Protein termini are not cuts and are never checked: a peptide starting at residue 1, or
        /// ending at the last residue, got that terminus from the sequence ending, not from the agent -- and
        /// neither is a truncation-product boundary the database annotates.</para>
        /// </remarks>
        protected bool IsUnreachableWithoutRequiredModification(Dictionary<int, Modification> variableModPattern,
            int productLength, DigestionAgent agent, IEnumerable<Modification> configuredModifications,
            IEnumerable<Modification> fixedModifications, IReadOnlyList<int> internalFeasibleSites,
            int maxMissedCleavagesAllowed, out int openMissedCleavages)
        {
            openMissedCleavages = MissedCleavages;

            // Nothing configured can require anything, so nothing can be unreachable. This gate is what
            // keeps an ordinary tryptic digest from paying for a feature it cannot use.
            if (agent is null || !agent.HasCleavageRequirement || Parent is null)
            {
                return false;
            }

            string parentSequence = Parent.BaseSequence;

            // Removing the initiator methionine is not a proteolytic cut, so the N-terminus of a peptide
            // that starts at residue 2 of a sequence beginning with Met was not produced by the agent and
            // must not be asked to justify itself. Without this exemption every initiator-cleaved form is
            // deleted -- for a glycoprotease that is HALF of all N-terminal peptidoforms, dropped silently
            // and with no read-through to replace them.
            bool startedAtInitiatorMethionineRemoval = OneBasedStartResidue == 2
                && parentSequence.Length > 0
                && parentSequence[0] == 'M';

            // A cut severs the bond AFTER the residue that names it, so the cut at this peptide's
            // N-terminus falls after the residue preceding it. A terminus the database annotates as a
            // truncation-product boundary (signal peptide, propeptide, chain) is a processing site, not
            // a cut this agent made, and is exempt for the same reason as initiator-Met removal.
            if (OneBasedStartResidue > 1
                && !startedAtInitiatorMethionineRemoval
                && !IsTruncationProductBoundary(OneBasedStartResidue - 1)
                && !AnyMotifJustifies(OneBasedStartResidue - 1, parentSequence, variableModPattern, productLength, agent,
                    configuredModifications, fixedModifications))
            {
                return true;
            }

            if (OneBasedEndResidue < parentSequence.Length
                && !IsTruncationProductBoundary(OneBasedEndResidue)
                && !AnyMotifJustifies(OneBasedEndResidue, parentSequence, variableModPattern, productLength, agent,
                    configuredModifications, fixedModifications))
            {
                return true;
            }

            // Both ends are real cuts, so this peptide exists. What remains is whether it exists WITHIN
            // the caller's missed-cleavage budget, and that is not the span's raw count: a feasible site
            // this peptidoform leaves unoccupied is a site the agent could not have cut, so skipping it
            // is not a missed cleavage. Discounting those is what lets the read-through across an
            // unoccupied site come back at the budget the caller actually set -- the peptide the drop
            // above used to remove with nothing to replace it.
            //
            // Clamped, so the count can never go negative, and compared afterwards so the generation
            // slack that bought this span cannot leak out as a peptide claiming more missed cleavages
            // than were asked for. Exactly the shape of the blocking mirror in ProteolyticPeptide.
            if (internalFeasibleSites is not null && internalFeasibleSites.Count > 0)
            {
                int unjustifiedInternalSites = 0;
                for (int i = 0; i < internalFeasibleSites.Count; i++)
                {
                    int site = internalFeasibleSites[i];
                    if (!AnyMotifJustifies(site, parentSequence, variableModPattern, productLength, agent,
                            configuredModifications, fixedModifications))
                    {
                        unjustifiedInternalSites++;
                    }
                }

                openMissedCleavages = MissedCleavages - Math.Min(unjustifiedInternalSites, MissedCleavages);
            }

            return openMissedCleavages > maxMissedCleavagesAllowed;
        }

        /// <summary>
        /// The internal positions of this product at which the agent could have cut -- the sites that
        /// survive <see cref="DigestionAgent.FilterToFeasibleCleavageSites"/> and so were in the list the
        /// span was enumerated from. One-based parent residues, each naming the bond after it.
        /// </summary>
        /// <remarks>
        /// Hoisted out of the peptidoform loop on purpose. Which positions are SITES depends only on the
        /// sequence and the configured modifications, not on where a particular peptidoform happens to put
        /// its glycans, so scanning per pattern would repeat identical work inside the hottest loop in the
        /// library. Only the occupancy question is per-peptidoform.
        /// </remarks>
        protected List<int> FindInternalFeasibleCleavageSites(DigestionAgent agent,
            IEnumerable<Modification> configuredModifications)
        {
            var sites = new List<int>();
            if (agent is null || !agent.HasCleavageRequirement || Parent is null)
            {
                return sites;
            }

            for (int residue = OneBasedStartResidue; residue < OneBasedEndResidue; residue++)
            {
                if (agent.IsFeasibleCleavageSite(residue, Parent, configuredModifications))
                {
                    sites.Add(residue);
                }
            }

            return sites;
        }

        /// <summary>
        /// True when at least one of the agent's motifs both MATCHES the parent sequence at this cut and
        /// has its modification requirement met there. A motif carrying no requirement justifies any cut
        /// it matches.
        /// </summary>
        /// <param name="cutAfterOneBasedResidue">
        /// The parent residue whose C-side bond is severed. The motif's recognition sequence therefore
        /// begins <see cref="DigestionMotif.CutIndex"/> residues earlier.
        /// </param>
        private bool AnyMotifJustifies(int cutAfterOneBasedResidue, string parentSequence,
            Dictionary<int, Modification> variableModPattern, int productLength, DigestionAgent agent,
            IEnumerable<Modification> configuredModifications, IEnumerable<Modification> fixedModifications)
        {
            foreach (DigestionMotif motif in agent.DigestionMotifs)
            {
                if (motif is null)
                {
                    continue;
                }

                // Where the motif's recognition sequence would have to start for its cut to land here.
                // Fits takes a ZERO-based index into the sequence it is handed, and that sequence is the
                // PARENT here -- the same API is used against a peptide's own sequence elsewhere in the
                // library, so the coordinate space has to be stated rather than assumed.
                int motifStartZeroBased = cutAfterOneBasedResidue - motif.CutIndex;
                if (motifStartZeroBased < 0
                    || motifStartZeroBased + motif.InducingCleavage.Length > parentSequence.Length)
                {
                    continue;
                }

                (bool fits, bool prevented) = motif.Fits(parentSequence, motifStartZeroBased);
                if (!fits || prevented)
                {
                    continue;
                }

                // ALL of this motif's conditions must hold for it to vouch for the cut. A motif with
                // none is satisfied by sequence alone.
                bool everyConditionMet = true;
                foreach (CleavageRequirement requirement in motif.CleavageRequirements)
                {
                    if (!RequirementIsMet(requirement, cutAfterOneBasedResidue, variableModPattern, productLength,
                            configuredModifications, fixedModifications))
                    {
                        everyConditionMet = false;
                        break;
                    }
                }

                if (everyConditionMet)
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// True when this condition holds at the residue it addresses: a REQUIRED modification present,
        /// or a FORBIDDEN one absent.
        /// </summary>
        /// <remarks>
        /// The two polarities fail in opposite directions, so they are given opposite benefit of the
        /// doubt wherever occupancy is unknowable. A required condition whose subsite lies in the
        /// neighbouring product falls back to whether the parent COULD carry the modification; a
        /// forbidden one is simply treated as met, because the alternative is to drop a peptide over a
        /// residue this peptidoform does not describe. Both choices only ever keep peptides.
        /// </remarks>
        private bool RequirementIsMet(CleavageRequirement requirement, int cutAfterOneBasedResidue,
            Dictionary<int, Modification> variableModPattern, int productLength,
            IEnumerable<Modification> configuredModifications, IEnumerable<Modification> fixedModifications)
        {
            // Subsites count outward from the severed bond: P1 is the residue before it, P1' the one
            // after. The bond falls after cutAfterOneBasedResidue, so Pk is (cut - k + 1) and Pk' is
            // (cut + k), both one-based in the parent.
            int constrainedResidue = requirement.IsPrimeSide
                ? cutAfterOneBasedResidue + requirement.Subsite
                : cutAfterOneBasedResidue - requirement.Subsite + 1;

            // No such residue: a forbidden modification cannot be there, so the condition holds
            // vacuously, while a required one can never be satisfied.
            if (constrainedResidue < 1 || constrainedResidue > Parent.BaseSequence.Length)
            {
                return requirement.IsForbidden;
            }

            // Inside this product, the pattern is authoritative: it says what this peptidoform carries.
            if (constrainedResidue >= OneBasedStartResidue && constrainedResidue <= OneBasedEndResidue)
            {
                int key = constrainedResidue - OneBasedStartResidue + 2;
                if (key < 2 || key > productLength + 1)
                {
                    return requirement.IsForbidden;
                }

                bool carriesMatchingModification = variableModPattern is not null
                    && variableModPattern.TryGetValue(key, out Modification placed)
                    && requirement.IsSatisfiedBy(placed);

                return requirement.IsConditionMet(carriesMatchingModification);
            }

            // Outside this product the peptidoform says nothing about the residue, so a forbidden
            // condition is normally given the benefit of the doubt -- guessing would delete a real peptide.
            //
            // A FIXED modification is the one case where the doubt is settled. Fixed means every copy of
            // the residue carries it, so if one satisfies this condition and fits there, the subsite IS
            // occupied and the cleavage was impossible. That is what lets a non-prime forbidden rule prune
            // from BOTH sides of its bond: the product ENDING at the cut sees the subsite in its own
            // pattern, and the product STARTING at it can only ever see it here. Without this, IMPa's
            // "no cleavage between two adjacent glycosites" removes the left-hand peptide and leaves the
            // right-hand one standing, which is half a rule.
            if (requirement.IsForbidden)
            {
                if (fixedModifications is null)
                {
                    return true;
                }

                string sequenceForFixed = Parent.BaseSequence;
                foreach (Modification fixedModification in fixedModifications)
                {
                    if (requirement.IsSatisfiedBy(fixedModification)
                        && ModificationLocalization.ModFits(fixedModification, sequenceForFixed, constrainedResidue,
                            sequenceForFixed.Length, constrainedResidue))
                    {
                        return false;
                    }
                }

                return true;
            }

            // Outside this product, in a neighbour: fall back to whether the parent could ever carry a
            // satisfying modification there. See the caller's remarks for why this is feasibility rather
            // than occupancy.
            //
            // BOTH sources are consulted, for the same reason the site-list filter consults both. A search
            // supplying the glycan as a variable modification against an unannotated database has nothing
            // in OneBasedPossibleLocalizedModifications, so an annotation-only test answers "impossible"
            // for every cut whose constrained residue sits in the neighbouring product -- and silently
            // deletes every peptide lying C-terminal to a glycan-justified cut.
            if (Parent.OneBasedPossibleLocalizedModifications is not null
                && Parent.OneBasedPossibleLocalizedModifications.TryGetValue(constrainedResidue, out var candidates)
                && candidates is not null
                && candidates.Any(requirement.IsSatisfiedBy))
            {
                return true;
            }

            if (configuredModifications is null)
            {
                return false;
            }

            string parentSequence = Parent.BaseSequence;
            foreach (Modification configured in configuredModifications)
            {
                if (requirement.IsSatisfiedBy(configured)
                    && ModificationLocalization.ModFits(configured, parentSequence, constrainedResidue,
                        parentSequence.Length, constrainedResidue))
                {
                    return true;
                }
            }

            return false;
        }


        /// <summary>
        /// The residues inside this product that MUST carry a modification satisfying a cleavage
        /// requirement, because the cuts that produced this product could not have happened otherwise.
        /// Empty for every ordinary protease. Keys are two-based, matching the modification-pattern
        /// convention: key 2 is the first residue, key <see cref="Length"/> + 1 the last.
        /// </summary>
        /// <remarks>
        /// <para><b>What this is for, and why it cannot be answered at digestion.</b> A glyco search never
        /// puts the glycan in the digestion: the peptide backbone is identified naked and the glycan
        /// arrives afterwards as a precursor-mass difference, to be localized against fragment ions. So
        /// the cleavage requirement cannot filter peptides there the way it does for a search that
        /// configures the glycan as a modification. What it CAN do is tell the localizer where the glycan
        /// must have been: if OpeRATOR produced this peptide, its first residue carries a glycan. That is
        /// not a hypothesis to be scored, it is a consequence of the peptide existing at all.</para>
        ///
        /// <para><b>Derived, never stored.</b> Everything needed is already reachable -- the agent's
        /// motifs, this product's coordinates and the parent sequence -- so this computes on demand and
        /// adds no field. That is not a style preference: MetaMorpheus caches the peptide index to disk
        /// with a field-based, schema-rigid serializer whose reuse check carries no assembly version, so a
        /// new serialized field silently changes the on-disk layout and corrupts indices users already
        /// have.</para>
        ///
        /// <para><b>Deliberately conservative in two places, because over-constraining DELETES correct
        /// answers.</b> A cut that any requirement-free motif explains is not obligated at all -- a
        /// tryptic cut inside a StcE-trypsin digest needs no glycan. And when several requirement-carrying
        /// motifs fit the same cut but point at DIFFERENT residues, the truth is a disjunction ("one of
        /// these must be occupied") which a set of forced positions cannot express, so nothing is emitted.
        /// Under-constraining only forgoes a refinement; over-constraining would forbid the right
        /// localization.</para>
        ///
        /// <para>A subsite lying outside this product is skipped: it belongs to the neighbouring peptide
        /// and is not a position this product can localize to. For a P2 requirement (StcE) the C-terminal
        /// cut's subsite is inside and the N-terminal cut's is not; for P1' (the OgpA family) it is the
        /// other way round.</para>
        /// </remarks>
        /// <summary>
        /// The obligated sites of <see cref="GetCleavageObligatedSites"/>, each with the conditions the
        /// protease imposes there -- so a caller can ask not only WHERE a modification must be but WHICH
        /// modifications would satisfy the rule.
        /// </summary>
        /// <remarks>
        /// This is what lets a glyco search use a composition-level rule. OpeRATOR does not merely require
        /// "a glycan" at P1'; it requires at least core 1, and is blocked again by core 2. The localizer
        /// holds real glycans, and a glycan IS a Modification, so it can put each candidate through
        /// <see cref="CleavageRequirement.IsSatisfiedBy"/> directly and refuse the arrangements that place
        /// the wrong KIND of glycan on a site the enzyme has already spoken for.
        /// </remarks>
        public Dictionary<int, List<CleavageRequirement>> GetCleavageObligations(DigestionAgent agent)
        {
            var obligations = new Dictionary<int, List<CleavageRequirement>>();
            foreach (int site in GetCleavageObligatedSites(agent, obligations))
            {
                // GetCleavageObligatedSites fills the dictionary as it goes; the loop is only to run it.
                _ = site;
            }

            return obligations;
        }

        public List<int> GetCleavageObligatedSites(DigestionAgent agent) =>
            GetCleavageObligatedSites(agent, null);

        private List<int> GetCleavageObligatedSites(DigestionAgent agent,
            Dictionary<int, List<CleavageRequirement>> conditionsBySite)
        {
            var obligated = new List<int>();
            if (agent is null || !agent.HasCleavageRequirement || Parent is null)
            {
                return obligated;
            }

            string parentSequence = Parent.BaseSequence;
            int productLength = OneBasedEndResidue - OneBasedStartResidue + 1;

            // Removing the initiator methionine is not a proteolytic cut, and neither is the sequence
            // simply beginning or ending, nor a truncation-product boundary the database annotates --
            // none of those needs a modification to explain it.
            bool startedAtInitiatorMethionineRemoval = OneBasedStartResidue == 2
                && parentSequence.Length > 0
                && parentSequence[0] == 'M';

            if (OneBasedStartResidue > 1 && !startedAtInitiatorMethionineRemoval
                && !IsTruncationProductBoundary(OneBasedStartResidue - 1))
            {
                AddObligationForCut(OneBasedStartResidue - 1, parentSequence, agent, productLength, obligated, conditionsBySite);
            }

            if (OneBasedEndResidue < parentSequence.Length && !IsTruncationProductBoundary(OneBasedEndResidue))
            {
                AddObligationForCut(OneBasedEndResidue, parentSequence, agent, productLength, obligated, conditionsBySite);
            }

            return obligated;
        }

        /// <summary>
        /// Adds the site this cut forces, if it forces exactly one and that one lies inside this product.
        /// </summary>
        private void AddObligationForCut(int cutAfterOneBasedResidue, string parentSequence, DigestionAgent agent,
            int productLength, List<int> obligated, Dictionary<int, List<CleavageRequirement>> conditionsBySite)
        {
            int forcedResidue = -1;

            // Each condition with the parent residue it addresses. A motif's conditions need not share a
            // residue -- IMPa is "P1':O-glycan;P1:!O-glycan", two residues either side of the bond -- so
            // only those addressing the forced residue may be reported against it.
            var conditionsHere = new List<(int Residue, CleavageRequirement Condition)>();

            foreach (DigestionMotif motif in agent.DigestionMotifs)
            {
                if (motif is null)
                {
                    continue;
                }

                int motifStartZeroBased = cutAfterOneBasedResidue - motif.CutIndex;
                if (motifStartZeroBased < 0
                    || motifStartZeroBased + motif.InducingCleavage.Length > parentSequence.Length)
                {
                    continue;
                }

                (bool fits, bool prevented) = motif.Fits(parentSequence, motifStartZeroBased);
                if (!fits || prevented)
                {
                    continue;
                }

                // Only a REQUIRED condition proves a residue was occupied. A forbidden one proves the
                // opposite and is not an obligation to localize anything, so it is skipped -- emitting it
                // would force a glycan onto the one residue the enzyme says cannot carry it. It is still
                // collected, because once a site IS obliged, a forbidden condition AT THAT SITE still
                // says what may not sit there.
                var required = new List<CleavageRequirement>();
                foreach (CleavageRequirement candidate in motif.CleavageRequirements)
                {
                    if (candidate.IsForbidden)
                    {
                        conditionsHere.Add((ConstrainedResidue(candidate, cutAfterOneBasedResidue), candidate));
                    }
                    else
                    {
                        required.Add(candidate);
                    }
                }

                // This motif explains the cut without requiring any modification, so the cut obliges
                // nothing at all.
                if (required.Count != 1)
                {
                    return;
                }

                CleavageRequirement requirement = required[0];
                int constrainedResidue = ConstrainedResidue(requirement, cutAfterOneBasedResidue);
                conditionsHere.Add((constrainedResidue, requirement));

                if (forcedResidue == -1)
                {
                    forcedResidue = constrainedResidue;
                }
                else if (forcedResidue != constrainedResidue)
                {
                    // Two motifs fit and disagree about which residue must be occupied. The obligation is
                    // a disjunction, which a forced-site set cannot say, so say nothing.
                    return;
                }
            }

            if (forcedResidue < OneBasedStartResidue || forcedResidue > OneBasedEndResidue)
            {
                return;
            }

            int twoBasedKey = forcedResidue - OneBasedStartResidue + 2;
            if (twoBasedKey >= 2 && twoBasedKey <= productLength + 1 && !obligated.Contains(twoBasedKey))
            {
                obligated.Add(twoBasedKey);

                if (conditionsBySite is not null)
                {
                    // Every condition, of every motif that vouched for this cut, that addresses THIS
                    // residue, so the caller can test a candidate modification against all of them.
                    // Forbidden conditions are carried too -- they do not oblige a site, but once a site
                    // IS obliged they still say what may not sit on it. A condition addressing another
                    // residue is left out: IMPa's P1 prohibition would otherwise refuse the very glycan
                    // its P1' requirement demands.
                    conditionsBySite[twoBasedKey] = conditionsHere
                        .Where(c => c.Residue == forcedResidue)
                        .Select(c => c.Condition)
                        .ToList();
                }
            }
        }

        /// <summary>
        /// True when the bond after <paramref name="cutAfterOneBasedResidue"/> is a boundary of one of the
        /// parent's annotated truncation products. Digestion yields products ending or starting there
        /// ("chain start", "chain end" and the intact products), so such a terminus came from the
        /// database, not from the agent.
        /// </summary>
        /// <remarks>
        /// When a truncation boundary happens to coincide with a real cleavage site this exempts that
        /// cut too, which only ever keeps a peptide -- the safe direction for this correction to err.
        /// </remarks>
        private bool IsTruncationProductBoundary(int cutAfterOneBasedResidue)
        {
            if (Parent?.TruncationProducts is null)
            {
                return false;
            }

            foreach (var truncationProduct in Parent.TruncationProducts)
            {
                if (truncationProduct.OneBasedBeginPosition - 1 == cutAfterOneBasedResidue
                    || truncationProduct.OneBasedEndPosition == cutAfterOneBasedResidue)
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// The one-based parent residue a condition addresses at the cut after
        /// <paramref name="cutAfterOneBasedResidue"/>. Subsites count outward from the severed bond: Pk is
        /// (cut - k + 1) and Pk' is (cut + k).
        /// </summary>
        private static int ConstrainedResidue(CleavageRequirement requirement, int cutAfterOneBasedResidue) =>
            requirement.IsPrimeSide
                ? cutAfterOneBasedResidue + requirement.Subsite
                : cutAfterOneBasedResidue - requirement.Subsite + 1;

        #endregion

        #region Digestion Helper Methods

        /// <summary>
        /// Generates all possible variable modification patterns for a peptide, which includes variable and localized modifications but excludes fixed mods
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
        protected static IEnumerable<Dictionary<int, Modification>> GetVariableModificationPatterns(Dictionary<int, SortedSet<Modification>> possibleVariableModifications, int maxModsForPeptide, int peptideLength)
        {
            if (possibleVariableModifications.Count <= 0) 
                yield break;

            int[] baseVariableModificationPattern = new int[peptideLength + 4];
            int totalAvailableMods = possibleVariableModifications.Values.Sum(modList => modList?.Count ?? 0);
            int maxVariableMods = Math.Min(totalAvailableMods, maxModsForPeptide);
            var variableModKvpList = possibleVariableModifications.ToList();

            for (int variable_modifications = 0; variable_modifications <= maxVariableMods; variable_modifications++)
            {
                foreach (int[] variable_modification_pattern in GetVariableModificationPatternsRecursive(variableModKvpList,
                             possibleVariableModifications.Count - variable_modifications, baseVariableModificationPattern, 0))
                {
                    // use modification pattern to construct a dictionary of modifications for the peptide
                    var modificationPattern = new Dictionary<int, Modification>(possibleVariableModifications.Count);

                    foreach (var variableModSet in possibleVariableModifications)
                    {
                        int modIndex = variable_modification_pattern[variableModSet.Key] - 1;
                        if (modIndex >= 0)
                        {
                            modificationPattern.Add(variableModSet.Key, variableModSet.Value.ElementAt(modIndex));
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
            IEnumerable<Modification> allKnownFixedModifications, in Dictionary<int, Modification> fixedModsOneIsNterminus)
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
                            else if (mod is CleavageModification)
                            {
                                if (OneBasedStartResidue != 1)
                                    fixedModsOneIsNterminus[1] = mod;
                            }
                            else if (OneBasedStartResidue == 1) // Modified BioPolymer Start Residue (e.g. Protein N-Terminal)
                            {
                                if (!fixedModsOneIsNterminus.TryAdd(1, mod)) // Check if a protein N-terminal mod is already present
                                {
                                    if (mod.LocationRestriction is "N-terminal." or "5'-terminal.") // Only overwrite if new mod is N-terminal, not peptide N-terminal
                                    {
                                        fixedModsOneIsNterminus[1] = mod;
                                    }
                                }
                            }
                            else //Normal N-terminal peptide modification
                            {
                                fixedModsOneIsNterminus[1] = mod;
                            }
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
                            else if (mod is CleavageModification)
                            {
                                if (OneBasedEndResidue != Parent.Length)
                                    fixedModsOneIsNterminus[length + 2] = mod;
                            }
                            else if (OneBasedEndResidue == Parent.Length) // Modified BioPolymer End Residue (e.g. Protein C-Terminal)
                            {
                                if (!fixedModsOneIsNterminus.TryAdd(length + 2, mod)) // Check if a protein C-terminal mod is already present
                                {
                                    if (mod.LocationRestriction is "C-terminal." or "3'-terminal.") // Only overwrite if new mod is C-terminal, not peptide C-terminal
                                    {
                                        fixedModsOneIsNterminus[length + 2] = mod;
                                    }
                                }
                            }
                            else //Normal C-terminal peptide modification 
                            {
                                fixedModsOneIsNterminus[length + 2] = mod;
                            }
                        }
                        break;

                    default:
                        throw new NotSupportedException("This terminus localization is not supported.");
                }
            }

            foreach (var entry in Parent.OneBasedFixedModifications)
            {
                if (entry.Key == 0 && OneBasedStartResidue == 1)
                    fixedModsOneIsNterminus[1] = entry.Value;
                else if (entry.Key == Parent.Length + 2 && OneBasedEndResidue == Parent.Length)
                    fixedModsOneIsNterminus[length + 2] = entry.Value;
                else if (entry.Key >= OneBasedStartResidue && entry.Key <= OneBasedEndResidue)
                    fixedModsOneIsNterminus[entry.Key - OneBasedStartResidue + 2] = entry.Value;
            }
        }

        /// <summary>
        /// Populates the variable modifications dictionary  from both the variable modifications and the localized mods from xml reading, 
        /// considering the N-terminal, C-terminal, and internal positions.
        /// </summary>
        /// <param name="allVariableMods">A list of all variable modifications.</param>
        /// <param name="twoBasedDictToPopulate">A reference to a dictionary that will hold the variable modifications, with the key representing the position.</param>
        /// <remarks>
        /// This method iterates through all variable modifications and assigns them to the appropriate positions in the peptide.
        /// It considers different location restrictions such as N-terminal, C-terminal, and anywhere within the peptide.
        /// </remarks>
        protected void PopulateVariableModifications(List<Modification> allVariableMods, in Dictionary<int, SortedSet<Modification>> twoBasedDictToPopulate)
        {
            int peptideLength = OneBasedEndResidue - OneBasedStartResidue + 1;
            var pepNTermVariableMods = new SortedSet<Modification>();
            twoBasedDictToPopulate.Add(1, pepNTermVariableMods);

            var pepCTermVariableMods = new SortedSet<Modification>();
            twoBasedDictToPopulate.Add(peptideLength + 2, pepCTermVariableMods);

            // VARIABLE MODS
            foreach (Modification variableModification in allVariableMods)
            {
                // Check if can be a n-term mod
                if (CanBeNTerminalOrFivePrime(variableModification, peptideLength)
                    && (!IsCleavageModification(variableModification) || OneBasedStartResidue != 1)
                    && !ModificationLocalization.UniprotModExists(Parent, 1, variableModification))
                {
                    pepNTermVariableMods.Add(variableModification);
                }

                for (int r = 0; r < peptideLength; r++)
                {
                    if (ModificationLocalization.ModFits(variableModification, Parent.BaseSequence, r + 1, peptideLength, OneBasedStartResidue + r)
                        && variableModification.LocationRestriction == "Anywhere." && !ModificationLocalization.UniprotModExists(Parent, r + 1, variableModification))
                    {
                        if (!twoBasedDictToPopulate.TryGetValue(r + 2, out var residueVariableMods))
                        {
                            residueVariableMods = new SortedSet<Modification>() { variableModification };
                            twoBasedDictToPopulate.Add(r + 2, residueVariableMods);
                        }
                        else
                        {
                            residueVariableMods.Add(variableModification);
                        }
                    }
                }
                // Check if can be a c-term mod
                if (CanBeCTerminalOrThreePrime(variableModification, peptideLength)
                    && (!IsCleavageModification(variableModification) || OneBasedEndResidue != Parent.Length)
                    && !ModificationLocalization.UniprotModExists(Parent, peptideLength, variableModification))
                {
                    pepCTermVariableMods.Add(variableModification);
                }
            }

            // LOCALIZED MODS
            foreach (var kvp in Parent.OneBasedPossibleLocalizedModifications)
            {
                bool inBounds = kvp.Key >= OneBasedStartResidue && kvp.Key <= OneBasedEndResidue;
                if (!inBounds)
                {
                    continue;
                }

                int locInPeptide = kvp.Key - OneBasedStartResidue + 1;
                foreach (Modification modWithMass in kvp.Value)
                {
                    if (modWithMass is not Modification variableModification)
                        continue;

                    // Check if can be a n-term mod
                    if (locInPeptide == 1 && CanBeNTerminalOrFivePrime(variableModification, peptideLength) && !Parent.IsDecoy)
                    {
                        pepNTermVariableMods.Add(variableModification);
                    }

                    int r = locInPeptide - 1;
                    if (r >= 0 && r < peptideLength
                               && (Parent.IsDecoy ||
                                   (ModificationLocalization.ModFits(variableModification, Parent.BaseSequence, r + 1, peptideLength, OneBasedStartResidue + r)
                                    && variableModification.LocationRestriction == "Anywhere.")))
                    {
                        if (!twoBasedDictToPopulate.TryGetValue(r + 2, out var residueVariableMods))
                        {
                            residueVariableMods = new SortedSet<Modification>() { variableModification };
                            twoBasedDictToPopulate.Add(r + 2, residueVariableMods);
                        }
                        else
                        {
                            residueVariableMods.Add(variableModification);
                        }
                    }

                    // Check if can be a c-term mod
                    if (locInPeptide == peptideLength && CanBeCTerminalOrThreePrime(variableModification, peptideLength) && !Parent.IsDecoy)
                    {
                        pepCTermVariableMods.Add(variableModification);
                    }
                }
            }
        }

        /// <summary>
        /// Appends fixed modifications to the variable modification pattern when no variable mod exists. 
        /// </summary>
        /// <param name="fixedModDict">The dictionary containing fixed modifications.</param>
        /// <param name="variableModPattern">The dictionary containing the variable modification pattern.</param>
        /// <param name="numFixedMods">The number of fixed modifications appended.</param>
        /// <remarks>
        /// This method iterates through the fixed modifications and adds them to the variable modification pattern
        /// if they are not already present. The number of fixed modifications appended is returned via the out parameter.
        /// </remarks>
        protected void AppendFixedModificationsToVariable(in Dictionary<int, Modification> fixedModDict, in Dictionary<int, Modification> variableModPattern, out int numFixedMods)
        {
            numFixedMods = 0;
            foreach (var fixedModPattern in fixedModDict)
            {
                if (variableModPattern.ContainsKey(fixedModPattern.Key))
                    continue;
                numFixedMods++;
                variableModPattern.Add(fixedModPattern.Key, fixedModPattern.Value);
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
        private static IEnumerable<int[]> GetVariableModificationPatternsRecursive(List<KeyValuePair<int, SortedSet<Modification>>> possibleVariableModifications,
            int unmodifiedResiduesDesired, int[] variableModificationPattern, int index)
        {
            if (index < possibleVariableModifications.Count - 1)
            {
                if (unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    foreach (int[] new_variable_modification_pattern in GetVariableModificationPatternsRecursive(possibleVariableModifications,
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
                        foreach (int[] new_variable_modification_pattern in GetVariableModificationPatternsRecursive(possibleVariableModifications,
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
        private bool CanBeNTerminalOrFivePrime(Modification mod, int peptideLength)
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
        private bool CanBeCTerminalOrThreePrime(Modification mod, int peptideLength)
        {
            return mod.LocationRestriction is "3'-terminal." or "Oligo 3'-terminal." or "C-terminal." or "Peptide C-terminal."
                   && ModificationLocalization.ModFits(mod, Parent.BaseSequence, peptideLength, peptideLength, OneBasedStartResidue + peptideLength - 1);
        }

        private static bool IsCleavageModification(Modification modification)
        {
            return modification.ModificationType == "Protease" || modification is CleavageModification;
        }

        #endregion
    }
}
