using System;
using System.Collections.Generic;
using System.Linq;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;

namespace Proteomics.ProteolyticDigestion
{
    /// <summary>
    /// The single decision behind the cleavage-blocking correction, so that its two halves cannot
    /// disagree.
    ///
    /// The correction is an exchange, not a filter. A peptidoform whose C-terminus is a cleavage the
    /// protease could not have made goes OUT (<see cref="Proteomics.ProteolyticDigestion.ProteolyticPeptide"/>
    /// drops it), and the read-through form that really would have been produced comes IN -- which it
    /// can only do if digestion was told to enumerate wider in the first place
    /// (<see cref="Protein.Digest(Omics.Digestion.IDigestionParams,List{Modification},List{Modification},List{SilacLabel},System.ValueTuple{SilacLabel,SilacLabel}?,bool)"/>
    /// adds <see cref="GenerationSlack"/>). Performing half the trade loses real peptides, so both
    /// halves are decided here, once, from the same <see cref="DigestionParams"/>, the same
    /// modification list and the same protein.
    /// </summary>
    /// <remarks>
    /// A default-valued instance is INERT: <see cref="IsActive"/> is false, <see cref="GenerationSlack"/>
    /// is zero and <see cref="Blocks"/> always returns false, so flag-off digestion is byte-for-byte the
    /// historical path.
    /// </remarks>
    public readonly struct CleavageBlockingPolicy
    {
        private readonly DigestionAgent _protease;

        /// <summary>
        /// Indexed by one-based residue, 0 through the protein length: true where the protease's site list
        /// has a cut after that residue, the protein's own termini included -- the list FullDigestion
        /// enumerates from. See <see cref="MissedCleavageBudgetFor"/>.
        /// </summary>
        private readonly bool[] _sites;

        /// <summary>
        /// Indexed by one-based residue: true when the protease really cuts after that residue in this
        /// protein AND only because of that residue's side chain, so a blocking modification on it
        /// abolishes the cut. See <see cref="ChargeDirectedSites"/>.
        /// </summary>
        private readonly bool[] _chargeDirectedSites;

        private readonly Protein _protein;

        // The variable and fixed modifications, and the fixed alone; null unless the cleavage-promoting correction is on.
        private readonly List<Modification> _configuredModifications;
        private readonly List<Modification> _fixedModifications;

        private CleavageBlockingPolicy(DigestionAgent protease, int maxMissedCleavages, int generationSlack,
            bool[] sites, bool[] chargeDirectedSites, Protein protein, List<Modification> configuredModifications,
            List<Modification> fixedModifications, List<int> cleavageSites)
        {
            _protease = protease;
            MaxMissedCleavages = maxMissedCleavages;
            GenerationSlack = generationSlack;
            _sites = sites;
            _chargeDirectedSites = chargeDirectedSites;
            _protein = protein;
            _configuredModifications = configuredModifications;
            _fixedModifications = fixedModifications;
            CleavageSites = cleavageSites;
        }

        /// <summary>
        /// The policy in force for this digestion of <paramref name="protein"/>, or the inert policy when
        /// the correction cannot fire. Built once per protein in Protein.Digest and handed to both ends of
        /// the exchange, which is what guarantees the slack and the drop agree about the protease and about
        /// the missed-cleavage budget.
        /// </summary>
        /// <param name="variableModifications">
        /// The VARIABLE modifications configured for the search. Fixed modifications are deliberately not
        /// consulted -- see the remarks.
        /// </param>
        /// <param name="protein">
        /// The protein being digested. The policy is per protein: it records which residues are real
        /// protease sites there, so a blocked residue that was never a site (a K before P under trypsin|P)
        /// cannot cancel a genuine missed cleavage, and it reads the modifications annotated on it.
        /// </param>
        /// <param name="fixedModifications">
        /// Read only by the cleavage-promoting correction, where a fixed glycan meets a requirement.
        /// </param>
        /// <remarks>
        /// Only the modifications digestion places OPTIONALLY activate the policy, and only they can block
        /// a cleavage under it: the configured variable modifications, and the protein's own annotated ones
        /// (<see cref="Protein.OneBasedPossibleLocalizedModifications"/> -- UniProt XML features, what
        /// G-PTM-D writes). Both enter the same per-peptide pattern under the same
        /// <see cref="DigestionParams.MaxMods"/> cap, so the slack below pays for either. An annotated
        /// blocking modification has to activate the policy by itself, or whether a database acetyllysine
        /// is corrected would depend on whether some unrelated blocking modification is in the variable list.
        ///
        /// A fixed blocking modification is unavoidable: it sits on every instance of its target residue,
        /// so the read-through that would replace a dropped peptidoform has to span every such residue in
        /// the protein, and no slack computed from MaxMods -- which does not bound fixed modifications --
        /// can pay for it. At MaxMods 0 the slack is zero while the drop would still fire, and the peptide
        /// is lost rather than corrected. Fixed blocking modifications therefore stay outside this
        /// correction entirely; modelling them means removing their sites from the protease's site list
        /// before enumeration, which is a different mechanism.
        /// </remarks>
        public static CleavageBlockingPolicy For(DigestionParams digestionParams, IEnumerable<Modification> variableModifications,
            Protein protein, IEnumerable<Modification> fixedModifications = null)
        {
            if (digestionParams is null || digestionParams.SearchModeType != CleavageSpecificity.Full)
            {
                return default;
            }

            DigestionAgent protease = digestionParams.Protease;
            bool blocking = digestionParams.RespectCleavageBlockingModifications
                && (AnyCanBlockCleavage(variableModifications, protease) || AnyAnnotatedCanBlockCleavage(protein, protease));
            bool promoting = digestionParams.RespectCleavagePromotingModifications
                && protease.DigestionMotifs.Any(motif => motif.CleavageRequirements.Count > 0);
            if (!blocking && !promoting)
            {
                return default;
            }

            string sequence = protein?.BaseSequence ?? string.Empty;
            List<int> cleavageSites = protease.GetDigestionSiteIndices(sequence);
            bool[] sites = Sites(cleavageSites, sequence.Length);
            bool[] chargeDirectedSites = blocking ? ChargeDirectedSites(protease, sequence, sites) : null;
            int generationSlack = blocking ? digestionParams.MaxMods : 0;
            List<Modification> configured = null;
            List<Modification> fixedMods = null;
            if (promoting)
            {
                // A site no required modification can reach, or where a forbidden one must sit, is not a site, so it leaves
                // the list before enumeration.
                // A site one can reach may still be unoccupied in a given peptidoform: the slack enumerates the
                // read-throughs across those, and the drop discounts them again.
                fixedMods = fixedModifications?.ToList() ?? new List<Modification>();
                configured = (variableModifications ?? Enumerable.Empty<Modification>()).Concat(fixedMods).ToList();
                cleavageSites = cleavageSites.Where(site => site == 0 || site == sequence.Length
                    || Justifies(protease, sequence, site, (requirement, residue) => AssumedToCarry(protein, configured, fixedMods, requirement, residue)))
                    .ToList();
                sites = Sites(cleavageSites, sequence.Length);
                // A read-through may also span blocked sites, which the blocking slack already pays for.
                generationSlack += ReadThroughSlack(protease, sequence, cleavageSites, digestionParams.MaxLength,
                    digestionParams.MaxMissedCleavages + generationSlack);
            }

            return new CleavageBlockingPolicy(protease, digestionParams.MaxMissedCleavages, generationSlack, sites,
                chargeDirectedSites, protein, configured, fixedMods, promoting ? cleavageSites : null);
        }

        private static bool AnyAnnotatedCanBlockCleavage(Protein protein, DigestionAgent protease) =>
            protein?.OneBasedPossibleLocalizedModifications is not null
            && protein.OneBasedPossibleLocalizedModifications.Values.Any(modifications => AnyCanBlockCleavage(modifications, protease));

        /// <summary>
        /// <see cref="DigestionAgent.GetDigestionSiteIndices"/> as a lookup, termini included.
        /// </summary>
        private static bool[] Sites(List<int> cleavageSites, int sequenceLength)
        {
            var sites = new bool[sequenceLength + 1];
            foreach (int site in cleavageSites)
            {
                sites[site] = true;
            }

            return sites;
        }

        /// <summary>
        /// The residues of <paramref name="sequence"/> after which <paramref name="protease"/> really cuts
        /// and could be stopped from cutting by a modification of that residue's side chain.
        /// </summary>
        /// <remarks>
        /// Starts from <paramref name="sites"/>, the same site list digestion enumerates from and counts
        /// missed cleavages against, so the discount and the count agree site for site. The protein's own
        /// termini are not cleavages and are excluded. A site is then removed if any motif that fits there
        /// does not name the residue literally at P1 -- StcE's "TX|T" cuts "TK|T" whatever the lysine
        /// carries, so in a StcE-trypsin co-digest that site survives an acetylated K even though trypsin's
        /// "K|" alone would not.
        /// </remarks>
        private static bool[] ChargeDirectedSites(DigestionAgent protease, string sequence, bool[] sites)
        {
            var chargeDirected = (bool[])sites.Clone();
            chargeDirected[0] = false;
            chargeDirected[sequence.Length] = false;

            for (int location = 0; location < sequence.Length; location++)
            {
                foreach (DigestionMotif motif in protease.DigestionMotifs)
                {
                    int site = location + motif.CutIndex;
                    if (site > 0 && site < sequence.Length && chargeDirected[site]
                        && !motif.CleavesCTerminalTo(sequence[site - 1])
                        && motif.Fits(sequence, location).Item1)
                    {
                        chargeDirected[site] = false;
                    }
                }
            }

            return chargeDirected;
        }

        /// <summary>
        /// True when at least one of <paramref name="variableModifications"/> could abolish a cleavage
        /// <paramref name="protease"/> would otherwise have performed. Exposed so the gate itself can be
        /// asserted: the enumeration it guards is trimmed again on the way out, so an output comparison
        /// passes whether or not the gate works.
        /// </summary>
        public static bool AnyCanBlockCleavage(IEnumerable<Modification> variableModifications, DigestionAgent protease) =>
            variableModifications is not null
            && variableModifications.Any(modification => CleavageBlockingModifications.BlocksCleavageBy(modification, protease));

        /// <summary>
        /// False for the inert policy, which is what flag-off digestion and any search whose variable
        /// modifications cannot block a cleavage of the configured protease both get.
        /// </summary>
        public bool IsActive => _protease is not null;

        /// <summary>
        /// Extra missed cleavages digestion must enumerate so that the read-through form of a blocked
        /// cleavage exists to replace the peptidoform this policy drops. Zero when inert, so a search
        /// that cannot use the slack does not pay for it.
        /// </summary>
        /// <remarks>
        /// <see cref="DigestionParams.MaxMods"/>: a peptidoform carries at most that many variable and
        /// annotated modifications and therefore at most that many blocked sites, so every read-through
        /// this policy can need is reachable, however many blocked sites co-occur.
        /// </remarks>
        public int GenerationSlack { get; }

        /// <summary>
        /// The budget the CALLER asked for. A surviving peptidoform's reported missed-cleavage count is
        /// bounded by it, through <see cref="MissedCleavageBudgetFor"/>, so the slack above stays an
        /// enumeration detail and is never observable in the result.
        /// </summary>
        public int MaxMissedCleavages { get; }

        /// <summary>
        /// The most open missed cleavages a surviving peptidoform spanning residues
        /// <paramref name="oneBasedStartResidue"/> to <paramref name="oneBasedEndResidue"/> may report: the
        /// caller's <see cref="MaxMissedCleavages"/>, or one fewer for a WHOLE proteolysis product.
        /// </summary>
        /// <remarks>
        /// One fewer because that is what flag-off digestion does. FullDigestion emits a whole proteolysis
        /// product -- a chain, signal peptide or truncation emitted entire, so that neither end is a
        /// cleavage site -- only while the sites inside it are strictly FEWER than the budget, and reports
        /// the sites inside it; every other product it emits may span as many sites as the budget. The
        /// slack enumerates under a larger budget and this bound trims back to the caller's, so it has to
        /// apply the same bound to each kind of product. With one bound for all, a whole product spanning
        /// exactly the budget's worth of sites came out with the flag on and not with it off, and the slack
        /// was visible in the result. If FullDigestion's bound changes, this must change with it;
        /// WholeProteolysisProducts_AreTrimmedToTheBudgetFlagOffDigestionUses asserts both together.
        ///
        /// A blocked site is not a site, so a whole product whose only sites are blocked counts as spanning
        /// none, and comes out exactly when flag-off digestion would emit a whole product spanning none.
        /// </remarks>
        public int MissedCleavageBudgetFor(int oneBasedStartResidue, int oneBasedEndResidue)
        {
            bool wholeProteolysisProduct = IsActive
                && oneBasedStartResidue >= 1 && oneBasedEndResidue < _sites.Length
                && !_sites[oneBasedStartResidue - 1] && !_sites[oneBasedEndResidue];
            return wholeProteolysisProduct ? MaxMissedCleavages - 1 : MaxMissedCleavages;
        }

        /// <summary>
        /// True when <paramref name="modification"/>, sitting on one-based residue
        /// <paramref name="oneBasedResidue"/> of the protein this policy was built for, abolishes a
        /// cleavage the configured protease would otherwise have performed there. Always false for the
        /// inert policy, and false wherever the protease never cut -- so the discount counts real sites only.
        /// </summary>
        public bool Blocks(Modification modification, int oneBasedResidue) =>
            _chargeDirectedSites is not null
            && oneBasedResidue > 0 && oneBasedResidue < _chargeDirectedSites.Length
            && _chargeDirectedSites[oneBasedResidue]
            && CleavageBlockingModifications.BlocksCleavageBy(modification, _protease);

        /// <summary>
        /// The sites FullDigestion enumerates when the cleavage-promoting correction is on: the protease's, less those
        /// no required modification can reach or a forbidden one must sit at. Null otherwise, and FullDigestion uses the protease's own.
        /// </summary>
        public List<int> CleavageSites { get; }

        /// <summary>False when a terminus of this peptidoform is a cut it lacks the required modification for.</summary>
        internal bool CutsHaveRequiredModifications(int start, int end, Dictionary<int, Modification> variableMods,
            Dictionary<int, Modification> fixedMods)
        {
            if (_configuredModifications is null)
            {
                return true;
            }

            bool nTerminusIsACut = start > 1 && !(start == 2 && _protein[0] == 'M') && !IsTruncationBoundary(start - 1);
            bool cTerminusIsACut = end < _protein.Length && !IsTruncationBoundary(end);
            return (!nTerminusIsACut || Justifies(start - 1, start, end, variableMods, fixedMods))
                && (!cTerminusIsACut || Justifies(end, start, end, variableMods, fixedMods));
        }

        /// <summary>The sites inside this peptidoform that lack their required modification, so are not missed cleavages.</summary>
        internal int UnoccupiedSites(int start, int end, Dictionary<int, Modification> variableMods,
            Dictionary<int, Modification> fixedMods)
        {
            int unoccupied = 0;
            for (int site = start; _configuredModifications is not null && site < end; site++)
            {
                if (_sites[site] && !Justifies(site, start, end, variableMods, fixedMods))
                {
                    unoccupied++;
                }
            }

            return unoccupied;
        }

        // Inside the peptidoform a subsite carries what the peptidoform says it carries.
        private bool Justifies(int cut, int start, int end, Dictionary<int, Modification> variableMods,
            Dictionary<int, Modification> fixedMods)
        {
            Protein protein = _protein;
            List<Modification> configured = _configuredModifications;
            List<Modification> fixedModifications = _fixedModifications;
            return Justifies(_protease, protein.BaseSequence, cut, (requirement, residue) =>
                residue >= start && residue <= end
                    ? (variableMods.TryGetValue(residue - start + 2, out Modification placed)
                        || fixedMods.TryGetValue(residue - start + 2, out placed))
                      && requirement.IsSatisfiedBy(placed)
                    : AssumedToCarry(protein, configured, fixedModifications, requirement, residue));
        }

        // True if a motif that cuts here has every glycan it requires and none it forbids.
        private static bool Justifies(DigestionAgent protease, string sequence, int cut, Func<CleavageRequirement, int, bool> carries) =>
            protease.DigestionMotifs.Any(motif => cut >= motif.CutIndex
                && motif.Fits(sequence, cut - motif.CutIndex) is (true, false)
                && motif.CleavageRequirements.All(requirement =>
                    carries(requirement, requirement.ConstrainedResidue(cut)) != requirement.IsForbidden));

        // Where no peptidoform says, assume what lets the protease cut: a required glycan wherever one could be, a forbidden
        // one only where a permanent or fixed one must be.
        private static bool AssumedToCarry(Protein protein, List<Modification> configured, List<Modification> fixedModifications,
            CleavageRequirement requirement, int residue) =>
            residue >= 1 && residue <= protein.Length
            && ((protein.OneBasedFixedModifications.TryGetValue(residue, out Modification permanent) && requirement.IsSatisfiedBy(permanent))
                || (!requirement.IsForbidden && protein.OneBasedPossibleLocalizedModifications.TryGetValue(residue, out List<Modification> annotated)
                    && annotated.Any(requirement.IsSatisfiedBy))
                || (requirement.IsForbidden ? fixedModifications : configured).Any(modification => requirement.IsSatisfiedBy(modification)
                    && ModificationLocalization.ModFits(modification, protein.BaseSequence, residue, protein.Length, residue)));

        private bool IsTruncationBoundary(int cut)
        {
            foreach (TruncationProduct product in _protein.TruncationProducts)
            {
                if (product.OneBasedBeginPosition - 1 == cut || product.OneBasedEndPosition == cut)
                {
                    return true;
                }
            }

            return false;
        }

        // The most sites only a required modification lets the protease cut that one read-through can span, so that
        // the unoccupied form across them is enumerated. Without a length limit a read-through is capped at 60 residues.
        private static int ReadThroughSlack(DigestionAgent protease, string sequence, List<int> sites, int maxLength,
            int maxMissedCleavages)
        {
            bool[] requirementOnly = sites
                .Select(site => site > 0 && site < sequence.Length && !Justifies(protease, sequence, site, (_, _) => false))
                .ToArray();
            int lengthLimit = maxLength == int.MaxValue ? 60 : maxLength;
            int most = 0, first = 0, ordinaryInside = 0, requirementOnlyInside = 0;
            for (int last = 1; last < sites.Count; last++)
            {
                if (last - 1 > first)
                {
                    if (requirementOnly[last - 1]) requirementOnlyInside++; else ordinaryInside++;
                }

                while (sites[last] - sites[first] > lengthLimit || ordinaryInside > maxMissedCleavages)
                {
                    first++;
                    if (first < last)
                    {
                        if (requirementOnly[first]) requirementOnlyInside--; else ordinaryInside--;
                    }
                }

                most = Math.Max(most, requirementOnlyInside);
            }

            return most;
        }
    }
}
