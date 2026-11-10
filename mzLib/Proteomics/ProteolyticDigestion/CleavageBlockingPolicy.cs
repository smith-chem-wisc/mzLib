using System.Collections.Generic;
using System.Linq;
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

        private CleavageBlockingPolicy(DigestionAgent protease, int maxMissedCleavages, int generationSlack,
            bool[] sites, bool[] chargeDirectedSites)
        {
            _protease = protease;
            MaxMissedCleavages = maxMissedCleavages;
            GenerationSlack = generationSlack;
            _sites = sites;
            _chargeDirectedSites = chargeDirectedSites;
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
            Protein protein)
        {
            if (digestionParams is null
                || !digestionParams.RespectCleavageBlockingModifications
                || digestionParams.SearchModeType != CleavageSpecificity.Full
                || !(AnyCanBlockCleavage(variableModifications, digestionParams.Protease)
                    || AnyAnnotatedCanBlockCleavage(protein, digestionParams.Protease)))
            {
                return default;
            }

            string sequence = protein?.BaseSequence ?? string.Empty;
            bool[] sites = Sites(digestionParams.Protease, sequence);
            return new CleavageBlockingPolicy(digestionParams.Protease, digestionParams.MaxMissedCleavages, digestionParams.MaxMods,
                sites, ChargeDirectedSites(digestionParams.Protease, sequence, sites));
        }

        private static bool AnyAnnotatedCanBlockCleavage(Protein protein, DigestionAgent protease) =>
            protein?.OneBasedPossibleLocalizedModifications is not null
            && protein.OneBasedPossibleLocalizedModifications.Values.Any(modifications => AnyCanBlockCleavage(modifications, protease));

        /// <summary>
        /// <see cref="DigestionAgent.GetDigestionSiteIndices"/> as a lookup, termini included.
        /// </summary>
        private static bool[] Sites(DigestionAgent protease, string sequence)
        {
            var sites = new bool[sequence.Length + 1];
            foreach (int site in protease.GetDigestionSiteIndices(sequence))
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
            IsActive
            && oneBasedResidue > 0 && oneBasedResidue < _chargeDirectedSites.Length
            && _chargeDirectedSites[oneBasedResidue]
            && CleavageBlockingModifications.BlocksCleavageBy(modification, _protease);
    }
}
