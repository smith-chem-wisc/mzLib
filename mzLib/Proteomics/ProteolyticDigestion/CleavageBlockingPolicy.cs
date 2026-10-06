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
    /// halves are decided here, once, from the same <see cref="DigestionParams"/> and the same
    /// modification list.
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
        /// Indexed by one-based residue: true when the protease really cuts after that residue in this
        /// protein AND only because of that residue's side chain, so a blocking modification on it
        /// abolishes the cut. See <see cref="ChargeDirectedSites"/>.
        /// </summary>
        private readonly bool[] _chargeDirectedSites;

        private CleavageBlockingPolicy(DigestionAgent protease, int maxMissedCleavages, int generationSlack,
            bool[] chargeDirectedSites)
        {
            _protease = protease;
            MaxMissedCleavages = maxMissedCleavages;
            GenerationSlack = generationSlack;
            _chargeDirectedSites = chargeDirectedSites;
        }

        /// <summary>
        /// The policy in force for this digestion, or the inert policy when the correction cannot fire.
        /// Called at both ends of the exchange with the same two arguments, which is what guarantees the
        /// slack and the drop agree about the protease and about the missed-cleavage budget.
        /// </summary>
        /// <param name="variableModifications">
        /// The VARIABLE modifications configured for the search. Fixed modifications are deliberately not
        /// consulted -- see the remarks.
        /// </param>
        /// <param name="sequence">
        /// The base sequence of the protein being digested. The policy is per protein: it records which
        /// residues are real protease sites there, so a blocked residue that was never a site (a K before
        /// P under trypsin|P) cannot cancel a genuine missed cleavage.
        /// </param>
        /// <remarks>
        /// Only variable modifications activate the policy, and only a variable modification can block a
        /// cleavage under it. A fixed blocking modification is unavoidable: it sits on every instance of
        /// its target residue, so the read-through that would replace a dropped peptidoform has to span
        /// every such residue in the protein, and no slack computed from <see cref="DigestionParams.MaxMods"/>
        /// -- which bounds variable modifications only -- can pay for it. At MaxMods 0 the slack is zero
        /// while the drop would still fire, and the peptide is lost rather than corrected. Fixed blocking
        /// modifications therefore stay outside this correction entirely; modelling them means removing
        /// their sites from the protease's site list before enumeration, which is a different mechanism.
        /// </remarks>
        public static CleavageBlockingPolicy For(DigestionParams digestionParams, IEnumerable<Modification> variableModifications,
            string sequence)
        {
            if (digestionParams is null
                || !digestionParams.RespectCleavageBlockingModifications
                || digestionParams.SearchModeType != CleavageSpecificity.Full
                || !AnyCanBlockCleavage(variableModifications, digestionParams.Protease))
            {
                return default;
            }

            return new CleavageBlockingPolicy(digestionParams.Protease, digestionParams.MaxMissedCleavages, digestionParams.MaxMods,
                ChargeDirectedSites(digestionParams.Protease, sequence ?? string.Empty));
        }

        /// <summary>
        /// The residues of <paramref name="sequence"/> after which <paramref name="protease"/> really cuts
        /// and could be stopped from cutting by a modification of that residue's side chain.
        /// </summary>
        /// <remarks>
        /// Starts from <see cref="DigestionAgent.GetDigestionSiteIndices"/>, the same site list digestion
        /// enumerates from and counts missed cleavages against, so the discount and the count agree site
        /// for site. The protein's own termini are not cleavages and are excluded. A site is then removed if
        /// any motif that fits there does not name the residue literally at P1 -- StcE's "TX|T" cuts
        /// "TK|T" whatever the lysine carries, so in a StcE-trypsin co-digest that site survives an
        /// acetylated K even though trypsin's "K|" alone would not.
        /// </remarks>
        private static bool[] ChargeDirectedSites(DigestionAgent protease, string sequence)
        {
            var sites = new bool[sequence.Length + 1];
            foreach (int site in protease.GetDigestionSiteIndices(sequence))
            {
                if (site > 0 && site < sequence.Length)
                {
                    sites[site] = true;
                }
            }

            for (int location = 0; location < sequence.Length; location++)
            {
                foreach (DigestionMotif motif in protease.DigestionMotifs)
                {
                    int site = location + motif.CutIndex;
                    if (site > 0 && site < sequence.Length && sites[site]
                        && !motif.CleavesCTerminalTo(sequence[site - 1])
                        && motif.Fits(sequence, location).Item1)
                    {
                        sites[site] = false;
                    }
                }
            }

            return sites;
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
        /// <see cref="DigestionParams.MaxMods"/>: a peptidoform carries at most that many variable
        /// modifications and therefore at most that many blocked sites, so every read-through this policy
        /// can need is reachable, however many blocked sites co-occur.
        /// </remarks>
        public int GenerationSlack { get; }

        /// <summary>
        /// The budget the CALLER asked for -- the bound a surviving peptidoform's reported missed-cleavage
        /// count must satisfy, so the slack above stays an enumeration detail and is never observable in
        /// the result.
        /// </summary>
        public int MaxMissedCleavages { get; }

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
