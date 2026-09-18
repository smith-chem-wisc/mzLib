using MzLibUtil;
using Omics.Modifications;

namespace Omics.Digestion
{
    public abstract class DigestionAgent
    {
        protected static readonly HashSetPool<int> HashSetPool = new HashSetPool<int>(8);

        protected DigestionAgent(string name, CleavageSpecificity cleavageSpecificity, List<DigestionMotif> motifList, Modification cleavageMod)
        {
            Name = name;
            CleavageSpecificity = cleavageSpecificity;
            DigestionMotifs = motifList ?? new List<DigestionMotif>();
            CleavageMod = cleavageMod;
        }

        public readonly string Name;
        public CleavageSpecificity CleavageSpecificity { get; init; }
        public List<DigestionMotif> DigestionMotifs { get; init; }
        public Modification CleavageMod { get; set; }

        public override string ToString()
        {
            return Name;
        }

        public override bool Equals(object? obj)
        {
            return obj is DigestionAgent agent && agent.Name == Name;
        }

        public override int GetHashCode()
        {
            return Name.GetHashCode();
        }

        /// <summary>
        /// True when this agent cleaves the bond C-TERMINAL to <paramref name="residue"/>, under any of
        /// its motifs. Trypsin reports true for K and R; Glu-C reports true for E alone; Asp-N and Lys-N
        /// report false for every residue, because they cut N-terminal to their recognition residue.
        ///
        /// This is what makes a cleavage-blocking modification a question about the PAIR rather than
        /// about the modification alone: an acetylated lysine abolishes a trypsin site and abolishes
        /// nothing at all in a Glu-C digest.
        /// </summary>
        /// <remarks>
        /// Residue-level, context-free -- see <see cref="DigestionMotif.CleavesCTerminalTo"/> for what
        /// that does and does not take into account.
        /// </remarks>
        /// <summary>
        /// True when at least one of this agent's motifs will not cleave unless a modification is present
        /// at one of its subsites -- that is, when this is a glycoprotease rather than an ordinary
        /// sequence-directed protease. False for every agent that ships today.
        /// </summary>
        /// <remarks>
        /// The inertness gate for the whole cleavage-promoting correction, and the counterpart of
        /// ProteinDigestion's AnyConfiguredModificationCanBlockCleavage. It matters because the discharge
        /// runs once per generated peptidoform, in the same loop as
        /// <see cref="Modifications.ModificationLocalization.ModFits"/> -- roughly 8.8 billion calls in a
        /// bottom-up run -- so a digest with no glycoprotease must pay nothing at all for a feature it
        /// cannot use.
        /// </remarks>
        public bool HasCleavageRequirement
        {
            get
            {
                if (DigestionMotifs is null)
                {
                    return false;
                }

                foreach (DigestionMotif motif in DigestionMotifs)
                {
                    if (motif?.CleavageRequirement is not null)
                    {
                        return true;
                    }
                }

                return false;
            }
        }

        /// <summary>
        /// Removes from <paramref name="oneBasedIndicesToCleaveAfter"/> every internal site that no motif
        /// can justify, because the modification a motif REQUIRES at one of its subsites cannot be
        /// present there at all. Returns the list unchanged when this agent requires nothing.
        /// </summary>
        /// <remarks>
        /// <para><b>Why the filter belongs here and not after digestion.</b> A site that is not a site
        /// must never enter the enumeration in the first place. Dropping the peptidoforms afterwards
        /// removes the two fragments either side of a bad cut but does NOT produce the read-through
        /// peptide that replaces them -- that peptide carries one more missed cleavage and, at
        /// MaxMissedCleavages = 0, was never generated. Filtering here instead means peptide spans,
        /// missed-cleavage counts and read-throughs all come out right with no generation slack at all.</para>
        ///
        /// <para><b>Feasibility, not occupancy.</b> This asks whether the parent COULD carry a satisfying
        /// modification at the constrained residue, using the localized modifications the database
        /// declares. It cannot ask whether a particular peptidoform DOES carry one, because peptidoforms
        /// do not exist yet. A site kept here may still be refused per-peptidoform later; a site dropped
        /// here could never have been real. Being wrong in that direction only ever keeps peptides.</para>
        ///
        /// <para>The first and last entries are the sequence's own termini rather than cleavage sites, so
        /// they are always kept: removing them would discard the peptide that runs to the end of the
        /// protein.</para>
        /// </remarks>
        public List<int> FilterToFeasibleCleavageSites(List<int> oneBasedIndicesToCleaveAfter, IBioPolymer parent,
            IEnumerable<Modification> configuredModifications = null)
        {
            if (!HasCleavageRequirement || oneBasedIndicesToCleaveAfter is null || parent is null)
            {
                return oneBasedIndicesToCleaveAfter;
            }

            string sequence = parent.BaseSequence;
            var feasible = new List<int>(oneBasedIndicesToCleaveAfter.Count);

            for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count; i++)
            {
                int site = oneBasedIndicesToCleaveAfter[i];

                // The injected termini are not cleavage events and are never filtered.
                if (i == 0 || i == oneBasedIndicesToCleaveAfter.Count - 1 || site <= 0 || site >= sequence.Length)
                {
                    feasible.Add(site);
                    continue;
                }

                if (AnyMotifCouldJustify(site, sequence, parent, configuredModifications))
                {
                    feasible.Add(site);
                }
            }

            return feasible;
        }

        /// <summary>
        /// True when some motif both matches the sequence at this cut and could have its requirement met
        /// there. A motif carrying no requirement justifies any cut it matches, which is what lets a
        /// composite agent keep cutting at its ordinary sequence motifs.
        /// </summary>
        private bool AnyMotifCouldJustify(int cutAfterOneBasedResidue, string sequence, IBioPolymer parent,
            IEnumerable<Modification> configuredModifications)
        {
            foreach (DigestionMotif motif in DigestionMotifs)
            {
                if (motif is null)
                {
                    continue;
                }

                // Fits takes a ZERO-based index into the sequence handed to it, and the motif's
                // recognition sequence begins CutIndex residues before the bond it severs.
                int motifStartZeroBased = cutAfterOneBasedResidue - motif.CutIndex;
                if (motifStartZeroBased < 0 || motifStartZeroBased + motif.InducingCleavage.Length > sequence.Length)
                {
                    continue;
                }

                (bool fits, bool prevented) = motif.Fits(sequence, motifStartZeroBased);
                if (!fits || prevented)
                {
                    continue;
                }

                if (motif.CleavageRequirement is null)
                {
                    return true;
                }

                // Subsites count outward from the bond, which falls after cutAfterOneBasedResidue:
                // Pk is (cut - k + 1) and Pk' is (cut + k), both one-based in the parent.
                CleavageRequirement requirement = motif.CleavageRequirement;
                int constrainedResidue = requirement.IsPrimeSide
                    ? cutAfterOneBasedResidue + requirement.Subsite
                    : cutAfterOneBasedResidue - requirement.Subsite + 1;

                if (constrainedResidue < 1 || constrainedResidue > sequence.Length)
                {
                    continue;
                }

                if (CouldCarrySatisfyingModification(requirement, constrainedResidue, sequence, parent, configuredModifications))
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// True when this agent could cut after <paramref name="cutAfterOneBasedResidue"/> in
        /// <paramref name="parent"/> -- that is, when the site would survive
        /// <see cref="FilterToFeasibleCleavageSites"/> and therefore appears in the site list digestion
        /// actually enumerated.
        /// </summary>
        /// <remarks>
        /// Needed by the missed-cleavage discount, which has to tell three kinds of internal position
        /// apart: a position that is no site at all (never in the list, so never a missed cleavage), a
        /// site that is real for this peptidoform (a genuine missed cleavage), and a site that is feasible
        /// but unoccupied here (in the list, but not a cleavage the agent could have made, so it must be
        /// discounted). Only the middle one may count against the caller's budget.
        /// </remarks>
        public bool IsFeasibleCleavageSite(int cutAfterOneBasedResidue, IBioPolymer parent,
            IEnumerable<Modification> configuredModifications = null)
        {
            if (parent is null || cutAfterOneBasedResidue < 1 || cutAfterOneBasedResidue >= parent.BaseSequence.Length)
            {
                return false;
            }

            return AnyMotifCouldJustify(cutAfterOneBasedResidue, parent.BaseSequence, parent, configuredModifications);
        }

        /// <summary>
        /// The most feasible cleavage sites that can fall strictly inside a single peptide of at most
        /// <paramref name="maxPeptideLength"/> residues, and therefore the generation slack the
        /// cleavage-promoting correction needs.
        /// </summary>
        /// <remarks>
        /// <para><b>Why the slack is not MaxMods, unlike the blocking mirror.</b> A blocked site is a site
        /// carrying a modification, and a peptidoform carries at most MaxMods of them, so MaxMods bounds
        /// the blocking slack exactly. An unjustified site is the COMPLEMENT -- a feasible site with NO
        /// modification on it -- and nothing about the modification budget limits how many of those a
        /// peptide may span. A peptidoform with no glycan at all has every feasible site it spans
        /// unjustified. So the bound has to come from the sequence and the length limit instead.</para>
        ///
        /// <para>This is an exact upper bound: any peptide the caller could legally receive is at most
        /// maxPeptideLength residues long, so it cannot contain more internal sites than the densest
        /// window of that length. Buying exactly that much slack means every read-through is reachable
        /// and none is merely hoped for.</para>
        ///
        /// <para>The slack is cheap even when it is large. It widens the outer loop of the span
        /// enumeration, but every extra span is rejected by the length check before a peptide object is
        /// built, so the cost is loop iterations rather than digestion work. A protein with few sites
        /// buys little slack; a mucin buys a lot and needs it.</para>
        /// </remarks>
        public static int MaximumInternalSitesInOnePeptide(List<int> oneBasedIndicesToCleaveAfter, int maxPeptideLength)
        {
            if (oneBasedIndicesToCleaveAfter is null || oneBasedIndicesToCleaveAfter.Count < 3)
            {
                return 0;
            }

            // Widest pair of sites still within the length limit; the internal sites are those between.
            int most = 0;
            int start = 0;
            for (int end = 1; end < oneBasedIndicesToCleaveAfter.Count; end++)
            {
                while (start < end
                       && oneBasedIndicesToCleaveAfter[end] - oneBasedIndicesToCleaveAfter[start] > maxPeptideLength)
                {
                    start++;
                }

                int internalSites = end - start - 1;
                if (internalSites > most)
                {
                    most = internalSites;
                }
            }

            return most;
        }

        /// <summary>
        /// True when a modification satisfying <paramref name="requirement"/> could occupy
        /// <paramref name="constrainedResidue"/> -- either because the database annotates one there, or
        /// because the search configures one that fits there.
        /// </summary>
        /// <remarks>
        /// <para><b>Both sources have to be consulted, and consulting only the first is a real defect.</b>
        /// Database annotations answer "this protein is known to be glycosylated here"; configured variable
        /// and fixed modifications answer "this search is willing to place a glycan wherever it fits". A
        /// search supplying an O-glycan as a variable modification against an unannotated database -- the
        /// ordinary MetaMorpheus configuration -- has every motif site feasible, and judging it on
        /// annotations alone found none of them, filtered away every site, and returned the undigested
        /// protein as the only product. The protease was silently switched off.</para>
        ///
        /// <para>Feasibility from a configured modification is deliberately weak: a variable modification
        /// that fits anywhere makes every site feasible, so the filter stops discriminating and the
        /// per-peptidoform occupancy check downstream does the work instead. That is the correct division --
        /// "could be anywhere" genuinely is no constraint on the site list -- and it is why this filter
        /// bites hardest exactly where the evidence is strongest, on a database with localized glycosites.</para>
        /// </remarks>
        private static bool CouldCarrySatisfyingModification(CleavageRequirement requirement, int constrainedResidue,
            string sequence, IBioPolymer parent, IEnumerable<Modification> configuredModifications)
        {
            if (parent.OneBasedPossibleLocalizedModifications is not null
                && parent.OneBasedPossibleLocalizedModifications.TryGetValue(constrainedResidue, out var candidates)
                && candidates is not null
                && candidates.Any(requirement.IsSatisfiedBy))
            {
                return true;
            }

            if (configuredModifications is null)
            {
                return false;
            }

            foreach (Modification configured in configuredModifications)
            {
                // The whole parent is passed as the "product" because this asks whether the modification
                // could ever sit at this residue of this sequence, not whether it sits on some peptide.
                if (requirement.IsSatisfiedBy(configured)
                    && ModificationLocalization.ModFits(configured, sequence, constrainedResidue, sequence.Length, constrainedResidue))
                {
                    return true;
                }
            }

            return false;
        }

        public bool CleavesCTerminalTo(char residue)
        {
            foreach (DigestionMotif motif in DigestionMotifs)
            {
                if (motif.CleavesCTerminalTo(residue))
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Is length of given peptide okay, given minimum and maximum?
        /// </summary>
        /// <param name="length"></param>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <returns></returns>
        protected static bool ValidLength(int length, int minLength, int maxLength)
        {
            return ValidMinLength(length, minLength) && ValidMaxLength(length, maxLength);
        }

        /// <summary>
        /// Is length of given peptide okay, given minimum?
        /// </summary>
        /// <param name="length"></param>
        /// <param name="minLength"></param>
        /// <returns></returns>
        protected static bool ValidMinLength(int length, int minLength)
        {
            return length >= minLength;
        }

        /// <summary>
        /// Is length of given peptide okay, given maximum?
        /// </summary>
        /// <param name="length"></param>
        /// <param name="maxLength"></param>
        /// <returns></returns>
        protected static bool ValidMaxLength(int? length, int maxLength)
        {
            return !length.HasValue || length <= maxLength;
        }

        /// <summary>
        /// Gets the indices after which this protease will cleave a given protein sequence
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        public List<int> GetDigestionSiteIndices(string sequence)
        {
            var indices = HashSetPool.Get(); // use hash set to ensure no duplicates
            try // Try block is to ensure that, even if an error gets thrown, the hashset is returned to the pool
            {
                indices.Add(0); // The start of the protein is treated as a cleavage site to retain the n-terminal peptide

                for (int r = 0; r < sequence.Length; r++)
                {
                    var cutSiteIndex = -1;
                    bool cleavagePrevented = false;

                    foreach (DigestionMotif motif in DigestionMotifs)
                    {
                        var motifResults = motif.Fits(sequence, r);
                        bool motifFits = motifResults.Item1;
                        bool motifPreventsCleavage = motifResults.Item2;

                        if (motifFits && r + motif.CutIndex < sequence.Length)
                        {
                            cutSiteIndex = Math.Max(r + motif.CutIndex, cutSiteIndex);
                        }

                        if (motifPreventsCleavage) // if any motif prevents cleave
                        {
                            cleavagePrevented = true;
                        }
                    }

                    // if no motif prevents cleave
                    if (!cleavagePrevented && cutSiteIndex != -1)
                    {
                        indices.Add(cutSiteIndex);
                    }
                }

                indices.Add(sequence.Length); // The end of the protein is treated as a cleavage site to retain the c-terminal peptide
                return indices.ToList(); // convert the hashset to a list for return. 
            }
            finally
            {
                // return hashset to pool. This clears it and gets it ready for the next time it is needed from the pool.
                HashSetPool.Return(indices);
            }
        }
    }
}
