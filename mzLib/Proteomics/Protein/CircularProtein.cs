using System;
using System.Collections.Generic;
using System.Linq;
using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;

namespace Proteomics
{
    /// <summary>
    /// Represents a head-to-tail cyclic protein or peptide.
    ///
    /// CANONICAL FORM AND NUMBERING
    /// ----------------------------
    /// The sequence is always stored as the lexicographically smallest rotation
    /// (canonical form), regardless of the input rotation. This is the shared
    /// numbering system for both the protein and all peptides derived from it:
    ///
    ///   - Residue 1 is the first character of <see cref="Protein.BaseSequence"/>.
    ///   - Residues are numbered 1..N continuously around the ring.
    ///   - For fragments that wrap past residue N, numbering continues into a
    ///     conceptual second copy of the ring: residue N+1 = residue 1, N+2 = 2,
    ///     etc. A wrap-around fragment starting at residue s of length L is
    ///     annotated [s, s+L-1] where s+L-1 may exceed N.
    ///
    /// Example: ring G-K-A-E-C, canonical form A-E-C-G-K.
    ///   Residue 1=A, 2=E, 3=C, 4=G, 5=K.
    ///   A fragment covering K,A,E (wrapping) starts at residue 5 with length 3:
    ///   annotation [5, 7].
    ///
    /// Because both <see cref="CircularProtein"/> and
    /// <see cref="CircularPeptideWithSetModifications"/> are always canonicalized,
    /// two circular molecules represent the same ring if and only if their
    /// <see cref="Protein.BaseSequence"/> strings are equal.
    ///
    /// DIGESTION
    /// ---------
    /// Digestion is performed on a doubled (2N−1) proxy sequence
    /// (<see cref="BaseSequence"/> + <see cref="BaseSequence"/>[..^1]), then filtered
    /// to peptides of length ≤ N and deduplicated by full sequence. Wrapping peptides
    /// have <c>OneBasedEndResidueInProtein &gt; N</c>, consistent with the annotation
    /// convention above.
    ///
    /// DISPATCH NOTE
    /// -------------
    /// Callers must hold a <see cref="CircularProtein"/> reference (not
    /// <see cref="Protein"/> or <see cref="IBioPolymer"/>) for the shadowed
    /// <see cref="Digest"/> to be dispatched correctly, since
    /// <see cref="Protein.Digest"/> is not virtual.
    /// </summary>
    public class CircularProtein : Protein
    {
        /// <summary>
        /// The canonical circular sequence (lexicographically smallest rotation of the input).
        /// Identical to <see cref="Protein.BaseSequence"/>; exposed under this name for clarity.
        /// </summary>
        public string CircularSequence => BaseSequence;

        /// <summary>
        /// Monoisotopic mass of the cyclic molecule: sum of residue masses with no H₂O addition,
        /// since there are no free N- or C-termini in a head-to-tail cyclized molecule.
        /// </summary>
        public double CyclicMonoisotopicMass { get; }

        // ── Constructors ──────────────────────────────────────────────────────

        /// <summary>
        /// Constructs a <see cref="CircularProtein"/>. The sequence is canonicalized
        /// to the lexicographically smallest rotation before being stored.
        /// </summary>
        public CircularProtein(
            string sequence,
            string accession,
            string organism = null,
            List<Tuple<string, string>> geneNames = null,
            IDictionary<int, List<Modification>> oneBasedModifications = null,
            string name = null,
            string fullName = null,
            bool isDecoy = false,
            bool isContaminant = false,
            List<DatabaseReference> databaseReferences = null,
            string databaseFilePath = null)
            : base(
                sequence: GetCanonicalRotation(sequence),
                accession: accession,
                organism: organism,
                geneNames: geneNames,
                oneBasedModifications: oneBasedModifications,
                name: name,
                fullName: fullName,
                isDecoy: isDecoy,
                isContaminant: isContaminant,
                databaseReferences: databaseReferences,
                databaseFilePath: databaseFilePath)
        {
            CyclicMonoisotopicMass = BaseSequence.Sum(aa => Residue.ResidueMonoisotopicMass[aa]);
        }

        /// <summary>
        /// Creates a <see cref="CircularProtein"/> from an existing <see cref="Protein"/>,
        /// canonicalizing its sequence. Use after loading entries from a FASTA or XML database.
        /// </summary>
        public static CircularProtein FromProtein(Protein source) =>
            new(source.BaseSequence,
                source.Accession,
                source.Organism,
                source.GeneNames,
                source.OneBasedPossibleLocalizedModifications
                      .ToDictionary(kv => kv.Key, kv => kv.Value),
                source.Name,
                source.FullName,
                source.IsDecoy,
                source.IsContaminant,
                source.DatabaseReferences,
                source.DatabaseFilePath);

        // ── Digestion ─────────────────────────────────────────────────────────

        /// <summary>
        /// Digests the circular sequence by running linear digestion on a doubled
        /// proxy sequence of length 2N−1, then filtering and deduplicating results.
        /// Wrapping peptides have <c>OneBasedEndResidueInProtein &gt; N</c>.
        /// </summary>
        public new IEnumerable<IBioPolymerWithSetMods> Digest(
    IDigestionParams digestionParams,
    List<Modification> allKnownFixedModifications,
    List<Modification> variableModifications,
    List<SilacLabel> silacLabels = null,
    (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null,
    bool topDownTruncationSearch = false)
        {
            int n = Length;

            // ── Step 1: Identify cleavage sites in the ring ───────────────────
            //
            // Find all 1-based positions in the canonical ring at which the
            // protease cuts after that residue.
            var cleavagePositions = GetCleavagePositionsInRing(digestionParams);
            int numCleavageSites = cleavagePositions.Count;
            int maxMissedCleavages = ((DigestionParams)digestionParams).MaxMissedCleavages;

            // ── Step 2: Emit the CircularPeptideWithSetModifications (1:1) ────
            //
            // The CircularPeptideWithSetModifications is ALWAYS a direct 1:1
            // conversion of the CircularProtein — the whole ring, intact.
            // It is NEVER produced by cutting: not even a single cut.
            // It is emitted here, independently of the proxy digestion below.
            //
            // It is only possible when the missed-cleavage budget is large enough
            // to absorb every cleavage site in the ring, meaning the entire ring
            // is one un-cut cyclic segment with no free N- or C-termini.
            //
            // A CircularPeptideWithSetModifications is ALWAYS length N.
            // It is NEVER shorter than N and NEVER starts at any position other
            // than position 1 (the canonical origin).
            if (maxMissedCleavages >= numCleavageSites)
            {
                // Unmodified full-ring product — always emitted when the budget allows.
                yield return new CircularPeptideWithSetModifications(
                    protein: this,
                    digestionParams: digestionParams,
                    oneBasedStartResidueInProtein: 1,
                    oneBasedEndResidueInProtein: n,
                    cleavageSpecificity: CleavageSpecificity.Full,
                    peptideDescription: null,
                    missedCleavages: numCleavageSites,
                    allModsOneIsNterminus: new Dictionary<int, Modification>(),
                    numFixedMods: 0,
                    baseSequence: BaseSequence);

                // Variable-modified full-ring products.
                // Circular peptides have no N- or C-termini, so only
                // "Anywhere." modifications are applicable.
                //
                // We enumerate all combinations of up to MaxMods modifications
                // across all matching (key, mod) pairs, mirroring the combinatorial
                // enumeration that linear peptide digestion performs.
                // MaxMods == 0 means no modified forms are produced.
                if (variableModifications != null && maxMissedCleavages >= numCleavageSites)
                {
                    // Build the list of all (key, mod) pairs that fit this ring.
                    // Key convention: side-chain at 1-based position pos → key = pos + 1.
                    // Circular peptides have no termini, so N-term (key=1) and
                    // C-term (key=n+2) restrictions are excluded.
                    var applicableSites = new List<(int Key, Modification Mod)>();
                    foreach (var mod in variableModifications)
                    {
                        if (mod.LocationRestriction != "Anywhere.") continue;

                        for (int pos = 1; pos <= n; pos++)
                        {
                            if (ModificationLocalization.ModFits(
                                    mod,
                                    BaseSequence,
                                    digestionProductOneBasedIndex: pos,
                                    digestionProductLength: n,
                                    bioPolymerOneBasedIndex: pos))
                            {
                                applicableSites.Add((pos + 1, mod));
                            }
                        }
                    }

                    // Enumerate all subsets of applicableSites of size 1..MaxMods.
                    // Each subset where all keys are distinct produces one modified peptide.
                    // (Two mods at the same site are excluded — one mod per residue.)
                    int maxMods = ((DigestionParams)digestionParams).MaxMods;
                    foreach (var subset in GetModSubsets(applicableSites, maxMods))
                    {
                        yield return new CircularPeptideWithSetModifications(
                            protein: this,
                            digestionParams: digestionParams,
                            oneBasedStartResidueInProtein: 1,
                            oneBasedEndResidueInProtein: n,
                            cleavageSpecificity: CleavageSpecificity.Full,
                            peptideDescription: null,
                            missedCleavages: numCleavageSites,
                            allModsOneIsNterminus: subset,
                            numFixedMods: 0,
                            baseSequence: BaseSequence);
                    }
                }
            }

            // ── Step 3: Emit linear PeptideWithSetModifications products ──────
            //
            // All linear products come from the proxy digestion of a doubled
            // sequence (BaseSequence + BaseSequence[..^1], length 2N-1), which
            // allows the linear digestion engine to discover wrap-around fragments.
            //
            // CRITICAL: Every product from the proxy digestion is ALWAYS a
            // PeptideWithSetModifications — linear, with free termini.
            // The CircularPeptideWithSetModifications was emitted above and must
            // NEVER be produced from the proxy.
            //
            // There are two sub-cases, both yielding PeptideWithSetModifications:
            //
            //   Single cut (Length == N):
            //     The ring was opened at exactly one cleavage site. The product
            //     spans the full canonical sequence and has the same length N as
            //     the circular product, but it is a completely different object:
            //     it has free N- and C-termini and carries the standard +H2O mass.
            //     Do NOT conflate length equality with type equality.
            //
            //   Two or more cuts (Length < N):
            //     The ring was opened and sub-divided. Always linear, always < N.
            if (numCleavageSites == 0)
                yield break; // No cuts possible; only the circular product exists.

            // Valid start positions: one position after each cleavage site in
            // 1-based canonical ring coordinates.
            var validStartPositions = new HashSet<int>(
                cleavagePositions.Select(p => p == n ? 1 : p + 1));

            var proxyProtein = new Protein(this, BaseSequence + BaseSequence[..^1]);
            var proxyPeptides = proxyProtein.Digest(digestionParams,
                allKnownFixedModifications, variableModifications,
                silacLabels, turnoverLabels, topDownTruncationSearch);

            var seen = new HashSet<string>();

            foreach (var peptide in proxyPeptides)
            {
                if (peptide is not PeptideWithSetModifications pwsm)
                    continue; // Non-linear proxy products are not expected; skip.

                // Discard fragments longer than N, or starting in the second copy
                // of the proxy (those duplicate first-copy fragments).
                if (pwsm.Length > n || pwsm.OneBasedStartResidueInProtein > n)
                    continue;

                // Discard proxy artifacts not starting at a genuine cut boundary.
                if (!validStartPositions.Contains(pwsm.OneBasedStartResidueInProtein))
                    continue;

                // Deduplicate by full sequence.
                if (!seen.Add(pwsm.FullSequence))
                    continue;

                // Every proxy product is linear — length N (single-cut) or < N (multi-cut).
                // Neither case is ever a CircularPeptideWithSetModifications.
                yield return new PeptideWithSetModifications(
                    protein: this,
                    digestionParams: pwsm.DigestionParams,
                    oneBasedStartResidueInProtein: pwsm.OneBasedStartResidueInProtein,
                    oneBasedEndResidueInProtein: pwsm.OneBasedEndResidueInProtein,
                    cleavageSpecificity: pwsm.CleavageSpecificityForFdrCategory,
                    peptideDescription: pwsm.PeptideDescription,
                    missedCleavages: pwsm.MissedCleavages,
                    allModsOneIsNterminus: pwsm.AllModsOneIsNterminus,
                    numFixedMods: pwsm.NumFixedMods,
                    baseSequence: pwsm.BaseSequence);
            }
        }

        /// <summary>
        /// Yields all non-empty subsets of <paramref name="sites"/> of size at most
        /// <paramref name="maxMods"/>, where every key in the subset is distinct
        /// (one modification per residue). Each subset is returned as a
        /// <see cref="Dictionary{TKey,TValue}"/> ready for use as
        /// <c>AllModsOneIsNterminus</c>.
        /// </summary>
        private static IEnumerable<Dictionary<int, Modification>> GetModSubsets(
            List<(int Key, Modification Mod)> sites,
            int maxMods)
        {
            if (maxMods <= 0 || sites.Count == 0)
                yield break;

            // Enumerate subsets via recursive backtracking.
            var current = new Dictionary<int, Modification>();

            foreach (var subset in Backtrack(sites, 0, current, maxMods))
                yield return subset;
        }

        private static IEnumerable<Dictionary<int, Modification>> Backtrack(
            List<(int Key, Modification Mod)> sites,
            int startIndex,
            Dictionary<int, Modification> current,
            int maxMods)
        {
            // Yield a copy of the current non-empty subset.
            if (current.Count > 0)
                yield return new Dictionary<int, Modification>(current);

            if (current.Count == maxMods)
                yield break;

            for (int i = startIndex; i < sites.Count; i++)
            {
                var (key, mod) = sites[i];

                // Skip if this residue already has a mod in the current subset.
                if (current.ContainsKey(key))
                    continue;

                current[key] = mod;

                foreach (var subset in Backtrack(sites, i + 1, current, maxMods))
                    yield return subset;

                current.Remove(key);
            }
        }

        private HashSet<int> GetCleavagePositionsInRing(IDigestionParams digestionParams)
        {
            // Proxy sequence: BaseSequence + BaseSequence[..^1], length 2N-1.
            // The last character is omitted deliberately — the proxy represents one
            // full traversal of the ring starting at the canonical origin, continuing
            // through a second partial copy. This means every cleavage residue in the
            // first copy [1..N] is followed by more sequence (from the second copy),
            // so the linear digestion engine will always produce a cut after it and
            // reveal its position as a peptide boundary.
            //
            // We collect the end position of every peptide whose end falls within
            // [1, N] — these are the genuine cleavage sites of the canonical ring.
            var proxyProtein = new Protein(
                BaseSequence + BaseSequence[..^1],
                Accession + "_cleavage");

            var fragments = proxyProtein.Digest(
                new DigestionParams(
                    protease: ((DigestionParams)digestionParams).Protease.Name,
                    maxMissedCleavages: 0,
                    minPeptideLength: 1),
                [], [])
                .Cast<PeptideWithSetModifications>()
                .Where(p => p.OneBasedEndResidueInProtein <= Length)
                .ToList();

            var positions = new HashSet<int>();
            foreach (var f in fragments)
                positions.Add(f.OneBasedEndResidueInProtein);

            return positions;
        }

        // ── Shared numbering system ───────────────────────────────────────────


        /// <summary>
        /// Returns the lexicographically smallest rotation of <paramref name="sequence"/>.
        /// This is the single shared implementation used by both
        /// <see cref="CircularProtein"/> and <see cref="CircularPeptideWithSetModifications"/>.
        ///
        /// The algorithm reads all N rotations, comparing character by character, and
        /// selects the one that sorts earliest under ordinal (byte-value) comparison —
        /// i.e., 'A' &lt; 'C' &lt; ... &lt; 'Y', with subsequent characters as tiebreakers.
        /// </summary>
        public static string GetCanonicalRotation(string sequence)
        {
            if (string.IsNullOrEmpty(sequence))
                return sequence;

            int n = sequence.Length;
            string best = sequence;

            for (int i = 1; i < n; i++)
            {
                // Compare rotation starting at i with the current best, character by character.
                for (int k = 0; k < n; k++)
                {
                    char ci = sequence[(i + k) % n];
                    char cb = best[k];               // best is already a rotation, read linearly

                    if (ci < cb) { best = (sequence + sequence).Substring(i, n); break; }
                    if (ci > cb) break;
                    // equal: continue to next character (tiebreaker)
                }
            }

            return best;
        }
    }
}