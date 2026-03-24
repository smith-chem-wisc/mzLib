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

            // Identify the 1-based positions of all cleavage sites in the
            // canonical ring. A cleavage site at position p means the protease
            // cuts after residue p, so genuine sub-peptides and ring-opening
            // products always start at position p+1 (mod N, 1-based).
            // Products that start at position 1 are only genuine if position N
            // (the last residue) is itself a cleavage site.
            var cleavagePositions = GetCleavagePositionsInRing(digestionParams);

            // Valid start positions for genuine cut products: one position after
            // each cleavage site, in 1-based canonical ring coordinates.
            // These are the only positions at which a cut can open a new fragment.
            var validStartPositions = new HashSet<int>(
                cleavagePositions.Select(p => p == n ? 1 : p + 1));

            var proxyProtein = new Protein(this, BaseSequence + BaseSequence[..^1]);
            var proxyPeptides = proxyProtein.Digest(digestionParams,
                allKnownFixedModifications, variableModifications,
                silacLabels, turnoverLabels, topDownTruncationSearch);

            var seen = new HashSet<string>();
            bool anyYielded = false;

            foreach (var peptide in proxyPeptides)
            {
                if (peptide is not PeptideWithSetModifications pwsm)
                {
                    yield return peptide;
                    anyYielded = true;
                    continue;
                }

                if (pwsm.Length > n || pwsm.OneBasedStartResidueInProtein > n)
                    continue;

                // Only accept products whose start position corresponds to a
                // genuine cut site in the canonical ring. This excludes proxy
                // boundary artifacts (products starting at the proxy N-terminus
                // when that position is not itself after a cleavage residue).
                int canonicalStart = pwsm.OneBasedStartResidueInProtein <= n
                    ? pwsm.OneBasedStartResidueInProtein
                    : pwsm.OneBasedStartResidueInProtein - n;

                if (!validStartPositions.Contains(canonicalStart))
                    continue;

                if (!seen.Add(pwsm.FullSequence))
                    continue;

                anyYielded = true;

                if (pwsm.Length == n)
                {
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
                else
                {
                    yield return new CircularPeptideWithSetModifications(
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

            if (!anyYielded)
            {
                // Unmodified full-ring peptide — always produced
                yield return new CircularPeptideWithSetModifications(
                    protein: this,
                    digestionParams: digestionParams,
                    oneBasedStartResidueInProtein: 1,
                    oneBasedEndResidueInProtein: n,
                    cleavageSpecificity: CleavageSpecificity.Full,
                    peptideDescription: null,
                    missedCleavages: 0,
                    allModsOneIsNterminus: new Dictionary<int, Modification>(),
                    numFixedMods: 0,
                    baseSequence: BaseSequence);

                // For each variable modification with LocationRestriction "Anywhere.",
                // scan the canonical sequence for matching residues and yield one
                // modified full-ring peptide per match site.
                // N-terminal and C-terminal restrictions are excluded — circular
                // peptides have no termini.
                if (variableModifications != null)
                {
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
                                int modKey = pos + 1;
                                var modDict = new Dictionary<int, Modification> { [modKey] = mod };

                                yield return new CircularPeptideWithSetModifications(
                                    protein: this,
                                    digestionParams: digestionParams,
                                    oneBasedStartResidueInProtein: 1,
                                    oneBasedEndResidueInProtein: n,
                                    cleavageSpecificity: CleavageSpecificity.Full,
                                    peptideDescription: null,
                                    missedCleavages: 0,
                                    allModsOneIsNterminus: modDict,
                                    numFixedMods: 0,
                                    baseSequence: BaseSequence);
                            }
                        }
                    }
                }
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