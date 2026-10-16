#nullable enable
using MzLibUtil;
using Omics.Modifications;
using Omics.Digestion;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace UsefulProteomicsDatabases.EntrapmentGeneration;

/// <summary>
/// Finds the target peptide an entrapment peptide was built from, using nothing but the target
/// protein and the peptide itself.
/// </summary>
/// <remarks>
/// <para>A rearrangement preserves the residue multiset and leaves every cleavage-site residue in
/// place, so those two things together are a key: the partner of an entrapment peptide is the
/// target peptide of the same protein with the same free-residue composition and the same pinned
/// pattern. That is guaranteed by the construction rather than recorded anywhere, which is why no
/// index file has to be written, shipped, or kept in step with the database. A search reports a
/// protein accession and a peptide sequence, and that is enough.</para>
/// <para>The key can collide: a protein may contain two different peptides sharing both, and
/// P57071 really does contain <c>LIHTGVK</c> and <c>LIHTVGK</c>. Those are reported in
/// <see cref="AmbiguousPeptides"/> and refuse to resolve, because a guess between them would be
/// wrong half the time and silently so.</para>
/// </remarks>
public sealed class EntrapmentPairing
{
    private readonly Dictionary<string, string> _byKey;
    private readonly HashSet<string> _ambiguous;
    private readonly List<DigestionMotif> _motifs;
    private readonly HashSet<string> _truncationProductPeptides = new();

    /// <summary>Indexes a target protein's peptides, missed cleavages included, by composition and pinned pattern.</summary>
    public EntrapmentPairing(Protein target, IDigestionParams digestionParams)
    {
        if (target is null)
        {
            throw new MzLibException("Cannot build a pairing without a target protein.");
        }
        EntrapmentAssembler.RefuseUnusableDigestionParams(digestionParams);

        _byKey = new Dictionary<string, string>();
        _ambiguous = new HashSet<string>();
        var collided = new HashSet<string>();
        var indexed = new HashSet<string>();

        _motifs = digestionParams.DigestionAgent.DigestionMotifs;
        List<int> sites = digestionParams.DigestionAgent.GetDigestionSiteIndices(target.BaseSequence);
        sites.Sort();

        // Index runs of adjacent base pieces, not just single ones: a search reports missed-cleavage
        // peptides too, and those are pairable for the same reason -- base-piece compositions add, so
        // a run is isomeric with the run its partner was built from.
        int maxPieces = digestionParams.MaxMissedCleavages + 1;

        // Only peptides a search could actually report. Without this the index holds every run
        // however short, and short peptides collide constantly -- "AK" shares its key with every
        // other "AK" -- which inflated the ambiguity count nearly tenfold with peptides nobody will
        // ever see (6,889 against 736 on the reviewed human proteome). Filtering is safe as well as
        // honest: a key begins with the peptide's length, so peptides of different lengths never
        // collide and removing the short ones cannot change how a searchable peptide resolves.
        int minLength = digestionParams.MinLength;
        int maxLength = digestionParams.MaxLength;

        // The initiator methionine, by the same rule DigestionAgent.FullDigestion applies: a run
        // beginning at the protein's first residue is also emitted from residue 2 unless the
        // behaviour is Retain, and is emitted only from residue 2 when it is Cleave. Nothing extra
        // when the opening piece is a lone M, since the run from residue 2 is then an ordinary one.
        string sequence = target.BaseSequence;
        InitiatorMethionineBehavior initiatorMethionine =
            (digestionParams as DigestionParams)?.InitiatorMethionineBehavior ?? InitiatorMethionineBehavior.Variable;
        bool startsWithMethionine = sequence.Length > 0 && sequence[0] == 'M';
        bool retainMethionine = initiatorMethionine != InitiatorMethionineBehavior.Cleave || !startsWithMethionine;
        bool cleaveMethionine = initiatorMethionine != InitiatorMethionineBehavior.Retain && startsWithMethionine
                                && sites.Count > 1 && sites[1] != 1;

        for (int first = 0; first < sites.Count - 1; first++)
        {
            for (int pieces = 1; pieces <= maxPieces && first + pieces < sites.Count; pieces++)
            {
                int start = sites[first];
                int end = sites[first + pieces];
                if (first != 0 || retainMethionine)
                {
                    Index(start, end);
                }
                if (first == 0 && cleaveMethionine)
                {
                    Index(1, end);
                }
            }
        }

        void Index(int start, int end)
        {
            int length = end - start;
            if (length < minLength || length > maxLength)
            {
                return;
            }

            string peptide = sequence.Substring(start, length);
            indexed.Add(peptide);
            string key = KeyOf(peptide);

            if (_byKey.TryGetValue(key, out string? existing))
            {
                if (existing != peptide)
                {
                    collided.Add(key);
                    _ambiguous.Add(existing);
                    _ambiguous.Add(peptide);
                }
                return;
            }

            _byKey[key] = peptide;
        }

        foreach (string key in collided)
        {
            _byKey.Remove(key);
        }

        // Peptides a search reports at signal-peptide, propeptide and chain boundaries. Digestion
        // emits them because the target carries those truncation products; a partner carries none
        // (positional annotations do not survive the rearrangement), so they have no partner and
        // never will. Counted and named rather than indexed, so the r computed from this class is
        // honest about the population it covers. Protein.Digest is the oracle, so this is exactly
        // what a search adds and nothing else.
        if (target.TruncationProducts.Any())
        {
            var noMods = new List<Modification>();
            foreach (var peptide in target.Digest(digestionParams, noMods, noMods))
            {
                if (!indexed.Contains(peptide.BaseSequence))
                {
                    _truncationProductPeptides.Add(peptide.BaseSequence);
                }
            }
        }
    }

    /// <summary>
    /// Target peptides a search reports only because of the target's truncation products -- they
    /// begin or end at a signal-peptide, propeptide or chain boundary rather than at a cleavage
    /// site. A partner carries no truncation products, so these have no partner, and they are left
    /// out of <see cref="SearchablePeptideCount"/>. A paired estimator should exclude them as it
    /// excludes <see cref="AmbiguousPeptides"/>.
    /// </summary>
    public IReadOnlyCollection<string> TruncationProductPeptides => _truncationProductPeptides;

    /// <summary>
    /// Target peptides sharing a composition-and-pinning key with another peptide of the same
    /// protein, and therefore not resolvable. Report these rather than pairing them.
    /// </summary>
    public IReadOnlyCollection<string> AmbiguousPeptides => _ambiguous;

    /// <summary>
    /// Distinct peptides of this protein a search could report -- runs of up to
    /// <c>MaxMissedCleavages + 1</c> base pieces, within the length bounds, plus the forms of the
    /// N-terminal runs a search reports without their initiator methionine. This is the population
    /// an FDP estimator's <c>r</c> is over, and the denominator an ambiguity rate needs: the
    /// report's own peptide counts are over <i>base pieces</i>, which is a different and smaller
    /// population, and dividing one by the other gives a rate of nothing.
    /// <para>Peptides at truncation-product boundaries are <b>not</b> included, although a search of
    /// an XML database reports them: they have no partner, so counting them would measure r over a
    /// population no partner can reach. They are in <see cref="TruncationProductPeptides"/>
    /// instead, so the two together are what the search actually covers.</para>
    /// </summary>
    public int SearchablePeptideCount => _byKey.Count + _ambiguous.Count;

    /// <summary>
    /// The same count for a bare sequence, so the <b>entrapment</b> side of a database can be
    /// measured on the same footing as the target side.
    /// </summary>
    /// <remarks>
    /// Without this a consumer has only <c>entrapmentPeptides</c>, which counts base pieces, so the
    /// achieved ratio is a base-piece ratio rather than the peptide-level <c>r</c> an FDP estimator
    /// is over -- and excision means the two search spaces differ by more than a fold factor, so it
    /// cannot be derived from the target side and a correction.
    /// </remarks>
    public static int CountSearchablePeptides(string sequence, IDigestionParams digestionParams)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            return 0;
        }

        return new EntrapmentPairing(new Protein(sequence, "counting-only"), digestionParams)
            .SearchablePeptideCount;
    }

    /// <summary>The target peptide <paramref name="entrapmentPeptide"/> was built from.</summary>
    /// <param name="entrapmentPeptide">The peptide's BASE sequence. A full sequence with modifications
    /// written into it matches nothing and returns false.</param>
    /// <returns>False when the peptide belongs to no target peptide of this protein, or when its
    /// key is ambiguous.</returns>
    public bool TryResolve(string entrapmentPeptide, out string targetPeptide)
    {
        targetPeptide = string.Empty;
        if (string.IsNullOrEmpty(entrapmentPeptide))
        {
            return false;
        }

        // Through a local, so a miss leaves the documented empty string rather than the null that
        // TryGetValue writes to its out parameter -- this one is declared non-nullable.
        if (!_byKey.TryGetValue(KeyOf(entrapmentPeptide), out string? found))
        {
            return false;
        }

        targetPeptide = found;
        return true;
    }

    /// <summary>
    /// Length, the residues held in place and where, and the sorted free residues -- exactly what a
    /// rearrangement preserves, and nothing it is free to change.
    /// </summary>
    private string KeyOf(string peptide)
    {
        HashSet<int> pinned = DecoySequenceValidator.CleavageSitePositions(peptide, _motifs);

        var pinnedPart = new List<string>();
        var free = new List<char>();
        for (int i = 0; i < peptide.Length; i++)
        {
            if (pinned.Contains(i))
            {
                pinnedPart.Add(i + ":" + peptide[i]);
            }
            else
            {
                free.Add(peptide[i]);
            }
        }

        free.Sort();
        return peptide.Length + "|" + string.Join(",", pinnedPart) + "|" + new string(free.ToArray());
    }
}
