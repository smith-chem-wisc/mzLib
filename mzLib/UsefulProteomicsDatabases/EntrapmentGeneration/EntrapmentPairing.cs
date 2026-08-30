#nullable enable
using MzLibUtil;
using Omics.Digestion;
using Proteomics;
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

    /// <summary>Indexes a target protein's peptides, missed cleavages included, by composition and pinned pattern.</summary>
    public EntrapmentPairing(Protein target, IDigestionParams digestionParams)
    {
        if (target is null)
        {
            throw new MzLibException("Cannot build a pairing without a target protein.");
        }

        _byKey = new Dictionary<string, string>();
        _ambiguous = new HashSet<string>();
        var collided = new HashSet<string>();

        _motifs = digestionParams.DigestionAgent.DigestionMotifs;
        List<int> sites = digestionParams.DigestionAgent.GetDigestionSiteIndices(target.BaseSequence);
        sites.Sort();

        // Index runs of adjacent base pieces, not just single ones: a search reports missed-cleavage
        // peptides too, and those are pairable for the same reason -- base-piece compositions add, so
        // a run is isomeric with the run its partner was built from.
        int maxPieces = digestionParams.MaxMissedCleavages + 1;

        for (int first = 0; first < sites.Count - 1; first++)
        {
            for (int pieces = 1; pieces <= maxPieces && first + pieces < sites.Count; pieces++)
            {
                int start = sites[first];
                string peptide = target.BaseSequence.Substring(start, sites[first + pieces] - start);
                string key = KeyOf(peptide);

                if (_byKey.TryGetValue(key, out string? existing))
                {
                    if (existing != peptide)
                    {
                        collided.Add(key);
                        _ambiguous.Add(existing);
                        _ambiguous.Add(peptide);
                    }
                    continue;
                }

                _byKey[key] = peptide;
            }
        }

        foreach (string key in collided)
        {
            _byKey.Remove(key);
        }
    }

    /// <summary>
    /// Target peptides sharing a composition-and-pinning key with another peptide of the same
    /// protein, and therefore not resolvable. Report these rather than pairing them.
    /// </summary>
    public IReadOnlyCollection<string> AmbiguousPeptides => _ambiguous;

    /// <summary>The target peptide <paramref name="entrapmentPeptide"/> was built from.</summary>
    /// <returns>False when the peptide belongs to no target peptide of this protein, or when its
    /// key is ambiguous.</returns>
    public bool TryResolve(string entrapmentPeptide, out string targetPeptide)
    {
        targetPeptide = string.Empty;
        if (string.IsNullOrEmpty(entrapmentPeptide))
        {
            return false;
        }

        return _byKey.TryGetValue(KeyOf(entrapmentPeptide), out targetPeptide!);
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
