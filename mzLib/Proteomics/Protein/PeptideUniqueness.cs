using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;

namespace Proteomics
{
    /// <summary>
    /// How widely a peptide's sequence is shared across a protein database, with I and L treated as
    /// the same residue. See <see cref="PeptideUniquenessClassifier"/>.
    /// </summary>
    public enum PeptideSharing
    {
        /// <summary>No target protein in the database contains the peptide.</summary>
        NotInDatabase,

        /// <summary>
        /// Every protein that contains the peptide has the same sequence. Two accessions with identical
        /// sequences cannot be told apart by any peptide, so they count as one sequence here.
        /// </summary>
        Unique,

        /// <summary>
        /// The peptide occurs in more than one distinct sequence, and all of them share a gene: the
        /// isoform-level case. Such a peptide supports the gene but not any one isoform.
        /// </summary>
        SharedWithinGene,

        /// <summary>The peptide occurs in sequences with no gene in common.</summary>
        SharedAcrossGenes
    }

    /// <summary>
    /// The classification of one peptide.
    /// </summary>
    /// <param name="Peptide">The peptide exactly as given.</param>
    /// <param name="Sharing">How widely it is shared.</param>
    /// <param name="Accessions">Every target protein containing it, distinct, in ordinal order.</param>
    /// <param name="SharedGeneKeys">
    /// The gene keys common to every protein in <paramref name="Accessions"/>, in ordinal order. Empty
    /// for <see cref="PeptideSharing.NotInDatabase"/> and <see cref="PeptideSharing.SharedAcrossGenes"/>.
    /// </param>
    public sealed record PeptideUniqueness(string Peptide, PeptideSharing Sharing,
        IReadOnlyList<string> Accessions, IReadOnlyList<string> SharedGeneKeys);

    /// <summary>
    /// Classifies peptides as unique to one sequence, shared among isoforms of one gene, or shared
    /// across genes, treating I and L as identical because a mass spectrometer cannot tell them apart.
    ///
    /// A peptide belongs to a protein when the protein's sequence contains it, whatever the protease.
    /// That is deliberately conservative: a peptide is only called unique when no other sequence in
    /// the database contains it anywhere, so an isoform claimed from a unique peptide cannot be
    /// explained by another entry at a site the search's protease rules happened to skip.
    ///
    /// Decoys are ignored. Contaminants are real sequences in the search space and are included, so
    /// a peptide shared with a contaminant is reported as shared.
    /// </summary>
    public static class PeptideUniquenessClassifier
    {
        /// <summary>
        /// The longest prefix packed into one index key: 12 residues at 5 bits each fit in 60 bits.
        /// </summary>
        private const int MaxKeyLength = 12;

        /// <summary>
        /// Classifies each peptide against the target proteins in <paramref name="proteins"/>.
        /// </summary>
        /// <param name="peptides">Unmodified base sequences in upper case. One result is returned per
        /// entry, in the same order, duplicates included.</param>
        /// <param name="proteins">The database. Decoys are skipped.</param>
        /// <param name="geneKeys">The gene keys of a protein; two sequences share a gene when their
        /// keys intersect. Defaults to <see cref="DefaultGeneKeys"/>.</param>
        /// <exception cref="ArgumentNullException">Either collection, or an entry of it, is null.</exception>
        /// <exception cref="ArgumentException">A peptide is empty or holds anything but the letters A-Z.</exception>
        public static IReadOnlyList<PeptideUniqueness> Classify(IEnumerable<string> peptides,
            IEnumerable<Protein> proteins, Func<Protein, IEnumerable<string>> geneKeys = null)
        {
            ArgumentNullException.ThrowIfNull(peptides);
            ArgumentNullException.ThrowIfNull(proteins);
            geneKeys ??= DefaultGeneKeys;

            var given = peptides.ToList();
            var targets = proteins.Select(p => p ?? throw new ArgumentNullException(nameof(proteins),
                    "The protein collection contains a null entry."))
                .Where(p => !p.IsDecoy)
                .ToList();

            // Distinct folded peptides; many inputs can share one (duplicates, or I/L variants).
            var foldedIndex = new Dictionary<string, int>(StringComparer.Ordinal);
            var folded = new List<string>();
            var inputToFolded = new int[given.Count];
            for (int i = 0; i < given.Count; i++)
            {
                string key = FoldPeptide(given[i]);
                if (!foldedIndex.TryGetValue(key, out int id))
                {
                    id = folded.Count;
                    foldedIndex.Add(key, id);
                    folded.Add(key);
                }
                inputToFolded[i] = id;
            }

            var hits = FindContainingProteins(folded, targets);

            var byFolded = new PeptideUniqueness[folded.Count];
            var results = new PeptideUniqueness[given.Count];
            for (int i = 0; i < given.Count; i++)
            {
                int id = inputToFolded[i];
                byFolded[id] ??= Classify(folded[id], hits[id], targets, geneKeys);
                results[i] = byFolded[id] with { Peptide = given[i] };
            }
            return results;
        }

        /// <summary>
        /// The default gene keys, namespaced so that keys of different kinds never collide. A protein
        /// gets every key it has, so an XML entry and a FASTA isoform of it still meet on one:
        /// <list type="bullet">
        /// <item>"ensembl:" + each of <see cref="Protein.EnsemblGeneIds"/> (never a pick among several);</item>
        /// <item>"gene:" + organism + ":" + the primary gene name, so bovine and human ALB stay apart;</item>
        /// <item>"entry:" + the UniProt entry or unversioned RefSeq accession, so P12345-2 meets P12345.
        /// An accession outside both grammars is used verbatim.</item>
        /// </list>
        /// </summary>
        public static IEnumerable<string> DefaultGeneKeys(Protein protein)
        {
            ArgumentNullException.ThrowIfNull(protein);

            var keys = protein.EnsemblGeneIds.Select(geneId => "ensembl:" + geneId).ToList();

            string primary = protein.GeneNames?.FirstOrDefault(n => n.Item1 == "primary")?.Item2;
            if (!string.IsNullOrEmpty(primary))
            {
                keys.Add("gene:" + (protein.Organism ?? "") + ":" + primary);
            }

            keys.Add("entry:" + ProteinAccession.Parse(protein.Accession).EntryAccession);
            return keys;
        }

        private static PeptideUniqueness Classify(string foldedPeptide, List<int> proteinIndexes,
            List<Protein> targets, Func<Protein, IEnumerable<string>> geneKeys)
        {
            var containing = proteinIndexes.Select(i => targets[i]).ToList();
            var accessions = containing.Select(p => p.Accession)
                .Distinct(StringComparer.Ordinal)
                .OrderBy(a => a, StringComparer.Ordinal)
                .ToList();

            if (containing.Count == 0)
            {
                return new PeptideUniqueness(foldedPeptide, PeptideSharing.NotInDatabase, accessions, Array.Empty<string>());
            }

            HashSet<string> shared = null;
            foreach (var protein in containing)
            {
                var keys = (geneKeys(protein) ?? Enumerable.Empty<string>()).Where(k => k != null);
                if (shared == null)
                {
                    shared = new HashSet<string>(keys, StringComparer.Ordinal);
                }
                else
                {
                    shared.IntersectWith(keys);
                }
            }
            var sharedKeys = shared.OrderBy(k => k, StringComparer.Ordinal).ToList();

            int distinctSequences = containing.Select(p => FoldSequence(p.BaseSequence))
                .Distinct(StringComparer.Ordinal)
                .Count();

            PeptideSharing sharing = distinctSequences == 1 ? PeptideSharing.Unique
                : sharedKeys.Count > 0 ? PeptideSharing.SharedWithinGene
                : PeptideSharing.SharedAcrossGenes;

            return new PeptideUniqueness(foldedPeptide, sharing, accessions,
                sharing == PeptideSharing.SharedAcrossGenes ? Array.Empty<string>() : sharedKeys);
        }

        /// <summary>
        /// For each folded peptide, the indexes of the target proteins containing it. Peptides are
        /// indexed by a packed prefix of the shortest peptide's length, so each protein position costs
        /// one rolling update and one lookup rather than a scan over every peptide.
        /// </summary>
        private static List<int>[] FindContainingProteins(List<string> folded, List<Protein> targets)
        {
            var hits = new List<int>[folded.Count];
            for (int i = 0; i < hits.Length; i++)
            {
                hits[i] = new List<int>();
            }
            if (folded.Count == 0)
            {
                return hits;
            }

            int keyLength = Math.Min(MaxKeyLength, folded.Min(p => p.Length));
            ulong mask = (1UL << (5 * keyLength)) - 1;

            var index = new Dictionary<ulong, List<int>>();
            for (int id = 0; id < folded.Count; id++)
            {
                ulong key = 0;
                for (int j = 0; j < keyLength; j++)
                {
                    key = (key << 5) | Code(folded[id][j]);
                }
                if (!index.TryGetValue(key, out var ids))
                {
                    index.Add(key, ids = new List<int>());
                }
                ids.Add(id);
            }

            for (int p = 0; p < targets.Count; p++)
            {
                string sequence = FoldSequence(targets[p].BaseSequence);
                ulong key = 0;
                int valid = 0; // residues since the last character that cannot be in a peptide
                for (int end = 0; end < sequence.Length; end++)
                {
                    ulong code = Code(sequence[end]);
                    if (code == 0)
                    {
                        valid = 0;
                        key = 0;
                        continue;
                    }
                    key = ((key << 5) | code) & mask;
                    if (++valid < keyLength || !index.TryGetValue(key, out var candidates))
                    {
                        continue;
                    }

                    int start = end - keyLength + 1;
                    foreach (int id in candidates)
                    {
                        string peptide = folded[id];
                        var peptideHits = hits[id];
                        if (start + peptide.Length <= sequence.Length
                            && (peptideHits.Count == 0 || peptideHits[^1] != p)
                            && sequence.AsSpan(start, peptide.Length).SequenceEqual(peptide))
                        {
                            peptideHits.Add(p);
                        }
                    }
                }
            }
            return hits;
        }

        /// <summary>A letter's 5-bit code (1-26), or 0 for anything a peptide cannot contain.</summary>
        private static ulong Code(char c) => c is >= 'A' and <= 'Z' ? (ulong)(c - 'A' + 1) : 0;

        private static string FoldSequence(string sequence) => (sequence ?? "").Replace('I', 'L');

        private static string FoldPeptide(string peptide)
        {
            if (peptide == null)
            {
                throw new ArgumentNullException(nameof(peptide), "The peptide collection contains a null entry.");
            }
            if (peptide.Length == 0 || peptide.Any(c => c is < 'A' or > 'Z'))
            {
                throw new ArgumentException(
                    $"Peptide '{peptide}' must be a non-empty unmodified base sequence of the letters A-Z.",
                    nameof(peptide));
            }
            return peptide.Replace('I', 'L');
        }
    }
}
