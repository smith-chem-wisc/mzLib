using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using MzLibUtil;
using Proteomics;
using UsefulProteomicsDatabases.Annotation;

namespace UsefulProteomicsDatabases.Pathways
{
    /// <summary>What one entry's Reactome rows say about it.</summary>
    public enum ReactomeMembershipStatus
    {
        /// <summary>The row carries a pathway, the entry's own or inherited.</summary>
        Annotated = 0,

        /// <summary>
        /// UniProt lists no Reactome reference for the entry, and it inherits none, while its organism is covered.
        /// This means "UniProt lists none", not "in no pathway".
        /// </summary>
        NoReference = 1,

        /// <summary>No entry of the entry's organism in the loaded list carries any Reactome reference.</summary>
        SpeciesNotCovered = 2
    }

    /// <summary>
    /// One row of a database's Reactome membership: an entry and one pathway it is a member of, or the entry's single
    /// pathway-less row.
    /// </summary>
    /// <param name="Accession">The database entry's accession, verbatim.</param>
    /// <param name="Organism">The entry's own NCBI taxon as a CURIE, "NCBITaxon:9606"; null when the entry names none.</param>
    /// <param name="PathwayId">The Reactome stable id, or null on a pathway-less row.</param>
    /// <param name="Molecule">The isoform the membership is scoped to (UniProt's &lt;molecule&gt;), or null when it covers
    /// the whole entry.</param>
    /// <param name="SequenceAccession">Null on a whole-entry row. On an isoform row: the entry's own sequence accession
    /// (<see cref="Annotation.MemberResolver.SequenceAccessionOf"/>) when the molecule is the entry's displayed isoform,
    /// otherwise the molecule. A member takes the isoform rows whose sequence_accession is its own.</param>
    /// <param name="Inherited">True when the entry is a UniProt isoform with no reference of its own and the row is
    /// its entry's whole-entry row; null on a pathway-less row.</param>
    /// <param name="Status">Annotated on every row that carries a pathway.</param>
    public sealed record ReactomeMembershipRow(string Accession, string Organism, string PathwayId, string Molecule,
        string SequenceAccession, bool? Inherited, ReactomeMembershipStatus Status);

    /// <summary>
    /// The Reactome membership of one annotation database, independent of any search
    /// (kegg:DEF-REACTOME-MEMBERSHIP v1). For every entry: the pathways UniProt lists on it, each scoped to the whole
    /// entry or to one isoform, after a species check. Ancestors are not here; they are a property of the Reactome
    /// release (the per-release Reactome tables). Which entry a protein-group member takes is
    /// <see cref="Annotation.MemberResolver"/>'s.
    ///
    /// The rules, in order:
    /// <list type="number">
    /// <item>Only references whose type is exactly "Reactome" are read (<see cref="Protein.ReactomePathways"/>).</item>
    /// <item>Species. An entry's organism is its NCBI taxon. An organism's code is the id prefix ("HSA" of R-HSA-70171)
    /// carried by a strict majority of its references, counted once per entry and pathway, when it has at least
    /// <see cref="MinimumReferencesForSpeciesCode"/>. With a code, a reference with another prefix is dropped and
    /// counted. Without one, the organism's references are kept unfiltered and it is listed as unvalidated. Entries
    /// that name no taxon are one organism, <see cref="UnknownOrganism"/>.</item>
    /// <item>Isoform scope. A reference scoped to an isoform gives one row per isoform; see
    /// <see cref="ReactomeMembershipRow.SequenceAccession"/>.</item>
    /// <item>Inheritance (as GO's). A UniProt isoform with no Reactome reference of its own, whose entry is in the
    /// database with at least one whole-entry row, takes those rows, flagged inherited. Isoform rows are never
    /// copied: they reach the isoform's members through sequence_accession.</item>
    /// <item>An accession present more than once has its references unioned.</item>
    /// <item>Every entry has at least one row: an entry with no pathway gets one pathway-less row, reading
    /// <see cref="ReactomeMembershipStatus.SpeciesNotCovered"/> when no entry of its organism carries any Reactome
    /// reference, otherwise <see cref="ReactomeMembershipStatus.NoReference"/>.</item>
    /// </list>
    /// </summary>
    public sealed class ReactomeAccessionAnnotation
    {
        /// <summary>A Reactome stable id, R-HSA-70171, and its species prefix.</summary>
        private static readonly Regex StableId = new(@"^R-(?<species>[A-Z]{3})-\d+$", RegexOptions.Compiled);

        private readonly ReactomeMembershipRow[] _rows;

        /// <summary>The definition every row follows, written in the header.</summary>
        public const string DefinitionId = "kegg:DEF-REACTOME-MEMBERSHIP";

        /// <summary>The definition's version. It changes only when the same inputs would give different rows.</summary>
        public const int DefinitionVersion = 1;

        /// <summary>An organism needs at least this many Reactome references for its species code to be read.</summary>
        public const int MinimumReferencesForSpeciesCode = 10;

        /// <summary>
        /// The name the species rule and the header give the entries that name no NCBI taxon (a FASTA entry without
        /// OX=). They are treated as one organism; their rows' organism is null.
        /// </summary>
        public const string UnknownOrganism = "unknown";

        /// <param name="annotationProteins">The target proteins of one search, every file it loaded. Decoys and nulls are
        /// skipped. An accession present more than once has its references unioned.</param>
        /// <param name="annotationDbSha256">The target-set hash the artefact is keyed by.</param>
        /// <exception cref="ArgumentNullException">annotationProteins is null.</exception>
        /// <exception cref="ArgumentException">annotationDbSha256 is null or empty: the table is per database.</exception>
        public ReactomeAccessionAnnotation(IEnumerable<Protein> annotationProteins, string annotationDbSha256)
        {
            ArgumentNullException.ThrowIfNull(annotationProteins);
            if (string.IsNullOrEmpty(annotationDbSha256))
            {
                throw new ArgumentException("The annotation database's sha256 is required: the table is per database.",
                    nameof(annotationDbSha256));
            }
            AnnotationDbSha256 = annotationDbSha256;

            var entries = ReadEntries(annotationProteins);
            Accessions = entries.Keys.ToArray();

            // Species: count each organism's references, read its code, then drop the references a code rules out.
            var counts = new Dictionary<string, Dictionary<string, int>>(StringComparer.Ordinal);
            foreach (var entry in entries.Values)
            {
                if (!counts.TryGetValue(entry.SpeciesKey, out var prefixes))
                {
                    prefixes = new Dictionary<string, int>(StringComparer.Ordinal);
                    counts.Add(entry.SpeciesKey, prefixes);
                }
                foreach (var pathway in entry.OwnPathways)
                {
                    string prefix = PrefixOf(pathway.Id);
                    prefixes[prefix] = prefixes.GetValueOrDefault(prefix) + 1;
                }
            }
            var codes = new SortedDictionary<string, string>(StringComparer.Ordinal);
            var unvalidated = new SortedDictionary<string, int>(StringComparer.Ordinal);
            foreach (var (organism, prefixes) in counts)
            {
                int total = prefixes.Values.Sum();
                if (total == 0)
                {
                    continue;
                }
                var leader = prefixes.OrderByDescending(p => p.Value).First();
                if (total >= MinimumReferencesForSpeciesCode && leader.Key.Length > 0 && 2 * leader.Value > total)
                {
                    codes.Add(organism, leader.Key);
                }
                else
                {
                    unvalidated.Add(organism, total);
                }
            }
            SpeciesCodes = codes;
            UnvalidatedOrganisms = unvalidated;

            int dropped = 0;
            foreach (var entry in entries.Values)
            {
                entry.KeptPathways = entry.OwnPathways;
                if (codes.TryGetValue(entry.SpeciesKey, out string code))
                {
                    entry.KeptPathways = entry.OwnPathways.Where(p => PrefixOf(p.Id) == code).ToList();
                    dropped += entry.OwnPathways.Count - entry.KeptPathways.Count;
                }
            }
            DroppedForeignReferences = dropped;

            _rows = entries.Values
                .SelectMany(entry => RowsOf(entry, entries, counts[entry.SpeciesKey].Count > 0))
                .ToArray();
        }

        /// <summary>The target-set hash, as the caller gave it.</summary>
        public string AnnotationDbSha256 { get; }

        /// <summary>Every entry, in ordinal order.</summary>
        public IReadOnlyList<string> Accessions { get; }

        /// <summary>Each organism whose code was read (a strict majority of at least
        /// <see cref="MinimumReferencesForSpeciesCode"/> references), and its code, e.g. "HSA". In ordinal order.</summary>
        public IReadOnlyDictionary<string, string> SpeciesCodes { get; }

        /// <summary>Each organism with references but no code, and its reference count. Its references are kept
        /// unfiltered. In ordinal order.</summary>
        public IReadOnlyDictionary<string, int> UnvalidatedOrganisms { get; }

        /// <summary>The references dropped because their prefix is not their organism's code.</summary>
        public int DroppedForeignReferences { get; }

        /// <summary>Every row: entries in ordinal order, then pathway id, then molecule (whole-entry first).</summary>
        public IReadOnlyList<ReactomeMembershipRow> Rows() => _rows;

        /// <summary>One accession of the database, its copies merged.</summary>
        private sealed class Entry
        {
            public Entry(string accession) => Accession = accession;

            public string Accession { get; }

            /// <summary>The first NCBI taxon its copies name, as a CURIE; null when none names one.</summary>
            public string Organism { get; set; }

            /// <summary>The organism the species rule counts it under.</summary>
            public string SpeciesKey => Organism ?? UnknownOrganism;

            /// <summary>The first displayed isoform its copies name.</summary>
            public string DisplayedIsoformId { get; set; }

            public List<DatabaseReference> References { get; } = new();

            /// <summary>Its own pathways, before the species check.</summary>
            public IReadOnlyList<ReactomePathway> OwnPathways { get; set; }

            /// <summary>Its own pathways after the species check.</summary>
            public IReadOnlyList<ReactomePathway> KeptPathways { get; set; }
        }

        private static SortedDictionary<string, Entry> ReadEntries(IEnumerable<Protein> proteins)
        {
            var entries = new SortedDictionary<string, Entry>(StringComparer.Ordinal);
            foreach (var protein in proteins)
            {
                if (protein == null || protein.IsDecoy)
                {
                    continue;
                }
                if (!entries.TryGetValue(protein.Accession, out var entry))
                {
                    entry = new Entry(protein.Accession);
                    entries.Add(protein.Accession, entry);
                }
                string taxon = protein.NcbiTaxonomyId;
                entry.Organism ??= string.IsNullOrEmpty(taxon) ? null : "NCBITaxon:" + taxon;
                entry.DisplayedIsoformId ??= protein.DisplayedIsoformId;
                entry.References.AddRange(protein.DatabaseReferences);
            }
            foreach (var entry in entries.Values)
            {
                entry.OwnPathways = ReactomePathway.FromDatabaseReferences(entry.References);
            }
            return entries;
        }

        /// <summary>The species prefix of a stable id ("HSA" of R-HSA-70171), or "" when the id is not one.</summary>
        private static string PrefixOf(string pathwayId)
        {
            var match = StableId.Match(pathwayId);
            return match.Success ? match.Groups["species"].Value : "";
        }

        private static IEnumerable<ReactomeMembershipRow> RowsOf(Entry entry, SortedDictionary<string, Entry> entries,
            bool organismCovered)
        {
            var rows = new List<ReactomeMembershipRow>();
            foreach (var pathway in entry.KeptPathways)
            {
                if (pathway.Molecules.Count == 0)
                {
                    rows.Add(new ReactomeMembershipRow(entry.Accession, entry.Organism, pathway.Id, null, null, false,
                        ReactomeMembershipStatus.Annotated));
                    continue;
                }
                foreach (string molecule in pathway.Molecules)
                {
                    string sequenceAccession = molecule == entry.DisplayedIsoformId
                        ? MemberResolver.SequenceAccessionOf(entry.Accession)
                        : molecule;
                    rows.Add(new ReactomeMembershipRow(entry.Accession, entry.Organism, pathway.Id, molecule,
                        sequenceAccession, false, ReactomeMembershipStatus.Annotated));
                }
            }

            if (entry.OwnPathways.Count == 0
                && ProteinAccession.Parse(entry.Accession) is { Namespace: AccessionNamespace.UniProt, Isoform: not null } isoform
                && entries.TryGetValue(isoform.EntryAccession, out var parent))
            {
                rows.AddRange(parent.KeptPathways
                    .Where(p => p.Molecules.Count == 0)
                    .Select(p => new ReactomeMembershipRow(entry.Accession, entry.Organism, p.Id, null, null, true,
                        ReactomeMembershipStatus.Annotated)));
            }

            if (rows.Count == 0)
            {
                rows.Add(new ReactomeMembershipRow(entry.Accession, entry.Organism, null, null, null, null,
                    organismCovered ? ReactomeMembershipStatus.NoReference : ReactomeMembershipStatus.SpeciesNotCovered));
            }
            return rows
                .OrderBy(r => r.PathwayId, StringComparer.Ordinal)
                .ThenBy(r => r.Molecule ?? "", StringComparer.Ordinal);
        }
    }
}
