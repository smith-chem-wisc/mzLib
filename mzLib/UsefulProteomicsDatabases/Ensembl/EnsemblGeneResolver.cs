using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MzLibUtil;
using Omics.BioPolymerGroup;
using Proteomics;

namespace UsefulProteomicsDatabases.Ensembl
{
    /// <summary>
    /// How an accession met the gene set. Every accession gets exactly one, so "we did not look",
    /// "the source has nothing" and "the source knows it only on an ALT haplotype" never share a null.
    /// </summary>
    public enum GeneResolutionOutcome
    {
        /// <summary>Exactly one gene on the assembly the gene set covers.</summary>
        Resolved,

        /// <summary>More than one gene on that assembly. One row per gene; never a pick.</summary>
        MultiGene,

        /// <summary>
        /// The source links the accession only to genes outside the gene set -- with the primary-assembly
        /// GTF, ALT haplotypes, patches and scaffolds. The source has an answer, just not one on the
        /// assembly every other row is counted against.
        /// </summary>
        OffPrimaryOnly,

        /// <summary>A well-formed accession the source links to no gene.</summary>
        NotInSource,

        /// <summary>
        /// Not a UniProt or RefSeq accession, and the source links it to nothing -- a decoy or
        /// entrapment prefix, or a mangled id. Kept apart from NotInSource so a data problem is not read
        /// as a biology result.
        /// </summary>
        UnrecognizedAccession,

        /// <summary>
        /// Flagged a contaminant. Never mapped, even when it would resolve, and never silently dropped.
        /// </summary>
        ContaminantNotMapped
    }

    /// <summary>
    /// One row of a resolution: one (accession, gene) pair, or one outcome row for an accession with
    /// no gene. Restrictions are columns, never modes -- a consumer filters on Outcome and can see what
    /// it dropped.
    /// </summary>
    public sealed record GeneResolution(
        string Accession,
        string EntryAccession,
        int? Isoform,
        AccessionNamespace Namespace,
        GeneResolutionOutcome Outcome,
        int GeneCount,
        string GeneId,
        string VersionedGeneId,
        string GeneSymbol,
        string GeneBiotype,
        int OffPrimaryGenes,
        string UniProtGeneName,
        string Source,
        string SearchDatabaseSha256,
        string GeneSetRelease,
        string GeneSetSha256);

    /// <summary>
    /// Resolves proteins to stable Ensembl gene ids, counted against one release's gene set.
    ///
    /// The gene links are read from the search database itself (Protein.EnsemblGeneReferences), so a
    /// resolution is keyed on (accession, search database sha256) and never on a display symbol. A
    /// protein-group TSV's Gene column is the first gene name per protein, '|'-joined; it can be empty
    /// for a protein, it cannot express several genes, and it differs by database format. This is the
    /// replacement for reading it.
    /// </summary>
    public sealed class EnsemblGeneResolver
    {
        /// <summary>The Source value for rows resolved from the search database's own dbReferences.</summary>
        public const string SearchDatabaseSource = "search_database_dbreference";

        public EnsemblGeneResolver(EnsemblGeneSet geneSet)
        {
            GeneSet = geneSet ?? throw new ArgumentNullException(nameof(geneSet));
        }

        /// <summary>The gene set every resolution is counted against.</summary>
        public EnsemblGeneSet GeneSet { get; }

        /// <summary>
        /// Every row for one protein. Never empty: a protein always gets at least one outcome.
        /// </summary>
        /// <param name="protein">A protein as loaded from the search database.</param>
        /// <param name="searchDatabaseSha256">The sha256 of the decompressed database the search read.</param>
        public IReadOnlyList<GeneResolution> Resolve(Protein protein, string searchDatabaseSha256)
        {
            ArgumentNullException.ThrowIfNull(protein);

            // A sequence-variant proteoform ("P12345_S70N", as LoadProteinXML names applied variants) is
            // not an accession grammar, but it is a known entry: the entry, its grammar and its gene
            // links come from the consensus, while the verbatim accession is kept as the search reported it.
            var entry = protein.ConsensusVariant as Protein ?? protein;
            var accession = entry.Accession.ParseProteinAccession() with { Verbatim = protein.Accession };
            string uniProtGeneName = PrimaryGeneName(protein) ?? PrimaryGeneName(entry);

            GeneResolution Row(GeneResolutionOutcome outcome, int geneCount = 0, string geneId = null,
                string versionedGeneId = null, int offPrimary = 0)
            {
                EnsemblGene gene = null;
                if (geneId != null) GeneSet.TryGetGene(geneId, out gene);
                return new GeneResolution(accession.Verbatim, accession.EntryAccession, accession.Isoform,
                    accession.Namespace, outcome, geneCount, geneId, versionedGeneId, gene?.Symbol, gene?.Biotype,
                    offPrimary, uniProtGeneName, SearchDatabaseSource, searchDatabaseSha256, GeneSet.Release,
                    GeneSet.SourceSha256);
            }

            if (protein.IsContaminant)
            {
                return new[] { Row(GeneResolutionOutcome.ContaminantNotMapped) };
            }

            // First versioned id seen per stable gene: one gene's transcripts share a gene version.
            var versionedByGene = new Dictionary<string, string>(StringComparer.Ordinal);
            var links = protein.EnsemblGeneReferences;
            if (links.Count == 0 && !ReferenceEquals(entry, protein))
            {
                links = entry.EnsemblGeneReferences;
            }

            foreach (var link in links)
            {
                versionedByGene.TryAdd(link.GeneId, link.VersionedGeneId);
            }

            if (versionedByGene.Count == 0)
            {
                return new[]
                {
                    Row(accession.Namespace == AccessionNamespace.Unrecognized
                        ? GeneResolutionOutcome.UnrecognizedAccession
                        : GeneResolutionOutcome.NotInSource)
                };
            }

            var onAssembly = versionedByGene.Keys.Where(GeneSet.Contains).OrderBy(g => g, StringComparer.Ordinal).ToList();
            int offAssembly = versionedByGene.Count - onAssembly.Count;
            if (onAssembly.Count == 0)
            {
                return new[] { Row(GeneResolutionOutcome.OffPrimaryOnly, offPrimary: offAssembly) };
            }

            var outcome = onAssembly.Count == 1 ? GeneResolutionOutcome.Resolved : GeneResolutionOutcome.MultiGene;
            return onAssembly
                .Select(g => Row(outcome, onAssembly.Count, g, versionedByGene[g], offAssembly))
                .ToList();
        }

        /// <summary>Resolves every protein, in order.</summary>
        public IEnumerable<GeneResolution> ResolveAll(IEnumerable<Protein> proteins, string searchDatabaseSha256)
        {
            ArgumentNullException.ThrowIfNull(proteins);
            foreach (var protein in proteins)
            {
                foreach (var row in Resolve(protein, searchDatabaseSha256))
                {
                    yield return row;
                }
            }
        }

        /// <summary>
        /// The entry's primary gene name, as the search database wrote it. A label only -- read once from
        /// bytes pinned by the database hash, so it cannot come out ragged across datasets.
        /// </summary>
        private static string PrimaryGeneName(Protein protein) =>
            protein.GeneNames?.FirstOrDefault(n => n.Item1 == "primary")?.Item2;
    }

    /// <summary>
    /// Writes resolutions as a long-format table: one row per (accession, gene), never a '|'-joined
    /// cell. Column names are snake_case because this is a data-interchange table whose names are a
    /// contract with its consumers. Nulls are empty cells; every row still carries an outcome saying why.
    /// </summary>
    public static class GeneResolutionTsv
    {
        public static readonly IReadOnlyList<TsvColumn<GeneResolution>> Schema = new List<TsvColumn<GeneResolution>>
        {
            new("accession", r => r.Accession),
            new("entry_accession", r => r.EntryAccession),
            new("isoform", r => r.Isoform?.ToString()),
            new("namespace", r => r.Namespace.ToString().ToLowerInvariant()),
            new("outcome", r => OutcomeName(r.Outcome)),
            new("n_genes", r => r.GeneCount.ToString()),
            new("gene_id", r => r.GeneId),
            new("versioned_gene_id", r => r.VersionedGeneId),
            new("gene_symbol", r => r.GeneSymbol),
            new("gene_biotype", r => r.GeneBiotype),
            new("off_primary_genes", r => r.OffPrimaryGenes.ToString()),
            new("uniprot_gene_name", r => r.UniProtGeneName),
            new("source", r => r.Source),
            new("search_database_sha256", r => r.SearchDatabaseSha256),
            new("gene_set_release", r => r.GeneSetRelease),
            new("gene_set_sha256", r => r.GeneSetSha256),
        };

        public static void Write(TextWriter output, IEnumerable<GeneResolution> rows) =>
            TsvWriter.Write(output, Schema, rows);

        /// <summary>The outcome as written in the table, e.g. "off_primary_only".</summary>
        public static string OutcomeName(GeneResolutionOutcome outcome) => outcome switch
        {
            GeneResolutionOutcome.Resolved => "resolved",
            GeneResolutionOutcome.MultiGene => "multi_gene",
            GeneResolutionOutcome.OffPrimaryOnly => "off_primary_only",
            GeneResolutionOutcome.NotInSource => "not_in_source",
            GeneResolutionOutcome.UnrecognizedAccession => "unrecognized_accession",
            GeneResolutionOutcome.ContaminantNotMapped => "contaminant_not_mapped",
            _ => throw new ArgumentOutOfRangeException(nameof(outcome), outcome, null)
        };
    }
}
