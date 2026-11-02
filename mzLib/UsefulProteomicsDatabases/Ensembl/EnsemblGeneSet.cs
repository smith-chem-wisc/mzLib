using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Security.Cryptography;
using System.Text.RegularExpressions;

namespace UsefulProteomicsDatabases.Ensembl
{
    /// <summary>
    /// One gene row of an Ensembl GTF.
    /// </summary>
    /// <param name="GeneId">The stable id, e.g. "ENSG00000111640". Ensembl GTFs write it unversioned.</param>
    /// <param name="Version">The gene_version attribute, or null when absent. Absent is not 0.</param>
    /// <param name="Biotype">The gene_biotype attribute, or "unknown" when absent.</param>
    /// <param name="Symbol">The gene_name attribute -- a display label, never a key. Null when absent,
    /// which is common for novel lncRNAs.</param>
    /// <param name="SeqRegion">The GTF's first column: the chromosome or scaffold.</param>
    public sealed record EnsemblGene(string GeneId, int? Version, string Biotype, string Symbol, string SeqRegion);

    /// <summary>
    /// The gene set of one Ensembl release, read from its GTF, together with what it was read from.
    ///
    /// It exists to answer one question the protein database cannot: is this gene on the assembly a
    /// resolution is being counted against? Load the PRIMARY-ASSEMBLY GTF
    /// (Species.Assembly.Release.gtf.gz, not the chr_patch_hapl_scaff one). Ensembl's cross-references
    /// also name genes on ALT haplotypes and patches, where one locus is described many times; counting
    /// those as separate genes made human release 116 look 6.99% multi-gene when it is 0.36%.
    ///
    /// Reference data this size (tens of MB) is never embedded in the assembly -- see
    /// ControlledVocabulary for the line -- so the caller supplies the file and the set records its name,
    /// sha256, release and genome build, which travel with anything resolved against it.
    ///
    /// A GTF is large (human release 116: 141 MB compressed, 4.66 GB unzipped) and only its gene rows
    /// are used. <see cref="EnsemblGeneSetWriter"/> writes those rows and the GTF's provenance to a
    /// compact table (about 0.5 MB for human) that <see cref="EnsemblGeneSetReader"/> reads back into
    /// an identical set.
    /// </summary>
    public sealed class EnsemblGeneSet
    {
        private static readonly Regex GeneIdAttribute = new(@"gene_id ""([^""]+)""", RegexOptions.Compiled);
        private static readonly Regex GeneVersionAttribute = new(@"gene_version ""(\d+)""", RegexOptions.Compiled);
        private static readonly Regex GeneBiotypeAttribute = new(@"gene_biotype ""([^""]+)""", RegexOptions.Compiled);
        private static readonly Regex GeneNameAttribute = new(@"gene_name ""([^""]+)""", RegexOptions.Compiled);

        /// <summary>Ensembl names GTFs Species.Assembly.Release.gtf[.gz]; the release is the number before ".gtf".</summary>
        private static readonly Regex ReleaseInFileName = new(@"\.(\d+)\.gtf(\.gz)?$", RegexOptions.Compiled | RegexOptions.IgnoreCase);

        private readonly Dictionary<string, EnsemblGene> _genes;

        private EnsemblGeneSet(Dictionary<string, EnsemblGene> genes, string sourceFileName, string sourceSha256,
            string release, string genomeBuild, string genebuildLastUpdated)
        {
            _genes = genes;
            SourceFileName = sourceFileName;
            SourceSha256 = sourceSha256;
            Release = release;
            GenomeBuild = genomeBuild;
            GenebuildLastUpdated = genebuildLastUpdated;
            GeneIds = genes.Keys.OrderBy(k => k, StringComparer.Ordinal).ToList();
        }

        /// <summary>The file name the set was read from.</summary>
        public string SourceFileName { get; }

        /// <summary>
        /// Lower-case hex sha256 of the file's bytes as read -- the compressed bytes for a .gz, which is
        /// what Ensembl publishes and checksums.
        /// </summary>
        public string SourceSha256 { get; }

        /// <summary>The Ensembl release parsed from the file name, or null when the name carries none.</summary>
        public string Release { get; }

        /// <summary>The "#!genome-build" header value, e.g. "GRCh38.p14", or null when absent.</summary>
        public string GenomeBuild { get; }

        /// <summary>The "#!genebuild-last-updated" header value, or null when absent.</summary>
        public string GenebuildLastUpdated { get; }

        /// <summary>Every stable gene id in the set, in ordinal order.</summary>
        public IReadOnlyList<string> GeneIds { get; }

        public int Count => _genes.Count;

        /// <summary>Every gene in the set, in ordinal order of stable id.</summary>
        public IEnumerable<EnsemblGene> Genes => GeneIds.Select(id => _genes[id]);

        /// <summary>True when the stable id is in the set. Versioned ids are not stripped here: that is the caller's call.</summary>
        public bool Contains(string geneId) => geneId != null && _genes.ContainsKey(geneId);

        public bool TryGetGene(string geneId, out EnsemblGene gene)
        {
            gene = null;
            return geneId != null && _genes.TryGetValue(geneId, out gene);
        }

        /// <summary>
        /// Reads the gene rows of a GTF (plain or .gz). Transcript, exon and other feature rows are
        /// skipped. A gene row without a gene_id is malformed and throws rather than being skipped,
        /// because a silently missing gene would later read as "not on the assembly".
        /// </summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">A gene row carries no gene_id.</exception>
        public static EnsemblGeneSet LoadGtf(string path)
        {
            if (!File.Exists(path))
            {
                throw new FileNotFoundException("Ensembl GTF not found.", path);
            }

            string sha256;
            using (var hashStream = File.OpenRead(path))
            {
                sha256 = Convert.ToHexString(SHA256.HashData(hashStream)).ToLowerInvariant();
            }

            var genes = new Dictionary<string, EnsemblGene>(StringComparer.Ordinal);
            string genomeBuild = null;
            string lastUpdated = null;

            using var file = File.OpenRead(path);
            using Stream content = path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase)
                ? new GZipStream(file, CompressionMode.Decompress)
                : file;
            using var reader = new StreamReader(content);

            int lineNumber = 0;
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                lineNumber++;
                if (line.StartsWith("#", StringComparison.Ordinal))
                {
                    genomeBuild ??= HeaderValue(line, "#!genome-build ");
                    lastUpdated ??= HeaderValue(line, "#!genebuild-last-updated ");
                    continue;
                }

                // Only the first 8 columns are fixed; attributes are the 9th.
                string[] columns = line.Split('\t', 9);
                if (columns.Length < 9 || columns[2] != "gene")
                {
                    continue;
                }

                string attributes = columns[8];
                var id = GeneIdAttribute.Match(attributes);
                if (!id.Success)
                {
                    throw new InvalidDataException($"{Path.GetFileName(path)} line {lineNumber}: gene row has no gene_id.");
                }

                var version = GeneVersionAttribute.Match(attributes);
                var biotype = GeneBiotypeAttribute.Match(attributes);
                var name = GeneNameAttribute.Match(attributes);

                genes[id.Groups[1].Value] = new EnsemblGene(
                    id.Groups[1].Value,
                    version.Success ? int.Parse(version.Groups[1].Value) : null,
                    biotype.Success ? biotype.Groups[1].Value : "unknown",
                    name.Success ? name.Groups[1].Value : null,
                    columns[0]);
            }

            var release = ReleaseInFileName.Match(Path.GetFileName(path));
            return new EnsemblGeneSet(genes, Path.GetFileName(path), sha256,
                release.Success ? release.Groups[1].Value : null, genomeBuild, lastUpdated);
        }

        /// <summary>
        /// Rebuilds a set from genes and the provenance of the GTF they were originally read from. Used by
        /// <see cref="EnsemblGeneSetReader"/>: the provenance stays the GTF's, so anything resolved against
        /// the rebuilt set is keyed exactly as it would be against the GTF.
        /// </summary>
        internal static EnsemblGeneSet FromGenes(Dictionary<string, EnsemblGene> genes, string sourceFileName,
            string sourceSha256, string release, string genomeBuild, string genebuildLastUpdated) =>
            new(genes, sourceFileName, sourceSha256, release, genomeBuild, genebuildLastUpdated);

        private static string HeaderValue(string line, string prefix) =>
            line.StartsWith(prefix, StringComparison.Ordinal) ? line.Substring(prefix.Length).Trim() : null;
    }
}
