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
    /// Ensembl's own external-accession -> gene cross-references for one release
    /// (Species.Assembly.Release.uniprot.tsv.gz), used as a second opinion beside the gene links a
    /// UniProt XML carries.
    ///
    /// The two are expected to disagree, and for a reason: UniProt links an accession to every Ensembl
    /// transcript that encodes it -- readthrough genes, unnamed novel genes and identical paralogs
    /// included -- while this table assigns it where Ensembl's own mapping puts it. Neither is dropped;
    /// the resolver reports per gene whether this table agrees.
    ///
    /// The dump has one row per TRANSCRIPT, so a pair can appear many times with different evidence;
    /// the strongest is kept (DIRECT, an assertion, over SEQUENCE_MATCH and INFERRED_PAIR, inferences).
    /// </summary>
    public sealed class EnsemblXrefTable
    {
        /// <summary>The column layout this reader accepts. Anything else is refused rather than guessed at.</summary>
        public static readonly IReadOnlyList<string> Columns = new[]
        {
            "gene_stable_id", "transcript_stable_id", "protein_stable_id", "xref", "db_name",
            "info_type", "source_identity", "xref_identity", "linkage_type"
        };

        private static readonly Dictionary<string, int> InfoTypeRank = new(StringComparer.Ordinal)
        {
            ["DIRECT"] = 0,
            ["SEQUENCE_MATCH"] = 1,
            ["INFERRED_PAIR"] = 2,
        };

        private static readonly Regex ReleaseInFileName =
            new(@"\.(\d+)\.[a-z]+\.tsv(\.gz)?$", RegexOptions.Compiled | RegexOptions.IgnoreCase);

        /// <summary>accession -> gene -> strongest info_type for that pair.</summary>
        private readonly Dictionary<string, Dictionary<string, string>> _links;

        private EnsemblXrefTable(Dictionary<string, Dictionary<string, string>> links, string sourceFileName,
            string sourceSha256, string release)
        {
            _links = links;
            SourceFileName = sourceFileName;
            SourceSha256 = sourceSha256;
            Release = release;
        }

        public string SourceFileName { get; }

        /// <summary>Lower-case hex sha256 of the file's bytes as read (the compressed bytes for a .gz).</summary>
        public string SourceSha256 { get; }

        /// <summary>The Ensembl release parsed from the file name, or null when it carries none.</summary>
        public string Release { get; }

        /// <summary>The number of distinct accessions in the table.</summary>
        public int AccessionCount => _links.Count;

        public bool ContainsAccession(string accession) => accession != null && _links.ContainsKey(accession);

        /// <summary>
        /// True when the table links <paramref name="accession"/> to <paramref name="geneId"/> (both exact,
        /// stable gene id), with the strongest evidence type Ensembl gave for that pair.
        /// </summary>
        public bool TryGetLink(string accession, string geneId, out string infoType)
        {
            infoType = null;
            return accession != null && geneId != null
                && _links.TryGetValue(accession, out var genes)
                && genes.TryGetValue(geneId, out infoType);
        }

        /// <summary>
        /// Every gene the table links <paramref name="accession"/> to (exact), with the strongest evidence
        /// type for each, in ordinal order of gene id. Empty when the accession is not in the table.
        /// </summary>
        public IReadOnlyList<(string GeneId, string InfoType)> GenesFor(string accession) =>
            accession != null && _links.TryGetValue(accession, out var genes)
                ? genes.OrderBy(g => g.Key, StringComparer.Ordinal).Select(g => (g.Key, g.Value)).ToList()
                : Array.Empty<(string, string)>();

        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header is not the expected layout.</exception>
        public static EnsemblXrefTable Load(string path)
        {
            if (!File.Exists(path))
            {
                throw new FileNotFoundException("Ensembl xref table not found.", path);
            }

            string sha256;
            using (var hashStream = File.OpenRead(path))
            {
                sha256 = Convert.ToHexString(SHA256.HashData(hashStream)).ToLowerInvariant();
            }

            using var file = File.OpenRead(path);
            using Stream content = path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase)
                ? new GZipStream(file, CompressionMode.Decompress)
                : file;
            using var reader = new StreamReader(content);

            string header = reader.ReadLine() ?? "";
            if (!header.Split('\t').SequenceEqual(Columns))
            {
                throw new InvalidDataException(
                    $"{Path.GetFileName(path)}: unexpected columns. Got [{header.Replace('\t', ',')}], " +
                    $"expected [{string.Join(",", Columns)}].");
            }

            var links = new Dictionary<string, Dictionary<string, string>>(StringComparer.Ordinal);
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                string[] cells = line.Split('\t');
                if (cells.Length != Columns.Count)
                {
                    continue;
                }

                string gene = cells[0], accession = cells[3], info = cells[5];
                if (!links.TryGetValue(accession, out var genes))
                {
                    genes = new Dictionary<string, string>(StringComparer.Ordinal);
                    links[accession] = genes;
                }

                if (!genes.TryGetValue(gene, out string existing) || Rank(info) < Rank(existing))
                {
                    genes[gene] = info;
                }
            }

            var release = ReleaseInFileName.Match(Path.GetFileName(path));
            return new EnsemblXrefTable(links, Path.GetFileName(path), sha256, release.Success ? release.Groups[1].Value : null);
        }

        private static int Rank(string infoType) => InfoTypeRank.TryGetValue(infoType ?? "", out int r) ? r : 99;
    }
}
