#nullable enable
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Text;
using UsefulProteomicsDatabases.Annotation;

namespace UsefulProteomicsDatabases.Taxonomy
{
    /// <summary>
    /// Reads and writes the two tables a taxonomy is stored as: the taxa of a <see cref="TaxonomyTree"/> and the
    /// proteomes of a <see cref="ProteomeTaxonMap"/>.
    /// </summary>
    /// <remarks>
    /// <para>Both are tab-separated text in the layout mzLib's other annotation tables use: "#!key value" header lines,
    /// then a column header, then one row per record, with "\n" line endings and rows in a fixed order so the same
    /// taxonomy always writes the same bytes.</para>
    /// <code>
    /// #!taxonomy_nodes_format 1
    /// #!mzlib_version 1.0.600+abc1234
    /// #!mzlib_release 1.0.600
    /// #!source https://rest.uniprot.org/proteomes/search
    /// #!source_release 2026_04
    /// #!source_query reference:true
    /// taxon_id	parent_taxon_id	rank	scientific_name	common_name	hidden
    /// 2		domain	Bacteria	bacteria	false
    /// 1224	3379134	phylum	Pseudomonadota	proteobacteria	false
    /// </code>
    /// <para>A file a user writes by hand to add taxa needs only the format line, a source line and the columns. An
    /// empty parent marks a root; an empty rank or common name means the source stated none.</para>
    /// <para>Reading refuses rather than guesses: a wrong header, a short row, an identifier that is not a whole
    /// number or that repeats are all <see cref="InvalidDataException"/>s naming the file and line.</para>
    /// </remarks>
    public static class TaxonomySnapshotTsv
    {
        /// <summary>The version both table formats are at.</summary>
        public const int FormatVersion = 1;

        private const string NodesFormatKey = "taxonomy_nodes_format";
        private const string ProteomesFormatKey = "proteome_taxa_format";
        private const string SourceKey = "source";
        private const string ReleaseKey = "source_release";
        private const string QueryKey = "source_query";

        private static readonly string[] RequiredKeys = { SourceKey };
        private static readonly string[] OptionalKeys = { ReleaseKey, QueryKey };

        /// <summary>The node table's columns, in order.</summary>
        public static readonly IReadOnlyList<string> NodeColumns = new[]
        {
            "taxon_id", "parent_taxon_id", "rank", "scientific_name", "common_name", "hidden"
        };

        /// <summary>The proteome table's columns, in order.</summary>
        public static readonly IReadOnlyList<string> ProteomeColumns = new[]
        {
            "proteome_id", "taxon_id", "organism", "protein_count", "proteome_type"
        };

        /// <summary>Writes a tree's taxa, in ascending order of identifier.</summary>
        /// <exception cref="ArgumentException">A name or rank contains a tab or line break.</exception>
        public static void WriteNodes(TextWriter output, TaxonomyTree tree)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(tree);
            WriteHeader(output, NodesFormatKey, NodeColumns, tree.Provenance);
            foreach (TaxonomyNode node in tree.Nodes)
            {
                TsvHeader.RejectSeparators(node.Rank, "rank");
                TsvHeader.RejectSeparators(node.ScientificName, "scientific_name");
                TsvHeader.RejectSeparators(node.CommonName, "common_name");
                output.Write(string.Join('\t',
                    Number(node.TaxonId),
                    node.ParentTaxonId is long parent ? Number(parent) : "",
                    node.Rank,
                    node.ScientificName,
                    node.CommonName ?? "",
                    node.IsHidden ? "true" : "false") + "\n");
            }
        }

        /// <summary>
        /// Writes a tree's taxa to a file; a path ending in ".gz" is written as gzip. An existing file is replaced only
        /// once the whole table is written, so a failed write leaves it as it was.
        /// </summary>
        public static void WriteNodes(string path, TaxonomyTree tree) =>
            WriteFile(path, output => WriteNodes(output, tree));

        /// <summary>Writes a map's proteomes, in ordinal order of identifier.</summary>
        /// <exception cref="ArgumentException">A proteome identifier, organism or type contains a tab or line break.</exception>
        public static void WriteProteomes(TextWriter output, ProteomeTaxonMap proteomes)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(proteomes);
            WriteHeader(output, ProteomesFormatKey, ProteomeColumns, proteomes.Provenance);
            foreach (ProteomeTaxon proteome in proteomes.Proteomes)
            {
                TsvHeader.RejectSeparators(proteome.ProteomeId, "proteome_id");
                TsvHeader.RejectSeparators(proteome.Organism, "organism");
                TsvHeader.RejectSeparators(proteome.ProteomeType, "proteome_type");
                output.Write(string.Join('\t',
                    proteome.ProteomeId,
                    Number(proteome.TaxonId),
                    proteome.Organism,
                    proteome.ProteinCount is long count ? Number(count) : "",
                    proteome.ProteomeType ?? "") + "\n");
            }
        }

        /// <summary>
        /// Writes a map's proteomes to a file; a path ending in ".gz" is written as gzip. An existing file is replaced
        /// only once the whole table is written, so a failed write leaves it as it was.
        /// </summary>
        public static void WriteProteomes(string path, ProteomeTaxonMap proteomes) =>
            WriteFile(path, output => WriteProteomes(output, proteomes));

        /// <summary>Reads a node table and builds its tree.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The file is not a node table, or its taxa do not form a tree.</exception>
        public static TaxonomyTree ReadNodes(string path)
        {
            AnnotationTable<TaxonomyNode> table = ReadNodeRows(path);
            return BuildTree(table);
        }

        /// <summary>
        /// Reads a node table's rows without requiring them to form a tree on their own, which an overlay file does
        /// not: its taxa name parents in the tree they are added to.
        /// </summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The file is not a node table.</exception>
        public static AnnotationTable<TaxonomyNode> ReadNodeRows(string path)
        {
            ArgumentNullException.ThrowIfNull(path);
            if (!File.Exists(path))
            {
                throw new FileNotFoundException("Taxonomy node table not found.", path);
            }
            using TextReader reader = OpenText(path);
            return ReadNodeRows(reader, Path.GetFileName(path));
        }

        /// <summary>Reads a proteome table.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The file is not a proteome table.</exception>
        public static ProteomeTaxonMap ReadProteomes(string path)
        {
            ArgumentNullException.ThrowIfNull(path);
            if (!File.Exists(path))
            {
                throw new FileNotFoundException("Proteome taxon table not found.", path);
            }
            using TextReader reader = OpenText(path);
            return ReadProteomes(reader, Path.GetFileName(path));
        }

        internal static TaxonomyTree ReadNodesResource(string resourceName)
        {
            using TextReader reader = OpenResource(resourceName);
            return BuildTree(ReadNodeRows(reader, resourceName));
        }

        internal static ProteomeTaxonMap ReadProteomesResource(string resourceName)
        {
            using TextReader reader = OpenResource(resourceName);
            return ReadProteomes(reader, resourceName);
        }

        private static TaxonomyTree BuildTree(AnnotationTable<TaxonomyNode> table)
        {
            try
            {
                return TaxonomyTree.FromNodes(table.Rows, ProvenanceOf(table.Header));
            }
            catch (InvalidDataException e)
            {
                throw new InvalidDataException($"{table.SourceFileName}: {e.Message}", e);
            }
        }

        private static AnnotationTable<TaxonomyNode> ReadNodeRows(TextReader reader, string fileName)
        {
            var seen = new HashSet<long>();
            return AnnotationTableReader.Read(reader, fileName, NodesFormatKey, FormatVersion, RequiredKeys, OptionalKeys, NodeColumns,
                (cells, file, lineNumber) =>
                {
                    long taxonId = TaxonId(cells[0], "taxon_id", file, lineNumber);
                    if (!seen.Add(taxonId))
                    {
                        throw new InvalidDataException($"{file} line {lineNumber}: taxon_id {cells[0]} appears twice.");
                    }
                    long? parent = cells[1].Length == 0 ? null : TaxonId(cells[1], "parent_taxon_id", file, lineNumber);
                    if (parent == taxonId)
                    {
                        throw new InvalidDataException($"{file} line {lineNumber}: taxon {cells[0]} names itself as its parent.");
                    }
                    if (cells[3].Length == 0)
                    {
                        throw new InvalidDataException($"{file} line {lineNumber}: scientific_name is empty.");
                    }
                    bool hidden = AnnotationTableReader.Bool(cells[5], "hidden", file, lineNumber)
                        ?? throw new InvalidDataException($"{file} line {lineNumber}: hidden is empty; expected true or false.");
                    return new TaxonomyNode(taxonId, parent, cells[2], cells[3], AnnotationTableReader.Cell(cells[4]), hidden);
                });
        }

        private static ProteomeTaxonMap ReadProteomes(TextReader reader, string fileName)
        {
            var seen = new HashSet<string>(StringComparer.Ordinal);
            AnnotationTable<ProteomeTaxon> table = AnnotationTableReader.Read(reader, fileName, ProteomesFormatKey, FormatVersion,
                RequiredKeys, OptionalKeys, ProteomeColumns,
                (cells, file, lineNumber) =>
                {
                    if (cells[0].Length == 0 || !seen.Add(cells[0]))
                    {
                        throw new InvalidDataException($"{file} line {lineNumber}: proteome_id '{cells[0]}' is empty or repeats.");
                    }
                    long taxonId = TaxonId(cells[1], "taxon_id", file, lineNumber);
                    long? proteinCount = null;
                    if (cells[3].Length != 0)
                    {
                        if (!long.TryParse(cells[3], NumberStyles.None, CultureInfo.InvariantCulture, out long count))
                        {
                            throw new InvalidDataException($"{file} line {lineNumber}: protein_count '{cells[3]}' is not a whole number.");
                        }
                        proteinCount = count;
                    }
                    return new ProteomeTaxon(cells[0], taxonId, cells[2], proteinCount, AnnotationTableReader.Cell(cells[4]));
                });
            return ProteomeTaxonMap.FromProteomes(table.Rows, ProvenanceOf(table.Header));
        }

        private static long TaxonId(string cell, string column, string fileName, int lineNumber)
        {
            if (!long.TryParse(cell, NumberStyles.None, CultureInfo.InvariantCulture, out long taxonId) || taxonId <= 0)
            {
                throw new InvalidDataException($"{fileName} line {lineNumber}: {column} '{cell}' is not a positive whole number.");
            }
            return taxonId;
        }

        private static TaxonomyProvenance ProvenanceOf(IReadOnlyDictionary<string, string> header) =>
            new(header[SourceKey],
                header.TryGetValue(ReleaseKey, out var release) && release.Length != 0 ? release : null,
                header.TryGetValue(QueryKey, out var query) && query.Length != 0 ? query : null);

        private static void WriteHeader(TextWriter output, string formatKey, IReadOnlyList<string> columns, TaxonomyProvenance provenance)
        {
            TsvHeader.Write(output, formatKey, FormatVersion.ToString(CultureInfo.InvariantCulture));
            TsvHeader.WriteProducer(output);
            TsvHeader.Write(output, SourceKey, provenance.Source);
            TsvHeader.Write(output, ReleaseKey, provenance.Release);
            TsvHeader.Write(output, QueryKey, provenance.Query);
            output.Write(string.Join('\t', columns) + "\n");
        }

        private static string Number(long value) => value.ToString(CultureInfo.InvariantCulture);

        private static TextReader OpenText(string path)
        {
            Stream file = File.OpenRead(path);
            Stream content = path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase)
                ? new GZipStream(file, CompressionMode.Decompress)
                : file;
            return new StreamReader(content);
        }

        /// <summary>Writes beside <paramref name="path"/> first, then moves the finished file into place.</summary>
        private static void WriteFile(string path, Action<TextWriter> write)
        {
            ArgumentNullException.ThrowIfNull(path);
            string partial = path + ".partial";
            try
            {
                using (TextWriter output = CreateText(partial, gzip: path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase)))
                {
                    write(output);
                }
                File.Move(partial, path, overwrite: true);
            }
            finally
            {
                if (File.Exists(partial))
                {
                    File.Delete(partial);
                }
            }
        }

        private static TextWriter CreateText(string path, bool gzip)
        {
            Stream file = File.Create(path);
            Stream content = gzip
                ? new GZipStream(file, CompressionLevel.SmallestSize)
                : file;
            return new StreamWriter(content, new UTF8Encoding(encoderShouldEmitUTF8Identifier: false)) { NewLine = "\n" };
        }

        private static TextReader OpenResource(string resourceName)
        {
            Stream? resource = typeof(TaxonomySnapshotTsv).Assembly.GetManifestResourceStream(resourceName)
                ?? throw new MzLibException($"The embedded taxonomy resource '{resourceName}' is missing from the assembly.");
            return new StreamReader(new GZipStream(resource, CompressionMode.Decompress));
        }
    }
}
