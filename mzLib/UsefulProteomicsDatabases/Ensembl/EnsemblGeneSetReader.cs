using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.IO.Compression;

namespace UsefulProteomicsDatabases.Ensembl
{
    /// <summary>
    /// Reads the compact gene table <see cref="EnsemblGeneSetWriter"/> writes back into an
    /// <see cref="EnsemblGeneSet"/>. The set carries the provenance of the GTF the table was made from --
    /// its file name, sha256, release and genome build -- not the table's own, so a resolution counted
    /// against it is keyed exactly as one counted against the GTF.
    ///
    /// The reader is strict: an unknown format version, an unknown header key, a missing source file or
    /// sha256, a different column layout, a short row or a repeated gene id is refused rather than guessed
    /// at. A table that cannot say which GTF it came from cannot key a resolution.
    /// </summary>
    public static class EnsemblGeneSetReader
    {
        /// <summary>Reads a gene table (plain, or gzip-compressed when the path ends in ".gz").</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The file is not a gene table this reader understands.</exception>
        public static EnsemblGeneSet Load(string filePath)
        {
            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException("Ensembl gene table not found.", filePath);
            }

            string name = Path.GetFileName(filePath);
            using var file = File.OpenRead(filePath);
            using Stream content = filePath.EndsWith(".gz", StringComparison.OrdinalIgnoreCase)
                ? new GZipStream(file, CompressionMode.Decompress)
                : file;
            using var reader = new StreamReader(content);

            var header = new Dictionary<string, string>(StringComparer.Ordinal);
            int lineNumber = 0;
            string line;
            while ((line = reader.ReadLine()) != null && line.StartsWith("#!", StringComparison.Ordinal))
            {
                lineNumber++;
                int space = line.IndexOf(' ');
                if (space < 0)
                {
                    throw new InvalidDataException($"{name} line {lineNumber}: header line has no value.");
                }
                string key = line.Substring(2, space - 2);
                if (!IsKnownKey(key))
                {
                    throw new InvalidDataException($"{name} line {lineNumber}: unknown header key '{key}'.");
                }
                if (!header.TryAdd(key, line.Substring(space + 1)))
                {
                    throw new InvalidDataException($"{name} line {lineNumber}: header key '{key}' is repeated.");
                }
                if (lineNumber == 1 && key != EnsemblGeneSetWriter.FormatKey)
                {
                    throw new InvalidDataException($"{name}: the first line must be #!{EnsemblGeneSetWriter.FormatKey}.");
                }
            }

            string expectedFormat = EnsemblGeneSetWriter.FormatVersion.ToString(CultureInfo.InvariantCulture);
            if (!header.TryGetValue(EnsemblGeneSetWriter.FormatKey, out var format) || format != expectedFormat)
            {
                throw new InvalidDataException(
                    $"{name}: not an Ensembl gene table of format {expectedFormat} (found '{format ?? "none"}').");
            }
            string sourceFile = Required(header, EnsemblGeneSetWriter.SourceFileKey, name);
            string sourceSha256 = Required(header, EnsemblGeneSetWriter.SourceSha256Key, name);

            lineNumber++;
            string expectedColumns = string.Join('\t', ColumnHeaders());
            if (line != expectedColumns)
            {
                throw new InvalidDataException(
                    $"{name} line {lineNumber}: expected the columns '{expectedColumns.Replace('\t', ' ')}'.");
            }

            var genes = new Dictionary<string, EnsemblGene>(StringComparer.Ordinal);
            while ((line = reader.ReadLine()) != null)
            {
                lineNumber++;
                string[] cells = line.Split('\t');
                if (cells.Length != EnsemblGeneSetWriter.Schema.Count)
                {
                    throw new InvalidDataException(
                        $"{name} line {lineNumber}: {cells.Length} cells, expected {EnsemblGeneSetWriter.Schema.Count}.");
                }
                if (cells[0].Length == 0 || cells[2].Length == 0 || cells[4].Length == 0)
                {
                    throw new InvalidDataException($"{name} line {lineNumber}: gene_id, gene_biotype and seq_region are required.");
                }

                int? version = null;
                if (cells[1].Length > 0)
                {
                    if (!int.TryParse(cells[1], NumberStyles.None, CultureInfo.InvariantCulture, out int v))
                    {
                        throw new InvalidDataException($"{name} line {lineNumber}: gene_version '{cells[1]}' is not a number.");
                    }
                    version = v;
                }

                var gene = new EnsemblGene(cells[0], version, cells[2], cells[3].Length == 0 ? null : cells[3], cells[4]);
                if (!genes.TryAdd(gene.GeneId, gene))
                {
                    throw new InvalidDataException($"{name} line {lineNumber}: {gene.GeneId} is repeated.");
                }
            }

            header.TryGetValue(EnsemblGeneSetWriter.ReleaseKey, out var release);
            header.TryGetValue(EnsemblGeneSetWriter.GenomeBuildKey, out var genomeBuild);
            header.TryGetValue(EnsemblGeneSetWriter.GenebuildLastUpdatedKey, out var lastUpdated);
            return EnsemblGeneSet.FromGenes(genes, sourceFile, sourceSha256, release, genomeBuild, lastUpdated);
        }

        private static bool IsKnownKey(string key) => key is EnsemblGeneSetWriter.FormatKey
            or EnsemblGeneSetWriter.SourceFileKey or EnsemblGeneSetWriter.SourceSha256Key
            or EnsemblGeneSetWriter.ReleaseKey or EnsemblGeneSetWriter.GenomeBuildKey
            or EnsemblGeneSetWriter.GenebuildLastUpdatedKey;

        private static string Required(Dictionary<string, string> header, string key, string name) =>
            header.TryGetValue(key, out var value) && value.Length > 0
                ? value
                : throw new InvalidDataException($"{name}: header #!{key} is required.");

        private static IEnumerable<string> ColumnHeaders()
        {
            foreach (var column in EnsemblGeneSetWriter.Schema)
            {
                yield return column.Header;
            }
        }
    }
}
