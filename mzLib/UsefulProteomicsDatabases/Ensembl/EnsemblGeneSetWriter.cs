using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using Omics.BioPolymerGroup;

namespace UsefulProteomicsDatabases.Ensembl
{
    /// <summary>
    /// Writes an <see cref="EnsemblGeneSet"/> as a compact gene table: the provenance of the GTF it was
    /// read from, as "#!key value" lines in the GTF's own header style, then one tab-separated row per
    /// gene. <see cref="EnsemblGeneSetReader"/> reads it back into an identical set.
    ///
    /// The provenance written is the GTF's, not this file's, so a resolution counted against the table is
    /// keyed exactly as one counted against the GTF. The output is deterministic -- genes in ordinal
    /// order, "\n" line endings on every platform -- so the same set always gives the same rows.
    /// </summary>
    public static class EnsemblGeneSetWriter
    {
        /// <summary>The format version written on the first line. The reader refuses any other.</summary>
        public const int FormatVersion = 1;

        internal const string FormatKey = "ensembl-gene-set-format";
        internal const string SourceFileKey = "source-file";
        internal const string SourceSha256Key = "source-sha256";
        internal const string ReleaseKey = "release";
        internal const string GenomeBuildKey = "genome-build";
        internal const string GenebuildLastUpdatedKey = "genebuild-last-updated";

        /// <summary>The gene row layout, in column order.</summary>
        public static readonly IReadOnlyList<TsvColumn<EnsemblGene>> Schema = new[]
        {
            new TsvColumn<EnsemblGene>("gene_id", g => g.GeneId),
            new TsvColumn<EnsemblGene>("gene_version", g => g.Version?.ToString(CultureInfo.InvariantCulture)),
            new TsvColumn<EnsemblGene>("gene_biotype", g => g.Biotype),
            new TsvColumn<EnsemblGene>("gene_name", g => g.Symbol),
            new TsvColumn<EnsemblGene>("seq_region", g => g.SeqRegion),
        };

        /// <summary>
        /// Writes <paramref name="geneSet"/> to <paramref name="outputPath"/>, gzip-compressed when the path
        /// ends in ".gz". A null gene version or symbol is written as an empty cell; a null release, genome
        /// build or genebuild date omits that header line.
        /// </summary>
        /// <exception cref="ArgumentException">A value contains a tab or line break, which would shift the
        /// row. A GTF read by <see cref="EnsemblGeneSet.LoadGtf"/> cannot produce one.</exception>
        public static void Write(string outputPath, EnsemblGeneSet geneSet)
        {
            ArgumentNullException.ThrowIfNull(outputPath);
            ArgumentNullException.ThrowIfNull(geneSet);

            using var file = File.Create(outputPath);
            using Stream content = outputPath.EndsWith(".gz", StringComparison.OrdinalIgnoreCase)
                ? new GZipStream(file, CompressionLevel.SmallestSize)
                : file;
            using var writer = new StreamWriter(content) { NewLine = "\n" };

            WriteHeaderLine(writer, FormatKey, FormatVersion.ToString(CultureInfo.InvariantCulture));
            WriteHeaderLine(writer, SourceFileKey, geneSet.SourceFileName);
            WriteHeaderLine(writer, SourceSha256Key, geneSet.SourceSha256);
            WriteHeaderLine(writer, ReleaseKey, geneSet.Release);
            WriteHeaderLine(writer, GenomeBuildKey, geneSet.GenomeBuild);
            WriteHeaderLine(writer, GenebuildLastUpdatedKey, geneSet.GenebuildLastUpdated);

            writer.WriteLine(TsvWriter.HeaderLine(Schema));
            foreach (var gene in geneSet.Genes)
            {
                foreach (var column in Schema)
                {
                    RejectSeparators(column.GetValue(gene), column.Header, gene.GeneId);
                }
                writer.WriteLine(TsvWriter.RowLine(Schema, gene));
            }
        }

        private static void WriteHeaderLine(TextWriter writer, string key, string value)
        {
            if (value == null)
            {
                return;
            }
            RejectSeparators(value, key, null);
            writer.WriteLine($"#!{key} {value}");
        }

        private static void RejectSeparators(string value, string field, string geneId)
        {
            if (value != null && value.IndexOfAny(new[] { '\t', '\r', '\n' }) >= 0)
            {
                string where = geneId == null ? field : $"{field} of {geneId}";
                throw new ArgumentException($"{where} contains a tab or line break and cannot be written as one cell.");
            }
        }
    }
}
