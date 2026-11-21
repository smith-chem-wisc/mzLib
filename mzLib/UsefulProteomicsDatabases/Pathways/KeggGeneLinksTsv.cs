using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Proteomics;
using UsefulProteomicsDatabases.Annotation;

namespace UsefulProteomicsDatabases.Pathways
{
    /// <summary>One KEGG GENES id UniProt lists on an entry, verbatim (e.g. "hsa:2597").</summary>
    /// <param name="Accession">The database entry's accession, verbatim.</param>
    /// <param name="KeggGeneId">The KEGG GENES id, verbatim.</param>
    public sealed record KeggGeneLinkRow(string Accession, string KeggGeneId);

    /// <summary>
    /// Writes and reads one database's KEGG gene links: the KEGG GENES ids UniProt distributes, for link-outs only. No
    /// KEGG pathway, orthology or module is stored, and nothing here calls KEGG. There are no status rows: for a
    /// database that has the table, an entry with no row is one UniProt lists no KEGG id for.
    /// <code>
    /// #!kegg_gene_links_format 1
    /// #!mzlib_version 1.0.597+...
    /// #!mzlib_release 1.0.597
    /// #!annotation_db_sha256 ...
    /// accession	kegg_gene_id
    /// P04406	hsa:2597
    /// </code>
    /// </summary>
    public static class KeggGeneLinksTsv
    {
        /// <summary>The format version written on the first line.</summary>
        public const int FormatVersion = 1;

        private const string FormatKey = "kegg_gene_links_format";

        /// <summary>The column names, in order.</summary>
        public static readonly IReadOnlyList<string> Columns = new[] { "accession", "kegg_gene_id" };

        /// <summary>
        /// Every (entry, KEGG GENES id) pair of the database, in ordinal order. Decoys and nulls are skipped, and an
        /// accession present more than once has its ids unioned.
        /// </summary>
        /// <exception cref="ArgumentNullException">annotationProteins is null.</exception>
        public static IReadOnlyList<KeggGeneLinkRow> Links(IEnumerable<Protein> annotationProteins)
        {
            ArgumentNullException.ThrowIfNull(annotationProteins);
            return annotationProteins
                .Where(p => p != null && !p.IsDecoy)
                .SelectMany(p => p.KeggGeneIds.Select(id => new KeggGeneLinkRow(p.Accession, id)))
                .Distinct()
                .OrderBy(r => r.Accession, StringComparer.Ordinal)
                .ThenBy(r => r.KeggGeneId, StringComparer.Ordinal)
                .ToList();
        }

        /// <summary>Writes the rows in the order given, with the header. Every row is checked before anything is written.</summary>
        /// <exception cref="ArgumentNullException">output or rows is null.</exception>
        /// <exception cref="ArgumentException">annotationDbSha256 is null or empty; a value is empty or contains a tab or
        /// a line break; or a row repeats.</exception>
        public static void Write(TextWriter output, IEnumerable<KeggGeneLinkRow> rows, string annotationDbSha256)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(rows);
            if (string.IsNullOrEmpty(annotationDbSha256))
            {
                throw new ArgumentException("The annotation database's sha256 is required: the table is per database.",
                    nameof(annotationDbSha256));
            }
            var checkedRows = rows.ToList();
            var seen = new HashSet<KeggGeneLinkRow>();
            foreach (var row in checkedRows)
            {
                if (string.IsNullOrEmpty(row.Accession) || string.IsNullOrEmpty(row.KeggGeneId))
                {
                    throw new ArgumentException("A KEGG gene link has an empty accession or id.", nameof(rows));
                }
                TsvHeader.RejectSeparators(row.Accession, "accession");
                TsvHeader.RejectSeparators(row.KeggGeneId, "kegg_gene_id");
                if (!seen.Add(row))
                {
                    throw new ArgumentException($"'{row.Accession}', {row.KeggGeneId} appears twice.", nameof(rows));
                }
            }

            TsvHeader.Write(output, FormatKey, FormatVersion.ToString(CultureInfo.InvariantCulture));
            TsvHeader.WriteProducer(output);
            TsvHeader.Write(output, "annotation_db_sha256", annotationDbSha256);
            output.Write(string.Join('\t', Columns) + "\n");
            foreach (var row in checkedRows)
            {
                output.Write($"{row.Accession}\t{row.KeggGeneId}" + "\n");
            }
        }

        /// <summary>Reads a database's KEGG gene links file.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's; a value is empty; or a
        /// row repeats.</exception>
        public static AnnotationTable<KeggGeneLinkRow> Read(string path)
        {
            var seen = new HashSet<KeggGeneLinkRow>();
            return AnnotationTableReader.Read(path, FormatKey, FormatVersion, new[] { "annotation_db_sha256" },
                Array.Empty<string>(), Columns,
                (cells, fileName, lineNumber) =>
                {
                    string where = $"{fileName} line {lineNumber}";
                    if (cells[0].Length == 0)
                    {
                        throw new InvalidDataException($"{where}: accession is empty.");
                    }
                    if (cells[1].Length == 0)
                    {
                        throw new InvalidDataException($"{where}: kegg_gene_id is empty.");
                    }
                    var row = new KeggGeneLinkRow(cells[0], cells[1]);
                    if (!seen.Add(row))
                    {
                        throw new InvalidDataException($"{where}: '{row.Accession}', {row.KeggGeneId} appears twice.");
                    }
                    return row;
                });
        }
    }
}
