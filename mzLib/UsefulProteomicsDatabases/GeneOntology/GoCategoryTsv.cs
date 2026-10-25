using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases.Annotation;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>One row of a category table: a term and one category (and subcategory) of a map it falls under.</summary>
    /// <param name="Subcategory">"category:subcategory", or null when the term reaches the category only through a
    /// category-level anchor.</param>
    public sealed record GoCategoryRow(string GoId, string Category, string Subcategory);

    /// <summary>
    /// Writes one consumer map's term-to-category table for a set of GO annotation rows:
    /// <code>
    /// #!go_category_format 1
    /// #!mzlib_version 1.0.593+...
    /// #!mzlib_release 1.0.593
    /// #!go_release releases/2026-07-26
    /// #!go_obo_sha256 ...
    /// #!category_map organelle 1 ...
    /// go_id	category	subcategory
    /// </code>
    /// One row per (term, category, subcategory), for every term in the annotation rows that is at or below
    /// one of the map's anchors; a term under no anchor gets no row, so absence means "outside the map".
    ///
    /// Categories are deliberately not columns of the annotation rows: a run may apply any number of consumer
    /// maps, and a changed map rewrites only this small table. The two files join on (go_id, go_release),
    /// so the writer refuses rows whose release differs from the resolver's -- a table from one release
    /// joined to rows from another would silently mis-assign.
    /// </summary>
    public static class GoCategoryTsv
    {
        /// <summary>The format version written on the first line.</summary>
        public const int FormatVersion = 1;

        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">A row's go_release or go.obo sha256 differs from the resolver's
        /// ontology.</exception>
        public static void Write(TextWriter output, GoCategoryResolver resolver, IEnumerable<GoAnnotationRow> rows)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(resolver);
            ArgumentNullException.ThrowIfNull(rows);

            var ontology = resolver.Ontology;
            var terms = new SortedSet<string>(StringComparer.Ordinal);
            foreach (var row in rows)
            {
                RequireSame(row.GoRelease, ontology.Release, "go_release", row);
                RequireSame(row.GoOboSha256, ontology.SourceSha256, "go_obo_sha256", row);
                if (row.GoId != null)
                {
                    terms.Add(row.GoId);
                }
            }

            WriteTable(output, resolver, terms);
        }

        /// <summary>
        /// Writes the map's table for every term of the resolver's ontology release, not only the terms one
        /// search's rows hold. The table then depends on the release and the map alone, so a store keeps one per
        /// (go.obo, map) and joins it to any database's annotation.
        /// </summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        public static void Write(TextWriter output, GoCategoryResolver resolver)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(resolver);
            WriteTable(output, resolver, resolver.Ontology.TermIds);
        }

        /// <summary>The column names, in order.</summary>
        public static readonly IReadOnlyList<string> Columns = new[] { "go_id", "category", "subcategory" };

        /// <summary>Reads a category table.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's; go_id or category
        /// is empty; a subcategory is not qualified by its row's category; or a row repeats.</exception>
        public static AnnotationTable<GoCategoryRow> Read(string path)
        {
            var seen = new HashSet<GoCategoryRow>();
            return AnnotationTableReader.Read(path, "go_category_format", FormatVersion,
                new[] { "go_obo_sha256", "category_map" }, new[] { "go_release" }, Columns,
                (cells, fileName, lineNumber) =>
                {
                    string where = $"{fileName} line {lineNumber}";
                    if (cells[0].Length == 0 || cells[1].Length == 0)
                    {
                        throw new InvalidDataException($"{where}: go_id or category is empty.");
                    }
                    var row = new GoCategoryRow(cells[0], cells[1], AnnotationTableReader.Cell(cells[2]));
                    if (row.Subcategory != null && !row.Subcategory.StartsWith(row.Category + ":", StringComparison.Ordinal))
                    {
                        throw new InvalidDataException($"{where}: subcategory '{row.Subcategory}' is not qualified by category '{row.Category}'.");
                    }
                    if (!seen.Add(row))
                    {
                        throw new InvalidDataException($"{where}: repeats an earlier row.");
                    }
                    return row;
                });
        }

        private static void WriteTable(TextWriter output, GoCategoryResolver resolver, IEnumerable<string> terms)
        {
            var ontology = resolver.Ontology;
            var map = resolver.Map;
            TsvHeader.Write(output, "go_category_format", FormatVersion.ToString(CultureInfo.InvariantCulture));
            TsvHeader.WriteProducer(output);
            TsvHeader.Write(output, "go_release", ontology.Release);
            TsvHeader.Write(output, "go_obo_sha256", ontology.SourceSha256);
            TsvHeader.Write(output, "category_map", $"{map.MapName} {map.MapVersion} {map.SourceSha256}");
            output.Write("go_id\tcategory\tsubcategory" + "\n");
            foreach (string goId in terms)
            {
                foreach (var category in resolver.Categorize(goId))
                {
                    output.Write($"{goId}\t{category.Category}\t{category.Subcategory}" + "\n");
                }
            }
        }

        private static void RequireSame(string actual, string expected, string field, GoAnnotationRow row)
        {
            if (!string.Equals(actual, expected, StringComparison.Ordinal))
            {
                throw new ArgumentException(
                    $"Row for group '{row.ProteinGroup}' has {field} '{actual}', but the category map is applied " +
                    $"against {field} '{expected}'.");
            }
        }
    }
}
