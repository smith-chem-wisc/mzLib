using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace UsefulProteomicsDatabases.GeneOntology
{
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
                if (!string.Equals(row.GoRelease, ontology.Release, StringComparison.Ordinal)
                    || !string.Equals(row.GoOboSha256, ontology.SourceSha256, StringComparison.Ordinal))
                {
                    throw new ArgumentException(
                        $"Row for group '{row.ProteinGroup}' comes from Gene Ontology release '{row.GoRelease}', " +
                        $"but the category map is applied against '{ontology.Release}'.");
                }
                if (row.GoId != null)
                {
                    terms.Add(row.GoId);
                }
            }

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
    }
}
