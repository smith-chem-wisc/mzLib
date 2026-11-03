using System;
using System.Collections.Generic;
using System.IO;
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

        /// <summary>The column names, in order.</summary>
        public static readonly IReadOnlyList<string> Columns = new[] { "accession", "kegg_gene_id" };

        /// <summary>
        /// Every (entry, KEGG GENES id) pair of the database, in ordinal order. Decoys and nulls are skipped, and an
        /// accession present more than once has its ids unioned.
        /// </summary>
        /// <exception cref="ArgumentNullException">annotationProteins is null.</exception>
        public static IReadOnlyList<KeggGeneLinkRow> Links(IEnumerable<Protein> annotationProteins) => throw new NotImplementedException();

        /// <summary>Writes the rows in the order given, with the header. Every row is checked before anything is written.</summary>
        /// <exception cref="ArgumentNullException">output or rows is null.</exception>
        /// <exception cref="ArgumentException">annotationDbSha256 is null or empty; a value is empty or contains a tab or
        /// a line break; or a row repeats.</exception>
        public static void Write(TextWriter output, IEnumerable<KeggGeneLinkRow> rows, string annotationDbSha256) =>
            throw new NotImplementedException();

        /// <summary>Reads a database's KEGG gene links file.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's; a value is empty; or a
        /// row repeats.</exception>
        public static AnnotationTable<KeggGeneLinkRow> Read(string path) => throw new NotImplementedException();
    }
}
