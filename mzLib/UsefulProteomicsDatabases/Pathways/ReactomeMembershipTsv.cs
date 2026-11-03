using System;
using System.Collections.Generic;
using System.IO;
using UsefulProteomicsDatabases.Annotation;

namespace UsefulProteomicsDatabases.Pathways
{
    /// <summary>
    /// Writes and reads one database's Reactome membership (<see cref="ReactomeAccessionAnnotation"/>).
    /// <code>
    /// #!reactome_membership_format 1
    /// #!mzlib_version 1.0.597+...
    /// #!mzlib_release 1.0.597
    /// #!definition kegg:DEF-REACTOME-MEMBERSHIP v1
    /// #!annotation_db_sha256 ...
    /// #!species_codes NCBITaxon:10090=MMU;NCBITaxon:9606=HSA
    /// #!unvalidated_organisms NCBITaxon:9031=3
    /// #!dropped_foreign_references 8
    /// accession	organism	pathway_id	molecule	sequence_accession	inherited	status
    /// P04406	NCBITaxon:9606	R-HSA-70171			false	annotated
    /// P04406-2	NCBITaxon:9606	R-HSA-70171			true	annotated
    /// Q9XXXX	NCBITaxon:9606					no_reference
    /// </code>
    /// species_codes, unvalidated_organisms and dropped_foreign_references are written only when not empty or zero.
    /// </summary>
    public static class ReactomeMembershipTsv
    {
        /// <summary>The format version written on the first line.</summary>
        public const int FormatVersion = 1;

        /// <summary>The column names, in order.</summary>
        public static readonly IReadOnlyList<string> Columns =
            new[] { "accession", "organism", "pathway_id", "molecule", "sequence_accession", "inherited", "status" };

        /// <summary>"annotated", "no_reference" or "species_not_covered".</summary>
        public static string StatusName(ReactomeMembershipStatus status) => throw new NotImplementedException();

        /// <summary>Writes every row of <paramref name="annotation"/>, with the header.</summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">A value contains a tab or a line break.</exception>
        public static void Write(TextWriter output, ReactomeAccessionAnnotation annotation) => throw new NotImplementedException();

        /// <summary>Reads a database's Reactome membership file.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's; a status is unknown; a row
        /// that carries a pathway is not annotated, or an annotated row carries none; an entry has a pathway-less row
        /// and another row; a sequence_accession disagrees with its rule; an inherited row is isoform-scoped; or a row
        /// repeats.</exception>
        public static AnnotationTable<ReactomeMembershipRow> Read(string path) => throw new NotImplementedException();
    }
}
