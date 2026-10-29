using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases.GeneOntology;

namespace UsefulProteomicsDatabases.Annotation
{
    /// <summary>
    /// Writes and reads a search's member resolution: one row per distinct protein-group member, saying which
    /// annotation-database entry it takes.
    /// <code>
    /// #!member_resolution_format 1
    /// #!mzlib_version 1.0.596+...
    /// #!mzlib_release 1.0.596
    /// #!annotation_db_sha256 ...
    /// #!source_file_sha256 ...
    /// accession	annotated_as	resolution	entrapment	sequence_accession
    /// P04406	P04406	self	false	P04406
    /// P04406-2	P04406	inherited	false	P04406-2
    /// P04406-2_A20T	P04406	inherited	false	P04406-2
    /// Random_P12345_f0		no_entry	true
    /// </code>
    /// It holds no annotation, and it does not restate group membership: the search's protein-group table
    /// already lists each group's members. A view joins the groups' members to this table on accession, and
    /// annotated_as to any per-database annotation table (GO, pathways) on accession. An annotation that names
    /// one isoform (a Reactome reference's molecule) joins on sequence_accession instead. The resolution
    /// depends on the annotation database, so the header names it.
    /// </summary>
    public static class MemberResolutionTsv
    {
        /// <summary>The format version written on the first line.</summary>
        public const int FormatVersion = 1;

        private const string FormatKey = "member_resolution_format";

        /// <summary>The column names, in order.</summary>
        public static readonly IReadOnlyList<string> Columns = new[] { "accession", "annotated_as", "resolution", "entrapment", "sequence_accession" };

        /// <summary>Writes the rows in the order given, with the header.</summary>
        /// <param name="annotationDbSha256">sha256 of the annotation database the members were resolved against.</param>
        /// <param name="sourceFileSha256">sha256 of the protein-group table the members came from, or null to omit the line.</param>
        /// <exception cref="ArgumentNullException">output or rows is null.</exception>
        /// <exception cref="ArgumentException">annotationDbSha256 is empty; an accession repeats or contains a
        /// tab or line break; or a row's annotated_as or sequence_accession contradicts its resolution or
        /// accession.</exception>
        public static void Write(TextWriter output, IEnumerable<MemberResolution> rows, string annotationDbSha256,
            string sourceFileSha256 = null)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(rows);
            if (string.IsNullOrEmpty(annotationDbSha256))
            {
                throw new ArgumentException("The annotation database's sha256 is required: the resolution depends on it.", nameof(annotationDbSha256));
            }
            var materialized = rows.ToList();
            var seen = new HashSet<string>(StringComparer.Ordinal);
            foreach (var row in materialized)
            {
                if (!seen.Add(row.Accession))
                {
                    throw new ArgumentException($"Member '{row.Accession}' appears twice.");
                }
                Check(row, row.Accession, message => new ArgumentException(message));
            }

            TsvHeader.Write(output, FormatKey, FormatVersion.ToString(CultureInfo.InvariantCulture));
            TsvHeader.WriteProducer(output);
            TsvHeader.Write(output, "annotation_db_sha256", annotationDbSha256);
            TsvHeader.Write(output, "source_file_sha256", sourceFileSha256);
            output.Write(string.Join('\t', Columns) + "\n");
            foreach (var row in materialized)
            {
                TsvHeader.RejectSeparators(row.Accession, "accession");
                TsvHeader.RejectSeparators(row.AnnotatedAs, "annotated_as");
                output.Write($"{row.Accession}\t{row.AnnotatedAs}\t{KindName(row.Kind)}\t{(row.IsEntrapment ? "true" : "false")}\t{row.SequenceAccession}" + "\n");
            }
        }

        /// <summary>Reads a member resolution file.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's; a row's
        /// resolution or entrapment is not a known value; annotated_as or sequence_accession contradicts the
        /// resolution or accession; or an accession repeats.</exception>
        public static AnnotationTable<MemberResolution> Read(string path)
        {
            var seen = new HashSet<string>(StringComparer.Ordinal);
            return AnnotationTableReader.Read(path, FormatKey, FormatVersion,
                new[] { "annotation_db_sha256" }, new[] { "source_file_sha256" }, Columns,
                (cells, fileName, lineNumber) =>
                {
                    string where = $"{fileName} line {lineNumber}";
                    if (cells[0].Length == 0)
                    {
                        throw new InvalidDataException($"{where}: accession is empty.");
                    }
                    if (!seen.Add(cells[0]))
                    {
                        throw new InvalidDataException($"{where}: member '{cells[0]}' appears twice.");
                    }
                    var kind = cells[2] switch
                    {
                        "self" => MemberResolutionKind.Self,
                        "inherited" => MemberResolutionKind.Inherited,
                        "no_entry" => MemberResolutionKind.NoEntry,
                        _ => throw new InvalidDataException($"{where}: resolution '{cells[2]}' is not self, inherited or no_entry.")
                    };
                    bool entrapment = AnnotationTableReader.Bool(cells[3], "entrapment", fileName, lineNumber)
                        ?? throw new InvalidDataException($"{where}: entrapment is empty.");
                    var row = new MemberResolution(cells[0], AnnotationTableReader.Cell(cells[1]), kind, entrapment,
                        AnnotationTableReader.Cell(cells[4]));
                    Check(row, where, message => new InvalidDataException(message));
                    return row;
                });
        }

        /// <summary>"self", "inherited" or "no_entry".</summary>
        public static string KindName(MemberResolutionKind kind) => kind switch
        {
            MemberResolutionKind.Self => "self",
            MemberResolutionKind.Inherited => "inherited",
            MemberResolutionKind.NoEntry => "no_entry",
            _ => throw new ArgumentOutOfRangeException(nameof(kind), kind, null)
        };

        /// <summary>
        /// self names the member itself, inherited another entry, no_entry none; and the sequence accession is
        /// the one the member's accession gives, or none for no_entry.
        /// </summary>
        private static void Check(MemberResolution row, string where, Func<string, Exception> fail)
        {
            bool consistent = row.Kind switch
            {
                MemberResolutionKind.Self => row.AnnotatedAs == row.Accession,
                MemberResolutionKind.Inherited => !string.IsNullOrEmpty(row.AnnotatedAs) && row.AnnotatedAs != row.Accession,
                MemberResolutionKind.NoEntry => string.IsNullOrEmpty(row.AnnotatedAs),
                _ => false
            };
            if (!consistent)
            {
                throw fail($"{where}: annotated_as '{row.AnnotatedAs}' contradicts resolution '{KindName(row.Kind)}'.");
            }
            string sequence = row.Kind == MemberResolutionKind.NoEntry ? null : MemberResolver.SequenceAccessionOf(row.Accession);
            if (row.SequenceAccession != sequence)
            {
                throw fail($"{where}: sequence_accession '{row.SequenceAccession}' contradicts accession '{row.Accession}' ({KindName(row.Kind)}).");
            }
        }
    }
}
