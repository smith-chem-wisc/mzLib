using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases.Annotation;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>
    /// Writes and reads one database's GO annotation (<see cref="GoAccessionAnnotation"/>): one row per
    /// (entry, directly annotated term), and one term-less row per entry with no term.
    /// <code>
    /// #!go_accession_terms_format 1
    /// #!mzlib_version 1.0.596+...
    /// #!mzlib_release 1.0.596
    /// #!go_release releases/2026-07-26
    /// #!go_obo_sha256 ...
    /// #!annotation_db_sha256 ...
    /// accession	go_id	aspect	evidence	inherited
    /// P04406	GO:0005634	cellular_component	ECO:0000314	false
    /// P04406-2	GO:0005634	cellular_component	ECO:0000314	true
    /// Q9XXXX
    /// </code>
    /// The values every row shares are in the header, once. Ancestors are not rows: join go_id to the
    /// release's closure (<see cref="GoOntologyTsv.WriteClosure"/>), and the term's name to its term table.
    /// Evidence is ';'-joined; an id or code containing ';' is refused.
    /// </summary>
    public static class GoAccessionTermsTsv
    {
        /// <summary>The format version written on the first line.</summary>
        public const int FormatVersion = 1;

        private const string FormatKey = "go_accession_terms_format";

        /// <summary>The column names, in order.</summary>
        public static readonly IReadOnlyList<string> Columns = new[] { "accession", "go_id", "aspect", "evidence", "inherited" };

        /// <summary>Writes every row of <paramref name="annotation"/>, with the header.</summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">A value contains a tab, a line break or, in a set, ';'.</exception>
        public static void Write(TextWriter output, GoAccessionAnnotation annotation)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(annotation);
            var ontology = annotation.Ontology;

            TsvHeader.Write(output, FormatKey, FormatVersion.ToString(CultureInfo.InvariantCulture));
            TsvHeader.WriteProducer(output);
            TsvHeader.Write(output, "go_release", ontology.Release);
            TsvHeader.Write(output, "go_obo_sha256", ontology.SourceSha256);
            TsvHeader.Write(output, "annotation_db_sha256", annotation.AnnotationDbSha256);
            TsvHeader.Write(output, "unresolved_go_ids", annotation.UnresolvedGoIds.Count > 0
                ? annotation.UnresolvedGoIds.Count.ToString(CultureInfo.InvariantCulture) : null);
            output.Write(string.Join('\t', Columns) + "\n");
            foreach (var row in annotation.Rows())
            {
                TsvHeader.RejectSeparators(row.Accession, "accession");
                if (row.Evidence.Any(e => e.Contains(';')))
                {
                    throw new ArgumentException($"An evidence code of '{row.Accession}' contains ';', the set separator.");
                }
                string evidence = string.Join(';', row.Evidence);
                TsvHeader.RejectSeparators(evidence, "evidence");
                string inherited = row.Inherited switch { true => "true", false => "false", null => "" };
                output.Write($"{row.Accession}\t{row.GoId}\t{GoAnnotationTsv.AspectName(row.Aspect)}\t{evidence}\t{inherited}" + "\n");
            }
        }

        /// <summary>Reads a database's GO annotation file.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's; a value is not
        /// valid for its column; a term-less row carries a value; a termed row lacks its aspect or inherited flag;
        /// an (accession, go_id) pair repeats; or an entry has both a term-less row and terms.</exception>
        public static AnnotationTable<GoAccessionTermRow> Read(string path)
        {
            var seen = new HashSet<(string, string)>();
            var termless = new HashSet<string>(StringComparer.Ordinal);
            var termed = new HashSet<string>(StringComparer.Ordinal);
            return AnnotationTableReader.Read(path, FormatKey, FormatVersion,
                new[] { "go_obo_sha256", "annotation_db_sha256" }, new[] { "go_release", "unresolved_go_ids" }, Columns,
                (cells, fileName, lineNumber) =>
                {
                    string where = $"{fileName} line {lineNumber}";
                    string accession = cells[0];
                    if (accession.Length == 0)
                    {
                        throw new InvalidDataException($"{where}: accession is empty.");
                    }
                    string goId = AnnotationTableReader.Cell(cells[1]);
                    var aspect = GoAnnotationTsv.ParseAspect(cells[2], where);
                    var evidence = AnnotationTableReader.Set(cells[3]);
                    bool? inherited = AnnotationTableReader.Bool(cells[4], "inherited", fileName, lineNumber);

                    if (goId == null)
                    {
                        if (aspect != null || evidence.Count > 0 || inherited != null)
                        {
                            throw new InvalidDataException($"{where}: a term-less row of '{accession}' carries a value.");
                        }
                        if (termed.Contains(accession) || !termless.Add(accession))
                        {
                            throw new InvalidDataException($"{where}: '{accession}' has a term-less row and another row.");
                        }
                    }
                    else
                    {
                        if (aspect == null || inherited == null)
                        {
                            throw new InvalidDataException($"{where}: the row for '{accession}', {goId} lacks its aspect or inherited flag.");
                        }
                        if (termless.Contains(accession))
                        {
                            throw new InvalidDataException($"{where}: '{accession}' has a term-less row and another row.");
                        }
                        if (!seen.Add((accession, goId)))
                        {
                            throw new InvalidDataException($"{where}: '{accession}', {goId} appears twice.");
                        }
                        termed.Add(accession);
                    }
                    return new GoAccessionTermRow(accession, goId, aspect, evidence, inherited);
                });
        }
    }
}
