using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Proteomics;
using UsefulProteomicsDatabases.Annotation;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>One term of an ontology release: its primary id, name, aspect and whether it is obsolete.</summary>
    public sealed record GoTermRow(string GoId, string GoName, GoAspect Aspect, bool IsObsolete);

    /// <summary>One pair of an ontology release's closure: <paramref name="AncestorGoId"/> is the term itself or one
    /// of its ancestors over is_a and part_of.</summary>
    public sealed record GoClosureRow(string GoId, string AncestorGoId);

    /// <summary>
    /// Writes and reads the two tables that depend on the ontology release alone, so a store keeps them once
    /// per go.obo rather than once per database or per search:
    /// <list type="bullet">
    /// <item><b>terms</b> (go_terms_format 1): every term, obsolete ones included -- go_id, go_name, aspect,
    /// obsolete.</item>
    /// <item><b>closure</b> (go_closure_format 1): every (go_id, ancestor_go_id) pair, where the ancestor is the
    /// term itself or a term above it over is_a and part_of. The term itself is always a pair, so a join of
    /// direct annotations through the closure keeps the direct terms; an obsolete term has only that pair, so
    /// nothing propagates from it.</item>
    /// </list>
    /// Joining a database's direct annotations (<see cref="GoAccessionTermsTsv"/>) to the closure on go_id gives
    /// exactly the terms <see cref="GoGroupAnnotator"/> propagates to: the pair's ancestor is the propagated term,
    /// and the direct term's evidence is carried to it.
    /// </summary>
    public static class GoOntologyTsv
    {
        /// <summary>The format version of both tables.</summary>
        public const int FormatVersion = 1;

        private const string TermsFormatKey = "go_terms_format";
        private const string ClosureFormatKey = "go_closure_format";

        /// <summary>The term table's column names, in order.</summary>
        public static readonly IReadOnlyList<string> TermColumns = new[] { "go_id", "go_name", "aspect", "obsolete" };

        /// <summary>The closure's column names, in order.</summary>
        public static readonly IReadOnlyList<string> ClosureColumns = new[] { "go_id", "ancestor_go_id" };

        /// <summary>Writes every term of <paramref name="ontology"/>, in ordinal order of id.</summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">A term's name contains a tab or line break.</exception>
        public static void WriteTerms(TextWriter output, GeneOntologyGraph ontology)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(ontology);
            WriteHeader(output, TermsFormatKey, ontology);
            output.Write(string.Join('\t', TermColumns) + "\n");
            foreach (string id in ontology.TermIds)
            {
                ontology.TryGetTerm(id, out var term);
                TsvHeader.RejectSeparators(term.Name, "go_name");
                output.Write($"{term.Id}\t{term.Name}\t{GoAnnotationTsv.AspectName(term.Aspect)}\t{(term.IsObsolete ? "true" : "false")}" + "\n");
            }
        }

        /// <summary>
        /// Writes every (term, ancestor) pair of <paramref name="ontology"/>, the term itself included, ordered by
        /// term then ancestor (ordinal).
        /// </summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        public static void WriteClosure(TextWriter output, GeneOntologyGraph ontology)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(ontology);
            WriteHeader(output, ClosureFormatKey, ontology);
            output.Write(string.Join('\t', ClosureColumns) + "\n");
            foreach (string id in ontology.TermIds)
            {
                foreach (string ancestor in ontology.Ancestors(id).Append(id).OrderBy(a => a, StringComparer.Ordinal))
                {
                    output.Write($"{id}\t{ancestor}" + "\n");
                }
            }
        }

        /// <summary>Reads a term table.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's; a cell is empty or
        /// invalid; or a term repeats.</exception>
        public static AnnotationTable<GoTermRow> ReadTerms(string path)
        {
            var seen = new HashSet<string>(StringComparer.Ordinal);
            return AnnotationTableReader.Read(path, TermsFormatKey, FormatVersion,
                new[] { "go_obo_sha256" }, new[] { "go_release" }, TermColumns,
                (cells, fileName, lineNumber) =>
                {
                    string where = $"{fileName} line {lineNumber}";
                    if (cells[0].Length == 0 || !seen.Add(cells[0]))
                    {
                        throw new InvalidDataException($"{where}: go_id '{cells[0]}' is empty or repeats.");
                    }
                    var aspect = GoAnnotationTsv.ParseAspect(cells[2], where)
                        ?? throw new InvalidDataException($"{where}: aspect is empty.");
                    bool obsolete = AnnotationTableReader.Bool(cells[3], "obsolete", fileName, lineNumber)
                        ?? throw new InvalidDataException($"{where}: obsolete is empty.");
                    return new GoTermRow(cells[0], cells[1], aspect, obsolete);
                });
        }

        /// <summary>Reads a closure.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's; a cell is empty;
        /// a pair repeats; or a term lacks its own (term, term) pair.</exception>
        public static AnnotationTable<GoClosureRow> ReadClosure(string path)
        {
            var seen = new HashSet<(string, string)>();
            var table = AnnotationTableReader.Read(path, ClosureFormatKey, FormatVersion,
                new[] { "go_obo_sha256" }, new[] { "go_release" }, ClosureColumns,
                (cells, fileName, lineNumber) =>
                {
                    string where = $"{fileName} line {lineNumber}";
                    if (cells[0].Length == 0 || cells[1].Length == 0)
                    {
                        throw new InvalidDataException($"{where}: a closure cell is empty.");
                    }
                    if (!seen.Add((cells[0], cells[1])))
                    {
                        throw new InvalidDataException($"{where}: {cells[0]}, {cells[1]} appears twice.");
                    }
                    return new GoClosureRow(cells[0], cells[1]);
                });
            var withoutSelf = table.Rows.Select(r => r.GoId).Distinct(StringComparer.Ordinal)
                .FirstOrDefault(id => !seen.Contains((id, id)));
            if (withoutSelf != null)
            {
                throw new InvalidDataException($"{table.SourceFileName}: {withoutSelf} has no ({withoutSelf}, {withoutSelf}) pair, so its direct annotations would be lost in a join.");
            }
            return table;
        }

        private static void WriteHeader(TextWriter output, string formatKey, GeneOntologyGraph ontology)
        {
            TsvHeader.Write(output, formatKey, FormatVersion.ToString(CultureInfo.InvariantCulture));
            TsvHeader.WriteProducer(output);
            TsvHeader.Write(output, "go_release", ontology.Release);
            TsvHeader.Write(output, "go_obo_sha256", ontology.SourceSha256);
        }
    }
}
