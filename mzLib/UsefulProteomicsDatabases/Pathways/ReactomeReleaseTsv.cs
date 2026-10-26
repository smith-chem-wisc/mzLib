using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases.Annotation;
using UsefulProteomicsDatabases.GeneOntology;

namespace UsefulProteomicsDatabases.Pathways
{
    /// <summary>One pathway of a Reactome release: its id, name and species as Reactome writes them, and the human
    /// pathway it is (<see cref="ReactomeRelease.HumanEquivalentId"/>), or null.</summary>
    public sealed record ReactomePathwayRow(string PathwayId, string Name, string SpeciesName, string HumanEquivalentId);

    /// <summary>One parent-child pair of a Reactome release's hierarchy.</summary>
    public sealed record ReactomeHierarchyRow(string ParentId, string ChildId);

    /// <summary>One pair of a Reactome release's ancestor closure: <paramref name="AncestorPathwayId"/> is the pathway
    /// itself or a pathway above it.</summary>
    public sealed record ReactomeAncestorRow(string PathwayId, string AncestorPathwayId);

    /// <summary>
    /// Writes and reads the three tables that depend on the Reactome release alone, so a store keeps them once per
    /// release rather than once per database or per search -- the Reactome counterpart of GoOntologyTsv:
    /// <list type="bullet">
    /// <item><b>pathways</b> (reactome_pathways_format 1): pathway_id, name, species_name, human_equivalent_id (empty
    /// when the pathway has no human equivalent).</item>
    /// <item><b>hierarchy</b> (reactome_hierarchy_format 1): parent_id, child_id.</item>
    /// <item><b>ancestors</b> (reactome_ancestors_format 1): every (pathway_id, ancestor_pathway_id) pair, the pathway
    /// itself included, so joining direct memberships through it keeps the direct pathway; named like GO's closure
    /// (go_id, ancestor_go_id).</item>
    /// </list>
    /// Every table's header carries its format, the mzLib build that wrote it, the release and both source files'
    /// sha256s. A value every row shares is a header line, never a column.
    /// </summary>
    public static class ReactomeReleaseTsv
    {
        /// <summary>The format version of all three tables.</summary>
        public const int FormatVersion = 1;

        private const string PathwaysFormatKey = "reactome_pathways_format";
        private const string HierarchyFormatKey = "reactome_hierarchy_format";
        private const string AncestorsFormatKey = "reactome_ancestors_format";

        private static readonly string[] RequiredKeys = { "reactome_release", "reactome_pathways_sha256", "reactome_relation_sha256" };

        /// <summary>The pathway table's column names, in order.</summary>
        public static readonly IReadOnlyList<string> PathwayColumns = new[] { "pathway_id", "name", "species_name", "human_equivalent_id" };

        /// <summary>The hierarchy's column names, in order.</summary>
        public static readonly IReadOnlyList<string> HierarchyColumns = new[] { "parent_id", "child_id" };

        /// <summary>The ancestor closure's column names, in order.</summary>
        public static readonly IReadOnlyList<string> AncestorColumns = new[] { "pathway_id", "ancestor_pathway_id" };

        /// <summary>Writes every pathway of <paramref name="release"/>, in ordinal order of id.</summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">A name contains a tab or line break.</exception>
        public static void WritePathways(TextWriter output, ReactomeRelease release)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(release);
            WriteHeader(output, PathwaysFormatKey, PathwayColumns, release);
            foreach (var pathway in release.Pathways)
            {
                TsvHeader.RejectSeparators(pathway.Name, "name");
                TsvHeader.RejectSeparators(pathway.SpeciesName, "species_name");
                output.Write($"{pathway.PathwayId}\t{pathway.Name}\t{pathway.SpeciesName}\t{pathway.HumanEquivalentId}" + "\n");
            }
        }

        /// <summary>Writes every parent-child pair of <paramref name="release"/>, ordered by parent then child.</summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        public static void WriteHierarchy(TextWriter output, ReactomeRelease release)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(release);
            WriteHeader(output, HierarchyFormatKey, HierarchyColumns, release);
            foreach (var edge in release.Hierarchy)
            {
                output.Write($"{edge.ParentId}\t{edge.ChildId}" + "\n");
            }
        }

        /// <summary>
        /// Writes every (pathway, ancestor) pair of <paramref name="release"/>, the pathway itself included, ordered by
        /// pathway then ancestor (ordinal).
        /// </summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        public static void WriteAncestors(TextWriter output, ReactomeRelease release)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(release);
            WriteHeader(output, AncestorsFormatKey, AncestorColumns, release);
            foreach (var pathway in release.Pathways)
            {
                foreach (string ancestor in release.Ancestors(pathway.PathwayId).Append(pathway.PathwayId).OrderBy(a => a, StringComparer.Ordinal))
                {
                    output.Write($"{pathway.PathwayId}\t{ancestor}" + "\n");
                }
            }
        }

        /// <summary>Reads a pathway table. An empty human_equivalent_id reads as null.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's, or a pathway id is empty
        /// or repeats.</exception>
        public static AnnotationTable<ReactomePathwayRow> ReadPathways(string path)
        {
            var seen = new HashSet<string>(StringComparer.Ordinal);
            return AnnotationTableReader.Read(path, PathwaysFormatKey, FormatVersion, RequiredKeys, Array.Empty<string>(), PathwayColumns,
                (cells, fileName, lineNumber) =>
                {
                    if (cells[0].Length == 0 || !seen.Add(cells[0]))
                    {
                        throw new InvalidDataException($"{fileName} line {lineNumber}: pathway_id '{cells[0]}' is empty or repeats.");
                    }
                    return new ReactomePathwayRow(cells[0], cells[1], cells[2], AnnotationTableReader.Cell(cells[3]));
                });
        }

        /// <summary>Reads a hierarchy.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's, a cell is empty, or a
        /// pair repeats.</exception>
        public static AnnotationTable<ReactomeHierarchyRow> ReadHierarchy(string path)
        {
            var seen = new HashSet<(string, string)>();
            return AnnotationTableReader.Read(path, HierarchyFormatKey, FormatVersion, RequiredKeys, Array.Empty<string>(), HierarchyColumns,
                (cells, fileName, lineNumber) =>
                {
                    RequirePair(cells, seen, fileName, lineNumber, "hierarchy");
                    return new ReactomeHierarchyRow(cells[0], cells[1]);
                });
        }

        /// <summary>Reads an ancestor closure.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's; a cell is empty; a pair
        /// repeats; or a pathway lacks its own (pathway, pathway) pair.</exception>
        public static AnnotationTable<ReactomeAncestorRow> ReadAncestors(string path)
        {
            var seen = new HashSet<(string, string)>();
            var table = AnnotationTableReader.Read(path, AncestorsFormatKey, FormatVersion, RequiredKeys, Array.Empty<string>(), AncestorColumns,
                (cells, fileName, lineNumber) =>
                {
                    RequirePair(cells, seen, fileName, lineNumber, "ancestors");
                    return new ReactomeAncestorRow(cells[0], cells[1]);
                });
            var withoutSelf = table.Rows.Select(r => r.PathwayId).Distinct(StringComparer.Ordinal)
                .FirstOrDefault(id => !seen.Contains((id, id)));
            if (withoutSelf != null)
            {
                throw new InvalidDataException($"{table.SourceFileName}: {withoutSelf} has no ({withoutSelf}, {withoutSelf}) pair, so its direct memberships would be lost in a join.");
            }
            return table;
        }

        private static void RequirePair(string[] cells, HashSet<(string, string)> seen, string fileName, int lineNumber, string table)
        {
            string where = $"{fileName} line {lineNumber}";
            if (cells[0].Length == 0 || cells[1].Length == 0)
            {
                throw new InvalidDataException($"{where}: a cell of the {table} table is empty.");
            }
            if (!seen.Add((cells[0], cells[1])))
            {
                throw new InvalidDataException($"{where}: {cells[0]}, {cells[1]} appears twice.");
            }
        }

        private static void WriteHeader(TextWriter output, string formatKey, IReadOnlyList<string> columns, ReactomeRelease release)
        {
            TsvHeader.Write(output, formatKey, FormatVersion.ToString(CultureInfo.InvariantCulture));
            TsvHeader.WriteProducer(output);
            TsvHeader.Write(output, "reactome_release", release.Release);
            TsvHeader.Write(output, "reactome_pathways_sha256", release.PathwaysSha256);
            TsvHeader.Write(output, "reactome_relation_sha256", release.RelationSha256);
            output.Write(string.Join('\t', columns) + "\n");
        }
    }
}
