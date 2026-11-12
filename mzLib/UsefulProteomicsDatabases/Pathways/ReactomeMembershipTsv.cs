using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
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
    /// Each map is one line, its entries ';'-joined as organism=value in ordinal order, because the reader refuses a
    /// repeated header key.
    /// </summary>
    public static class ReactomeMembershipTsv
    {
        /// <summary>The format version written on the first line.</summary>
        public const int FormatVersion = 1;

        private const string FormatKey = "reactome_membership_format";

        private static readonly Regex OrganismCurie = new(@"^NCBITaxon:\d+$", RegexOptions.Compiled);

        /// <summary>The column names, in order.</summary>
        public static readonly IReadOnlyList<string> Columns =
            new[] { "accession", "organism", "pathway_id", "molecule", "sequence_accession", "inherited", "status" };

        private static string Definition =>
            $"{ReactomeAccessionAnnotation.DefinitionId} v{ReactomeAccessionAnnotation.DefinitionVersion}";

        /// <summary>"annotated", "no_reference" or "species_not_covered".</summary>
        /// <exception cref="ArgumentOutOfRangeException">The value is not a defined status.</exception>
        public static string StatusName(ReactomeMembershipStatus status) => status switch
        {
            ReactomeMembershipStatus.Annotated => "annotated",
            ReactomeMembershipStatus.NoReference => "no_reference",
            ReactomeMembershipStatus.SpeciesNotCovered => "species_not_covered",
            _ => throw new ArgumentOutOfRangeException(nameof(status), status, "Not a Reactome membership status.")
        };

        /// <summary>
        /// Writes every row of <paramref name="annotation"/>, with the header. Every row is checked before anything is
        /// written.
        /// </summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">A value contains a tab or a line break, an accession is empty, or an
        /// organism is not an NCBITaxon CURIE: <see cref="Read"/> would refuse the file.</exception>
        public static void Write(TextWriter output, ReactomeAccessionAnnotation annotation)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(annotation);
            var rows = annotation.Rows();
            foreach (var row in rows)
            {
                if (string.IsNullOrEmpty(row.Accession))
                {
                    throw new ArgumentException("An accession is empty.", nameof(annotation));
                }
                if (row.Organism != null && !OrganismCurie.IsMatch(row.Organism))
                {
                    throw new ArgumentException($"The organism '{row.Organism}' of '{row.Accession}' is not an NCBITaxon CURIE.",
                        nameof(annotation));
                }
                TsvHeader.RejectSeparators(row.Accession, "accession");
                TsvHeader.RejectSeparators(row.PathwayId, "pathway_id");
                TsvHeader.RejectSeparators(row.Molecule, "molecule");
                TsvHeader.RejectSeparators(row.SequenceAccession, "sequence_accession");
            }

            TsvHeader.Write(output, FormatKey, FormatVersion.ToString(CultureInfo.InvariantCulture));
            TsvHeader.WriteProducer(output);
            TsvHeader.Write(output, "definition", Definition);
            TsvHeader.Write(output, "annotation_db_sha256", annotation.AnnotationDbSha256);
            TsvHeader.Write(output, "species_codes", Map(annotation.SpeciesCodes));
            TsvHeader.Write(output, "unvalidated_organisms",
                Map(annotation.UnvalidatedOrganisms.ToDictionary(kv => kv.Key, kv => kv.Value.ToString(CultureInfo.InvariantCulture))));
            TsvHeader.Write(output, "dropped_foreign_references", annotation.DroppedForeignReferences > 0
                ? annotation.DroppedForeignReferences.ToString(CultureInfo.InvariantCulture) : null);
            output.Write(string.Join('\t', Columns) + "\n");
            foreach (var row in rows)
            {
                string inherited = row.Inherited switch { true => "true", false => "false", null => "" };
                output.Write($"{row.Accession}\t{row.Organism}\t{row.PathwayId}\t{row.Molecule}\t{row.SequenceAccession}\t" +
                    $"{inherited}\t{StatusName(row.Status)}" + "\n");
            }
        }

        /// <summary>organism=value pairs in ordinal order, ';'-joined; null when there are none.</summary>
        private static string Map(IEnumerable<KeyValuePair<string, string>> pairs)
        {
            var entries = pairs.OrderBy(kv => kv.Key, StringComparer.Ordinal).Select(kv => $"{kv.Key}={kv.Value}").ToList();
            return entries.Count == 0 ? null : string.Join(';', entries);
        }

        /// <summary>Reads a database's Reactome membership file.</summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The header or columns are not this format's, or the definition is
        /// another; a status is unknown; a row that carries a pathway is not annotated, or an annotated row carries
        /// none; an entry has a pathway-less row and another row; a sequence_accession disagrees with its rule; an
        /// inherited row is isoform-scoped; or a row repeats.</exception>
        public static AnnotationTable<ReactomeMembershipRow> Read(string path)
        {
            var seen = new HashSet<(string, string, string)>();
            var pathwayless = new HashSet<string>(StringComparer.Ordinal);
            var withPathway = new HashSet<string>(StringComparer.Ordinal);
            var table = AnnotationTableReader.Read(path, FormatKey, FormatVersion,
                new[] { "definition", "annotation_db_sha256" },
                new[] { "species_codes", "unvalidated_organisms", "dropped_foreign_references" }, Columns,
                (cells, fileName, lineNumber) =>
                {
                    string where = $"{fileName} line {lineNumber}";
                    string accession = cells[0];
                    if (accession.Length == 0)
                    {
                        throw new InvalidDataException($"{where}: accession is empty.");
                    }
                    string organism = AnnotationTableReader.Cell(cells[1]);
                    if (organism != null && !OrganismCurie.IsMatch(organism))
                    {
                        throw new InvalidDataException($"{where}: organism '{organism}' is not an NCBITaxon CURIE (NCBITaxon:9606).");
                    }
                    string pathwayId = AnnotationTableReader.Cell(cells[2]);
                    string molecule = AnnotationTableReader.Cell(cells[3]);
                    string sequenceAccession = AnnotationTableReader.Cell(cells[4]);
                    bool? inherited = AnnotationTableReader.Bool(cells[5], "inherited", fileName, lineNumber);
                    var status = ParseStatus(cells[6], where);

                    if (pathwayId == null)
                    {
                        if (status == ReactomeMembershipStatus.Annotated)
                        {
                            throw new InvalidDataException($"{where}: an annotated row of '{accession}' carries no pathway.");
                        }
                        if (molecule != null || sequenceAccession != null || inherited != null)
                        {
                            throw new InvalidDataException($"{where}: a pathway-less row of '{accession}' carries a value.");
                        }
                        if (withPathway.Contains(accession) || !pathwayless.Add(accession))
                        {
                            throw new InvalidDataException($"{where}: '{accession}' has a pathway-less row and another row.");
                        }
                    }
                    else
                    {
                        if (status != ReactomeMembershipStatus.Annotated)
                        {
                            throw new InvalidDataException(
                                $"{where}: the row for '{accession}', {pathwayId} carries a pathway but reads '{cells[6]}'.");
                        }
                        if (inherited == null)
                        {
                            throw new InvalidDataException($"{where}: inherited is empty on the row for '{accession}', {pathwayId}.");
                        }
                        if (inherited == true && molecule != null)
                        {
                            throw new InvalidDataException($"{where}: an inherited row is isoform-scoped ('{accession}', {pathwayId}).");
                        }
                        bool sequenceFollowsRule = molecule == null
                            ? sequenceAccession == null
                            : sequenceAccession == molecule || sequenceAccession == MemberResolver.SequenceAccessionOf(accession);
                        if (!sequenceFollowsRule)
                        {
                            throw new InvalidDataException(
                                $"{where}: sequence_accession '{cells[4]}' is not the molecule or the entry's own sequence " +
                                $"accession, or is set on a whole-entry row ('{accession}', {pathwayId}).");
                        }
                        if (pathwayless.Contains(accession))
                        {
                            throw new InvalidDataException($"{where}: '{accession}' has a pathway-less row and another row.");
                        }
                        if (!seen.Add((accession, pathwayId, molecule)))
                        {
                            throw new InvalidDataException($"{where}: '{accession}', {pathwayId}, {molecule ?? "whole entry"} appears twice.");
                        }
                        withPathway.Add(accession);
                    }
                    return new ReactomeMembershipRow(accession, organism, pathwayId, molecule, sequenceAccession, inherited, status);
                });

            if (table.Header["definition"] != Definition)
            {
                throw new InvalidDataException(
                    $"{table.SourceFileName}: definition must be '{Definition}', found '{table.Header["definition"]}'.");
            }
            return table;
        }

        private static ReactomeMembershipStatus ParseStatus(string value, string where) => value switch
        {
            "annotated" => ReactomeMembershipStatus.Annotated,
            "no_reference" => ReactomeMembershipStatus.NoReference,
            "species_not_covered" => ReactomeMembershipStatus.SpeciesNotCovered,
            _ => throw new InvalidDataException($"{where}: status '{value}' is not annotated, no_reference or species_not_covered.")
        };
    }
}
