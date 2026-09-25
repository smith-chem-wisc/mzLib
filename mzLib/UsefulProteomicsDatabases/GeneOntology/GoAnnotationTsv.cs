using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using Omics.BioPolymerGroup;
using Proteomics;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>
    /// Writes GO annotation rows as a tab-separated file: "#!key value" header lines, a column header, then
    /// one line per (protein group, term) row.
    ///
    /// The header states the provenance every row shares -- format version first, then the mzLib build that
    /// wrote the file and its release ("none" for a local build), then the go.obo release
    /// and sha256, the annotation database's sha256 and, when there is one, the sha256 of the results file
    /// the groups were read from -- and five run counters: the number of multi-member groups and the number
    /// of groups in each annotation status. The counters count GROUPS, not rows, and only groups at
    /// q_value &lt;= <see cref="CounterQValueMax"/>; the header says so on its own line. The rows are not
    /// filtered: every group is in the file, and the consumer filters.
    ///
    /// This is a data-interchange table, so its column names and header keys are a contract with whoever
    /// ingests it. Long format: one term per row. Set-valued cells (accession_used, evidence) are
    /// ';'-joined, and a member containing ';' is refused rather than written ambiguously. Groups keep input
    /// order; terms within a group are in ordinal order; lines end in "\n" on every platform.
    /// </summary>
    public static class GoAnnotationTsv
    {
        /// <summary>The format version written on the first line.</summary>
        public const int FormatVersion = 1;

        /// <summary>The q-value at or below which a group counts toward the header counters.</summary>
        public const double CounterQValueMax = 0.01;

        /// <summary>The row layout, in column order.</summary>
        public static readonly IReadOnlyList<TsvColumn<GoAnnotationRow>> Schema = new[]
        {
            new TsvColumn<GoAnnotationRow>("protein_group", r => r.ProteinGroup),
            new TsvColumn<GoAnnotationRow>("accession_used", r => JoinSet(r.AccessionUsed)),
            new TsvColumn<GoAnnotationRow>("go_id", r => r.GoId),
            new TsvColumn<GoAnnotationRow>("go_name", r => r.GoName),
            new TsvColumn<GoAnnotationRow>("aspect", r => AspectName(r.Aspect)),
            new TsvColumn<GoAnnotationRow>("evidence", r => JoinSet(r.Evidence)),
            new TsvColumn<GoAnnotationRow>("inherited", r => Bool(r.Inherited)),
            new TsvColumn<GoAnnotationRow>("propagated", r => Bool(r.Propagated)),
            new TsvColumn<GoAnnotationRow>("n_members", r => r.NMembers.ToString(CultureInfo.InvariantCulture)),
            new TsvColumn<GoAnnotationRow>("n_with", r => r.NWith.ToString(CultureInfo.InvariantCulture)),
            new TsvColumn<GoAnnotationRow>("annotation_status", r => StatusName(r.Status)),
            new TsvColumn<GoAnnotationRow>("q_value", r => r.QValue.ToString("R", CultureInfo.InvariantCulture)),
            new TsvColumn<GoAnnotationRow>("go_release", r => r.GoRelease),
            new TsvColumn<GoAnnotationRow>("go_obo_sha256", r => r.GoOboSha256),
            new TsvColumn<GoAnnotationRow>("annotation_db_sha256", r => r.AnnotationDbSha256),
        };

        /// <summary>
        /// Writes <paramref name="rows"/> with the header. Every row must come from <paramref name="ontology"/>
        /// and the annotation database <paramref name="annotationDbSha256"/>: a file mixing releases would
        /// describe itself falsely.
        /// </summary>
        /// <param name="sourceFileSha256">sha256 of the results file the groups came from, or null to omit the line.</param>
        /// <exception cref="ArgumentNullException">output, rows or ontology is null.</exception>
        /// <exception cref="ArgumentException">A row's release, go.obo sha256 or annotation database sha256
        /// differs from the header's; a value contains a tab or line break; or a set member contains ';'.</exception>
        public static void Write(TextWriter output, IEnumerable<GoAnnotationRow> rows, GeneOntologyGraph ontology,
            string annotationDbSha256, string sourceFileSha256 = null)
        {
            ArgumentNullException.ThrowIfNull(output);
            ArgumentNullException.ThrowIfNull(rows);
            ArgumentNullException.ThrowIfNull(ontology);

            var materialized = rows.ToList();
            foreach (var row in materialized)
            {
                RequireSame(row.GoRelease, ontology.Release, "go_release", row);
                RequireSame(row.GoOboSha256, ontology.SourceSha256, "go_obo_sha256", row);
                RequireSame(row.AnnotationDbSha256, annotationDbSha256, "annotation_db_sha256", row);
            }

            var counted = materialized
                .GroupBy(r => r.ProteinGroup, StringComparer.Ordinal)
                .Select(g => g.First())
                .Where(r => r.QValue <= CounterQValueMax)
                .ToList();

            TsvHeader.Write(output, "go_annotation_format", FormatVersion.ToString(CultureInfo.InvariantCulture));
            TsvHeader.WriteProducer(output);
            TsvHeader.Write(output, "go_release", ontology.Release);
            TsvHeader.Write(output, "go_obo_sha256", ontology.SourceSha256);
            TsvHeader.Write(output, "annotation_db_sha256", annotationDbSha256);
            TsvHeader.Write(output, "source_file_sha256", sourceFileSha256);
            TsvHeader.Write(output, "counter_q_value_max", CounterQValueMax.ToString("R", CultureInfo.InvariantCulture));
            TsvHeader.Write(output, "n_multi_member_groups", Count(counted.Count(r => r.NMembers > 1)));
            foreach (var status in new[] { GoAnnotationStatus.Annotated, GoAnnotationStatus.NoGoTerms,
                         GoAnnotationStatus.NoEntry, GoAnnotationStatus.Contaminant })
            {
                TsvHeader.Write(output, "status_" + StatusName(status), Count(counted.Count(r => r.Status == status)));
            }

            output.Write(TsvWriter.HeaderLine(Schema) + "\n");
            foreach (var row in materialized)
            {
                RejectSetSeparators(row.AccessionUsed, "accession_used", row);
                RejectSetSeparators(row.Evidence, "evidence", row);
                foreach (var column in Schema)
                {
                    TsvHeader.RejectSeparators(column.GetValue(row), column.Header);
                }
                output.Write(TsvWriter.RowLine(Schema, row) + "\n");
            }
        }

        /// <summary>"annotated", "no_go_terms", "no_entry" or "contaminant".</summary>
        public static string StatusName(GoAnnotationStatus status) => status switch
        {
            GoAnnotationStatus.Annotated => "annotated",
            GoAnnotationStatus.NoGoTerms => "no_go_terms",
            GoAnnotationStatus.NoEntry => "no_entry",
            GoAnnotationStatus.Contaminant => "contaminant",
            _ => throw new ArgumentOutOfRangeException(nameof(status), status, null)
        };

        /// <summary>GO's own namespace names; "unknown" when the ontology did not say; null for no term.</summary>
        private static string AspectName(GoAspect? aspect) => aspect switch
        {
            null => null,
            GoAspect.BiologicalProcess => "biological_process",
            GoAspect.CellularComponent => "cellular_component",
            GoAspect.MolecularFunction => "molecular_function",
            GoAspect.Unknown => "unknown",
            _ => throw new ArgumentOutOfRangeException(nameof(aspect), aspect, null)
        };

        private static string Bool(bool? value) => value switch
        {
            true => "true",
            false => "false",
            null => null
        };

        private static string JoinSet(IReadOnlyList<string> members) => members == null ? null : string.Join(';', members);

        private static string Count(int n) => n.ToString(CultureInfo.InvariantCulture);

        private static void RequireSame(string actual, string expected, string field, GoAnnotationRow row)
        {
            if (!string.Equals(actual, expected, StringComparison.Ordinal))
            {
                throw new ArgumentException(
                    $"Row for group '{row.ProteinGroup}' has {field} '{actual}', but the file's is '{expected}'.");
            }
        }

        private static void RejectSetSeparators(IReadOnlyList<string> members, string field, GoAnnotationRow row)
        {
            if (members != null && members.Any(m => m != null && m.Contains(';')))
            {
                throw new ArgumentException($"A {field} member of group '{row.ProteinGroup}' contains ';', the set separator.");
            }
        }
    }

    /// <summary>The "#!key value" header lines the GO tables share.</summary>
    internal static class TsvHeader
    {
        /// <summary>Writes "#!key value"; a null value omits the line.</summary>
        public static void Write(TextWriter output, string key, string value)
        {
            if (value == null)
            {
                return;
            }
            RejectSeparators(value, key);
            output.Write($"#!{key} {value}" + "\n");
        }

        /// <summary>
        /// Writes "#!mzlib_version" (this assembly's informational version, commit included) and
        /// "#!mzlib_release" (the released version, or "none" for a local build), so an ingester can refuse
        /// a file no release produced without parsing the version itself.
        /// </summary>
        public static void WriteProducer(TextWriter output)
        {
            string version = typeof(TsvHeader).Assembly
                .GetCustomAttribute<AssemblyInformationalVersionAttribute>()?.InformationalVersion;
            Write(output, "mzlib_version", version);
            Write(output, "mzlib_release", ReleaseOf(version));
        }

        /// <summary>
        /// The release an informational version names, or "none". release.yml builds with
        /// /p:Version=(git tag); a build without it gets the SDK default 1.0.0, which no release carries.
        /// </summary>
        public static string ReleaseOf(string informationalVersion)
        {
            string version = informationalVersion?.Split('+')[0];
            return string.IsNullOrEmpty(version) || version == "1.0.0" ? "none" : version;
        }

        public static void RejectSeparators(string value, string field)
        {
            if (value != null && value.IndexOfAny(new[] { '\t', '\n', '\r' }) >= 0)
            {
                throw new ArgumentException($"The value of '{field}' contains a tab or line break, which would shift the row.");
            }
        }
    }
}
