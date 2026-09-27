using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using System.Threading;
using System.Threading.Tasks;
using Parquet;
using Parquet.Schema;
using UsefulProteomicsDatabases.Ensembl;

namespace OrthologyStore
{
    /// <summary>One file of a written snapshot. Rows is null for a file that is not a table.</summary>
    public sealed record SnapshotFile(string Path, long? Rows, long Bytes, string Sha256);

    /// <summary>What a write produced: the snapshot's id and every file but the manifest itself.</summary>
    public sealed record SnapshotManifest(string SnapshotId, IReadOnlyList<SnapshotFile> Files);

    /// <summary>
    /// Writes an <see cref="OrthologySnapshot"/> to a directory of Parquet tables that can be queried
    /// with no server (DuckDB, Arrow, pandas), plus views.sql and a manifest:
    ///
    ///   manifest.json                    format, source, release, species, every input and file with sha256
    ///   genes/&lt;species&gt;.parquet        the primary-assembly gene set, one row per gene
    ///   members/&lt;species&gt;.parquet      gene -> Compara gene tree, one row per treed gene
    ///   pairs/&lt;a&gt;__&lt;b&gt;.parquet         every relationship between a and b (a &lt;= b; a == b holds paralogs)
    ///   views.sql                        orthologs, pair_status and species_set3 as DuckDB table macros
    ///
    /// Values are Compara's, verbatim, and a NULL stays null. The same snapshot always gives the same
    /// bytes, so a snapshot is identified by its manifest and can be compared with a later one rather
    /// than overwritten: the writer refuses a directory that is not empty.
    /// </summary>
    public static class OrthologySnapshotWriter
    {
        public const string FormatName = "ensembl-orthology-snapshot";
        public const int FormatVersion = 1;

        private static readonly ParquetSchema GenesSchema = new(
            new DataField<string>("gene_id"), new DataField<int?>("gene_version"), new DataField<string>("gene_biotype"),
            new DataField<string>("gene_name"), new DataField<string>("seq_region"), new DataField<string>("species"));

        private static readonly ParquetSchema MembersSchema = new(
            new DataField<string>("gene_id"), new DataField<string>("species"), new DataField<string>("group_id"),
            new DataField<string>("source_group_id"), new DataField<string>("canonical_protein_id"));

        private static readonly ParquetSchema PairsSchema = new(
            new DataField<string>("homology_id"), new DataField<string>("relationship_type"),
            new DataField<string>("relationship_class"), new DataField<string>("species_a"), new DataField<string>("species_b"),
            new DataField<string>("gene_a"), new DataField<string>("protein_a"), new DataField<double?>("identity_a"),
            new DataField<string>("gene_b"), new DataField<string>("protein_b"), new DataField<double?>("identity_b"),
            new DataField<double?>("dn"), new DataField<double?>("ds"), new DataField<int?>("goc_score"),
            new DataField<double?>("wga_coverage"), new DataField<bool?>("is_high_confidence"),
            new DataField<string>("source_dump"));

        /// <summary>The column names of a table, in order: "genes", "members" or "pairs".</summary>
        public static IReadOnlyList<string> Columns(string table) => (table switch
        {
            "genes" => GenesSchema,
            "members" => MembersSchema,
            "pairs" => PairsSchema,
            _ => throw new ArgumentOutOfRangeException(nameof(table), table, "Not a snapshot table.")
        }).GetDataFields().Select(f => f.Name).ToList();

        /// <summary>
        /// The snapshot's id: a sha256 over the format, release, species and every input's sha256. Two
        /// snapshots with the same id were built from the same bytes by the same format.
        /// </summary>
        public static string SnapshotId(OrthologySnapshot snapshot)
        {
            ArgumentNullException.ThrowIfNull(snapshot);
            string text = string.Join("\n", new[] { $"{FormatName} {FormatVersion}", snapshot.Release }
                .Concat(snapshot.Species)
                .Concat(Inputs(snapshot).Select(i => i.Sha256).OrderBy(s => s, StringComparer.Ordinal)));
            return Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(text))).ToLowerInvariant();
        }

        /// <exception cref="IOException">The directory exists and is not empty.</exception>
        public static async Task<SnapshotManifest> WriteAsync(OrthologySnapshot snapshot, string directory,
            CancellationToken cancellationToken = default)
        {
            ArgumentNullException.ThrowIfNull(snapshot);
            ArgumentException.ThrowIfNullOrEmpty(directory);
            if (Directory.Exists(directory) && Directory.EnumerateFileSystemEntries(directory).Any())
            {
                throw new IOException($"{directory} is not empty. A snapshot is written once, never updated in place.");
            }

            string id = SnapshotId(snapshot);
            var metadata = new Dictionary<string, string>
            {
                ["orthology.format"] = $"{FormatName} {FormatVersion}",
                ["orthology.snapshot_id"] = id,
                ["orthology.source"] = $"ensembl_compara {snapshot.Release}",
            };
            var files = new List<SnapshotFile>();

            foreach (string species in snapshot.Species)
            {
                var genes = snapshot.GeneSets[species].Genes.ToList();
                files.Add(await WriteTable(directory, $"genes/{species}.parquet", GenesSchema, metadata, genes.Count, async rg =>
                {
                    await Strings(rg, GenesSchema, 0, genes.Select(g => g.GeneId));
                    await Values(rg, GenesSchema, 1, genes.Select(g => g.Version), cancellationToken);
                    await Strings(rg, GenesSchema, 2, genes.Select(g => g.Biotype));
                    await Strings(rg, GenesSchema, 3, genes.Select(g => g.Symbol));
                    await Strings(rg, GenesSchema, 4, genes.Select(g => g.SeqRegion));
                    await Strings(rg, GenesSchema, 5, genes.Select(_ => species));
                }, cancellationToken));

                var members = genes.Select(g => snapshot.GeneTrees.TryGetMember(g.GeneId, out var m) ? m : null)
                    .Where(m => m != null).ToList();
                files.Add(await WriteTable(directory, $"members/{species}.parquet", MembersSchema, metadata, members.Count, async rg =>
                {
                    await Strings(rg, MembersSchema, 0, members.Select(m => m.GeneId));
                    await Strings(rg, MembersSchema, 1, members.Select(_ => species));
                    await Strings(rg, MembersSchema, 2, members.Select(m => $"compara-{snapshot.Release}:{m.TreeId}"));
                    await Strings(rg, MembersSchema, 3, members.Select(m => m.TreeId));
                    await Strings(rg, MembersSchema, 4, members.Select(m => m.CanonicalProteinId));
                }, cancellationToken));
            }

            foreach (var (a, b) in snapshot.Pairs)
            {
                var rows = snapshot.Pair(a, b);
                files.Add(await WriteTable(directory, $"pairs/{a}__{b}.parquet", PairsSchema, metadata, rows.Count, async rg =>
                {
                    await Strings(rg, PairsSchema, 0, rows.Select(r => r.HomologyId));
                    await Strings(rg, PairsSchema, 1, rows.Select(r => r.HomologyType));
                    await Strings(rg, PairsSchema, 2, rows.Select(r => ComparaHomologyDump.ClassName(r.Class)));
                    await Strings(rg, PairsSchema, 3, rows.Select(r => r.SpeciesA));
                    await Strings(rg, PairsSchema, 4, rows.Select(r => r.SpeciesB));
                    await Strings(rg, PairsSchema, 5, rows.Select(r => r.GeneA));
                    await Strings(rg, PairsSchema, 6, rows.Select(r => r.ProteinA));
                    await Values(rg, PairsSchema, 7, rows.Select(r => r.IdentityA), cancellationToken);
                    await Strings(rg, PairsSchema, 8, rows.Select(r => r.GeneB));
                    await Strings(rg, PairsSchema, 9, rows.Select(r => r.ProteinB));
                    await Values(rg, PairsSchema, 10, rows.Select(r => r.IdentityB), cancellationToken);
                    await Values(rg, PairsSchema, 11, rows.Select(r => r.Dn), cancellationToken);
                    await Values(rg, PairsSchema, 12, rows.Select(r => r.Ds), cancellationToken);
                    await Values(rg, PairsSchema, 13, rows.Select(r => r.GocScore), cancellationToken);
                    await Values(rg, PairsSchema, 14, rows.Select(r => r.WgaCoverage), cancellationToken);
                    await Values(rg, PairsSchema, 15, rows.Select(r => r.IsHighConfidence), cancellationToken);
                    await Strings(rg, PairsSchema, 16, rows.Select(r => r.SourceFile));
                }, cancellationToken));
            }

            files.Add(WriteText(directory, "views.sql", ViewsSql()));
            files.Sort((x, y) => string.CompareOrdinal(x.Path, y.Path));
            WriteManifest(directory, snapshot, id, files);
            return new SnapshotManifest(id, files);
        }

        /// <summary>The DuckDB table macros written beside every snapshot.</summary>
        public static string ViewsSql()
        {
            using var stream = typeof(OrthologySnapshotWriter).Assembly.GetManifestResourceStream("OrthologyStore.views.sql")
                ?? throw new InvalidOperationException("views.sql is not embedded in the assembly.");
            using var reader = new StreamReader(stream);
            return reader.ReadToEnd().Replace("\r\n", "\n");
        }

        private sealed record Input(string Role, string Species, string File, string Sha256);

        private static IEnumerable<Input> Inputs(OrthologySnapshot snapshot) =>
            snapshot.Dumps.Select(d => new Input("homologies", d.Genome, d.SourceFileName, d.SourceSha256))
                .Append(new Input("gene_tree_content", null, snapshot.GeneTrees.SourceFileName, snapshot.GeneTrees.SourceSha256))
                .Concat(snapshot.Species.Select(s =>
                    new Input("gene_set", s, snapshot.GeneSets[s].SourceFileName, snapshot.GeneSets[s].SourceSha256)));

        private static async Task<SnapshotFile> WriteTable(string directory, string relative, ParquetSchema schema,
            Dictionary<string, string> metadata, int rows, Func<ParquetRowGroupWriter, Task> columns,
            CancellationToken cancellationToken)
        {
            string path = Path.Combine(directory, relative);
            Directory.CreateDirectory(Path.GetDirectoryName(path));
            await using (var file = File.Create(path))
            {
                var options = new ParquetOptions { CompressionMethod = CompressionMethod.Zstd };
                await using var writer = await ParquetWriter.CreateAsync(schema, file, options, false, cancellationToken);
                writer.CustomMetadata = metadata;
                using var rowGroup = writer.CreateRowGroup();
                if (rows > 0)
                {
                    await columns(rowGroup);
                }
            }
            return Describe(directory, relative, rows);
        }

        private static Task Strings(ParquetRowGroupWriter rg, ParquetSchema schema, int column, IEnumerable<string> values) =>
            rg.WriteAsync(schema.GetDataFields()[column], values.ToList());

        private static Task Values<T>(ParquetRowGroupWriter rg, ParquetSchema schema, int column, IEnumerable<T?> values,
            CancellationToken cancellationToken) where T : struct =>
            rg.WriteAsync<T>(schema.GetDataFields()[column], new ReadOnlyMemory<T?>(values.ToArray()), null, null, cancellationToken);

        private static SnapshotFile WriteText(string directory, string relative, string text)
        {
            File.WriteAllBytes(Path.Combine(directory, relative), new UTF8Encoding(false).GetBytes(text));
            return Describe(directory, relative, null);
        }

        private static SnapshotFile Describe(string directory, string relative, long? rows)
        {
            byte[] bytes = File.ReadAllBytes(Path.Combine(directory, relative));
            return new SnapshotFile(relative, rows, bytes.Length, Convert.ToHexString(SHA256.HashData(bytes)).ToLowerInvariant());
        }

        private static void WriteManifest(string directory, OrthologySnapshot snapshot, string id, IReadOnlyList<SnapshotFile> files)
        {
            using var stream = new MemoryStream();
            using (var json = new Utf8JsonWriter(stream, new JsonWriterOptions { Indented = true, NewLine = "\n" }))
            {
                json.WriteStartObject();
                json.WriteString("format", FormatName);
                json.WriteNumber("format_version", FormatVersion);
                json.WriteString("snapshot_id", id);
                json.WriteString("source", "ensembl_compara");
                json.WriteString("release", snapshot.Release);
                json.WriteString("collection", snapshot.GeneTrees.Collection);
                json.WriteStartArray("species");
                foreach (string s in snapshot.Species) json.WriteStringValue(s);
                json.WriteEndArray();

                // What a gene-resolution row carries as gene_set_sha256, so a consumer can check that a
                // resolution and this snapshot were counted against the same gene set.
                json.WriteStartObject("gene_set_sha256");
                foreach (string s in snapshot.Species) json.WriteString(s, snapshot.GeneSets[s].SourceSha256);
                json.WriteEndObject();

                json.WriteStartArray("inputs");
                foreach (var input in Inputs(snapshot))
                {
                    json.WriteStartObject();
                    json.WriteString("role", input.Role);
                    if (input.Species != null) json.WriteString("species", input.Species);
                    json.WriteString("file", input.File);
                    json.WriteString("sha256", input.Sha256);
                    json.WriteEndObject();
                }
                json.WriteEndArray();

                json.WriteStartObject("build_checks");
                json.WriteNumber("identical_duplicate_rows_dropped", snapshot.IdenticalDuplicatesDropped);
                json.WriteEndObject();

                json.WriteStartArray("files");
                foreach (var f in files)
                {
                    json.WriteStartObject();
                    json.WriteString("path", f.Path);
                    if (f.Rows.HasValue) json.WriteNumber("rows", f.Rows.Value);
                    json.WriteNumber("bytes", f.Bytes);
                    json.WriteString("sha256", f.Sha256);
                    json.WriteEndObject();
                }
                json.WriteEndArray();
                json.WriteEndObject();
            }
            stream.WriteByte((byte)'\n');
            File.WriteAllBytes(Path.Combine(directory, "manifest.json"), stream.ToArray());
        }
    }
}
