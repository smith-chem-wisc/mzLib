using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Omics.SpectrumMatch;
using Parquet;
using Parquet.Data;
using Parquet.Schema;

namespace Readers.SpectralLibrary.DiaNNSpectralLibrary
{
    /// <summary>
    /// Reads and writes DIA-NN spectral libraries in Apache Parquet format (.parquet).
    ///
    /// DIA-NN 1.9+ uses Parquet as its primary empirical library format.  It is the preferred
    /// interchange format over TSV because it is:
    ///   • ~10x smaller on disk (columnar + Snappy compression)
    ///   • Strongly typed (no parsing ambiguity)
    ///   • Directly readable from Python (pandas), R (arrow), and Spark
    ///
    /// ⚠️  CRITICAL TYPE REQUIREMENTS (DIA-NN will reject the library if violated):
    ///   • All floating-point columns MUST be FLOAT (single precision, 32-bit).
    ///     Using DOUBLE (64-bit) causes DIA-NN to abort with a schema error.
    ///   • All integer/boolean columns MUST be INT64 (64-bit).
    ///     Using INT32 causes the same rejection.
    ///   • String columns use BYTE_ARRAY with UTF-8 logical type.
    ///
    /// Row structure: one row per fragment ion.  Precursor-level fields are repeated on every
    /// fragment row belonging to that precursor (denormalized / wide format).
    ///
    /// NuGet dependency:  Parquet.Net (Aloneguid)
    ///   dotnet add package Parquet.Net
    /// </summary>
    public static class DiaNNParquetSpectralLibrary
    {
        // ─────────────────────────────────────────────────────────────────────────────────────
        // Column name constants
        //
        // These match DIA-NN's canonical column names.  The reader additionally accepts the
        // aliases listed in ColumnAliases below.
        // ─────────────────────────────────────────────────────────────────────────────────────

        private static class Columns
        {
            // Precursor-level
            public const string ModifiedPeptide   = "ModifiedPeptide";
            public const string StrippedPeptide   = "StrippedPeptide";
            public const string PrecursorCharge   = "PrecursorCharge";
            public const string PrecursorMz       = "PrecursorMz";
            public const string Tr_recalibrated   = "Tr_recalibrated";
            public const string IonMobility       = "IonMobility";
            public const string ProteinId         = "ProteinId";
            public const string ProteinName       = "ProteinName";
            public const string Genes             = "Genes";
            public const string Proteotypic       = "Proteotypic";
            public const string Decoy             = "Decoy";

            // Fragment-level
            public const string FragmentMz        = "FragmentMz";
            public const string RelativeIntensity = "RelativeIntensity";
            public const string FragmentType      = "FragmentType";
            public const string FragmentNumber    = "FragmentNumber";
            public const string FragmentCharge    = "FragmentCharge";
            public const string FragmentLossType  = "FragmentLossType";
            public const string ExcludeFromAssay  = "ExcludeFromAssay";
        }

        /// <summary>
        /// Column aliases accepted on read.  Maps alias → canonical Columns.* name.
        /// Covers Spectronaut-style and OpenSWATH-style variants.
        /// </summary>
        private static readonly Dictionary<string, string> ColumnAliases =
            new(StringComparer.OrdinalIgnoreCase)
            {
                // ModifiedPeptide
                { "ModifiedPeptideSequence",    Columns.ModifiedPeptide   },
                { "FullUniModPeptideName",       Columns.ModifiedPeptide   },
                { "LabeledSequence",             Columns.ModifiedPeptide   },

                // PrecursorCharge
                { "Charge",                      Columns.PrecursorCharge   },
                { "PrecursorZ",                  Columns.PrecursorCharge   },

                // Tr_recalibrated
                { "iRT",                         Columns.Tr_recalibrated   },
                { "RetentionTime",               Columns.Tr_recalibrated   },
                { "NormalizedRetentionTime",     Columns.Tr_recalibrated   },
                { "RT",                          Columns.Tr_recalibrated   },

                // RelativeIntensity
                { "LibraryIntensity",            Columns.RelativeIntensity },
                { "Intensity",                   Columns.RelativeIntensity },

                // ProteinId
                { "UniprotID",                   Columns.ProteinId         },
                { "ProteinAccession",            Columns.ProteinId         },

                // ProteinName
                { "Protein",                     Columns.ProteinName       },
                { "ProteinDescription",          Columns.ProteinName       },
            };

        /// <summary>Required columns — read fails with a clear message if any are absent.</summary>
        private static readonly HashSet<string> RequiredColumns = new(StringComparer.OrdinalIgnoreCase)
        {
            Columns.ModifiedPeptide,
            Columns.PrecursorCharge,
            Columns.PrecursorMz,
            Columns.Tr_recalibrated,
            Columns.FragmentMz,
            Columns.RelativeIntensity,
            Columns.FragmentType,
            Columns.FragmentNumber,
            Columns.FragmentCharge,
        };

        // ─────────────────────────────────────────────────────────────────────────────────────
        // Schema definition
        // ─────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Returns the canonical DIA-NN Parquet schema for a spectral library.
        ///
        /// Type constraints enforced here (DIA-NN is strict):
        ///   • floats  → DataType.Float    (FLOAT / single-precision)
        ///   • ints    → DataType.Int64    (INT64)
        ///   • strings → DataType.String   (BYTE_ARRAY / UTF-8)
        /// </summary>
        public static ParquetSchema GetSchema() => new ParquetSchema(
            // ── Precursor columns ──────────────────────────────────────────────────────────
            new DataField<string>(Columns.ModifiedPeptide),
            new DataField<string>(Columns.StrippedPeptide),
            new DataField<long>(Columns.PrecursorCharge),         // INT64, not INT32
            new DataField<float>(Columns.PrecursorMz),            // FLOAT, not DOUBLE
            new DataField<float>(Columns.Tr_recalibrated),        // FLOAT, not DOUBLE
            new DataField<float>(Columns.IonMobility),            // FLOAT, not DOUBLE
            new DataField<string>(Columns.ProteinId),
            new DataField<string>(Columns.ProteinName),
            new DataField<string>(Columns.Genes),
            new DataField<long>(Columns.Proteotypic),             // INT64
            new DataField<long>(Columns.Decoy),                   // INT64

            // ── Fragment columns ───────────────────────────────────────────────────────────
            new DataField<float>(Columns.FragmentMz),             // FLOAT, not DOUBLE
            new DataField<float>(Columns.RelativeIntensity),      // FLOAT, not DOUBLE
            new DataField<string>(Columns.FragmentType),
            new DataField<long>(Columns.FragmentNumber),          // INT64
            new DataField<long>(Columns.FragmentCharge),          // INT64
            new DataField<string>(Columns.FragmentLossType),
            new DataField<long>(Columns.ExcludeFromAssay)         // INT64
        );

        // ─────────────────────────────────────────────────────────────────────────────────────
        // Write
        // ─────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Writes a list of <see cref="DiaNNLibraryEntry"/> objects to a Parquet file.
        ///
        /// The output is written in a single row group (simplest and fastest for libraries
        /// up to ~100M rows on modern hardware).  For libraries with billions of rows,
        /// use <see cref="WriteInRowGroups"/> to split into manageable chunks.
        /// </summary>
        /// <param name="filePath">Output .parquet file path.</param>
        /// <param name="entries">Library entries to write.</param>
        /// <param name="compressionMethod">
        /// Compression codec.  Snappy (default) gives the best read speed.
        /// Gzip gives ~20% better compression at the cost of ~5x slower write.
        /// </param>
        public static void Write(
            string filePath,
            IReadOnlyList<DiaNNLibraryEntry> entries,
            CompressionMethod compressionMethod = CompressionMethod.Snappy)
        {
            if (entries == null) throw new ArgumentNullException(nameof(entries));
            if (entries.Count == 0) throw new ArgumentException("Cannot write an empty spectral library.", nameof(entries));

            // Expand entries → flat fragment rows
            var rows = ExpandToRows(entries);
            WriteRows(filePath, rows, compressionMethod, append: false);
        }

        /// <summary>
        /// Convenience overload that converts mzLib <see cref="LibrarySpectrum"/> objects first.
        /// </summary>
        public static void WriteFromLibrarySpectra(
            string filePath,
            IReadOnlyList<LibrarySpectrum> spectra,
            CompressionMethod compressionMethod = CompressionMethod.Snappy)
        {
            var entries = spectra.Select(DiaNNLibraryEntry.FromLibrarySpectrum).ToList();
            Write(filePath, entries, compressionMethod);
        }

        /// <summary>
        /// Writes entries in multiple row groups, each containing at most
        /// <paramref name="rowsPerGroup"/> fragment rows.
        ///
        /// Prefer this for very large libraries (>50M fragment rows) to avoid allocating
        /// the entire column arrays in memory at once.
        /// </summary>
        public static void WriteInRowGroups(
            string filePath,
            IReadOnlyList<DiaNNLibraryEntry> entries,
            int rowsPerGroup = 5_000_000,
            CompressionMethod compressionMethod = CompressionMethod.Snappy)
        {
            if (entries == null) throw new ArgumentNullException(nameof(entries));
            if (rowsPerGroup <= 0) throw new ArgumentOutOfRangeException(nameof(rowsPerGroup));

            var rows = ExpandToRows(entries);
            bool firstGroup = true;

            for (int offset = 0; offset < rows.Count; offset += rowsPerGroup)
            {
                var chunk = rows.Skip(offset).Take(rowsPerGroup).ToList();
                WriteRows(filePath, chunk, compressionMethod, append: !firstGroup);
                firstGroup = false;
            }
        }

        // ─────────────────────────────────────────────────────────────────────────────────────
        // Read
        // ─────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Reads a DIA-NN Parquet spectral library from disk.
        /// Groups fragment rows by precursor key (ModifiedPeptide/PrecursorCharge).
        /// </summary>
        /// <param name="filePath">Path to the .parquet file.</param>
        /// <returns>List of <see cref="DiaNNLibraryEntry"/> objects, one per unique precursor.</returns>
        /// <exception cref="InvalidDataException">
        /// Thrown if required columns are missing from the file schema.
        /// </exception>
        public static List<DiaNNLibraryEntry> Read(string filePath)
        {
            if (!File.Exists(filePath))
                throw new FileNotFoundException($"Parquet library not found: {filePath}");

            using var stream = File.OpenRead(filePath);
            return ReadFromStream(stream, filePath);
        }

        /// <summary>
        /// Reads a Parquet library from a stream (useful for testing with MemoryStream).
        /// </summary>
        public static List<DiaNNLibraryEntry> ReadFromStream(Stream stream, string? diagnosticPath = null)
        {
            using var reader = new ParquetReader(stream);
            return ReadFromParquetReader(reader, diagnosticPath ?? "<stream>");
        }

        /// <summary>
        /// Reads a Parquet library and returns mzLib <see cref="LibrarySpectrum"/> objects directly.
        /// </summary>
        public static List<LibrarySpectrum> ReadAsLibrarySpectra(string filePath) =>
            Read(filePath).Select(e => e.ToLibrarySpectrum()).ToList();

        // ─────────────────────────────────────────────────────────────────────────────────────
        // Schema inspection
        // ─────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Inspects the schema of a Parquet file without loading its data.
        /// Returns a summary string with column names and types — useful for diagnosing
        /// DIA-NN type mismatch errors.
        /// </summary>
        public static string InspectSchema(string filePath)
        {
            using var stream = File.OpenRead(filePath);
            using var reader = new ParquetReader(stream);

            var sb = new System.Text.StringBuilder();
            sb.AppendLine($"File: {filePath}");
            sb.AppendLine($"Row groups: {reader.RowGroupCount}");
            sb.AppendLine("Schema:");

            foreach (var field in reader.Schema.Fields)
            {
                string required = RequiredColumns.Contains(field.Name) ? " [required]" : "";
                sb.AppendLine($"  {field.Name}  ({field.GetType().Name}){required}");
            }

            return sb.ToString();
        }

        /// <summary>
        /// Validates that a Parquet file's column types conform to DIA-NN's strict requirements.
        /// Returns a list of type violations (empty = valid for DIA-NN).
        /// </summary>
        public static List<string> ValidateSchema(string filePath)
        {
            var violations = new List<string>();
            using var stream = File.OpenRead(filePath);
            using var reader = new ParquetReader(stream);

            var floatColumns = new HashSet<string>(StringComparer.OrdinalIgnoreCase)
            {
                Columns.PrecursorMz, Columns.Tr_recalibrated, Columns.IonMobility,
                Columns.FragmentMz, Columns.RelativeIntensity
            };
            var int64Columns = new HashSet<string>(StringComparer.OrdinalIgnoreCase)
            {
                Columns.PrecursorCharge, Columns.FragmentNumber, Columns.FragmentCharge,
                Columns.Proteotypic, Columns.Decoy, Columns.ExcludeFromAssay
            };

            var schemaFieldsByName = reader.Schema.Fields
                .OfType<DataField>()
                .ToDictionary(f => f.Name, f => (Field)f, StringComparer.OrdinalIgnoreCase);

            foreach (var col in floatColumns)
            {
                if (!schemaFieldsByName.TryGetValue(col, out var field)) continue;
                if (field is not DataField<float>)
                    violations.Add($"Column '{col}' must be FLOAT (single-precision) but is {field.GetType().Name}. DIA-NN rejects DOUBLE.");
            }

            foreach (var col in int64Columns)
            {
                if (!schemaFieldsByName.TryGetValue(col, out var field)) continue;
                if (field is not DataField<long>)
                    violations.Add($"Column '{col}' must be INT64 but is {field.GetType().Name}. DIA-NN rejects INT32.");
            }

            return violations;
        }

        // ─────────────────────────────────────────────────────────────────────────────────────
        // Internal: flat row model
        // ─────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Intermediate flat representation of one Parquet row (one fragment ion).
        /// Avoids re-allocating strings for every column during columnar assembly.
        /// </summary>
        private sealed class ParquetFragmentRow
        {
            // Precursor-level (repeated on every fragment row)
            public string ModifiedPeptide   = string.Empty;
            public string StrippedPeptide   = string.Empty;
            public long   PrecursorCharge;
            public float  PrecursorMz;
            public float  Tr_recalibrated;
            public float  IonMobility;
            public string ProteinId         = string.Empty;
            public string ProteinName       = string.Empty;
            public string Genes             = string.Empty;
            public long   Proteotypic;
            public long   Decoy;

            // Fragment-level
            public float  FragmentMz;
            public float  RelativeIntensity;
            public string FragmentType      = string.Empty;
            public long   FragmentNumber;
            public long   FragmentCharge;
            public string FragmentLossType  = "noloss";
            public long   ExcludeFromAssay;
        }

        /// <summary>
        /// Expands a list of library entries to a flat list of per-fragment rows.
        /// Fragments within each entry are emitted in m/z order (ascending).
        /// </summary>
        private static List<ParquetFragmentRow> ExpandToRows(IReadOnlyList<DiaNNLibraryEntry> entries)
        {
            int totalRows = entries.Sum(e => e.Fragments?.Count ?? 0);
            var rows = new List<ParquetFragmentRow>(totalRows);

            foreach (var entry in entries)
            {
                if (entry.Fragments == null || entry.Fragments.Count == 0) continue;

                string modSeq      = entry.ModifiedSequence ?? string.Empty;
                string stripped    = entry.StrippedSequence  ?? DiaNNModificationMapping.GetStrippedSequence(modSeq);
                long   charge      = entry.PrecursorCharge;
                float  mz          = (float)entry.PrecursorMz;
                float  rt          = (float)entry.RetentionTime;
                float  im          = (float)entry.IonMobility;
                string proteinId   = entry.ProteinId   ?? string.Empty;
                string proteinName = entry.ProteinName ?? string.Empty;
                string genes       = entry.GeneName    ?? string.Empty;
                long   proteotypic = entry.IsProteotypic ? 1L : 0L;
                long   decoy       = entry.IsDecoy       ? 1L : 0L;

                // Sort fragments by m/z before writing (good practice; DIA-NN may expect this)
                var sortedFrags = entry.Fragments.OrderBy(f => f.Mz);

                foreach (var frag in sortedFrags)
                {
                    rows.Add(new ParquetFragmentRow
                    {
                        ModifiedPeptide   = modSeq,
                        StrippedPeptide   = stripped,
                        PrecursorCharge   = charge,
                        PrecursorMz       = mz,
                        Tr_recalibrated   = rt,
                        IonMobility       = im,
                        ProteinId         = proteinId,
                        ProteinName       = proteinName,
                        Genes             = genes,
                        Proteotypic       = proteotypic,
                        Decoy             = decoy,

                        FragmentMz        = (float)frag.Mz,
                        RelativeIntensity = (float)Math.Clamp(frag.Intensity, 0.0, 1.0),
                        FragmentType      = frag.IonType.ToString(),
                        FragmentNumber    = frag.SeriesNumber,
                        FragmentCharge    = frag.Charge,
                        FragmentLossType  = frag.LossType ?? "noloss",
                        ExcludeFromAssay  = 0L
                    });
                }
            }

            return rows;
        }

        // ─────────────────────────────────────────────────────────────────────────────────────
        // Internal: columnar write
        // ─────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Writes a flat list of rows as one Parquet row group.
        /// Uses the canonical DIA-NN schema from <see cref="GetSchema"/>.
        /// </summary>
        private static void WriteRows(
            string filePath,
            List<ParquetFragmentRow> rows,
            CompressionMethod compression,
            bool append)
        {
            int n = rows.Count;
            if (n == 0) return;

            // ── Allocate column arrays ─────────────────────────────────────────────────────
            var col_ModifiedPeptide   = new string[n];
            var col_StrippedPeptide   = new string[n];
            var col_PrecursorCharge   = new long[n];
            var col_PrecursorMz       = new float[n];
            var col_Tr_recalibrated   = new float[n];
            var col_IonMobility       = new float[n];
            var col_ProteinId         = new string[n];
            var col_ProteinName       = new string[n];
            var col_Genes             = new string[n];
            var col_Proteotypic       = new long[n];
            var col_Decoy             = new long[n];
            var col_FragmentMz        = new float[n];
            var col_RelativeIntensity = new float[n];
            var col_FragmentType      = new string[n];
            var col_FragmentNumber    = new long[n];
            var col_FragmentCharge    = new long[n];
            var col_FragmentLossType  = new string[n];
            var col_ExcludeFromAssay  = new long[n];

            // ── Fill column arrays ─────────────────────────────────────────────────────────
            for (int i = 0; i < n; i++)
            {
                var r = rows[i];
                col_ModifiedPeptide[i]   = r.ModifiedPeptide;
                col_StrippedPeptide[i]   = r.StrippedPeptide;
                col_PrecursorCharge[i]   = r.PrecursorCharge;
                col_PrecursorMz[i]       = r.PrecursorMz;
                col_Tr_recalibrated[i]   = r.Tr_recalibrated;
                col_IonMobility[i]       = r.IonMobility;
                col_ProteinId[i]         = r.ProteinId;
                col_ProteinName[i]       = r.ProteinName;
                col_Genes[i]             = r.Genes;
                col_Proteotypic[i]       = r.Proteotypic;
                col_Decoy[i]             = r.Decoy;
                col_FragmentMz[i]        = r.FragmentMz;
                col_RelativeIntensity[i] = r.RelativeIntensity;
                col_FragmentType[i]      = r.FragmentType;
                col_FragmentNumber[i]    = r.FragmentNumber;
                col_FragmentCharge[i]    = r.FragmentCharge;
                col_FragmentLossType[i]  = r.FragmentLossType;
                col_ExcludeFromAssay[i]  = r.ExcludeFromAssay;
            }

            // ── Build strongly-typed field references ──────────────────────────────────────
            // We create the schema once and cache the DataField<T> references so that each
            // WriteColumn call uses the correct generic type.  This avoids the need for a
            // string-keyed indexer on ParquetSchema (which Parquet.Net does not provide).
            var f_ModifiedPeptide   = new DataField<string>(Columns.ModifiedPeptide);
            var f_StrippedPeptide   = new DataField<string>(Columns.StrippedPeptide);
            var f_PrecursorCharge   = new DataField<long>  (Columns.PrecursorCharge);
            var f_PrecursorMz       = new DataField<float> (Columns.PrecursorMz);
            var f_Tr_recalibrated   = new DataField<float> (Columns.Tr_recalibrated);
            var f_IonMobility       = new DataField<float> (Columns.IonMobility);
            var f_ProteinId         = new DataField<string>(Columns.ProteinId);
            var f_ProteinName       = new DataField<string>(Columns.ProteinName);
            var f_Genes             = new DataField<string>(Columns.Genes);
            var f_Proteotypic       = new DataField<long>  (Columns.Proteotypic);
            var f_Decoy             = new DataField<long>  (Columns.Decoy);
            var f_FragmentMz        = new DataField<float> (Columns.FragmentMz);
            var f_RelativeIntensity = new DataField<float> (Columns.RelativeIntensity);
            var f_FragmentType      = new DataField<string>(Columns.FragmentType);
            var f_FragmentNumber    = new DataField<long>  (Columns.FragmentNumber);
            var f_FragmentCharge    = new DataField<long>  (Columns.FragmentCharge);
            var f_FragmentLossType  = new DataField<string>(Columns.FragmentLossType);
            var f_ExcludeFromAssay  = new DataField<long>  (Columns.ExcludeFromAssay);

            var schema = new ParquetSchema(
                f_ModifiedPeptide, f_StrippedPeptide, f_PrecursorCharge,
                f_PrecursorMz, f_Tr_recalibrated, f_IonMobility,
                f_ProteinId, f_ProteinName, f_Genes,
                f_Proteotypic, f_Decoy,
                f_FragmentMz, f_RelativeIntensity, f_FragmentType,
                f_FragmentNumber, f_FragmentCharge, f_FragmentLossType, f_ExcludeFromAssay);

            // ── Write ──────────────────────────────────────────────────────────────────────
            var mode = append ? FileMode.Append : FileMode.Create;

            using var fileStream = new FileStream(filePath, mode, FileAccess.Write, FileShare.None);
            using var writer     = new ParquetWriter(schema, fileStream, append: append);
            writer.CompressionMethod = compression;

            using var groupWriter = writer.CreateRowGroup();

            groupWriter.WriteColumn(new DataColumn(f_ModifiedPeptide,   col_ModifiedPeptide));
            groupWriter.WriteColumn(new DataColumn(f_StrippedPeptide,   col_StrippedPeptide));
            groupWriter.WriteColumn(new DataColumn(f_PrecursorCharge,   col_PrecursorCharge));
            groupWriter.WriteColumn(new DataColumn(f_PrecursorMz,       col_PrecursorMz));
            groupWriter.WriteColumn(new DataColumn(f_Tr_recalibrated,   col_Tr_recalibrated));
            groupWriter.WriteColumn(new DataColumn(f_IonMobility,       col_IonMobility));
            groupWriter.WriteColumn(new DataColumn(f_ProteinId,         col_ProteinId));
            groupWriter.WriteColumn(new DataColumn(f_ProteinName,       col_ProteinName));
            groupWriter.WriteColumn(new DataColumn(f_Genes,             col_Genes));
            groupWriter.WriteColumn(new DataColumn(f_Proteotypic,       col_Proteotypic));
            groupWriter.WriteColumn(new DataColumn(f_Decoy,             col_Decoy));
            groupWriter.WriteColumn(new DataColumn(f_FragmentMz,        col_FragmentMz));
            groupWriter.WriteColumn(new DataColumn(f_RelativeIntensity, col_RelativeIntensity));
            groupWriter.WriteColumn(new DataColumn(f_FragmentType,      col_FragmentType));
            groupWriter.WriteColumn(new DataColumn(f_FragmentNumber,    col_FragmentNumber));
            groupWriter.WriteColumn(new DataColumn(f_FragmentCharge,    col_FragmentCharge));
            groupWriter.WriteColumn(new DataColumn(f_FragmentLossType,  col_FragmentLossType));
            groupWriter.WriteColumn(new DataColumn(f_ExcludeFromAssay,  col_ExcludeFromAssay));
        }

        // ─────────────────────────────────────────────────────────────────────────────────────
        // Internal: columnar read
        // ─────────────────────────────────────────────────────────────────────────────────────

        private static List<DiaNNLibraryEntry> ReadFromParquetReader(
            ParquetReader reader,
            string diagnosticPath)
        {
            // ── 1. Resolve column names (including aliases) ────────────────────────────────
            var colMap = BuildColumnMap(reader.Schema, diagnosticPath);

            // ── 2. Read all row groups into flat arrays ────────────────────────────────────
            // Accumulate across all row groups then group precursors.
            // For most DIA-NN libraries this is a single row group.
            var allRows = new List<ParquetFragmentRow>(capacity: 1_000_000);

            for (int rg = 0; rg < reader.RowGroupCount; rg++)
            {
                using var groupReader = reader.OpenRowGroupReader(rg);
                ReadRowGroup(groupReader, colMap, allRows);
            }

            // ── 3. Group flat rows by (ModifiedPeptide, PrecursorCharge) ──────────────────
            return GroupRowsToEntries(allRows);
        }

        /// <summary>
        /// Reads one row group from the Parquet file into flat <see cref="ParquetFragmentRow"/> objects.
        /// Handles column aliases transparently.
        /// </summary>
        private static void ReadRowGroup(
            ParquetRowGroupReader groupReader,
            Dictionary<string, string> colMap,
            List<ParquetFragmentRow> destination)
        {
            // Read all columns we care about.  Missing optional columns fall back to defaults.
            string[] col_ModifiedPeptide   = ReadStringColumn(groupReader,  colMap, Columns.ModifiedPeptide);
            string[] col_StrippedPeptide   = ReadStringColumn(groupReader,  colMap, Columns.StrippedPeptide,  optional: true);
            long[]   col_PrecursorCharge   = ReadInt64Column(groupReader,   colMap, Columns.PrecursorCharge);
            float[]  col_PrecursorMz       = ReadFloatColumn(groupReader,   colMap, Columns.PrecursorMz);
            float[]  col_Tr_recalibrated   = ReadFloatColumn(groupReader,   colMap, Columns.Tr_recalibrated);
            float[]  col_IonMobility       = ReadFloatColumn(groupReader,   colMap, Columns.IonMobility,      optional: true);
            string[] col_ProteinId         = ReadStringColumn(groupReader,  colMap, Columns.ProteinId,        optional: true);
            string[] col_ProteinName       = ReadStringColumn(groupReader,  colMap, Columns.ProteinName,      optional: true);
            string[] col_Genes             = ReadStringColumn(groupReader,  colMap, Columns.Genes,            optional: true);
            long[]   col_Proteotypic       = ReadInt64Column(groupReader,   colMap, Columns.Proteotypic,      optional: true);
            long[]   col_Decoy             = ReadInt64Column(groupReader,   colMap, Columns.Decoy,            optional: true);
            float[]  col_FragmentMz        = ReadFloatColumn(groupReader,   colMap, Columns.FragmentMz);
            float[]  col_RelativeIntensity = ReadFloatColumn(groupReader,   colMap, Columns.RelativeIntensity);
            string[] col_FragmentType      = ReadStringColumn(groupReader,  colMap, Columns.FragmentType);
            long[]   col_FragmentNumber    = ReadInt64Column(groupReader,   colMap, Columns.FragmentNumber);
            long[]   col_FragmentCharge    = ReadInt64Column(groupReader,   colMap, Columns.FragmentCharge);
            string[] col_FragmentLossType  = ReadStringColumn(groupReader,  colMap, Columns.FragmentLossType, optional: true);
            long[]   col_ExcludeFromAssay  = ReadInt64Column(groupReader,   colMap, Columns.ExcludeFromAssay, optional: true);

            int rowCount = col_ModifiedPeptide.Length;

            for (int i = 0; i < rowCount; i++)
            {
                destination.Add(new ParquetFragmentRow
                {
                    ModifiedPeptide   = col_ModifiedPeptide[i],
                    StrippedPeptide   = col_StrippedPeptide.Length > i ? col_StrippedPeptide[i] : string.Empty,
                    PrecursorCharge   = col_PrecursorCharge[i],
                    PrecursorMz       = col_PrecursorMz[i],
                    Tr_recalibrated   = col_Tr_recalibrated[i],
                    IonMobility       = col_IonMobility.Length > i ? col_IonMobility[i] : 0f,
                    ProteinId         = col_ProteinId.Length   > i ? col_ProteinId[i]   : string.Empty,
                    ProteinName       = col_ProteinName.Length > i ? col_ProteinName[i] : string.Empty,
                    Genes             = col_Genes.Length       > i ? col_Genes[i]       : string.Empty,
                    Proteotypic       = col_Proteotypic.Length > i ? col_Proteotypic[i] : 1L,
                    Decoy             = col_Decoy.Length       > i ? col_Decoy[i]       : 0L,

                    FragmentMz        = col_FragmentMz[i],
                    RelativeIntensity = col_RelativeIntensity[i],
                    FragmentType      = col_FragmentType[i],
                    FragmentNumber    = col_FragmentNumber[i],
                    FragmentCharge    = col_FragmentCharge[i],
                    FragmentLossType  = col_FragmentLossType.Length > i ? col_FragmentLossType[i] : "noloss",
                    ExcludeFromAssay  = col_ExcludeFromAssay.Length > i ? col_ExcludeFromAssay[i] : 0L,
                });
            }
        }

        /// <summary>
        /// Groups flat fragment rows by their precursor key and builds <see cref="DiaNNLibraryEntry"/> objects.
        /// Preserves insertion order of first occurrence of each precursor.
        /// </summary>
        private static List<DiaNNLibraryEntry> GroupRowsToEntries(List<ParquetFragmentRow> rows)
        {
            // Use an ordered dictionary (insertion-ordered in .NET 5+) to preserve precursor order
            var precursorOrder = new List<string>();
            var precursorMap   = new Dictionary<string, DiaNNLibraryEntry>(StringComparer.Ordinal);

            foreach (var row in rows)
            {
                // Normalize intensity: DIA-NN sometimes stores 0–10000 scale
                float intensity = row.RelativeIntensity;
                if (intensity > 1.0f) intensity /= 10000.0f;

                string key = row.ModifiedPeptide + "/" + row.PrecursorCharge.ToString();

                if (!precursorMap.TryGetValue(key, out var entry))
                {
                    entry = new DiaNNLibraryEntry
                    {
                        ModifiedSequence = row.ModifiedPeptide,
                        StrippedSequence = row.StrippedPeptide.Length > 0
                            ? row.StrippedPeptide
                            : DiaNNModificationMapping.GetStrippedSequence(row.ModifiedPeptide),
                        PrecursorCharge   = (int)row.PrecursorCharge,
                        PrecursorMz       = row.PrecursorMz,
                        RetentionTime     = row.Tr_recalibrated,
                        IonMobility       = row.IonMobility,
                        ProteinId         = row.ProteinId,
                        ProteinName       = row.ProteinName,
                        GeneName          = row.Genes,
                        IsProteotypic     = row.Proteotypic != 0,
                        IsDecoy           = row.Decoy != 0,
                        Fragments         = new List<DiaNNFragmentIon>()
                    };
                    precursorMap[key] = entry;
                    precursorOrder.Add(key);
                }

                // Only include fragment if it is not flagged for exclusion
                if (row.ExcludeFromAssay == 0)
                {
                    char ionType = row.FragmentType.Length > 0 ? char.ToUpperInvariant(row.FragmentType[0]) : 'y';
                    entry.Fragments.Add(new DiaNNFragmentIon
                    {
                        Mz           = row.FragmentMz,
                        Intensity    = intensity,
                        IonType      = ionType,
                        SeriesNumber = (int)row.FragmentNumber,
                        Charge       = (int)row.FragmentCharge,
                        LossType     = string.IsNullOrWhiteSpace(row.FragmentLossType) ? "noloss" : row.FragmentLossType
                    });
                }
            }

            // Return in original precursor insertion order
            return precursorOrder.Select(k => precursorMap[k]).ToList();
        }

        // ─────────────────────────────────────────────────────────────────────────────────────
        // Column-read helpers
        // ─────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Builds a canonical-name → actual-file-column-name mapping, resolving aliases.
        /// Validates that all required columns are present.
        /// </summary>
        private static Dictionary<string, string> BuildColumnMap(
            ParquetSchema schema,
            string diagnosticPath)
        {
            // Build set of actual file column names
            var fileColumns = new HashSet<string>(
                schema.Fields.Select(f => f.Name),
                StringComparer.OrdinalIgnoreCase);

            var map = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);

            // For each canonical name, check if the file has it directly or via alias
            var allCanonical = new[]
            {
                Columns.ModifiedPeptide,  Columns.StrippedPeptide,   Columns.PrecursorCharge,
                Columns.PrecursorMz,      Columns.Tr_recalibrated,   Columns.IonMobility,
                Columns.ProteinId,        Columns.ProteinName,        Columns.Genes,
                Columns.Proteotypic,      Columns.Decoy,
                Columns.FragmentMz,       Columns.RelativeIntensity,  Columns.FragmentType,
                Columns.FragmentNumber,   Columns.FragmentCharge,     Columns.FragmentLossType,
                Columns.ExcludeFromAssay
            };

            foreach (var canonical in allCanonical)
            {
                if (fileColumns.Contains(canonical))
                {
                    map[canonical] = canonical;
                    continue;
                }

                // Check aliases
                bool found = false;
                foreach (var (alias, target) in ColumnAliases)
                {
                    if (string.Equals(target, canonical, StringComparison.OrdinalIgnoreCase)
                        && fileColumns.Contains(alias))
                    {
                        map[canonical] = alias;
                        found = true;
                        break;
                    }
                }

                if (!found && !RequiredColumns.Contains(canonical))
                {
                    // Optional column absent — map to empty string as sentinel
                    map[canonical] = string.Empty;
                }
            }

            // Validate required columns
            var missing = RequiredColumns
                .Where(req => !map.ContainsKey(req) || map[req] == string.Empty)
                .ToList();

            if (missing.Count > 0)
            {
                throw new InvalidDataException(
                    $"Parquet file '{diagnosticPath}' is missing required columns: {string.Join(", ", missing)}.\n" +
                    $"Available columns: {string.Join(", ", fileColumns)}");
            }

            return map;
        }

        // ── Typed column readers ───────────────────────────────────────────────────────────

        /// <summary>
        /// Looks up a <see cref="DataField"/> in a schema by column name (case-insensitive).
        /// Returns null if the column is not present and it's optional.
        /// </summary>
        private static DataField? FindField(ParquetSchema schema, string columnName)
        {
            return schema.Fields
                .OfType<DataField>()
                .FirstOrDefault(f => string.Equals(f.Name, columnName, StringComparison.OrdinalIgnoreCase));
        }

        private static string[] ReadStringColumn(
            ParquetRowGroupReader groupReader,
            Dictionary<string, string> colMap,
            string canonicalName,
            bool optional = false)
        {
            string? actualName = ResolveColumnName(colMap, canonicalName, optional);
            if (actualName == null) return Array.Empty<string>();

            var field = FindField(groupReader.Schema, actualName);
            if (field == null) return Array.Empty<string>();

            var dataColumn = groupReader.ReadColumn(field);
            return dataColumn.Data as string[] ?? Array.Empty<string>();
        }

        private static float[] ReadFloatColumn(
            ParquetRowGroupReader groupReader,
            Dictionary<string, string> colMap,
            string canonicalName,
            bool optional = false)
        {
            string? actualName = ResolveColumnName(colMap, canonicalName, optional);
            if (actualName == null) return Array.Empty<float>();

            var field = FindField(groupReader.Schema, actualName);
            if (field == null) return Array.Empty<float>();

            var dataColumn = groupReader.ReadColumn(field);

            // Handle both FLOAT and DOUBLE gracefully on read (we enforce FLOAT on write)
            if (dataColumn.Data is float[]  floats)  return floats;
            if (dataColumn.Data is double[] doubles) return doubles.Select(d => (float)d).ToArray();

            return Array.Empty<float>();
        }

        private static long[] ReadInt64Column(
            ParquetRowGroupReader groupReader,
            Dictionary<string, string> colMap,
            string canonicalName,
            bool optional = false)
        {
            string? actualName = ResolveColumnName(colMap, canonicalName, optional);
            if (actualName == null) return Array.Empty<long>();

            var field = FindField(groupReader.Schema, actualName);
            if (field == null) return Array.Empty<long>();

            var dataColumn = groupReader.ReadColumn(field);

            // Handle both INT64 and INT32 gracefully on read (we enforce INT64 on write)
            if (dataColumn.Data is long[] longs) return longs;
            if (dataColumn.Data is int[]  ints)  return ints.Select(x => (long)x).ToArray();

            return Array.Empty<long>();
        }

        private static string? ResolveColumnName(
            Dictionary<string, string> colMap,
            string canonicalName,
            bool optional)
        {
            if (!colMap.TryGetValue(canonicalName, out string? actualName))
                return optional ? null : throw new InvalidOperationException($"Column '{canonicalName}' not in column map.");

            return string.IsNullOrEmpty(actualName) ? null : actualName;
        }
    }
}
