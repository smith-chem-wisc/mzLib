using System;
using System.Buffers;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using Omics.SpectrumMatch;

namespace Readers.SpectralLibrary.DiaNNSpectralLibrary
{
    /// <summary>
    /// Writes spectral libraries to the DIA-NN binary .speclib format.
    ///
    /// Output format overview (matching the reader's expectations):
    ///
    ///   [4  bytes]  int32  version   (written as the SpecLibVersion enum value)
    ///   [4  bytes]  int32  gen_decoys flag
    ///   [4  bytes]  int32  precursor count
    ///   [variable]  precursor records  (count × PrecursorRecord.Size bytes, sequential)
    ///   [variable]  fragment blocks    (one contiguous FragmentRecord array per precursor,
    ///                                   at the file offset stored in PrecursorRecord.FragmentOffset)
    ///   [variable]  string table       (null-terminated UTF-8 strings; precursor records store
    ///                                   byte offsets into this section relative to its start)
    ///
    /// All numeric values are written in little-endian byte order (native on x86/x64, matching DIA-NN).
    ///
    /// SCAFFOLD NOTE:  The exact byte layout of PrecursorRecord and FragmentRecord are estimates
    /// pending validation against real .speclib files (see DiaNNBinaryStructs.cs).  The writer
    /// produces output that is structurally consistent with DiaNNSpecLibReader and will round-trip
    /// correctly through that reader.  Compatibility with the actual DIA-NN executable requires
    /// Prompt 8 binary validation.
    /// </summary>
    public class DiaNNSpecLibWriter : IDisposable
    {
        // -----------------------------------------------------------------------------------------
        // Constants
        // -----------------------------------------------------------------------------------------

        /// <summary>
        /// Version written by default.  Use SpecLibVersion.V3 for broad compatibility with
        /// DIA-NN 1.8.x readers while still including ion mobility fields.
        /// Change to SpecLibVersion.V8 once the v8 struct layout is validated.
        /// </summary>
        public const SpecLibVersion DefaultWriteVersion = SpecLibVersion.V3;

        // -----------------------------------------------------------------------------------------
        // Public API
        // -----------------------------------------------------------------------------------------

        /// <summary>
        /// Writes a list of <see cref="DiaNNLibraryEntry"/> objects to a binary .speclib file.
        /// </summary>
        /// <param name="filePath">Destination file path (will be created or overwritten).</param>
        /// <param name="entries">Library entries to write.  Must not be null or empty.</param>
        /// <param name="version">Binary format version to emit.  Default: V3.</param>
        /// <param name="progress">Optional progress callback — invoked with (written, total) counts.</param>
        /// <exception cref="ArgumentNullException">entries is null.</exception>
        /// <exception cref="ArgumentException">entries is empty.</exception>
        /// <exception cref="IOException">File cannot be created or written.</exception>
        public static void Write(
            string filePath,
            IReadOnlyList<DiaNNLibraryEntry> entries,
            SpecLibVersion version = DefaultWriteVersion,
            Action<int, int>? progress = null)
        {
            if (entries == null) throw new ArgumentNullException(nameof(entries));
            if (entries.Count == 0) throw new ArgumentException("Cannot write an empty spectral library.", nameof(entries));

            using var stream = new FileStream(filePath, FileMode.Create, FileAccess.Write, FileShare.None,
                bufferSize: 1 << 20, // 1 MB write buffer
                useAsync: false);
            using var writer = new BinaryWriter(stream, Encoding.UTF8, leaveOpen: false);

            WriteInternal(writer, entries, version, progress);
        }

        /// <summary>
        /// Convenience overload — converts <see cref="LibrarySpectrum"/> objects to
        /// <see cref="DiaNNLibraryEntry"/> wrappers and writes them.
        /// Useful for exporting MetaMorpheus-generated libraries to DIA-NN format.
        /// </summary>
        public static void WriteFromLibrarySpectra(
            string filePath,
            IReadOnlyList<LibrarySpectrum> spectra,
            SpecLibVersion version = DefaultWriteVersion,
            Action<int, int>? progress = null)
        {
            if (spectra == null) throw new ArgumentNullException(nameof(spectra));
            var entries = spectra.Select(DiaNNLibraryEntry.FromLibrarySpectrum).ToList();
            Write(filePath, entries, version, progress);
        }

        // -----------------------------------------------------------------------------------------
        // Core write pipeline
        // -----------------------------------------------------------------------------------------

        private static void WriteInternal(
            BinaryWriter writer,
            IReadOnlyList<DiaNNLibraryEntry> entries,
            SpecLibVersion version,
            Action<int, int>? progress)
        {
            int precursorCount = entries.Count;
            bool hasDecoys = entries.Any(e => e.IsDecoy);

            // ── Step 1: Build the string table ──────────────────────────────────────────────────
            // We need string offsets before we can fill in PrecursorRecord, so build it first.
            var (stringTable, sequenceOffsets, proteinOffsets, geneOffsets) =
                BuildStringTable(entries);

            // ── Step 2: Compute fragment block offsets ────────────────────────────────────────
            // Fragment blocks start immediately after the precursor array.
            // Precursor array starts after the 12-byte header (version + gen_decoys + count).
            const int headerBytes = 12;
            long precursorArrayBytes = (long)precursorCount * DiaNNBinaryStructs.SizeOf.PrecursorRecord;
            long fragmentBlocksStart = headerBytes + precursorArrayBytes;

            // Walk entries to compute each precursor's fragment data offset.
            var fragmentOffsets = new long[precursorCount];
            long runningFragmentOffset = fragmentBlocksStart;
            for (int i = 0; i < precursorCount; i++)
            {
                fragmentOffsets[i] = runningFragmentOffset;
                runningFragmentOffset += (long)entries[i].Fragments.Count * DiaNNBinaryStructs.SizeOf.FragmentRecord;
            }

            long stringTableStart = runningFragmentOffset;

            // ── Step 3: Write header ──────────────────────────────────────────────────────────
            writer.Write((int)version);
            writer.Write(hasDecoys ? 1 : 0);
            writer.Write(precursorCount);

            // ── Step 4: Write precursor records ───────────────────────────────────────────────
            // String table offsets are relative to stringTableStart so that they stay
            // meaningful regardless of where in the file the string table lands.
            // The reader (DiaNNSpecLibReader) applies the same convention.
            for (int i = 0; i < precursorCount; i++)
            {
                var entry = entries[i];
                var record = new DiaNNBinaryStructs.PrecursorRecord
                {
                    PrecursorMz       = (float)entry.PrecursorMz,
                    Charge            = (short)entry.PrecursorCharge,
                    RetentionTime     = (float)entry.RetentionTime,
                    IonMobility       = (float)entry.IonMobility,
                    SequenceOffset    = (int)sequenceOffsets[i],
                    SequenceLength    = (short)Math.Min(entry.ModifiedSequence?.Length ?? 0, short.MaxValue),
                    ProteinOffset     = (int)proteinOffsets[i],
                    GeneOffset        = (int)geneOffsets[i],
                    FragmentOffset    = fragmentOffsets[i],
                    FragmentCount     = (short)Math.Min(entry.Fragments.Count, short.MaxValue),
                    IsDecoy           = entry.IsDecoy ? (byte)1 : (byte)0,
                    IsProteotypic     = entry.IsProteotypic ? (byte)1 : (byte)0,
                    QValue            = (float)(entry.QValue ?? 0.0)
                };
                WritePrecursorRecord(writer, record);
                progress?.Invoke(i + 1, precursorCount);
            }

            // ── Step 5: Write fragment blocks ─────────────────────────────────────────────────
            // Fragments are written in the same order as the precursor array.
            // Within each precursor block, fragments are sorted by m/z (ascending) so the
            // reader can binary-search them during scoring.
            byte[] fragmentBuffer = ArrayPool<byte>.Shared.Rent(DiaNNBinaryStructs.SizeOf.FragmentRecord);
            try
            {
                for (int i = 0; i < precursorCount; i++)
                {
                    var sortedFragments = entries[i].Fragments
                        .OrderBy(f => f.Mz)
                        .ToList();

                    foreach (var frag in sortedFragments)
                    {
                        var record = new DiaNNBinaryStructs.FragmentRecord
                        {
                            Mz           = (float)frag.Mz,
                            Intensity    = (float)Math.Clamp(frag.Intensity, 0f, 1f),
                            IonType      = DiaNNBinaryStructs.IonTypeEncoding.ToByteFromChar(frag.IonType),
                            SeriesNumber = (byte)Math.Min(frag.SeriesNumber, byte.MaxValue),
                            Charge       = (byte)Math.Min(frag.Charge, byte.MaxValue),
                            LossType     = DiaNNBinaryStructs.LossTypeEncoding.ToByteFromString(frag.LossType)
                        };
                        WriteFragmentRecord(writer, record, fragmentBuffer);
                    }
                }
            }
            finally
            {
                ArrayPool<byte>.Shared.Return(fragmentBuffer);
            }

            // ── Step 6: Write string table ────────────────────────────────────────────────────
            writer.Write(stringTable);

            // ── Step 7: Flush ─────────────────────────────────────────────────────────────────
            writer.Flush();
        }

        // -----------------------------------------------------------------------------------------
        // String table construction
        // -----------------------------------------------------------------------------------------

        /// <summary>
        /// Builds a flat byte array of null-terminated UTF-8 strings and returns per-entry offsets.
        ///
        /// Deduplication: identical strings (same sequence, protein, gene) share a single entry in
        /// the table, saving significant space in large predicted libraries where the same protein
        /// group appears thousands of times.
        /// </summary>
        private static (byte[] table, int[] seqOffsets, int[] protOffsets, int[] geneOffsets)
            BuildStringTable(IReadOnlyList<DiaNNLibraryEntry> entries)
        {
            // Map from string → byte offset in table (for deduplication)
            var internTable = new Dictionary<string, int>(StringComparer.Ordinal);
            var tableBytes = new List<byte>(capacity: entries.Count * 20);

            int[] seqOffsets  = new int[entries.Count];
            int[] protOffsets = new int[entries.Count];
            int[] geneOffsets = new int[entries.Count];

            for (int i = 0; i < entries.Count; i++)
            {
                seqOffsets[i]  = InternString(entries[i].ModifiedSequence ?? string.Empty, internTable, tableBytes);
                protOffsets[i] = InternString(entries[i].ProteinId        ?? string.Empty, internTable, tableBytes);
                geneOffsets[i] = InternString(entries[i].GeneName          ?? string.Empty, internTable, tableBytes);
            }

            return (tableBytes.ToArray(), seqOffsets, protOffsets, geneOffsets);
        }

        /// <summary>
        /// Adds <paramref name="s"/> to the string table if not already present (deduplication),
        /// then returns its byte offset within the table.
        /// </summary>
        private static int InternString(string s, Dictionary<string, int> internTable, List<byte> tableBytes)
        {
            if (internTable.TryGetValue(s, out int existingOffset))
                return existingOffset;

            int offset = tableBytes.Count;
            internTable[s] = offset;

            // UTF-8 encode + null terminator
            byte[] encoded = Encoding.UTF8.GetBytes(s);
            tableBytes.AddRange(encoded);
            tableBytes.Add(0); // null terminator
            return offset;
        }

        // -----------------------------------------------------------------------------------------
        // Low-level struct writing (using MemoryMarshal for zero-copy struct → bytes)
        // -----------------------------------------------------------------------------------------

        private static void WritePrecursorRecord(BinaryWriter writer, DiaNNBinaryStructs.PrecursorRecord record)
        {
            int size = DiaNNBinaryStructs.SizeOf.PrecursorRecord;
            byte[] buffer = ArrayPool<byte>.Shared.Rent(size);
            try
            {
                MemoryMarshal.Write(buffer.AsSpan(0, size), ref record);
                writer.Write(buffer, 0, size);
            }
            finally
            {
                ArrayPool<byte>.Shared.Return(buffer);
            }
        }

        private static void WriteFragmentRecord(
            BinaryWriter writer,
            DiaNNBinaryStructs.FragmentRecord record,
            byte[] reusableBuffer)
        {
            int size = DiaNNBinaryStructs.SizeOf.FragmentRecord;
            MemoryMarshal.Write(reusableBuffer.AsSpan(0, size), ref record);
            writer.Write(reusableBuffer, 0, size);
        }

        // -----------------------------------------------------------------------------------------
        // IDisposable
        // -----------------------------------------------------------------------------------------

        private bool _disposed;

        public void Dispose()
        {
            if (_disposed) return;
            _disposed = true;
        }

        // -----------------------------------------------------------------------------------------
        // Diagnostic helpers
        // -----------------------------------------------------------------------------------------

        /// <summary>
        /// Estimates the output file size in bytes for a given list of entries.
        /// Useful for pre-flight checks before writing to disk.
        /// </summary>
        public static long EstimateFileSize(IReadOnlyList<DiaNNLibraryEntry> entries)
        {
            const int headerBytes = 12;
            long precursorBytes  = (long)entries.Count * DiaNNBinaryStructs.SizeOf.PrecursorRecord;
            long fragmentBytes   = entries.Sum(e => (long)e.Fragments.Count * DiaNNBinaryStructs.SizeOf.FragmentRecord);
            long stringBytes     = EstimateStringTableSize(entries);
            return headerBytes + precursorBytes + fragmentBytes + stringBytes;
        }

        private static long EstimateStringTableSize(IReadOnlyList<DiaNNLibraryEntry> entries)
        {
            // Estimate without full deduplication — worst case is sum of all string lengths + null bytes.
            long total = 0;
            foreach (var e in entries)
            {
                total += (e.ModifiedSequence?.Length ?? 0) + 1;
                total += (e.ProteinId?.Length        ?? 0) + 1;
                total += (e.GeneName?.Length          ?? 0) + 1;
            }
            return total;
        }

        /// <summary>
        /// Validates that a list of entries is well-formed before writing.
        /// Returns a list of validation errors (empty list = valid).
        /// </summary>
        public static List<string> ValidateEntries(IReadOnlyList<DiaNNLibraryEntry> entries)
        {
            var errors = new List<string>();
            var seen = new HashSet<string>(StringComparer.Ordinal);

            for (int i = 0; i < entries.Count; i++)
            {
                var e = entries[i];

                if (string.IsNullOrWhiteSpace(e.ModifiedSequence))
                    errors.Add($"Entry {i}: ModifiedSequence is null or empty.");

                if (e.PrecursorCharge <= 0)
                    errors.Add($"Entry {i} ({e.ModifiedSequence}): PrecursorCharge must be > 0 (got {e.PrecursorCharge}).");

                if (e.PrecursorMz <= 0)
                    errors.Add($"Entry {i} ({e.ModifiedSequence}): PrecursorMz must be > 0 (got {e.PrecursorMz:F4}).");

                if (e.Fragments == null || e.Fragments.Count == 0)
                    errors.Add($"Entry {i} ({e.ModifiedSequence}/{e.PrecursorCharge}): Has no fragment ions.");

                string key = $"{e.ModifiedSequence}/{e.PrecursorCharge}";
                if (!seen.Add(key))
                    errors.Add($"Entry {i}: Duplicate precursor key '{key}'.");

                if (e.Fragments != null)
                {
                    for (int j = 0; j < e.Fragments.Count; j++)
                    {
                        var f = e.Fragments[j];
                        if (f.Mz <= 0)
                            errors.Add($"Entry {i}, Fragment {j}: Mz must be > 0 (got {f.Mz:F4}).");
                        if (f.Intensity < 0 || f.Intensity > 1.0)
                            errors.Add($"Entry {i}, Fragment {j}: Intensity {f.Intensity:F4} is outside [0, 1]. Intensities should be normalized.");
                        if (f.SeriesNumber <= 0)
                            errors.Add($"Entry {i}, Fragment {j}: SeriesNumber must be > 0 (got {f.SeriesNumber}).");
                    }
                }
            }

            return errors;
        }
    }
}
