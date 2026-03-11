using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace Readers.SpectralLibrary.DiaNNSpectralLibrary
{
    /// <summary>
    /// Writes DIA-NN 2.3.2 .speclib binary spectral library files (format version −10).
    ///
    /// The binary specification is fully documented in Session7_SpecLibWrapup.md.
    ///
    /// Key design decisions:
    ///   - Always writes exactly 12 fragment slots per precursor (up to 11 real + 1 null).
    ///   - Fragment intensities are normalized so the maximum = 1.0 before writing.
    ///   - Fragment bytes 0–2 are written as safe placeholder values [1, 2, 1, 0] for real
    ///     records and [1, 1, 1, 0] for null slots (DIA-NN internal indices — not interpretable).
    ///   - The two shared-peptide markers are written as raw bytes, NOT float literals,
    ///     because their subnormal float values cannot be reliably round-tripped through
    ///     float arithmetic.
    ///   - A minimal secondary table stub is written (one group per protein, no joint groups).
    ///   - The file is written atomically: to a temp file first, then renamed.
    /// </summary>
    public static class DiaNNSpecLibWriter
    {
        // ─── Public API ──────────────────────────────────────────────────────────────

        /// <summary>
        /// Writes a list of <see cref="DiaNNLibraryEntry"/> objects to a DIA-NN .speclib file.
        /// </summary>
        /// <param name="entries">
        /// Library entries to write. Must not be null or empty.
        /// Entries are sorted by GlobalIndex before writing; if GlobalIndex is 0 for all,
        /// the current list order is used and indices are assigned sequentially.
        /// </param>
        /// <param name="outputPath">Destination file path (created or overwritten).</param>
        /// <param name="fastaPath">
        /// FASTA path string embedded in the header. Use the source FASTA path if known,
        /// or an empty string if not available.
        /// </param>
        /// <exception cref="ArgumentNullException">entries is null.</exception>
        /// <exception cref="ArgumentException">entries is empty.</exception>
        public static void Write(
            IReadOnlyList<DiaNNLibraryEntry> entries,
            string outputPath,
            string fastaPath = "")
        {
            if (entries == null) throw new ArgumentNullException(nameof(entries));
            if (entries.Count == 0)
                throw new ArgumentException("Cannot write an empty spectral library.", nameof(entries));

            // Sort by GlobalIndex; assign sequential indices if all are 0
            var sorted = entries
                .OrderBy(e => e.GlobalIndex)
                .ToList();

            bool needsIndexAssignment = sorted.All(e => e.GlobalIndex == 0);
            if (needsIndexAssignment)
            {
                for (int i = 0; i < sorted.Count; i++)
                    sorted[i].GlobalIndex = i;
            }

            // Collect unique proteins in the order they first appear
            var proteins = BuildProteinList(sorted);

            // Write to temp file, then rename atomically
            string tempPath = outputPath + ".tmp";
            try
            {
                using var stream = new FileStream(tempPath, FileMode.Create,
                    FileAccess.Write, FileShare.None, bufferSize: 1 << 20);
                using var writer = new BinaryWriter(stream, Encoding.UTF8, leaveOpen: false);

                WriteHeader(writer, fastaPath, proteins.Count);
                WriteProteinTable(writer, proteins);
                WriteSecondaryTable(writer, proteins, sorted);
                WritePrecursors(writer, sorted, proteins);
                WriteTrailingSection(writer, sorted.Count);
            }
            catch
            {
                if (File.Exists(tempPath)) File.Delete(tempPath);
                throw;
            }

            // Atomic rename
            if (File.Exists(outputPath)) File.Delete(outputPath);
            File.Move(tempPath, outputPath);
        }

        // ─── Header ───────────────────────────────────────────────────────────────────

        private static void WriteHeader(BinaryWriter w, string fastaPath, int nProteins)
        {
            // int32 = −10  (version)
            w.Write(DiaNNBinaryStructs.ExpectedVersion);

            // Four opaque flag int32 fields: 1, 1, 1, 0
            w.Write(1); w.Write(1); w.Write(1); w.Write(0);

            // FASTA path: length-prefixed UTF-8, NOT null-terminated
            byte[] fastaBytes = Encoding.UTF8.GetBytes(fastaPath);
            w.Write(fastaBytes.Length);
            w.Write(fastaBytes);

            // n_proteins
            w.Write(nProteins);
        }

        // ─── Protein table ────────────────────────────────────────────────────────────

        private static void WriteProteinTable(BinaryWriter w, List<ProteinSlot> proteins)
        {
            for (int i = 0; i < proteins.Count; i++)
            {
                var p = proteins[i];

                if (i == 0)
                {
                    // Protein 0: 2 opaque prefix int32 fields
                    w.Write(1); w.Write(0);
                }
                else
                {
                    // Proteins 1..N-1: 4 opaque prefix int32 fields
                    // Pattern observed: [protein_index, protein_index, 1, 0]
                    w.Write(i); w.Write(i); w.Write(1); w.Write(0);
                }

                WriteLengthPrefixedString(w, p.Accession);
                WriteLengthPrefixedString(w, p.Name);
                WriteLengthPrefixedString(w, p.Gene);
            }
        }

        // ─── Secondary table (minimal stub) ──────────────────────────────────────────

        private static void WriteSecondaryTable(
            BinaryWriter w, List<ProteinSlot> proteins, List<DiaNNLibraryEntry> entries)
        {
            // Minimal valid stub: one individual group per protein, no joint groups.
            //
            // Structure (from Session7_SpecLibWrapup.md §1.3):
            //   int32 = 1  (flag)
            //   int32 = 1  (flag)
            //   int32 = n_groups
            //   [n_groups group records]
            //   [pre-precursor preamble]

            int nGroups = proteins.Count;
            w.Write(1); w.Write(1);
            w.Write(nGroups);

            for (int gi = 0; gi < proteins.Count; gi++)
            {
                var p = proteins[gi];

                // 2 × int32 opaque fields (use zeros)
                w.Write(0); w.Write(0);

                // Accession / name / gene
                WriteLengthPrefixedString(w, p.Accession);
                WriteLengthPrefixedString(w, p.Name);
                WriteLengthPrefixedString(w, p.Gene);

                // Precursor index list: empty (int32 count = 0)
                w.Write(0);
            }

            // Pre-precursor preamble:
            //   float 0.0, float iRT_max, float 0.0, float −iRT_max, int32 n_precursors, int32 0
            float irtMax = entries.Count > 0
                ? entries.Max(e => Math.Abs((float)e.RetentionTime))
                : 0f;

            w.Write(0.0f);
            w.Write(irtMax);
            w.Write(0.0f);
            w.Write(-irtMax);
            w.Write(entries.Count);
            w.Write(0);
        }

        // ─── Precursor records ────────────────────────────────────────────────────────

        private static void WritePrecursors(
            BinaryWriter w, List<DiaNNLibraryEntry> entries, List<ProteinSlot> proteins)
        {
            var proteinIndexMap = proteins
                .Select((p, i) => (p.Accession, i))
                .ToDictionary(x => x.Accession, x => x.i);

            for (int seq = 0; seq < entries.Count; seq++)
            {
                var e = entries[seq];

                // Normalize fragment intensities
                var realFragments = e.Fragments
                    .Where(f => f.Mz > 0)
                    .ToList();

                float maxInty = realFragments.Count > 0
                    ? (float)realFragments.Max(f => f.Intensity)
                    : 1.0f;
                if (maxInty <= 0) maxInty = 1.0f;

                // Sort by intensity descending; take at most 11
                var topFragments = realFragments
                    .OrderByDescending(f => f.Intensity)
                    .Take(DiaNNBinaryStructs.FragmentSlotsPerPrecursor - 1)
                    .ToList();

                // Resolve protein group index
                int protGroupIndex = e.ProteinGroupIndex;
                if (protGroupIndex == 0 && !string.IsNullOrEmpty(e.ProteinId))
                {
                    // Re-derive from our protein table if available
                    if (proteinIndexMap.TryGetValue(e.ProteinId, out int pIdx))
                        protGroupIndex = pIdx;
                }

                int charge           = e.PrecursorCharge;
                int strippedLen      = e.StrippedSequence?.Length ?? 0;
                float precursorMz    = (float)e.PrecursorMz;
                float irt            = (float)e.RetentionTime;
                float ionMobility    = (float)e.IonMobility;

                // Build name string: ModifiedSequence + charge digit
                string nameStr = BuildNameString(e.ModifiedSequence, charge);
                byte[] nameBytes = Encoding.ASCII.GetBytes(nameStr);

                // ── Pre-name block (48 bytes) ──────────────────────────────────────
                // Written in field order matching nlp offsets
                w.Write(seq);            // global_index    (nlp-48)
                w.Write(charge);         // charge          (nlp-44)
                w.Write(strippedLen);    // stripped_len    (nlp-40)
                w.Write(precursorMz);    // precursor_mz    (nlp-36)
                w.Write(irt);            // iRT             (nlp-32)
                w.Write(ionMobility);    // ion_mobility    (nlp-28)
                w.Write(1.0f);           // placeholder     (nlp-24)
                w.Write(1.0f);           // placeholder     (nlp-20)
                w.Write(1.0f);           // placeholder     (nlp-16)
                w.Write(1.0f);           // placeholder     (nlp-12)
                WriteMarker(w, e.TerminusType); // marker (nlp-8)
                w.Write(protGroupIndex); // protein_group_index (nlp-4)

                // ── Name field ────────────────────────────────────────────────────
                w.Write(nameBytes.Length);  // name_length (nlp+0)
                w.Write(nameBytes);          // name        (nlp+4)

                // ── Post-name block (28 bytes) ─────────────────────────────────────
                float topMz = topFragments.Count > 0 ? (float)topFragments[0].Mz : 0f;
                float topInty = topFragments.Count > 0 ? (float)topFragments[0].Intensity / maxInty : 0f;

                w.Write(0.0f);                                  // +0
                w.Write(1.0f);                                  // +4
                w.Write(0.0f);                                  // +8
                w.Write(0.0f);                                  // +12
                w.Write(DiaNNBinaryStructs.FragmentSlotsPerPrecursor); // n_fragments +16
                w.Write(topMz);                                 // top_fragment_mz +20
                w.Write(topInty);                               // top_fragment_intensity +24

                // ── Fragment records (12 × 12 bytes) ──────────────────────────────
                // Write up to 11 real records
                foreach (var frag in topFragments)
                {
                    w.Write(DiaNNBinaryStructs.RealFragmentOpaqueBytes);
                    w.Write((float)frag.Mz);
                    w.Write((float)frag.Intensity / maxInty);
                }

                // Pad with null slots to reach exactly 12 total
                int nullCount = DiaNNBinaryStructs.FragmentSlotsPerPrecursor - topFragments.Count;
                for (int n = 0; n < nullCount; n++)
                {
                    w.Write(DiaNNBinaryStructs.NullFragmentOpaqueBytes);
                    w.Write(0.0f); // mz = 0 → null placeholder
                    w.Write(0.0f);
                }
            }
        }

        // ─── Trailing section ─────────────────────────────────────────────────────────

        private static void WriteTrailingSection(BinaryWriter w, int nPrecursors)
        {
            // n_precursors + 1 int32 values.
            // A simple monotonic sequence suffices; the exact semantic is not confirmed.
            for (int i = 0; i < nPrecursors; i++)
                w.Write(i);
            // Sentinel
            w.Write(0);
        }

        // ─── Marker writer ────────────────────────────────────────────────────────────

        private static void WriteMarker(BinaryWriter w, DiaNNBinaryStructs.TerminusType t)
        {
            switch (t)
            {
                case DiaNNBinaryStructs.TerminusType.Internal:
                    w.Write(2.9375f); break;
                case DiaNNBinaryStructs.TerminusType.NTerminal:
                    w.Write(2.8906f); break;
                case DiaNNBinaryStructs.TerminusType.CTerminal:
                    w.Write(2.8438f); break;
                case DiaNNBinaryStructs.TerminusType.NAndCTerminal:
                    w.Write(2.7969f); break;
                case DiaNNBinaryStructs.TerminusType.SharedConflict:
                    // Subnormal float — MUST be written as raw bytes, not a float literal.
                    w.Write(new byte[] { 0x07, 0x00, 0x3D, 0x00 }); break;
                case DiaNNBinaryStructs.TerminusType.SharedSameType:
                    // Subnormal float — MUST be written as raw bytes.
                    w.Write(new byte[] { 0x07, 0x00, 0x39, 0x00 }); break;
                default:
                    w.Write(2.9375f); break;
            }
        }

        // ─── Helpers ──────────────────────────────────────────────────────────────────

        /// <summary>
        /// Constructs the DIA-NN name string: ModifiedSequence (with (UniMod:N) tags) + charge digit.
        /// The modified sequence is already in parenthesis notation (as stored in DiaNNLibraryEntry).
        /// </summary>
        private static string BuildNameString(string modifiedSequence, int charge)
        {
            // modifiedSequence is stored without flanking underscores and without charge digit
            return modifiedSequence + charge.ToString();
        }

        private static void WriteLengthPrefixedString(BinaryWriter w, string s)
        {
            byte[] bytes = string.IsNullOrEmpty(s)
                ? Array.Empty<byte>()
                : Encoding.UTF8.GetBytes(s);
            w.Write(bytes.Length);
            if (bytes.Length > 0) w.Write(bytes);
        }

        /// <summary>
        /// Collects the ordered list of unique proteins from the library entries,
        /// preserving first-appearance order.
        /// </summary>
        private static List<ProteinSlot> BuildProteinList(List<DiaNNLibraryEntry> entries)
        {
            var seen = new HashSet<string>(StringComparer.Ordinal);
            var result = new List<ProteinSlot>();

            foreach (var e in entries)
            {
                string acc = e.ProteinId ?? string.Empty;
                if (seen.Add(acc))
                {
                    result.Add(new ProteinSlot(acc, e.ProteinName ?? string.Empty, e.GeneName ?? string.Empty));
                }
            }

            // Ensure at least one protein entry even for empty metadata
            if (result.Count == 0)
                result.Add(new ProteinSlot("unknown", "unknown", "unknown"));

            return result;
        }

        private readonly struct ProteinSlot
        {
            public readonly string Accession;
            public readonly string Name;
            public readonly string Gene;

            public ProteinSlot(string accession, string name, string gene)
            {
                Accession = accession; Name = name; Gene = gene;
            }
        }
    }
}
