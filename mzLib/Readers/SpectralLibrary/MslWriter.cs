using MassSpectrometry;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using System.Buffers;
using System.Runtime.InteropServices;
using System.Text;

namespace Readers.SpectralLibrary;

/// <summary>
/// Static writer for the .msl (mzLib Spectral Library) binary format.
///
/// Writing is a strict two-pass algorithm:
///   Pass 1 — <see cref="MslWriteLayout"/> computes every byte offset, string index,
///             protein deduplication table, and elution group assignment arithmetically
///             before a single byte is written.
///   Pass 2 — the file is written sequentially from byte 0 to EOF with no Seek() calls;
///             <c>MslPrecursorRecord.FragmentBlockOffset</c> is filled from the pre-computed
///             layout rather than from a back-patch.
///
/// The output is written atomically: bytes go to a temp file first, then the temp file
/// is renamed to the final path with <c>File.Move(overwrite: true)</c> so that a failed
/// write never leaves a partially written file at the destination path.
/// </summary>
public static class MslWriter
{
	// ── Static constructor ───────────────────────────────────────────────────

	/// <summary>
	/// Runs <see cref="MslStructs.SizeCheck"/> once when the class is first used so that
	/// a Pack-setting mistake is caught immediately rather than silently corrupting files.
	/// </summary>
	static MslWriter()
	{
		MslStructs.SizeCheck();
	}

	// ────────────────────────────────────────────────────────────────────────
	// Magic endianness constant
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// The file magic stored as a uint for use inside <c>MslFileHeader.Magic</c> and
	/// <c>MslFooter.TrailingMagic</c> struct fields, which are serialized to disk by
	/// <see cref="MemoryMarshal.Write{T}"/> in native (little-endian) byte order on
	/// x86/x64 platforms.
	///
	/// The canonical on-disk magic is the four ASCII bytes "MZLB": 0x4D 0x5A 0x4C 0x42.
	/// When MemoryMarshal serializes a uint little-endian, it writes the least-significant
	/// byte first. To produce 0x4D 0x5A 0x4C 0x42 on disk the uint must therefore hold the
	/// byte-swapped value 0x424C5A4Du (decimal 1112756813), so that the four bytes written
	/// are reversed back to the canonical order.
	///
	/// Do NOT use <see cref="MslFormat.MagicAsUInt32"/> (= 0x4D5A4C42u) directly inside
	/// a struct field written via MemoryMarshal; that value is the big-endian form and
	/// would produce the wrong on-disk bytes 0x42 0x4C 0x5A 0x4D.
	/// </summary>
	private const uint MagicForLEStruct = 0x424C_5A4Du;

	// ────────────────────────────────────────────────────────────────────────
	// Public API
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Writes a list of <see cref="MslLibraryEntry"/> objects to a .msl binary file.
	///
	/// The output path is written atomically: bytes go to a temp file first, then the
	/// temp file is renamed to the final path so that a failed write never leaves a
	/// partially written file at the destination.
	///
	/// Fragment intensities are normalized per precursor before writing (max = 1.0).
	/// Fragment records within each precursor are sorted by m/z ascending to support
	/// binary search during scoring.
	/// </summary>
	/// <param name="outputPath">
	///   Destination file path. The file is created or overwritten. The directory must
	///   already exist; this method does not create parent directories.
	/// </param>
	/// <param name="entries">
	///   Ordered list of library entries to write. The on-disk precursor order mirrors
	///   the list order. Must not be null; may be empty (produces a valid zero-precursor file).
	/// </param>
	/// <exception cref="ArgumentNullException">
	///   Thrown when <paramref name="outputPath"/> or <paramref name="entries"/> is null.
	/// </exception>
	/// <exception cref="NotSupportedException">
	///   Thrown when any fragment carries <see cref="MslFormat.NeutralLossCode.Custom"/>,
	///   which requires the extended annotation table not yet implemented.
	/// </exception>
	public static void Write(string outputPath, IReadOnlyList<MslLibraryEntry> entries)
	{
		if (outputPath is null) throw new ArgumentNullException(nameof(outputPath));
		if (entries is null) throw new ArgumentNullException(nameof(entries));

		// Derive a temp path in the same directory so the final rename is an atomic
		// same-volume move rather than a cross-device copy.
		string tempPath = outputPath + ".tmp~";

		try
		{
			// ── Pass 1: compute layout (all offsets, indices, deduplication) ─────
			var layout = new MslWriteLayout(entries);

			// ── Pass 2a: write all sections except the footer ────────────────────
			// The writer is fully closed before the CRC read-back so that the OS
			// releases its exclusive write lock on the temp file. Attempting to
			// open the same file for reading while a BinaryWriter still holds it
			// open with FileShare.None raises an IOException on Windows.
			using (var fs = new FileStream(tempPath, FileMode.Create, FileAccess.Write, FileShare.None,
											   bufferSize: 65536, FileOptions.SequentialScan))
			using (var writer = new BinaryWriter(fs, Encoding.UTF8, leaveOpen: false))
			{
				WriteHeader(writer, layout);
				WriteProteinTable(writer, layout);
				WriteStringTable(writer, layout);
				WritePrecursorArray(writer, layout, entries);
				WriteFragmentBlocks(writer, layout, entries);
				WriteOffsetTable(writer, layout);
				// writer and fs are disposed here; OS write lock is released
			}

			// ── Pass 2b: compute CRC over the closed file, then append footer ─────
			// A fresh read stream is safe because the write stream above has been
			// fully disposed. The footer is appended in a new write stream.
			uint crc = ComputeCrc32(tempPath, layout.DataEndOffset);

			using (var fs = new FileStream(tempPath, FileMode.Append, FileAccess.Write, FileShare.None,
											   bufferSize: 128))
			using (var writer = new BinaryWriter(fs, Encoding.UTF8, leaveOpen: false))
			{
				WriteFooter(writer, layout, crc);
			}

			// Atomic rename: replaces any existing file at outputPath
			File.Move(tempPath, outputPath, overwrite: true);
		}
		catch
		{
			// Clean up the temp file on any failure so stale partials don't accumulate
			if (File.Exists(tempPath))
				File.Delete(tempPath);
			throw;
		}
	}

	/// <summary>
	/// Convenience overload that converts <see cref="LibrarySpectrum"/> objects to
	/// <see cref="MslLibraryEntry"/> via <see cref="MslLibraryEntry.FromLibrarySpectrum"/>
	/// before writing. All extended fields (IonMobility, QValue, protein metadata, etc.)
	/// are set to their documented defaults by the conversion method.
	/// </summary>
	/// <param name="outputPath">
	///   Destination file path. Follows the same atomicity guarantee as
	///   <see cref="Write(string, IReadOnlyList{MslLibraryEntry})"/>.
	/// </param>
	/// <param name="spectra">
	///   Source library spectra. Must not be null. Each spectrum must have a non-null,
	///   non-empty <c>MatchedFragmentIons</c> list or the resulting entry will fail validation.
	/// </param>
	/// <exception cref="ArgumentNullException">
	///   Thrown when <paramref name="outputPath"/> or <paramref name="spectra"/> is null.
	/// </exception>
	public static void WriteFromLibrarySpectra(string outputPath,
		IReadOnlyList<LibrarySpectrum> spectra)
	{
		if (outputPath is null) throw new ArgumentNullException(nameof(outputPath));
		if (spectra is null) throw new ArgumentNullException(nameof(spectra));

		// Convert each LibrarySpectrum to an MslLibraryEntry; preserve list order
		var entries = new List<MslLibraryEntry>(spectra.Count);
		foreach (LibrarySpectrum spectrum in spectra)
			entries.Add(MslLibraryEntry.FromLibrarySpectrum(spectrum));

		Write(outputPath, entries);
	}

	/// <summary>
	/// Returns the estimated file size in bytes without writing any data.
	/// Useful for pre-flight disk-space checks before committing to a write.
	///
	/// The estimate uses the formula:
	///   HeaderSize + NProteins × ProteinRecordSize
	///   + StringTable(approx) + NPrecursors × PrecursorRecordSize
	///   + NPrecursors × AvgFragments × FragmentRecordSize
	///   + NPrecursors × 8 (offset table) + FooterSize
	///
	/// Actual size may vary by up to ±5% depending on UTF-8 encoding lengths and
	/// the number of unique strings vs total string references.
	/// </summary>
	/// <param name="nPrecursors">
	///   Number of precursor entries that will be written.
	/// </param>
	/// <param name="avgFragmentsPerPrecursor">
	///   Average number of fragment records per precursor.
	/// </param>
	/// <param name="avgStringBytes">
	///   Average UTF-8 byte length of each unique string (sequence, protein accession, etc.).
	///   A value of 20 is a reasonable default for tryptic peptide libraries.
	/// </param>
	/// <returns>Estimated total file size in bytes.</returns>
	public static long EstimateFileSize(int nPrecursors, int avgFragmentsPerPrecursor,
		int avgStringBytes)
	{
		// Estimate unique strings in the string table.
		//
		// Each precursor contributes up to 3 string fields that are likely to be unique
		// per-precursor (modified sequence, stripped sequence, protein accession). The
		// protein name and gene name are typically shared across many precursors within a
		// protein, so they are not counted per-precursor. A deduplication factor of 0.6
		// is applied to account for charge-state variants sharing the same sequences.
		// The +1 accounts for the mandatory empty-string entry at index 0.
		long estimatedUniqueStrings = (long)(nPrecursors * 3 * 0.6) + 1;

		// String table on disk: 4-byte NStrings header + 4-byte TotalBytes header
		// + per-entry: 4-byte length prefix + UTF-8 body
		long stringTableBytes = 8 + estimatedUniqueStrings * (4 + avgStringBytes);

		// Fragment section: one MslFragmentRecord (20 bytes) per fragment per precursor
		long fragmentBytes = (long)nPrecursors * avgFragmentsPerPrecursor * MslFormat.FragmentRecordSize;

		return MslFormat.HeaderSize
			 + (long)nPrecursors * MslFormat.ProteinRecordSize   // ~1 protein record per precursor
			 + stringTableBytes
			 + (long)nPrecursors * MslFormat.PrecursorRecordSize
			 + fragmentBytes
			 + (long)nPrecursors * sizeof(long)                   // offset table: 8 bytes per precursor
			 + MslFormat.FooterSize;
	}

	/// <summary>
	/// Validates a list of <see cref="MslLibraryEntry"/> objects before writing.
	/// Returns a list of human-readable error strings. An empty return list means all
	/// entries passed every check and the list is safe to pass to <see cref="Write"/>.
	///
	/// Checks performed:
	/// <list type="bullet">
	///   <item>Modified sequence is not null or empty.</item>
	///   <item>Charge is greater than zero.</item>
	///   <item>PrecursorMz is greater than zero.</item>
	///   <item>At least one fragment is present (FragmentCount > 0).</item>
	///   <item>Each fragment has Mz > 0.</item>
	///   <item>Each fragment has Charge > 0.</item>
	///   <item>For internal ions: SecondaryFragmentNumber > FragmentNumber.</item>
	/// </list>
	/// </summary>
	/// <param name="entries">The entries to validate. Must not be null.</param>
	/// <returns>A list of error description strings; empty when all entries are valid.</returns>
	public static List<string> ValidateEntries(IReadOnlyList<MslLibraryEntry> entries)
	{
		if (entries is null) throw new ArgumentNullException(nameof(entries));

		var errors = new List<string>();

		for (int i = 0; i < entries.Count; i++)
		{
			MslLibraryEntry entry = entries[i];

			// Modified sequence must be a non-empty string
			if (string.IsNullOrEmpty(entry.ModifiedSequence))
				errors.Add($"Entry [{i}]: ModifiedSequence is null or empty.");

			// Precursor charge must be positive
			if (entry.Charge <= 0)
				errors.Add($"Entry [{i}] '{entry.ModifiedSequence}': Charge must be > 0 (got {entry.Charge}).");

			// Precursor m/z must be a positive, finite value
			if (entry.PrecursorMz <= 0.0 || !double.IsFinite(entry.PrecursorMz))
				errors.Add($"Entry [{i}] '{entry.ModifiedSequence}': PrecursorMz must be > 0 (got {entry.PrecursorMz}).");

			// At least one fragment ion is required for a useful library entry
			if (entry.Fragments == null || entry.Fragments.Count == 0)
			{
				errors.Add($"Entry [{i}] '{entry.ModifiedSequence}': FragmentCount must be > 0.");
				continue; // Skip per-fragment checks; Fragments list may be null
			}

			// Per-fragment checks
			for (int j = 0; j < entry.Fragments.Count; j++)
			{
				MslFragmentIon frag = entry.Fragments[j];

				// Fragment m/z must be a positive, finite value
				if (frag.Mz <= 0f || !float.IsFinite(frag.Mz))
					errors.Add($"Entry [{i}] '{entry.ModifiedSequence}', fragment [{j}]: " +
							   $"Mz must be > 0 (got {frag.Mz}).");

				// Fragment charge must be positive
				if (frag.Charge <= 0)
					errors.Add($"Entry [{i}] '{entry.ModifiedSequence}', fragment [{j}]: " +
							   $"Charge must be > 0 (got {frag.Charge}).");

				// Internal ion residue range must be valid (end > start)
				if (frag.IsInternalFragment && frag.SecondaryFragmentNumber <= frag.FragmentNumber)
					errors.Add($"Entry [{i}] '{entry.ModifiedSequence}', fragment [{j}]: " +
							   $"Internal ion SecondaryFragmentNumber ({frag.SecondaryFragmentNumber}) " +
							   $"must be > FragmentNumber ({frag.FragmentNumber}).");
			}
		}

		return errors;
	}

	// ────────────────────────────────────────────────────────────────────────
	// Pass-2 write helpers  (called in strict file-section order)
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Writes the 64-byte <see cref="MslFileHeader"/> to the current stream position.
	/// All offset and count fields are taken from the pre-computed <paramref name="layout"/>.
	/// </summary>
	/// <param name="writer">Open binary writer positioned at byte 0.</param>
	/// <param name="layout">Pre-computed write layout containing all offset values.</param>
	private static void WriteHeader(BinaryWriter writer, MslWriteLayout layout)
	{
		// Populate every field from the pre-computed layout
		var header = new MslFileHeader
		{
			// MslFileHeader.Magic is a uint field serialized by MemoryMarshal.Write, which
			// writes the value in native (little-endian) byte order on x86/x64. To produce
			// the canonical on-disk bytes 0x4D 0x5A 0x4C 0x42 ("MZLB") we store the
			// byte-swapped value 0x424C5A4Du so that the LE serialization reverses it back.
			// MslFormat.MagicAsUInt32 (0x4D5A4C42u) is the big-endian representation and
			// must NOT be used directly here.
			Magic = MagicForLEStruct,
			FormatVersion = MslFormat.CurrentVersion,
			FileFlags = layout.FileFlags,
			NPrecursors = layout.NPrecursors,
			NProteins = layout.NProteins,
			NElutionGroups = layout.NElutionGroups,
			NStrings = layout.NStrings,
			Reserved = 0,
			ProteinTableOffset = layout.ProteinTableOffset,
			StringTableOffset = layout.StringTableOffset,
			PrecursorSectionOffset = layout.PrecursorSectionOffset,
			FragmentSectionOffset = layout.FragmentSectionOffset
		};

		WriteStruct(writer, header);
	}

	/// <summary>
	/// Writes one <see cref="MslProteinRecord"/> per unique protein in the layout's protein
	/// slot list. Fields are taken from the <see cref="MslWriteLayout.ProteinSlots"/> list
	/// and the pre-computed string-table indices.
	/// </summary>
	/// <param name="writer">Open binary writer positioned at the start of the protein table.</param>
	/// <param name="layout">Pre-computed write layout.</param>
	private static void WriteProteinTable(BinaryWriter writer, MslWriteLayout layout)
	{
		foreach (ProteinSlot slot in layout.ProteinSlots)
		{
			var record = new MslProteinRecord
			{
				AccessionStringIdx = slot.AccessionStringIdx,
				NameStringIdx = slot.NameStringIdx,
				GeneStringIdx = slot.GeneStringIdx,
				ProteinGroupId = 0,                   // not yet assigned; default to 0
				NPrecursors = slot.PrecursorCount,
				ProteinFlags = 0                    // is_reviewed not determined here
			};

			WriteStruct(writer, record);
		}
	}

	/// <summary>
	/// Writes the string table section.
	///
	/// On-disk layout:
	/// <code>
	///   int32   NStrings       — total number of string entries
	///   int32   TotalBytes     — total byte count of all encoded string bodies
	///   For each string:
	///     int32  length        — byte count of the UTF-8 body (no null terminator)
	///     bytes  body          — UTF-8 encoded characters
	/// </code>
	/// The empty string "" is always at index 0.
	/// </summary>
	/// <param name="writer">Open binary writer positioned at the start of the string table.</param>
	/// <param name="layout">Pre-computed write layout containing the ordered string list.</param>
	private static void WriteStringTable(BinaryWriter writer, MslWriteLayout layout)
	{
		// Header: entry count and total body byte size
		writer.Write(layout.StringList.Count);
		writer.Write(layout.StringTableBodyBytes);

		// Each string: 4-byte length prefix + UTF-8 body (no null terminator)
		foreach (string s in layout.StringList)
		{
			byte[] encoded = Encoding.UTF8.GetBytes(s);
			writer.Write(encoded.Length);
			writer.Write(encoded);
		}
	}

	/// <summary>
	/// Writes one <see cref="MslPrecursorRecord"/> per entry in <paramref name="entries"/>,
	/// preserving the caller's list order. All integer indices and offsets come from the
	/// pre-computed <paramref name="layout"/>.
	/// </summary>
	/// <param name="writer">Open binary writer positioned at the start of the precursor array.</param>
	/// <param name="layout">Pre-computed write layout.</param>
	/// <param name="entries">Source library entries in write order.</param>
	private static void WritePrecursorArray(BinaryWriter writer, MslWriteLayout layout,
		IReadOnlyList<MslLibraryEntry> entries)
	{
		for (int i = 0; i < entries.Count; i++)
		{
			MslLibraryEntry entry = entries[i];
			MslPrecursorLayout pl = layout.PrecursorLayouts[i];

			var record = new MslPrecursorRecord
			{
				PrecursorMz = (float)entry.PrecursorMz,
				Irt = (float)entry.Irt,
				IonMobility = (float)entry.IonMobility,
				Charge = (short)entry.Charge,
				FragmentCount = (short)entry.Fragments.Count,
				ElutionGroupId = entry.ElutionGroupId,
				ProteinIdx = pl.ProteinIdx,
				ModifiedSeqStringIdx = pl.ModifiedSeqStringIdx,
				StrippedSeqStringIdx = pl.StrippedSeqStringIdx,
				FragmentBlockOffset = pl.FragmentBlockOffset,
				QValue = entry.QValue,
				StrippedSeqLength = string.IsNullOrEmpty(entry.StrippedSequence)
										   ? 0
										   : entry.StrippedSequence.Length,
				MoleculeType = (short)entry.MoleculeType,
				DissociationType = (short)entry.DissociationType,
				Nce = (short)(entry.Nce * 10),   // stored as NCE × 10
				PrecursorFlags = MslFormat.EncodePrecursorFlags(
										   isDecoy: entry.IsDecoy,
										   isProteotypic: entry.IsProteotypic,
										   rtCalibrated: false),
				SourceType = (byte)entry.Source
			};

			WriteStruct(writer, record);
		}
	}

	/// <summary>
	/// Writes the fragment blocks section: for each precursor, writes its fragment ions as
	/// contiguous <see cref="MslFragmentRecord"/> structs in m/z ascending order (sorted
	/// during pass 1). Intensities must already be normalized before this is called (done
	/// in <see cref="MslWriteLayout"/> constructor).
	/// </summary>
	/// <param name="writer">Open binary writer positioned at the start of the fragment section.</param>
	/// <param name="layout">Pre-computed write layout (used for NeutralLossCode per fragment).</param>
	/// <param name="entries">Source library entries; fragment lists are already sorted and normalized.</param>
	private static void WriteFragmentBlocks(BinaryWriter writer, MslWriteLayout layout,
		IReadOnlyList<MslLibraryEntry> entries)
	{
		foreach (MslLibraryEntry entry in entries)
		{
			foreach (MslFragmentIon frag in entry.Fragments)
			{
				// Map the neutral-loss double to the 3-bit NeutralLossCode enum value.
				// Custom losses are not yet supported and were already checked during layout.
				MslFormat.NeutralLossCode lossCode = ClassifyNeutralLoss(frag.NeutralLoss);

				// Build the packed flags byte from the four independent flag properties
				byte flags = MslFormat.EncodeFragmentFlags(
					isInternal: frag.IsInternalFragment,
					isDiagnostic: frag.IsDiagnosticIon,
					lossCode: lossCode,
					excludeFromQuant: frag.ExcludeFromQuant);

				var record = new MslFragmentRecord
				{
					Mz = frag.Mz,
					Intensity = frag.Intensity,
					ProductType = (short)(int)frag.ProductType,
					SecondaryProductType = frag.SecondaryProductType.HasValue
											   ? (short)(int)frag.SecondaryProductType.Value
											   : (short)-1,
					FragmentNumber = (short)frag.FragmentNumber,
					SecondaryFragmentNumber = (short)frag.SecondaryFragmentNumber,
					ResiduePosition = (short)frag.ResiduePosition,
					Charge = (byte)frag.Charge,
					Flags = flags
				};

				WriteStruct(writer, record);
			}
		}
	}

	/// <summary>
	/// Writes the offset table: one int64 per precursor giving the absolute file byte offset
	/// of that precursor's <see cref="MslFragmentRecord"/> block. The offsets come directly
	/// from the pre-computed <see cref="MslWriteLayout.PrecursorLayouts"/> list.
	/// </summary>
	/// <param name="writer">Open binary writer positioned at the start of the offset table.</param>
	/// <param name="layout">Pre-computed write layout.</param>
	private static void WriteOffsetTable(BinaryWriter writer, MslWriteLayout layout)
	{
		foreach (MslPrecursorLayout pl in layout.PrecursorLayouts)
			writer.Write(pl.FragmentBlockOffset);
	}

	/// <summary>
	/// Writes the 20-byte <see cref="MslFooter"/> as the final bytes of the file.
	/// </summary>
	/// <param name="writer">Open binary writer positioned at the start of the footer.</param>
	/// <param name="layout">Pre-computed write layout (supplies OffsetTableOffset and NPrecursors).</param>
	/// <param name="crc32">CRC-32 computed over bytes 0..(OffsetTableOffset−1) by the caller.</param>
	private static void WriteFooter(BinaryWriter writer, MslWriteLayout layout, uint crc32)
	{
		var footer = new MslFooter
		{
			OffsetTableOffset = layout.OffsetTableOffset,
			NPrecursors = layout.NPrecursors,
			DataCrc32 = crc32,
			TrailingMagic = MagicForLEStruct
		};

		WriteStruct(writer, footer);
	}

	// ────────────────────────────────────────────────────────────────────────
	// Struct serialization helper
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Serializes a blittable struct of type <typeparamref name="T"/> to the
	/// <paramref name="writer"/> stream using <see cref="MemoryMarshal.Write{T}"/> so that
	/// the entire struct is written in one span operation rather than one field at a time.
	/// A byte buffer of the correct size is rented from <see cref="ArrayPool{T}"/> to avoid
	/// per-call heap allocations on the hot path.
	/// </summary>
	/// <typeparam name="T">An unmanaged (blittable) struct type with a known Marshal size.</typeparam>
	/// <param name="writer">The <see cref="BinaryWriter"/> to write into.</param>
	/// <param name="value">The struct value to serialize.</param>
	private static void WriteStruct<T>(BinaryWriter writer, T value) where T : unmanaged
	{
		// Determine the exact byte count the OS will see for this struct
		int size = Marshal.SizeOf<T>();

		// Rent a buffer to avoid a heap allocation per struct write
		byte[] buffer = ArrayPool<byte>.Shared.Rent(size);
		try
		{
			// Write the struct's raw bytes into the buffer via MemoryMarshal
			MemoryMarshal.Write(buffer.AsSpan(0, size), ref value);
			writer.Write(buffer, 0, size);
		}
		finally
		{
			// Always return the rented buffer, even on exception
			ArrayPool<byte>.Shared.Return(buffer);
		}
	}

	// ────────────────────────────────────────────────────────────────────────
	// CRC-32 computation
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Computes the CRC-32/ISO-HDLC checksum over bytes 0..(dataEndOffset−1) of the file
	/// at <paramref name="filePath"/>. The file must be flushed to disk before calling this.
	///
	/// Uses a self-contained lookup-table implementation (standard 0xEDB88320 reflected
	/// polynomial) to avoid any dependency on <c>System.IO.Hashing</c>, which is only
	/// available in .NET 6+ and may not be present in all target frameworks.
	/// Reads the file in 64 KiB chunks to bound peak memory usage.
	/// </summary>
	/// <param name="filePath">Path of the file to hash.</param>
	/// <param name="dataEndOffset">
	///   Number of bytes to include in the hash (exclusive upper bound).
	///   Bytes at or beyond this offset (the offset table and footer) are excluded.
	/// </param>
	/// <returns>The CRC-32 checksum as a uint32.</returns>
	private static uint ComputeCrc32(string filePath, long dataEndOffset)
	{
		// 64 KiB read buffer — large enough to amortise I/O overhead without wasting memory
		const int ChunkSize = 65536;

		uint crc = 0xFFFF_FFFFu; // standard CRC-32 initial value
		byte[] buffer = ArrayPool<byte>.Shared.Rent(ChunkSize);

		try
		{
			using var fs = new FileStream(filePath, FileMode.Open, FileAccess.Read,
										  FileShare.None, bufferSize: ChunkSize);

			long remaining = dataEndOffset;

			while (remaining > 0)
			{
				// Read up to ChunkSize bytes but never past dataEndOffset
				int toRead = (int)Math.Min(remaining, ChunkSize);
				int read = fs.Read(buffer, 0, toRead);

				if (read == 0)
					break; // Unexpected EOF — CRC will be over whatever was read

				// Feed each byte through the CRC-32 lookup table
				for (int i = 0; i < read; i++)
					crc = (crc >> 8) ^ Crc32Table[(crc ^ buffer[i]) & 0xFF];

				remaining -= read;
			}
		}
		finally
		{
			ArrayPool<byte>.Shared.Return(buffer);
		}

		// Final XOR with 0xFFFFFFFF is part of the CRC-32/ISO-HDLC specification
		return crc ^ 0xFFFF_FFFFu;
	}

	/// <summary>
	/// Precomputed CRC-32/ISO-HDLC lookup table generated from the standard reflected
	/// polynomial 0xEDB88320. Each entry gives the CRC contribution of one input byte value.
	/// Computed once at class load time; stored as a static readonly array to avoid
	/// regeneration overhead on every Write() call.
	/// </summary>
	private static readonly uint[] Crc32Table = BuildCrc32Table();

	/// <summary>
	/// Builds the 256-entry CRC-32 lookup table using the standard reflected polynomial
	/// 0xEDB88320 (the bit-reversed form of the IEEE 802.3 / PKZIP polynomial 0x04C11DB7).
	/// </summary>
	/// <returns>A 256-element uint array where index <c>i</c> is the CRC remainder for byte value <c>i</c>.</returns>
	private static uint[] BuildCrc32Table()
	{
		// Standard reflected CRC-32 polynomial (IEEE 802.3 / PKZIP / zlib convention)
		const uint Polynomial = 0xEDB8_8320u;

		var table = new uint[256];

		for (uint i = 0; i < 256; i++)
		{
			// Process the 8 bits of each possible byte value
			uint entry = i;
			for (int bit = 0; bit < 8; bit++)
				entry = (entry & 1u) != 0
					? (entry >> 1) ^ Polynomial
					: entry >> 1;
			table[i] = entry;
		}

		return table;
	}

	/// <summary>
	/// Computes the CRC-32/ISO-HDLC checksum over a byte array segment.
	/// Used by the test suite to independently verify the footer CRC without re-reading a file.
	/// </summary>
	/// <param name="data">The byte array to hash.</param>
	/// <param name="length">Number of bytes to include, starting from index 0.</param>
	/// <returns>The CRC-32 checksum as a uint32.</returns>
	internal static uint ComputeCrc32OfArray(byte[] data, int length)
	{
		uint crc = 0xFFFF_FFFFu;

		for (int i = 0; i < length; i++)
			crc = (crc >> 8) ^ Crc32Table[(crc ^ data[i]) & 0xFF];

		return crc ^ 0xFFFF_FFFFu;
	}

	// ────────────────────────────────────────────────────────────────────────
	// Intensity normalization
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Normalizes the fragment intensities of a single precursor entry in-place so that
	/// the most abundant fragment has intensity 1.0. All other intensities are scaled
	/// proportionally.
	///
	/// If the fragment list is empty, or if all intensities are zero or non-finite,
	/// normalization is skipped and the list is left unchanged.
	/// </summary>
	/// <param name="fragments">
	///   The mutable fragment list belonging to a single precursor. Modified in place.
	/// </param>
	private static void NormalizeIntensities(List<MslFragmentIon> fragments)
	{
		if (fragments.Count == 0)
			return;

		// Find the maximum finite, positive intensity across all fragments
		float maxIntensity = 0f;
		foreach (MslFragmentIon f in fragments)
		{
			if (float.IsFinite(f.Intensity) && f.Intensity > maxIntensity)
				maxIntensity = f.Intensity;
		}

		// Skip normalization when max is zero (all-zero or all-NaN intensities)
		if (maxIntensity == 0f)
			return;

		// Divide every intensity by the maximum; result lies in [0, 1]
		float invMax = 1.0f / maxIntensity;
		for (int i = 0; i < fragments.Count; i++)
		{
			MslFragmentIon f = fragments[i];
			f.Intensity = f.Intensity * invMax;
		}
	}

	// ────────────────────────────────────────────────────────────────────────
	// Neutral-loss classification (shared by layout and writer)
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Maps a neutral-loss mass in Daltons to the nearest named
	/// <see cref="MslFormat.NeutralLossCode"/>, or
	/// <see cref="MslFormat.NeutralLossCode.Custom"/> when the mass does not correspond to
	/// any of the seven defined codes.
	///
	/// Tolerance is ±0.01 Da to accommodate minor floating-point variations between tools.
	/// </summary>
	/// <param name="neutralLoss">
	///   Neutral-loss mass in Daltons. 0.0 = no loss. Negative = mass removed from fragment.
	/// </param>
	/// <returns>The best-matching <see cref="MslFormat.NeutralLossCode"/>.</returns>
	internal static MslFormat.NeutralLossCode ClassifyNeutralLoss(double neutralLoss)
	{
		// Named neutral-loss mass constants (negative = lost from the fragment ion)
		const double MassH2O = -18.010565;
		const double MassNH3 = -17.026549;
		const double MassH3PO4 = -97.976895;
		const double MassHPO3 = -79.966331;
		const double MassPlusH2O = MassH3PO4 + MassH2O; // combined loss
		const double Tolerance = 0.01; // Da — intentionally loose for cross-tool compatibility

		if (neutralLoss == 0.0) return MslFormat.NeutralLossCode.None;
		if (Math.Abs(neutralLoss - MassH2O) < Tolerance) return MslFormat.NeutralLossCode.H2O;
		if (Math.Abs(neutralLoss - MassNH3) < Tolerance) return MslFormat.NeutralLossCode.NH3;
		if (Math.Abs(neutralLoss - MassH3PO4) < Tolerance) return MslFormat.NeutralLossCode.H3PO4;
		if (Math.Abs(neutralLoss - MassHPO3) < Tolerance) return MslFormat.NeutralLossCode.HPO3;
		if (Math.Abs(neutralLoss - MassPlusH2O) < Tolerance) return MslFormat.NeutralLossCode.PlusH2O;
		return MslFormat.NeutralLossCode.Custom;
	}

	// ────────────────────────────────────────────────────────────────────────
	// Pass-1 layout engine
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Encapsulates all pre-computed layout information for a single write operation:
	/// section offsets, string table, protein deduplication, elution group assignments,
	/// per-precursor fragment-block offsets, and file-level flags.
	///
	/// Constructing this object constitutes Pass 1 of the two-pass algorithm. After
	/// construction every offset is finalized and can be written sequentially in Pass 2
	/// without any back-patching or Seek() calls.
	/// </summary>
	private sealed class MslWriteLayout
	{
		// ── Counts ───────────────────────────────────────────────────────────

		/// <summary>Total number of precursor entries to write.</summary>
		public int NPrecursors { get; }

		/// <summary>Number of unique proteins in the protein table.</summary>
		public int NProteins { get; }

		/// <summary>Number of unique elution groups (distinct stripped sequences).</summary>
		public int NElutionGroups { get; }

		/// <summary>Total number of unique strings in the string table.</summary>
		public int NStrings => StringList.Count;

		// ── Section offsets (absolute file byte positions) ───────────────────

		/// <summary>Byte offset of the first MslProteinRecord (= MslFormat.HeaderSize).</summary>
		public long ProteinTableOffset { get; }

		/// <summary>Byte offset of the first string table entry.</summary>
		public long StringTableOffset { get; }

		/// <summary>Byte offset of the first MslPrecursorRecord.</summary>
		public long PrecursorSectionOffset { get; }

		/// <summary>Byte offset of the first MslFragmentRecord.</summary>
		public long FragmentSectionOffset { get; }

		/// <summary>Byte offset of the per-precursor offset table (one int64 per precursor).</summary>
		public long OffsetTableOffset { get; }

		/// <summary>
		/// Byte offset immediately after all data to be covered by the CRC-32 checksum.
		/// Equals <see cref="OffsetTableOffset"/>; the offset table and footer are excluded
		/// from the CRC.
		/// </summary>
		public long DataEndOffset => OffsetTableOffset;

		// ── String table ─────────────────────────────────────────────────────

		/// <summary>
		/// Ordered list of all unique strings; index 0 is always the empty string "".
		/// The integer index of each string in this list is what is stored in every
		/// string-index field of the precursor and protein records.
		/// </summary>
		public List<string> StringList { get; }

		/// <summary>
		/// Total number of UTF-8 encoded bytes in all string bodies (excluding the
		/// 4-byte length prefix of each entry). Written as the TotalBytes header in
		/// <see cref="MslWriter.WriteStringTable"/>.
		/// </summary>
		public int StringTableBodyBytes { get; }

		// ── Protein table ────────────────────────────────────────────────────

		/// <summary>
		/// One slot per unique protein accession, in first-occurrence order.
		/// Written verbatim as the protein table section.
		/// </summary>
		public List<ProteinSlot> ProteinSlots { get; }

		// ── Per-precursor layout ─────────────────────────────────────────────

		/// <summary>
		/// Per-precursor layout data (string indices, protein index, fragment-block offset)
		/// in the same order as the original entries list. Indexed by precursor position.
		/// </summary>
		public List<MslPrecursorLayout> PrecursorLayouts { get; }

		// ── File-level flags ─────────────────────────────────────────────────

		/// <summary>
		/// Packed int32 file-level flag field for the header. Computed by inspecting all
		/// entries for ion-mobility, protein data, gene data, and predicted-source presence.
		/// </summary>
		public int FileFlags { get; }

		// ── Constructor (Pass 1) ─────────────────────────────────────────────

		/// <summary>
		/// Performs the complete Pass-1 layout computation over <paramref name="entries"/>:
		/// string internment, protein deduplication, elution-group assignment, per-precursor
		/// m/z sort + intensity normalization, fragment-block offset computation, and
		/// file-level flag evaluation. All results are stored in the public properties above.
		/// </summary>
		/// <param name="entries">
		///   The library entries to be written. Modified in place (fragment lists are sorted
		///   and intensities are normalized). Must not be null.
		/// </param>
		/// <exception cref="NotSupportedException">
		///   Thrown when any fragment's neutral loss maps to
		///   <see cref="MslFormat.NeutralLossCode.Custom"/>, which is not yet supported.
		/// </exception>
		public MslWriteLayout(IReadOnlyList<MslLibraryEntry> entries)
		{
			NPrecursors = entries.Count;

			// ── Step 1: String internment ─────────────────────────────────────
			// Dictionary maps string value → integer index in StringList.
			// Index 0 is reserved for the empty string so that absent optional fields
			// (protein name, gene name, etc.) resolve to a defined empty-string slot.
			var stringIndex = new Dictionary<string, int>(StringComparer.Ordinal);
			StringList = new List<string>();

			// Pre-register the empty string at index 0
			InternString(string.Empty, stringIndex, StringList);

			// Walk every entry and intern all string fields
			foreach (MslLibraryEntry entry in entries)
			{
				InternString(entry.ModifiedSequence ?? string.Empty, stringIndex, StringList);
				InternString(entry.StrippedSequence ?? string.Empty, stringIndex, StringList);
				InternString(entry.ProteinAccession ?? string.Empty, stringIndex, StringList);
				InternString(entry.ProteinName ?? string.Empty, stringIndex, StringList);
				InternString(entry.GeneName ?? string.Empty, stringIndex, StringList);
			}

			// Compute total UTF-8 body bytes (sum of per-string UTF-8 byte counts)
			int bodyBytes = 0;
			foreach (string s in StringList)
				bodyBytes += Encoding.UTF8.GetByteCount(s);
			StringTableBodyBytes = bodyBytes;

			// ── Step 2: Protein deduplication ─────────────────────────────────
			// Map from accession string → protein slot index (zero-based)
			var proteinSlotIndex = new Dictionary<string, int>(StringComparer.Ordinal);
			ProteinSlots = new List<ProteinSlot>();

			foreach (MslLibraryEntry entry in entries)
			{
				string acc = entry.ProteinAccession ?? string.Empty;

				if (!proteinSlotIndex.TryGetValue(acc, out int slotIdx))
				{
					// First time we've seen this accession: create a new slot
					slotIdx = ProteinSlots.Count;
					proteinSlotIndex[acc] = slotIdx;
					ProteinSlots.Add(new ProteinSlot
					{
						Accession = acc,
						AccessionStringIdx = stringIndex[acc],
						NameStringIdx = stringIndex[entry.ProteinName ?? string.Empty],
						GeneStringIdx = stringIndex[entry.GeneName ?? string.Empty],
						PrecursorCount = 0
					});
				}

				// Increment the precursor count for this protein slot
				ProteinSlots[slotIdx].PrecursorCount++;
			}

			NProteins = ProteinSlots.Count;

			// ── Step 3: Elution group assignment ──────────────────────────────
			// Groups entries by StrippedSequence; groups appear in first-occurrence order.
			var elutionGroupIndex = new Dictionary<string, int>(StringComparer.Ordinal);
			int nextGroupId = 0;

			foreach (MslLibraryEntry entry in entries)
			{
				string stripped = entry.StrippedSequence ?? string.Empty;

				if (!elutionGroupIndex.TryGetValue(stripped, out int groupId))
				{
					groupId = nextGroupId++;
					elutionGroupIndex[stripped] = groupId;
				}

				// Assign the computed elution group ID back to the entry
				entry.ElutionGroupId = groupId;
			}

			NElutionGroups = nextGroupId;

			// ── Step 4: Fragment ordering + intensity normalization ─────────────
			foreach (MslLibraryEntry entry in entries)
			{
				if (entry.Fragments == null || entry.Fragments.Count == 0)
					continue;

				// Validate no Custom neutral-loss codes before sorting
				foreach (MslFragmentIon frag in entry.Fragments)
				{
					if (ClassifyNeutralLoss(frag.NeutralLoss) == MslFormat.NeutralLossCode.Custom)
						throw new NotSupportedException(
							$"Fragment in '{entry.ModifiedSequence}' has a custom neutral loss " +
							$"({frag.NeutralLoss:F4} Da). Custom neutral-loss codes are not yet " +
							"supported; the extended annotation table mechanism has not been designed.");
				}

				// Sort fragments by m/z ascending to support binary search during scoring
				entry.Fragments.Sort((a, b) => a.Mz.CompareTo(b.Mz));

				// Normalize intensities so the most abundant fragment = 1.0
				NormalizeIntensities(entry.Fragments);
			}

			// ── Step 5: Offset computation ────────────────────────────────────
			// Compute every absolute file offset in dependency order.

			// Protein table starts immediately after the fixed-size header
			ProteinTableOffset = MslFormat.HeaderSize;

			// String table starts after the protein table
			StringTableOffset = ProteinTableOffset + (long)NProteins * MslFormat.ProteinRecordSize;

			// String table on-disk size: 4-byte NStrings + 4-byte TotalBytes + per-entry data
			// Each entry: 4-byte length prefix + UTF-8 body bytes
			long stringTableTotalBytes = 4L + 4L +  // NStrings header + TotalBytes header
										 (long)StringList.Count * 4 + // length prefixes
										 StringTableBodyBytes;         // body bytes

			// Precursor section starts after the string table
			PrecursorSectionOffset = StringTableOffset + stringTableTotalBytes;

			// Fragment section starts after all precursor records
			FragmentSectionOffset = PrecursorSectionOffset +
									(long)NPrecursors * MslFormat.PrecursorRecordSize;

			// ── Step 6: Per-precursor layout (fragment-block offsets) ──────────
			PrecursorLayouts = new List<MslPrecursorLayout>(NPrecursors);

			// Accumulate the running fragment section offset for each precursor
			long runningFragmentOffset = FragmentSectionOffset;

			for (int i = 0; i < entries.Count; i++)
			{
				MslLibraryEntry entry = entries[i];

				string acc = entry.ProteinAccession ?? string.Empty;
				int proteinIdx = proteinSlotIndex.TryGetValue(acc, out int pIdx) ? pIdx : -1;

				var pl = new MslPrecursorLayout
				{
					ModifiedSeqStringIdx = stringIndex[entry.ModifiedSequence ?? string.Empty],
					StrippedSeqStringIdx = stringIndex[entry.StrippedSequence ?? string.Empty],
					ProteinIdx = proteinIdx,
					FragmentBlockOffset = runningFragmentOffset
				};

				PrecursorLayouts.Add(pl);

				// Advance past this precursor's fragment block
				int fragCount = entry.Fragments?.Count ?? 0;
				runningFragmentOffset += (long)fragCount * MslFormat.FragmentRecordSize;
			}

			// Offset table starts after all fragment blocks
			OffsetTableOffset = runningFragmentOffset;

			// ── Step 7: File-level flags ───────────────────────────────────────
			int flags = 0;

			bool anyIonMobility = false;
			bool anyProteinData = false;
			bool anyGeneData = false;
			bool allPredicted = NPrecursors > 0; // assume true until falsified

			foreach (MslLibraryEntry entry in entries)
			{
				if (entry.IonMobility != 0.0) anyIonMobility = true;
				if (!string.IsNullOrEmpty(entry.ProteinAccession)) anyProteinData = true;
				if (!string.IsNullOrEmpty(entry.GeneName)) anyGeneData = true;
				if (entry.Source != MslFormat.SourceType.Predicted) allPredicted = false;
			}

			if (anyIonMobility) flags |= MslFormat.FileFlagHasIonMobility;
			if (anyProteinData) flags |= MslFormat.FileFlagHasProteinData;
			if (anyGeneData) flags |= MslFormat.FileFlagHasGeneData;
			if (allPredicted) flags |= MslFormat.FileFlagIsPredicted;

			FileFlags = flags;
		}

		// ── String internment helper ─────────────────────────────────────────

		/// <summary>
		/// Adds <paramref name="value"/> to the string table if it is not already present
		/// and returns its assigned integer index. If the string is already interned its
		/// existing index is returned without modification to the table.
		/// </summary>
		/// <param name="value">The string to intern. Must not be null.</param>
		/// <param name="index">The lookup dictionary mapping string → int index.</param>
		/// <param name="list">The ordered string list; the string is appended when new.</param>
		/// <returns>The zero-based integer index assigned to <paramref name="value"/>.</returns>
		private static int InternString(string value,
			Dictionary<string, int> index,
			List<string> list)
		{
			// Fast path: string is already interned
			if (index.TryGetValue(value, out int existingIdx))
				return existingIdx;

			// Slow path: assign the next sequential index
			int newIdx = list.Count;
			list.Add(value);
			index[value] = newIdx;
			return newIdx;
		}
	}
}

// ────────────────────────────────────────────────────────────────────────────
// Supporting data types (file-private)
// ────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Mutable holder for one unique protein in the write layout's protein deduplication table.
/// Created during Pass 1; read during the protein-table write in Pass 2.
/// </summary>
internal sealed class ProteinSlot
{
	/// <summary>UniProt-style accession string, used as the deduplication key.</summary>
	public string Accession { get; init; } = string.Empty;

	/// <summary>String-table index for the protein accession string.</summary>
	public int AccessionStringIdx { get; init; }

	/// <summary>String-table index for the human-readable protein name.</summary>
	public int NameStringIdx { get; init; }

	/// <summary>String-table index for the gene symbol.</summary>
	public int GeneStringIdx { get; init; }

	/// <summary>
	/// Number of precursor entries that reference this protein. Populated during Pass 1
	/// and written to <see cref="MslProteinRecord.NPrecursors"/> in Pass 2.
	/// </summary>
	public int PrecursorCount { get; set; }
}

/// <summary>
/// Pre-computed layout data for a single precursor entry. All values are resolved during
/// Pass 1 so that Pass 2 can write the <see cref="MslPrecursorRecord"/> without any
/// additional computation or dictionary lookups.
/// </summary>
internal sealed class MslPrecursorLayout
{
	/// <summary>
	/// String-table index for the modified sequence (mzLib bracket notation).
	/// Written to <see cref="MslPrecursorRecord.ModifiedSeqStringIdx"/>.
	/// </summary>
	public int ModifiedSeqStringIdx { get; init; }

	/// <summary>
	/// String-table index for the stripped (unmodified) amino-acid sequence.
	/// Written to <see cref="MslPrecursorRecord.StrippedSeqStringIdx"/>.
	/// </summary>
	public int StrippedSeqStringIdx { get; init; }

	/// <summary>
	/// Zero-based index into the protein table for the owning protein.
	/// −1 when the entry has no protein assignment.
	/// Written to <see cref="MslPrecursorRecord.ProteinIdx"/>.
	/// </summary>
	public int ProteinIdx { get; init; }

	/// <summary>
	/// Absolute file byte offset of the first <see cref="MslFragmentRecord"/> for this
	/// precursor. Computed arithmetically in Pass 1 by accumulating fragment block sizes.
	/// Written to <see cref="MslPrecursorRecord.FragmentBlockOffset"/> and to the offset table.
	/// </summary>
	public long FragmentBlockOffset { get; init; }
}