using MassSpectrometry;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using System.Buffers;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using ZstdSharp;
using static UsefulProteomicsDatabases.ProteinDbRetriever;

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
/// When compression is enabled (<c>compressionLevel &gt; 0</c>), all fragment records are
/// first serialized to an in-memory buffer, compressed with zstd, and the compressed size
/// is stored in a 16-byte compression descriptor written immediately after the file header.
/// This pre-compression in memory means the descriptor's sizes are known before the first
/// byte is written to disk, preserving the no-Seek invariant.
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
	// Constants
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// The file magic stored as a uint for use inside struct fields serialized by
	/// <see cref="MemoryMarshal.Write{T}"/> in native (little-endian) byte order.
	/// Produces the canonical on-disk bytes 0x4D 0x5A 0x4C 0x42 ("MZLB").
	/// Do NOT use <see cref="MslFormat.MagicAsUInt32"/> directly in a struct field.
	/// </summary>
	private const uint MagicForLEStruct = 0x424C_5A4Du;

	/// <summary>
	/// Size in bytes of the compression descriptor block written immediately after the
	/// file header when <see cref="MslFormat.FileFlagIsCompressed"/> is set.
	/// Layout: int64 CompressedFragmentSize + int64 UncompressedFragmentSize = 16 bytes.
	/// </summary>
	private const int CompressionDescriptorSize = 16;

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
	///   Ordered list of library entries to write. Must not be null; may be empty
	///   (produces a valid zero-precursor file).
	/// </param>
	/// <param name="compressionLevel">
	///   zstd compression level for the fragment section. 0 = no compression (default);
	///   1–22 = zstd levels in ascending compression ratio / decreasing speed order.
	///   Level 3 is recommended for balanced throughput. When &gt; 0,
	///   <see cref="MslFormat.FileFlagIsCompressed"/> is set in the file header, a
	///   16-byte compression descriptor is inserted at file offset 64, and
	///   <c>MslPrecursorRecord.FragmentBlockOffset</c> values are stored as offsets
	///   into the decompressed fragment buffer rather than absolute file positions.
	///   Note: <see cref="MslLibrary.LoadIndexOnly"/> falls back to full-load for
	///   compressed files; index-only mode is unavailable when the file is compressed.
	/// </param>
	/// <exception cref="ArgumentNullException">
	///   Thrown when <paramref name="outputPath"/> or <paramref name="entries"/> is null.
	/// </exception>
	/// <exception cref="ArgumentOutOfRangeException">
	///   Thrown when <paramref name="compressionLevel"/> is outside [0, 22].
	/// </exception>
	public static void Write(string outputPath, IReadOnlyList<MslLibraryEntry> entries,
		int compressionLevel = 0)
	{
		if (outputPath is null) throw new ArgumentNullException(nameof(outputPath));
		if (entries is null) throw new ArgumentNullException(nameof(entries));
		if (compressionLevel < 0 || compressionLevel > 22)
			throw new ArgumentOutOfRangeException(nameof(compressionLevel),
				"compressionLevel must be in [0, 22]. Use 0 for no compression.");

		string tempPath = outputPath + ".tmp~";

		try
		{
			// ── Pass 1: compute layout ────────────────────────────────────────
			// For compressed files, the fragment records are also serialized and
			// compressed here so that CompressedFragmentSize is known before we
			// open the output file (preserving the no-Seek-in-Pass-2 invariant).
			var layout = new MslWriteLayout(entries, compressionLevel);

			// Pre-compress the fragment data when compression is requested.
			// This gives us the exact compressed size before writing the descriptor.
			byte[]? compressedFragmentData = null;
			if (layout.IsCompressed)
			{
				compressedFragmentData = BuildAndCompressFragmentBuffer(layout, entries);
				// Store sizes in layout so WriteHeader / WriteCompressionDescriptor can use them
				layout.CompressedFragmentSize = compressedFragmentData.Length;
				// UncompressedFragmentSize was already set during Pass 1
			}

			// ── Pass 2a: write all sections except the footer ─────────────────
			using (var fs = new FileStream(tempPath, FileMode.Create, FileAccess.Write,
										   FileShare.None, bufferSize: 65536,
										   FileOptions.SequentialScan))
			using (var writer = new BinaryWriter(fs, Encoding.UTF8, leaveOpen: false))
			{
				WriteHeader(writer, layout);
				if (layout.IsCompressed)
					WriteCompressionDescriptor(writer, layout);
				WriteProteinTable(writer, layout);
				WriteStringTable(writer, layout);
				WritePrecursorArray(writer, layout, entries);

				if (layout.IsCompressed)
					writer.Write(compressedFragmentData!);   // already compressed above
				else
					WriteFragmentRecords(writer, layout, entries);  // uncompressed path

				WriteExtAnnotationTable(writer, layout);
				WriteOffsetTable(writer, layout);
			}

			// ── Pass 2b: CRC + footer ─────────────────────────────────────────
			uint crc = ComputeCrc32(tempPath, layout.DataEndOffset);

			using (var fs = new FileStream(tempPath, FileMode.Append, FileAccess.Write,
										   FileShare.None, bufferSize: 128))
			using (var writer = new BinaryWriter(fs, Encoding.UTF8, leaveOpen: false))
			{
				WriteFooter(writer, layout, crc);
			}

			File.Move(tempPath, outputPath, overwrite: true);
		}
		catch
		{
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
	///   <see cref="Write(string, IReadOnlyList{MslLibraryEntry}, int)"/>.
	/// </param>
	/// <param name="spectra">
	///   Source library spectra. Must not be null.
	/// </param>
	/// <param name="compressionLevel">
	///   zstd compression level. 0 = no compression (default). See
	///   <see cref="Write(string, IReadOnlyList{MslLibraryEntry}, int)"/> for details.
	/// </param>
	/// <exception cref="ArgumentNullException">
	///   Thrown when <paramref name="outputPath"/> or <paramref name="spectra"/> is null.
	/// </exception>
	public static void WriteFromLibrarySpectra(string outputPath,
		IReadOnlyList<LibrarySpectrum> spectra, int compressionLevel = 0)
	{
		if (outputPath is null) throw new ArgumentNullException(nameof(outputPath));
		if (spectra is null) throw new ArgumentNullException(nameof(spectra));

		var entries = new List<MslLibraryEntry>(spectra.Count);
		foreach (LibrarySpectrum spectrum in spectra)
			entries.Add(MslLibraryEntry.FromLibrarySpectrum(spectrum));

		Write(outputPath, entries, compressionLevel);
	}

	/// <summary>
	/// Returns the estimated file size in bytes without writing any data.
	/// Useful for pre-flight disk-space checks before committing to a write.
	/// </summary>
	/// <param name="nPrecursors">Number of precursor entries that will be written.</param>
	/// <param name="avgFragmentsPerPrecursor">Average number of fragment records per precursor.</param>
	/// <param name="avgStringBytes">
	///   Average UTF-8 byte length of each unique string. 20 is a reasonable default.
	/// </param>
	/// <returns>Estimated total file size in bytes.</returns>
	public static long EstimateFileSize(int nPrecursors, int avgFragmentsPerPrecursor,
		int avgStringBytes)
	{
		long estimatedUniqueStrings = (long)(nPrecursors * 3 * 0.6) + 1;
		long stringTableBytes = 8 + estimatedUniqueStrings * (4 + avgStringBytes);
		long fragmentBytes = (long)nPrecursors * avgFragmentsPerPrecursor
										  * MslFormat.FragmentRecordSize;

		return MslFormat.HeaderSize
			 + (long)nPrecursors * MslFormat.ProteinRecordSize
			 + stringTableBytes
			 + (long)nPrecursors * MslFormat.PrecursorRecordSize
			 + fragmentBytes
			 + (long)nPrecursors * sizeof(long)
			 + MslFormat.FooterSize;
	}

	/// <summary>
	/// Validates a list of <see cref="MslLibraryEntry"/> objects before writing.
	/// Returns a list of human-readable error strings. An empty list means all entries
	/// are safe to pass to <see cref="Write"/>.
	/// </summary>
	public static List<string> ValidateEntries(IReadOnlyList<MslLibraryEntry> entries)
	{
		if (entries is null) throw new ArgumentNullException(nameof(entries));

		var errors = new List<string>();

		for (int i = 0; i < entries.Count; i++)
		{
			MslLibraryEntry entry = entries[i];

			if (string.IsNullOrEmpty(entry.ModifiedSequence))
				errors.Add($"Entry [{i}]: ModifiedSequence is null or empty.");

			if (entry.Charge <= 0)
				errors.Add($"Entry [{i}] '{entry.ModifiedSequence}': Charge must be > 0 (got {entry.Charge}).");

			if (entry.PrecursorMz <= 0.0 || !double.IsFinite(entry.PrecursorMz))
				errors.Add($"Entry [{i}] '{entry.ModifiedSequence}': PrecursorMz must be > 0 (got {entry.PrecursorMz}).");

			if (entry.Fragments == null || entry.Fragments.Count == 0)
			{
				errors.Add($"Entry [{i}] '{entry.ModifiedSequence}': FragmentCount must be > 0.");
				continue;
			}

			for (int j = 0; j < entry.Fragments.Count; j++)
			{
				MslFragmentIon frag = entry.Fragments[j];

				if (frag.Mz <= 0f || !float.IsFinite(frag.Mz))
					errors.Add($"Entry [{i}] '{entry.ModifiedSequence}', fragment [{j}]: " +
							   $"Mz must be > 0 (got {frag.Mz}).");

				if (frag.Charge <= 0)
					errors.Add($"Entry [{i}] '{entry.ModifiedSequence}', fragment [{j}]: " +
							   $"Charge must be > 0 (got {frag.Charge}).");

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
	/// Sets <see cref="MslFormat.FileFlagIsCompressed"/> in <c>FileFlags</c> when the
	/// layout uses compression.
	/// </summary>
	private static void WriteHeader(BinaryWriter writer, MslWriteLayout layout)
	{
		var header = new MslFileHeader
		{
			Magic = MagicForLEStruct,
			FormatVersion = MslFormat.CurrentVersion,
			FileFlags = layout.FileFlags,
			NPrecursors = layout.NPrecursors,
			NProteins = layout.NProteins,
			NElutionGroups = layout.NElutionGroups,
			NStrings = layout.NStrings,
			ExtAnnotationTableOffset = layout.HasCustomLosses
										   ? (int)layout.ExtAnnotationTableOffset : 0,
			ProteinTableOffset = layout.ProteinTableOffset,
			StringTableOffset = layout.StringTableOffset,
			PrecursorSectionOffset = layout.PrecursorSectionOffset,
			FragmentSectionOffset = layout.FragmentSectionOffset
		};

		WriteStruct(writer, header);
	}

	/// <summary>
	/// Writes the 16-byte compression descriptor immediately after the file header.
	/// Present only when <see cref="MslFormat.FileFlagIsCompressed"/> is set.
	/// <code>
	///   int64  CompressedFragmentSize    — byte count of the compressed zstd frame
	///   int64  UncompressedFragmentSize  — byte count after decompression
	/// </code>
	/// Both values are known after pre-compression in Pass 1.
	/// </summary>
	private static void WriteCompressionDescriptor(BinaryWriter writer, MslWriteLayout layout)
	{
		writer.Write(layout.CompressedFragmentSize);    // int64
		writer.Write(layout.UncompressedFragmentSize);  // int64
	}

	/// <summary>
	/// Writes one <see cref="MslProteinRecord"/> per unique protein.
	/// </summary>
	private static void WriteProteinTable(BinaryWriter writer, MslWriteLayout layout)
	{
		foreach (ProteinSlot slot in layout.ProteinSlots)
		{
			var record = new MslProteinRecord
			{
				AccessionStringIdx = slot.AccessionStringIdx,
				NameStringIdx = slot.NameStringIdx,
				GeneStringIdx = slot.GeneStringIdx,
				ProteinGroupId = 0,
				NPrecursors = slot.PrecursorCount,
				ProteinFlags = 0
			};

			WriteStruct(writer, record);
		}
	}

	/// <summary>
	/// Writes the string table section.
	/// On-disk: int32 NStrings, int32 TotalBytes, then per-entry: int32 length + UTF-8 body.
	/// The empty string "" is always at index 0.
	/// </summary>
	private static void WriteStringTable(BinaryWriter writer, MslWriteLayout layout)
	{
		writer.Write(layout.StringList.Count);
		writer.Write(layout.StringTableBodyBytes);

		foreach (string s in layout.StringList)
		{
			byte[] encoded = Encoding.UTF8.GetBytes(s);
			writer.Write(encoded.Length);
			writer.Write(encoded);
		}
	}

	/// <summary>
	/// Writes one <see cref="MslPrecursorRecord"/> per entry.
	///
	/// When compressed, <c>FragmentBlockOffset</c> values are decompressed-buffer-relative
	/// (offsets from 0 into the uncompressed fragment buffer, not absolute file positions).
	/// The reader detects compression via <see cref="MslFormat.FileFlagIsCompressed"/> and
	/// interprets offsets accordingly.
	/// </summary>
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
										   ? 0 : entry.StrippedSequence.Length,
				MoleculeType = (short)entry.MoleculeType,
				DissociationType = (short)entry.DissociationType,
				Nce = (short)(entry.Nce * 10),
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
	/// Writes raw <see cref="MslFragmentRecord"/> structs directly to the output stream
	/// (uncompressed path only).
	/// For the compressed path, <see cref="BuildAndCompressFragmentBuffer"/> is used instead.
	/// </summary>
	private static void WriteFragmentRecords(BinaryWriter writer, MslWriteLayout layout,
		IReadOnlyList<MslLibraryEntry> entries)
	{
		foreach (MslLibraryEntry entry in entries)
		{
			foreach (MslFragmentIon frag in entry.Fragments)
			{
				WriteFragmentRecord(writer, layout, frag);
			}
		}
	}

	/// <summary>
	/// Serializes all fragment records into an in-memory buffer, then compresses the buffer
	/// with zstd at <see cref="MslWriteLayout.CompressionLevel"/>. Returns the compressed bytes.
	/// Called during Pass 1 (before the output file is opened) so that
	/// <see cref="MslWriteLayout.CompressedFragmentSize"/> is known when writing the header.
	/// </summary>
	private static byte[] BuildAndCompressFragmentBuffer(MslWriteLayout layout,
		IReadOnlyList<MslLibraryEntry> entries)
	{
		int bufferSize = (int)layout.UncompressedFragmentSize;
		using var ms = new MemoryStream(bufferSize);
		using var writer = new BinaryWriter(ms, Encoding.UTF8, leaveOpen: false);

		foreach (MslLibraryEntry entry in entries)
		{
			foreach (MslFragmentIon frag in entry.Fragments)
			{
				WriteFragmentRecord(writer, layout, frag);
			}
		}

		writer.Flush();
		byte[] uncompressed = ms.ToArray();

		using var compressor = new Compressor(layout.CompressionLevel);
		return compressor.Wrap(uncompressed).ToArray();
	}

	/// <summary>
	/// Writes a single <see cref="MslFragmentRecord"/> to any <see cref="BinaryWriter"/>.
	/// Shared by the uncompressed direct-write path and the compressed buffer-build path.
	/// </summary>
	private static void WriteFragmentRecord(BinaryWriter writer, MslWriteLayout layout,
		MslFragmentIon frag)
	{
		MslFormat.NeutralLossCode lossCode = ClassifyNeutralLoss(frag.NeutralLoss);

		byte flags = MslFormat.EncodeFragmentFlags(
			isInternal: frag.IsInternalFragment,
			isDiagnostic: frag.IsDiagnosticIon,
			lossCode: lossCode,
			excludeFromQuant: frag.ExcludeFromQuant);

		// When neutral_loss_code == Custom, repurpose ResiduePosition as ExtAnnotationIdx
		// (1-based index into the extended annotation table).
		// For all other loss codes, ResiduePosition retains its normal meaning.
		short residuePos = lossCode == MslFormat.NeutralLossCode.Custom
			? (short)layout.CustomLossTable[frag.NeutralLoss]
			: (short)frag.ResiduePosition;

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
			ResiduePosition = residuePos,
			Charge = (byte)frag.Charge,
			Flags = flags
		};

		WriteStruct(writer, record);
	}

	/// <summary>
	/// Writes the extended annotation table section when the library contains at least one
	/// fragment with a custom neutral-loss mass. Does nothing when
	/// <see cref="MslWriteLayout.HasCustomLosses"/> is false.
	///
	/// On-disk format:
	/// <code>
	///   int32    NCustomLosses     — count (≥ 2 when section present: sentinel + entries)
	///   double[] CustomLossMasses  — NCustomLosses × 8 bytes, little-endian IEEE 754
	/// </code>
	/// Index 0 is the reserved sentinel (0.0). Valid custom-loss indices start at 1.
	/// </summary>
	private static void WriteExtAnnotationTable(BinaryWriter writer, MslWriteLayout layout)
	{
		if (!layout.HasCustomLosses)
			return;

		double[] sorted = layout.CustomLossTable
			.OrderBy(kv => kv.Value)
			.Select(kv => kv.Key)
			.ToArray();

		int totalCount = sorted.Length + 1;   // +1 for index-0 sentinel
		writer.Write(totalCount);             // int32: NCustomLosses
		writer.Write(0.0);                    // index 0: reserved sentinel
		foreach (double mass in sorted)
			writer.Write(mass);               // indices 1..N: actual custom masses
	}

	/// <summary>
	/// Writes the offset table: one int64 per precursor giving the fragment-block offset.
	/// For uncompressed files this is the absolute file byte offset.
	/// For compressed files it is the byte offset within the decompressed fragment buffer.
	/// </summary>
	private static void WriteOffsetTable(BinaryWriter writer, MslWriteLayout layout)
	{
		foreach (MslPrecursorLayout pl in layout.PrecursorLayouts)
			writer.Write(pl.FragmentBlockOffset);
	}

	/// <summary>
	/// Writes the 20-byte <see cref="MslFooter"/> as the final bytes of the file.
	/// </summary>
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
	/// Serializes a blittable struct to the writer stream via <see cref="MemoryMarshal.Write{T}"/>.
	/// Rents a buffer from <see cref="ArrayPool{T}"/> to avoid per-call heap allocations.
	/// </summary>
	private static void WriteStruct<T>(BinaryWriter writer, T value) where T : unmanaged
	{
		int size = Marshal.SizeOf<T>();
		byte[] buffer = ArrayPool<byte>.Shared.Rent(size);
		try
		{
			MemoryMarshal.Write(buffer.AsSpan(0, size), ref value);
			writer.Write(buffer, 0, size);
		}
		finally
		{
			ArrayPool<byte>.Shared.Return(buffer);
		}
	}

	// ────────────────────────────────────────────────────────────────────────
	// CRC-32 computation
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Computes the CRC-32/ISO-HDLC checksum over bytes 0..(dataEndOffset−1) of the file.
	/// Reads in 64 KiB chunks. Uses the standard reflected polynomial 0xEDB88320.
	/// </summary>
	private static uint ComputeCrc32(string filePath, long dataEndOffset)
	{
		const int ChunkSize = 65536;
		uint crc = 0xFFFF_FFFFu;
		byte[] buffer = ArrayPool<byte>.Shared.Rent(ChunkSize);

		try
		{
			using var fs = new FileStream(filePath, FileMode.Open, FileAccess.Read,
										  FileShare.None, bufferSize: ChunkSize);
			long remaining = dataEndOffset;

			while (remaining > 0)
			{
				int toRead = (int)Math.Min(remaining, ChunkSize);
				int read = fs.Read(buffer, 0, toRead);
				if (read == 0) break;

				for (int i = 0; i < read; i++)
					crc = (crc >> 8) ^ Crc32Table[(crc ^ buffer[i]) & 0xFF];

				remaining -= read;
			}
		}
		finally
		{
			ArrayPool<byte>.Shared.Return(buffer);
		}

		return crc ^ 0xFFFF_FFFFu;
	}

	private static readonly uint[] Crc32Table = BuildCrc32Table();

	private static uint[] BuildCrc32Table()
	{
		const uint Polynomial = 0xEDB8_8320u;
		var table = new uint[256];

		for (uint i = 0; i < 256; i++)
		{
			uint entry = i;
			for (int bit = 0; bit < 8; bit++)
				entry = (entry & 1u) != 0 ? (entry >> 1) ^ Polynomial : entry >> 1;
			table[i] = entry;
		}

		return table;
	}

	/// <summary>
	/// Computes CRC-32 over a byte array segment.
	/// Used by the test suite to independently verify the footer CRC.
	/// </summary>
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

	private static void NormalizeIntensities(List<MslFragmentIon> fragments)
	{
		if (fragments.Count == 0) return;

		float maxIntensity = 0f;
		foreach (MslFragmentIon f in fragments)
		{
			if (float.IsFinite(f.Intensity) && f.Intensity > maxIntensity)
				maxIntensity = f.Intensity;
		}

		if (maxIntensity == 0f) return;

		float invMax = 1.0f / maxIntensity;
		for (int i = 0; i < fragments.Count; i++)
			fragments[i].Intensity = fragments[i].Intensity * invMax;
	}

	// ────────────────────────────────────────────────────────────────────────
	// Neutral-loss classification
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Maps a neutral-loss mass in Daltons to the nearest named
	/// <see cref="MslFormat.NeutralLossCode"/>, or
	/// <see cref="MslFormat.NeutralLossCode.Custom"/> when the mass does not correspond
	/// to any of the defined codes. Tolerance is ±0.01 Da.
	/// </summary>
	internal static MslFormat.NeutralLossCode ClassifyNeutralLoss(double neutralLoss)
	{
		const double MassH2O = -18.010565;
		const double MassNH3 = -17.026549;
		const double MassH3PO4 = -97.976895;
		const double MassHPO3 = -79.966331;
		const double MassPlusH2O = MassH3PO4 + MassH2O;
		const double Tolerance = 0.01;

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
	/// per-precursor fragment-block offsets, custom neutral-loss table, compression state,
	/// and file-level flags.
	///
	/// Constructing this object constitutes Pass 1 of the two-pass algorithm. After
	/// construction every offset is finalized and can be written sequentially in Pass 2
	/// without any back-patching or Seek() calls.
	/// </summary>
	private sealed class MslWriteLayout
	{
		// ── Counts ───────────────────────────────────────────────────────────
		public int NPrecursors { get; }
		public int NProteins { get; private set; }
		public int NElutionGroups { get; private set; }
		public int NStrings => StringList.Count;

		// ── Section offsets (absolute file byte positions) ───────────────────

		/// <summary>Absolute byte offset of the first MslProteinRecord.</summary>
		public long ProteinTableOffset { get; private set; }

		/// <summary>Absolute byte offset of the first string table entry.</summary>
		public long StringTableOffset { get; private set; }

		/// <summary>Absolute byte offset of the first MslPrecursorRecord.</summary>
		public long PrecursorSectionOffset { get; private set; }

		/// <summary>
		/// Absolute byte offset in the file where fragment data begins.
		/// For uncompressed files: offset of the first MslFragmentRecord.
		/// For compressed files: offset of the start of the zstd frame.
		/// Written into MslFileHeader.FragmentSectionOffset so the reader can locate it.
		/// </summary>
		public long FragmentSectionOffset { get; private set; }

		/// <summary>
		/// Absolute byte offset of the per-precursor offset table.
		/// For uncompressed files: immediately after the fragment section (+ ext table if present).
		/// For compressed files: computed lazily as FragmentSectionOffset + CompressedFragmentSize
		/// + ext annotation table size, once CompressedFragmentSize is set.
		/// </summary>
		public long OffsetTableOffset { get; private set; }

		/// <summary>CRC covers bytes [0, DataEndOffset).</summary>
		public long DataEndOffset => OffsetTableOffset;

		// ── String table ─────────────────────────────────────────────────────
		public List<string> StringList { get; private set; }
		public int StringTableBodyBytes { get; private set; }

		// ── Protein table ────────────────────────────────────────────────────
		public List<ProteinSlot> ProteinSlots { get; private set; }

		// ── Per-precursor layout ─────────────────────────────────────────────
		public List<MslPrecursorLayout> PrecursorLayouts { get; private set; }

		// ── File-level flags ─────────────────────────────────────────────────
		public int FileFlags { get; private set; }

		// ── Extended annotation table ────────────────────────────────────────
		public bool HasCustomLosses { get; private set; }
		public Dictionary<double, int> CustomLossTable { get; private set; } = new();

		/// <summary>
		/// Absolute byte offset of the extended annotation table in the file.
		/// 0 when <see cref="HasCustomLosses"/> is false.
		/// For compressed files this is after the compressed fragment data.
		/// </summary>
		public long ExtAnnotationTableOffset { get; private set; }

		// ── Compression ──────────────────────────────────────────────────────

		/// <summary>zstd compression level supplied to the constructor; 0 = no compression.</summary>
		public int CompressionLevel { get; }

		/// <summary>True when <see cref="CompressionLevel"/> is greater than 0.</summary>
		public bool IsCompressed => CompressionLevel > 0;

		/// <summary>
		/// Byte count of the compressed zstd frame.
		/// Set by <see cref="MslWriter.Write"/> after <see cref="BuildAndCompressFragmentBuffer"/>
		/// returns; used by <see cref="WriteCompressionDescriptor"/> and, for compressed files,
		/// to compute <see cref="OffsetTableOffset"/>. 0 until set.
		/// </summary>
		public long CompressedFragmentSize
		{
			get => _compressedFragmentSize;
			set
			{
				_compressedFragmentSize = value;
				if (IsCompressed)
					UpdateCompressedOffsets();
			}
		}
		private long _compressedFragmentSize;

		/// <summary>
		/// Total byte count of the uncompressed fragment records.
		/// Computed during Pass 1 from total fragment counts × FragmentRecordSize.
		/// Written into the compression descriptor alongside CompressedFragmentSize.
		/// </summary>
		public long UncompressedFragmentSize { get; private set; }

		// ── Constructor (Pass 1) ─────────────────────────────────────────────

		/// <summary>
		/// Performs the complete Pass-1 layout computation: string internment, protein
		/// deduplication, elution-group assignment, per-precursor m/z sort + intensity
		/// normalization, custom neutral-loss table construction, section offset computation,
		/// and file-level flag evaluation.
		/// </summary>
		/// <param name="entries">
		///   Library entries to be written. Modified in place (fragment lists are sorted
		///   and intensities are normalized). Must not be null.
		/// </param>
		/// <param name="compressionLevel">
		///   zstd level; 0 = no compression. When &gt; 0 the compression descriptor (16 bytes)
		///   is included in section offset calculations.
		/// </param>
		public MslWriteLayout(IReadOnlyList<MslLibraryEntry> entries, int compressionLevel)
		{
			NPrecursors = entries.Count;
			CompressionLevel = compressionLevel;

			// ── Step 1: String internment ─────────────────────────────────────
			var stringIndex = new Dictionary<string, int>(StringComparer.Ordinal);
			StringList = new List<string>();
			InternString(string.Empty, stringIndex, StringList);

			foreach (MslLibraryEntry entry in entries)
			{
				InternString(entry.ModifiedSequence ?? string.Empty, stringIndex, StringList);
				InternString(entry.StrippedSequence ?? string.Empty, stringIndex, StringList);
				InternString(entry.ProteinAccession ?? string.Empty, stringIndex, StringList);
				InternString(entry.ProteinName ?? string.Empty, stringIndex, StringList);
				InternString(entry.GeneName ?? string.Empty, stringIndex, StringList);
			}

			int bodyBytes = 0;
			foreach (string s in StringList)
				bodyBytes += Encoding.UTF8.GetByteCount(s);
			StringTableBodyBytes = bodyBytes;

			// ── Step 2: Protein deduplication ─────────────────────────────────
			var proteinSlotIndex = new Dictionary<string, int>(StringComparer.Ordinal);
			ProteinSlots = new List<ProteinSlot>();

			foreach (MslLibraryEntry entry in entries)
			{
				string acc = entry.ProteinAccession ?? string.Empty;

				if (!proteinSlotIndex.TryGetValue(acc, out int slotIdx))
				{
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

				ProteinSlots[slotIdx].PrecursorCount++;
			}

			NProteins = ProteinSlots.Count;

			// ── Step 3: Elution group assignment ──────────────────────────────
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

				entry.ElutionGroupId = groupId;
			}

			NElutionGroups = nextGroupId;

			// ── Step 4: Fragment ordering + normalization + custom-loss table ──
			CustomLossTable = new Dictionary<double, int>();
			int nextCustomIdx = 1; // index 0 is reserved as the 0.0 sentinel

			foreach (MslLibraryEntry entry in entries)
			{
				if (entry.Fragments == null || entry.Fragments.Count == 0)
					continue;

				foreach (MslFragmentIon frag in entry.Fragments)
				{
					if (ClassifyNeutralLoss(frag.NeutralLoss) == MslFormat.NeutralLossCode.Custom)
					{
						if (!CustomLossTable.ContainsKey(frag.NeutralLoss))
							CustomLossTable[frag.NeutralLoss] = nextCustomIdx++;
						HasCustomLosses = true;
					}
				}

				entry.Fragments.Sort((a, b) => a.Mz.CompareTo(b.Mz));
				NormalizeIntensities(entry.Fragments);
			}

			// ── Step 5: Offset computation ────────────────────────────────────
			// When compressed, the 16-byte descriptor is inserted at offset 64, pushing
			// the protein table (and everything after it) up by CompressionDescriptorSize.
			long headerEnd = MslFormat.HeaderSize
				+ (IsCompressed ? CompressionDescriptorSize : 0L);

			ProteinTableOffset = headerEnd;

			StringTableOffset = ProteinTableOffset
				+ (long)NProteins * MslFormat.ProteinRecordSize;

			long stringTableTotalBytes = 4L + 4L
				+ (long)StringList.Count * 4
				+ StringTableBodyBytes;

			PrecursorSectionOffset = StringTableOffset + stringTableTotalBytes;

			FragmentSectionOffset = PrecursorSectionOffset
				+ (long)NPrecursors * MslFormat.PrecursorRecordSize;

			// ── Step 6: Per-precursor layout (fragment-block offsets) ──────────
			PrecursorLayouts = new List<MslPrecursorLayout>(NPrecursors);

			// Uncompressed: FragmentBlockOffset = absolute file position.
			// Compressed:   FragmentBlockOffset = offset within decompressed buffer (starts at 0).
			long runningOffset = IsCompressed ? 0L : FragmentSectionOffset;

			for (int i = 0; i < entries.Count; i++)
			{
				MslLibraryEntry entry = entries[i];
				string acc = entry.ProteinAccession ?? string.Empty;
				int proteinIdx = proteinSlotIndex.TryGetValue(acc, out int pIdx) ? pIdx : -1;

				PrecursorLayouts.Add(new MslPrecursorLayout
				{
					ModifiedSeqStringIdx = stringIndex[entry.ModifiedSequence ?? string.Empty],
					StrippedSeqStringIdx = stringIndex[entry.StrippedSequence ?? string.Empty],
					ProteinIdx = proteinIdx,
					FragmentBlockOffset = runningOffset
				});

				int fragCount = entry.Fragments?.Count ?? 0;
				runningOffset += (long)fragCount * MslFormat.FragmentRecordSize;
			}

			// Total uncompressed fragment bytes
			UncompressedFragmentSize = IsCompressed
				? runningOffset                          // runningOffset started at 0
				: runningOffset - FragmentSectionOffset;

			// For uncompressed files, compute OffsetTableOffset now (it's fully known).
			// For compressed files, OffsetTableOffset is computed lazily in
			// UpdateCompressedOffsets() once CompressedFragmentSize is set.
			if (!IsCompressed)
			{
				long fragmentSectionEnd = runningOffset; // = absolute end of fragment section
				if (HasCustomLosses)
				{
					ExtAnnotationTableOffset = fragmentSectionEnd;
					long extAnnotSize = 4L + (long)(CustomLossTable.Count + 1) * sizeof(double);
					OffsetTableOffset = ExtAnnotationTableOffset + extAnnotSize;
				}
				else
				{
					ExtAnnotationTableOffset = 0;
					OffsetTableOffset = fragmentSectionEnd;
				}
			}
			// (For compressed files, OffsetTableOffset remains 0 until UpdateCompressedOffsets)

			// ── Step 7: File-level flags ───────────────────────────────────────
			int flags = 0;
			bool anyIonMobility = false;
			bool anyProteinData = false;
			bool anyGeneData = false;
			bool allPredicted = NPrecursors > 0;

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
			if (HasCustomLosses) flags |= MslFormat.FileFlagHasExtAnnotations;
			if (IsCompressed) flags |= MslFormat.FileFlagIsCompressed;

			FileFlags = flags;
		}

		/// <summary>
		/// Computes <see cref="OffsetTableOffset"/> and <see cref="ExtAnnotationTableOffset"/>
		/// for compressed files, once <see cref="CompressedFragmentSize"/> is known.
		/// Called automatically when the <see cref="CompressedFragmentSize"/> property is set.
		/// </summary>
		private void UpdateCompressedOffsets()
		{
			// For compressed files, the on-disk layout after the precursor section is:
			//   [Compressed Fragment Data]  CompressedFragmentSize bytes
			//   [Extended Annotation Table] present when HasCustomLosses
			//   [Offset Table]
			long compressedEnd = FragmentSectionOffset + _compressedFragmentSize;

			if (HasCustomLosses)
			{
				ExtAnnotationTableOffset = compressedEnd;
				long extAnnotSize = 4L + (long)(CustomLossTable.Count + 1) * sizeof(double);
				OffsetTableOffset = ExtAnnotationTableOffset + extAnnotSize;
			}
			else
			{
				ExtAnnotationTableOffset = 0;
				OffsetTableOffset = compressedEnd;
			}
		}

		// ── String internment helper ─────────────────────────────────────────

		private static int InternString(string value,
			Dictionary<string, int> index, List<string> list)
		{
			if (index.TryGetValue(value, out int existingIdx))
				return existingIdx;

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
	public string Accession { get; init; } = string.Empty;
	public int AccessionStringIdx { get; init; }
	public int NameStringIdx { get; init; }
	public int GeneStringIdx { get; init; }
	public int PrecursorCount { get; set; }
}

/// <summary>
/// Pre-computed layout data for a single precursor entry. All values are resolved during
/// Pass 1 so that Pass 2 can write the <see cref="MslPrecursorRecord"/> without any
/// additional computation or dictionary lookups.
/// </summary>
internal sealed class MslPrecursorLayout
{
	public int ModifiedSeqStringIdx { get; init; }
	public int StrippedSeqStringIdx { get; init; }
	public int ProteinIdx { get; init; }

	/// <summary>
	/// Absolute file byte offset (uncompressed) or decompressed-buffer offset (compressed)
	/// of the first <see cref="MslFragmentRecord"/> for this precursor.
	/// </summary>
	public long FragmentBlockOffset { get; init; }
}