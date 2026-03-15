using System.Buffers;
using System.Runtime.InteropServices;
using System.Text;
using MassSpectrometry;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using ZstdSharp;

namespace Readers.SpectralLibrary;

/// <summary>
/// Static reader for the .msl (mzLib Spectral Library) binary format.
///
/// Two read modes are supported:
///
///   <see cref="Load"/> — full load: reads the entire file into memory. All precursors,
///   fragment blocks, and strings are deserialized into <see cref="MslLibraryEntry"/>
///   objects. Best for libraries that fit comfortably in RAM (typically &lt;~2 GB). No
///   <see cref="System.IO.FileStream"/> is held open after this method returns.
///
///   <see cref="LoadIndexOnly"/> — index-only load: reads only the precursor records and
///   string table. Fragment blocks remain on disk and are fetched lazily via seeks into
///   the kept-open <see cref="System.IO.FileStream"/>. Best for multi-gigabyte libraries
///   where loading all fragments would exhaust available RAM.
///   Note: compressed files always fall back to full-load regardless of which method is
///   called — index-only mode is unavailable when <see cref="MslFormat.FileFlagIsCompressed"/>
///   is set.
///
/// Both modes return an <see cref="MslLibrary"/> with the same public API; the difference
/// is purely internal.
///
/// Validation sequence on every open:
///   1. Leading magic bytes verified against "MZLB".
///   2. <c>FormatVersion</c> accepted when in range [1, <see cref="MslFormat.CurrentVersion"/>].
///   3. Trailing footer magic verified.
///   4. Footer <c>NPrecursors</c> cross-checked against header <c>NPrecursors</c>.
///   5. CRC-32/ISO-HDLC checksum verified over all data bytes before the offset table.
/// </summary>
public static class MslReader
{
	// ── Static constructor ────────────────────────────────────────────────────

	/// <summary>
	/// Runs <see cref="MslStructs.SizeCheck"/> and builds the CRC-32 lookup table once
	/// when the class is first used. A Pack-setting mistake causes an immediate throw rather
	/// than a silent misread of every file.
	/// </summary>
	static MslReader()
	{
		MslStructs.SizeCheck();
		BuildCrc32Table();
	}

	// ── CRC-32/ISO-HDLC ──────────────────────────────────────────────────────

	/// <summary>
	/// Precomputed 256-entry lookup table for CRC-32/ISO-HDLC (reflected polynomial
	/// 0xEDB88320). Matches the output of zlib <c>crc32()</c> and PKZIP.
	/// Populated once in the static constructor by <see cref="BuildCrc32Table"/>.
	/// </summary>
	private static readonly uint[] Crc32Table = new uint[256];

	/// <summary>
	/// Fills <see cref="Crc32Table"/> using the reflected CRC-32/ISO-HDLC polynomial
	/// (0xEDB88320). Called exactly once from the static constructor.
	/// </summary>
	private static void BuildCrc32Table()
	{
		// Reflected CRC-32 polynomial: same constant as zlib, PKZIP, and Ethernet.
		const uint Poly = 0xEDB8_8320u;

		for (uint i = 0; i < 256; i++)
		{
			uint entry = i;

			// Shift 8 times, XOR with polynomial whenever the output LSB is 1
			for (int bit = 0; bit < 8; bit++)
				entry = (entry & 1u) != 0 ? (entry >> 1) ^ Poly : entry >> 1;

			Crc32Table[i] = entry;
		}
	}

	/// <summary>
	/// Computes the CRC-32/ISO-HDLC checksum over the first <paramref name="length"/> bytes
	/// of <paramref name="data"/>. The result matches <c>zlib crc32()</c>.
	/// </summary>
	/// <param name="data">Source byte array. Must not be null.</param>
	/// <param name="length">
	/// Number of bytes to include starting from index 0.
	/// Must satisfy 0 &lt;= length &lt;= data.Length.
	/// </param>
	/// <returns>The 32-bit CRC checksum.</returns>
	internal static uint ComputeCrc32OfArray(byte[] data, int length)
	{
		// CRC register starts all-ones; final result is the one's complement.
		uint crc = 0xFFFF_FFFFu;

		for (int i = 0; i < length; i++)
			crc = (crc >> 8) ^ Crc32Table[(crc ^ data[i]) & 0xFFu];

		return crc ^ 0xFFFF_FFFFu;
	}

	// ── Neutral-loss decoding ─────────────────────────────────────────────────

	/// <summary>
	/// Maps a <see cref="MslFormat.NeutralLossCode"/> to the corresponding neutral-loss mass
	/// in daltons. Returns 0.0 for <see cref="MslFormat.NeutralLossCode.None"/> and for any
	/// unrecognised code. Custom losses are not decoded here; callers must handle them via the
	/// extended annotation table.
	/// </summary>
	/// <param name="code">Neutral-loss code extracted from the fragment flags byte.</param>
	/// <returns>Neutral-loss mass in daltons (negative = loss); 0.0 when absent or unknown.</returns>
	private static double DecodeNeutralLoss(MslFormat.NeutralLossCode code) => code switch
	{
		MslFormat.NeutralLossCode.None => 0.0,
		MslFormat.NeutralLossCode.H2O => -18.010565,
		MslFormat.NeutralLossCode.NH3 => -17.026549,
		MslFormat.NeutralLossCode.H3PO4 => -97.976895,
		MslFormat.NeutralLossCode.HPO3 => -79.966331,
		MslFormat.NeutralLossCode.PlusH2O => -18.010565 + -97.976895,
		_ => 0.0
	};

	// ── Public API ────────────────────────────────────────────────────────────

	/// <summary>
	/// Reads only the <see cref="MslFileHeader"/> (64 bytes at offset 0) without loading any
	/// entries, string table, or fragments. Useful for rapid metadata inspection.
	///
	/// The leading magic is verified; the full five-step validation sequence is NOT performed.
	/// </summary>
	/// <param name="filePath">Absolute or relative path to the .msl file. Must not be null.</param>
	/// <returns>The deserialized <see cref="MslFileHeader"/> struct.</returns>
	/// <exception cref="FileNotFoundException">File does not exist at <paramref name="filePath"/>.</exception>
	/// <exception cref="FormatException">
	/// The first four bytes do not match the MSL magic, or the file is too short to hold a header.
	/// </exception>
	public static MslFileHeader ReadHeaderOnly(string filePath)
	{
		if (!File.Exists(filePath))
			throw new FileNotFoundException($"MSL file not found: '{filePath}'.", filePath);

		// Only read the header bytes — no need to load the whole file for metadata inspection
		using var fs = new FileStream(filePath, FileMode.Open, FileAccess.Read,
									  FileShare.Read, bufferSize: MslFormat.HeaderSize);

		byte[] headerBytes = new byte[MslFormat.HeaderSize];
		int n = fs.Read(headerBytes, 0, MslFormat.HeaderSize);

		if (n < 4)
			throw new FormatException($"File '{filePath}' is too short to contain a valid MSL header.");

		// Verify magic against raw bytes — MagicMatches handles the LE struct byte-swap correctly
		if (!MslFormat.MagicMatches(headerBytes.AsSpan(0, 4)))
			throw new FormatException(
				$"Magic mismatch in '{filePath}': not an MSL file " +
				$"(got 0x{headerBytes[0]:X2}{headerBytes[1]:X2}{headerBytes[2]:X2}{headerBytes[3]:X2}).");

		if (n < MslFormat.HeaderSize)
			throw new FormatException(
				$"File '{filePath}' is too short to contain a complete 64-byte MSL header.");

		return MemoryMarshal.Read<MslFileHeader>(headerBytes.AsSpan());
	}

	/// <summary>
	/// Reads the entire .msl file into memory and returns an <see cref="MslLibraryData"/> with
	/// all precursor entries and fragment ions fully loaded. No file handle is held open after
	/// this method returns.
	///
	/// When <see cref="MslFormat.FileFlagIsCompressed"/> is set the fragment section is
	/// decompressed from the zstd frame before fragment blocks are read. All other read paths
	/// are unchanged.
	///
	/// Optimal for libraries that fit comfortably in RAM (typically &lt;~2 GB). For larger
	/// libraries use <see cref="LoadIndexOnly"/> to avoid exhausting available RAM.
	/// </summary>
	/// <param name="filePath">Path to the .msl file. Must not be null.</param>
	/// <returns>
	/// A fully populated <see cref="MslLibraryData"/> whose <c>Entries</c> list contains one
	/// <see cref="MslLibraryEntry"/> per precursor, all with fragment ions loaded.
	/// </returns>
	/// <exception cref="FileNotFoundException">File does not exist.</exception>
	/// <exception cref="FormatException">
	/// Magic mismatch, unsupported version, trailing magic mismatch, or NPrecursors mismatch.
	/// </exception>
	/// <exception cref="InvalidDataException">CRC-32 checksum mismatch (data corruption).</exception>
	public static MslLibraryData Load(string filePath)
	{
		if (!File.Exists(filePath))
			throw new FileNotFoundException($"MSL file not found: '{filePath}'.", filePath);

		// Read the entire file into memory for CRC validation and zero-copy struct casting
		byte[] fileBytes = File.ReadAllBytes(filePath);

		ValidateFileBytes(fileBytes, filePath, out MslFileHeader header, out _);

		// Deserialise all three lookup tables
		string[] strings = ReadStringTable(fileBytes, header);
		MslProteinRecord[] proteins = ReadProteinTable(fileBytes, header);
		MslPrecursorRecord[] precursors = ReadPrecursorArray(fileBytes, header);

		// Read extended annotation table (custom neutral-loss masses) when present
		double[] customLossMasses = ReadExtAnnotationTable(fileBytes, header);

		// When compressed, decompress the fragment section into a separate buffer.
		// FragmentBlockOffset values in precursor records are then relative to this buffer.
		bool isCompressed = (header.FileFlags & MslFormat.FileFlagIsCompressed) != 0;
		byte[] fragmentBuffer = isCompressed
			? DecompressFragmentSection(fileBytes, header)
			: fileBytes;

		// Build the fully-loaded entry list — each entry includes its fragment ions
		var entries = new List<MslLibraryEntry>(precursors.Length);

		for (int i = 0; i < precursors.Length; i++)
		{
			MslPrecursorRecord p = precursors[i];
			List<MslFragmentIon> fragments = ReadFragmentBlockFromBytes(fragmentBuffer, p, customLossMasses);
			entries.Add(ConvertPrecursor(p, strings, proteins, fragments));
		}

		return new MslLibraryData(entries, header);
	}

	/// <summary>
	/// Reads only the precursor records and string table into memory; fragment blocks remain
	/// on disk and are fetched lazily via seeks into the kept-open <see cref="System.IO.FileStream"/>.
	///
	/// The returned <see cref="MslLibraryData"/> holds an open file handle until its
	/// <see cref="MslLibrary.Dispose"/> method is called. Callers must dispose the library
	/// when finished (critical on Windows, which uses mandatory file locking).
	///
	/// All five validation checks are performed on open before any entries are returned.
	///
	/// <para>
	/// <b>Compressed files:</b> when <see cref="MslFormat.FileFlagIsCompressed"/> is set,
	/// index-only mode is not available because fragment block offsets are relative to the
	/// decompressed buffer rather than the file. This method transparently falls back to full
	/// decompression and returns an <see cref="MslLibraryData"/> with <c>IsIndexOnly = false</c>.
	/// No exception is thrown; callers can detect this via <see cref="MslLibrary.IsIndexOnly"/>.
	/// </para>
	/// </summary>
	/// <param name="filePath">Path to the .msl file. Must not be null.</param>
	/// <returns>
	/// An <see cref="MslLibraryData"/> in index-only mode for uncompressed files, or in
	/// full-load mode for compressed files.
	/// </returns>
	/// <exception cref="FileNotFoundException">File does not exist.</exception>
	/// <exception cref="FormatException">Structural or version validation failed.</exception>
	/// <exception cref="InvalidDataException">CRC-32 checksum mismatch.</exception>
	public static MslLibraryData LoadIndexOnly(string filePath)
	{
		if (!File.Exists(filePath))
			throw new FileNotFoundException($"MSL file not found: '{filePath}'.", filePath);

		// Read the whole file for CRC validation and struct deserialization;
		// a separate persistent FileStream is then opened for on-demand fragment reads.
		byte[] fileBytes = File.ReadAllBytes(filePath);

		ValidateFileBytes(fileBytes, filePath, out MslFileHeader header, out _);

		// Compressed files cannot use index-only mode: fragment offsets are decompressed-
		// buffer-relative and there is no persistent decompressed buffer to seek into.
		// Fall back to full-load transparently.
		bool isCompressed = (header.FileFlags & MslFormat.FileFlagIsCompressed) != 0;
		if (isCompressed)
		{
			System.Diagnostics.Debug.WriteLine(
				$"[MslReader] LoadIndexOnly called on compressed file '{filePath}'; " +
				"falling back to full-load (index-only mode unavailable for compressed files).");

			string[] stringsC = ReadStringTable(fileBytes, header);
			MslProteinRecord[] proteinsC = ReadProteinTable(fileBytes, header);
			MslPrecursorRecord[] precursorsC = ReadPrecursorArray(fileBytes, header);
			double[] customLossMassesC = ReadExtAnnotationTable(fileBytes, header);
			byte[] fragmentBuffer = DecompressFragmentSection(fileBytes, header);

			var entriesC = new List<MslLibraryEntry>(precursorsC.Length);
			for (int i = 0; i < precursorsC.Length; i++)
			{
				MslPrecursorRecord p = precursorsC[i];
				List<MslFragmentIon> fragments =
					ReadFragmentBlockFromBytes(fragmentBuffer, p, customLossMassesC);
				entriesC.Add(ConvertPrecursor(p, stringsC, proteinsC, fragments));
			}

			// Return a full-load library (no open stream, IsIndexOnly = false)
			return new MslLibraryData(entriesC, header);
		}

		// Uncompressed: standard index-only path
		string[] strings = ReadStringTable(fileBytes, header);
		MslProteinRecord[] proteins = ReadProteinTable(fileBytes, header);
		MslPrecursorRecord[] precursors = ReadPrecursorArray(fileBytes, header);

		// Read the extended annotation table now so on-demand fragment reads can decode
		// custom neutral-loss indices without re-reading the file header each time.
		double[] customLossMasses = ReadExtAnnotationTable(fileBytes, header);

		// Build skeleton entries — fragment lists are intentionally left empty
		var entries = new List<MslLibraryEntry>(precursors.Length);

		for (int i = 0; i < precursors.Length; i++)
			entries.Add(ConvertPrecursor(precursors[i], strings, proteins, new List<MslFragmentIon>()));

		// Open a persistent read stream for on-demand fragment reads.
		// FileShare.Read allows concurrent readers; prevents writers from modifying the open file.
		var onDemandStream = new FileStream(
			filePath, FileMode.Open, FileAccess.Read, FileShare.Read, bufferSize: 4096);

		return new MslLibraryData(entries, header, precursors, strings, proteins, onDemandStream, customLossMasses);
	}

	// ── Validation ────────────────────────────────────────────────────────────

	/// <summary>
	/// Performs all five mandatory validation checks on the in-memory file bytes and
	/// outputs the deserialized header and footer on success. Throws on the first failure.
	///
	/// Version acceptance policy:
	/// <list type="bullet">
	///   <item>Versions 1 through <see cref="MslFormat.CurrentVersion"/> — accepted.</item>
	///   <item>Any other version — <see cref="FormatException"/> thrown.</item>
	/// </list>
	/// </summary>
	/// <param name="fileBytes">Complete file content already loaded into memory.</param>
	/// <param name="filePath">Original file path, used only for exception messages.</param>
	/// <param name="header">Output: deserialized <see cref="MslFileHeader"/> on success.</param>
	/// <param name="footer">Output: deserialized <see cref="MslFooter"/> on success.</param>
	/// <exception cref="FormatException">Any structural or version check fails.</exception>
	/// <exception cref="InvalidDataException">CRC-32 does not match the stored value.</exception>
	private static void ValidateFileBytes(
		byte[] fileBytes,
		string filePath,
		out MslFileHeader header,
		out MslFooter footer)
	{
		// ── Check 1: minimum size ─────────────────────────────────────────────
		int minimumSize = MslFormat.HeaderSize + MslFormat.FooterSize;

		if (fileBytes.Length < minimumSize)
			throw new FormatException(
				$"File '{filePath}' is too short ({fileBytes.Length} bytes) " +
				$"to be a valid .msl file (minimum {minimumSize} bytes).");

		// ── Check 2: leading magic ────────────────────────────────────────────
		// MagicMatches() operates on raw bytes; it is endian-safe and handles the
		// LE struct byte-swap applied by the writer (see Prompt 2 Handoff §5.1).
		if (!MslFormat.MagicMatches(fileBytes.AsSpan(0, 4)))
			throw new FormatException(
				$"Magic mismatch in '{filePath}': not an MSL file " +
				$"(got 0x{fileBytes[0]:X2}{fileBytes[1]:X2}{fileBytes[2]:X2}{fileBytes[3]:X2}).");

		// ── Check 3: format version ───────────────────────────────────────────
		// Deserialise the header now; re-use it for all subsequent field reads.
		header = MemoryMarshal.Read<MslFileHeader>(fileBytes.AsSpan(0, MslFormat.HeaderSize));

		// Accept versions 1 through CurrentVersion; reject anything outside that range.
		if (header.FormatVersion < 1 || header.FormatVersion > MslFormat.CurrentVersion)
			throw new FormatException(
				$"Unsupported version: {header.FormatVersion} in '{filePath}'. " +
				$"This reader supports versions 1–{MslFormat.CurrentVersion}.");

		// ── Check 4: trailing footer magic ────────────────────────────────────
		int footerStart = fileBytes.Length - MslFormat.FooterSize;
		footer = MemoryMarshal.Read<MslFooter>(fileBytes.AsSpan(footerStart, MslFormat.FooterSize));

		// The trailing magic occupies the last 4 bytes of the file.
		// Verify using MagicMatches() on the raw bytes for endian-safe comparison.
		int trailingMagicStart = fileBytes.Length - 4;

		if (!MslFormat.MagicMatches(fileBytes.AsSpan(trailingMagicStart, 4)))
			throw new FormatException(
				$"Trailing magic mismatch in '{filePath}': file is truncated or not a valid .msl file.");

		// ── Check 5: NPrecursors cross-check ──────────────────────────────────
		if (footer.NPrecursors != header.NPrecursors)
			throw new FormatException(
				$"NPrecursors mismatch in '{filePath}': " +
				$"header says {header.NPrecursors}, footer says {footer.NPrecursors}. " +
				"File may be truncated or corrupt.");

		// ── Check 6: CRC-32 checksum ──────────────────────────────────────────
		// Coverage: bytes 0..(OffsetTableOffset - 1); excludes the offset table and footer.
		long crcEndOffset = footer.OffsetTableOffset;

		if (crcEndOffset < 0 || crcEndOffset > fileBytes.Length)
			throw new FormatException(
				$"Invalid OffsetTableOffset ({crcEndOffset}) in footer of '{filePath}'.");

		uint computedCrc = ComputeCrc32OfArray(fileBytes, (int)crcEndOffset);

		if (computedCrc != footer.DataCrc32)
			throw new InvalidDataException(
				$"CRC32 mismatch: file may be corrupted. " +
				$"Stored: 0x{footer.DataCrc32:X8}, Computed: 0x{computedCrc:X8}.");
	}

	// ── Table readers ─────────────────────────────────────────────────────────

	/// <summary>
	/// Deserialises the string table from the in-memory file bytes starting at
	/// <c>header.StringTableOffset</c>. Returns a string array where index 0 is always
	/// the empty string.
	///
	/// On-disk layout per entry: [int32 length][UTF-8 body, no null terminator].
	/// The section header contains [int32 NStrings][int32 TotalBodyBytes].
	/// </summary>
	/// <param name="fileBytes">Complete file content.</param>
	/// <param name="header">Deserialized file header providing offset and count.</param>
	/// <returns>
	/// All interned strings as a zero-based string array; length equals the NStrings
	/// field found at the start of the string table section.
	/// </returns>
	/// <exception cref="FormatException">Index 0 is not the empty string.</exception>
	private static string[] ReadStringTable(byte[] fileBytes, MslFileHeader header)
	{
		// Walking position within fileBytes
		int pos = (int)header.StringTableOffset;

		// Two leading int32 fields: NStrings, then TotalBodyBytes (informational)
		int nStrings = ReadInt32LE(fileBytes, pos); pos += 4;
		// TotalBodyBytes is informational only; we don't need it for parsing
		pos += 4;

		var strings = new string[nStrings];

		for (int i = 0; i < nStrings; i++)
		{
			// Each entry: 4-byte length prefix + UTF-8 body bytes (no null terminator)
			int len = ReadInt32LE(fileBytes, pos); pos += 4;
			strings[i] = len > 0 ? Encoding.UTF8.GetString(fileBytes, pos, len) : string.Empty;
			pos += len;
		}

		// Invariant guaranteed by the writer: index 0 must always be the empty string
		if (nStrings > 0 && strings[0] != string.Empty)
			throw new FormatException(
				$"String table invariant violated: index 0 must be the empty string, " +
				$"but got '{strings[0]}'.");

		return strings;
	}

	/// <summary>
	/// Deserialises the protein table from the in-memory file bytes starting at
	/// <c>header.ProteinTableOffset</c>. Uses zero-copy batch casting via
	/// <see cref="MemoryMarshal.Cast{TFrom,TTo}"/>. Returns an empty array when
	/// <c>header.NProteins</c> is zero.
	/// </summary>
	/// <param name="fileBytes">Complete file content.</param>
	/// <param name="header">Deserialized file header providing the offset and count.</param>
	/// <returns>
	/// Array of <see cref="MslProteinRecord"/> structs; may be empty.
	/// </returns>
	private static MslProteinRecord[] ReadProteinTable(byte[] fileBytes, MslFileHeader header)
	{
		int nProteins = header.NProteins;

		if (nProteins == 0)
			return Array.Empty<MslProteinRecord>();

		int byteCount = nProteins * MslFormat.ProteinRecordSize;
		int startPos = (int)header.ProteinTableOffset;

		var proteins = new MslProteinRecord[nProteins];
		ReadOnlySpan<byte> span = fileBytes.AsSpan(startPos, byteCount);
		MemoryMarshal.Cast<byte, MslProteinRecord>(span).CopyTo(proteins);

		return proteins;
	}

	/// <summary>
	/// Deserialises the entire precursor array from the in-memory file bytes starting at
	/// <c>header.PrecursorSectionOffset</c>. Uses zero-copy batch casting for performance
	/// (hot path for large libraries). Returns an empty array for a zero-precursor library.
	/// </summary>
	/// <param name="fileBytes">Complete file content.</param>
	/// <param name="header">Deserialized file header providing the offset and count.</param>
	/// <returns>
	/// Array of <see cref="MslPrecursorRecord"/> structs in file order.
	/// </returns>
	private static MslPrecursorRecord[] ReadPrecursorArray(byte[] fileBytes, MslFileHeader header)
	{
		int nPrecursors = header.NPrecursors;

		if (nPrecursors == 0)
			return Array.Empty<MslPrecursorRecord>();

		int byteCount = nPrecursors * MslFormat.PrecursorRecordSize;
		int startPos = (int)header.PrecursorSectionOffset;

		var precursors = new MslPrecursorRecord[nPrecursors];
		ReadOnlySpan<byte> span = fileBytes.AsSpan(startPos, byteCount);
		MemoryMarshal.Cast<byte, MslPrecursorRecord>(span).CopyTo(precursors);

		return precursors;
	}

	/// <summary>
	/// Reads the extended annotation table section when
	/// <see cref="MslFormat.FileFlagHasExtAnnotations"/> is set in the file header.
	/// Returns an empty array for version-1 files or any file without the flag.
	///
	/// On-disk format:
	/// <code>
	///   int32    NCustomLosses
	///   double[] CustomLossMasses  (NCustomLosses × 8 bytes)
	/// </code>
	///
	/// Index 0 is the reserved sentinel (0.0 = no loss). Valid custom entries start at 1.
	/// </summary>
	private static double[] ReadExtAnnotationTable(byte[] fileBytes, MslFileHeader header)
	{
		// Flag absent or version 1: no extended annotation table in this file
		if ((header.FileFlags & MslFormat.FileFlagHasExtAnnotations) == 0
			|| header.ExtAnnotationTableOffset <= 0)
			return Array.Empty<double>();

		int pos = header.ExtAnnotationTableOffset;
		int count = ReadInt32LE(fileBytes, pos); pos += 4;

		if (count <= 0)
			return Array.Empty<double>();

		var masses = new double[count];
		for (int i = 0; i < count; i++)
		{
			masses[i] = BitConverter.ToDouble(fileBytes, pos);
			pos += 8;
		}

		return masses;
	}

	/// <summary>
	/// Reads the 16-byte compression descriptor at file offset 64 and decompresses the
	/// zstd fragment frame into a new <c>byte[]</c>. Called when
	/// <see cref="MslFormat.FileFlagIsCompressed"/> is set.
	///
	/// After this call, <c>MslPrecursorRecord.FragmentBlockOffset</c> values are treated as
	/// offsets into the returned buffer rather than absolute file positions.
	/// </summary>
	/// <param name="fileBytes">Complete file content.</param>
	/// <param name="header">Deserialized file header (provides <c>FragmentSectionOffset</c>).</param>
	/// <returns>The decompressed fragment section as a new byte array.</returns>
	private static byte[] DecompressFragmentSection(byte[] fileBytes, MslFileHeader header)
	{
		// Compression descriptor is at offset 64 (immediately after the 64-byte header):
		//   int64 CompressedFragmentSize   (offset 64)
		//   int64 UncompressedFragmentSize (offset 72)
		const int DescriptorOffset = MslFormat.HeaderSize; // = 64
		long compressedSize = BitConverter.ToInt64(fileBytes, DescriptorOffset);
		long uncompressedSize = BitConverter.ToInt64(fileBytes, DescriptorOffset + 8);

		// The compressed zstd frame starts at FragmentSectionOffset
		int frameStart = (int)header.FragmentSectionOffset;

		ReadOnlySpan<byte> compressedSpan = fileBytes.AsSpan(frameStart, (int)compressedSize);

		using var decompressor = new Decompressor();
		byte[] decompressed = new byte[uncompressedSize];
		decompressor.Unwrap(compressedSpan, decompressed);

		return decompressed;
	}

	/// <summary>
	/// Deserialises one precursor's fragment block from a byte buffer.
	/// Used in full-load mode (<see cref="Load"/>).
	///
	/// For uncompressed files <paramref name="buffer"/> is the full file content and
	/// <c>precursor.FragmentBlockOffset</c> is an absolute file position.
	/// For compressed files <paramref name="buffer"/> is the decompressed fragment buffer
	/// and the offset is relative to its start.
	///
	/// Returns an empty list when <c>precursor.FragmentCount</c> is zero.
	/// </summary>
	/// <param name="buffer">Byte buffer to read from.</param>
	/// <param name="precursor">
	/// The precursor record whose <c>FragmentBlockOffset</c> and <c>FragmentCount</c>
	/// identify the fragment block.
	/// </param>
	/// <param name="customLossMasses">
	/// Extended annotation table masses. Index 0 is the sentinel (0.0). Pass
	/// <see cref="Array.Empty{T}"/> for files with no custom neutral losses.
	/// </param>
	/// <returns>
	/// List of <see cref="MslFragmentIon"/> objects in m/z ascending order as written.
	/// </returns>
	private static List<MslFragmentIon> ReadFragmentBlockFromBytes(
		byte[] buffer,
		MslPrecursorRecord precursor,
		double[] customLossMasses)
	{
		int fragmentCount = precursor.FragmentCount;

		if (fragmentCount == 0)
			return new List<MslFragmentIon>(0);

		int byteCount = fragmentCount * MslFormat.FragmentRecordSize;
		int startPos = (int)precursor.FragmentBlockOffset;

		var records = new MslFragmentRecord[fragmentCount];
		ReadOnlySpan<byte> span = buffer.AsSpan(startPos, byteCount);
		MemoryMarshal.Cast<byte, MslFragmentRecord>(span).CopyTo(records);

		var ions = new List<MslFragmentIon>(fragmentCount);
		foreach (ref readonly MslFragmentRecord r in records.AsSpan())
			ions.Add(ConvertFragment(in r, customLossMasses));

		return ions;
	}

	/// <summary>
	/// Deserialises one precursor's fragment block from the open on-demand
	/// <see cref="System.IO.FileStream"/> by seeking to the stored offset. Used in
	/// index-only mode; called from <see cref="MslLibrary.LoadFragmentsOnDemand"/>.
	///
	/// Thread-safety: the caller (<see cref="MslLibrary"/>) must hold the library's internal
	/// stream lock around the entire call to prevent concurrent Seek + Read interleaving.
	/// </summary>
	/// <param name="stream">
	/// An open <see cref="System.IO.FileStream"/>; will be seeked to
	/// <c>precursor.FragmentBlockOffset</c> before reading.
	/// </param>
	/// <param name="precursor">The precursor record identifying the fragment block.</param>
	/// <param name="customLossMasses">
	/// Extended annotation table masses passed from the library's cached copy.
	/// Pass <see cref="Array.Empty{T}"/> for files with no custom neutral losses.
	/// </param>
	/// <returns>
	/// List of <see cref="MslFragmentIon"/> objects. Empty when <c>FragmentCount</c> is zero.
	/// </returns>
	internal static List<MslFragmentIon> ReadFragmentBlockFromStream(
		FileStream stream,
		MslPrecursorRecord precursor,
		double[] customLossMasses)
	{
		int fragmentCount = precursor.FragmentCount;

		if (fragmentCount == 0)
			return new List<MslFragmentIon>(0);

		int byteCount = fragmentCount * MslFormat.FragmentRecordSize;

		// Seek to the fragment block; the caller must hold the stream lock
		stream.Seek(precursor.FragmentBlockOffset, SeekOrigin.Begin);

		// Rent a pooled buffer to avoid per-call heap allocation on the hot path
		byte[] buffer = ArrayPool<byte>.Shared.Rent(byteCount);

		try
		{
			// ReadExactly guarantees a complete read; avoids silent short-read bugs
			stream.ReadExactly(buffer, 0, byteCount);

			var records = new MslFragmentRecord[fragmentCount];
			ReadOnlySpan<byte> span = buffer.AsSpan(0, byteCount);
			MemoryMarshal.Cast<byte, MslFragmentRecord>(span).CopyTo(records);

			var ions = new List<MslFragmentIon>(fragmentCount);
			foreach (ref readonly MslFragmentRecord r in records.AsSpan())
				ions.Add(ConvertFragment(in r, customLossMasses));

			return ions;
		}
		finally
		{
			// Return the rented buffer even if an exception occurs during conversion
			ArrayPool<byte>.Shared.Return(buffer);
		}
	}

	// ── Conversion helpers ─────────────────────────────────────────────────────

	/// <summary>
	/// Converts a single <see cref="MslFragmentRecord"/> (20-byte binary representation)
	/// into an <see cref="MslFragmentIon"/> (rich in-memory type).
	///
	/// Type-widening performed:
	///   <c>short</c> FragmentNumber, SecondaryFragmentNumber, ResiduePosition → <c>int</c>
	///   <c>byte</c>  Charge → <c>int</c>
	///   SecondaryProductType == -1 in the record → null in the output (terminal ion sentinel)
	///
	/// Custom neutral-loss decoding: when <c>neutral_loss_code == Custom</c>, the
	/// <c>ResiduePosition</c> field is repurposed as a 1-based index into
	/// <paramref name="customLossMasses"/>. <c>ResiduePosition</c> is set to 0 for such
	/// fragments (documented trade-off; see Prompt 11 design rationale).
	/// </summary>
	/// <param name="r">Raw fragment record, passed by read-only reference to avoid a copy.</param>
	/// <param name="customLossMasses">
	/// Extended annotation table masses. Index 0 is the sentinel (0.0 = no loss).
	/// Pass <see cref="Array.Empty{T}"/> for files with no custom neutral losses.
	/// </param>
	/// <returns>A fully populated <see cref="MslFragmentIon"/>.</returns>
	private static MslFragmentIon ConvertFragment(in MslFragmentRecord r, double[] customLossMasses)
	{
		// Decode the packed flags byte into its four named components
		var (_, _, lossCode, excludeFromQuant) = MslFormat.DecodeFragmentFlags(r.Flags);

		double neutralLoss;
		int residuePosition;

		if (lossCode == MslFormat.NeutralLossCode.Custom)
		{
			// ResiduePosition is repurposed as the 1-based index into the ext annotation table
			int extIdx = r.ResiduePosition;
			neutralLoss = (extIdx > 0 && extIdx < customLossMasses.Length)
				? customLossMasses[extIdx]
				: 0.0;  // defensive fallback; should not occur in a well-formed file
			residuePosition = 0; // not available for custom-loss fragments (documented trade-off)
		}
		else
		{
			neutralLoss = DecodeNeutralLoss(lossCode);
			residuePosition = r.ResiduePosition;
		}

		return new MslFragmentIon
		{
			Mz = r.Mz,
			Intensity = r.Intensity,
			ProductType = (ProductType)r.ProductType,
			// -1 on disk is the "not an internal ion" sentinel; map to null in the rich type
			SecondaryProductType = r.SecondaryProductType == -1
									? (ProductType?)null
									: (ProductType)r.SecondaryProductType,
			// short → int widening required for all three residue-number fields
			FragmentNumber = (int)r.FragmentNumber,
			SecondaryFragmentNumber = (int)r.SecondaryFragmentNumber,
			ResiduePosition = residuePosition,
			// byte → int widening required for charge
			Charge = (int)r.Charge,
			NeutralLoss = neutralLoss,
			ExcludeFromQuant = excludeFromQuant
		};
	}

	/// <summary>
	/// Converts a <see cref="MslPrecursorRecord"/> plus resolved strings, proteins, and
	/// fragment ions into a fully populated <see cref="MslLibraryEntry"/>.
	///
	/// Type-widening and scaling performed:
	///   <c>short</c>  Charge          → <c>int</c>  (explicit cast)
	///   <c>float</c>  PrecursorMz     → <c>double</c> (implicit widening)
	///   <c>float</c>  Irt, IonMobility → <c>double</c> (implicit widening)
	///   <c>short</c>  Nce (stored as NCE × 10) → <c>int</c> NCE via integer division by 10
	/// </summary>
	/// <param name="p">The raw precursor record.</param>
	/// <param name="strings">
	/// Complete string table. Every string-index field in <paramref name="p"/> and in the
	/// protein record indexes into this array.
	/// </param>
	/// <param name="proteins">
	/// Complete protein table. <paramref name="p"/>.ProteinIdx is a zero-based index;
	/// -1 means no protein assignment.
	/// </param>
	/// <param name="fragments">
	/// Pre-loaded fragment ions. Pass a populated list for full-load mode; pass an empty list
	/// for index-only mode (fragments loaded later via <see cref="MslLibrary.LoadFragmentsOnDemand"/>).
	/// </param>
	/// <returns>A fully populated <see cref="MslLibraryEntry"/>.</returns>
	private static MslLibraryEntry ConvertPrecursor(
		MslPrecursorRecord p,
		string[] strings,
		MslProteinRecord[] proteins,
		List<MslFragmentIon> fragments)
	{
		// Decode the single flags byte into three named Booleans
		var (isDecoy, isProteotypic, _) = MslFormat.DecodePrecursorFlags(p.PrecursorFlags);

		// Resolve the optional protein record (ProteinIdx == -1 means no protein)
		bool hasProtein = p.ProteinIdx >= 0 && p.ProteinIdx < proteins.Length;
		MslProteinRecord protein = hasProtein ? proteins[p.ProteinIdx] : default;

		return new MslLibraryEntry
		{
			ModifiedSequence = strings[p.ModifiedSeqStringIdx],
			StrippedSequence = strings[p.StrippedSeqStringIdx],
			// float → double: implicit widening, no precision loss beyond float32 representation
			PrecursorMz = p.PrecursorMz,
			// short → int: explicit cast required (no implicit narrowing in C#)
			Charge = (int)p.Charge,
			Irt = p.Irt,
			IonMobility = p.IonMobility,
			IsDecoy = isDecoy,
			IsProteotypic = isProteotypic,
			QValue = p.QValue,
			ElutionGroupId = p.ElutionGroupId,
			MoleculeType = (MslFormat.MoleculeType)p.MoleculeType,
			DissociationType = (DissociationType)p.DissociationType,
			// Nce stored on disk as NCE × 10 (short); divide by 10 to recover the actual NCE int
			Nce = (int)p.Nce / 10,
			Source = (MslFormat.SourceType)p.SourceType,
			ProteinAccession = hasProtein ? strings[protein.AccessionStringIdx] : string.Empty,
			ProteinName = hasProtein ? strings[protein.NameStringIdx] : string.Empty,
			GeneName = hasProtein ? strings[protein.GeneStringIdx] : string.Empty,
			Fragments = fragments
		};
	}

	// ── Low-level byte helpers ────────────────────────────────────────────────

	/// <summary>
	/// Reads a little-endian int32 from <paramref name="buf"/> at zero-based byte offset
	/// <paramref name="pos"/>. Equivalent to <see cref="System.IO.BinaryReader.ReadInt32"/>
	/// but operates directly on an in-memory byte array without stream overhead.
	/// </summary>
	/// <param name="buf">Source byte array. Must have at least <c>pos + 4</c> bytes.</param>
	/// <param name="pos">Zero-based start offset of the int32 within <paramref name="buf"/>.</param>
	/// <returns>The decoded int32 value in host (little-endian) byte order.</returns>
	private static int ReadInt32LE(byte[] buf, int pos) =>
		buf[pos] | (buf[pos + 1] << 8) | (buf[pos + 2] << 16) | (buf[pos + 3] << 24);
}