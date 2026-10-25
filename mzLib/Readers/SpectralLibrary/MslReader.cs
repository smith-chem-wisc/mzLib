using MassSpectrometry;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using System.Buffers;
using System.Buffers.Binary;
using System.Runtime.InteropServices;
using Microsoft.Win32.SafeHandles;
using System.IO.MemoryMappedFiles;
using System.Text;
using ZstdSharp;

namespace Readers.SpectralLibrary;

/// <summary>
/// Static reader for the .msl (mzLib Spectral Library) binary format.
///
/// Two read modes are supported:
///
///   <see cref="Load"/> — full load: streams the whole file. All precursors,
///   fragment blocks, and strings are deserialized into <see cref="MslLibraryEntry"/>
///   objects. Best for libraries whose entries fit in RAM. No
///   <see cref="System.IO.FileStream"/> is held open after this method returns.
///
///   <see cref="LoadIndexOnly"/> — index-only load: reads only the precursor records and
///   string table. Fragment blocks remain on disk and are fetched lazily from a read-only
///   memory map of the kept-open file. Best for multi-gigabyte libraries
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
	/// Runs <see cref="MslStructs.SizeCheck"/> once when the class is first used.
	/// A Pack-setting mistake causes an immediate throw rather than a silent misread of
	/// every file.
	/// </summary>
	static MslReader()
	{
		MslStructs.SizeCheck();
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
		MslFormat.NeutralLossCode.H3PO4AndH2O => -18.010565 + -97.976895,
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
	/// Reads the .msl file and returns an <see cref="MslLibraryData"/> with all precursor
	/// entries and fragment ions fully loaded. The file is streamed section by section, so its
	/// size is not limited to 2 GB. No file handle is held open after this method returns.
	///
	/// When <see cref="MslFormat.FileFlagIsCompressed"/> is set the fragment section is
	/// decoded from the zstd frame as a stream while fragment blocks are read.
	///
	/// Memory use is proportional to the library's in-memory entries. When those would not fit
	/// in RAM, use <see cref="LoadIndexOnly"/>.
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

		// Stream the file section by section: one byte[] for the whole file would stop at 2 GB
		using var fs = new FileStream(
			filePath, FileMode.Open, FileAccess.Read, FileShare.Read,
			bufferSize: 1 << 20, FileOptions.SequentialScan);

		StreamingValidateAndReadHeader(filePath, fs, out MslFileHeader header, out MslFooter footer);
		return LoadAllEntries(fs, header, footer);
	}

	/// <summary>
	/// Reads only the precursor records and string table into memory; fragment blocks remain
	/// on disk and are fetched lazily from a read-only memory map of the kept-open file.
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
	public static MslLibraryData LoadIndexOnly(string filePath) =>
		LoadIndexOnly(filePath, memoryMapFragments: IsOnLocalFixedDrive(filePath));

	/// <summary>
	/// <see cref="LoadIndexOnly(string)"/> with the fragment read strategy chosen by the caller.
	/// </summary>
	/// <param name="filePath">Path to the .msl file.</param>
	/// <param name="memoryMapFragments">
	/// True to read fragment blocks from a memory map of the file (fastest, lock-free); false to
	/// use positional reads, which surface I/O failures as catchable <see cref="IOException"/>s.
	/// </param>
	internal static MslLibraryData LoadIndexOnly(string filePath, bool memoryMapFragments)
	{
		if (!File.Exists(filePath))
			throw new FileNotFoundException($"MSL file not found: '{filePath}'.", filePath);

		// Open a single FileStream that serves dual purpose:
		//   (a) streaming CRC validation and targeted section reads during open
		//   (b) on-demand fragment reads after the method returns (index-only mode)
		// FileShare.Read allows concurrent readers; prevents writers from modifying the open file.
		var fs = new FileStream(
			filePath, FileMode.Open, FileAccess.Read, FileShare.Read,
			bufferSize: 65536, FileOptions.RandomAccess);

		try
		{
			// Validate magic, version, footer, NPrecursors, and CRC-32 — all streaming.
			// No full-file allocation occurs here; fragment bytes are never read.
			StreamingValidateAndReadHeader(filePath, fs, out MslFileHeader header, out MslFooter footer);

			// Compressed files cannot use index-only mode: fragment offsets are decompressed-
			// buffer-relative and there is no persistent decompressed buffer to seek into.
			// Fall back to full-load transparently.
			bool isCompressed = (header.FileFlags & MslFormat.FileFlagIsCompressed) != 0;
			if (isCompressed)
			{
				System.Diagnostics.Debug.WriteLine(
					$"[MslReader] LoadIndexOnly called on compressed file '{filePath}'; " +
					"falling back to full-load (index-only mode unavailable for compressed files).");

				// Full-load path for compressed files, read from the stream already open
				MslLibraryData full = LoadAllEntries(fs, header, footer);
				fs.Dispose();
				return full;
			}

			// Uncompressed: targeted reads — only the sections the index needs.
			// Fragment bytes are never read; the FileStream stays open for on-demand reads.
			string[] strings = ReadStringTableFromStream(fs, header);
			MslProteinRecord[] proteins = ReadProteinTableFromStream(fs, header);
			MslPrecursorRecord[] precursors = ReadPrecursorArrayFromStream(fs, header);
			double[] customLossMasses = ReadExtAnnotationTableFromStream(fs, header, footer, precursors);

			// Build skeleton entries — fragment lists are intentionally left empty.
			var entries = new List<MslLibraryEntry>(precursors.Length);
			for (int i = 0; i < precursors.Length; i++)
				entries.Add(ConvertPrecursor(precursors[i], strings, proteins, new List<MslFragmentIon>()));

			// Transfer stream ownership to MslLibraryData.
			// fs must NOT be disposed here on the success path.
			return new MslLibraryData(entries, header, precursors, strings, proteins, fs, customLossMasses,
				memoryMapFragments);
		}
		catch
		{
			// Dispose the stream on any failure before the ownership transfer,
			// so the file handle is never leaked on exception.
			fs.Dispose();
			throw;
		}
	}
	/// <summary>
	/// Reads the string table, which runs from <c>StringTableOffset</c> to the precursor section.
	/// Strings are parsed one after another from their length prefixes; the section header's
	/// TotalBodyBytes is informational and is not relied on.
	/// </summary>
	/// <exception cref="FormatException">A length runs past the section, or index 0 is not empty.</exception>
	private static string[] ReadStringTableFromStream(FileStream fs, MslFileHeader header)
	{
		long sectionLength = header.PrecursorSectionOffset - header.StringTableOffset;
		if (sectionLength < 8 || sectionLength > Array.MaxLength)
			throw new FormatException(
				$"Invalid string table bounds ({header.StringTableOffset}..{header.PrecursorSectionOffset}).");

		byte[] table = new byte[sectionLength];
		fs.Seek(header.StringTableOffset, SeekOrigin.Begin);
		fs.ReadExactly(table);

		// Section header: NStrings, then TotalBodyBytes (informational)
		int nStrings = ReadInt32LE(table, 0);
		if (nStrings < 0)
			throw new FormatException($"Invalid string count {nStrings} in the string table.");

		var strings = new string[nStrings];
		long pos = 8;
		for (int i = 0; i < nStrings; i++)
		{
			if (pos + 4 > table.Length)
				throw new FormatException("String table is truncated: it ends inside a length prefix.");
			int len = ReadInt32LE(table, (int)pos); pos += 4;
			if (len < 0 || pos + len > table.Length)
				throw new FormatException($"String {i} (length {len}) runs past the end of the string table.");
			strings[i] = len > 0 ? Encoding.UTF8.GetString(table, (int)pos, len) : string.Empty;
			pos += len;
		}

		if (nStrings > 0 && strings[0] != string.Empty)
			throw new FormatException("String table invariant violated: index 0 must be the empty string.");

		return strings;
	}
	private static MslProteinRecord[] ReadProteinTableFromStream(FileStream fs, MslFileHeader header)
	{
		int nProteins = header.NProteins;
		if (nProteins == 0) return Array.Empty<MslProteinRecord>();

		int byteCount = nProteins * MslFormat.ProteinRecordSize;
		fs.Seek(header.ProteinTableOffset, SeekOrigin.Begin);

		byte[] buf = new byte[byteCount];
		fs.ReadExactly(buf, 0, byteCount);

		var proteins = new MslProteinRecord[nProteins];
		MemoryMarshal.Cast<byte, MslProteinRecord>(buf.AsSpan()).CopyTo(proteins);
		return proteins;
	}

	/// <summary>
	/// Precursor records read per call in <see cref="ReadPrecursorArrayFromStream"/>: 56 MB.
	/// </summary>
	internal const int PrecursorReadChunkRecords = 1 << 20;

	/// <summary>
	/// Reads the precursor array by seeking to <c>header.PrecursorSectionOffset</c>
	/// and reading exactly <c>NPrecursors × PrecursorRecordSize</c> bytes.
	/// <para>
	/// The read is done <paramref name="chunkRecords"/> records at a time: a byte span over the
	/// whole array overflows its int32 length past about 38 M precursors (2^31 / 56).
	/// </para>
	/// </summary>
	internal static MslPrecursorRecord[] ReadPrecursorArrayFromStream(FileStream fs, MslFileHeader header,
		int chunkRecords = PrecursorReadChunkRecords)
	{
		int nPrecursors = header.NPrecursors;
		if (nPrecursors == 0) return Array.Empty<MslPrecursorRecord>();

		fs.Seek(header.PrecursorSectionOffset, SeekOrigin.Begin);

		var precursors = new MslPrecursorRecord[nPrecursors];
		for (int start = 0; start < nPrecursors;)
		{
			int count = Math.Min(chunkRecords, nPrecursors - start);
			fs.ReadExactly(MemoryMarshal.AsBytes(precursors.AsSpan(start, count)));
			start += count;
		}
		return precursors;
	}

	/// <summary>
	/// Reads the extended annotation table (custom neutral-loss masses) when the flag is set.
	/// Returns an empty array for files without custom neutral losses.
	/// <para>
	/// The table is located from the file layout (it follows the fragment section and ends where
	/// the offset table begins), not from <c>MslFileHeader.ExtAnnotationTableOffset</c>: that
	/// field is an int32 and cannot hold the position in a file over 2 GB.
	/// </para>
	/// </summary>
	/// <exception cref="FormatException">The table does not end where the offset table begins.</exception>
	private static double[] ReadExtAnnotationTableFromStream(
		FileStream fs, MslFileHeader header, MslFooter footer, MslPrecursorRecord[] precursors)
	{
		if ((header.FileFlags & MslFormat.FileFlagHasExtAnnotations) == 0)
			return Array.Empty<double>();

		long fragmentSectionLength;
		if ((header.FileFlags & MslFormat.FileFlagIsCompressed) != 0)
		{
			fragmentSectionLength = ReadCompressionDescriptor(fs).CompressedSize;
		}
		else
		{
			fragmentSectionLength = 0;
			foreach (MslPrecursorRecord p in precursors)
				fragmentSectionLength += (long)p.FragmentCount * MslFormat.FragmentRecordSize;
		}

		long tableStart = header.FragmentSectionOffset + fragmentSectionLength;
		fs.Seek(tableStart, SeekOrigin.Begin);

		Span<byte> countBuf = stackalloc byte[4];
		fs.ReadExactly(countBuf);
		int count = BinaryPrimitives.ReadInt32LittleEndian(countBuf);

		if (count < 0 || tableStart + 4 + 8L * count != footer.OffsetTableOffset)
			throw new FormatException(
				$"Extended annotation table at offset {tableStart} (count {count}) does not end where " +
				$"the offset table begins ({footer.OffsetTableOffset}); the file is corrupt.");

		if (count == 0) return Array.Empty<double>();

		var masses = new double[count];
		fs.ReadExactly(MemoryMarshal.AsBytes(masses.AsSpan()));
		return masses;
	}

	/// <summary>
	/// Reads the 16-byte compression descriptor that follows the header in compressed files.
	/// </summary>
	private static (long CompressedSize, long UncompressedSize) ReadCompressionDescriptor(FileStream fs)
	{
		Span<byte> descriptor = stackalloc byte[16];
		fs.Seek(MslFormat.HeaderSize, SeekOrigin.Begin);
		fs.ReadExactly(descriptor);
		return (BinaryPrimitives.ReadInt64LittleEndian(descriptor),
				BinaryPrimitives.ReadInt64LittleEndian(descriptor[8..]));
	}

	// ── Full load ─────────────────────────────────────────────────────────────

	/// <summary>
	/// Reads every section of an already-validated file from <paramref name="fs"/> and builds
	/// the fully-loaded entry list. The file is read in bounded pieces, never as one buffer, so
	/// there is no 2 GB limit (a single <c>byte[]</c> cannot exceed <see cref="Array.MaxLength"/>).
	/// </summary>
	/// <param name="fs">Open stream positioned anywhere; it is seeked as needed.</param>
	/// <param name="header">Header returned by <see cref="StreamingValidateAndReadHeader"/>.</param>
	/// <param name="footer">Footer returned by <see cref="StreamingValidateAndReadHeader"/>.</param>
	/// <returns>A full-load <see cref="MslLibraryData"/> (no open stream).</returns>
	private static MslLibraryData LoadAllEntries(FileStream fs, MslFileHeader header, MslFooter footer)
	{
		string[] strings = ReadStringTableFromStream(fs, header);
		MslProteinRecord[] proteins = ReadProteinTableFromStream(fs, header);
		MslPrecursorRecord[] precursors = ReadPrecursorArrayFromStream(fs, header);
		double[] customLossMasses = ReadExtAnnotationTableFromStream(fs, header, footer, precursors);

		var fragments = new List<MslFragmentIon>[precursors.Length];

		if ((header.FileFlags & MslFormat.FileFlagIsCompressed) != 0)
		{
			// Compressed: FragmentBlockOffset values are offsets into the decompressed section,
			// which is decoded as a stream rather than into one buffer.
			var (compressedSize, uncompressedSize) = ReadCompressionDescriptor(fs);

			fs.Seek(header.FragmentSectionOffset, SeekOrigin.Begin);
			using var frame = new BoundedReadStream(fs, compressedSize);
			using var decompressed = new DecompressionStream(frame, bufferSize: 1 << 20, leaveOpen: true);
			long position = ReadAllFragmentBlocks(
				decompressed, 0, uncompressedSize, precursors, customLossMasses, fragments);

			// Decode the rest of the frame: it must end exactly at the declared uncompressed size.
			// Reaching the end also makes the decoder check that the frame is complete.
			SkipForward(decompressed, uncompressedSize - position);
			if (decompressed.ReadByte() != -1)
				throw new FormatException(
					$"The compressed fragment section decodes to more than its declared {uncompressedSize} bytes.");
		}
		else
		{
			// Uncompressed: FragmentBlockOffset values are absolute file positions
			fs.Seek(header.FragmentSectionOffset, SeekOrigin.Begin);
			ReadAllFragmentBlocks(
				fs, header.FragmentSectionOffset, footer.OffsetTableOffset, precursors, customLossMasses, fragments);
		}

		var entries = new List<MslLibraryEntry>(precursors.Length);
		for (int i = 0; i < precursors.Length; i++)
			entries.Add(ConvertPrecursor(precursors[i], strings, proteins, fragments[i]));

		return new MslLibraryData(entries, header);
	}

	/// <summary>
	/// Reads every precursor's fragment block from <paramref name="source"/> in ascending
	/// offset order, so the source is consumed front to back (the writer lays blocks out in
	/// precursor order, so this is normally a single forward pass with no seeks).
	/// </summary>
	/// <param name="source">
	/// The fragment bytes: the file itself (seekable) or a decompression stream (forward-only).
	/// </param>
	/// <param name="position">
	/// The offset, in the same space as <c>FragmentBlockOffset</c>, that <paramref name="source"/>
	/// is currently positioned at.
	/// </param>
	/// <param name="end">Exclusive upper bound, in the same space, that no block may cross.</param>
	/// <param name="precursors">Precursor records providing each block's offset and count.</param>
	/// <param name="customLossMasses">Extended annotation table masses.</param>
	/// <param name="fragments">Output: the fragment list for each precursor, by precursor index.</param>
	/// <returns>The position after the last block read.</returns>
	/// <exception cref="FormatException">
	/// A block crosses <paramref name="end"/>, or blocks overlap in a forward-only source; a valid
	/// writer produces neither.
	/// </exception>
	private static long ReadAllFragmentBlocks(
		Stream source,
		long position,
		long end,
		MslPrecursorRecord[] precursors,
		double[] customLossMasses,
		List<MslFragmentIon>[] fragments)
	{
		int n = precursors.Length;

		// Visit blocks in offset order; skip the sort when they are already ascending
		int[] order = new int[n];
		for (int i = 0; i < n; i++) order[i] = i;
		bool ascending = true;
		for (int i = 1; i < n && ascending; i++)
			ascending = precursors[i].FragmentBlockOffset >= precursors[i - 1].FragmentBlockOffset;
		if (!ascending)
		{
			long[] offsets = new long[n];
			for (int i = 0; i < n; i++) offsets[i] = precursors[i].FragmentBlockOffset;
			Array.Sort(offsets, order);
		}

		foreach (int idx in order)
		{
			MslPrecursorRecord p = precursors[idx];
			int fragmentCount = p.FragmentCount;

			if (fragmentCount == 0)
			{
				fragments[idx] = new List<MslFragmentIon>(0);
				continue;
			}

			long offset = p.FragmentBlockOffset;
			long blockEnd = offset + (long)fragmentCount * MslFormat.FragmentRecordSize;
			if (offset < 0 || blockEnd > end)
				throw new FormatException(
					$"Fragment block at offset {offset} runs past the end of the fragment section ({end}).");

			if (offset != position)
			{
				if (source.CanSeek)
				{
					source.Seek(offset, SeekOrigin.Begin);
				}
				else if (offset > position)
				{
					// Forward-only source: read and discard the gap
					SkipForward(source, offset - position);
				}
				else
				{
					throw new FormatException(
						$"Fragment blocks overlap at offset {offset}; the fragment section is corrupt.");
				}
			}

			var records = new MslFragmentRecord[fragmentCount];
			source.ReadExactly(MemoryMarshal.AsBytes(records.AsSpan()));
			position = blockEnd;

			var ions = new List<MslFragmentIon>(fragmentCount);
			foreach (ref readonly MslFragmentRecord r in records.AsSpan())
				ions.Add(ConvertFragment(in r, customLossMasses));
			fragments[idx] = ions;
		}

		return position;
	}

	/// <summary>Reads and discards <paramref name="count"/> bytes from a forward-only stream.</summary>
	private static void SkipForward(Stream source, long count)
	{
		if (count <= 0) return;
		byte[] buffer = ArrayPool<byte>.Shared.Rent(81920);
		try
		{
			while (count > 0)
			{
				int chunk = (int)Math.Min(count, buffer.Length);
				source.ReadExactly(buffer, 0, chunk);
				count -= chunk;
			}
		}
		finally
		{
			ArrayPool<byte>.Shared.Return(buffer);
		}
	}

	/// <summary>
	/// Read-only view of the next <c>length</c> bytes of an underlying stream. Bounds the zstd
	/// decoder to the compressed frame so it cannot read into the sections that follow it.
	/// </summary>
	private sealed class BoundedReadStream : Stream
	{
		private readonly Stream _inner;
		private long _remaining;

		public BoundedReadStream(Stream inner, long length)
		{
			_inner = inner;
			_remaining = length;
		}

		public override int Read(byte[] buffer, int offset, int count) =>
			Read(buffer.AsSpan(offset, count));

		public override int Read(Span<byte> buffer)
		{
			if (_remaining <= 0) return 0;
			if (buffer.Length > _remaining) buffer = buffer[..(int)_remaining];
			int read = _inner.Read(buffer);
			_remaining -= read;
			return read;
		}

		public override bool CanRead => true;
		public override bool CanSeek => false;
		public override bool CanWrite => false;
		public override long Length => throw new NotSupportedException();
		public override long Position
		{
			get => throw new NotSupportedException();
			set => throw new NotSupportedException();
		}
		public override void Flush() { }
		public override long Seek(long offset, SeekOrigin origin) => throw new NotSupportedException();
		public override void SetLength(long value) => throw new NotSupportedException();
		public override void Write(byte[] buffer, int offset, int count) => throw new NotSupportedException();
	}

	/// <summary>
	/// Deserialises one precursor's fragment block from a memory-mapped view of the file.
	/// Used in index-only mode; called from <c>MslLibraryData.LoadFragmentsOnDemand</c>.
	///
	/// Thread-safe without a lock: the view is read-only and has no shared position.
	/// </summary>
	/// <param name="view">Read-only view over the whole file, starting at offset 0.</param>
	/// <param name="fileLength">File length; the view's capacity is page-rounded, so it bounds the read.</param>
	/// <param name="precursor">The precursor record identifying the fragment block.</param>
	/// <param name="customLossMasses">
	/// Extended annotation table masses passed from the library's cached copy.
	/// Pass <see cref="Array.Empty{T}"/> for files with no custom neutral losses.
	/// </param>
	/// <returns>
	/// List of <see cref="MslFragmentIon"/> objects. Empty when <c>FragmentCount</c> is zero.
	/// </returns>
	/// <exception cref="EndOfStreamException">The block extends past the end of the file.</exception>
	internal static List<MslFragmentIon> ReadFragmentBlockAt(
		MemoryMappedViewAccessor view,
		long fileLength,
		MslPrecursorRecord precursor,
		double[] customLossMasses)
	{
		int fragmentCount = precursor.FragmentCount;

		if (fragmentCount == 0)
			return new List<MslFragmentIon>(0);

		long end = precursor.FragmentBlockOffset + (long)fragmentCount * MslFormat.FragmentRecordSize;
		if (precursor.FragmentBlockOffset < 0 || end > fileLength)
			throw new EndOfStreamException(
				$"Fragment block at offset {precursor.FragmentBlockOffset} extends past the end of the file.");

		var records = new MslFragmentRecord[fragmentCount];
		view.ReadArray(precursor.FragmentBlockOffset, records, 0, fragmentCount);
		return ToIons(records, customLossMasses);
	}

	/// <summary>
	/// Deserialises one precursor's fragment block with positional reads from
	/// <paramref name="handle"/>. Used in index-only mode for files that are not memory-mapped
	/// (see <see cref="IsOnLocalFixedDrive"/>). Thread-safe without a lock: positional reads do
	/// not use a shared file position.
	/// </summary>
	/// <exception cref="EndOfStreamException">The block extends past the end of the file.</exception>
	internal static List<MslFragmentIon> ReadFragmentBlockAt(
		SafeFileHandle handle,
		MslPrecursorRecord precursor,
		double[] customLossMasses)
	{
		int fragmentCount = precursor.FragmentCount;

		if (fragmentCount == 0)
			return new List<MslFragmentIon>(0);

		var records = new MslFragmentRecord[fragmentCount];
		Span<byte> target = MemoryMarshal.AsBytes(records.AsSpan());

		// RandomAccess.Read may return fewer bytes than asked; loop until the block is complete
		int filled = 0;
		while (filled < target.Length)
		{
			int read = RandomAccess.Read(handle, target[filled..], precursor.FragmentBlockOffset + filled);
			if (read == 0)
				throw new EndOfStreamException(
					$"Fragment block at offset {precursor.FragmentBlockOffset} extends past the end of the file.");
			filled += read;
		}

		return ToIons(records, customLossMasses);
	}

	private static List<MslFragmentIon> ToIons(MslFragmentRecord[] records, double[] customLossMasses)
	{
		var ions = new List<MslFragmentIon>(records.Length);
		foreach (ref readonly MslFragmentRecord r in records.AsSpan())
			ions.Add(ConvertFragment(in r, customLossMasses));
		return ions;
	}

	/// <summary>
	/// True when <paramref name="filePath"/> is on a local fixed (or RAM) drive. Index-only mode
	/// memory-maps only such files: if a network or removable drive fails while a mapped page is
	/// read, the fault ends the process instead of raising a catchable <see cref="IOException"/>.
	/// Any doubt (UNC path, unknown drive, lookup failure) answers false.
	/// </summary>
	internal static bool IsOnLocalFixedDrive(string filePath)
	{
		try
		{
			string fullPath = Path.GetFullPath(filePath);
			if (fullPath.StartsWith(@"\\", StringComparison.Ordinal))
				return false;  // UNC path: a network share

			// The drive (mount) holding the file is the one with the longest matching root
			StringComparison comparison = OperatingSystem.IsWindows()
				? StringComparison.OrdinalIgnoreCase
				: StringComparison.Ordinal;
			DriveInfo? drive = null;
			foreach (DriveInfo candidate in DriveInfo.GetDrives())
			{
				string root = candidate.RootDirectory.FullName;
				if (fullPath.StartsWith(root, comparison)
					&& (drive is null || root.Length > drive.RootDirectory.FullName.Length))
					drive = candidate;
			}

			return drive?.DriveType is DriveType.Fixed or DriveType.Ram;
		}
		catch (Exception)
		{
			return false;
		}
	}

	// ── Conversion helpers ─────────────────────────────────────────────────────

	/// <summary>
	/// Converts a single <see cref="MslFragmentRecord"/> (20-byte binary representation)
	/// into an <see cref="MslFragmentIon"/> (rich in-memory type).
	///
	/// Type-widening performed:
	///   <c>short</c> FragmentNumber, SecondaryFragmentNumber, ResiduePosition → <c>int</c>
	///   <c>byte</c>  ChargeState → <c>int</c>
	///   SecondaryProductType == -1 in the record → null in the output (terminal ion sentinel)
	///
	/// Custom neutral-loss decoding: when <c>neutral_loss_code == Custom</c>, the
	/// <c>ResiduePosition</c> field is repurposed as a 1-based index into
	/// <paramref name="customLossMasses"/>. <c>ResiduePosition</c> is set to 0 for such
	/// fragments (documented trade-off).
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
	///   <c>short</c>  ChargeState          → <c>int</c>  (explicit cast)
	///   <c>float</c>  PrecursorMz     → <c>double</c> (implicit widening)
	///   <c>float</c>  RetentionTime, IonMobility → <c>double</c> (implicit widening)
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
	/// for index-only mode (fragments loaded later via <c>MslLibrary.LoadFragmentsOnDemand</c>).
	/// </param>
	/// <returns>A fully populated <see cref="MslLibraryEntry"/>.</returns>
	private static MslLibraryEntry ConvertPrecursor(
		MslPrecursorRecord p,
		string[] strings,
		MslProteinRecord[] proteins,
		List<MslFragmentIon> fragments)
	{
		// Decode the single flags byte into its named Booleans
		var (isDecoy, isProteotypic, rtCalibrated) = MslFormat.DecodePrecursorFlags(p.PrecursorFlags);
		bool isEntrapment = MslFormat.DecodeIsEntrapment(p.PrecursorFlags);

		// Resolve the optional protein record (ProteinIdx == -1 means no protein)
		bool hasProtein = p.ProteinIdx >= 0 && p.ProteinIdx < proteins.Length;
		MslProteinRecord protein = hasProtein ? proteins[p.ProteinIdx] : default;

		return new MslLibraryEntry
		{
			FullSequence = strings[p.ModifiedSeqStringIdx],
			BaseSequence = strings[p.StrippedSeqStringIdx],
			// float → double: implicit widening, no precision loss beyond float32 representation
			PrecursorMz = p.PrecursorMz,
			// short → int: explicit cast required (no implicit narrowing in C#)
			ChargeState = (int)p.Charge,
			RetentionTime = p.Irt,
			RtIsCalibrated = rtCalibrated,
			IonMobility = p.IonMobility,
			IsDecoy = isDecoy,
			IsEntrapment = isEntrapment,
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
			MatchedFragmentIons = fragments
		};
	}
	private static void StreamingValidateAndReadHeader(
	string filePath,
	FileStream fs,
	out MslFileHeader header,
	out MslFooter footer,
	int chunkSize = 1 << 20)
	{
		long fileLength = fs.Length;
		int minimumSize = MslFormat.HeaderSize + MslFormat.FooterSize;

		if (fileLength < minimumSize)
			throw new FormatException(
				$"File '{filePath}' is too short ({fileLength} bytes) " +
				$"to be a valid .msl file (minimum {minimumSize} bytes).");

		// ── Step 1: read header (64 bytes) ───────────────────────────────
		byte[] headerBytes = new byte[MslFormat.HeaderSize];
		fs.Seek(0, SeekOrigin.Begin);
		fs.ReadExactly(headerBytes, 0, MslFormat.HeaderSize);

		if (!MslFormat.MagicMatches(headerBytes.AsSpan(0, 4)))
			throw new FormatException(
				$"Magic mismatch in '{filePath}': not an MSL file " +
				$"(got 0x{headerBytes[0]:X2}{headerBytes[1]:X2}" +
				$"{headerBytes[2]:X2}{headerBytes[3]:X2}).");

		header = MemoryMarshal.Read<MslFileHeader>(headerBytes.AsSpan());

		if (header.FormatVersion < 1 || header.FormatVersion > MslFormat.CurrentVersion)
			throw new FormatException(
				$"Unsupported version: {header.FormatVersion} in '{filePath}'. " +
				$"This reader supports versions 1–{MslFormat.CurrentVersion}.");

		// ── Step 2: read footer (last 20 bytes) ──────────────────────────
		byte[] footerBytes = new byte[MslFormat.FooterSize];
		fs.Seek(-MslFormat.FooterSize, SeekOrigin.End);
		fs.ReadExactly(footerBytes, 0, MslFormat.FooterSize);

		footer = MemoryMarshal.Read<MslFooter>(footerBytes.AsSpan());

		// Trailing magic check (last 4 bytes of file)
		if (!MslFormat.MagicMatches(footerBytes.AsSpan(MslFormat.FooterSize - 4, 4)))
			throw new FormatException(
				$"Trailing magic mismatch in '{filePath}': file is truncated or not a valid .msl file.");

		// NPrecursors cross-check
		if (footer.NPrecursors != header.NPrecursors)
			throw new FormatException(
				$"NPrecursors mismatch in '{filePath}': " +
				$"header says {header.NPrecursors}, footer says {footer.NPrecursors}. " +
				"File may be truncated or corrupt.");

		// ── Step 3: streaming CRC-32 ──────────────────────────────────────
		long crcEndOffset = footer.OffsetTableOffset;

		if (crcEndOffset < 0 || crcEndOffset > fileLength)
			throw new FormatException(
				$"Invalid OffsetTableOffset ({crcEndOffset}) in footer of '{filePath}'.");

		uint crc = MslCrc32.Initial;
		byte[] chunk = ArrayPool<byte>.Shared.Rent(chunkSize);

		try
		{
			fs.Seek(0, SeekOrigin.Begin);
			long remaining = crcEndOffset;

			while (remaining > 0)
			{
				int toRead = (int)Math.Min(remaining, chunkSize);
				int read = fs.Read(chunk, 0, toRead);
				if (read == 0) break;

				crc = MslCrc32.Update(crc, chunk.AsSpan(0, read));

				remaining -= read;
			}
		}
		finally
		{
			ArrayPool<byte>.Shared.Return(chunk);
		}

		uint computedCrc = MslCrc32.Finish(crc);

		if (computedCrc != footer.DataCrc32)
			throw new InvalidDataException(
				$"CRC32 mismatch: file may be corrupted. " +
				$"Stored: 0x{footer.DataCrc32:X8}, Computed: 0x{computedCrc:X8}.");
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