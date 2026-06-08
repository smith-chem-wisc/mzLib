using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Targeted coverage for every untested exception path in <see cref="MslReader"/>.
///
/// Each test names the method and the check number from ValidateFileBytes / 
/// StreamingValidateAndReadHeader so the mapping back to source is unambiguous.
///
/// Technique: write a valid .msl file, read it to bytes, corrupt the specific field
/// that guards the throw, write the bytes back (without re-computing CRC where the
/// corrupt field is itself what we are testing), then assert the expected exception.
/// Where the corrupted field is inside the CRC-covered region, the CRC is re-computed
/// so we reach the intended check rather than failing earlier at CRC.
///
/// Helper used throughout: MslWriter.ComputeCrc32OfArray (internal, same assembly).
/// </summary>
[TestFixture]
public sealed class TestMslReaderCoverage
{
	// ── Fixture setup ─────────────────────────────────────────────────────────

	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslReaderCoverageTests");

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDirectory);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
			Directory.Delete(OutputDirectory, recursive: true);
	}

	private static string TempPath(string name) =>
		Path.Combine(OutputDirectory, name + ".msl");

	// ── Canonical one-entry library used by most tests ────────────────────────

	private static MslLibraryEntry MakeEntry() => new MslLibraryEntry
	{
		FullSequence = "PEPTIDE",
		BaseSequence = "PEPTIDE",
		PrecursorMz = 449.74,
		ChargeState = 2,
		RetentionTime = 35.4,
		MoleculeType = MslFormat.MoleculeType.Peptide,
		DissociationType = DissociationType.HCD,
		MatchedFragmentIons = new List<MslFragmentIon>
		{
			new MslFragmentIon
			{
				Mz = 175.119f, Intensity = 1.0f,
				ProductType = ProductType.y, FragmentNumber = 1,
				ResiduePosition = 6, Charge = 1
			}
		}
	};

	/// <summary>
	/// Write a valid library and return its bytes. All patch helpers below call this.
	/// </summary>
	private static byte[] WriteAndReadBytes(string name)
	{
		string path = TempPath(name);
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() });
		return File.ReadAllBytes(path);
	}

	/// <summary>
	/// Re-computes the CRC over bytes[0..offsetTableOffset) and patches it into the
	/// footer so that a byte change inside the CRC-covered region does not cause the
	/// test to fail at the CRC check before reaching the intended check.
	/// </summary>
	private static void RecomputeCrc(byte[] bytes)
	{
		int footerStart = bytes.Length - MslFormat.FooterSize;
		long offsetTableOffset = BitConverter.ToInt64(bytes, footerStart);
		uint newCrc = MslWriter.ComputeCrc32OfArray(bytes, (int)offsetTableOffset);
		byte[] crcBytes = BitConverter.GetBytes(newCrc);
		// DataCrc32 is at footer offset +12
		Array.Copy(crcBytes, 0, bytes, footerStart + 12, 4);
	}

	// ═════════════════════════════════════════════════════════════════════════
	// ReadHeaderOnly — 4 throw sites
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// ReadHeaderOnly: FileNotFoundException when file does not exist.
	/// </summary>
	[Test]
	public void ReadHeaderOnly_MissingFile_ThrowsFileNotFoundException()
	{
		Assert.Throws<FileNotFoundException>(() =>
			MslReader.ReadHeaderOnly(Path.Combine(OutputDirectory, "no_such_file.msl")));
	}

	/// <summary>
	/// ReadHeaderOnly: FormatException when the first 4 bytes are not "MZLB".
	/// </summary>
	[Test]
	public void ReadHeaderOnly_BadMagic_ThrowsFormatException()
	{
		string path = TempPath("rheader_bad_magic");
		byte[] bytes = WriteAndReadBytes("rheader_bad_magic_src");
		// Corrupt the first 4 magic bytes
		bytes[0] = 0xFF; bytes[1] = 0xFF; bytes[2] = 0xFF; bytes[3] = 0xFF;
		File.WriteAllBytes(path, bytes);

		var ex = Assert.Throws<FormatException>(() => MslReader.ReadHeaderOnly(path));
		Assert.That(ex!.Message, Does.Contain("Magic mismatch"));
	}

	/// <summary>
	/// ReadHeaderOnly: FormatException "too short to contain a complete 64-byte MSL header"
	/// when the file has valid magic but fewer than 64 bytes total.
	/// We write exactly 4 valid magic bytes followed by padding to produce a 10-byte file.
	/// </summary>
	[Test]
	public void ReadHeaderOnly_FileShorterThanHeader_ThrowsFormatException()
	{
		string path = TempPath("rheader_short");
		// Write valid magic bytes then stop — file is only 10 bytes
		byte[] tiny = new byte[10];
		// MZLB magic in little-endian struct representation used by the writer
		tiny[0] = 0x4D; tiny[1] = 0x5A; tiny[2] = 0x4C; tiny[3] = 0x42;
		File.WriteAllBytes(path, tiny);

		var ex = Assert.Throws<FormatException>(() => MslReader.ReadHeaderOnly(path));
		Assert.That(ex!.Message, Does.Contain("too short"));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Load — FileNotFoundException
	// ═════════════════════════════════════════════════════════════════════════

	[Test]
	public void Load_MissingFile_ThrowsFileNotFoundException()
	{
		Assert.Throws<FileNotFoundException>(() =>
			MslReader.Load(Path.Combine(OutputDirectory, "no_such_file.msl")));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// LoadIndexOnly — FileNotFoundException
	// ═════════════════════════════════════════════════════════════════════════

	[Test]
	public void LoadIndexOnly_MissingFile_ThrowsFileNotFoundException()
	{
		Assert.Throws<FileNotFoundException>(() =>
			MslReader.LoadIndexOnly(Path.Combine(OutputDirectory, "no_such_file.msl")));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// ValidateFileBytes (reached via Load) — 6 throw sites
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Check 1 — FormatException when file is smaller than HeaderSize + FooterSize.
	/// We write a file of only 4 bytes (less than the 84-byte minimum).
	/// </summary>
	[Test]
	public void Load_FileTooShort_ThrowsFormatException()
	{
		string path = TempPath("load_too_short");
		File.WriteAllBytes(path, new byte[4]);

		var ex = Assert.Throws<FormatException>(() => MslReader.Load(path));
		Assert.That(ex!.Message, Does.Contain("too short"));
	}

	/// <summary>
	/// Check 2 — FormatException when the leading magic is wrong.
	/// CRC is re-computed so we reach the magic check cleanly.
	/// </summary>
	[Test]
	public void Load_BadLeadingMagic_ThrowsFormatException()
	{
		string path = TempPath("load_bad_magic");
		byte[] bytes = WriteAndReadBytes("load_bad_magic_src");
		bytes[0] = 0xDE; bytes[1] = 0xAD; bytes[2] = 0xBE; bytes[3] = 0xEF;
		// Magic is checked before CRC — no need to recompute CRC here
		File.WriteAllBytes(path, bytes);

		var ex = Assert.Throws<FormatException>(() => MslReader.Load(path));
		Assert.That(ex!.Message, Does.Contain("Magic mismatch"));
	}

	/// <summary>
	/// Check 3 — FormatException for unsupported format version.
	/// Patch FormatVersion to 99, recompute CRC so validation reaches the version check.
	/// </summary>
	[Test]
	public void Load_UnsupportedVersion_ThrowsFormatException()
	{
		string path = TempPath("load_bad_version");
		byte[] bytes = WriteAndReadBytes("load_bad_version_src");
		// FormatVersion is a little-endian int32 at header offset 4
		byte[] v99 = BitConverter.GetBytes(99);
		Array.Copy(v99, 0, bytes, 4, 4);
		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		var ex = Assert.Throws<FormatException>(() => MslReader.Load(path));
		Assert.That(ex!.Message, Does.Contain("Unsupported version"));
		Assert.That(ex.Message, Does.Contain("99"));
	}

	/// <summary>
	/// Check 3 — FormatException for version 0 (below the minimum of 1).
	/// </summary>
	[Test]
	public void Load_VersionZero_ThrowsFormatException()
	{
		string path = TempPath("load_version_zero");
		byte[] bytes = WriteAndReadBytes("load_version_zero_src");
		byte[] v0 = BitConverter.GetBytes(0);
		Array.Copy(v0, 0, bytes, 4, 4);
		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		var ex = Assert.Throws<FormatException>(() => MslReader.Load(path));
		Assert.That(ex!.Message, Does.Contain("Unsupported version"));
	}

	/// <summary>
	/// Check 4 — FormatException when the trailing magic (last 4 bytes) is wrong.
	/// The trailing magic is NOT inside the CRC-covered region, so no CRC patch needed.
	/// </summary>
	[Test]
	public void Load_BadTrailingMagic_ThrowsFormatException()
	{
		string path = TempPath("load_bad_trailing_magic");
		byte[] bytes = WriteAndReadBytes("load_bad_trailing_magic_src");
		// Corrupt the last 4 bytes (trailing magic)
		int tail = bytes.Length - 4;
		bytes[tail] = 0x00; bytes[tail + 1] = 0x00;
		bytes[tail + 2] = 0x00; bytes[tail + 3] = 0x00;
		File.WriteAllBytes(path, bytes);

		var ex = Assert.Throws<FormatException>(() => MslReader.Load(path));
		Assert.That(ex!.Message, Does.Contain("Trailing magic mismatch"));
	}

	/// <summary>
	/// Check 5 — FormatException when header NPrecursors != footer NPrecursors.
	/// Patch header NPrecursors (offset 12) to a wrong value and recompute CRC.
	/// </summary>
	[Test]
	public void Load_NPrecursorsMismatch_ThrowsFormatException()
	{
		string path = TempPath("load_nprecursors_mismatch");
		byte[] bytes = WriteAndReadBytes("load_nprecursors_mismatch_src");
		// NPrecursors is int32 at header offset 12 — set it to 99
		byte[] n99 = BitConverter.GetBytes(99);
		Array.Copy(n99, 0, bytes, 12, 4);
		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		var ex = Assert.Throws<FormatException>(() => MslReader.Load(path));
		Assert.That(ex!.Message, Does.Contain("NPrecursors mismatch"));
	}

	/// <summary>
	/// Check 6 — InvalidDataException on CRC mismatch.
	/// Flip one byte in the middle of the file (inside the CRC-covered region)
	/// without recomputing the CRC.
	/// </summary>
	[Test]
	public void Load_CrcMismatch_ThrowsInvalidDataException()
	{
		string path = TempPath("load_crc_mismatch");
		byte[] bytes = WriteAndReadBytes("load_crc_mismatch_src");
		// Corrupt a byte well inside the fragment section (not header/footer)
		bytes[MslFormat.HeaderSize + 4] ^= 0xFF;
		File.WriteAllBytes(path, bytes);

		Assert.Throws<InvalidDataException>(() => MslReader.Load(path));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// StreamingValidateAndReadHeader (reached via LoadIndexOnly) — same checks
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Streaming Check 1 — FormatException when file is too short.
	/// </summary>
	[Test]
	public void LoadIndexOnly_FileTooShort_ThrowsFormatException()
	{
		string path = TempPath("indexonly_too_short");
		File.WriteAllBytes(path, new byte[4]);

		var ex = Assert.Throws<FormatException>(() => MslReader.LoadIndexOnly(path));
		Assert.That(ex!.Message, Does.Contain("too short"));
	}

	/// <summary>
	/// Streaming Check 2 — FormatException when the leading magic is wrong.
	/// </summary>
	[Test]
	public void LoadIndexOnly_BadLeadingMagic_ThrowsFormatException()
	{
		string path = TempPath("indexonly_bad_magic");
		byte[] bytes = WriteAndReadBytes("indexonly_bad_magic_src");
		bytes[0] = 0xDE; bytes[1] = 0xAD; bytes[2] = 0xBE; bytes[3] = 0xEF;
		File.WriteAllBytes(path, bytes);

		var ex = Assert.Throws<FormatException>(() => MslReader.LoadIndexOnly(path));
		Assert.That(ex!.Message, Does.Contain("Magic mismatch"));
	}

	/// <summary>
	/// Streaming Check 3 — FormatException for unsupported version.
	/// </summary>
	[Test]
	public void LoadIndexOnly_UnsupportedVersion_ThrowsFormatException()
	{
		string path = TempPath("indexonly_bad_version");
		byte[] bytes = WriteAndReadBytes("indexonly_bad_version_src");
		byte[] v99 = BitConverter.GetBytes(99);
		Array.Copy(v99, 0, bytes, 4, 4);
		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		var ex = Assert.Throws<FormatException>(() => MslReader.LoadIndexOnly(path));
		Assert.That(ex!.Message, Does.Contain("Unsupported version"));
	}

	/// <summary>
	/// Streaming Check 4 — FormatException when trailing magic is wrong.
	/// Also verifies the FileStream is not leaked (file can be deleted immediately after).
	/// </summary>
	[Test]
	public void LoadIndexOnly_BadTrailingMagic_ThrowsAndDoesNotLeakStream()
	{
		string path = TempPath("indexonly_bad_trailing");
		byte[] bytes = WriteAndReadBytes("indexonly_bad_trailing_src");
		int tail = bytes.Length - 4;
		bytes[tail] = 0x00; bytes[tail + 1] = 0x00;
		bytes[tail + 2] = 0x00; bytes[tail + 3] = 0x00;
		File.WriteAllBytes(path, bytes);

		Assert.Throws<FormatException>(() => MslReader.LoadIndexOnly(path));

		// If the stream were leaked, this Delete would throw on Windows
		Assert.DoesNotThrow(() => File.Delete(path),
			"FileStream must be disposed in the catch block — file must be deletable after throw");
	}

	/// <summary>
	/// Streaming Check 5 — FormatException when header NPrecursors != footer NPrecursors.
	/// </summary>
	[Test]
	public void LoadIndexOnly_NPrecursorsMismatch_ThrowsFormatException()
	{
		string path = TempPath("indexonly_nprecursors");
		byte[] bytes = WriteAndReadBytes("indexonly_nprecursors_src");
		byte[] n99 = BitConverter.GetBytes(99);
		Array.Copy(n99, 0, bytes, 12, 4);
		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		var ex = Assert.Throws<FormatException>(() => MslReader.LoadIndexOnly(path));
		Assert.That(ex!.Message, Does.Contain("NPrecursors mismatch"));
	}

	/// <summary>
	/// Streaming CRC — InvalidDataException on checksum mismatch via LoadIndexOnly.
	/// Also verifies the stream is disposed on failure (file deletable after throw).
	/// </summary>
	[Test]
	public void LoadIndexOnly_CrcMismatch_ThrowsAndDoesNotLeakStream()
	{
		string path = TempPath("indexonly_crc");
		byte[] bytes = WriteAndReadBytes("indexonly_crc_src");
		bytes[MslFormat.HeaderSize + 4] ^= 0xFF;
		File.WriteAllBytes(path, bytes);

		Assert.Throws<InvalidDataException>(() => MslReader.LoadIndexOnly(path));

		Assert.DoesNotThrow(() => File.Delete(path),
			"FileStream must be disposed in the catch block — file must be deletable after throw");
	}

	// ═════════════════════════════════════════════════════════════════════════
	// ReadStringTable / ReadStringTableFromStream — string table invariant
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// FormatException when string table index 0 is not the empty string.
	///
	/// Strategy: patch NStrings to 1 so the reader loop exits after reading only
	/// entry 0, then set entry 0's length prefix to 1. With NStrings=1 the reader
	/// reads exactly one byte of body (whatever happens to follow in the file — a
	/// valid byte, no overrun) and decodes a non-empty strings[0], which trips the
	/// invariant check. NStrings is in the CRC-covered region so we recompute.
	/// </summary>
	[Test]
	public void Load_StringTableIndex0NotEmpty_ThrowsFormatException()
	{
		string path = TempPath("stringtable_invariant");
		byte[] bytes = WriteAndReadBytes("stringtable_invariant_src");

		MslFileHeader header = MemoryMarshal.Read<MslFileHeader>(
			bytes.AsSpan(0, MslFormat.HeaderSize));

		int tableBase = (int)header.StringTableOffset;

		// Patch NStrings to 1 — loop reads only entry 0, no overrun from later entries
		bytes[tableBase] = 1;
		bytes[tableBase + 1] = 0;
		bytes[tableBase + 2] = 0;
		bytes[tableBase + 3] = 0;

		// Patch entry 0 length prefix (at tableBase+8, after NStrings:4 + TotalBodyBytes:4) to 1
		// The one body byte read will be whatever follows — a valid file byte, no overrun
		int entry0LenOffset = tableBase + 8;
		bytes[entry0LenOffset] = 1;
		bytes[entry0LenOffset + 1] = 0;
		bytes[entry0LenOffset + 2] = 0;
		bytes[entry0LenOffset + 3] = 0;

		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		var ex = Assert.Throws<FormatException>(() => MslReader.Load(path));
		Assert.That(ex!.Message, Does.Contain("String table invariant"));
	}

	/// <summary>
	/// Same invariant check reached via LoadIndexOnly (ReadStringTableFromStream path).
	///
	/// ReadStringTableFromStream allocates tableData = new byte[nStrings*4 + totalBodyBytes].
	/// We must also patch TotalBodyBytes (+1) so the allocated buffer is large enough
	/// to hold the 1-byte body we are injecting for entry 0.
	/// </summary>
	[Test]
	public void LoadIndexOnly_StringTableIndex0NotEmpty_ThrowsFormatException()
	{
		string path = TempPath("stringtable_stream_invariant");
		byte[] bytes = WriteAndReadBytes("stringtable_stream_invariant_src");

		MslFileHeader header = MemoryMarshal.Read<MslFileHeader>(
			bytes.AsSpan(0, MslFormat.HeaderSize));

		int tableBase = (int)header.StringTableOffset;

		// Patch NStrings to 1 so the stream reader only processes entry 0
		bytes[tableBase] = 1;
		bytes[tableBase + 1] = 0;
		bytes[tableBase + 2] = 0;
		bytes[tableBase + 3] = 0;

		// Patch TotalBodyBytes (tableBase+4) from 0 to 1 so tableData allocation includes
		// the 1 body byte we need for entry 0
		int totalBodyOffset = tableBase + 4;
		int originalTotal = BitConverter.ToInt32(bytes, totalBodyOffset);
		byte[] newTotal = BitConverter.GetBytes(originalTotal + 1);
		Array.Copy(newTotal, 0, bytes, totalBodyOffset, 4);

		// Patch entry 0 length prefix to 1
		int entry0LenOffset = tableBase + 8;
		bytes[entry0LenOffset] = 1;
		bytes[entry0LenOffset + 1] = 0;
		bytes[entry0LenOffset + 2] = 0;
		bytes[entry0LenOffset + 3] = 0;

		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		Assert.Throws<FormatException>(() => MslReader.LoadIndexOnly(path));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// ValidateFileBytes — invalid OffsetTableOffset check
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// FormatException when OffsetTableOffset in the footer is negative or beyond EOF.
	///
	/// OffsetTableOffset is the first 8 bytes of the footer. We set it to a value
	/// larger than the file so the bounds check fires. The trailing magic must remain
	/// intact (last 4 bytes), and the NPrecursors cross-check must still pass,
	/// so we only change OffsetTableOffset and leave everything else alone.
	/// The CRC check uses OffsetTableOffset itself as the end boundary — a value
	/// beyond EOF triggers the bounds check before CRC is attempted.
	/// </summary>
	[Test]
	public void Load_InvalidOffsetTableOffset_ThrowsFormatException()
	{
		string path = TempPath("bad_offset_table");
		byte[] bytes = WriteAndReadBytes("bad_offset_table_src");

		int footerStart = bytes.Length - MslFormat.FooterSize;
		// Write a huge OffsetTableOffset that is beyond EOF
		byte[] bigOffset = BitConverter.GetBytes((long)(bytes.Length + 99999));
		Array.Copy(bigOffset, 0, bytes, footerStart, 8);
		File.WriteAllBytes(path, bytes);

		var ex = Assert.Throws<FormatException>(() => MslReader.Load(path));
		Assert.That(ex!.Message, Does.Contain("OffsetTableOffset"));
	}
}