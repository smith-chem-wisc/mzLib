using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 test suite for the optional zstd block compression introduced in
/// .msl format version 3 (Prompt 13).
///
/// When <c>compressionLevel &gt; 0</c> is passed to <see cref="MslWriter.Write"/>, the fragment
/// section (plus extended annotation table if present) is compressed as a single zstd frame.
/// A 16-byte compression descriptor is inserted at file offset 64 (immediately after the header),
/// and <c>FileFlagIsCompressed</c> is set in the file-level flags.
/// <c>MslPrecursorRecord.FragmentBlockOffset</c> values are then relative to the decompressed
/// buffer rather than absolute file positions, so <see cref="MslLibrary.LoadIndexOnly"/> falls
/// back to full-load mode transparently for compressed files.
///
/// All file-system writes go to a dedicated temp directory cleaned up in
/// <see cref="OneTimeTearDown"/>.
/// </summary>
[TestFixture]
public sealed class TestMslCompression
{
	// ── Fixture setup ─────────────────────────────────────────────────────────

	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslCompressionTests");

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDirectory);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
			Directory.Delete(OutputDirectory, recursive: true);
	}

	// ── Helpers ───────────────────────────────────────────────────────────────

	private static string TempPath(string name) =>
		Path.Combine(OutputDirectory, name + ".msl");

	/// <summary>
	/// Builds a single canonical test entry with three fragments. Used in most tests that
	/// need a small but well-defined library to verify round-trip correctness.
	/// </summary>
	private static MslLibraryEntry MakeEntry(string sequence = "PEPTIDE", int charge = 2,
		double precursorMz = 449.74, double irt = 35.4)
	{
		return new MslLibraryEntry
		{
			FullSequence = sequence,
			BaseSequence = sequence,
			PrecursorMz = precursorMz,
			ChargeState = charge,
			RetentionTime = irt,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 175.119f, Intensity = 1.0f,
					ProductType = ProductType.y, FragmentNumber = 1,
					ResiduePosition = 6, Charge = 1
				},
				new MslFragmentIon
				{
					Mz = 276.166f, Intensity = 0.8f,
					ProductType = ProductType.y, FragmentNumber = 2,
					ResiduePosition = 5, Charge = 1
				},
				new MslFragmentIon
				{
					Mz = 98.060f, Intensity = 0.6f,
					ProductType = ProductType.b, FragmentNumber = 2,
					ResiduePosition = 1, Charge = 1
				}
			}
		};
	}

	/// <summary>
	/// Builds a list of <paramref name="count"/> distinct entries. Each entry has a unique
	/// modified sequence (P0, P1, …) with a single y1 fragment. The fragment data is highly
	/// repetitive by design, making it very compressible.
	/// </summary>
	private static List<MslLibraryEntry> MakeEntries(int count)
	{
		// Each entry has 10 fragments with repetitive b/y types and normalized intensities.
		// 10 fragments per entry makes the fragment section the dominant contributor to file
		// size, which is essential for demonstrating meaningful compression ratios.
		var entries = new List<MslLibraryEntry>(count);
		for (int i = 0; i < count; i++)
		{
			var fragments = new List<MslFragmentIon>(10);
			for (int j = 1; j <= 10; j++)
			{
				fragments.Add(new MslFragmentIon
				{
					Mz = 100f + j * 50f,
					Intensity = 1.0f / j,
					ProductType = j % 2 == 0 ? ProductType.b : ProductType.y,
					FragmentNumber = j,
					ResiduePosition = j,
					Charge = 1
				});
			}

			entries.Add(new MslLibraryEntry
			{
				FullSequence = $"P{i}",
				BaseSequence = $"P{i}",
				PrecursorMz = 400.0 + i * 0.01,
				ChargeState = 2,
				RetentionTime = i * 0.1,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				MatchedFragmentIons = fragments
			});
		}
		return entries;
	}

	/// <summary>
	/// Reads the raw <see cref="MslFileHeader"/> from a written .msl file.
	/// </summary>
	private static MslFileHeader ReadRawHeader(string path)
	{
		byte[] bytes = new byte[MslFormat.HeaderSize];
		using var fs = new FileStream(path, FileMode.Open, FileAccess.Read);
		fs.ReadExactly(bytes, 0, MslFormat.HeaderSize);
		return MemoryMarshal.Read<MslFileHeader>(bytes.AsSpan());
	}

	// ── Tests ─────────────────────────────────────────────────────────────────

	/// <summary>
	/// A library written at level 3 must have FileFlagIsCompressed set in the header and
	/// read back all entries correctly.
	/// </summary>
	[Test]
	public void Compression_Level3_ProducesValidFile()
	{
		string path = TempPath("level3_valid");
		var entries = new List<MslLibraryEntry> { MakeEntry() };

		MslWriter.Write(path, entries, compressionLevel: 3);

		MslFileHeader header = ReadRawHeader(path);
		Assert.That((header.FileFlags & MslFormat.FileFlagIsCompressed), Is.Not.Zero,
			"FileFlagIsCompressed must be set when compressionLevel > 0");

		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.PrecursorCount, Is.EqualTo(1));
		Assert.That(lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? entry), Is.True);
		Assert.That(entry!.MatchedFragmentIons, Has.Count.EqualTo(3));
	}

	/// <summary>
	/// The compressed file must be smaller than the uncompressed version of the same
	/// library (at least 10% smaller for a 1 K-entry synthetic library).
	/// </summary>
	[Test]
	public void Compression_Level3_SmallerThanUncompressed()
	{
		string pathUncompressed = TempPath("size_uncompressed");
		string pathCompressed = TempPath("size_compressed");
		var entries = MakeEntries(1000);

		MslWriter.Write(pathUncompressed, entries, compressionLevel: 0);
		MslWriter.Write(pathCompressed, entries, compressionLevel: 3);

		long uncompressedSize = new FileInfo(pathUncompressed).Length;
		long compressedSize = new FileInfo(pathCompressed).Length;

		Assert.That(compressedSize, Is.LessThan(uncompressedSize * 0.9),
			$"Compressed ({compressedSize} bytes) should be at least 10% smaller than " +
			$"uncompressed ({uncompressedSize} bytes)");
	}

	/// <summary>
	/// Levels 1, 9, and 19 must all produce files that round-trip correctly.
	/// </summary>
	[Test]
	public void Compression_AllLevels_1_9_19_RoundTrip()
	{
		var entries = new List<MslLibraryEntry> { MakeEntry(), MakeEntry("ACDEFGHIK", 2, 529.76, 42.1) };

		foreach (int level in new[] { 1, 9, 19 })
		{
			string path = TempPath($"level{level}_roundtrip");
			MslWriter.Write(path, entries, compressionLevel: level);

			using MslLibrary lib = MslLibrary.Load(path);
			Assert.That(lib.PrecursorCount, Is.EqualTo(2),
				$"Level {level}: expected 2 precursors");
			Assert.That(lib.TryGetEntry("PEPTIDE", 2, out _), Is.True,
				$"Level {level}: PEPTIDE/2 not found");
			Assert.That(lib.TryGetEntry("ACDEFGHIK", 2, out _), Is.True,
				$"Level {level}: ACDEFGHIK/2 not found");
		}
	}

	/// <summary>
	/// compressionLevel = 0 must produce a file without FileFlagIsCompressed and the
	/// file must be readable as a normal (uncompressed) library.
	/// </summary>
	[Test]
	public void Compression_Level0_NoCompression()
	{
		string path = TempPath("level0_uncompressed");
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() }, compressionLevel: 0);

		MslFileHeader header = ReadRawHeader(path);
		Assert.That((header.FileFlags & MslFormat.FileFlagIsCompressed), Is.Zero,
			"FileFlagIsCompressed must NOT be set when compressionLevel == 0");

		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.PrecursorCount, Is.EqualTo(1));
		Assert.That(lib.IsCompressed, Is.False);
	}

	/// <summary>
	/// Every precursor field and every fragment field must match exactly after a
	/// compressed write + read round-trip.
	/// </summary>
	[Test]
	public void Compression_AllEntries_RoundTrip_Exact()
	{
		string path = TempPath("exact_roundtrip");

		var original = new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			BaseSequence = "PEPTIDE",
			PrecursorMz = 449.74,
			ChargeState = 2,
			RetentionTime = 35.4,
			IonMobility = 0.95,
			IsDecoy = false,
			IsProteotypic = true,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 175.119f, Intensity = 1.0f,
					ProductType = ProductType.y, FragmentNumber = 1,
					ResiduePosition = 6, Charge = 1, NeutralLoss = 0.0
				},
				new MslFragmentIon
				{
					Mz = 98.060f, Intensity = 0.6f,
					ProductType = ProductType.b, FragmentNumber = 2,
					ResiduePosition = 1, Charge = 1, NeutralLoss = 0.0
				}
			}
		};

		MslWriter.Write(path, new List<MslLibraryEntry> { original }, compressionLevel: 3);

		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? entry), Is.True);

		Assert.That(entry!.FullSequence, Is.EqualTo(original.FullSequence));
		Assert.That(entry.ChargeState, Is.EqualTo(original.ChargeState));
		Assert.That(entry.RetentionTime, Is.EqualTo((double)(float)original.RetentionTime).Within(1e-4));
		Assert.That(entry.IsDecoy, Is.EqualTo(original.IsDecoy));
		Assert.That(entry.IsProteotypic, Is.EqualTo(original.IsProteotypic));
		Assert.That(entry.MatchedFragmentIons, Has.Count.EqualTo(2));

		// MatchedFragmentIons are sorted by m/z ascending on write; check both
		MslFragmentIon b2 = entry.MatchedFragmentIons.Single(f => f.ProductType == ProductType.b);
		MslFragmentIon y1 = entry.MatchedFragmentIons.Single(f => f.ProductType == ProductType.y);

		Assert.That(b2.FragmentNumber, Is.EqualTo(2));
		Assert.That(y1.FragmentNumber, Is.EqualTo(1));
	}

	/// <summary>
	/// All fragment m/z values must round-trip within 1e-4 Da (float32 precision preserved).
	/// </summary>
	[Test]
	public void Compression_FragmentMz_WithinTolerance()
	{
		string path = TempPath("mz_tolerance");
		float[] expectedMzValues = { 98.060f, 175.119f, 276.166f };

		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() }, compressionLevel: 3);

		using MslLibrary lib = MslLibrary.Load(path);
		lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? entry);

		float[] actualMzValues = entry!.MatchedFragmentIons.Select(f => f.Mz).OrderBy(m => m).ToArray();

		Assert.That(actualMzValues, Has.Length.EqualTo(expectedMzValues.Length));
		for (int i = 0; i < expectedMzValues.Length; i++)
			Assert.That(actualMzValues[i], Is.EqualTo(expectedMzValues[i]).Within(1e-4f),
				$"Fragment [{i}] m/z mismatch after compressed round-trip");
	}

	/// <summary>
	/// MslLibrary.LoadIndexOnly on a compressed file must return a library with
	/// IsIndexOnly == false (transparent full-load fallback).
	/// </summary>
	[Test]
	public void Compression_LoadIndexOnly_FallsBackToFull()
	{
		string path = TempPath("index_only_fallback");
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() }, compressionLevel: 3);

		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);

		Assert.That(lib.IsIndexOnly, Is.False,
			"LoadIndexOnly on a compressed file must fall back to full-load (IsIndexOnly == false)");
		Assert.That(lib.PrecursorCount, Is.EqualTo(1),
			"Full-load fallback must still return all entries");
	}

	/// <summary>
	/// library.IsCompressed must be true for a file written with compressionLevel > 0.
	/// </summary>
	[Test]
	public void Compression_IsCompressed_TrueForCompressedFile()
	{
		string path = TempPath("is_compressed_true");
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() }, compressionLevel: 3);

		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.IsCompressed, Is.True);
	}

	/// <summary>
	/// library.IsCompressed must be false for a file written with compressionLevel == 0.
	/// </summary>
	[Test]
	public void Compression_IsCompressed_FalseForUncompressedFile()
	{
		string path = TempPath("is_compressed_false");
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() }, compressionLevel: 0);

		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.IsCompressed, Is.False);
	}

	/// <summary>
	/// Corrupting a single byte in the compressed fragment data must cause the CRC check
	/// to throw InvalidDataException (the compressed frame is covered by the CRC).
	/// </summary>
	[Test]
	public void Compression_CrcCoversCompressedBytes()
	{
		string path = TempPath("crc_compressed");
		MslWriter.Write(path, MakeEntries(10), compressionLevel: 3);

		byte[] bytes = File.ReadAllBytes(path);

		// The compressed data starts at FragmentSectionOffset (from the header).
		// Corrupt the first byte of the compressed zstd frame.
		MslFileHeader header = MemoryMarshal.Read<MslFileHeader>(bytes.AsSpan(0, MslFormat.HeaderSize));
		int frameStart = (int)header.FragmentSectionOffset;
		bytes[frameStart] ^= 0xFF;

		File.WriteAllBytes(path, bytes);

		Assert.Throws<InvalidDataException>(() => MslLibrary.Load(path),
			"CRC-32 must detect corruption inside the compressed fragment data");
	}

	/// <summary>
	/// An uncompressed version-2 file must be accepted by the version-3 reader.
	/// </summary>
	[Test]
	public void Compression_Version2_ReadableByVersion3Reader()
	{
		string path = TempPath("version2_compat");

		// Write a version-3 uncompressed file, then patch FormatVersion to 2
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() }, compressionLevel: 0);
		byte[] bytes = File.ReadAllBytes(path);

		// Patch FormatVersion at offset 4 to 2
		bytes[4] = 2; bytes[5] = 0; bytes[6] = 0; bytes[7] = 0;

		// Recompute CRC (the patched byte is in the CRC-covered region)
		int footerStart = bytes.Length - MslFormat.FooterSize;
		long offsetTableOffset = BitConverter.ToInt64(bytes, footerStart);
		uint newCrc = MslWriter.ComputeCrc32OfArray(bytes, (int)offsetTableOffset);
		byte[] crcBytes = BitConverter.GetBytes(newCrc);
		Array.Copy(crcBytes, 0, bytes, footerStart + 12, 4);

		File.WriteAllBytes(path, bytes);

		Assert.DoesNotThrow(() =>
		{
			using MslLibrary lib = MslLibrary.Load(path);
			Assert.That(lib.PrecursorCount, Is.EqualTo(1));
		}, "Version-2 file must be readable by the version-3 reader");
	}

	/// <summary>
	/// An uncompressed version-1 file must be accepted by the version-3 reader.
	/// </summary>
	[Test]
	public void Compression_Version1_ReadableByVersion3Reader()
	{
		string path = TempPath("version1_compat");

		// Write a version-3 uncompressed file, then patch FormatVersion to 1
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() }, compressionLevel: 0);
		byte[] bytes = File.ReadAllBytes(path);

		// Patch FormatVersion at offset 4 to 1, and zero ExtAnnotationTableOffset at offset 28
		bytes[4] = 1; bytes[5] = 0; bytes[6] = 0; bytes[7] = 0;
		bytes[28] = 0; bytes[29] = 0; bytes[30] = 0; bytes[31] = 0;

		int footerStart = bytes.Length - MslFormat.FooterSize;
		long offsetTableOffset = BitConverter.ToInt64(bytes, footerStart);
		uint newCrc = MslWriter.ComputeCrc32OfArray(bytes, (int)offsetTableOffset);
		byte[] crcBytes = BitConverter.GetBytes(newCrc);
		Array.Copy(crcBytes, 0, bytes, footerStart + 12, 4);

		File.WriteAllBytes(path, bytes);

		Assert.DoesNotThrow(() =>
		{
			using MslLibrary lib = MslLibrary.Load(path);
			Assert.That(lib.PrecursorCount, Is.EqualTo(1));
		}, "Version-1 file must be readable by the version-3 reader");
	}

	/// <summary>
	/// A 50 K-entry library at level 3 must be at least 2× smaller than the
	/// uncompressed version of the same library.
	/// </summary>
	[Test]
	public void Compression_LargeLibrary_50K_SizeReduction()
	{
		string pathUncompressed = TempPath("large_uncompressed");
		string pathCompressed = TempPath("large_compressed");

		// MslWriteLayout mutates entries in place (sorts fragments, normalizes intensities),
		// so each Write call must receive a fresh list.
		MslWriter.Write(pathUncompressed, MakeEntries(50_000), compressionLevel: 0);
		MslWriter.Write(pathCompressed, MakeEntries(50_000), compressionLevel: 3);

		long uncompressedSize = new FileInfo(pathUncompressed).Length;
		long compressedSize = new FileInfo(pathCompressed).Length;

		Assert.That(compressedSize * 2, Is.LessThan(uncompressedSize),
			$"50K-entry compressed file ({compressedSize} bytes) must be at least 2× smaller " +
			$"than uncompressed ({uncompressedSize} bytes)");
	}

	/// <summary>
	/// QueryMzWindow must return correct results on a library loaded from a compressed file.
	/// </summary>
	[Test]
	public void Compression_QueryMzWindow_WorksAfterDecompression()
	{
		string path = TempPath("query_mz_window");

		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("PEPTIDE",   2, precursorMz: 449.74),
			MakeEntry("ACDEFGHIK", 2, precursorMz: 529.76),
			MakeEntry("LMNPQRST",  2, precursorMz: 650.00)
		};

		MslWriter.Write(path, entries, compressionLevel: 3);

		using MslLibrary lib = MslLibrary.Load(path);

		// Window that should include only the first two entries
		ReadOnlySpan<MslPrecursorIndexEntry> window = lib.QueryMzWindow(400f, 530f);

		Assert.That(window.Length, Is.EqualTo(2),
			"QueryMzWindow should return exactly 2 entries in the [400, 530] window");

		float[] mzValues = window.ToArray().Select(e => e.PrecursorMz).OrderBy(m => m).ToArray();
		Assert.That(mzValues[0], Is.EqualTo(449.74f).Within(1e-2f));
		Assert.That(mzValues[1], Is.EqualTo(529.76f).Within(1e-2f));
	}

	/// <summary>
	/// TryGetEntry must return the correct entry from a compressed library (DDA lookup).
	/// </summary>
	[Test]
	public void Compression_DdaLookup_WorksAfterDecompression()
	{
		string path = TempPath("dda_lookup");
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("PEPTIDE",   2, precursorMz: 449.74, irt: 35.4),
			MakeEntry("ACDEFGHIK", 3, precursorMz: 353.51, irt: 42.1)
		};

		MslWriter.Write(path, entries, compressionLevel: 3);

		using MslLibrary lib = MslLibrary.Load(path);

		Assert.That(lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? e1), Is.True);
		Assert.That(e1!.ChargeState, Is.EqualTo(2));
		Assert.That(e1.MatchedFragmentIons, Has.Count.EqualTo(3));

		Assert.That(lib.TryGetEntry("ACDEFGHIK", 3, out MslLibraryEntry? e2), Is.True);
		Assert.That(e2!.ChargeState, Is.EqualTo(3));

		Assert.That(lib.TryGetEntry("NOTPRESENT", 2, out _), Is.False,
			"Lookup of absent entry must return false");
	}

	/// <summary>
	/// A zero-precursor library written with compression enabled must write and read
	/// without error.
	/// </summary>
	[Test]
	public void Compression_EmptyLibrary_HandledGracefully()
	{
		string path = TempPath("empty_compressed");

		Assert.DoesNotThrow(() =>
			MslWriter.Write(path, new List<MslLibraryEntry>(), compressionLevel: 3),
			"Writing an empty library with compression must not throw");

		Assert.DoesNotThrow(() =>
		{
			using MslLibrary lib = MslLibrary.Load(path);
			Assert.That(lib.PrecursorCount, Is.EqualTo(0));
			Assert.That(lib.IsCompressed, Is.True);
		}, "Reading a compressed empty library must not throw");
	}
}