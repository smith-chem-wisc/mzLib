using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 test class for <see cref="MslWriter"/>.
///
/// All tests that inspect file content read raw bytes and parse manually; no
/// <c>MslReader</c> is used (Prompt 3 deliverable). Each test that needs a file calls
/// <see cref="BuildTestEntries"/> or one of its siblings, writes to a per-test temp path,
/// and deletes the file in a teardown to keep the working directory clean.
///
/// Fixture layout:
/// <list type="bullet">
///   <item>Write produces a file — basic existence and magic checks</item>
///   <item>File structure — section offsets read from raw header bytes</item>
///   <item>Content correctness — first precursor fields parsed manually</item>
///   <item>String table — deduplication, unicode, index-0 empty string</item>
///   <item>Elution groups — same/different stripped sequence → same/different ID</item>
///   <item>Internal ions — secondary product type, residue numbers, flags byte</item>
///   <item>Validation — error conditions surfaced by ValidateEntries</item>
///   <item>CRC32 — footer checksum and corruption detection</item>
///   <item>WriteFromLibrarySpectra — round-trip through LibrarySpectrum</item>
/// </list>
/// </summary>
[TestFixture]
public sealed class TestMslWriter
{
	// ────────────────────────────────────────────────────────────────────────
	// Test fixture constants
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Directory into which all test output files are written. Using a dedicated subfolder
	/// avoids polluting the working directory and makes bulk teardown trivial.
	/// </summary>
	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslWriterTests");

	// ────────────────────────────────────────────────────────────────────────
	// Per-run setup / teardown
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>Creates the test output directory if it does not already exist.</summary>
	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		Directory.CreateDirectory(OutputDirectory);
	}

	/// <summary>Deletes every .msl file written by this test run.</summary>
	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
		{
			foreach (string file in Directory.GetFiles(OutputDirectory, "*.msl"))
				File.Delete(file);
		}
	}

	// ────────────────────────────────────────────────────────────────────────
	// Test data factories
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Builds a minimal list of two valid <see cref="MslLibraryEntry"/> objects suitable for
	/// use as a write fixture throughout the test suite.
	///
	/// Entry 0: PEPTIDE/2 — unmodified, two b/y fragments, protein P12345.
	/// Entry 1: PEPTIDE/3 — same stripped sequence as entry 0 (same elution group),
	///          one y fragment, no protein.
	///
	/// Both entries share the stripped sequence "PEPTIDE", so they must receive the same
	/// ElutionGroupId after the writer's pass-1 layout computation.
	/// </summary>
	/// <returns>A list of two fully populated MslLibraryEntry objects.</returns>
	private static List<MslLibraryEntry> BuildTestEntries()
	{
		// ── Entry 0: PEPTIDE/2 ───────────────────────────────────────────────
		var entry0 = new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			BaseSequence = "PEPTIDE",
			PrecursorMz = 449.7358,
			ChargeState = 2,
			RetentionTime = 35.4,
			IonMobility = 0.0,
			ProteinAccession = "P12345",
			ProteinName = "Test protein alpha",
			GeneName = "TESTA",
			IsDecoy = false,
			IsProteotypic = true,
			QValue = 0.01f,
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz                   = 175.119f,   // y1 of PEPTIDE
                    Intensity            = 8000f,
					ProductType          = ProductType.y,
					SecondaryProductType = null,
					FragmentNumber       = 1,
					SecondaryFragmentNumber = 0,
					ResiduePosition      = 6,
					Charge               = 1,
					NeutralLoss          = 0.0,
					ExcludeFromQuant     = false
				},
				new MslFragmentIon
				{
					Mz                   = 98.060f,    // b2 of PEPTIDE
                    Intensity            = 10000f,    // highest — will become 1.0 after normalize
                    ProductType          = ProductType.b,
					SecondaryProductType = null,
					FragmentNumber       = 2,
					SecondaryFragmentNumber = 0,
					ResiduePosition      = 1,
					Charge               = 1,
					NeutralLoss          = 0.0,
					ExcludeFromQuant     = false
				}
			}
		};

		// ── Entry 1: PEPTIDE/3 — same stripped sequence, different charge ────
		var entry1 = new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			BaseSequence = "PEPTIDE",
			PrecursorMz = 300.160,
			ChargeState = 3,
			RetentionTime = 35.4,
			IonMobility = 0.0,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsDecoy = false,
			IsProteotypic = false,
			QValue = float.NaN,
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz                   = 175.119f,
					Intensity            = 5000f,
					ProductType          = ProductType.y,
					SecondaryProductType = null,
					FragmentNumber       = 1,
					SecondaryFragmentNumber = 0,
					ResiduePosition      = 6,
					Charge               = 1,
					NeutralLoss          = 0.0,
					ExcludeFromQuant     = false
				}
			}
		};

		return new List<MslLibraryEntry> { entry0, entry1 };
	}

	/// <summary>
	/// Builds a single entry that contains one internal fragment ion (bIy type) for tests
	/// that specifically exercise the internal-ion code path.
	/// </summary>
	/// <returns>A list containing exactly one entry with one internal fragment.</returns>
	private static List<MslLibraryEntry> BuildInternalIonEntry()
	{
		var entry = new MslLibraryEntry
		{
			FullSequence = "ACDEFGHIK",
			BaseSequence = "ACDEFGHIK",
			PrecursorMz = 529.760,
			ChargeState = 2,
			RetentionTime = 42.1,
			IonMobility = 0.0,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsDecoy = false,
			IsProteotypic = true,
			QValue = float.NaN,
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
                // Terminal ion: y3 of ACDEFGHIK
                new MslFragmentIon
				{
					Mz                   = 364.197f,
					Intensity            = 9000f,
					ProductType          = ProductType.y,
					SecondaryProductType = null,
					FragmentNumber       = 3,
					SecondaryFragmentNumber = 0,
					ResiduePosition      = 6,
					Charge               = 1,
					NeutralLoss          = 0.0,
					ExcludeFromQuant     = false
				},
                // Internal ion: bIy spanning residues 2–5 (start=2, end=5)
                new MslFragmentIon
				{
					Mz                   = 485.215f,
					Intensity            = 10000f,
					ProductType          = ProductType.b,
					SecondaryProductType = ProductType.y,     // C-terminal type of the internal pair
                    FragmentNumber       = 2,                 // start residue index (0-based)
                    SecondaryFragmentNumber = 5,              // end residue index (0-based, exclusive)
                    ResiduePosition      = 2,
					Charge               = 1,
					NeutralLoss          = 0.0,
					ExcludeFromQuant     = false
				}
			}
		};

		return new List<MslLibraryEntry> { entry };
	}

	/// <summary>
	/// Builds a list of two entries with different stripped sequences ("PEPTIDE" and
	/// "ACDEFGHIK") to verify that different stripped sequences produce different elution
	/// group IDs.
	/// </summary>
	/// <returns>A list of two entries with distinct stripped sequences.</returns>
	private static List<MslLibraryEntry> BuildTwoDistinctSequenceEntries()
	{
		var entries = BuildTestEntries();            // entries[0] and [1] share "PEPTIDE"

		// Add a third entry with a different stripped sequence
		entries.Add(new MslLibraryEntry
		{
			FullSequence = "ACDEFGHIK",
			BaseSequence = "ACDEFGHIK",
			PrecursorMz = 529.760,
			ChargeState = 2,
			RetentionTime = 42.1,
			IonMobility = 0.0,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsDecoy = false,
			IsProteotypic = false,
			QValue = float.NaN,
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Nce = 0,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 175.119f, Intensity = 5000f,
					ProductType = ProductType.y, FragmentNumber = 1,
					ResiduePosition = 8, Charge = 1
				}
			}
		});

		return entries;
	}

	/// <summary>
	/// Returns a unique temp file path for a given test method name, rooted in
	/// <see cref="OutputDirectory"/>. The file will have the .msl extension.
	/// </summary>
	/// <param name="testName">Name of the calling test method (used as filename stem).</param>
	/// <returns>A full path to a non-existent .msl temp file.</returns>
	private static string TempPath(string testName) =>
		Path.Combine(OutputDirectory, testName + ".msl");

	// ────────────────────────────────────────────────────────────────────────
	// Raw byte-reading helpers
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Reads all bytes from a file at <paramref name="path"/> into a byte array.
	/// </summary>
	private static byte[] ReadAllBytes(string path) => File.ReadAllBytes(path);

	/// <summary>
	/// Reads a little-endian int32 from <paramref name="data"/> at the given byte
	/// <paramref name="offset"/>. Used to parse raw header / footer fields without an MslReader.
	/// </summary>
	private static int ReadInt32LE(byte[] data, int offset) =>
		data[offset]
		| (data[offset + 1] << 8)
		| (data[offset + 2] << 16)
		| (data[offset + 3] << 24);

	/// <summary>
	/// Reads a little-endian int64 from <paramref name="data"/> at the given byte
	/// <paramref name="offset"/>. Used to parse absolute file offsets stored in the header.
	/// </summary>
	private static long ReadInt64LE(byte[] data, int offset) =>
		(long)ReadInt32LE(data, offset)
		| ((long)ReadInt32LE(data, offset + 4) << 32);

	/// <summary>
	/// Reads a little-endian float32 from <paramref name="data"/> at the given byte offset.
	/// Used to parse PrecursorMz from the raw precursor record bytes.
	/// </summary>
	private static float ReadFloat32LE(byte[] data, int offset) =>
		BitConverter.ToSingle(data, offset);

	/// <summary>
	/// Reads a little-endian uint32 from <paramref name="data"/> at the given byte offset.
	/// Used to read magic and CRC fields.
	/// </summary>
	private static uint ReadUInt32LE(byte[] data, int offset) =>
		(uint)ReadInt32LE(data, offset);

	// ────────────────────────────────────────────────────────────────────────
	// Header field byte offsets (matching MslFileHeader layout)
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>Byte offset of MslFileHeader.Magic (uint, 4 bytes).</summary>
	private const int HdrOffMagic = 0;
	/// <summary>Byte offset of MslFileHeader.FormatVersion (int, 4 bytes).</summary>
	private const int HdrOffVersion = 4;
	/// <summary>Byte offset of MslFileHeader.NPrecursors (int, 4 bytes).</summary>
	private const int HdrOffNPrecursors = 12;
	/// <summary>Byte offset of MslFileHeader.ProteinTableOffset (long, 8 bytes).</summary>
	private const int HdrOffProteinTableOffset = 32;
	/// <summary>Byte offset of MslFileHeader.StringTableOffset (long, 8 bytes).</summary>
	private const int HdrOffStringTableOffset = 40;
	/// <summary>Byte offset of MslFileHeader.PrecursorSectionOffset (long, 8 bytes).</summary>
	private const int HdrOffPrecursorSectionOffset = 48;
	/// <summary>Byte offset of MslFileHeader.FragmentSectionOffset (long, 8 bytes).</summary>
	private const int HdrOffFragmentSectionOffset = 56;

	// ────────────────────────────────────────────────────────────────────────
	// Group: Write produces a file
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>Verifies that Write() creates a file at the specified path.</summary>
	[Test]
	public void Write_CreatesFile_WhenGivenValidEntries()
	{
		string path = TempPath(nameof(Write_CreatesFile_WhenGivenValidEntries));

		MslWriter.Write(path, BuildTestEntries());

		Assert.That(File.Exists(path), Is.True,
			"Write() must create a file at the output path.");
	}

	/// <summary>
	/// Verifies that the first four bytes of the written file equal the MZLB magic sequence.
	/// </summary>
	[Test]
	public void Write_FileHasCorrectMagicBytes()
	{
		string path = TempPath(nameof(Write_FileHasCorrectMagicBytes));
		MslWriter.Write(path, BuildTestEntries());

		byte[] data = ReadAllBytes(path);

		// Magic is written as a uint via BinaryWriter on a little-endian system;
		// on-disk bytes will be the reversed byte order of MagicAsUInt32 = 0x4D5A4C42u.
		// Use MslFormat.MagicMatches() which handles the byte-order question correctly.
		Assert.That(MslFormat.MagicMatches(data.AsSpan(0, 4)), Is.True,
			"First four bytes must match the MZLB magic sequence.");
	}

	/// <summary>
	/// Verifies that the last four bytes of the written file equal the trailing magic
	/// stored in MslFooter.TrailingMagic (same value as the header magic).
	/// </summary>
	[Test]
	public void Write_FileHasCorrectFooterMagic()
	{
		string path = TempPath(nameof(Write_FileHasCorrectFooterMagic));
		MslWriter.Write(path, BuildTestEntries());

		byte[] data = ReadAllBytes(path);

		// Footer layout: 8-byte OffsetTableOffset + 4-byte NPrecursors + 4-byte CRC + 4-byte magic
		// TrailingMagic is at the last 4 bytes (offset data.Length - 4)
		int footerMagicOffset = data.Length - 4;
		Assert.That(MslFormat.MagicMatches(data.AsSpan(footerMagicOffset, 4)), Is.True,
			"Last four bytes of the file must contain the MZLB trailing magic.");
	}

	/// <summary>
	/// Verifies that MslFooter.NPrecursors equals the number of entries passed to Write().
	/// </summary>
	[Test]
	public void Write_FooterNPrecursors_MatchesEntryCount()
	{
		var entries = BuildTestEntries();
		string path = TempPath(nameof(Write_FooterNPrecursors_MatchesEntryCount));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);

		// Footer layout offsets from EOF:
		//   -20: OffsetTableOffset (8 bytes)
		//   -12: NPrecursors       (4 bytes)
		//   -8:  DataCrc32         (4 bytes)
		//   -4:  TrailingMagic     (4 bytes)
		int nPrecursorsFooterOffset = data.Length - MslFormat.FooterSize + 8; // offset of NPrecursors in footer
		int storedNPrecursors = ReadInt32LE(data, nPrecursorsFooterOffset);

		Assert.That(storedNPrecursors, Is.EqualTo(entries.Count),
			"Footer NPrecursors must equal the number of entries written.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group: File structure tests (raw byte checks, no MslReader)
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Verifies that the PrecursorSectionOffset stored in the header points to a byte
	/// position consistent with the expected layout formula:
	///   HeaderSize + NProteins × ProteinRecordSize + StringTableSize.
	/// </summary>
	[Test]
	public void Write_PrecursorSection_StartsAtCorrectOffset()
	{
		var entries = BuildTestEntries();
		string path = TempPath(nameof(Write_PrecursorSection_StartsAtCorrectOffset));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);

		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);
		long stringTableOffset = ReadInt64LE(data, HdrOffStringTableOffset);

		// Precursor section must start after the string table
		Assert.That(precursorSectionOffset, Is.GreaterThan(stringTableOffset),
			"PrecursorSectionOffset must be after StringTableOffset.");

		// And must be consistent with header size + protein table
		long proteinTableOffset = ReadInt64LE(data, HdrOffProteinTableOffset);
		Assert.That(proteinTableOffset, Is.EqualTo(MslFormat.HeaderSize),
			"ProteinTableOffset must equal HeaderSize.");
	}

	/// <summary>
	/// Verifies that the FragmentSectionOffset in the header equals
	/// PrecursorSectionOffset + NPrecursors × PrecursorRecordSize.
	/// </summary>
	[Test]
	public void Write_FragmentSection_ImmediatelyAfterPrecursorSection()
	{
		var entries = BuildTestEntries();
		string path = TempPath(nameof(Write_FragmentSection_ImmediatelyAfterPrecursorSection));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);

		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);
		long fragmentSectionOffset = ReadInt64LE(data, HdrOffFragmentSectionOffset);
		int nPrecursors = ReadInt32LE(data, HdrOffNPrecursors);

		long expectedFragmentStart = precursorSectionOffset
									 + (long)nPrecursors * MslFormat.PrecursorRecordSize;

		Assert.That(fragmentSectionOffset, Is.EqualTo(expectedFragmentStart),
			"FragmentSectionOffset must be PrecursorSectionOffset + NPrecursors × PrecursorRecordSize.");
	}

	/// <summary>
	/// Verifies that MslFooter.OffsetTableOffset equals the byte position immediately
	/// following all fragment records, i.e. it is placed directly after the fragment section.
	/// </summary>
	[Test]
	public void Write_OffsetTable_ImmediatelyAfterFragmentSection()
	{
		var entries = BuildTestEntries();
		string path = TempPath(nameof(Write_OffsetTable_ImmediatelyAfterFragmentSection));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);

		long fragmentSectionOffset = ReadInt64LE(data, HdrOffFragmentSectionOffset);

		// Sum all fragment counts to compute total fragment bytes
		long totalFragmentBytes = 0;
		foreach (MslLibraryEntry e in entries)
			totalFragmentBytes += (e.MatchedFragmentIons?.Count ?? 0) * MslFormat.FragmentRecordSize;

		long expectedOffsetTableStart = fragmentSectionOffset + totalFragmentBytes;

		// Read OffsetTableOffset from footer (first 8 bytes of footer)
		long offsetTableOffset = ReadInt64LE(data, data.Length - MslFormat.FooterSize);

		Assert.That(offsetTableOffset, Is.EqualTo(expectedOffsetTableStart),
			"OffsetTableOffset must be immediately after all fragment records.");
	}

	/// <summary>
	/// Verifies that EstimateFileSize() returns a value within ±5% of the actual file size
	/// for the standard two-entry test fixture.
	/// </summary>
	[Test]
	public void Write_FileSizeMatchesEstimate()
	{
		var entries = BuildTestEntries();
		string path = TempPath(nameof(Write_FileSizeMatchesEstimate));
		MslWriter.Write(path, entries);

		long actualSize = new FileInfo(path).Length;
		int avgFragments = entries.Sum(e => e.MatchedFragmentIons?.Count ?? 0) / entries.Count;
		long estimate = MslWriter.EstimateFileSize(entries.Count, avgFragments, 20);

		double ratio = (double)actualSize / estimate;

		// EstimateFileSize is a heuristic pre-flight check documented as ±5% for 1,000-entry
		// libraries. For a two-entry fixture with heavy string deduplication the variance is
		// higher, so we allow ±15% here. The important property is that the estimate is
		// in the right order of magnitude — it should never be off by a factor of 2×.
		Assert.That(ratio, Is.InRange(0.85, 1.15),
			$"EstimateFileSize should be within ±15% of actual size for small fixtures. " +
			$"Actual={actualSize}, Estimate={estimate}, ratio={ratio:F3}.");
	}

	// 

	// ────────────────────────────────────────────────────────────────────────
	// Group: Content correctness (raw byte parsing)
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Verifies that the PrecursorMz stored in the first precursor record on disk matches
	/// the value supplied in the first entry (within float32 precision).
	/// </summary>
	[Test]
	public void Write_FirstPrecursor_HasCorrectMz()
	{
		var entries = BuildTestEntries();
		string path = TempPath(nameof(Write_FirstPrecursor_HasCorrectMz));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);

		// MslPrecursorRecord.PrecursorMz is at offset 0 of the first record (float32)
		float storedMz = ReadFloat32LE(data, (int)precursorSectionOffset);

		Assert.That(storedMz, Is.EqualTo((float)entries[0].PrecursorMz).Within(1e-4f),
			"First precursor PrecursorMz on disk must match the input entry.");
	}

	/// <summary>
	/// Verifies that the ChargeState stored in the first precursor record matches the input entry.
	/// ChargeState is at offset 12 of MslPrecursorRecord as an int16.
	/// </summary>
	[Test]
	public void Write_FirstPrecursor_HasCorrectCharge()
	{
		var entries = BuildTestEntries();
		string path = TempPath(nameof(Write_FirstPrecursor_HasCorrectCharge));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);

		// MslPrecursorRecord.ChargeState is at offset 12 from the start of the record (int16, 2 bytes)
		const int ChargeOffsetInRecord = 12;
		short storedCharge = (short)(data[(int)precursorSectionOffset + ChargeOffsetInRecord]
									 | (data[(int)precursorSectionOffset + ChargeOffsetInRecord + 1] << 8));

		Assert.That(storedCharge, Is.EqualTo(entries[0].ChargeState),
			"First precursor ChargeState on disk must match the input entry.");
	}

	/// <summary>
	/// Verifies that the FragmentCount stored in the first precursor record matches the
	/// fragment count of the first input entry.
	/// FragmentCount is at offset 14 of MslPrecursorRecord as an int16.
	/// </summary>
	[Test]
	public void Write_FirstPrecursor_HasCorrectFragmentCount()
	{
		var entries = BuildTestEntries();
		string path = TempPath(nameof(Write_FirstPrecursor_HasCorrectFragmentCount));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);

		// MslPrecursorRecord.FragmentCount is at offset 14 from the start of the record (int16)
		const int FragCountOffsetInRecord = 14;
		short storedFragCount = (short)(
			data[(int)precursorSectionOffset + FragCountOffsetInRecord]
			| (data[(int)precursorSectionOffset + FragCountOffsetInRecord + 1] << 8));

		Assert.That(storedFragCount, Is.EqualTo(entries[0].MatchedFragmentIons.Count),
			"First precursor FragmentCount on disk must equal the number of fragments in the entry.");
	}

	/// <summary>
	/// Verifies that the FragmentBlockOffset stored in the first precursor record points to
	/// a byte position that contains a valid fragment m/z value (the first fragment of
	/// entry 0, which after m/z sort should be the lower of the two m/z values).
	/// FragmentBlockOffset is at offset 32 of MslPrecursorRecord as an int64.
	/// </summary>
	[Test]
	public void Write_FragmentBlockOffset_PointsToCorrectFragmentData()
	{
		var entries = BuildTestEntries();
		string path = TempPath(nameof(Write_FragmentBlockOffset_PointsToCorrectFragmentData));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);

		// MslPrecursorRecord.FragmentBlockOffset is at offset 32 of the first precursor record
		const int FragBlockOffsetInRecord = 32;
		long fragmentBlockOffset = ReadInt64LE(data, (int)precursorSectionOffset + FragBlockOffsetInRecord);

		// First fragment m/z stored at byte 0 of MslFragmentRecord is a float32
		float firstFragMz = ReadFloat32LE(data, (int)fragmentBlockOffset);

		// After m/z sort the smallest Mz in entry 0 is 98.060 (b2); verify it is close
		float expectedMinMz = entry0MinMz(entries);
		Assert.That(firstFragMz, Is.EqualTo(expectedMinMz).Within(0.01f),
			"FragmentBlockOffset must point to the first fragment (m/z sorted) of the precursor.");
	}

	/// <summary>Returns the minimum Mz value among all fragments in the first entry.</summary>
	private static float entry0MinMz(IReadOnlyList<MslLibraryEntry> entries) =>
		entries[0].MatchedFragmentIons.Min(f => f.Mz);

	/// <summary>
	/// Verifies that the fragment records on disk for the first precursor are ordered by
	/// m/z ascending (pass-1 sorts fragments before writing).
	/// </summary>
	[Test]
	public void Write_Fragments_AreSortedByMz()
	{
		var entries = BuildTestEntries();
		string path = TempPath(nameof(Write_Fragments_AreSortedByMz));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);
		const int FragBlockOffsetInRecord = 32;
		long fragmentBlockOffset = ReadInt64LE(data, (int)precursorSectionOffset + FragBlockOffsetInRecord);

		int fragCount = entries[0].MatchedFragmentIons.Count;
		float prevMz = float.MinValue;

		for (int i = 0; i < fragCount; i++)
		{
			// MslFragmentRecord.Mz is at byte 0 of each 20-byte record
			long recordStart = fragmentBlockOffset + (long)i * MslFormat.FragmentRecordSize;
			float mz = ReadFloat32LE(data, (int)recordStart);

			Assert.That(mz, Is.GreaterThanOrEqualTo(prevMz),
				$"Fragment [{i}] m/z ({mz}) must be ≥ previous m/z ({prevMz}). " +
				"MatchedFragmentIons must be sorted by m/z ascending.");
			prevMz = mz;
		}
	}

	/// <summary>
	/// Verifies that after writing, the maximum intensity in the first precursor's
	/// fragment block is exactly 1.0 (the result of pass-1 normalization).
	/// Intensity is at byte offset 4 of each MslFragmentRecord (float32).
	/// </summary>
	[Test]
	public void Write_IntensitiesAreNormalized_MaxIsOne()
	{
		var entries = BuildTestEntries();
		string path = TempPath(nameof(Write_IntensitiesAreNormalized_MaxIsOne));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);
		const int FragBlockOffsetInRecord = 32;
		long fragmentBlockOffset = ReadInt64LE(data, (int)precursorSectionOffset + FragBlockOffsetInRecord);

		int fragCount = entries[0].MatchedFragmentIons.Count;
		float maxIntensity = float.MinValue;

		for (int i = 0; i < fragCount; i++)
		{
			// MslFragmentRecord.Intensity is at byte offset 4 within the record (float32)
			long recordStart = fragmentBlockOffset + (long)i * MslFormat.FragmentRecordSize;
			float intensity = ReadFloat32LE(data, (int)recordStart + 4);
			if (intensity > maxIntensity) maxIntensity = intensity;
		}

		Assert.That(maxIntensity, Is.EqualTo(1.0f).Within(1e-6f),
			"The maximum fragment intensity within a precursor must be 1.0 after normalization.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group: String table tests
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Verifies that the string table on disk contains the empty string at index 0.
	/// Index 0 corresponds to the first entry after the two 4-byte header ints; its
	/// length prefix should be 0.
	/// </summary>
	[Test]
	public void Write_StringTable_ContainsEmptyStringAtIndex0()
	{
		string path = TempPath(nameof(Write_StringTable_ContainsEmptyStringAtIndex0));
		MslWriter.Write(path, BuildTestEntries());

		byte[] data = ReadAllBytes(path);
		long stringTableOffset = ReadInt64LE(data, HdrOffStringTableOffset);

		// String table layout: int32 NStrings, int32 TotalBytes, then entries
		// Entry 0 starts at stringTableOffset + 8 (after the two header ints)
		int entry0LengthOffset = (int)stringTableOffset + 8;
		int entry0Length = ReadInt32LE(data, entry0LengthOffset);

		Assert.That(entry0Length, Is.EqualTo(0),
			"String table entry at index 0 must be the empty string (length prefix = 0).");
	}

	/// <summary>
	/// Verifies that identical strings in two different entries are stored only once in the
	/// string table (deduplication). The test fixture has both entries sharing "PEPTIDE"
	/// as stripped and modified sequence; the string table should contain "PEPTIDE" only once.
	/// </summary>
	[Test]
	public void Write_StringTable_DeduplicatesIdenticalStrings()
	{
		string path = TempPath(nameof(Write_StringTable_DeduplicatesIdenticalStrings));
		MslWriter.Write(path, BuildTestEntries());

		byte[] data = ReadAllBytes(path);
		long stringTableOffset = ReadInt64LE(data, HdrOffStringTableOffset);

		int nStrings = ReadInt32LE(data, (int)stringTableOffset);

		// Parse all strings from the string table
		int pos = (int)stringTableOffset + 8; // skip NStrings + TotalBytes headers
		var strings = new List<string>(nStrings);

		for (int i = 0; i < nStrings; i++)
		{
			int len = ReadInt32LE(data, pos);
			pos += 4;
			string s = Encoding.UTF8.GetString(data, pos, len);
			strings.Add(s);
			pos += len;
		}

		int peptideCount = strings.Count(s => s == "PEPTIDE");
		Assert.That(peptideCount, Is.EqualTo(1),
			"The string 'PEPTIDE' must appear exactly once in the string table (deduplication).");
	}

	/// <summary>
	/// Verifies that a modified sequence containing non-ASCII characters (e.g. a Unicode
	/// sequence label) round-trips correctly through the string table UTF-8 encoding.
	/// </summary>
	[Test]
	public void Write_StringTable_UnicodeSequenceRoundTrips()
	{
		// Construct an entry with a Unicode-containing modified sequence
		const string UnicodeSequence = "PΕPTIDΕ"; // two Greek epsilon characters

		var entries = new List<MslLibraryEntry>
		{
			new MslLibraryEntry
			{
				FullSequence = UnicodeSequence,
				BaseSequence = "PEPTIDE",
				PrecursorMz      = 449.7,
				ChargeState           = 2,
				Source           = MslFormat.SourceType.Predicted,
				MoleculeType     = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				MatchedFragmentIons        = new List<MslFragmentIon>
				{
					new MslFragmentIon
					{
						Mz = 175.119f, Intensity = 1.0f,
						ProductType = ProductType.y, FragmentNumber = 1,
						ResiduePosition = 6, Charge = 1
					}
				}
			}
		};

		string path = TempPath(nameof(Write_StringTable_UnicodeSequenceRoundTrips));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long stringTableOffset = ReadInt64LE(data, HdrOffStringTableOffset);

		int nStrings = ReadInt32LE(data, (int)stringTableOffset);
		int pos = (int)stringTableOffset + 8;
		bool found = false;

		for (int i = 0; i < nStrings; i++)
		{
			int len = ReadInt32LE(data, pos);
			pos += 4;
			string s = Encoding.UTF8.GetString(data, pos, len);
			pos += len;

			if (s == UnicodeSequence) { found = true; break; }
		}

		Assert.That(found, Is.True,
			$"The Unicode sequence '{UnicodeSequence}' must be recoverable from the string table.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group: Elution group tests
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Verifies that two precursors with the same stripped sequence receive the same
	/// ElutionGroupId in their on-disk precursor records.
	/// ElutionGroupId is at offset 16 of MslPrecursorRecord (int32).
	/// </summary>
	[Test]
	public void Write_SameStrippedSequence_SameElutionGroupId()
	{
		var entries = BuildTestEntries(); // both entries have BaseSequence = "PEPTIDE"
		string path = TempPath(nameof(Write_SameStrippedSequence_SameElutionGroupId));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);

		// ElutionGroupId is at byte offset 16 within MslPrecursorRecord (int32)
		const int ElutionGroupIdOffset = 16;

		int egId0 = ReadInt32LE(data, (int)precursorSectionOffset + ElutionGroupIdOffset);
		int egId1 = ReadInt32LE(data, (int)precursorSectionOffset
									   + MslFormat.PrecursorRecordSize
									   + ElutionGroupIdOffset);

		Assert.That(egId0, Is.EqualTo(egId1),
			"Precursors with the same BaseSequence must have the same ElutionGroupId.");
	}

	/// <summary>
	/// Verifies that two precursors with different stripped sequences receive different
	/// ElutionGroupId values in their on-disk precursor records.
	/// </summary>
	[Test]
	public void Write_DifferentStrippedSequences_DifferentElutionGroupIds()
	{
		var entries = BuildTwoDistinctSequenceEntries(); // [0] and [1]: "PEPTIDE", [2]: "ACDEFGHIK"
		string path = TempPath(nameof(Write_DifferentStrippedSequences_DifferentElutionGroupIds));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);

		const int ElutionGroupIdOffset = 16;

		// Entry 0 ("PEPTIDE") and entry 2 ("ACDEFGHIK") should have different group IDs
		int egId0 = ReadInt32LE(data, (int)precursorSectionOffset + ElutionGroupIdOffset);
		int egId2 = ReadInt32LE(data, (int)precursorSectionOffset
									   + 2 * MslFormat.PrecursorRecordSize
									   + ElutionGroupIdOffset);

		Assert.That(egId0, Is.Not.EqualTo(egId2),
			"Precursors with different StrippedSequences must have different ElutionGroupIds.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group: Internal ion tests
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Verifies that when an internal ion is written, the SecondaryProductType field in the
	/// fragment record is not −1 (i.e. a real product type code is stored).
	/// SecondaryProductType is at byte offset 10 of MslFragmentRecord (int16).
	/// </summary>
	[Test]
	public void Write_InternalIon_SecondaryProductType_IsStored()
	{
		var entries = BuildInternalIonEntry();
		string path = TempPath(nameof(Write_InternalIon_SecondaryProductType_IsStored));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);

		const int FragBlockOffsetInRecord = 32;
		long fragmentBlockOffset = ReadInt64LE(data, (int)precursorSectionOffset + FragBlockOffsetInRecord);

		// After m/z sort the internal ion (Mz=485.215) is the second record (index 1)
		long internalRecordStart = fragmentBlockOffset + MslFormat.FragmentRecordSize; // record [1]

		// SecondaryProductType is at offset 10 within MslFragmentRecord (int16)
		const int SecondaryProductTypeOffsetInRecord = 10;
		short secondaryType = (short)(
			data[(int)internalRecordStart + SecondaryProductTypeOffsetInRecord]
			| (data[(int)internalRecordStart + SecondaryProductTypeOffsetInRecord + 1] << 8));

		Assert.That(secondaryType, Is.Not.EqualTo(-1),
			"An internal ion must have a valid SecondaryProductType (not -1) stored on disk.");
	}

	/// <summary>
	/// Verifies that the start residue (FragmentNumber) and end residue (SecondaryFragmentNumber)
	/// of an internal ion are stored correctly in the on-disk fragment record.
	/// FragmentNumber at offset 12, SecondaryFragmentNumber at offset 14 (both int16).
	/// </summary>
	[Test]
	public void Write_InternalIon_StartAndEndResidue_AreStored()
	{
		var entries = BuildInternalIonEntry();
		string path = TempPath(nameof(Write_InternalIon_StartAndEndResidue_AreStored));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);

		const int FragBlockOffsetInRecord = 32;
		long fragmentBlockOffset = ReadInt64LE(data, (int)precursorSectionOffset + FragBlockOffsetInRecord);

		// After m/z sort internal ion (Mz=485.215) is at index 1
		long internalRecordStart = fragmentBlockOffset + MslFormat.FragmentRecordSize;

		const int FragmentNumberOffset = 12; // int16
		const int SecondaryFragmentNumberOffset = 14; // int16

		short storedStart = (short)(data[(int)internalRecordStart + FragmentNumberOffset]
									| (data[(int)internalRecordStart + FragmentNumberOffset + 1] << 8));
		short storedEnd = (short)(data[(int)internalRecordStart + SecondaryFragmentNumberOffset]
									| (data[(int)internalRecordStart + SecondaryFragmentNumberOffset + 1] << 8));

		// The internal ion was created with FragmentNumber=2 and SecondaryFragmentNumber=5
		Assert.That(storedStart, Is.EqualTo(2),
			"Internal ion FragmentNumber (start residue) must match the input.");
		Assert.That(storedEnd, Is.EqualTo(5),
			"Internal ion SecondaryFragmentNumber (end residue) must match the input.");
	}

	/// <summary>
	/// Verifies that the flags byte (offset 19 of MslFragmentRecord) has bit 0 set for
	/// an internal fragment ion.
	/// </summary>
	[Test]
	public void Write_InternalIon_FlagsByte_IsInternalBitSet()
	{
		var entries = BuildInternalIonEntry();
		string path = TempPath(nameof(Write_InternalIon_FlagsByte_IsInternalBitSet));
		MslWriter.Write(path, entries);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);

		const int FragBlockOffsetInRecord = 32;
		long fragmentBlockOffset = ReadInt64LE(data, (int)precursorSectionOffset + FragBlockOffsetInRecord);

		// After m/z sort, internal ion is at index 1
		long internalRecordStart = fragmentBlockOffset + MslFormat.FragmentRecordSize;

		// Flags byte is at offset 19 within MslFragmentRecord
		const int FlagsOffset = 19;
		const byte IsInternalBit = 0b_0000_0001; // bit 0

		byte flags = data[(int)internalRecordStart + FlagsOffset];

		Assert.That((flags & IsInternalBit), Is.Not.Zero,
			"Internal fragment ion must have bit 0 (is_internal) set in the flags byte.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group: Validation tests
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Verifies that ValidateEntries() returns at least one error when an entry has an
	/// empty FullSequence.
	/// </summary>
	[Test]
	public void ValidateEntries_EmptySequence_ReturnsError()
	{
		var entries = new List<MslLibraryEntry>
		{
			new MslLibraryEntry
			{
				FullSequence = string.Empty, // invalid
                ChargeState           = 2,
				PrecursorMz      = 449.7,
				MatchedFragmentIons        = new List<MslFragmentIon>
				{
					new MslFragmentIon { Mz = 175f, Intensity = 1f, ProductType = ProductType.y,
										 FragmentNumber = 1, Charge = 1 }
				}
			}
		};

		List<string> errors = MslWriter.ValidateEntries(entries);

		Assert.That(errors, Is.Not.Empty,
			"An entry with an empty FullSequence must produce a validation error.");
	}

	/// <summary>
	/// Verifies that ValidateEntries() returns at least one error when an entry has ChargeState = 0.
	/// </summary>
	[Test]
	public void ValidateEntries_ZeroCharge_ReturnsError()
	{
		var entries = new List<MslLibraryEntry>
		{
			new MslLibraryEntry
			{
				FullSequence = "PEPTIDE",
				ChargeState           = 0,          // invalid
                PrecursorMz      = 449.7,
				MatchedFragmentIons        = new List<MslFragmentIon>
				{
					new MslFragmentIon { Mz = 175f, Intensity = 1f, ProductType = ProductType.y,
										 FragmentNumber = 1, Charge = 1 }
				}
			}
		};

		List<string> errors = MslWriter.ValidateEntries(entries);

		Assert.That(errors, Is.Not.Empty,
			"An entry with ChargeState = 0 must produce a validation error.");
	}

	/// <summary>
	/// Verifies that ValidateEntries() returns an error when an internal ion has
	/// SecondaryFragmentNumber ≤ FragmentNumber (invalid residue range).
	/// </summary>
	[Test]
	public void ValidateEntries_InternalIon_EndLessThanStart_ReturnsError()
	{
		var entries = new List<MslLibraryEntry>
		{
			new MslLibraryEntry
			{
				FullSequence = "ACDEFGHIK",
				ChargeState           = 2,
				PrecursorMz      = 529.76,
				MatchedFragmentIons        = new List<MslFragmentIon>
				{
					new MslFragmentIon
					{
						Mz                      = 300f,
						Intensity               = 1f,
						ProductType             = ProductType.b,
						SecondaryProductType    = ProductType.y, // marks as internal
                        FragmentNumber          = 5,
						SecondaryFragmentNumber = 3,             // invalid: end <= start
                        Charge                  = 1
					}
				}
			}
		};

		List<string> errors = MslWriter.ValidateEntries(entries);

		Assert.That(errors, Is.Not.Empty,
			"An internal ion with SecondaryFragmentNumber ≤ FragmentNumber must produce a validation error.");
	}

	/// <summary>
	/// Verifies that ValidateEntries() returns an empty list for the canonical test entries.
	/// </summary>
	[Test]
	public void ValidateEntries_ValidEntries_ReturnsEmptyList()
	{
		List<string> errors = MslWriter.ValidateEntries(BuildTestEntries());

		Assert.That(errors, Is.Empty,
			"ValidateEntries() must return an empty list for a fully valid entry set.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group: CRC32 tests
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Verifies that the CRC-32 stored in the footer matches an independently computed
	/// CRC-32 over bytes 0..(OffsetTableOffset−1) of the written file.
	/// </summary>
	[Test]
	public void Write_CRC32_IsCorrectInFooter()
	{
		string path = TempPath(nameof(Write_CRC32_IsCorrectInFooter));
		MslWriter.Write(path, BuildTestEntries());

		byte[] data = ReadAllBytes(path);

		// Read OffsetTableOffset from footer (first 8 bytes of footer)
		long offsetTableOffset = ReadInt64LE(data, data.Length - MslFormat.FooterSize);

		// Independently compute CRC-32 over data[0..(offsetTableOffset-1)]
		// Uses the same self-contained CRC-32/ISO-HDLC table as the writer.
		uint expectedCrc = MslWriter.ComputeCrc32OfArray(data, (int)offsetTableOffset);

		// Read the stored CRC from footer: at offset (FooterSize - 8) from start of footer
		// Footer layout: 8 OffsetTableOffset + 4 NPrecursors + 4 CRC + 4 TrailingMagic
		int storedCrcOffset = data.Length - MslFormat.FooterSize + 8 + 4; // = offset of DataCrc32
		uint storedCrc = ReadUInt32LE(data, storedCrcOffset);

		Assert.That(storedCrc, Is.EqualTo(expectedCrc),
			"The CRC-32 in the footer must equal an independently computed CRC of the data bytes.");
	}

	/// <summary>
	/// Verifies that corrupting a single data byte causes the stored CRC-32 to no longer
	/// match an independently computed CRC over the corrupted data.
	/// This does not test the reader — it simply confirms that the CRC is sensitive to changes.
	/// </summary>
	[Test]
	public void Write_CorruptedByte_WouldFailCRC32()
	{
		string path = TempPath(nameof(Write_CorruptedByte_WouldFailCRC32));
		MslWriter.Write(path, BuildTestEntries());

		byte[] data = ReadAllBytes(path);

		// Read OffsetTableOffset to know the CRC boundary
		long offsetTableOffset = ReadInt64LE(data, data.Length - MslFormat.FooterSize);

		// Read the stored CRC
		int storedCrcOffset = data.Length - MslFormat.FooterSize + 12;
		uint storedCrc = ReadUInt32LE(data, storedCrcOffset);

		// Corrupt a byte in the middle of the data section (byte 32, inside the header)
		byte[] corrupted = (byte[])data.Clone();
		corrupted[32] ^= 0xFF; // flip all bits of one byte

		// Recompute CRC over corrupted data using the same algorithm as the writer
		uint corruptedCrc = MslWriter.ComputeCrc32OfArray(corrupted, (int)offsetTableOffset);

		Assert.That(corruptedCrc, Is.Not.EqualTo(storedCrc),
			"Corrupting a data byte must cause the CRC to differ from the stored footer value.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group: WriteFromLibrarySpectra
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Verifies that WriteFromLibrarySpectra() produces a valid .msl file: magic bytes are
	/// correct and the file is non-empty.
	/// </summary>
	[Test]
	public void WriteFromLibrarySpectra_ProducesReadableFile()
	{
		// Build LibrarySpectrum objects using the mzLib constructor
		var matchedIons = new List<MatchedFragmentIon>
		{
			new MatchedFragmentIon(
				new Product(ProductType.y, FragmentationTerminus.C, 0.0, 1, 6, 0.0),
				experMz: 175.119, experIntensity: 1.0, charge: 1),
			new MatchedFragmentIon(
				new Product(ProductType.b, FragmentationTerminus.N, 0.0, 2, 1, 0.0),
				experMz: 98.060, experIntensity: 0.8, charge: 1)
		};

		var spectra = new List<LibrarySpectrum>
		{
			new LibrarySpectrum("PEPTIDE", 449.74, 2, matchedIons, rt: 35.4)
		};

		string path = TempPath(nameof(WriteFromLibrarySpectra_ProducesReadableFile));
		MslWriter.WriteFromLibrarySpectra(path, spectra);

		Assert.That(File.Exists(path), Is.True, "WriteFromLibrarySpectra must create a file.");

		byte[] data = ReadAllBytes(path);
		Assert.That(MslFormat.MagicMatches(data.AsSpan(0, 4)), Is.True,
			"File written by WriteFromLibrarySpectra must have correct magic bytes.");
	}

	/// <summary>
	/// Verifies that WriteFromLibrarySpectra() preserves fragment m/z values to within
	/// float32 precision (conversion path: double MatchedFragmentIon.Mz → float MslFragmentIon.Mz
	/// → float32 on disk).
	/// </summary>
	[Test]
	public void WriteFromLibrarySpectra_PreservesFragmentMzValues()
	{
		const double TargetMz = 175.119;

		var matchedIons = new List<MatchedFragmentIon>
		{
			new MatchedFragmentIon(
				new Product(ProductType.y, FragmentationTerminus.C, 0.0, 1, 6, 0.0),
				experMz: TargetMz, experIntensity: 1.0, charge: 1)
		};

		var spectra = new List<LibrarySpectrum>
		{
			new LibrarySpectrum("PEPTIDE", 449.74, 2, matchedIons, rt: 35.4)
		};

		string path = TempPath(nameof(WriteFromLibrarySpectra_PreservesFragmentMzValues));
		MslWriter.WriteFromLibrarySpectra(path, spectra);

		byte[] data = ReadAllBytes(path);
		long precursorSectionOffset = ReadInt64LE(data, HdrOffPrecursorSectionOffset);

		const int FragBlockOffsetInRecord = 32;
		long fragmentBlockOffset = ReadInt64LE(data, (int)precursorSectionOffset + FragBlockOffsetInRecord);

		// MslFragmentRecord.Mz is at byte offset 0 (float32)
		float storedMz = ReadFloat32LE(data, (int)fragmentBlockOffset);

		Assert.That(storedMz, Is.EqualTo((float)TargetMz).Within(1e-3f),
			"Fragment m/z must be preserved to within float32 precision through the conversion path.");
	}
}