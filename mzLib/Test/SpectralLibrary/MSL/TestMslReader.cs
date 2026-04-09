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
using System.Threading.Tasks;

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 test class for <see cref="MslReader"/>.
///
/// Tests are organized into groups:
/// <list type="bullet">
///   <item>Validation — wrong magic, wrong version, truncated footer, NPrecursors mismatch, CRC corruption</item>
///   <item>Header-only read — ReadHeaderOnly returns correct metadata</item>
///   <item>Full load — precursor fields, fragment fields, protein data round-trip</item>
///   <item>Index-only load — lazy fragment retrieval, thread safety, Dispose behavior</item>
///   <item>Round-trip fidelity — all fields identical after Write → Load cycle</item>
///   <item>Edge cases — empty library, zero fragments, large fragment counts</item>
/// </list>
///
/// All test files are written to a dedicated temp directory using <see cref="MslWriter.Write"/>
/// on the same synthetic fixtures used in <c>TestMslWriter</c>. The writer is treated as a
/// black box; its correctness is assumed from the Prompt 2 test suite.
/// </summary>
[TestFixture]
public sealed class TestMslReader
{
	// ── Output directory ──────────────────────────────────────────────────────

	/// <summary>
	/// Temporary directory for all test output files. Using a dedicated subfolder
	/// keeps the system temp directory tidy and makes bulk teardown easy.
	/// </summary>
	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslReaderTests");

	// ── Fixture lifecycle ─────────────────────────────────────────────────────

	/// <summary>Creates the test output directory once before any test in this fixture runs.</summary>
	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		Directory.CreateDirectory(OutputDirectory);
	}

	/// <summary>Deletes every .msl test file after all tests in this fixture have run.</summary>
	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
		{
			foreach (string file in Directory.GetFiles(OutputDirectory, "*.msl"))
				File.Delete(file);
			foreach (string file in Directory.GetFiles(OutputDirectory, "*.tmp"))
				File.Delete(file);
		}
	}

	// ── Path helper ───────────────────────────────────────────────────────────

	/// <summary>
	/// Returns a unique temp file path inside <see cref="OutputDirectory"/> for a test.
	/// Using the test method name as the filename makes failed-test artifacts easy to locate.
	/// </summary>
	/// <param name="testName">The test method name (typically <c>nameof(...)</c>).</param>
	/// <returns>A fully qualified file path with a .msl extension.</returns>
	private static string TempPath(string testName) =>
		Path.Combine(OutputDirectory, testName + ".msl");

	// ── Test data factories ───────────────────────────────────────────────────

	/// <summary>
	/// Builds the primary two-entry test fixture used throughout this file. Identical to the
	/// fixture from TestMslWriter so that writer and reader tests share the same ground truth.
	///
	/// Entry 0: PEPTIDE/2, protein P12345, two fragments (b2 and y1).
	/// Entry 1: PEPTIDE/3, no protein, one fragment (y1).
	/// Both share stripped sequence "PEPTIDE" → same ElutionGroupId after write.
	/// </summary>
	/// <returns>A list of two fully populated <see cref="MslLibraryEntry"/> objects.</returns>
	private static List<MslLibraryEntry> BuildTestEntries()
	{
		// ── Entry 0: PEPTIDE/2 ────────────────────────────────────────────────
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
				new()
				{
					Mz                      = 175.119f,
					Intensity               = 8000f,
					ProductType             = ProductType.y,
					SecondaryProductType    = null,
					FragmentNumber          = 1,
					SecondaryFragmentNumber = 0,
					ResiduePosition         = 6,
					Charge                  = 1,
					NeutralLoss             = 0.0,
					ExcludeFromQuant        = false
				},
				new()
				{
					Mz                      = 98.060f,
					Intensity               = 10000f,   // highest → 1.0 after normalize
                    ProductType             = ProductType.b,
					SecondaryProductType    = null,
					FragmentNumber          = 2,
					SecondaryFragmentNumber = 0,
					ResiduePosition         = 1,
					Charge                  = 1,
					NeutralLoss             = 0.0,
					ExcludeFromQuant        = false
				}
			}
		};

		// ── Entry 1: PEPTIDE/3 — same stripped seq, different charge ──────────
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
				new()
				{
					Mz                      = 175.119f,
					Intensity               = 5000f,
					ProductType             = ProductType.y,
					SecondaryProductType    = null,
					FragmentNumber          = 1,
					SecondaryFragmentNumber = 0,
					ResiduePosition         = 6,
					Charge                  = 1,
					NeutralLoss             = 0.0,
					ExcludeFromQuant        = false
				}
			}
		};

		return new List<MslLibraryEntry> { entry0, entry1 };
	}

	/// <summary>
	/// Builds a single-entry fixture containing one internal fragment ion (bIy type).
	/// Used by internal-ion round-trip tests.
	/// </summary>
	/// <returns>A list containing exactly one <see cref="MslLibraryEntry"/>.</returns>
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
			IsProteotypic = false,
			QValue = float.NaN,
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new()
				{
					Mz                      = 364.197f,
					Intensity               = 9000f,
					ProductType             = ProductType.y,
					SecondaryProductType    = null,
					FragmentNumber          = 3,
					SecondaryFragmentNumber = 0,
					ResiduePosition         = 6,
					Charge                  = 1,
					NeutralLoss             = 0.0,
					ExcludeFromQuant        = false
				},
				new()
				{
					Mz                      = 485.215f,
					Intensity               = 10000f,
					ProductType             = ProductType.b,
					SecondaryProductType    = ProductType.y,  // internal ion: bIy
                    FragmentNumber          = 2,              // start residue (0-based)
                    SecondaryFragmentNumber = 5,              // end residue (0-based)
                    ResiduePosition         = 2,
					Charge                  = 1,
					NeutralLoss             = 0.0,
					ExcludeFromQuant        = false
				}
			}
		};

		return new List<MslLibraryEntry> { entry };
	}

	// ── Group: Validation tests ───────────────────────────────────────────────

	/// <summary>
	/// Verifies that <see cref="MslReader.Load"/> throws <see cref="FileNotFoundException"/>
	/// when the requested path does not exist on disk.
	/// </summary>
	[Test]
	public void Load_ThrowsFileNotFoundException_WhenFileDoesNotExist()
	{
		string nonExistent = TempPath("_does_not_exist_");

		// Remove any leftover from a prior run
		if (File.Exists(nonExistent)) File.Delete(nonExistent);

		Assert.That(() => MslReader.Load(nonExistent),
			Throws.TypeOf<FileNotFoundException>(),
			"Load must throw FileNotFoundException for a path that does not exist.");
	}

	/// <summary>
	/// Verifies that <see cref="MslReader.Load"/> throws <see cref="FormatException"/>
	/// when the first four bytes of the file do not match the MSL magic.
	/// </summary>
	[Test]
	public void Load_ThrowsFormatException_WhenMagicIsWrong()
	{
		string path = TempPath(nameof(Load_ThrowsFormatException_WhenMagicIsWrong));

		// Write a valid file, then corrupt the first four bytes
		MslWriter.Write(path, BuildTestEntries());
		byte[] data = File.ReadAllBytes(path);
		data[0] = 0xFF; data[1] = 0xFF; data[2] = 0xFF; data[3] = 0xFF;
		File.WriteAllBytes(path, data);

		Assert.That(() => MslReader.Load(path),
			Throws.TypeOf<FormatException>().With.Message.Contains("Magic"),
			"Load must throw FormatException when magic bytes are wrong.");
	}

	/// <summary>
	/// Verifies that <see cref="MslReader.Load"/> throws <see cref="FormatException"/>
	/// when the format version field in the header does not equal
	/// <see cref="MslFormat.CurrentVersion"/>.
	/// </summary>
	[Test]
	public void Load_ThrowsFormatException_WhenVersionIsNotOne()
	{
		string path = TempPath(nameof(Load_ThrowsFormatException_WhenVersionIsNotOne));

		MslWriter.Write(path, BuildTestEntries());
		byte[] data = File.ReadAllBytes(path);

		// FormatVersion is at bytes 4–7 (little-endian int32); write version 99
		data[4] = 99; data[5] = 0; data[6] = 0; data[7] = 0;
		File.WriteAllBytes(path, data);

		Assert.That(() => MslReader.Load(path),
			Throws.TypeOf<FormatException>().With.Message.Contains("version").Or.With.Message.Contains("Version"),
			"Load must throw FormatException when format version is unsupported.");
	}

	/// <summary>
	/// Verifies that <see cref="MslReader.Load"/> throws <see cref="FormatException"/>
	/// when the trailing magic in the footer does not match, indicating the file is
	/// truncated or otherwise corrupt.
	/// </summary>
	[Test]
	public void Load_ThrowsFormatException_WhenFooterMagicIsWrong()
	{
		string path = TempPath(nameof(Load_ThrowsFormatException_WhenFooterMagicIsWrong));

		MslWriter.Write(path, BuildTestEntries());
		byte[] data = File.ReadAllBytes(path);

		// Corrupt the last 4 bytes (trailing magic)
		data[^1] = 0xDE; data[^2] = 0xAD; data[^3] = 0xBE; data[^4] = 0xEF;
		File.WriteAllBytes(path, data);

		Assert.That(() => MslReader.Load(path),
			Throws.TypeOf<FormatException>().With.Message.Contains("magic").Or.With.Message.Contains("Magic").Or.With.Message.Contains("trailing").Or.With.Message.Contains("Trailing"),
			"Load must throw FormatException when trailing footer magic is corrupt.");
	}

	/// <summary>
	/// Verifies that <see cref="MslReader.Load"/> throws <see cref="FormatException"/>
	/// when the NPrecursors count stored in the footer differs from the one in the header.
	/// </summary>
	[Test]
	public void Load_ThrowsFormatException_WhenNPrecursorsMismatch()
	{
		string path = TempPath(nameof(Load_ThrowsFormatException_WhenNPrecursorsMismatch));

		MslWriter.Write(path, BuildTestEntries());
		byte[] data = File.ReadAllBytes(path);

		// Footer NPrecursors is at offset (fileLength - FooterSize + 8): OffsetTableOffset=8B, NPrecursors=4B
		// Layout: [OffsetTableOffset:8][NPrecursors:4][DataCrc32:4][TrailingMagic:4]
		int footerNPrecursorsOffset = data.Length - MslFormat.FooterSize + 8;

		// Set footer NPrecursors to a value that cannot match the header (100)
		data[footerNPrecursorsOffset] = 100;
		data[footerNPrecursorsOffset + 1] = 0;
		data[footerNPrecursorsOffset + 2] = 0;
		data[footerNPrecursorsOffset + 3] = 0;
		File.WriteAllBytes(path, data);

		Assert.That(() => MslReader.Load(path),
			Throws.TypeOf<FormatException>().With.Message.Contains("NPrecursors").Or.With.Message.Contains("mismatch").Or.With.Message.Contains("Mismatch"),
			"Load must throw FormatException when footer NPrecursors doesn't match header.");
	}

	/// <summary>
	/// Verifies that <see cref="MslReader.Load"/> throws <see cref="InvalidDataException"/>
	/// when a data byte is corrupted, causing the CRC-32 check to fail.
	/// </summary>
	[Test]
	public void Load_ThrowsInvalidDataException_WhenCRC32IsWrong()
	{
		string path = TempPath(nameof(Load_ThrowsInvalidDataException_WhenCRC32IsWrong));

		MslWriter.Write(path, BuildTestEntries());
		byte[] data = File.ReadAllBytes(path);

		// Corrupt a byte in the middle of the header section (byte 32 is inside the header)
		// but do NOT corrupt the footer itself (which is excluded from CRC coverage)
		data[32] ^= 0xFF;
		File.WriteAllBytes(path, data);

		Assert.That(() => MslReader.Load(path),
			Throws.TypeOf<InvalidDataException>().With.Message.Contains("CRC"),
			"Load must throw InvalidDataException when data bytes fail the CRC-32 check.");
	}

	// ── Group: ReadHeaderOnly ─────────────────────────────────────────────────

	/// <summary>
	/// Verifies that <see cref="MslReader.ReadHeaderOnly"/> returns the correct NPrecursors
	/// count without loading entries or fragments.
	/// </summary>
	[Test]
	public void ReadHeaderOnly_ReturnsCorrectNPrecursors()
	{
		string path = TempPath(nameof(ReadHeaderOnly_ReturnsCorrectNPrecursors));
		List<MslLibraryEntry> entries = BuildTestEntries();
		MslWriter.Write(path, entries);

		MslFileHeader header = MslReader.ReadHeaderOnly(path);

		Assert.That(header.NPrecursors, Is.EqualTo(entries.Count),
			"ReadHeaderOnly must return the correct precursor count from the file header.");
	}

	/// <summary>
	/// Verifies that <see cref="MslReader.ReadHeaderOnly"/> reports format version 1.
	/// </summary>
	[Test]
	public void ReadHeaderOnly_ReturnsCorrectVersion()
	{
		string path = TempPath(nameof(ReadHeaderOnly_ReturnsCorrectVersion));
		MslWriter.Write(path, BuildTestEntries());

		MslFileHeader header = MslReader.ReadHeaderOnly(path);

		Assert.That(header.FormatVersion, Is.EqualTo(MslFormat.CurrentVersion),
			"ReadHeaderOnly must return FormatVersion == MslFormat.CurrentVersion.");
	}

	// ── Group: Full load — precursor fields ────────────────────────────────────

	/// <summary>
	/// Verifies that the number of entries returned by <see cref="MslReader.Load"/>
	/// equals the number written.
	/// </summary>
	[Test]
	public void Load_FullLoad_PrecursorCount_MatchesWrittenCount()
	{
		string path = TempPath(nameof(Load_FullLoad_PrecursorCount_MatchesWrittenCount));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Count, Is.EqualTo(written.Count),
			"Full load must return exactly as many entries as were written.");
	}

	/// <summary>
	/// Verifies that precursor m/z values survive write → read to within float32 precision.
	/// </summary>
	[Test]
	public void Load_FullLoad_PrecursorMz_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_PrecursorMz_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		// PrecursorMz is stored as float32; compare to float32 precision
		Assert.That((float)lib.Entries[0].PrecursorMz,
			Is.EqualTo((float)written[0].PrecursorMz).Within(1e-4f),
			"PrecursorMz must round-trip to float32 precision.");

		Assert.That((float)lib.Entries[1].PrecursorMz,
			Is.EqualTo((float)written[1].PrecursorMz).Within(1e-4f),
			"PrecursorMz for entry 1 must round-trip.");
	}

	/// <summary>
	/// Verifies that precursor charge states survive write → read exactly.
	/// </summary>
	[Test]
	public void Load_FullLoad_Charge_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_Charge_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].ChargeState, Is.EqualTo(written[0].ChargeState));
		Assert.That(lib.Entries[1].ChargeState, Is.EqualTo(written[1].ChargeState));
	}

	/// <summary>
	/// Verifies that iRT values survive write → read to float32 precision.
	/// </summary>
	[Test]
	public void Load_FullLoad_Irt_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_Irt_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That((float)lib.Entries[0].RetentionTime,
			Is.EqualTo((float)written[0].RetentionTime).Within(1e-4f),
			"RetentionTime must round-trip to float32 precision.");
	}

	/// <summary>
	/// Verifies that ion mobility values survive write → read to float32 precision.
	/// </summary>
	[Test]
	public void Load_FullLoad_IonMobility_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_IonMobility_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That((float)lib.Entries[0].IonMobility,
			Is.EqualTo((float)written[0].IonMobility).Within(1e-6f),
			"IonMobility must round-trip.");
	}

	/// <summary>
	/// Verifies that modified sequences survive write → read as exact string matches.
	/// </summary>
	[Test]
	public void Load_FullLoad_ModifiedSequence_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_ModifiedSequence_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].FullSequence, Is.EqualTo(written[0].FullSequence));
		Assert.That(lib.Entries[1].FullSequence, Is.EqualTo(written[1].FullSequence));
	}

	/// <summary>
	/// Verifies that stripped sequences survive write → read as exact string matches.
	/// </summary>
	[Test]
	public void Load_FullLoad_StrippedSequence_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_StrippedSequence_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].BaseSequence, Is.EqualTo(written[0].BaseSequence));
	}

	/// <summary>
	/// Verifies that the IsDecoy flag survives write → read.
	/// </summary>
	[Test]
	public void Load_FullLoad_IsDecoy_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_IsDecoy_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].IsDecoy, Is.EqualTo(written[0].IsDecoy));
	}

	/// <summary>
	/// Verifies that entries with the same stripped sequence receive the same ElutionGroupId
	/// after write → read.
	/// </summary>
	[Test]
	public void Load_FullLoad_ElutionGroupId_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_ElutionGroupId_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		// Both entries share "PEPTIDE" → same elution group ID
		Assert.That(lib.Entries[0].ElutionGroupId, Is.EqualTo(lib.Entries[1].ElutionGroupId),
			"Entries sharing the same stripped sequence must have the same ElutionGroupId.");
	}

	/// <summary>
	/// Verifies that the molecule type enum value survives write → read.
	/// </summary>
	[Test]
	public void Load_FullLoad_MoleculeType_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_MoleculeType_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].MoleculeType, Is.EqualTo(written[0].MoleculeType));
	}

	/// <summary>
	/// Verifies that the dissociation type enum value survives write → read.
	/// </summary>
	[Test]
	public void Load_FullLoad_DissociationType_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_DissociationType_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].DissociationType, Is.EqualTo(written[0].DissociationType));
	}

	// ── Group: Full load — fragment fields ────────────────────────────────────

	/// <summary>
	/// Verifies that the fragment count per precursor survives write → read.
	/// </summary>
	[Test]
	public void Load_FullLoad_FragmentCount_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_FragmentCount_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		// Writer sorts + normalizes in-place; fragment counts are unchanged
		Assert.That(lib.Entries[0].MatchedFragmentIons.Count, Is.EqualTo(written[0].MatchedFragmentIons.Count),
			"Entry 0 fragment count must match.");
		Assert.That(lib.Entries[1].MatchedFragmentIons.Count, Is.EqualTo(written[1].MatchedFragmentIons.Count),
			"Entry 1 fragment count must match.");
	}

	/// <summary>
	/// Verifies that fragment m/z values survive write → read to float32 precision.
	/// </summary>
	[Test]
	public void Load_FullLoad_FragmentMz_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_FragmentMz_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		// MatchedFragmentIons are sorted by m/z ascending; entry 0 has b2 first, then y1
		// After write: [b2@98.060, y1@175.119]
		Assert.That(lib.Entries[0].MatchedFragmentIons[0].Mz,
			Is.EqualTo(98.060f).Within(1e-3f),
			"b2 fragment m/z (lowest) must be first after m/z sort.");

		Assert.That(lib.Entries[0].MatchedFragmentIons[1].Mz,
			Is.EqualTo(175.119f).Within(1e-3f),
			"y1 fragment m/z must survive round-trip.");
	}

	/// <summary>
	/// Verifies that fragment intensities are normalized to [0, 1] and survive write → read.
	/// The highest-intensity fragment in each precursor must have intensity 1.0 on disk.
	/// </summary>
	[Test]
	public void Load_FullLoad_FragmentIntensity_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_FragmentIntensity_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		// Max intensity fragment for entry 0 is b2 (10000 raw → 1.0 normalized), m/z=98.060
		float maxIntensity = lib.Entries[0].MatchedFragmentIons.Max(f => f.Intensity);
		Assert.That(maxIntensity, Is.EqualTo(1.0f).Within(1e-5f),
			"The most-abundant fragment must have normalized intensity 1.0.");

		// y1 was 8000 raw out of 10000 max → 0.8 after normalization
		float y1Intensity = lib.Entries[0].MatchedFragmentIons.First(f => f.ProductType == ProductType.y).Intensity;
		Assert.That(y1Intensity, Is.EqualTo(0.8f).Within(1e-4f),
			"y1 intensity must be 0.8 after normalization (8000/10000).");
	}

	/// <summary>
	/// Verifies that fragment product types survive write → read exactly.
	/// </summary>
	[Test]
	public void Load_FullLoad_FragmentProductType_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_FragmentProductType_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		// Sorted ascending: [b2, y1]
		Assert.That(lib.Entries[0].MatchedFragmentIons[0].ProductType, Is.EqualTo(ProductType.b));
		Assert.That(lib.Entries[0].MatchedFragmentIons[1].ProductType, Is.EqualTo(ProductType.y));
	}

	/// <summary>
	/// Verifies that fragment numbers (ion series indices) survive write → read exactly.
	/// </summary>
	[Test]
	public void Load_FullLoad_FragmentNumber_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_FragmentNumber_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		// Sorted: [b2@98.060 → FragmentNumber=2, y1@175.119 → FragmentNumber=1]
		Assert.That(lib.Entries[0].MatchedFragmentIons[0].FragmentNumber, Is.EqualTo(2));
		Assert.That(lib.Entries[0].MatchedFragmentIons[1].FragmentNumber, Is.EqualTo(1));
	}

	/// <summary>
	/// Verifies that fragment charge states survive write → read.
	/// </summary>
	[Test]
	public void Load_FullLoad_FragmentCharge_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_FragmentCharge_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].MatchedFragmentIons[0].Charge, Is.EqualTo(1));
		Assert.That(lib.Entries[0].MatchedFragmentIons[1].Charge, Is.EqualTo(1));
	}

	/// <summary>
	/// Verifies that the SecondaryProductType of an internal fragment ion is correctly
	/// reconstructed after write → read (must be non-null).
	/// </summary>
	[Test]
	public void Load_FullLoad_InternalIon_SecondaryProductType_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_InternalIon_SecondaryProductType_RoundTrips));
		List<MslLibraryEntry> written = BuildInternalIonEntry();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		// The bIy ion is at higher m/z (485.215), so it should be at index 1 after sort
		MslFragmentIon internalIon = lib.Entries[0].MatchedFragmentIons
			.FirstOrDefault(f => f.SecondaryProductType != null)!;

		Assert.That(internalIon, Is.Not.Null,
			"At least one fragment with SecondaryProductType != null must exist.");
		Assert.That(internalIon.SecondaryProductType, Is.EqualTo(ProductType.y),
			"Internal ion SecondaryProductType must be y.");
		Assert.That(internalIon.ProductType, Is.EqualTo(ProductType.b),
			"Internal ion primary ProductType must be b.");
	}

	/// <summary>
	/// Verifies that the start residue (FragmentNumber) of an internal fragment ion
	/// survives write → read.
	/// </summary>
	[Test]
	public void Load_FullLoad_InternalIon_StartResidue_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_InternalIon_StartResidue_RoundTrips));
		List<MslLibraryEntry> written = BuildInternalIonEntry();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		MslFragmentIon internalIon = lib.Entries[0].MatchedFragmentIons
			.First(f => f.SecondaryProductType != null);

		Assert.That(internalIon.FragmentNumber, Is.EqualTo(2),
			"Internal ion start residue (FragmentNumber) must be 2.");
	}

	/// <summary>
	/// Verifies that the end residue (SecondaryFragmentNumber) of an internal fragment ion
	/// survives write → read.
	/// </summary>
	[Test]
	public void Load_FullLoad_InternalIon_EndResidue_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_InternalIon_EndResidue_RoundTrips));
		List<MslLibraryEntry> written = BuildInternalIonEntry();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		MslFragmentIon internalIon = lib.Entries[0].MatchedFragmentIons
			.First(f => f.SecondaryProductType != null);

		Assert.That(internalIon.SecondaryFragmentNumber, Is.EqualTo(5),
			"Internal ion end residue (SecondaryFragmentNumber) must be 5.");
	}

	/// <summary>
	/// Verifies that neutral-loss masses for named loss codes (H2O) round-trip correctly.
	/// Creates a custom entry with a -H2O fragment, writes, reads, and checks the loss mass.
	/// </summary>
	[Test]
	public void Load_FullLoad_NeutralLoss_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_NeutralLoss_RoundTrips));

		// Build a single entry with a fragment carrying a -H2O neutral loss
		var entries = new List<MslLibraryEntry>
		{
			new()
			{
				FullSequence = "PEPTIDE",
				BaseSequence = "PEPTIDE",
				PrecursorMz      = 449.74,
				ChargeState           = 2,
				RetentionTime              = 30.0,
				Source           = MslFormat.SourceType.Predicted,
				MoleculeType     = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				MatchedFragmentIons        = new List<MslFragmentIon>
				{
					new()
					{
						Mz            = 157.108f, // y1 - H2O approx
                        Intensity     = 1.0f,
						ProductType   = ProductType.y,
						FragmentNumber = 1,
						Charge        = 1,
						NeutralLoss   = -18.010565  // H2O loss
                    }
				}
			}
		};

		MslWriter.Write(path, entries);

		using MslLibraryData lib = MslReader.Load(path);

		double readLoss = lib.Entries[0].MatchedFragmentIons[0].NeutralLoss;
		Assert.That(readLoss, Is.EqualTo(-18.010565).Within(0.001),
			"H2O neutral-loss mass must round-trip from NeutralLossCode.H2O.");
	}

	/// <summary>
	/// Verifies that fragments within each precursor are in m/z ascending order after
	/// full load, confirming the writer's sort is preserved on read.
	/// </summary>
	[Test]
	public void Load_FullLoad_Fragments_AreSortedByMz()
	{
		string path = TempPath(nameof(Load_FullLoad_Fragments_AreSortedByMz));
		MslWriter.Write(path, BuildTestEntries());

		using MslLibraryData lib = MslReader.Load(path);

		foreach (MslLibraryEntry entry in lib.Entries)
		{
			var mzValues = entry.MatchedFragmentIons.Select(f => f.Mz).ToList();
			var sorted = mzValues.OrderBy(mz => mz).ToList();

			Assert.That(mzValues, Is.EqualTo(sorted),
				$"MatchedFragmentIons for entry '{entry.FullSequence}/{entry.ChargeState}' must be m/z ascending.");
		}
	}

	// ── Group: Full load — protein data ──────────────────────────────────────

	/// <summary>
	/// Verifies that the protein accession string survives write → read.
	/// </summary>
	[Test]
	public void Load_FullLoad_ProteinAccession_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_ProteinAccession_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].ProteinAccession, Is.EqualTo("P12345"),
			"ProteinAccession must survive write → read.");
	}

	/// <summary>
	/// Verifies that the gene name string survives write → read.
	/// </summary>
	[Test]
	public void Load_FullLoad_GeneName_RoundTrips()
	{
		string path = TempPath(nameof(Load_FullLoad_GeneName_RoundTrips));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].GeneName, Is.EqualTo("TESTA"),
			"GeneName must survive write → read.");
	}

	/// <summary>
	/// Verifies that entries with no protein assignment return an empty string for
	/// ProteinAccession (not null and not some garbage index-zero value).
	/// </summary>
	[Test]
	public void Load_FullLoad_NoProtein_ProteinAccessionIsEmpty()
	{
		string path = TempPath(nameof(Load_FullLoad_NoProtein_ProteinAccessionIsEmpty));
		MslWriter.Write(path, BuildTestEntries());

		using MslLibraryData lib = MslReader.Load(path);

		// Entry 1 has no protein in the fixture
		Assert.That(lib.Entries[1].ProteinAccession, Is.EqualTo(string.Empty),
			"An entry with no protein must have an empty ProteinAccession after load.");
	}

	// ── Group: Index-only load ────────────────────────────────────────────────

	/// <summary>
	/// Verifies that <see cref="MslReader.LoadIndexOnly"/> returns the correct precursor
	/// count without loading any fragments.
	/// </summary>
	[Test]
	public void LoadIndexOnly_PrecursorCount_MatchesWrittenCount()
	{
		string path = TempPath(nameof(LoadIndexOnly_PrecursorCount_MatchesWrittenCount));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.LoadIndexOnly(path);

		Assert.That(lib.Count, Is.EqualTo(written.Count),
			"LoadIndexOnly must report the correct precursor count.");
		Assert.That(lib.IsIndexOnly, Is.True,
			"IsIndexOnly must be true for a library opened in index-only mode.");
	}

	/// <summary>
	/// Verifies that precursor m/z values are accessible in index-only mode (they are
	/// loaded with the precursor records, not deferred to on-demand reads).
	/// </summary>
	[Test]
	public void LoadIndexOnly_PrecursorMz_IsAvailableWithoutFragments()
	{
		string path = TempPath(nameof(LoadIndexOnly_PrecursorMz_IsAvailableWithoutFragments));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.LoadIndexOnly(path);

		// Fragment lists must be empty in index-only mode
		Assert.That(lib.Entries[0].MatchedFragmentIons, Is.Empty,
			"Fragment list must be empty before on-demand loading.");

		// But PrecursorMz must already be populated
		Assert.That((float)lib.Entries[0].PrecursorMz,
			Is.EqualTo((float)written[0].PrecursorMz).Within(1e-4f),
			"PrecursorMz must be available in index-only mode even though fragments are not loaded.");
	}

	/// <summary>
	/// Verifies that fragments loaded on demand via <see cref="MslLibraryData.LoadFragmentsOnDemand"/>
	/// match the fragments returned by a full load of the same file.
	/// </summary>
	[Test]
	public void LoadIndexOnly_FragmentsLoadedOnDemand_MatchFullLoad()
	{
		string path = TempPath(nameof(LoadIndexOnly_FragmentsLoadedOnDemand_MatchFullLoad));
		MslWriter.Write(path, BuildTestEntries());

		// Full load for reference
		using MslLibraryData fullLib = MslReader.Load(path);
		// Index-only load for on-demand comparison
		using MslLibraryData indexLib = MslReader.LoadIndexOnly(path);

		for (int i = 0; i < fullLib.Count; i++)
		{
			List<MslFragmentIon> onDemand = indexLib.LoadFragmentsOnDemand(i);
			List<MslFragmentIon> preLoaded = fullLib.Entries[i].MatchedFragmentIons;

			Assert.That(onDemand.Count, Is.EqualTo(preLoaded.Count),
				$"Fragment count for entry {i} must match between full-load and on-demand.");

			for (int j = 0; j < preLoaded.Count; j++)
			{
				Assert.That(onDemand[j].Mz,
					Is.EqualTo(preLoaded[j].Mz).Within(1e-5f),
					$"Fragment[{j}].Mz for entry {i} must match.");

				Assert.That(onDemand[j].ProductType,
					Is.EqualTo(preLoaded[j].ProductType),
					$"Fragment[{j}].ProductType for entry {i} must match.");

				Assert.That(onDemand[j].FragmentNumber,
					Is.EqualTo(preLoaded[j].FragmentNumber),
					$"Fragment[{j}].FragmentNumber for entry {i} must match.");
			}
		}
	}

	/// <summary>
	/// Verifies that <see cref="MslLibraryData.Dispose"/> closes the file handle held in
	/// index-only mode. After disposal, subsequent calls to
	/// <see cref="MslLibraryData.LoadFragmentsOnDemand"/> must throw
	/// <see cref="ObjectDisposedException"/>.
	/// </summary>
	[Test]
	public void LoadIndexOnly_DisposesFileStream_OnDispose()
	{
		string path = TempPath(nameof(LoadIndexOnly_DisposesFileStream_OnDispose));
		MslWriter.Write(path, BuildTestEntries());

		MslLibraryData lib = MslReader.LoadIndexOnly(path);
		lib.Dispose();

		Assert.That(() => lib.LoadFragmentsOnDemand(0),
			Throws.TypeOf<ObjectDisposedException>(),
			"LoadFragmentsOnDemand must throw ObjectDisposedException after Dispose().");

		// Also verify the file can now be opened exclusively (stream is truly closed)
		Assert.That(() =>
		{
			using var exclusive = new FileStream(path, FileMode.Open, FileAccess.ReadWrite, FileShare.None);
		}, Throws.Nothing,
		"File must be openable exclusively after the library is disposed.");
	}

	/// <summary>
	/// Verifies that multiple threads can call <see cref="MslLibraryData.LoadFragmentsOnDemand"/>
	/// concurrently without data corruption or exceptions. Each thread reads the same
	/// precursor repeatedly and checks the fragment count.
	/// </summary>
	[Test]
	public void LoadIndexOnly_IsThreadSafe_ConcurrentFragmentReads()
	{
		string path = TempPath(nameof(LoadIndexOnly_IsThreadSafe_ConcurrentFragmentReads));
		MslWriter.Write(path, BuildTestEntries());

		using MslLibraryData lib = MslReader.LoadIndexOnly(path);

		// 8 threads, each reading all entries 25 times = 400 reads total
		const int ThreadCount = 8;
		const int ReadsPerThread = 25;

		var errors = new System.Collections.Concurrent.ConcurrentBag<string>();

		Parallel.For(0, ThreadCount, _ =>
		{
			for (int iter = 0; iter < ReadsPerThread; iter++)
			{
				for (int i = 0; i < lib.Count; i++)
				{
					try
					{
						List<MslFragmentIon> frags = lib.LoadFragmentsOnDemand(i);

						// We know the expected count for each entry in the two-entry fixture
						// (writer sorts in-place; entry 0 has 2 fragments, entry 1 has 1)
						int expectedCount = i == 0 ? 2 : 1;
						if (frags.Count != expectedCount)
						{
							errors.Add(
								$"Thread: entry {i} had {frags.Count} fragments, expected {expectedCount}.");
						}
					}
					catch (Exception ex)
					{
						errors.Add($"Exception reading entry {i}: {ex.Message}");
					}
				}
			}
		});

		Assert.That(errors, Is.Empty,
			"Concurrent LoadFragmentsOnDemand calls must not produce errors:\n" +
			string.Join("\n", errors));
	}

	// ── Group: Round-trip fidelity ────────────────────────────────────────────

	/// <summary>
	/// Verifies that all precursor-level fields for all fixture entries are identical
	/// after a write → read cycle (to the precision of their stored types).
	/// </summary>
	[Test]
	public void WriteRead_RoundTrip_AllPrecursorFields_Identical()
	{
		string path = TempPath(nameof(WriteRead_RoundTrip_AllPrecursorFields_Identical));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		for (int i = 0; i < written.Count; i++)
		{
			MslLibraryEntry w = written[i];
			MslLibraryEntry r = lib.Entries[i];

			Assert.That((float)r.PrecursorMz, Is.EqualTo((float)w.PrecursorMz).Within(1e-4f), $"[{i}] PrecursorMz");
			Assert.That(r.ChargeState, Is.EqualTo(w.ChargeState), $"[{i}] ChargeState");
			Assert.That((float)r.RetentionTime, Is.EqualTo((float)w.RetentionTime).Within(1e-4f), $"[{i}] RetentionTime");
			Assert.That(r.FullSequence, Is.EqualTo(w.FullSequence), $"[{i}] FullSequence");
			Assert.That(r.BaseSequence, Is.EqualTo(w.BaseSequence), $"[{i}] BaseSequence");
			Assert.That(r.IsDecoy, Is.EqualTo(w.IsDecoy), $"[{i}] IsDecoy");
			Assert.That(r.MoleculeType, Is.EqualTo(w.MoleculeType), $"[{i}] MoleculeType");
			Assert.That(r.DissociationType, Is.EqualTo(w.DissociationType), $"[{i}] DissociationType");
			Assert.That(r.Source, Is.EqualTo(w.Source), $"[{i}] Source");
		}
	}

	/// <summary>
	/// Verifies that all fragment-level fields for all fixture entries are identical
	/// after a write → read cycle.
	/// </summary>
	[Test]
	public void WriteRead_RoundTrip_AllFragmentFields_Identical()
	{
		string path = TempPath(nameof(WriteRead_RoundTrip_AllFragmentFields_Identical));
		List<MslLibraryEntry> written = BuildTestEntries();
		MslWriter.Write(path, written);

		using MslLibraryData lib = MslReader.Load(path);

		for (int i = 0; i < written.Count; i++)
		{
			// After write the writer has sorted and normalized in-place;
			// compare against the (now mutated) written lists
			List<MslFragmentIon> wFrags = written[i].MatchedFragmentIons;
			List<MslFragmentIon> rFrags = lib.Entries[i].MatchedFragmentIons;

			Assert.That(rFrags.Count, Is.EqualTo(wFrags.Count), $"[{i}] fragment count");

			for (int j = 0; j < wFrags.Count; j++)
			{
				Assert.That(rFrags[j].Mz, Is.EqualTo(wFrags[j].Mz).Within(1e-4f), $"[{i}][{j}] Mz");
				Assert.That(rFrags[j].Intensity, Is.EqualTo(wFrags[j].Intensity).Within(1e-5f), $"[{i}][{j}] Intensity");
				Assert.That(rFrags[j].ProductType, Is.EqualTo(wFrags[j].ProductType), $"[{i}][{j}] ProductType");
				Assert.That(rFrags[j].FragmentNumber, Is.EqualTo(wFrags[j].FragmentNumber), $"[{i}][{j}] FragmentNumber");
				Assert.That(rFrags[j].Charge, Is.EqualTo(wFrags[j].Charge), $"[{i}][{j}] ChargeState");
			}
		}
	}

	/// <summary>
	/// Verifies that calling <see cref="MslLibraryEntry.ToLibrarySpectrum"/> on a read entry
	/// produces a <see cref="LibrarySpectrum"/> whose sequence matches the written
	/// FullSequence (the "Name" in the DDA dictionary).
	/// </summary>
	[Test]
	public void WriteRead_RoundTrip_ToLibrarySpectrum_NameMatchesExpected()
	{
		string path = TempPath(nameof(WriteRead_RoundTrip_ToLibrarySpectrum_NameMatchesExpected));
		MslWriter.Write(path, BuildTestEntries());

		using MslLibraryData lib = MslReader.Load(path);

		LibrarySpectrum spectrum = lib.Entries[0].ToLibrarySpectrum();

		Assert.That(spectrum.Sequence, Is.EqualTo("PEPTIDE"),
			"ToLibrarySpectrum().Sequence must match the original FullSequence.");
	}

	/// <summary>
	/// Verifies that fragment ion annotations produced by <see cref="MslLibraryEntry.ToLibrarySpectrum"/>
	/// are correctly typed (b or y) after a write → read cycle.
	/// </summary>
	[Test]
	public void WriteRead_RoundTrip_ToLibrarySpectrum_FragmentIonAnnotationsCorrect()
	{
		string path = TempPath(nameof(WriteRead_RoundTrip_ToLibrarySpectrum_FragmentIonAnnotationsCorrect));
		MslWriter.Write(path, BuildTestEntries());

		using MslLibraryData lib = MslReader.Load(path);

		LibrarySpectrum spectrum = lib.Entries[0].ToLibrarySpectrum();

		// Entry 0 has [b2, y1] after sort; both should survive the ToLibrarySpectrum conversion
		var productTypes = spectrum.MatchedFragmentIons
			.Select(ion => ion.NeutralTheoreticalProduct.ProductType)
			.ToHashSet();

		Assert.That(productTypes, Contains.Item(ProductType.b),
			"ToLibrarySpectrum must include b-type fragments.");
		Assert.That(productTypes, Contains.Item(ProductType.y),
			"ToLibrarySpectrum must include y-type fragments.");
	}

	// ── Group: Edge cases ─────────────────────────────────────────────────────

	/// <summary>
	/// Verifies that a file written with zero precursors loads successfully and returns
	/// an empty Entries list.
	/// </summary>
	[Test]
	public void Load_EmptyLibrary_ZeroPrecursors_LoadsSuccessfully()
	{
		string path = TempPath(nameof(Load_EmptyLibrary_ZeroPrecursors_LoadsSuccessfully));
		MslWriter.Write(path, new List<MslLibraryEntry>());

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Count, Is.EqualTo(0),
			"An empty library must load with zero entries.");
	}

	/// <summary>
	/// Verifies that a precursor with no protein data (ProteinAccession = "") loads
	/// without error and returns empty protein strings.
	/// </summary>
	[Test]
	public void Load_SinglePrecursor_NoProtein_LoadsSuccessfully()
	{
		string path = TempPath(nameof(Load_SinglePrecursor_NoProtein_LoadsSuccessfully));

		var entries = new List<MslLibraryEntry>
		{
			new()
			{
				FullSequence = "PEPTIDE",
				BaseSequence = "PEPTIDE",
				PrecursorMz      = 449.74,
				ChargeState           = 2,
				Source           = MslFormat.SourceType.Predicted,
				MoleculeType     = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = string.Empty,
				ProteinName      = string.Empty,
				GeneName         = string.Empty,
				MatchedFragmentIons        = new List<MslFragmentIon>
				{
					new()
					{
						Mz            = 175.119f,
						Intensity     = 1.0f,
						ProductType   = ProductType.y,
						FragmentNumber = 1,
						Charge        = 1
					}
				}
			}
		};

		MslWriter.Write(path, entries);

		Assert.That(() =>
		{
			using MslLibraryData lib = MslReader.Load(path);
			Assert.That(lib.Count, Is.EqualTo(1));
			Assert.That(lib.Entries[0].ProteinAccession, Is.EqualTo(string.Empty));
		}, Throws.Nothing,
		"A single no-protein precursor must load without errors.");
	}

	/// <summary>
	/// Verifies that a precursor with a large number of fragments (50) loads correctly
	/// with all fragments intact.
	/// </summary>
	[Test]
	public void Load_PrecursorWithMaxFragmentCount_LoadsCorrectly()
	{
		string path = TempPath(nameof(Load_PrecursorWithMaxFragmentCount_LoadsCorrectly));

		// Build 50 y-ion fragments with distinct m/z values
		const int FragmentCount = 50;
		var fragments = new List<MslFragmentIon>(FragmentCount);
		for (int i = 1; i <= FragmentCount; i++)
		{
			fragments.Add(new MslFragmentIon
			{
				Mz = 100.0f + i * 10.0f,
				Intensity = 1.0f / i,  // decreasing intensity to avoid all-equal
				ProductType = ProductType.y,
				FragmentNumber = i,
				Charge = 1
			});
		}

		var entries = new List<MslLibraryEntry>
		{
			new()
			{
				FullSequence = "ACDEFGHIKLMNPQRSTVWY",
				BaseSequence = "ACDEFGHIKLMNPQRSTVWY",
				PrecursorMz      = 1000.0,
				ChargeState           = 2,
				Source           = MslFormat.SourceType.Predicted,
				MoleculeType     = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				MatchedFragmentIons        = fragments
			}
		};

		MslWriter.Write(path, entries);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].MatchedFragmentIons.Count, Is.EqualTo(FragmentCount),
			$"All {FragmentCount} fragments must survive write → read.");
	}

	/// <summary>
	/// Verifies that an oligonucleotide precursor with
	/// <see cref="MslFormat.MoleculeType.Oligonucleotide"/> round-trips correctly, and that
	/// calling <see cref="MslLibraryEntry.ToLibrarySpectrum"/> on the loaded entry does not
	/// throw. This exercises the MoleculeType-branching terminus-reconstruction code path
	/// that is distinct from the peptide path.
	/// </summary>
	[Test]
	public void Load_OligonucleotidePrecursor_FivePrimeTerminus_ReconstructedCorrectly()
	{
		string path = TempPath(nameof(Load_OligonucleotidePrecursor_FivePrimeTerminus_ReconstructedCorrectly));

		// Build a minimal oligonucleotide entry with a w-type fragment (5'-terminus series).
		// ProductType.w is the 5'-terminus ion class used in the Omics.Fragmentation.Oligo namespace.
		var entries = new List<MslLibraryEntry>
		{
			new()
			{
				FullSequence = "ACGU",
				BaseSequence = "ACGU",
				PrecursorMz      = 600.0,
				ChargeState           = 2,
				Source           = MslFormat.SourceType.Predicted,
				MoleculeType     = MslFormat.MoleculeType.Oligonucleotide,
				DissociationType = DissociationType.HCD,
				MatchedFragmentIons        = new List<MslFragmentIon>
				{
					new()
					{
						Mz             = 309.0f,  // w1 of ACGU (approximate)
                        Intensity      = 1.0f,
						ProductType    = ProductType.w,
						FragmentNumber = 1,
						Charge         = 1
					}
				}
			}
		};

		MslWriter.Write(path, entries);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Count, Is.EqualTo(1),
			"Oligonucleotide entry must load successfully.");
		Assert.That(lib.Entries[0].MoleculeType,
			Is.EqualTo(MslFormat.MoleculeType.Oligonucleotide),
			"MoleculeType must be Oligonucleotide after round-trip.");
		Assert.That(lib.Entries[0].MatchedFragmentIons.Count, Is.EqualTo(1),
			"Fragment count must be 1 after round-trip.");

		// Call ToLibrarySpectrum() to exercise the terminus-reconstruction branch for oligos.
		// This must not throw even though ProductType.w has a different FragmentationTerminus
		// mapping than standard peptide ions.
		Assert.That(() => lib.Entries[0].ToLibrarySpectrum(), Throws.Nothing,
			"ToLibrarySpectrum() must not throw for an oligonucleotide entry.");
	}

	/// <summary>
	/// Verifies that a proteoform (large protein sequence) precursor loads without error.
	/// This exercises the code path that handles large StrippedSeqLength values.
	/// </summary>
	[Test]
	public void Load_ProteoformPrecursor_LargeStrippedSeqLength_LoadsCorrectly()
	{
		string path = TempPath(nameof(Load_ProteoformPrecursor_LargeStrippedSeqLength_LoadsCorrectly));

		// Build a synthetic proteoform with a 200-residue sequence
		string longSeq = new string('A', 200);

		var entries = new List<MslLibraryEntry>
		{
			new()
			{
				FullSequence = longSeq,
				BaseSequence = longSeq,
				PrecursorMz      = 1200.0,
				ChargeState           = 10,
				Source           = MslFormat.SourceType.Predicted,
				MoleculeType     = MslFormat.MoleculeType.Proteoform,
				DissociationType = DissociationType.HCD,
				MatchedFragmentIons        = new List<MslFragmentIon>
				{
					new()
					{
						Mz             = 500.0f,
						Intensity      = 1.0f,
						ProductType    = ProductType.y,
						FragmentNumber = 10,
						Charge         = 2
					}
				}
			}
		};

		MslWriter.Write(path, entries);

		using MslLibraryData lib = MslReader.Load(path);

		Assert.That(lib.Count, Is.EqualTo(1), "Proteoform entry must load.");
		Assert.That(lib.Entries[0].BaseSequence.Length, Is.EqualTo(200),
			"Proteoform stripped sequence length must survive write → read.");
		Assert.That(lib.Entries[0].MoleculeType, Is.EqualTo(MslFormat.MoleculeType.Proteoform),
			"MoleculeType must be Proteoform.");
	}
}