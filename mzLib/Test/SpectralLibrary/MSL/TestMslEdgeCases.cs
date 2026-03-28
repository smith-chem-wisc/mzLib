using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Boundary conditions, corruption detection, and unusual-input tests for the MSL binary
/// spectral library format.
///
/// Coverage areas:
/// <list type="number">
///   <item>Empty and minimal libraries (zero precursors, single precursor, no fragments)</item>
///   <item>Long sequences (proteoforms with hundreds of residues)</item>
///   <item>Unicode in string-table fields (protein names, gene names)</item>
///   <item>File-level corruption detection (magic bytes, truncation, bit flip, version mismatch)</item>
///   <item>Fragmentation edge cases (all neutral-loss codes, diagnostic ions, M/Ycore ions)</item>
///   <item>All ProductType values round-tripping through the binary format</item>
///   <item>Float precision guarantees for m/z, iRT, and q-value storage</item>
/// </list>
///
/// Every test is self-contained: it writes its own temp file in <see cref="SetUp"/> and
/// cleans up in <see cref="TearDown"/>.
/// </summary>
[TestFixture]
public sealed class TestMslEdgeCases
{
	// ── Temp-file management ──────────────────────────────────────────────────

	/// <summary>
	/// Path of the primary .msl temp file for the current test. Set in <see cref="SetUp"/>
	/// and deleted (if present) in <see cref="TearDown"/>.
	/// </summary>
	private string _tempMslPath = string.Empty;

	/// <summary>
	/// Optional second temp path for tests that need two .msl files simultaneously.
	/// Deleted in <see cref="TearDown"/> if non-empty.
	/// </summary>
	private string _tempMslPath2 = string.Empty;

	/// <summary>Allocates unique .msl-extension temp paths before each test.</summary>
	[SetUp]
	public void SetUp()
	{
		string raw = Path.GetTempFileName();
		_tempMslPath = Path.ChangeExtension(raw, ".msl");
		File.Delete(raw);

		string raw2 = Path.GetTempFileName();
		_tempMslPath2 = Path.ChangeExtension(raw2, ".msl");
		File.Delete(raw2);
	}

	/// <summary>Deletes all temp files created for the current test.</summary>
	[TearDown]
	public void TearDown()
	{
		if (File.Exists(_tempMslPath)) File.Delete(_tempMslPath);
		if (File.Exists(_tempMslPath2)) File.Delete(_tempMslPath2);
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Empty and minimal libraries
	// ─────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Writing an empty entry list and then loading it must not throw any exception.
	/// A zero-precursor .msl file is a valid file that can arise during incremental library
	/// construction workflows.
	/// </summary>
	[Test]
	public void EdgeCase_EmptyLibrary_WritesAndLoads_WithoutError()
	{
		// Act / Assert — no exception expected
		Assert.DoesNotThrow(() =>
		{
			MslLibrary.Save(_tempMslPath, new List<MslLibraryEntry>());
			using MslLibrary lib = MslLibrary.Load(_tempMslPath);
			_ = lib.PrecursorCount; // force property access to confirm library is usable
		});
	}

	/// <summary>
	/// After writing and loading an empty library, <see cref="MslLibrary.PrecursorCount"/>
	/// must return exactly zero.
	/// </summary>
	[Test]
	public void EdgeCase_EmptyLibrary_PrecursorCount_IsZero()
	{
		MslLibrary.Save(_tempMslPath, new List<MslLibraryEntry>());
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);

		Assert.That(lib.PrecursorCount, Is.EqualTo(0));
	}

	/// <summary>
	/// A single precursor entry with an empty fragment list (no peaks) must write and load
	/// without error. The loaded entry must have zero fragments.
	/// </summary>
	[Test]
	public void EdgeCase_SinglePrecursor_NoFragments_Loads()
	{
		// Arrange
		var entry = new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			BaseSequence = "PEPTIDE",
			PrecursorMz = 400.21f,
			ChargeState = 2,
			RetentionTime = 35.4,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			QValue = float.NaN,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			MatchedFragmentIons = new List<MslFragmentIon>()
		};

		// Act
		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? loaded);

		// Assert
		Assert.That(found, Is.True);
		Assert.That(loaded!.MatchedFragmentIons.Count, Is.EqualTo(0));
	}

	/// <summary>
	/// A single precursor with 1,000 fragment ions must write and load without error or data
	/// truncation. The test verifies count equality and spot-checks several fragment m/z values
	/// and metadata fields to detect silent data corruption.
	///
	/// 1,000 fragments exercises the large-list code path (well above typical tryptic spectra)
	/// while keeping I/O and allocation suitable for CI.
	/// </summary>
	[Test]
	public void EdgeCase_SinglePrecursor_ManyFragments_WritesAndLoads()
	{
		// Arrange — generate 1 000 fragments with distinct m/z values
		const int FragmentCount = 1_000;
		var frags = new List<MslFragmentIon>(FragmentCount);
		for (int i = 0; i < FragmentCount; i++)
		{
			frags.Add(new MslFragmentIon
			{
				Mz = 100.0f + i * 0.01f,
				Intensity = 1.0f,
				ProductType = ProductType.y,
				FragmentNumber = (i % 10000) + 1,
				ResiduePosition = (i % 10000) + 1,
				Charge = 1,
				NeutralLoss = 0.0
			});
		}

		var entry = new MslLibraryEntry
		{
			FullSequence = "BIGPEPTIDE",
			BaseSequence = "BIGPEPTIDE",
			PrecursorMz = 500.0f,
			ChargeState = 2,
			RetentionTime = 50.0,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			QValue = float.NaN,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			MatchedFragmentIons = frags
		};

		// Act
		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry("BIGPEPTIDE", 2, out MslLibraryEntry? loaded);

		// Assert — count equality
		Assert.That(found, Is.True);
		Assert.That(loaded!.MatchedFragmentIons.Count, Is.EqualTo(FragmentCount));

		// Assert — spot-check content correctness on first, middle, and last fragments (sorted by m/z)
		List<MslFragmentIon> sortedLoaded = loaded.MatchedFragmentIons.OrderBy(f => f.Mz).ToList();

		// First fragment: m/z = 100.0 + 0 * 0.01 = 100.0
		Assert.That((double)sortedLoaded[0].Mz, Is.EqualTo(100.0).Within(1e-4),
			"First fragment m/z mismatch");
		Assert.That(sortedLoaded[0].ProductType, Is.EqualTo(ProductType.y),
			"First fragment ProductType mismatch");

		// Middle fragment: m/z = 100.0 + 500 * 0.01 = 105.0
		Assert.That((double)sortedLoaded[500].Mz, Is.EqualTo(105.0).Within(1e-4),
			"Middle fragment m/z mismatch");
		Assert.That(sortedLoaded[500].FragmentNumber, Is.EqualTo(501),
			"Middle fragment FragmentNumber mismatch");

		// Last fragment: m/z = 100.0 + 999 * 0.01 = 109.99
		Assert.That((double)sortedLoaded[999].Mz, Is.EqualTo(109.99).Within(1e-2),
			"Last fragment m/z mismatch");
		Assert.That(sortedLoaded[999].Charge, Is.EqualTo(1),
			"Last fragment Charge mismatch");
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Long sequences
	// ─────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// A proteoform entry whose stripped sequence is 500 residues long must be written and
	/// loaded without error. The stripped sequence length field in the binary record is int32,
	/// so this exercises values far above typical tryptic peptide lengths (7–30 residues).
	/// </summary>
	[Test]
	public void EdgeCase_LongSequence_500Residues_StrippedSeqLengthFitsInInt32()
	{
		// Arrange — create a 500-residue alanine sequence
		string longSeq = new string('A', 500);

		var entry = new MslLibraryEntry
		{
			FullSequence = longSeq,
			BaseSequence = longSeq,
			PrecursorMz = 900.0f,
			ChargeState = 10,
			RetentionTime = 120.0,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Proteoform,
			DissociationType = DissociationType.ECD,
			QValue = float.NaN,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 200.0f, Intensity = 1.0f, ProductType = ProductType.b,
					FragmentNumber = 5, ResiduePosition = 5, Charge = 1, NeutralLoss = 0.0
				}
			}
		};

		// Act
		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry(longSeq, 10, out MslLibraryEntry? loaded);

		// Assert
		Assert.That(found, Is.True);
		Assert.That(loaded!.BaseSequence.Length, Is.EqualTo(500));
	}

	/// <summary>
	/// A modified sequence containing many bracketed modification tags (e.g. oxidations on
	/// every Met) must survive the string-table round-trip with the notation intact.
	/// </summary>
	[Test]
	public void EdgeCase_LongModifiedSequence_WithManyMods_PreservesNotation()
	{
		// Arrange — 10 Met residues each with oxidation notation
		const string modTag = "[Common Variable:Oxidation on M]";
		string modSeq = string.Concat(Enumerable.Repeat("M" + modTag, 10));
		string stripped = new string('M', 10);

		var entry = new MslLibraryEntry
		{
			FullSequence = modSeq,
			BaseSequence = stripped,
			PrecursorMz = 700.0f,
			ChargeState = 3,
			RetentionTime = 40.0,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			QValue = float.NaN,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 300.0f, Intensity = 1.0f, ProductType = ProductType.y,
					FragmentNumber = 2, ResiduePosition = 2, Charge = 1, NeutralLoss = 0.0
				}
			}
		};

		// Act
		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry(modSeq, 3, out MslLibraryEntry? loaded);

		// Assert
		Assert.That(found, Is.True);
		Assert.That(loaded!.FullSequence, Is.EqualTo(modSeq));
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Unicode in string fields
	// ─────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// A protein name containing Unicode characters (e.g. Greek letters, accented characters)
	/// must round-trip through the UTF-8 string table without mojibake or truncation.
	/// </summary>
	[Test]
	public void EdgeCase_UnicodeProteinName_RoundTrips()
	{
		// Arrange — Unicode protein name with Greek and accented characters
		const string unicodeName = "Hémoprotéine α-chaîne (Ω-oxidase) — résidu β";

		var entry = new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			BaseSequence = "PEPTIDE",
			PrecursorMz = 400.21f,
			ChargeState = 2,
			RetentionTime = 35.4,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			QValue = float.NaN,
			ProteinAccession = "P99999",
			ProteinName = unicodeName,
			GeneName = "HEMO",
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 200.0f, Intensity = 1.0f, ProductType = ProductType.b,
					FragmentNumber = 2, ResiduePosition = 2, Charge = 1, NeutralLoss = 0.0
				}
			}
		};

		// Act
		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? loaded);

		// Assert
		Assert.That(found, Is.True);
		Assert.That(loaded!.ProteinName, Is.EqualTo(unicodeName));
	}

	/// <summary>
	/// An entry with an empty GeneName must be written and loaded without error, and the
	/// loaded GeneName must be an empty string (not null). The string table stores empty
	/// strings at index 0 by convention.
	/// </summary>
	[Test]
	public void EdgeCase_EmptyGeneName_StoredAsIndex0()
	{
		// Arrange
		var entry = new MslLibraryEntry
		{
			FullSequence = "KVFGR",
			BaseSequence = "KVFGR",
			PrecursorMz = 306.18f,
			ChargeState = 2,
			RetentionTime = 28.3,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			QValue = float.NaN,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,  // explicitly empty
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 319.2f, Intensity = 1.0f, ProductType = ProductType.b,
					FragmentNumber = 3, ResiduePosition = 3, Charge = 1, NeutralLoss = 0.0
				}
			}
		};

		// Act
		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry("KVFGR", 2, out MslLibraryEntry? loaded);

		// Assert
		Assert.That(found, Is.True);
		Assert.That(loaded!.GeneName, Is.Not.Null);
		Assert.That(loaded.GeneName, Is.EqualTo(string.Empty));
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Corruption detection
	// ─────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Overwriting the first four bytes of a valid .msl file with arbitrary bytes (not the
	/// MZLB magic) must cause <see cref="MslLibrary.Load"/> to throw
	/// <see cref="FormatException"/>. This prevents silent acceptance of wrong-format files.
	/// </summary>
	[Test]
	public void EdgeCase_CorruptedMagicBytes_ThrowsFormatException()
	{
		// Arrange — write a valid file first, then corrupt the magic bytes
		MslLibrary.Save(_tempMslPath, new[] { MakeMinimalEntry() });

		byte[] bytes = File.ReadAllBytes(_tempMslPath);
		bytes[0] = 0xFF; bytes[1] = 0xFF; bytes[2] = 0xFF; bytes[3] = 0xFF;
		File.WriteAllBytes(_tempMslPath, bytes);

		// Act / Assert
		Assert.Throws<FormatException>(() => MslLibrary.Load(_tempMslPath));
	}

	/// <summary>
	/// A file that is truncated mid-stream (all bytes after the header are removed) must
	/// cause <see cref="MslLibrary.Load"/> to throw <see cref="FormatException"/> or
	/// <see cref="InvalidDataException"/>. Either exception is acceptable; the test simply
	/// confirms no silent partial load occurs.
	/// </summary>
	[Test]
	public void EdgeCase_TruncatedFile_ThrowsFormatException()
	{
		// Arrange — write a valid file, then truncate it to just the header (64 bytes)
		MslLibrary.Save(_tempMslPath, new[] { MakeMinimalEntry() });

		byte[] allBytes = File.ReadAllBytes(_tempMslPath);
		// Truncate to header size only (MslFormat.HeaderSize = 64)
		byte[] truncated = allBytes.Take(MslFormat.HeaderSize).ToArray();
		File.WriteAllBytes(_tempMslPath, truncated);

		// Act / Assert — either FormatException or InvalidDataException is acceptable
		Assert.That(() => MslLibrary.Load(_tempMslPath),
			Throws.InstanceOf<FormatException>().Or.InstanceOf<InvalidDataException>());
	}

	/// <summary>
	/// Flipping a single bit in the fragment data section of an otherwise valid .msl file must
	/// cause <see cref="MslLibrary.Load"/> to throw <see cref="InvalidDataException"/> (CRC-32
	/// mismatch). This confirms that the CRC check catches single-bit corruption.
	/// </summary>
	[Test]
	public void EdgeCase_SingleFlippedDataBit_ThrowsInvalidDataException()
	{
		// Arrange — write a valid file
		MslLibrary.Save(_tempMslPath, new[] { MakeMinimalEntry() });

		byte[] bytes = File.ReadAllBytes(_tempMslPath);

		// Flip one bit in the middle of the file (well past the header, before the footer).
		// The footer occupies the last MslFormat.FooterSize (20) bytes.
		int corruptionOffset = bytes.Length / 2;
		bytes[corruptionOffset] ^= 0x01;
		File.WriteAllBytes(_tempMslPath, bytes);

		// Act / Assert
		Assert.Throws<InvalidDataException>(() => MslLibrary.Load(_tempMslPath));
	}

	/// <summary>
	/// Writing a version number other than a supported version into the file header (by direct
	/// byte manipulation) must cause <see cref="MslLibrary.Load"/> to throw
	/// <see cref="FormatException"/> with a message that mentions the version number.
	/// The reader supports versions 1 and <see cref="MslFormat.CurrentVersion"/> (2); version
	/// 99 is unsupported and must be rejected.
	/// </summary>
	[Test]
	public void EdgeCase_WrongVersion_ThrowsFormatExceptionWithVersionInMessage()
	{
		// Arrange — write a valid file, then overwrite the version int32 at offset 4
		MslLibrary.Save(_tempMslPath, new[] { MakeMinimalEntry() });

		byte[] bytes = File.ReadAllBytes(_tempMslPath);

		// Write version 99 (little-endian int32) at offset 4
		const int wrongVersion = 99;
		BitConverter.GetBytes(wrongVersion).CopyTo(bytes, 4);
		File.WriteAllBytes(_tempMslPath, bytes);

		// Act — load should throw FormatException
		var ex = Assert.Throws<FormatException>(() => MslLibrary.Load(_tempMslPath));

		// Assert — message should mention the bad version number
		Assert.That(ex!.Message, Does.Contain(wrongVersion.ToString()),
			"FormatException message should include the unexpected version number");
	}

	/// <summary>
	/// Overwriting the footer's NPrecursors field with a value that does not match the header's
	/// NPrecursors must cause <see cref="MslLibrary.Load"/> to throw
	/// <see cref="FormatException"/>. The footer/header cross-check exists specifically to
	/// detect truncated writes.
	/// </summary>
	[Test]
	public void EdgeCase_NprecursorsMismatch_ThrowsFormatException()
	{
		// Arrange — write a valid one-entry file
		MslLibrary.Save(_tempMslPath, new[] { MakeMinimalEntry() });

		byte[] bytes = File.ReadAllBytes(_tempMslPath);

		// The footer is the last MslFormat.FooterSize (20) bytes.
		// NPrecursors is at footer offset +8 (after OffsetTableOffset = 8 bytes).
		// We write 999 there to create a mismatch with the header's NPrecursors = 1.
		int footerStart = bytes.Length - MslFormat.FooterSize;
		int nPrecursorsAt = footerStart + 8; // OffsetTableOffset (8 bytes) comes first
		BitConverter.GetBytes(999).CopyTo(bytes, nPrecursorsAt);
		File.WriteAllBytes(_tempMslPath, bytes);

		// Act / Assert
		Assert.Throws<FormatException>(() => MslLibrary.Load(_tempMslPath));
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Fragmentation edge cases
	// ─────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Every named <see cref="MslFormat.NeutralLossCode"/> value (None, H2O, NH3, H3PO4,
	/// HPO3, PlusH2O) must survive the fragment-flags byte round-trip. The test writes one
	/// fragment per named loss type and verifies that the loaded NeutralLoss mass matches the
	/// written value within 1 × 10⁻⁴ Da.
	///
	/// Custom (arbitrary) neutral losses are handled via the extended annotation table
	/// introduced in format version 2; they are covered by <see cref="TestMslCustomNeutralLoss"/>.
	/// </summary>
	[Test]
	public void EdgeCase_AllNeutralLossValues_RoundTrip()
	{
		// All six named neutral-loss masses (None through PlusH2O).
		var lossValues = new (string name, double mass)[]
		{
			("None",     0.0),
			("H2O",    -18.010565),
			("NH3",    -17.026549),
			("H3PO4",  -97.976895),
			("HPO3",   -79.966331),
			("PlusH2O", -97.976895 + -18.010565),
		};

		// Build one entry with one fragment per named loss type at distinct m/z values
		var frags = new List<MslFragmentIon>();
		for (int i = 0; i < lossValues.Length; i++)
		{
			frags.Add(new MslFragmentIon
			{
				Mz = 200.0f + i * 50.0f,
				Intensity = 1.0f,
				ProductType = ProductType.y,
				FragmentNumber = i + 2,
				ResiduePosition = i + 2,
				Charge = 1,
				NeutralLoss = lossValues[i].mass
			});
		}

		var entry = new MslLibraryEntry
		{
			FullSequence = "NEUTRALLOSS",
			BaseSequence = "NEUTRALLOSS",
			PrecursorMz = 600.0f,
			ChargeState = 2,
			RetentionTime = 50.0,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			QValue = float.NaN,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			MatchedFragmentIons = frags
		};

		// Act
		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry("NEUTRALLOSS", 2, out MslLibraryEntry? loaded);

		// Assert
		Assert.That(found, Is.True);
		Assert.That(loaded!.MatchedFragmentIons.Count, Is.EqualTo(lossValues.Length));

		// Sort both by m/z ascending for positional comparison
		List<MslFragmentIon> origSorted = frags.OrderBy(f => f.Mz).ToList();
		List<MslFragmentIon> loadedSorted = loaded.MatchedFragmentIons.OrderBy(f => f.Mz).ToList();

		for (int i = 0; i < lossValues.Length; i++)
		{
			Assert.That(loadedSorted[i].NeutralLoss,
				Is.EqualTo(origSorted[i].NeutralLoss).Within(1e-4),
				$"NeutralLoss mismatch for named loss '{lossValues[i].name}'");
		}
	}

	/// <summary>
	/// A fragment ion whose neutral-loss mass does not match any of the five named codes
	/// (i.e. <see cref="MslFormat.NeutralLossCode.Custom"/> would be assigned) must now write
	/// and read back successfully via the extended annotation table introduced in format
	/// version 2. The custom mass must round-trip within 1 × 10⁻⁹ Da (double precision).
	///
	/// This test replaces the former <c>EdgeCase_CustomNeutralLoss_ThrowsNotSupportedException</c>
	/// which documented the pre-v2 limitation. Custom neutral losses are now fully supported;
	/// the dedicated test suite is in <see cref="TestMslCustomNeutralLoss"/>.
	/// </summary>
	[Test]
	public void EdgeCase_CustomNeutralLoss_WritesAndReadsBack()
	{
		// Arrange — arbitrary loss mass that ClassifyNeutralLoss maps to Custom
		const double customLoss = -42.010565; // acetyl loss, not a named code

		var entry = new MslLibraryEntry
		{
			FullSequence = "CUSTOMLOSS",
			BaseSequence = "CUSTOMLOSS",
			PrecursorMz = 400.0f,
			ChargeState = 2,
			RetentionTime = 30.0,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			QValue = float.NaN,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz              = 250.0f,
					Intensity       = 1.0f,
					ProductType     = ProductType.y,
					FragmentNumber  = 3,
					ResiduePosition = 3,
					Charge          = 1,
					NeutralLoss     = customLoss
				}
			}
		};

		// Act — must not throw
		Assert.DoesNotThrow(() => MslLibrary.Save(_tempMslPath, new[] { entry }),
			"Save must not throw for a custom neutral-loss fragment in format version 2");

		// Load and verify round-trip
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry("CUSTOMLOSS", 2, out MslLibraryEntry? loaded);

		Assert.That(found, Is.True);
		Assert.That(loaded!.MatchedFragmentIons.Count, Is.EqualTo(1));
		Assert.That(loaded.MatchedFragmentIons[0].NeutralLoss,
			Is.EqualTo(customLoss).Within(1e-9),
			"Custom neutral-loss mass must round-trip within double precision");
	}

	/// <summary>
	/// A fragment with <see cref="ProductType.D"/> (diagnostic ion) must have
	/// <see cref="MslFragmentIon.IsDiagnosticIon"/> equal to true after the round-trip.
	/// The is_diagnostic bit (bit 1 of the flags byte) must be set correctly by the writer and
	/// decoded by the reader.
	/// </summary>
	[Test]
	public void EdgeCase_DiagnosticIon_PreservedWithCorrectFlag()
	{
		// Arrange
		var entry = new MslLibraryEntry
		{
			FullSequence = "GLYCOPEP",
			BaseSequence = "GLYCOPEP",
			PrecursorMz = 550.0f,
			ChargeState = 2,
			RetentionTime = 45.0,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Glycopeptide,
			DissociationType = DissociationType.HCD,
			QValue = float.NaN,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 204.07f, Intensity = 1.0f,
					ProductType = ProductType.D,
					FragmentNumber = 1, ResiduePosition = 0, Charge = 1, NeutralLoss = 0.0
				}
			}
		};

		// Act
		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry("GLYCOPEP", 2, out MslLibraryEntry? loaded);

		// Assert
		Assert.That(found, Is.True);
		MslFragmentIon dIon = loaded!.MatchedFragmentIons.Single();
		Assert.That(dIon.IsDiagnosticIon, Is.True);
		Assert.That(dIon.ProductType, Is.EqualTo(ProductType.D));
	}

	/// <summary>
	/// A fragment with <see cref="ProductType.M"/> (immonium / precursor-related ion) must
	/// survive the round-trip with its ProductType preserved exactly.
	/// </summary>
	[Test]
	public void EdgeCase_MIon_Preserved()
	{
		var entry = MakeEntryWithSingleFragment(ProductType.M, 120.08f);

		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry(entry.FullSequence, entry.ChargeState, out MslLibraryEntry? loaded);

		Assert.That(found, Is.True);
		Assert.That(loaded!.MatchedFragmentIons.Single().ProductType, Is.EqualTo(ProductType.M));
	}

	/// <summary>
	/// A fragment with <see cref="ProductType.Y"/> (glycan Ycore / oxonium ion) must survive
	/// the round-trip with its ProductType preserved exactly.
	/// </summary>
	[Test]
	public void EdgeCase_GlycanYcoreIon_Preserved()
	{
		var entry = MakeEntryWithSingleFragment(ProductType.Y, 366.14f);

		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry(entry.FullSequence, entry.ChargeState, out MslLibraryEntry? loaded);

		Assert.That(found, Is.True);
		Assert.That(loaded!.MatchedFragmentIons.Single().ProductType, Is.EqualTo(ProductType.Y));
	}

	/// <summary>
	/// Writes one entry per ProductType value from 1 through 33 (skipping 0 / None) and verifies
	/// that each ProductType survives the binary round-trip exactly.
	///
	/// ProductType values are cast from int; any value that is not defined in the enum will
	/// still be stored and retrieved as the raw int16 — the test documents the expected
	/// round-trip behavior for all 33 assigned values.
	/// </summary>
	[Test]
	public void EdgeCase_AllProductTypes_1Through33_RoundTrip()
	{
		// Build one entry per ProductType value, each with a unique sequence to avoid key collisions
		var entries = new List<MslLibraryEntry>();
		for (int pt = 1; pt <= 33; pt++)
		{
			var productType = (ProductType)pt;
			string seq = $"PT{pt:D2}PEPTIDE"; // unique per value

			entries.Add(new MslLibraryEntry
			{
				FullSequence = seq,
				BaseSequence = seq,
				PrecursorMz = 300f + pt * 5f,
				ChargeState = 2,
				RetentionTime = 30.0 + pt,
				IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				QValue = float.NaN,
				ProteinAccession = string.Empty,
				ProteinName = string.Empty,
				GeneName = string.Empty,
				MatchedFragmentIons = new List<MslFragmentIon>
				{
					new MslFragmentIon
					{
						Mz              = 200f + pt,
						Intensity       = 1.0f,
						ProductType     = productType,
						FragmentNumber  = 2,
						ResiduePosition = 2,
						Charge          = 1,
						NeutralLoss     = 0.0
					}
				}
			});
		}

		// Act
		MslLibrary.Save(_tempMslPath, entries);
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);

		// Assert — each ProductType must round-trip correctly
		for (int pt = 1; pt <= 33; pt++)
		{
			var productType = (ProductType)pt;
			string seq = $"PT{pt:D2}PEPTIDE";

			bool found = lib.TryGetEntry(seq, 2, out MslLibraryEntry? loaded);
			Assert.That(found, Is.True, $"Entry for ProductType {pt} ({productType}) not found");
			Assert.That(loaded!.MatchedFragmentIons.Single().ProductType, Is.EqualTo(productType),
				$"ProductType mismatch for value {pt}");
		}
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Float precision
	// ─────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Precursor m/z is stored as float32 in the binary record. The round-trip value must
	/// equal the written value within single-precision float tolerance (1 × 10⁻⁴ Th).
	/// This test uses a m/z value with significant sub-float64 precision to confirm that
	/// only float32 rounding (not additional error) is introduced.
	/// </summary>
	[Test]
	public void EdgeCase_PrecursorMz_StoredAsSinglePrecisionFloat()
	{
		// double precision value that differs from its float32 representation
		const double highPrecisionMz = 500.123456789;
		float expectedFloat = (float)highPrecisionMz;

		var entry = MakeMinimalEntryWithMz((float)highPrecisionMz);

		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry(entry.FullSequence, entry.ChargeState, out MslLibraryEntry? loaded);

		Assert.That(found, Is.True);
		Assert.That((double)loaded!.PrecursorMz, Is.EqualTo(expectedFloat).Within(1e-4));
	}

	/// <summary>
	/// Fragment m/z is stored as float32. The round-trip value must equal the written float32
	/// value within 1 × 10⁻⁴ Th.
	/// </summary>
	[Test]
	public void EdgeCase_FragmentMz_StoredAsSinglePrecisionFloat()
	{
		const float fragMz = 357.1873f;
		var entry = MakeEntryWithSingleFragment(ProductType.y, fragMz);

		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry(entry.FullSequence, entry.ChargeState, out MslLibraryEntry? loaded);

		Assert.That(found, Is.True);
		Assert.That((double)loaded!.MatchedFragmentIons.Single().Mz, Is.EqualTo(fragMz).Within(1e-4));
	}

	/// <summary>
	/// iRT values can be negative (e.g. for hydrophilic peptides that elute before the iRT
	/// standard ladder anchors). A negative iRT must survive the float32 round-trip unchanged.
	/// </summary>
	[Test]
	public void EdgeCase_Irt_NegativeValue_Preserved()
	{
		const float negativeIrt = -12.34f;
		var entry = MakeMinimalEntry();
		entry.RetentionTime = negativeIrt;

		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry(entry.FullSequence, entry.ChargeState, out MslLibraryEntry? loaded);

		Assert.That(found, Is.True);
		Assert.That((double)loaded!.RetentionTime, Is.EqualTo(negativeIrt).Within(1e-4));
	}

	/// <summary>
	/// A q-value of <see cref="float.NaN"/> (indicating "not available") must survive the
	/// float32 round-trip as NaN. The test uses <see cref="float.IsNaN"/> to verify because
	/// NaN != NaN under IEEE 754 semantics.
	/// </summary>
	[Test]
	public void EdgeCase_QValue_NaN_Preserved()
	{
		var entry = MakeMinimalEntry();
		entry.QValue = float.NaN;

		MslLibrary.Save(_tempMslPath, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		bool found = lib.TryGetEntry(entry.FullSequence, entry.ChargeState, out MslLibraryEntry? loaded);

		Assert.That(found, Is.True);
		Assert.That(float.IsNaN(loaded!.QValue), Is.True,
			"QValue written as NaN must be read back as NaN");
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Private helpers
	// ─────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Constructs the simplest possible valid <see cref="MslLibraryEntry"/> (sequence "MINPEP",
	/// charge 2, one b2 fragment) for use in tests that only need a file to exist on disk.
	/// </summary>
	private static MslLibraryEntry MakeMinimalEntry() =>
		MakeMinimalEntryWithMz(300.0f);

	/// <summary>
	/// Constructs a minimal <see cref="MslLibraryEntry"/> with the specified precursor m/z.
	/// </summary>
	private static MslLibraryEntry MakeMinimalEntryWithMz(float precursorMz) =>
		new MslLibraryEntry
		{
			FullSequence = "MINPEP",
			BaseSequence = "MINPEP",
			PrecursorMz = precursorMz,
			ChargeState = 2,
			RetentionTime = 30.0,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			QValue = float.NaN,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 200.0f, Intensity = 1.0f, ProductType = ProductType.b,
					FragmentNumber = 2, ResiduePosition = 2, Charge = 1, NeutralLoss = 0.0
				}
			}
		};

	/// <summary>
	/// Constructs an <see cref="MslLibraryEntry"/> with a single fragment of the specified
	/// <see cref="ProductType"/> at the specified m/z. The entry's sequence is derived from
	/// the product type to avoid collisions when multiple entries are written in the same test.
	/// </summary>
	private static MslLibraryEntry MakeEntryWithSingleFragment(ProductType productType, float mz) =>
		new MslLibraryEntry
		{
			FullSequence = $"FRAG{(int)productType}",
			BaseSequence = $"FRAG{(int)productType}",
			PrecursorMz = 300f + (int)productType,
			ChargeState = 2,
			RetentionTime = 30.0,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			QValue = float.NaN,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz              = mz,
					Intensity       = 1.0f,
					ProductType     = productType,
					FragmentNumber  = 2,
					ResiduePosition = 2,
					Charge          = 1,
					NeutralLoss     = 0.0
				}
			}
		};
}