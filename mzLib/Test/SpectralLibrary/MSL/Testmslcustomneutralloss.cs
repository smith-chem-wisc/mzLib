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
/// NUnit 4 test suite for the custom neutral-loss extension table introduced in
/// .msl format version 2 (Prompt 11).
///
/// The extended annotation table stores arbitrary neutral-loss masses as a flat
/// <c>double[]</c> array at a new optional file section, referenced from fragment
/// records via the repurposed <c>ResiduePosition</c> field when
/// <c>NeutralLossCode.Custom</c> is stored in the flags byte.
///
/// All tests write to a dedicated temp directory and clean up in
/// <see cref="OneTimeTearDown"/>.
/// </summary>
[TestFixture]
public sealed class TestMslCustomNeutralLoss
{
	// ── Fixture setup ─────────────────────────────────────────────────────────

	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslCustomNeutralLossTests");

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
	/// Builds a minimal single-entry library with a single fragment that has the
	/// given neutral-loss mass.
	/// </summary>
	private static MslLibraryEntry EntryWithLoss(double neutralLoss,
		ProductType productType = ProductType.b,
		int fragmentNumber = 3,
		int residuePosition = 2,
		int charge = 1) =>
		new()
		{
			ModifiedSequence = "PEPTIDE",
			StrippedSequence = "PEPTIDE",
			PrecursorMz = 450.25,
			Charge = 2,
			Irt = 30.0,
			Fragments = new List<MslFragmentIon>
			{
				new()
				{
					Mz              = 250.13f,
					Intensity       = 1.0f,
					ProductType     = productType,
					FragmentNumber  = fragmentNumber,
					ResiduePosition = residuePosition,
					Charge          = charge,
					NeutralLoss     = neutralLoss
				}
			}
		};

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
	/// A single fragment with a custom glycan loss (HexNAc, −203.0794 Da) round-trips
	/// with the exact double-precision mass preserved within 1e-9 Da.
	/// </summary>
	[Test]
	public void CustomLoss_SingleEntry_RoundTrips_ExactMass()
	{
		const double hexNacLoss = -203.0794;
		string path = TempPath("single_custom");

		MslWriter.Write(path, new[] { EntryWithLoss(hexNacLoss) });
		using MslLibrary lib = MslLibrary.Load(path);

		Assert.That(lib.PrecursorCount, Is.EqualTo(1));
		bool found = lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? entry);
		Assert.That(found, Is.True);
		Assert.That(entry!.Fragments, Has.Count.EqualTo(1));
		Assert.That(entry.Fragments[0].NeutralLoss,
			Is.EqualTo(hexNacLoss).Within(1e-9));
	}

	/// <summary>
	/// Five fragments each with a distinct custom mass all survive round-trip;
	/// each mass is correct within 1e-9 Da.
	/// </summary>
	[Test]
	public void CustomLoss_MultipleDistinctMasses_AllRoundTrip()
	{
		double[] losses = { -203.0794, -162.0528, -291.0954, -365.1322, -111.0320 };
		string path = TempPath("multi_custom");

		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "GLYCOPEPTIDE",
			StrippedSequence = "GLYCOPEPTIDE",
			PrecursorMz = 700.35,
			Charge = 2,
			Irt = 45.0,
			Fragments = losses.Select((loss, i) => new MslFragmentIon
			{
				Mz = 200f + i * 50f,
				Intensity = 1.0f,
				ProductType = ProductType.b,
				FragmentNumber = i + 1,
				ResiduePosition = i,
				Charge = 1,
				NeutralLoss = loss
			}).ToList()
		};

		MslWriter.Write(path, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(path);

		bool found = lib.TryGetEntry("GLYCOPEPTIDE", 2, out MslLibraryEntry? readBack);
		Assert.That(found, Is.True);
		Assert.That(readBack!.Fragments, Has.Count.EqualTo(losses.Length));

		// Sort both by Mz to compare in a deterministic order (writer sorts by m/z)
		var readLosses = readBack.Fragments
			.OrderBy(f => f.Mz)
			.Select(f => f.NeutralLoss)
			.ToArray();
		var sortedExpected = losses.OrderByDescending(l => l)  // m/z ascending = loss ascending (least negative first)
			.ToArray();

		// Verify each mass matches within 1e-9 regardless of order — use a set check
		foreach (double expectedLoss in losses)
		{
			Assert.That(readLosses.Any(rl => Math.Abs(rl - expectedLoss) < 1e-9),
				Is.True,
				$"Expected loss {expectedLoss} not found in round-trip data.");
		}
	}

	/// <summary>
	/// Two fragments with identical custom masses produce an extended annotation table
	/// with only one entry (deduplicated), not two.
	/// </summary>
	[Test]
	public void CustomLoss_DeduplicatedInTable()
	{
		const double sharedLoss = -203.0794;
		string path = TempPath("dedup_custom");

		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "PEPTIDE",
			StrippedSequence = "PEPTIDE",
			PrecursorMz = 450.25,
			Charge = 2,
			Irt = 30.0,
			Fragments = new List<MslFragmentIon>
			{
				new() { Mz = 200f, Intensity = 1.0f, ProductType = ProductType.b,
					FragmentNumber = 2, ResiduePosition = 1, Charge = 1, NeutralLoss = sharedLoss },
				new() { Mz = 350f, Intensity = 0.8f, ProductType = ProductType.b,
					FragmentNumber = 4, ResiduePosition = 3, Charge = 1, NeutralLoss = sharedLoss }
			}
		};

		MslWriter.Write(path, new[] { entry });

		// Read the raw header to find ExtAnnotationTableOffset, then parse manually
		MslFileHeader header = ReadRawHeader(path);
		Assert.That((header.FileFlags & MslFormat.FileFlagHasExtAnnotations), Is.Not.Zero,
			"FileFlagHasExtAnnotations should be set");

		byte[] fileBytes = File.ReadAllBytes(path);
		int pos = header.ExtAnnotationTableOffset;
		int count = BitConverter.ToInt32(fileBytes, pos);

		// Exactly one unique custom loss — deduplication must have occurred.
		// NCustomLosses on disk includes the index-0 sentinel, so count == 2 (sentinel + 1 unique mass).
		Assert.That(count, Is.EqualTo(2),
			"Extended annotation table should have 2 entries (sentinel + 1 deduplicated custom loss) for 2 fragments with identical custom loss");
	}

	/// <summary>
	/// A library with H2O loss, NH3 loss, and a custom glycan loss all reads back correctly.
	/// </summary>
	[Test]
	public void CustomLoss_MixedNamedAndCustom_AllCorrect()
	{
		const double h2oLoss = -18.010565;
		const double nh3Loss = -17.026549;
		const double glycanLoss = -203.0794;

		string path = TempPath("mixed_named_custom");

		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "PEPTIDE",
			StrippedSequence = "PEPTIDE",
			PrecursorMz = 450.25,
			Charge = 2,
			Irt = 30.0,
			Fragments = new List<MslFragmentIon>
			{
				new() { Mz = 200f, Intensity = 1.0f, ProductType = ProductType.b,
					FragmentNumber = 2, ResiduePosition = 1, Charge = 1, NeutralLoss = h2oLoss },
				new() { Mz = 250f, Intensity = 0.9f, ProductType = ProductType.y,
					FragmentNumber = 3, ResiduePosition = 4, Charge = 1, NeutralLoss = nh3Loss },
				new() { Mz = 400f, Intensity = 0.7f, ProductType = ProductType.b,
					FragmentNumber = 5, ResiduePosition = 4, Charge = 1, NeutralLoss = glycanLoss }
			}
		};

		MslWriter.Write(path, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(path);

		bool found = lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? readBack);
		Assert.That(found, Is.True);

		var frags = readBack!.Fragments.OrderBy(f => f.Mz).ToList();
		Assert.That(frags, Has.Count.EqualTo(3));

		Assert.That(frags[0].NeutralLoss, Is.EqualTo(h2oLoss).Within(1e-6),
			"H2O loss fragment should decode correctly");
		Assert.That(frags[1].NeutralLoss, Is.EqualTo(nh3Loss).Within(1e-6),
			"NH3 loss fragment should decode correctly");
		Assert.That(frags[2].NeutralLoss, Is.EqualTo(glycanLoss).Within(1e-9),
			"Custom glycan loss fragment should decode correctly");
	}

	/// <summary>
	/// A fragment with NeutralLoss = 0.0 is stored as NeutralLossCode.None, not Custom;
	/// the extended annotation table section is absent from the file.
	/// </summary>
	[Test]
	public void CustomLoss_ZeroMass_NotCustom()
	{
		string path = TempPath("zero_loss");

		MslWriter.Write(path, new[] { EntryWithLoss(0.0) });

		MslFileHeader header = ReadRawHeader(path);
		Assert.That((header.FileFlags & MslFormat.FileFlagHasExtAnnotations), Is.Zero,
			"FileFlagHasExtAnnotations must not be set when NeutralLoss is 0.0");
		Assert.That(header.ExtAnnotationTableOffset, Is.Zero,
			"ExtAnnotationTableOffset must be 0 when no custom losses are present");
	}

	/// <summary>
	/// H3PO4 loss (−97.976895 Da, a named code) is stored in the 3-bit field, not in the
	/// extended annotation table; FileFlagHasExtAnnotations is 0.
	/// </summary>
	[Test]
	public void CustomLoss_NamedLoss_NotInTable()
	{
		const double h3po4Loss = -97.976895;
		string path = TempPath("named_phospho");

		MslWriter.Write(path, new[] { EntryWithLoss(h3po4Loss) });

		MslFileHeader header = ReadRawHeader(path);
		Assert.That((header.FileFlags & MslFormat.FileFlagHasExtAnnotations), Is.Zero,
			"Named H3PO4 loss must NOT set FileFlagHasExtAnnotations");
	}

	/// <summary>
	/// A file containing a custom-loss fragment has FileFlagHasExtAnnotations set in the
	/// header flags that the library exposes.
	/// </summary>
	[Test]
	public void CustomLoss_FileFlag_SetWhenPresent()
	{
		string path = TempPath("flag_set");

		MslWriter.Write(path, new[] { EntryWithLoss(-203.0794) });
		using MslLibrary lib = MslLibrary.Load(path);

		Assert.That((lib.Header.FileFlags & MslFormat.FileFlagHasExtAnnotations), Is.Not.Zero,
			"FileFlagHasExtAnnotations must be set in loaded header when custom losses present");
	}

	/// <summary>
	/// A file containing only named losses has FileFlagHasExtAnnotations clear in the
	/// loaded header.
	/// </summary>
	[Test]
	public void CustomLoss_FileFlag_ClearWhenAbsent()
	{
		string path = TempPath("flag_clear");

		MslWriter.Write(path, new[] { EntryWithLoss(-18.010565) }); // H2O — named loss
		using MslLibrary lib = MslLibrary.Load(path);

		Assert.That((lib.Header.FileFlags & MslFormat.FileFlagHasExtAnnotations), Is.Zero,
			"FileFlagHasExtAnnotations must be clear when only named losses are present");
	}

	/// <summary>Written files carry FormatVersion == 2.</summary>
	[Test]
	public void CustomLoss_Version2_WrittenCorrectly()
	{
		string path = TempPath("version2");

		// Write any library — version 2 is always used by the updated writer
		MslWriter.Write(path, new[] { EntryWithLoss(0.0) });

		MslFileHeader header = ReadRawHeader(path);
		Assert.That(header.FormatVersion, Is.EqualTo(2),
			"All files written by the updated writer must carry FormatVersion 2");
	}

	/// <summary>
	/// A version-1 file (no custom losses, FormatVersion 1) is accepted by the updated
	/// reader without throwing an exception.
	///
	/// Implementation: write a normal v2 file with no custom losses, then patch the
	/// FormatVersion field to 1. Because the patched byte is inside the CRC-covered region,
	/// the CRC is recomputed and re-written to the footer before attempting to read.
	/// </summary>
	[Test]
	public void CustomLoss_Version1_ReadableByVersion2Reader()
	{
		string path = TempPath("version1_compat");

		// Write a normal file with no custom losses (FileFlagHasExtAnnotations will be 0,
		// ExtAnnotationTableOffset will be 0 — already looks like a v1 file except version).
		MslWriter.Write(path, new[] { EntryWithLoss(0.0) });
		byte[] bytes = File.ReadAllBytes(path);

		// Patch FormatVersion at offset 4 (little-endian int32) from 2 → 1
		bytes[4] = 1; bytes[5] = 0; bytes[6] = 0; bytes[7] = 0;

		// The FormatVersion field is inside the CRC-covered region (before OffsetTableOffset),
		// so we must recompute and patch the CRC in the footer.
		int footerStart = bytes.Length - MslFormat.FooterSize;
		// OffsetTableOffset is the first 8 bytes of the footer
		long offsetTableOffset = BitConverter.ToInt64(bytes, footerStart);
		uint newCrc = MslWriter.ComputeCrc32OfArray(bytes, (int)offsetTableOffset);
		// DataCrc32 sits at footer offset +12 (after OffsetTableOffset:8 + NPrecursors:4)
		byte[] crcBytes = BitConverter.GetBytes(newCrc);
		Array.Copy(crcBytes, 0, bytes, footerStart + 12, 4);

		File.WriteAllBytes(path, bytes);

		// The reader must accept version-1 files without throwing
		Assert.DoesNotThrow(() =>
		{
			using MslLibrary lib = MslLibrary.Load(path);
			Assert.That(lib.PrecursorCount, Is.EqualTo(1));
		});
	}

	/// <summary>
	/// <c>header.ExtAnnotationTableOffset</c> points to exactly the right byte in the file.
	/// Verified by seeking there and reading NCustomLosses, which must equal 1.
	/// </summary>
	[Test]
	public void CustomLoss_TableOffset_MatchesActualPosition()
	{
		string path = TempPath("offset_check");

		MslWriter.Write(path, new[] { EntryWithLoss(-203.0794) });

		MslFileHeader header = ReadRawHeader(path);
		Assert.That((header.FileFlags & MslFormat.FileFlagHasExtAnnotations), Is.Not.Zero);
		Assert.That(header.ExtAnnotationTableOffset, Is.GreaterThan(0));

		byte[] fileBytes = File.ReadAllBytes(path);
		int pos = header.ExtAnnotationTableOffset;
		int count = BitConverter.ToInt32(fileBytes, pos);

		// NCustomLosses includes the index-0 sentinel, so count == 2 for one unique custom loss
		Assert.That(count, Is.EqualTo(2),
			"Reading NCustomLosses at the header-declared offset must yield 2 (sentinel + 1 custom loss)");

		// Index 0 is the sentinel (0.0); the real custom mass is at index 1 (offset +4 +8)
		double sentinel = BitConverter.ToDouble(fileBytes, pos + 4);
		Assert.That(sentinel, Is.EqualTo(0.0),
			"Index 0 must be the reserved sentinel (0.0)");

		double storedMass = BitConverter.ToDouble(fileBytes, pos + 4 + 8);
		Assert.That(storedMass, Is.EqualTo(-203.0794).Within(1e-9));
	}

	/// <summary>
	/// For a non-custom fragment, ResiduePosition round-trips correctly after the format change.
	/// </summary>
	[Test]
	public void CustomLoss_ResiduePosition_PreservedForNonCustom()
	{
		const int expectedResiduePos = 7;
		string path = TempPath("residue_pos_noncustom");

		MslWriter.Write(path, new[] { EntryWithLoss(0.0, residuePosition: expectedResiduePos) });
		using MslLibrary lib = MslLibrary.Load(path);

		lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.Fragments[0].ResiduePosition, Is.EqualTo(expectedResiduePos),
			"ResiduePosition must be preserved for fragments with no custom neutral loss");
	}

	/// <summary>
	/// For a custom-loss fragment, ResiduePosition is 0 after read-back (documented trade-off:
	/// the field is repurposed as ExtAnnotationIdx in the binary record).
	/// </summary>
	[Test]
	public void CustomLoss_ResiduePosition_ZeroForCustom()
	{
		string path = TempPath("residue_pos_custom");

		// Write with a non-zero ResiduePosition — it will be lost for custom-loss fragments
		MslWriter.Write(path, new[] { EntryWithLoss(-203.0794, residuePosition: 5) });
		using MslLibrary lib = MslLibrary.Load(path);

		lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.Fragments[0].ResiduePosition, Is.EqualTo(0),
			"ResiduePosition must be 0 for custom-loss fragments (repurposed as ExtAnnotationIdx in binary)");
	}

	/// <summary>
	/// The call that previously threw <c>NotSupportedException</c> for a custom neutral loss
	/// now completes without error.
	/// </summary>
	[Test]
	public void CustomLoss_FormerlyThrows_NowSucceeds()
	{
		string path = TempPath("formerly_throws");

		Assert.DoesNotThrow(() =>
			MslWriter.Write(path, new[] { EntryWithLoss(-203.0794) }),
			"Writing a custom neutral-loss fragment must no longer throw NotSupportedException");

		Assert.That(File.Exists(path), Is.True,
			"File must exist after successful write");
	}

	/// <summary>
	/// Corrupting a byte inside the extended annotation table causes the CRC check to throw
	/// <see cref="InvalidDataException"/> (the table is covered by the CRC).
	/// </summary>
	[Test]
	public void CustomLoss_CrcCoversExtAnnotationTable()
	{
		string path = TempPath("crc_ext_table");

		MslWriter.Write(path, new[] { EntryWithLoss(-203.0794) });

		// Find the start of the extended annotation table and corrupt the first double byte
		MslFileHeader header = ReadRawHeader(path);
		byte[] bytes = File.ReadAllBytes(path);

		// Corrupt the first byte of the stored custom-loss double (at offset +4 past NCustomLosses)
		int corruptOffset = header.ExtAnnotationTableOffset + 4;
		bytes[corruptOffset] ^= 0xFF;

		File.WriteAllBytes(path, bytes);

		// The reader must detect the corruption via CRC and throw InvalidDataException
		Assert.Throws<InvalidDataException>(() => MslLibrary.Load(path),
			"CRC-32 must detect corruption inside the extended annotation table");
	}
}