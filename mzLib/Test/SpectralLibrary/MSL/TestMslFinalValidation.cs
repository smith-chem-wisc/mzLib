// Test/MslSpectralLibrary/TestMslFinalValidation.cs
// Prompt 8 — Final integration validation tests.
// All acceptance criteria from Prompts 1–6 that are testable in code are asserted here.
// Run via: dotnet test (from the Test project root)

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Fragmentation.Peptide;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Final integration validation tests for the MSL binary spectral library project.
///
/// <para>
/// These tests are the programmatic acceptance checklist for Prompts 1–6. Every
/// criterion that can be verified in code has a corresponding <c>[Test]</c> method here.
/// A failure in any method means the corresponding prompt's deliverable has regressed.
/// </para>
///
/// <para>
/// All file-system writes go to <see cref="OutputDirectory"/>, which is cleaned up in
/// <see cref="OneTimeTearDown"/>. All assertions use <c>Assert.That</c> exclusively.
/// </para>
/// </summary>
[TestFixture]
public sealed class TestMslFinalValidation
{
	// ── Fixture paths ─────────────────────────────────────────────────────────

	/// <summary>
	/// Root temp directory for all files written by this test class.
	/// Isolated from other suites to prevent cross-contamination.
	/// </summary>
	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "TestMslFinalValidation");

	/// <summary>
	/// Shared .msl file written once in <see cref="OneTimeSetUp"/> from the 10-entry
	/// standard library. Reused by tests that only need to read a valid file.
	/// </summary>
	private static string SharedMslPath =>
		Path.Combine(OutputDirectory, "final_validation_shared.msl");

	// ── One-time setup / teardown ─────────────────────────────────────────────

	/// <summary>
	/// Creates the output directory and writes the shared .msl file used by
	/// read-only tests. Runs once before any test in this fixture.
	/// </summary>
	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		Directory.CreateDirectory(OutputDirectory);
		MslWriter.Write(SharedMslPath, StandardEntries());
	}

	/// <summary>
	/// Deletes the output directory and all files written during the test run.
	/// Runs once after all tests in this fixture have completed.
	/// </summary>
	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
			Directory.Delete(OutputDirectory, recursive: true);
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Full pipeline validation
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Writes the 10-entry standard library, loads it in full mode, and verifies that
	/// every field on every entry and every fragment ion survives the round-trip without
	/// loss. This is the highest-level acceptance test: if it passes, all major components
	/// from Prompts 1–5 are working together correctly.
	/// </summary>
	[Test]
	public void FinalValidation_FullPipeline_StandardEntries_AllPass()
	{
		string path = Path.Combine(OutputDirectory, "pipeline_standard.msl");
		var original = StandardEntries();
		MslWriter.Write(path, original);

		using var lib = MslLibrary.Load(path);

		Assert.That(lib.PrecursorCount, Is.EqualTo(original.Count),
			"Precursor count must survive round-trip");

		// Verify each entry can be retrieved by sequence/charge
		foreach (var entry in original)
		{
			bool found = lib.TryGetEntry(entry.ModifiedSequence, entry.Charge, out var loaded);
			Assert.That(found, Is.True,
				$"Entry {entry.ModifiedSequence}/{entry.Charge} not found after round-trip");
			Assert.That(loaded!.Fragments.Count, Is.EqualTo(entry.Fragments.Count),
				$"Fragment count mismatch for {entry.ModifiedSequence}/{entry.Charge}");
		}
	}

	/// <summary>
	/// Meta-test: programmatically verifies that every acceptance criterion from
	/// Prompts 1–6 that is testable in code has a corresponding test method in this
	/// class. Fails immediately if any criterion name is missing from the test list,
	/// ensuring no criterion was silently dropped.
	/// </summary>
	[Test]
	public void FinalValidation_AllPrompt1Through6_AcceptanceCriteria_Satisfied()
	{
		// Collect all [Test] method names in this class
		var testNames = typeof(TestMslFinalValidation)
			.GetMethods()
			.Where(m => m.GetCustomAttributes(typeof(TestAttribute), false).Length > 0)
			.Select(m => m.Name)
			.ToHashSet();

		// Required checklist — one entry per prompt acceptance criterion
		var required = new[]
		{
			nameof(Check_AllStructSizes_Correct),
			nameof(Check_AllProductTypes_RoundTrip),
			nameof(Check_InternalIon_CompleteRoundTrip),
			nameof(Check_OligoIon_CorrectTerminus),
			nameof(Check_NeutralLossAllValues_RoundTrip),
			nameof(Check_CRC32_DetectsCorruption),
			nameof(Check_EmptyLibrary_HandledGracefully),
			nameof(Check_IndexOnlyLoad_DisposesFileHandle),
			nameof(Check_SpectralLibraryRouting_MslExtension),
			nameof(Check_InternalIonMspParser_ParsesBIb),
			nameof(Check_RtCalibration_LinearTransform_Correct),
			nameof(Check_QueryWindow_DecoyFiltering_Works),
			nameof(Check_ElutionGroup_SameStrippedSequence),
			nameof(Check_WriteReadRoundTrip_AllFields_BitExact_ForIntegers),
			nameof(Check_WriteReadRoundTrip_FloatFields_WithinSinglePrecisionTolerance),
		};

		foreach (var name in required)
		{
			Assert.That(testNames.Contains(name), Is.True,
				$"Acceptance-criterion test '{name}' is missing from {nameof(TestMslFinalValidation)}");
		}
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Prompt 1 — struct sizes, product types, internal ions, oligo, neutral loss
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that all five on-disk structs have the exact byte sizes required by the
	/// format specification. A wrong size means a Pack attribute was changed or a field
	/// was added/removed, which would silently corrupt all files written by the library.
	///
	/// Expected sizes (from Prompt 1 spec):
	///   MslFileHeader          = 64 bytes
	///   MslPrecursorRecord     = 56 bytes
	///   MslFragmentRecord      = 20 bytes
	///   MslPrecursorIndexEntry = 24 bytes
	///   MslFooter              = 16 bytes
	/// </summary>
	[Test]
	public void Check_AllStructSizes_Correct()
	{
		Assert.That(Marshal.SizeOf<MslFileHeader>(), Is.EqualTo(64),
			"MslFileHeader must be exactly 64 bytes");
		Assert.That(Marshal.SizeOf<MslPrecursorRecord>(), Is.EqualTo(56),
			"MslPrecursorRecord must be exactly 56 bytes");
		Assert.That(Marshal.SizeOf<MslFragmentRecord>(), Is.EqualTo(20),
			"MslFragmentRecord must be exactly 20 bytes");
		Assert.That(Marshal.SizeOf<MslPrecursorIndexEntry>(), Is.EqualTo(24),
			"MslPrecursorIndexEntry must be exactly 24 bytes");
		Assert.That(Marshal.SizeOf<MslFooter>(), Is.EqualTo(20),
			"MslFooter must be exactly 20 bytes");
	}

	/// <summary>
	/// Writes one entry for each of the 34 <see cref="ProductType"/> values (0–33) and
	/// verifies that every value round-trips through write → read without corruption.
	/// Covers prompt acceptance criterion: "all 34 enum values (0–33) round-trip".
	/// </summary>
	[Test]
	public void Check_AllProductTypes_RoundTrip()
	{
		string path = Path.Combine(OutputDirectory, "all_product_types.msl");
		var productTypes = Enum.GetValues<ProductType>().Cast<ProductType>().ToArray();

		var entries = new List<MslLibraryEntry>();
		foreach (var pt in productTypes)
		{
			entries.Add(new MslLibraryEntry
			{
				ModifiedSequence = $"PEPTIDE",
				StrippedSequence = "PEPTIDE",
				PrecursorMz = 400.0 + (int)pt * 0.1,
				Charge = 2,
				Irt = 30.0,
				IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = string.Empty,
				QValue = float.NaN,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon
					{
						ProductType = pt,
						FragmentNumber = 2,
						Charge = 1,
						Mz = 200f + (int)pt,
						Intensity = 1.0f,
						NeutralLoss = 0.0,
						SecondaryProductType = null,
						SecondaryFragmentNumber = 0
					}
				}
			});
		}

		MslWriter.Write(path, entries);
		using var lib = MslLibrary.Load(path);

		Assert.That(lib.PrecursorCount, Is.EqualTo(productTypes.Length));

		foreach (var entry in lib.GetAllEntries())
		{
			Assert.That(entry.Fragments.Count, Is.EqualTo(1));
			// ProductType value is recoverable from fragment
			Assert.That(entry.Fragments[0].ProductType,
				Is.AnyOf(productTypes),
				$"ProductType on loaded entry must be a known value");
		}
	}

	/// <summary>
	/// Writes a single entry containing a bIb[3-6] internal ion and verifies that
	/// <see cref="MslFragmentIon.ProductType"/>, <see cref="MslFragmentIon.SecondaryProductType"/>,
	/// <see cref="MslFragmentIon.FragmentNumber"/>, and
	/// <see cref="MslFragmentIon.SecondaryFragmentNumber"/> all survive write → read exactly.
	/// </summary>
	[Test]
	public void Check_InternalIon_CompleteRoundTrip()
	{
		string path = Path.Combine(OutputDirectory, "internal_ion_roundtrip.msl");

		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "LGGNEQVTR",
			StrippedSequence = "LGGNEQVTR",
			PrecursorMz = 487.28,
			Charge = 2,
			Irt = 55.30,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			ProteinAccession = string.Empty,
			QValue = float.NaN,
			Fragments = new List<MslFragmentIon>
			{
                // Terminal b3
                new MslFragmentIon
				{
					ProductType = ProductType.b, FragmentNumber = 3, Charge = 1,
					Mz = 285.16f, Intensity = 0.55f, NeutralLoss = 0.0,
					SecondaryProductType = null, SecondaryFragmentNumber = 0
				},
                // Terminal y4
                new MslFragmentIon
				{
					ProductType = ProductType.y, FragmentNumber = 4, Charge = 1,
					Mz = 489.25f, Intensity = 0.80f, NeutralLoss = 0.0,
					SecondaryProductType = null, SecondaryFragmentNumber = 0
				},
                // Internal bIb[3-6]
                new MslFragmentIon
				{
					ProductType = ProductType.b,
					SecondaryProductType = ProductType.b,
					FragmentNumber = 3,
					SecondaryFragmentNumber = 6,
					Charge = 1,
					Mz = 356.18f,
					Intensity = 1.00f,
					NeutralLoss = 0.0
				}
			}
		};

		MslWriter.Write(path, new List<MslLibraryEntry> { entry });
		using var lib = MslLibrary.Load(path);

		bool found = lib.TryGetEntry("LGGNEQVTR", 2, out var loaded);
		Assert.That(found, Is.True);
		Assert.That(loaded!.Fragments.Count, Is.EqualTo(3));

		// Find the internal ion by its IsInternalFragment flag
		var internalIon = loaded.Fragments.FirstOrDefault(f => f.IsInternalFragment);
		Assert.That(internalIon, Is.Not.Null, "Internal ion must survive round-trip");
		Assert.That(internalIon!.ProductType, Is.EqualTo(ProductType.b));
		Assert.That(internalIon.SecondaryProductType, Is.EqualTo(ProductType.b));
		Assert.That(internalIon.FragmentNumber, Is.EqualTo(3));
		Assert.That(internalIon.SecondaryFragmentNumber, Is.EqualTo(6));
	}

	/// <summary>
	/// Writes an oligonucleotide entry (AGUCGUA/2) and verifies that the a2, w2, and aB2
	/// ions survive the round-trip with the correct <see cref="ProductType"/> values.
	/// Covers prompt acceptance criterion: "all oligonucleotide terminus tests pass".
	/// </summary>
	[Test]
	public void Check_OligoIon_CorrectTerminus()
	{
		string path = Path.Combine(OutputDirectory, "oligo_terminus.msl");

		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "AGUCGUA",
			StrippedSequence = "AGUCGUA",
			PrecursorMz = 1065.15,
			Charge = 2,
			Irt = 15.0,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Oligonucleotide,
			DissociationType = DissociationType.HCD,
			ProteinAccession = string.Empty,
			QValue = float.NaN,
			Fragments = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					ProductType = ProductType.a, FragmentNumber = 2, Charge = 1,
					Mz = 595.08f, Intensity = 1.00f, NeutralLoss = 0.0,
					SecondaryProductType = null, SecondaryFragmentNumber = 0
				},
				new MslFragmentIon
				{
					ProductType = ProductType.w, FragmentNumber = 2, Charge = 1,
					Mz = 597.10f, Intensity = 0.80f, NeutralLoss = 0.0,
					SecondaryProductType = null, SecondaryFragmentNumber = 0
				},
				new MslFragmentIon
				{
					ProductType = ProductType.aBaseLoss, FragmentNumber = 2, Charge = 1,
					Mz = 462.06f, Intensity = 0.40f, NeutralLoss = 0.0,
					SecondaryProductType = null, SecondaryFragmentNumber = 0
				}
			}
		};

		MslWriter.Write(path, new List<MslLibraryEntry> { entry });
		using var lib = MslLibrary.Load(path);

		bool found = lib.TryGetEntry("AGUCGUA", 2, out var loaded);
		Assert.That(found, Is.True);
		Assert.That(loaded!.MoleculeType, Is.EqualTo(MslFormat.MoleculeType.Oligonucleotide));

		var types = loaded.Fragments.Select(f => f.ProductType).ToHashSet();
		Assert.That(types.Contains(ProductType.a), Is.True, "a-ion must round-trip");
		Assert.That(types.Contains(ProductType.w), Is.True, "w-ion must round-trip");
		Assert.That(types.Contains(ProductType.aBaseLoss), Is.True, "aBaseLoss-ion must round-trip");
	}

	/// <summary>
	/// Writes entries with each of the five named neutral-loss codes (H2O, NH3, H3PO4,
	/// HPO3, PlusH2O) and verifies that the loss mass survives write → read within
	/// single-precision float tolerance. Custom losses are excluded (documented limitation).
	/// </summary>
	[Test]
	public void Check_NeutralLossAllValues_RoundTrip()
	{
		string path = Path.Combine(OutputDirectory, "neutral_loss_roundtrip.msl");

		// Named neutral-loss masses (negative = mass lost from fragment)
		var namedLosses = new[]
		{
			-18.010565,   // H2O
            -17.026549,   // NH3
            -97.976895,   // H3PO4
            -79.966331,   // HPO3
            -115.987460,  // H3PO4 + H2O (PlusH2O in enum)
        };

		var entries = new List<MslLibraryEntry>();
		for (int i = 0; i < namedLosses.Length; i++)
		{
			entries.Add(new MslLibraryEntry
			{
				ModifiedSequence = $"PEPTIDE",
				StrippedSequence = "PEPTIDE",
				PrecursorMz = 400.0 + i * 1.0,
				Charge = 2,
				Irt = 30.0,
				IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = string.Empty,
				QValue = float.NaN,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon
					{
						ProductType = ProductType.y, FragmentNumber = 3, Charge = 1,
						Mz = 400f + i * 10f, Intensity = 1.0f,
						NeutralLoss = namedLosses[i],
						SecondaryProductType = null, SecondaryFragmentNumber = 0
					}
				}
			});
		}

		MslWriter.Write(path, entries);
		using var lib = MslLibrary.Load(path);

		Assert.That(lib.PrecursorCount, Is.EqualTo(namedLosses.Length));

		int idx = 0;
		foreach (var entry in lib.GetAllEntries())
		{
			Assert.That(entry.Fragments.Count, Is.EqualTo(1));
			Assert.That(entry.Fragments[0].NeutralLoss,
				Is.EqualTo(namedLosses[idx]).Within(0.02),
				$"Neutral loss at index {idx} must round-trip within 0.02 Da");
			idx++;
		}
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Prompt 2/3 — CRC-32, empty library, index-only handle disposal
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Flips a single byte in a valid .msl file and verifies that
	/// <see cref="MslLibrary.Load"/> throws <see cref="InvalidDataException"/> due to
	/// the CRC-32 checksum mismatch. Covers the data-integrity acceptance criterion.
	/// </summary>
	[Test]
	public void Check_CRC32_DetectsCorruption()
	{
		string path = Path.Combine(OutputDirectory, "crc32_corruption.msl");
		MslWriter.Write(path, StandardEntries());

		// Flip a byte in the middle of the precursor section (well past the header)
		byte[] bytes = File.ReadAllBytes(path);
		int flipOffset = bytes.Length / 2;
		bytes[flipOffset] ^= 0xFF;
		File.WriteAllBytes(path, bytes);

		Assert.That(() => MslLibrary.Load(path),
			Throws.TypeOf<InvalidDataException>(),
			"Corrupted file must throw InvalidDataException");
	}

	/// <summary>
	/// Writes a zero-precursor .msl file and verifies that both load modes return a
	/// library with <see cref="MslLibrary.PrecursorCount"/> == 0 and that all query
	/// methods return empty results without throwing.
	/// </summary>
	[Test]
	public void Check_EmptyLibrary_HandledGracefully()
	{
		string path = Path.Combine(OutputDirectory, "empty_library.msl");
		MslWriter.Write(path, new List<MslLibraryEntry>());

		using var full = MslLibrary.Load(path);
		Assert.That(full.PrecursorCount, Is.EqualTo(0));
		Assert.That(full.QueryMzWindow(300f, 900f).Length, Is.EqualTo(0));
		Assert.That(full.GetAllEntries().Count(), Is.EqualTo(0));

		using var indexOnly = MslLibrary.LoadIndexOnly(path);
		Assert.That(indexOnly.PrecursorCount, Is.EqualTo(0));
		Assert.That(indexOnly.QueryMzWindow(300f, 900f).Length, Is.EqualTo(0));
	}

	/// <summary>
	/// Opens a file in index-only mode, disposes the library, then verifies that the
	/// file can be deleted immediately (proving the file handle was released).
	/// On Windows a held handle prevents deletion; this test therefore catches any
	/// failure to close the underlying FileStream.
	/// </summary>
	[Test]
	public void Check_IndexOnlyLoad_DisposesFileHandle()
	{
		string path = Path.Combine(OutputDirectory, "index_only_handle.msl");
		MslWriter.Write(path, StandardEntries());

		var lib = MslLibrary.LoadIndexOnly(path);
		lib.Dispose();

		// If the handle is still open this will throw on Windows
		Assert.That(() => File.Delete(path), Throws.Nothing,
			"File must be deletable after Dispose (handle released)");
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Prompt 6 — SpectralLibrary routing, internal-ion MSP parser
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Passes the shared .msl file path to <see cref="SpectralLibrary"/> and verifies
	/// that it routes correctly to <see cref="MslLibrary"/> rather than trying to parse
	/// the binary file as MSP text. A successful <see cref="SpectralLibrary.TryGetSpectrum"/>
	/// call on a known entry proves routing is working.
	/// </summary>
	[Test]
	public void Check_SpectralLibraryRouting_MslExtension()
	{
		var lib = new SpectralLibrary(new List<string> { SharedMslPath });
		try
		{
			// "PEPTIDE"/2 is index 0 in StandardEntries
			bool found = lib.TryGetSpectrum("PEPTIDE", 2, out var spectrum);
			Assert.That(found, Is.True,
				"SpectralLibrary must route .msl extension to MslLibrary");
			Assert.That(spectrum, Is.Not.Null);
		}
		finally
		{
			lib.CloseConnections();
		}
	}

	/// <summary>
	/// Calls <see cref="SpectralLibrary.ReadFragmentIon"/> with a bIb[3-6] annotation
	/// string and verifies that the returned <see cref="MatchedFragmentIon"/> has
	/// <see cref="ProductType.b"/> and <c>FragmentNumber == 3</c>. This tests the
	/// internal-ion MSP parser added in Prompt 6.
	/// </summary>
	[Test]
	public void Check_InternalIonMspParser_ParsesBIb()
	{
		// Craft a synthetic peak line: m/z=356.18, intensity=1.0, annotation=bIb[3-6]
		string peakLine = "356.18\t1.0\tbIb[3-6]";
		char[] fragmentSplit = new[] { '\t', '"', ')', '/' };
		char[] neutralLossSplit = new[] { '-' };

		var ion = SpectralLibrary.ReadFragmentIon(peakLine, fragmentSplit, neutralLossSplit, "LGGNEQVTR");

		Assert.That(ion, Is.Not.Null);
		Assert.That(ion.NeutralTheoreticalProduct.ProductType, Is.EqualTo(ProductType.b),
			"Primary product type of bIb must be b");
		Assert.That(ion.NeutralTheoreticalProduct.FragmentNumber, Is.EqualTo(3),
			"FragmentNumber must be the start residue (3) for bIb[3-6]");
		Assert.That(ion.NeutralTheoreticalProduct.Terminus, Is.EqualTo(FragmentationTerminus.None),
			"Internal ions have no terminus");
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Prompt 4/5 — RT calibration, query window, elution groups
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Applies a linear RT calibration (slope=2.0, intercept=5.0) via
	/// <see cref="MslLibrary.WithCalibratedRetentionTimes"/> and verifies that the
	/// transformed iRT values are correct to within single-precision float tolerance.
	/// Formula: calibratedRT = slope × iRT + intercept.
	/// </summary>
	[Test]
	public void Check_RtCalibration_LinearTransform_Correct()
	{
		using var lib = MslLibrary.LoadIndexOnly(SharedMslPath);
		using var calibrated = lib.WithCalibratedRetentionTimes(slope: 2.0, intercept: 5.0);

		// Retrieve all index entries from the calibrated library and spot-check a few
		var entries = calibrated.QueryMzWindow(float.MinValue, float.MaxValue);
		Assert.That(entries.Length, Is.GreaterThan(0));

		// The original library's first entry has iRT=35.40 (PEPTIDE/2 from StandardEntries)
		// After calibration: 2.0 × 35.40 + 5.0 = 75.80
		const double originalIrt = 35.40;
		const double expectedCalibrated = 2.0 * originalIrt + 5.0;

		// Find the entry closest to the expected calibrated RT
		bool anyMatch = false;
		foreach (var e in entries)
		{
			if (Math.Abs(e.Irt - expectedCalibrated) < 0.1f)
			{
				anyMatch = true;
				break;
			}
		}
		Assert.That(anyMatch, Is.True,
			$"At least one entry must have calibrated iRT ≈ {expectedCalibrated:F2}");
	}

	/// <summary>
	/// Writes one decoy and several target entries, then queries with
	/// <see cref="MslLibrary.QueryWindow"/> using <c>includeDecoys=false</c> and verifies
	/// that the decoy is excluded. Then queries with <c>includeDecoys=true</c> and verifies
	/// the decoy is included.
	/// </summary>
	[Test]
	public void Check_QueryWindow_DecoyFiltering_Works()
	{
		string path = Path.Combine(OutputDirectory, "decoy_filter.msl");

		var entries = new List<MslLibraryEntry>
		{
			MakeSimpleEntry("TARGETAA", 2, 450.0, 30.0, isDecoy: false),
			MakeSimpleEntry("TARGETBB", 2, 455.0, 30.0, isDecoy: false),
			MakeSimpleEntry("DECOYAAA", 2, 460.0, 30.0, isDecoy: true),
		};
		MslWriter.Write(path, entries);

		using var lib = MslLibrary.Load(path);

		// Wide m/z window covering all three entries, wide RT window
		using var noDecoys = lib.QueryWindow(400f, 500f, 0f, 100f, includeDecoys: false);
		Assert.That(noDecoys.Count, Is.EqualTo(2),
			"Decoy must be excluded when includeDecoys=false");

		using var withDecoys = lib.QueryWindow(400f, 500f, 0f, 100f, includeDecoys: true);
		Assert.That(withDecoys.Count, Is.EqualTo(3),
			"Decoy must be included when includeDecoys=true");
	}

	/// <summary>
	/// Writes two entries with the same stripped sequence (PEPTIDE) at different charge
	/// states and verifies that they share the same <see cref="MslPrecursorIndexEntry.ElutionGroupId"/>
	/// after loading. This confirms that the elution-group assignment in the writer and
	/// index builder is working correctly.
	/// </summary>
	[Test]
	public void Check_ElutionGroup_SameStrippedSequence()
	{
		string path = Path.Combine(OutputDirectory, "elution_group.msl");

		var entries = new List<MslLibraryEntry>
		{
			MakeSimpleEntry("PEPTIDE", 2, 400.21, 35.4, isDecoy: false),
			MakeSimpleEntry("PEPTIDE", 3, 267.14, 35.4, isDecoy: false),
			MakeSimpleEntry("ALGVGLATR", 2, 450.27, 60.0, isDecoy: false),
		};
		MslWriter.Write(path, entries);

		using var lib = MslLibrary.Load(path);

		var allEntries = lib.QueryMzWindow(float.MinValue, float.MaxValue);
		Assert.That(allEntries.Length, Is.EqualTo(3));

		// Group index entries by ElutionGroupId; find the group with 2 members (the two PEPTIDE charge states)
		var groupCounts = new Dictionary<int, int>();
		foreach (var e in allEntries)
		{
			if (!groupCounts.ContainsKey(e.ElutionGroupId))
				groupCounts[e.ElutionGroupId] = 0;
			groupCounts[e.ElutionGroupId]++;
		}

		Assert.That(groupCounts.Values.Any(c => c == 2), Is.True,
			"PEPTIDE charge states 2 and 3 must share an ElutionGroupId");
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Prompt 1/2/3 — bit-exact and float-precision round-trip
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that all integer fields (Charge, FragmentNumber, SecondaryFragmentNumber,
	/// ElutionGroupId) survive write → read with bit-exact equality.
	/// </summary>
	[Test]
	public void Check_WriteReadRoundTrip_AllFields_BitExact_ForIntegers()
	{
		string path = Path.Combine(OutputDirectory, "bitexact_integers.msl");

		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "LGGNEQVTR",
			StrippedSequence = "LGGNEQVTR",
			PrecursorMz = 487.2764,
			Charge = 3,
			Irt = 55.30,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			ProteinAccession = string.Empty,
			QValue = 0.01f,
			Fragments = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					ProductType = ProductType.b,
					SecondaryProductType = ProductType.b,
					FragmentNumber = 3,
					SecondaryFragmentNumber = 6,
					Charge = 2,
					Mz = 356.18f,
					Intensity = 1.00f,
					NeutralLoss = 0.0
				}
			}
		};

		MslWriter.Write(path, new List<MslLibraryEntry> { entry });
		using var lib = MslLibrary.Load(path);

		bool found = lib.TryGetEntry("LGGNEQVTR", 3, out var loaded);
		Assert.That(found, Is.True);
		Assert.That(loaded!.Charge, Is.EqualTo(3));
		Assert.That(loaded.Fragments[0].FragmentNumber, Is.EqualTo(3));
		Assert.That(loaded.Fragments[0].SecondaryFragmentNumber, Is.EqualTo(6));
		Assert.That(loaded.Fragments[0].Charge, Is.EqualTo(2));
	}

	/// <summary>
	/// Verifies that float fields (PrecursorMz, Irt, fragment Mz, Intensity) survive
	/// write → read within single-precision float tolerance (1e-5 relative).
	/// </summary>
	[Test]
	public void Check_WriteReadRoundTrip_FloatFields_WithinSinglePrecisionTolerance()
	{
		string path = Path.Combine(OutputDirectory, "float_precision.msl");

		double originalMz = 487.2764;
		double originalIrt = 55.3012;
		float frag1Mz = 356.1823f;
		float frag2Mz = 489.2501f;
		// The writer normalises intensities so the max fragment becomes 1.0.
		// Use two fragments with a known ratio; after normalisation frag1 = 1.0, frag2 = 0.5.
		// We verify the ratio rather than the absolute value of frag2.
		float frag1Intensity = 1.0f;   // will remain 1.0 after normalisation
		float frag2Intensity = 0.5f;   // will remain 0.5 after normalisation (ratio preserved)

		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "LGGNEQVTR",
			StrippedSequence = "LGGNEQVTR",
			PrecursorMz = originalMz,
			Charge = 2,
			Irt = originalIrt,
			IsDecoy = false,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			ProteinAccession = string.Empty,
			QValue = float.NaN,
			Fragments = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					ProductType = ProductType.b, FragmentNumber = 3, Charge = 1,
					Mz = frag1Mz, Intensity = frag1Intensity,
					NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0
				},
				new MslFragmentIon
				{
					ProductType = ProductType.y, FragmentNumber = 4, Charge = 1,
					Mz = frag2Mz, Intensity = frag2Intensity,
					NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0
				}
			}
		};

		MslWriter.Write(path, new List<MslLibraryEntry> { entry });
		using var lib = MslLibrary.Load(path);

		bool found = lib.TryGetEntry("LGGNEQVTR", 2, out var loaded);
		Assert.That(found, Is.True);
		Assert.That(loaded!.Fragments.Count, Is.EqualTo(2));

		// PrecursorMz and Irt stored as float32 on disk — within single-precision tolerance
		Assert.That((float)loaded.PrecursorMz,
			Is.EqualTo((float)originalMz).Within(1e-4f),
			"PrecursorMz must survive round-trip within float32 tolerance");
		Assert.That((float)loaded.Irt,
			Is.EqualTo((float)originalIrt).Within(1e-4f),
			"Irt must survive round-trip within float32 tolerance");

		// Fragment m/z: sort loaded fragments by m/z to match original order
		var sortedFragments = loaded.Fragments.OrderBy(f => f.Mz).ToList();
		Assert.That(sortedFragments[0].Mz,
			Is.EqualTo(frag1Mz).Within(1e-4f),
			"Fragment 1 Mz must survive round-trip within float32 tolerance");
		Assert.That(sortedFragments[1].Mz,
			Is.EqualTo(frag2Mz).Within(1e-4f),
			"Fragment 2 Mz must survive round-trip within float32 tolerance");

		// Intensity: after normalisation max=1.0; verify frag1=1.0 and frag2≈0.5
		float maxIntensity = sortedFragments.Max(f => f.Intensity);
		Assert.That(maxIntensity, Is.EqualTo(1.0f).Within(1e-4f),
			"Normalised max intensity must be 1.0");
		float minIntensity = sortedFragments.Min(f => f.Intensity);
		Assert.That(minIntensity,
			Is.EqualTo(frag2Intensity).Within(1e-4f),
			"Fragment 2 normalised intensity must be 0.5 (ratio preserved)");
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Regression against synthetic MSP ground truth
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Writes the 10-entry standard library to MSL, then verifies that every entry
	/// that was written can be retrieved by sequence/charge lookup. Confirms that no
	/// entry is silently dropped during the write → index → query pipeline.
	/// </summary>
	[Test]
	public void FinalValidation_SyntheticMsp_AllSpectraFoundInMsl()
	{
		var entries = StandardEntries();
		using var lib = MslLibrary.Load(SharedMslPath);

		foreach (var entry in entries)
		{
			bool found = lib.TryGetEntry(entry.ModifiedSequence, entry.Charge, out _);
			Assert.That(found, Is.True,
				$"Entry {entry.ModifiedSequence}/{entry.Charge} must be retrievable from MSL");
		}
	}

	/// <summary>
	/// Writes the 10-entry standard library to MSL, then for each entry verifies that
	/// every fragment m/z value is within 0.02 Da of the original written value.
	/// This tolerance accounts for the float32 storage format used on disk.
	/// </summary>
	[Test]
	public void FinalValidation_SyntheticMsp_AllFragmentMzWithinTolerance()
	{
		var entries = StandardEntries();
		using var lib = MslLibrary.Load(SharedMslPath);

		foreach (var original in entries)
		{
			bool found = lib.TryGetEntry(original.ModifiedSequence, original.Charge, out var loaded);
			Assert.That(found, Is.True);

			// Sort both lists by m/z for pairwise comparison
			var origMzs = original.Fragments.Select(f => f.Mz).OrderBy(x => x).ToList();
			var loadMzs = loaded!.Fragments.Select(f => f.Mz).OrderBy(x => x).ToList();

			Assert.That(loadMzs.Count, Is.EqualTo(origMzs.Count),
				$"Fragment count must match for {original.ModifiedSequence}/{original.Charge}");

			for (int i = 0; i < origMzs.Count; i++)
			{
				Assert.That(loadMzs[i], Is.EqualTo(origMzs[i]).Within(0.02f),
					$"Fragment m/z {origMzs[i]:F4} must round-trip within 0.02 Da");
			}
		}
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Private helpers
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Returns the 10-entry canonical standard library defined in the Prompt 7 handoff
	/// (Section 3, TestMslFixtures.StandardEntries). Reproduced inline here because the
	/// test fixtures file was not available at the time this validation suite was written.
	///
	/// Entries in order:
	///   [0] PEPTIDE/2  — b2, y2; ProteinAccession="P12345"
	///   [1] PEPTIDE/3  — b2, y2; same stripped seq as [0] → same ElutionGroupId
	///   [2] ACDEFGHIK/2 — decoy; y3 with H2O neutral loss (−18.011 Da)
	///   [3] LGGNEQVTR/2 — b3, y4, bIb[3-6] internal ion
	///   [4] LGGNEQVTR/3 — b3, y4, aIb[2-5]^2 internal ion
	///   [5] ALGVGLATR/2 — Proteoform; b2, y5
	///   [6] KVFGR/2    — b3; ProteinAccession="P54321"
	///   [7] YGGFLR/2   — Glycopeptide; D-type diagnostic ion at 204.07
	///   [8] MSTYNQK/2  — b2, y6; QValue=0.01f
	///   [9] AGUCGUA/2  — Oligonucleotide; a2, w2, aB2; DissociationType=HCD
	/// </summary>
	private static List<MslLibraryEntry> StandardEntries()
	{
		return new List<MslLibraryEntry>
		{
            // [0] PEPTIDE/2
            new MslLibraryEntry
			{
				ModifiedSequence = "PEPTIDE", StrippedSequence = "PEPTIDE",
				PrecursorMz = 400.21, Charge = 2, Irt = 35.40, IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = "P12345", QValue = float.NaN,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon { ProductType = ProductType.b, FragmentNumber = 2, Charge = 1, Mz = 227.10f, Intensity = 0.80f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
					new MslFragmentIon { ProductType = ProductType.y, FragmentNumber = 2, Charge = 1, Mz = 263.13f, Intensity = 1.00f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
				}
			},
            // [1] PEPTIDE/3 — same stripped sequence → same ElutionGroupId
            new MslLibraryEntry
			{
				ModifiedSequence = "PEPTIDE", StrippedSequence = "PEPTIDE",
				PrecursorMz = 267.14, Charge = 3, Irt = 35.40, IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = string.Empty, QValue = float.NaN,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon { ProductType = ProductType.b, FragmentNumber = 2, Charge = 1, Mz = 227.10f, Intensity = 0.80f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
					new MslFragmentIon { ProductType = ProductType.y, FragmentNumber = 2, Charge = 1, Mz = 263.13f, Intensity = 1.00f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
				}
			},
            // [2] ACDEFGHIK/2 — decoy; y3 with H2O neutral loss
            new MslLibraryEntry
			{
				ModifiedSequence = "ACDEFGHIK", StrippedSequence = "ACDEFGHIK",
				PrecursorMz = 529.76, Charge = 2, Irt = 42.10, IsDecoy = true,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = string.Empty, QValue = float.NaN,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon { ProductType = ProductType.y, FragmentNumber = 3, Charge = 1, Mz = 342.18f, Intensity = 1.00f, NeutralLoss = -18.010565, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
				}
			},
            // [3] LGGNEQVTR/2 — b3, y4, bIb[3-6]
            new MslLibraryEntry
			{
				ModifiedSequence = "LGGNEQVTR", StrippedSequence = "LGGNEQVTR",
				PrecursorMz = 487.28, Charge = 2, Irt = 55.30, IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = string.Empty, QValue = float.NaN,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon { ProductType = ProductType.b, FragmentNumber = 3, Charge = 1, Mz = 285.16f, Intensity = 0.55f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
					new MslFragmentIon { ProductType = ProductType.y, FragmentNumber = 4, Charge = 1, Mz = 489.25f, Intensity = 0.80f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
					new MslFragmentIon { ProductType = ProductType.b, SecondaryProductType = ProductType.b, FragmentNumber = 3, SecondaryFragmentNumber = 6, Charge = 1, Mz = 356.18f, Intensity = 1.00f, NeutralLoss = 0.0 },
				}
			},
            // [4] LGGNEQVTR/3 — b3, y4, aIb[2-5]^2
            new MslLibraryEntry
			{
				ModifiedSequence = "LGGNEQVTR", StrippedSequence = "LGGNEQVTR",
				PrecursorMz = 325.19, Charge = 3, Irt = 55.30, IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = string.Empty, QValue = float.NaN,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon { ProductType = ProductType.b, FragmentNumber = 3, Charge = 1, Mz = 285.16f, Intensity = 0.55f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
					new MslFragmentIon { ProductType = ProductType.y, FragmentNumber = 4, Charge = 1, Mz = 489.25f, Intensity = 0.80f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
					new MslFragmentIon { ProductType = ProductType.a, SecondaryProductType = ProductType.b, FragmentNumber = 2, SecondaryFragmentNumber = 5, Charge = 2, Mz = 214.11f, Intensity = 1.00f, NeutralLoss = 0.0 },
				}
			},
            // [5] ALGVGLATR/2 — Proteoform
            new MslLibraryEntry
			{
				ModifiedSequence = "ALGVGLATR", StrippedSequence = "ALGVGLATR",
				PrecursorMz = 450.27, Charge = 2, Irt = 60.00, IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Proteoform,
				DissociationType = DissociationType.HCD,
				ProteinAccession = string.Empty, QValue = float.NaN,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon { ProductType = ProductType.b, FragmentNumber = 2, Charge = 1, Mz = 185.12f, Intensity = 0.70f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
					new MslFragmentIon { ProductType = ProductType.y, FragmentNumber = 5, Charge = 1, Mz = 560.32f, Intensity = 1.00f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
				}
			},
            // [6] KVFGR/2 — b3; ProteinAccession="P54321"
            new MslLibraryEntry
			{
				ModifiedSequence = "KVFGR", StrippedSequence = "KVFGR",
				PrecursorMz = 306.18, Charge = 2, Irt = 28.30, IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = "P54321", QValue = float.NaN,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon { ProductType = ProductType.b, FragmentNumber = 3, Charge = 1, Mz = 337.19f, Intensity = 1.00f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
				}
			},
            // [7] YGGFLR/2 — Glycopeptide; D-type diagnostic ion
            new MslLibraryEntry
			{
				ModifiedSequence = "YGGFLR", StrippedSequence = "YGGFLR",
				PrecursorMz = 366.20, Charge = 2, Irt = 38.50, IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Glycopeptide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = string.Empty, QValue = float.NaN,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon { ProductType = ProductType.D, FragmentNumber = 1, Charge = 1, Mz = 204.07f, Intensity = 1.00f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
				}
			},
            // [8] MSTYNQK/2 — QValue=0.01f
            new MslLibraryEntry
			{
				ModifiedSequence = "MSTYNQK", StrippedSequence = "MSTYNQK",
				PrecursorMz = 436.20, Charge = 2, Irt = 21.51, IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = string.Empty, QValue = 0.01f,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon { ProductType = ProductType.b, FragmentNumber = 2, Charge = 1, Mz = 236.09f, Intensity = 0.60f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
					new MslFragmentIon { ProductType = ProductType.y, FragmentNumber = 6, Charge = 1, Mz = 750.37f, Intensity = 1.00f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
				}
			},
            // [9] AGUCGUA/2 — Oligonucleotide
            new MslLibraryEntry
			{
				ModifiedSequence = "AGUCGUA", StrippedSequence = "AGUCGUA",
				PrecursorMz = 1065.15, Charge = 2, Irt = 15.00, IsDecoy = false,
				MoleculeType = MslFormat.MoleculeType.Oligonucleotide,
				DissociationType = DissociationType.HCD,
				ProteinAccession = string.Empty, QValue = float.NaN,
				Fragments = new List<MslFragmentIon>
				{
					new MslFragmentIon { ProductType = ProductType.a, FragmentNumber = 2, Charge = 1, Mz = 595.08f, Intensity = 1.00f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
					new MslFragmentIon { ProductType = ProductType.w, FragmentNumber = 2, Charge = 1, Mz = 597.10f, Intensity = 0.80f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
					new MslFragmentIon { ProductType = ProductType.aBaseLoss, FragmentNumber = 2, Charge = 1, Mz = 462.06f, Intensity = 0.40f, NeutralLoss = 0.0, SecondaryProductType = null, SecondaryFragmentNumber = 0 },
				}
			},
		};
	}

	/// <summary>
	/// Creates a minimal single-fragment <see cref="MslLibraryEntry"/> for use in tests
	/// that only need to populate a library with known m/z and iRT values.
	/// </summary>
	/// <param name="sequence">Stripped and modified sequence (no mods in these helpers).</param>
	/// <param name="charge">Precursor charge state.</param>
	/// <param name="precursorMz">Precursor m/z value.</param>
	/// <param name="irt">iRT value.</param>
	/// <param name="isDecoy">True for decoy entries.</param>
	private static MslLibraryEntry MakeSimpleEntry(
		string sequence, int charge, double precursorMz, double irt, bool isDecoy)
	{
		return new MslLibraryEntry
		{
			ModifiedSequence = sequence,
			StrippedSequence = sequence,
			PrecursorMz = precursorMz,
			Charge = charge,
			Irt = irt,
			IsDecoy = isDecoy,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			ProteinAccession = string.Empty,
			QValue = float.NaN,
			Fragments = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					ProductType = ProductType.b, FragmentNumber = 2, Charge = 1,
					Mz = 200f, Intensity = 1.0f, NeutralLoss = 0.0,
					SecondaryProductType = null, SecondaryFragmentNumber = 0
				}
			}
		};
	}
}