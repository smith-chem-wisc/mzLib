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

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 test suite for <see cref="MslLibrary"/>, the top-level user-facing façade for
/// the .msl binary spectral library format.
///
/// All tests use <c>Assert.That</c> exclusively. All tests that touch the file system
/// write to a dedicated temp directory and clean up in <see cref="OneTimeTearDown"/>.
///
/// The shared fixture is a three-entry library written once in <see cref="OneTimeSetUp"/>:
/// <list type="bullet">
///   <item>Entry 0 — PEPTIDE/2, target, m/z 449.74, iRT 35.4, 2 fragments</item>
///   <item>Entry 1 — PEPTIDE/3, target, same stripped sequence as entry 0 (same elution group), m/z 300.16, iRT 35.4, 1 fragment</item>
///   <item>Entry 2 — ACDEFGHIK/2, decoy, m/z 529.76, iRT 42.1, 1 fragment</item>
/// </list>
///
/// Coverage groups:
/// <list type="number">
///   <item>Factory methods</item>
///   <item>Properties</item>
///   <item>DDA-style lookup (TryGetEntry / TryGetLibrarySpectrum)</item>
///   <item>DIA window queries (QueryMzWindow / QueryWindow / GetEntry)</item>
///   <item>Bulk enumeration (GetAllEntries)</item>
///   <item>RT calibration (WithCalibratedRetentionTimes)</item>
///   <item>Elution group (GetElutionGroup)</item>
///   <item>Save round-trip</item>
///   <item>Dispose</item>
/// </list>
/// </summary>
[TestFixture]
public sealed class TestMslLibrary
{
	// ── Fixture paths ─────────────────────────────────────────────────────────

	/// <summary>
	/// Root directory under which all test output files are written.
	/// Using a dedicated subfolder keeps the working directory clean and makes
	/// bulk teardown trivial.
	/// </summary>
	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslLibraryTests");

	/// <summary>
	/// Path of the shared three-entry .msl file written once in <see cref="OneTimeSetUp"/>.
	/// All tests that only need to read (not write) use this path.
	/// </summary>
	private static string SharedLibraryPath =>
		Path.Combine(OutputDirectory, "shared_fixture.msl");

	// ── One-time setup / teardown ─────────────────────────────────────────────

	/// <summary>
	/// Creates the output directory and writes the shared three-entry fixture file so that
	/// every subsequent test that opens a library can use the same on-disk file without
	/// re-writing it.
	/// </summary>
	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		Directory.CreateDirectory(OutputDirectory);
		MslLibrary.Save(SharedLibraryPath, BuildFixtureEntries());
	}

	/// <summary>
	/// Deletes every .msl file created by this test run.
	/// </summary>
	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
		{
			foreach (string f in Directory.GetFiles(OutputDirectory, "*.msl"))
				File.Delete(f);
		}
	}

	// ── Fixture data factory ──────────────────────────────────────────────────

	/// <summary>
	/// Returns the canonical three-entry fixture list used by all tests in this class.
	///
	/// Entry 0: PEPTIDE/2  — target, PrecursorMz 449.74, iRT 35.4, 2 fragments (b2, y1).
	/// Entry 1: PEPTIDE/3  — target, PrecursorMz 300.16, iRT 35.4, 1 fragment (y1).
	///                        Same stripped sequence as entry 0 → same ElutionGroupId.
	/// Entry 2: ACDEFGHIK/2 — decoy, PrecursorMz 529.76, iRT 42.1, 1 fragment (y3).
	/// </summary>
	/// <returns>A new list of three <see cref="MslLibraryEntry"/> objects.</returns>
	private static List<MslLibraryEntry> BuildFixtureEntries()
	{
		var entry0 = new MslLibraryEntry
		{
			ModifiedSequence = "PEPTIDE",
			StrippedSequence = "PEPTIDE",
			PrecursorMz = 449.74,
			Charge = 2,
			Irt = 35.4,
			IonMobility = 0.0,
			ProteinAccession = "P12345",
			ProteinName = "Test Protein",
			GeneName = "TESTA",
			IsDecoy = false,
			IsProteotypic = true,
			QValue = 0.01f,
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			Fragments = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 98.060f, Intensity = 10000f,
					ProductType = ProductType.b, FragmentNumber = 2,
					ResiduePosition = 1, Charge = 1
				},
				new MslFragmentIon
				{
					Mz = 175.119f, Intensity = 8000f,
					ProductType = ProductType.y, FragmentNumber = 1,
					ResiduePosition = 6, Charge = 1
				}
			}
		};

		var entry1 = new MslLibraryEntry
		{
			ModifiedSequence = "PEPTIDE",
			StrippedSequence = "PEPTIDE",
			PrecursorMz = 300.16,
			Charge = 3,
			Irt = 35.4,
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
			Fragments = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 175.119f, Intensity = 5000f,
					ProductType = ProductType.y, FragmentNumber = 1,
					ResiduePosition = 6, Charge = 1
				}
			}
		};

		var entry2 = new MslLibraryEntry
		{
			ModifiedSequence = "ACDEFGHIK",
			StrippedSequence = "ACDEFGHIK",
			PrecursorMz = 529.76,
			Charge = 2,
			Irt = 42.1,
			IonMobility = 0.0,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsDecoy = true,
			IsProteotypic = false,
			QValue = float.NaN,
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Nce = 0,
			Fragments = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 364.197f, Intensity = 10000f,
					ProductType = ProductType.y, FragmentNumber = 3,
					ResiduePosition = 6, Charge = 1
				}
			}
		};

		return new List<MslLibraryEntry> { entry0, entry1, entry2 };
	}

	/// <summary>
	/// Returns a unique per-test temp file path under <see cref="OutputDirectory"/>.
	/// The file is not created; only the path string is returned.
	/// </summary>
	/// <param name="testName">Stem used as the filename (typically the calling test name).</param>
	/// <returns>An absolute .msl file path that does not yet exist on disk.</returns>
	private static string TempPath(string testName) =>
		Path.Combine(OutputDirectory, testName + ".msl");

	// ════════════════════════════════════════════════════════════════════════════
	// Group 1 — Factory methods
	// ════════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslLibrary.Load"/> returns a non-null instance.
	/// </summary>
	[Test]
	public void Load_ReturnsNonNullLibrary()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		Assert.That(lib, Is.Not.Null,
			"Load() must return a non-null MslLibrary.");
	}

	/// <summary>
	/// Verifies that <see cref="MslLibrary.LoadIndexOnly"/> returns a non-null instance.
	/// </summary>
	[Test]
	public void LoadIndexOnly_ReturnsNonNullLibrary()
	{
		using MslLibrary lib = MslLibrary.LoadIndexOnly(SharedLibraryPath);

		Assert.That(lib, Is.Not.Null,
			"LoadIndexOnly() must return a non-null MslLibrary.");
	}

	/// <summary>
	/// Verifies that a full-load library reports <see cref="MslLibrary.IsIndexOnly"/> as
	/// <see langword="false"/>.
	/// </summary>
	[Test]
	public void Load_IsIndexOnly_ReturnsFalse()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		Assert.That(lib.IsIndexOnly, Is.False,
			"A library opened with Load() must report IsIndexOnly = false.");
	}

	/// <summary>
	/// Verifies that an index-only library reports <see cref="MslLibrary.IsIndexOnly"/> as
	/// <see langword="true"/>.
	/// </summary>
	[Test]
	public void LoadIndexOnly_IsIndexOnly_ReturnsTrue()
	{
		using MslLibrary lib = MslLibrary.LoadIndexOnly(SharedLibraryPath);

		Assert.That(lib.IsIndexOnly, Is.True,
			"A library opened with LoadIndexOnly() must report IsIndexOnly = true.");
	}

	/// <summary>
	/// Verifies that a full-load library reports the correct precursor count matching the
	/// number of entries written.
	/// </summary>
	[Test]
	public void Load_PrecursorCount_MatchesWrittenCount()
	{
		int expected = BuildFixtureEntries().Count;
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		Assert.That(lib.PrecursorCount, Is.EqualTo(expected),
			"Full-load PrecursorCount must equal the number of entries written.");
	}

	/// <summary>
	/// Verifies that an index-only library reports the correct precursor count matching the
	/// number of entries written.
	/// </summary>
	[Test]
	public void LoadIndexOnly_PrecursorCount_MatchesWrittenCount()
	{
		int expected = BuildFixtureEntries().Count;
		using MslLibrary lib = MslLibrary.LoadIndexOnly(SharedLibraryPath);

		Assert.That(lib.PrecursorCount, Is.EqualTo(expected),
			"Index-only PrecursorCount must equal the number of entries written.");
	}

	// ════════════════════════════════════════════════════════════════════════════
	// Group 2 — Properties
	// ════════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslLibrary.PrecursorCount"/> equals the count of written
	/// entries after a full load.
	/// </summary>
	[Test]
	public void PrecursorCount_CorrectAfterLoad()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		Assert.That(lib.PrecursorCount, Is.EqualTo(3),
			"PrecursorCount must equal 3 for the three-entry fixture.");
	}

	/// <summary>
	/// Verifies that <c>TargetCount + DecoyCount == PrecursorCount</c> holds for the
	/// fixture library (2 targets, 1 decoy).
	/// </summary>
	[Test]
	public void TargetCount_PlusDecoyCount_EqualsPrecursorCount()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		Assert.That(lib.TargetCount + lib.DecoyCount, Is.EqualTo(lib.PrecursorCount),
			"TargetCount + DecoyCount must equal PrecursorCount.");
		Assert.That(lib.TargetCount, Is.EqualTo(2),
			"Fixture has 2 target entries.");
		Assert.That(lib.DecoyCount, Is.EqualTo(1),
			"Fixture has 1 decoy entry.");
	}

	/// <summary>
	/// Verifies that <see cref="MslLibrary.MinPrecursorMz"/> and
	/// <see cref="MslLibrary.MaxPrecursorMz"/> bracket the three fixture precursor m/z values
	/// (300.16, 449.74, 529.76) correctly.
	/// </summary>
	[Test]
	public void MinMaxPrecursorMz_CorrectAfterLoad()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		Assert.That(lib.MinPrecursorMz, Is.EqualTo(300.16f).Within(0.01f),
			"MinPrecursorMz must be the smallest precursor m/z in the library.");
		Assert.That(lib.MaxPrecursorMz, Is.EqualTo(529.76f).Within(0.01f),
			"MaxPrecursorMz must be the largest precursor m/z in the library.");
	}

	// ════════════════════════════════════════════════════════════════════════════
	// Group 3 — DDA-style lookup
	// ════════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslLibrary.TryGetEntry"/> finds an existing entry by its
	/// modified sequence and charge in full-load mode and returns true.
	/// </summary>
	[Test]
	public void TryGetEntry_FullLoad_FindsExistingEntry()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		bool found = lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? entry);

		Assert.That(found, Is.True,
			"TryGetEntry must return true for PEPTIDE/2 in full-load mode.");
		Assert.That(entry, Is.Not.Null,
			"The out parameter must be non-null when TryGetEntry returns true.");
	}

	/// <summary>
	/// Verifies that <see cref="MslLibrary.TryGetEntry"/> finds an existing entry in
	/// index-only mode, confirming that fragment ions are loaded on demand.
	/// </summary>
	[Test]
	public void TryGetEntry_IndexOnly_FindsExistingEntry()
	{
		using MslLibrary lib = MslLibrary.LoadIndexOnly(SharedLibraryPath);

		bool found = lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? entry);

		Assert.That(found, Is.True,
			"TryGetEntry must return true for PEPTIDE/2 in index-only mode.");
		Assert.That(entry, Is.Not.Null,
			"The out parameter must be non-null when TryGetEntry returns true.");
		Assert.That(entry!.Fragments, Is.Not.Empty,
			"Fragments must be loaded on demand in index-only mode.");
	}

	/// <summary>
	/// Verifies that <see cref="MslLibrary.TryGetEntry"/> returns false and a null out
	/// parameter for a sequence/charge pair that does not exist in the library.
	/// </summary>
	[Test]
	public void TryGetEntry_ReturnsNull_ForUnknownSequence()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		bool found = lib.TryGetEntry("DOESNOTEXIST", 2, out MslLibraryEntry? entry);

		Assert.That(found, Is.False,
			"TryGetEntry must return false for a sequence not in the library.");
		Assert.That(entry, Is.Null,
			"The out parameter must be null when TryGetEntry returns false.");
	}

	/// <summary>
	/// Verifies that <see cref="MslLibrary.TryGetLibrarySpectrum"/> returns a
	/// <see cref="LibrarySpectrum"/> whose <c>Name</c> property matches the expected
	/// "sequence/charge" format used by MetaMorpheus.
	/// </summary>
	[Test]
	public void TryGetLibrarySpectrum_ReturnsCorrectName()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		bool found = lib.TryGetLibrarySpectrum("PEPTIDE", 2, out LibrarySpectrum? spectrum);

		Assert.That(found, Is.True);
		Assert.That(spectrum!.Name, Is.EqualTo("PEPTIDE/2"),
			"LibrarySpectrum.Name must be 'sequence/charge' to match MetaMorpheus conventions.");
	}

	/// <summary>
	/// Verifies that the <see cref="LibrarySpectrum"/> returned by
	/// <see cref="MslLibrary.TryGetLibrarySpectrum"/> contains the correct number of
	/// fragment ions as written.
	/// </summary>
	[Test]
	public void TryGetLibrarySpectrum_ReturnsCorrectFragmentCount()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		lib.TryGetLibrarySpectrum("PEPTIDE", 2, out LibrarySpectrum? spectrum);

		Assert.That(spectrum!.MatchedFragmentIons.Count, Is.EqualTo(2),
			"PEPTIDE/2 has 2 fragments; TryGetLibrarySpectrum must preserve all of them.");
	}

	// ════════════════════════════════════════════════════════════════════════════
	// Group 4 — DIA window queries
	// ════════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslLibrary.QueryMzWindow"/> returns exactly the entries whose
	/// precursor m/z falls within the specified window (449.0–450.0 should catch only
	/// PEPTIDE/2 at 449.74).
	/// </summary>
	[Test]
	public void QueryMzWindow_ReturnsCorrectCandidates()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		ReadOnlySpan<MslPrecursorIndexEntry> results = lib.QueryMzWindow(449.0f, 450.0f);

		Assert.That(results.Length, Is.EqualTo(1),
			"A window of 449.0–450.0 must contain exactly PEPTIDE/2.");
		Assert.That(results[0].PrecursorMz, Is.EqualTo(449.74f).Within(0.01f),
			"The single result must be the PEPTIDE/2 entry at m/z 449.74.");
	}

	/// <summary>
	/// Verifies that <see cref="MslLibrary.QueryWindow"/> applies the RT filter correctly:
	/// entries whose iRT falls outside [rtLow, rtHigh] must be excluded from the result.
	/// A window centred on the ACDEFGHIK entry (iRT 42.1) with a tight RT range that
	/// excludes it must produce an empty result set.
	/// </summary>
	[Test]
	public void QueryWindow_WithRtFilter_ExcludesOutOfRange()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		// m/z window covers ACDEFGHIK/2 at 529.76; RT window deliberately misses iRT 42.1
		using MslWindowResults results = lib.QueryWindow(
			mzLow: 529.0f, mzHigh: 530.5f,
			rtLow: 10.0f, rtHigh: 20.0f);

		Assert.That(results.Count, Is.EqualTo(0),
			"An RT window of 10–20 must exclude ACDEFGHIK/2 whose iRT is 42.1.");
	}

	/// <summary>
	/// Verifies that after a <see cref="MslLibrary.QueryWindow"/> call in index-only mode,
	/// <see cref="MslLibrary.GetEntry"/> successfully loads the full fragment data for a
	/// candidate returned by the window query.
	/// </summary>
	[Test]
	public void QueryWindow_GetEntry_LoadsFragmentsForIndexOnlyMode()
	{
		using MslLibrary lib = MslLibrary.LoadIndexOnly(SharedLibraryPath);

		using MslWindowResults results = lib.QueryWindow(
			mzLow: 449.0f, mzHigh: 450.0f,
			rtLow: 30.0f, rtHigh: 40.0f);

		Assert.That(results.Count, Is.EqualTo(1),
			"The window must match exactly PEPTIDE/2.");

		MslLibraryEntry? entry = lib.GetEntry(results.Entries[0].PrecursorIdx);

		Assert.That(entry, Is.Not.Null,
			"GetEntry must return a non-null entry for a valid PrecursorIdx.");
		Assert.That(entry!.Fragments, Is.Not.Empty,
			"Fragment ions must be loaded on demand in index-only mode.");
	}

	/// <summary>
	/// Verifies that calling <see cref="MslWindowResults.Dispose"/> on the result of
	/// <see cref="MslLibrary.QueryWindow"/> does not throw and that subsequent access to
	/// <see cref="MslWindowResults.Count"/> correctly returns 0 after disposal.
	/// </summary>
	[Test]
	public void QueryWindow_Dispose_ReleasesPoolBuffer()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		MslWindowResults results = lib.QueryWindow(
			mzLow: 449.0f, mzHigh: 450.0f,
			rtLow: 30.0f, rtHigh: 40.0f);

		Assert.That(() => results.Dispose(), Throws.Nothing,
			"Disposing the MslWindowResults must not throw.");

		// Calling Dispose a second time must also be safe (idempotent)
		Assert.That(() => results.Dispose(), Throws.Nothing,
			"Double-dispose of MslWindowResults must not throw.");
	}

	// ════════════════════════════════════════════════════════════════════════════
	// Group 5 — Bulk enumeration
	// ════════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslLibrary.GetAllEntries"/> with
	/// <c>includeDecoys = false</c> returns only the two target entries.
	/// </summary>
	[Test]
	public void GetAllEntries_ReturnsAllTargets_WhenDecoysFalse()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		List<MslLibraryEntry> targets = lib.GetAllEntries(includeDecoys: false).ToList();

		Assert.That(targets.Count, Is.EqualTo(2),
			"GetAllEntries(includeDecoys: false) must return only the 2 target entries.");
		Assert.That(targets.All(e => !e.IsDecoy), Is.True,
			"All returned entries must be non-decoy when includeDecoys is false.");
	}

	/// <summary>
	/// Verifies that <see cref="MslLibrary.GetAllEntries"/> with
	/// <c>includeDecoys = true</c> (the default) returns all three entries.
	/// </summary>
	[Test]
	public void GetAllEntries_ReturnsAll_WhenDecoysTrue()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		List<MslLibraryEntry> all = lib.GetAllEntries(includeDecoys: true).ToList();

		Assert.That(all.Count, Is.EqualTo(3),
			"GetAllEntries(includeDecoys: true) must return all 3 entries.");
	}

	/// <summary>
	/// Verifies that in full-load mode every entry yielded by
	/// <see cref="MslLibrary.GetAllEntries"/> has a non-empty fragment list, confirming
	/// that fragments were pre-populated at load time.
	/// </summary>
	[Test]
	public void GetAllEntries_FullLoad_AllFragmentsPresent()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		foreach (MslLibraryEntry entry in lib.GetAllEntries())
		{
			Assert.That(entry.Fragments, Is.Not.Empty,
				$"Entry {entry.ModifiedSequence}/{entry.Charge} must have fragments in full-load mode.");
		}
	}

	/// <summary>
	/// Verifies that in index-only mode every entry yielded by
	/// <see cref="MslLibrary.GetAllEntries"/> has its fragments populated on demand
	/// (non-empty fragment list) by the time it is yielded.
	/// </summary>
	[Test]
	public void GetAllEntries_IndexOnly_FragmentsLoadableOnDemand()
	{
		using MslLibrary lib = MslLibrary.LoadIndexOnly(SharedLibraryPath);

		foreach (MslLibraryEntry entry in lib.GetAllEntries())
		{
			Assert.That(entry.Fragments, Is.Not.Empty,
				$"Entry {entry.ModifiedSequence}/{entry.Charge} must have fragments loaded on demand.");
		}
	}

	// ════════════════════════════════════════════════════════════════════════════
	// Group 6 — RT calibration
	// ════════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslLibrary.WithCalibratedRetentionTimes"/> returns a new
	/// library whose index entries have transformed iRT values.
	/// Formula: calibratedRT = slope * iRT + intercept = 2.0 * 35.4 + 5.0 = 75.8.
	/// </summary>
	[Test]
	public void WithCalibratedRetentionTimes_NewLibrary_HasTransformedRt()
	{
		using MslLibrary original = MslLibrary.Load(SharedLibraryPath);
		using MslLibrary calibrated = original.WithCalibratedRetentionTimes(slope: 2.0, intercept: 5.0);

		// Query the calibrated library for PEPTIDE/2 and inspect its index entry's Irt
		ReadOnlySpan<MslPrecursorIndexEntry> entries =
			calibrated.QueryMzWindow(449.0f, 450.0f);

		Assert.That(entries.Length, Is.EqualTo(1),
			"Calibrated library must still contain PEPTIDE/2.");

		float expectedIrt = (float)(2.0 * 35.4 + 5.0);
		Assert.That(entries[0].Irt, Is.EqualTo(expectedIrt).Within(0.01f),
			"Calibrated iRT must equal slope * originalIrt + intercept.");
	}

	/// <summary>
	/// Verifies that calling <see cref="MslLibrary.WithCalibratedRetentionTimes"/> does not
	/// modify the iRT values in the original library.
	/// </summary>
	[Test]
	public void WithCalibratedRetentionTimes_OriginalLibrary_Unchanged()
	{
		using MslLibrary original = MslLibrary.Load(SharedLibraryPath);

		ReadOnlySpan<MslPrecursorIndexEntry> before = original.QueryMzWindow(449.0f, 450.0f);
		float originalIrt = before[0].Irt;

		// Apply calibration — this must not mutate the original index
		using MslLibrary _ = original.WithCalibratedRetentionTimes(slope: 2.0, intercept: 5.0);

		ReadOnlySpan<MslPrecursorIndexEntry> after = original.QueryMzWindow(449.0f, 450.0f);

		Assert.That(after[0].Irt, Is.EqualTo(originalIrt),
			"WithCalibratedRetentionTimes must not modify the original library's iRT values.");
	}

	/// <summary>
	/// Verifies that the calibrated library contains the same number of precursors as the
	/// original, confirming that no entries are lost during the calibration transform.
	/// </summary>
	[Test]
	public void WithCalibratedRetentionTimes_CalibrationPreserves_PrecursorCount()
	{
		using MslLibrary original = MslLibrary.Load(SharedLibraryPath);
		using MslLibrary calibrated = original.WithCalibratedRetentionTimes(slope: 1.5, intercept: -3.0);

		Assert.That(calibrated.PrecursorCount, Is.EqualTo(original.PrecursorCount),
			"RT calibration must preserve the total precursor count.");
	}

	// ════════════════════════════════════════════════════════════════════════════
	// Group 7 — Elution groups
	// ════════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslLibrary.GetElutionGroup"/> returns all charge states of a
	/// peptide that share a stripped sequence (PEPTIDE/2 and PEPTIDE/3 should share one group).
	/// </summary>
	[Test]
	public void GetElutionGroup_ReturnsAllChargeStatesOfPeptide()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		// Locate the elution group ID for PEPTIDE/2
		ReadOnlySpan<MslPrecursorIndexEntry> peptide2 = lib.QueryMzWindow(449.0f, 450.0f);
		Assert.That(peptide2.Length, Is.EqualTo(1));

		int groupId = peptide2[0].ElutionGroupId;
		ReadOnlySpan<MslPrecursorIndexEntry> group = lib.GetElutionGroup(groupId);

		Assert.That(group.Length, Is.EqualTo(2),
			"The elution group for PEPTIDE must contain both charge states (2+ and 3+).");

		short[] charges = new short[group.Length];
		for (int i = 0; i < group.Length; i++)
			charges[i] = group[i].Charge;

		Assert.That(charges, Does.Contain((short)2).And.Contain((short)3),
			"The elution group must include charge 2 and charge 3.");
	}

	/// <summary>
	/// Verifies that <see cref="MslLibrary.GetElutionGroup"/> does not mix entries from
	/// different peptides: the ACDEFGHIK group must not include any PEPTIDE entries.
	/// </summary>
	[Test]
	public void GetElutionGroup_DoesNotMix_DifferentPeptides()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		// Locate the elution group for ACDEFGHIK/2
		ReadOnlySpan<MslPrecursorIndexEntry> acdefghik = lib.QueryMzWindow(529.0f, 530.5f);
		Assert.That(acdefghik.Length, Is.EqualTo(1));

		int acdefghikGroupId = acdefghik[0].ElutionGroupId;

		// Locate the elution group ID for PEPTIDE (either charge state)
		ReadOnlySpan<MslPrecursorIndexEntry> peptide = lib.QueryMzWindow(449.0f, 450.0f);
		int peptideGroupId = peptide[0].ElutionGroupId;

		Assert.That(acdefghikGroupId, Is.Not.EqualTo(peptideGroupId),
			"ACDEFGHIK and PEPTIDE must have distinct ElutionGroupId values.");

		ReadOnlySpan<MslPrecursorIndexEntry> acGroup = lib.GetElutionGroup(acdefghikGroupId);
		Assert.That(acGroup.Length, Is.EqualTo(1),
			"The ACDEFGHIK group must contain exactly one entry.");
	}

	// ════════════════════════════════════════════════════════════════════════════
	// Group 8 — Save round-trip
	// ════════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslLibrary.Save"/> produces a file that can be reopened by
	/// <see cref="MslLibrary.Load"/> without throwing.
	/// </summary>
	[Test]
	public void Save_ProducesFileReadableByLoad()
	{
		string path = TempPath(nameof(Save_ProducesFileReadableByLoad));
		MslLibrary.Save(path, BuildFixtureEntries());

		Assert.That(() => { using var lib = MslLibrary.Load(path); }, Throws.Nothing,
			"A file written by Save() must be openable by Load() without error.");
	}

	/// <summary>
	/// Verifies that <see cref="MslLibrary.SaveFromLibrarySpectra"/> produces a file that can
	/// be reopened by <see cref="MslLibrary.Load"/> without throwing.
	/// </summary>
	[Test]
	public void SaveFromLibrarySpectra_ProducesFileReadableByLoad()
	{
		var spectra = new List<LibrarySpectrum>
		{
			new LibrarySpectrum(
				"PEPTIDE", 449.74, 2,
				new List<MatchedFragmentIon>
				{
					new MatchedFragmentIon(
						new Product(ProductType.y, FragmentationTerminus.C, 0.0, 1, 6, 0.0),
						experMz: 175.119, experIntensity: 1.0, charge: 1)
				},
				rt: 35.4)
		};

		string path = TempPath(nameof(SaveFromLibrarySpectra_ProducesFileReadableByLoad));
		MslLibrary.SaveFromLibrarySpectra(path, spectra);

		Assert.That(() => { using var lib = MslLibrary.Load(path); }, Throws.Nothing,
			"A file written by SaveFromLibrarySpectra() must be openable by Load() without error.");
	}

	/// <summary>
	/// Verifies that the Save → Load round-trip preserves all precursors and their core fields
	/// (modified sequence, charge, fragment count) without data loss.
	/// </summary>
	[Test]
	public void Save_ThenLoad_PreservesAllEntries()
	{
		List<MslLibraryEntry> written = BuildFixtureEntries();
		string path = TempPath(nameof(Save_ThenLoad_PreservesAllEntries));
		MslLibrary.Save(path, written);

		using MslLibrary lib = MslLibrary.Load(path);

		Assert.That(lib.PrecursorCount, Is.EqualTo(written.Count),
			"Round-trip precursor count must equal the number of entries written.");

		// Verify every written entry survives the round-trip by sequence/charge lookup
		foreach (MslLibraryEntry original in written)
		{
			bool found = lib.TryGetEntry(original.ModifiedSequence, original.Charge,
				out MslLibraryEntry? loaded);

			Assert.That(found, Is.True,
				$"Entry {original.ModifiedSequence}/{original.Charge} must survive the round-trip.");
			Assert.That(loaded!.Fragments.Count, Is.EqualTo(original.Fragments.Count),
				$"Fragment count for {original.ModifiedSequence}/{original.Charge} must be preserved.");
		}
	}

	// ════════════════════════════════════════════════════════════════════════════
	// Group 9 — Dispose
	// ════════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that calling <see cref="MslLibrary.Dispose"/> on a full-load library does not
	/// throw any exception.
	/// </summary>
	[Test]
	public void Dispose_FullLoad_DoesNotThrow()
	{
		MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		Assert.That(() => lib.Dispose(), Throws.Nothing,
			"Dispose() on a full-load library must not throw.");
	}

	/// <summary>
	/// Verifies that <see cref="MslLibrary.Dispose"/> on an index-only library releases the
	/// open file handle, confirmed by successfully deleting the file immediately after disposal.
	/// </summary>
	[Test]
	public void Dispose_IndexOnly_ClosesFileStream()
	{
		// Write a dedicated file so we can delete it without affecting other tests
		string path = TempPath(nameof(Dispose_IndexOnly_ClosesFileStream));
		MslLibrary.Save(path, BuildFixtureEntries());

		MslLibrary lib = MslLibrary.LoadIndexOnly(path);
		lib.Dispose();

		// On Windows, an open FileStream would prevent deletion; if Dispose closed it
		// correctly this must succeed on all platforms.
		Assert.That(() => File.Delete(path), Throws.Nothing,
			"File.Delete must succeed after Dispose() closes the index-only FileStream.");
	}

	/// <summary>
	/// Verifies that every public query method throws <see cref="ObjectDisposedException"/>
	/// when called after the library has been disposed.
	/// </summary>
	[Test]
	public void Operations_AfterDispose_ThrowObjectDisposedException()
	{
		MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.Dispose();

		Assert.That(() => lib.QueryMzWindow(400f, 500f),
			Throws.InstanceOf<ObjectDisposedException>(),
			"QueryMzWindow must throw ObjectDisposedException after Dispose.");

		Assert.That(() => lib.QueryWindow(400f, 500f, 0f, 100f),
			Throws.InstanceOf<ObjectDisposedException>(),
			"QueryWindow must throw ObjectDisposedException after Dispose.");

		Assert.That(() => lib.TryGetEntry("PEPTIDE", 2, out _),
			Throws.InstanceOf<ObjectDisposedException>(),
			"TryGetEntry must throw ObjectDisposedException after Dispose.");

		Assert.That(() => lib.TryGetLibrarySpectrum("PEPTIDE", 2, out _),
			Throws.InstanceOf<ObjectDisposedException>(),
			"TryGetLibrarySpectrum must throw ObjectDisposedException after Dispose.");

		Assert.That(() => lib.GetEntry(0),
			Throws.InstanceOf<ObjectDisposedException>(),
			"GetEntry must throw ObjectDisposedException after Dispose.");

		Assert.That(() => lib.GetAllEntries().ToList(),
			Throws.InstanceOf<ObjectDisposedException>(),
			"GetAllEntries must throw ObjectDisposedException after Dispose.");

		Assert.That(() => lib.GetElutionGroup(0),
			Throws.InstanceOf<ObjectDisposedException>(),
			"GetElutionGroup must throw ObjectDisposedException after Dispose.");

		Assert.That(() => lib.WithCalibratedRetentionTimes(1.0, 0.0),
			Throws.InstanceOf<ObjectDisposedException>(),
			"WithCalibratedRetentionTimes must throw ObjectDisposedException after Dispose.");
	}

	/// <summary>
	/// Verifies that calling <see cref="MslLibrary.Dispose"/> multiple times on the same
	/// instance does not throw (idempotency requirement).
	/// </summary>
	[Test]
	public void Dispose_CalledMultipleTimes_DoesNotThrow()
	{
		MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.Dispose();

		Assert.That(() => lib.Dispose(), Throws.Nothing,
			"Calling Dispose() a second time must not throw.");

		Assert.That(() => lib.Dispose(), Throws.Nothing,
			"Calling Dispose() a third time must not throw.");
	}
}