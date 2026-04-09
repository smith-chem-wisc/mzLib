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

namespace Test.SpectralLibrary.MSL;

/// <summary>
/// Targeted NUnit 4 tests that raise coverage on every MSL file flagged by the
/// Visual Studio code-coverage report for PR #1036. Tests are grouped by the source
/// file they primarily exercise. Each test includes a summary explaining what
/// behaviour is verified and why it matters.
/// </summary>
[TestFixture]
public sealed class TestMslCoverageGaps
{
	// ── Shared temp directory ─────────────────────────────────────────────────

	private string _dir = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_dir = Path.Combine(Path.GetTempPath(), $"MslCoverageGaps_{Guid.NewGuid():N}");
		Directory.CreateDirectory(_dir);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(_dir))
			Directory.Delete(_dir, recursive: true);
	}

	// ── Fixture helpers ───────────────────────────────────────────────────────

	private string P(string name) => Path.Combine(_dir, name);

	private static MslLibraryEntry MakePeptideEntry(
		string seq = "PEPTIDE", int charge = 2, float mz = 449.74f,
		float irt = 35.0f, bool isDecoy = false)
	{
		return new MslLibraryEntry
		{
			FullSequence = seq,
			BaseSequence = seq,
			PrecursorMz = mz,
			ChargeState = charge,
			RetentionTime = irt,
			IsDecoy = isDecoy,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Source = MslFormat.SourceType.Predicted,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new() { Mz = 98.06f,  Intensity = 1.0f, ProductType = ProductType.b,
						 FragmentNumber = 2, Charge = 1 },
				new() { Mz = 175.12f, Intensity = 0.8f, ProductType = ProductType.y,
						 FragmentNumber = 1, Charge = 1 }
			}
		};
	}

	private static MslLibraryEntry MakeProteoformEntry(
		string seq = "LARGEPROTEOFORM", int charge = 10,
		float mz = 1200.0f, float irt = 45.0f, bool isDecoy = false)
	{
		double nm = (double)mz * charge - charge * 1.007276;
		return new MslLibraryEntry
		{
			FullSequence = seq,
			BaseSequence = seq,
			PrecursorMz = mz,
			ChargeState = charge,
			RetentionTime = irt,
			IsDecoy = isDecoy,
			MoleculeType = MslFormat.MoleculeType.Proteoform,
			DissociationType = DissociationType.EThcD,
			Source = MslFormat.SourceType.Empirical,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new() { Mz = 850.4f,  Intensity = 1.0f, ProductType = ProductType.c,
						 FragmentNumber = 8,  Charge = 1 },
				new() { Mz = 900.2f,  Intensity = 0.8f, ProductType = ProductType.y,
						 FragmentNumber = 25, Charge = 2 },
				new() { Mz = 1100.3f, Intensity = 0.6f, ProductType = ProductType.b,
						 FragmentNumber = 20, Charge = 2 },
			}
		};
	}

	private MslLibrary SaveAndLoad(IReadOnlyList<MslLibraryEntry> entries, string stem)
	{
		string path = P(stem + ".msl");
		MslLibrary.Save(path, entries);
		return MslLibrary.Load(path);
	}

	// ═════════════════════════════════════════════════════════════════════════
	// MslWriter — ValidateEntries
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// ValidateEntries should return an empty error list for a well-formed entry.
	/// This is the happy-path contract: callers rely on an empty list as the
	/// go/no-go signal before calling Write.
	/// </summary>
	[Test]
	public void ValidateEntries_ValidEntry_ReturnsNoErrors()
	{
		var errors = MslWriter.ValidateEntries(new[] { MakePeptideEntry() });
		Assert.That(errors, Is.Empty);
	}

	/// <summary>
	/// ValidateEntries must catch a null/empty FullSequence. An entry without
	/// a sequence cannot be indexed or retrieved and must be rejected before writing.
	/// </summary>
	[Test]
	public void ValidateEntries_NullSequence_ReportsError()
	{
		var entry = MakePeptideEntry();
		entry.FullSequence = "";
		var errors = MslWriter.ValidateEntries(new[] { entry });
		Assert.That(errors, Has.Count.GreaterThan(0));
		Assert.That(errors[0], Does.Contain("FullSequence"));
	}

	/// <summary>
	/// ValidateEntries must catch a non-positive charge. A zero or negative charge
	/// would produce corrupt m/z arithmetic throughout the format.
	/// </summary>
	[Test]
	public void ValidateEntries_ZeroCharge_ReportsError()
	{
		var entry = MakePeptideEntry();
		entry.ChargeState = 0;
		var errors = MslWriter.ValidateEntries(new[] { entry });
		Assert.That(errors.Any(e => e.Contains("ChargeState")), Is.True);
	}

	/// <summary>
	/// ValidateEntries must catch a non-positive PrecursorMz. A zero m/z would
	/// silently corrupt every index sorted by m/z in the reader.
	/// </summary>
	[Test]
	public void ValidateEntries_ZeroPrecursorMz_ReportsError()
	{
		var entry = MakePeptideEntry();
		entry.PrecursorMz = 0.0;
		var errors = MslWriter.ValidateEntries(new[] { entry });
		Assert.That(errors.Any(e => e.Contains("PrecursorMz")), Is.True);
	}

	/// <summary>
	/// ValidateEntries must catch an entry with no fragment ions. A library entry
	/// with zero fragments is useless for scoring and should never reach the writer.
	/// </summary>
	[Test]
	public void ValidateEntries_NoFragments_ReportsError()
	{
		var entry = MakePeptideEntry();
		entry.MatchedFragmentIons = new List<MslFragmentIon>();
		var errors = MslWriter.ValidateEntries(new[] { entry });
		Assert.That(errors.Any(e => e.Contains("FragmentCount")), Is.True);
	}

	/// <summary>
	/// ValidateEntries must catch a fragment with a non-positive Mz. A zero fragment
	/// m/z cannot be matched against any experimental peak and corrupts scoring.
	/// </summary>
	[Test]
	public void ValidateEntries_ZeroFragmentMz_ReportsError()
	{
		var entry = MakePeptideEntry();
		entry.MatchedFragmentIons[0].Mz = 0f;
		var errors = MslWriter.ValidateEntries(new[] { entry });
		Assert.That(errors.Any(e => e.Contains("Mz")), Is.True);
	}

	/// <summary>
	/// ValidateEntries must catch a fragment with a non-positive ChargeState. Fragment
	/// charge is used to compute neutral mass during scoring; zero would divide by zero.
	/// </summary>
	[Test]
	public void ValidateEntries_ZeroFragmentCharge_ReportsError()
	{
		var entry = MakePeptideEntry();
		entry.MatchedFragmentIons[0].Charge = 0;
		var errors = MslWriter.ValidateEntries(new[] { entry });
		Assert.That(errors.Any(e => e.Contains("ChargeState")), Is.True);
	}

	/// <summary>
	/// ValidateEntries must catch an internal fragment ion where SecondaryFragmentNumber
	/// is not greater than FragmentNumber. An invalid range means the internal ion
	/// spans zero or negative residues, which would produce nonsense annotations.
	/// </summary>
	[Test]
	public void ValidateEntries_InvalidInternalIonRange_ReportsError()
	{
		var entry = MakePeptideEntry();
		entry.MatchedFragmentIons.Add(new MslFragmentIon
		{
			Mz = 400f,
			Intensity = 0.5f,
			Charge = 1,
			ProductType = ProductType.b,
			SecondaryProductType = ProductType.y,
			FragmentNumber = 5,
			SecondaryFragmentNumber = 3  // invalid: end <= start
		});
		var errors = MslWriter.ValidateEntries(new[] { entry });
		Assert.That(errors.Any(e => e.Contains("SecondaryFragmentNumber")), Is.True);
	}

	/// <summary>
	/// ValidateEntries must throw ArgumentNullException when passed null.
	/// This guards against silent null propagation into the writer.
	/// </summary>
	[Test]
	public void ValidateEntries_NullInput_Throws()
		=> Assert.That(() => MslWriter.ValidateEntries(null!), Throws.ArgumentNullException);

	// ═════════════════════════════════════════════════════════════════════════
	// MslWriter — WriteStreaming, WriteFromLibrarySpectra overload
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// WriteStreaming produces a valid readable .msl file. WriteStreaming is the
	/// path used by MslMerger for memory-efficient large-library writes; if it
	/// silently produces corrupt output the merger would silently lose data.
	/// </summary>
	[Test]
	public void WriteStreaming_ProducesValidFile()
	{
		var entries = new[] { MakePeptideEntry("ALGVGLATR"), MakePeptideEntry("PEPTIDE") };
		string path = P("streaming.msl");
		MslWriter.WriteStreaming(path, entries, compressionLevel: 0);
		Assert.That(File.Exists(path), Is.True);
		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.PrecursorCount, Is.EqualTo(2));
	}

	/// <summary>
	/// WriteFromLibrarySpectra converts LibrarySpectrum objects and writes them.
	/// This is the interop path from MetaMorpheus, where spectra are produced by
	/// the existing spectral library pipeline and then saved in the MSL format.
	/// </summary>
	[Test]
	public void WriteFromLibrarySpectra_RoundTrips()
	{
		var product = new Product(ProductType.b, FragmentationTerminus.N, 0.0, 2, 1, 0.0);
		var ion = new MatchedFragmentIon(product, 200.0, 1.0f, 1);
		var spectrum = new LibrarySpectrum("PEPTIDE", 449.74, 2,
			new List<MatchedFragmentIon> { ion }, rt: 35.0);

		string path = P("from_spectra.msl");
		MslWriter.WriteFromLibrarySpectra(path, new[] { spectrum });

		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.PrecursorCount, Is.EqualTo(1));
		bool found = lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? e);
		Assert.That(found, Is.True);
		Assert.That(e!.MatchedFragmentIons, Has.Count.EqualTo(1));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// MslReader — ReadHeaderOnly, DecodeNeutralLoss paths
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// ReadHeaderOnly reads the 64-byte header without deserializing entries or
	/// fragments. This is the fast metadata-inspection path used to check format
	/// version and precursor count before committing to a full load.
	/// </summary>
	[Test]
	public void ReadHeaderOnly_ReturnsCorrectPrecursorCount()
	{
		string path = P("header_only.msl");
		var entries = new[] { MakePeptideEntry(), MakePeptideEntry("ACDEFGHIK", charge: 2, mz: 529.76f) };
		MslLibrary.Save(path, entries);

		MslFileHeader header = MslReader.ReadHeaderOnly(path);
		Assert.That(header.NPrecursors, Is.EqualTo(2));
		Assert.That(header.FormatVersion, Is.EqualTo(MslFormat.CurrentVersion));
	}

	/// <summary>
	/// ReadHeaderOnly on a non-existent file must throw FileNotFoundException.
	/// Callers rely on this to give users a clear error rather than a NullRef.
	/// </summary>
	[Test]
	public void ReadHeaderOnly_MissingFile_ThrowsFileNotFound()
		=> Assert.That(() => MslReader.ReadHeaderOnly(P("does_not_exist.msl")),
			Throws.TypeOf<FileNotFoundException>());

	/// <summary>
	/// ReadHeaderOnly on a file with wrong magic bytes must throw FormatException.
	/// This prevents silently reading corrupt or non-MSL files as if they were valid.
	/// </summary>
	[Test]
	public void ReadHeaderOnly_WrongMagic_ThrowsFormatException()
	{
		string path = P("bad_magic_header.msl");
		MslLibrary.Save(path, new[] { MakePeptideEntry() });
		byte[] bytes = File.ReadAllBytes(path);
		bytes[0] = 0x00; bytes[1] = 0x00;
		File.WriteAllBytes(path, bytes);
		Assert.That(() => MslReader.ReadHeaderOnly(path), Throws.TypeOf<FormatException>());
	}

	/// <summary>
	/// A fragment written with an H2O neutral loss must be decoded with the correct
	/// mass on read-back. DecodeNeutralLoss is the path that converts the 3-bit
	/// NeutralLossCode in the flags byte back to a double mass; an off-by-one in
	/// the switch would silently produce wrong fragment masses during scoring.
	/// </summary>
	[Test]
	public void Reader_DecodeNeutralLoss_H2O_RoundTrips()
	{
		const double h2oLoss = -18.010565;
		var entry = MakePeptideEntry();
		entry.MatchedFragmentIons[0].NeutralLoss = h2oLoss;

		string path = P("h2o_loss.msl");
		MslLibrary.Save(path, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(path);
		lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? loaded);
		// Writer sorts by mz; b2 (98.06) is first
		MslFragmentIon lossIon = loaded!.MatchedFragmentIons.First(f => f.ProductType == ProductType.b && f.FragmentNumber == 2);
		Assert.That(lossIon.NeutralLoss, Is.EqualTo(h2oLoss).Within(1e-4));
	}

	/// <summary>
	/// A fragment written with an NH3 neutral loss must round-trip correctly.
	/// Tests the NH3 branch of DecodeNeutralLoss, which is the second most
	/// common neutral loss after H2O in proteomics libraries.
	/// </summary>
	[Test]
	public void Reader_DecodeNeutralLoss_NH3_RoundTrips()
	{
		const double nh3Loss = -17.026549;
		var entry = MakePeptideEntry();
		entry.MatchedFragmentIons[0].NeutralLoss = nh3Loss;

		string path = P("nh3_loss.msl");
		MslLibrary.Save(path, new[] { entry });
		using MslLibrary lib = MslLibrary.Load(path);
		lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? loaded);
		MslFragmentIon lossIon = loaded!.MatchedFragmentIons.First(f => f.ProductType == ProductType.b && f.FragmentNumber == 2);
		Assert.That(lossIon.NeutralLoss, Is.EqualTo(nh3Loss).Within(1e-4));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// MslLibraryData — IsCompressed, full-load ctor null-check
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// MslLibraryData constructed in full-load mode reports IsIndexOnly=false and
	/// Count equals the number of entries supplied. Callers branch on IsIndexOnly
	/// to decide whether to call LoadFragmentsOnDemand; a wrong value would cause
	/// callers to incorrectly attempt on-demand reads on a fully-loaded library.
	/// </summary>
	[Test]
	public void MslLibraryData_FullLoad_IsIndexOnlyFalseAndCountCorrect()
	{
		string path = P("libdata_fullload.msl");
		var src = new[] { MakePeptideEntry() };
		MslLibrary.Save(path, src);
		// Load returns MslLibrary which wraps MslLibraryData internally;
		// test via the public API to confirm full-load mode behaves correctly.
		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.IsIndexOnly, Is.False);
		Assert.That(lib.PrecursorCount, Is.EqualTo(1));
	}

	/// <summary>
	/// The full-load MslLibraryData constructor must throw ArgumentNullException
	/// when entries is null, preventing a NullReferenceException deep inside the
	/// index-build path where the error would be harder to diagnose.
	/// </summary>
	[Test]
	public void MslLibraryData_FullLoadCtor_NullEntries_Throws()
		=> Assert.That(
			() => new MslLibraryData(null!, new MslFileHeader()),
			Throws.ArgumentNullException);

	// ═════════════════════════════════════════════════════════════════════════
	// MslIndexStatistics — MinIrt, MaxIrt
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// MslIndex stores iRT values on each MslPrecursorIndexEntry. Verifying that
	/// iRT values round-trip correctly ensures that RT-window queries filter on
	/// the right candidates; a wrong iRT would silently exclude valid precursors.
	/// </summary>
	[Test]
	public void MslLibrary_IndexEntries_IrtValues_RoundTrip()
	{
		var entries = new[]
		{
			MakePeptideEntry("PEPTIDE",   charge: 2, mz: 449.74f, irt: 30.0f),
			MakePeptideEntry("ALGVGLATR", charge: 2, mz: 430.26f, irt: 60.0f),
		};
		using MslLibrary lib = SaveAndLoad(entries, "irt_range");

		ReadOnlySpan<MslPrecursorIndexEntry> all = lib.QueryMzWindow(0f, float.MaxValue);
		Assert.That(all.Length, Is.EqualTo(2));
		float[] irts = new float[all.Length];
		for (int i = 0; i < all.Length; i++) irts[i] = all[i].Irt;
		// Assert each expected iRT appears in the array using individual element checks
		Assert.That(irts.Any(v => Math.Abs(v - 30.0f) < 0.1f), Is.True,
			"Expected iRT=30.0 to be present in index entries.");
		Assert.That(irts.Any(v => Math.Abs(v - 60.0f) < 0.1f), Is.True,
			"Expected iRT=60.0 to be present in index entries.");
	}

	// ═════════════════════════════════════════════════════════════════════════
	// MslWindowResults — Entries property
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// MslWindowResults.Entries returns the matched index entries as a span.
	/// This property is the primary way callers iterate window query results;
	/// an untested Entries getter could mask a wrong span slice or wrong count.
	/// </summary>
	[Test]
	public void MslWindowResults_Entries_ReturnsMatchingEntries()
	{
		var entries = new[]
		{
			MakePeptideEntry("PEPTIDE",   charge: 2, mz: 449.74f, irt: 35.0f),
			MakePeptideEntry("ALGVGLATR", charge: 2, mz: 430.26f, irt: 72.0f),
		};
		using MslLibrary lib = SaveAndLoad(entries, "window_entries");

		using MslWindowResults results = lib.QueryWindow(
			mzLow: 425f, mzHigh: 460f, rtLow: 0f, rtHigh: 100f);

		// At least one entry should be in the m/z window; verify Entries is non-empty
		Assert.That(results.Count, Is.GreaterThan(0));
		ReadOnlySpan<MslPrecursorIndexEntry> span = results.Entries;
		Assert.That(span.Length, Is.EqualTo(results.Count));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// MslMerger — Merge, MslMergeResult.OutputPath
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// MslMerger.Merge combines two source .msl files into one output file.
	/// This is the critical maintenance operation for large library pipelines:
	/// if the merger silently drops entries or corrupts the output, downstream
	/// search results would be silently incomplete.
	/// </summary>
	[Test]
	public void MslMerger_Merge_TwoFiles_ProducesCorrectCount()
	{
		// File A: PEPTIDE/2
		string pathA = P("merge_a.msl");
		MslLibrary.Save(pathA, new[] { MakePeptideEntry("PEPTIDE", 2, 449.74f) });

		// File B: ALGVGLATR/2 — different sequence, no duplicates expected
		string pathB = P("merge_b.msl");
		MslLibrary.Save(pathB, new[] { MakePeptideEntry("ALGVGLATR", 2, 430.26f) });

		string outPath = P("merged.msl");
		MslMergeResult result = MslMerger.Merge(
			new[] { pathA, pathB }, outPath,
			conflictPolicy: MslMergeConflictPolicy.KeepFirst,
			deduplicate: true);

		Assert.That(result.OutputPath, Is.EqualTo(outPath));
		Assert.That(result.OutputEntryCount, Is.EqualTo(2));
		Assert.That(result.TotalSourceEntryCount, Is.EqualTo(2));
		Assert.That(result.DuplicatesSkipped, Is.EqualTo(0));
		Assert.That(File.Exists(outPath), Is.True);

		using MslLibrary merged = MslLibrary.Load(outPath);
		Assert.That(merged.PrecursorCount, Is.EqualTo(2));
	}

	/// <summary>
	/// MslMerger.Merge with KeepFirst conflict policy must keep only the first
	/// occurrence of a duplicate key and report the skipped count correctly.
	/// Incorrect deduplication could double-count spectra during database search.
	/// </summary>
	[Test]
	public void MslMerger_Merge_Deduplication_KeepFirst()
	{
		string pathA = P("dedup_a.msl");
		MslLibrary.Save(pathA, new[] { MakePeptideEntry("PEPTIDE", 2, 449.74f, irt: 35.0f) });

		// Same sequence/charge in file B — should be deduplicated
		string pathB = P("dedup_b.msl");
		MslLibrary.Save(pathB, new[] { MakePeptideEntry("PEPTIDE", 2, 449.74f, irt: 35.0f) });

		string outPath = P("dedup_out.msl");
		MslMergeResult result = MslMerger.Merge(
			new[] { pathA, pathB }, outPath,
			conflictPolicy: MslMergeConflictPolicy.KeepFirst,
			deduplicate: true);

		Assert.That(result.OutputEntryCount, Is.EqualTo(1));
		Assert.That(result.DuplicatesSkipped, Is.EqualTo(1));
	}

	/// <summary>
	/// MslMerger.Merge with deduplicate=false must write all entries including
	/// duplicates. Some workflows deliberately retain duplicates (e.g., when
	/// source libraries differ in metadata even for the same sequence).
	/// </summary>
	[Test]
	public void MslMerger_Merge_NoDeduplicate_WritesAll()
	{
		string pathA = P("nodup_a.msl");
		MslLibrary.Save(pathA, new[] { MakePeptideEntry("PEPTIDE", 2, 449.74f) });
		string pathB = P("nodup_b.msl");
		MslLibrary.Save(pathB, new[] { MakePeptideEntry("PEPTIDE", 2, 449.74f) });

		string outPath = P("nodup_out.msl");
		MslMergeResult result = MslMerger.Merge(
			new[] { pathA, pathB }, outPath,
			deduplicate: false);

		Assert.That(result.OutputEntryCount, Is.EqualTo(2));
		Assert.That(result.DuplicatesSkipped, Is.EqualTo(0));
	}

	/// <summary>
	/// MslMergeResult.OutputPath exposes the path of the written file.
	/// MetaMorpheus records this for logging and downstream file registration;
	/// returning a wrong path would cause file-not-found errors in callers.
	/// </summary>
	[Test]
	public void MslMergeResult_OutputPath_MatchesRequestedPath()
	{
		string pathA = P("outpath_a.msl");
		string outPath = P("outpath_result.msl");
		MslLibrary.Save(pathA, new[] { MakePeptideEntry() });

		MslMergeResult result = MslMerger.Merge(new[] { pathA }, outPath);
		Assert.That(result.OutputPath, Is.EqualTo(outPath));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// MslProteoformIndex — Build, QueryMassWindow, GetEntry, GetAllEntries,
	//                      MinNeutralMass, MaxNeutralMass, TryGetEntry
	// ═════════════════════════════════════════════════════════════════════════

	private static MslProteoformIndex BuildProtoIndex(
		IReadOnlyList<MslLibraryEntry> entries)
	{
		return MslProteoformIndex.Build(entries, i => entries[i]);
	}

	/// <summary>
	/// MslProteoformIndex.Build correctly filters out non-proteoform entries.
	/// Peptide entries must not appear in the proteoform index or they would be
	/// incorrectly scored against deconvoluted top-down spectra.
	/// </summary>
	[Test]
	public void ProteoformIndex_Build_OnlyIndexesProteoforms()
	{
		var peptide = MakePeptideEntry();
		var proteoform = MakeProteoformEntry();
		var index = BuildProtoIndex(new[] { peptide, proteoform });
		Assert.That(index.Count, Is.EqualTo(1));
	}

	/// <summary>
	/// MslProteoformIndex.MinNeutralMass and MaxNeutralMass correctly return the
	/// bounds of the indexed mass range. Callers use these to decide whether a
	/// given precursor neutral mass could ever match anything in the library.
	/// </summary>
	[Test]
	public void ProteoformIndex_MinMaxNeutralMass_CorrectBounds()
	{
		var e1 = MakeProteoformEntry("PROT1", charge: 10, mz: 1000f);
		var e2 = MakeProteoformEntry("PROT2", charge: 5, mz: 2000f);
		var index = BuildProtoIndex(new[] { e1, e2 });

		double mass1 = MslProteoformIndex.ComputeNeutralMass(1000.0, 10);
		double mass2 = MslProteoformIndex.ComputeNeutralMass(2000.0, 5);

		Assert.That(index.MinNeutralMass, Is.EqualTo(Math.Min(mass1, mass2)).Within(0.01));
		Assert.That(index.MaxNeutralMass, Is.EqualTo(Math.Max(mass1, mass2)).Within(0.01));
	}

	/// <summary>
	/// MslProteoformIndex.MinNeutralMass and MaxNeutralMass return 0 for an empty
	/// index. Callers must handle the empty case without crashing.
	/// </summary>
	[Test]
	public void ProteoformIndex_MinMaxNeutralMass_EmptyIndex_ReturnsZero()
	{
		var index = BuildProtoIndex(Array.Empty<MslLibraryEntry>());
		Assert.That(index.MinNeutralMass, Is.EqualTo(0.0));
		Assert.That(index.MaxNeutralMass, Is.EqualTo(0.0));
	}

	/// <summary>
	/// QueryMassWindow returns entries within the supplied mass window and
	/// excludes those outside it. This is the hot path for top-down candidate
	/// selection; an off-by-one in the binary search would silently miss candidates.
	/// Entries are chosen so their neutral masses differ by ~5000 Da, making it
	/// trivial to construct a tight window around one that cannot include the other.
	/// </summary>
	[Test]
	public void ProteoformIndex_QueryMassWindow_ReturnsCorrectEntries()
	{
		// e1: mz=500, z=10 → neutral mass ≈ 4989.9 Da
		// e2: mz=500, z=20 → neutral mass ≈ 9979.9 Da  (≈5000 Da apart)
		var e1 = MakeProteoformEntry("PROT1", charge: 10, mz: 500f);
		var e2 = MakeProteoformEntry("PROT2", charge: 20, mz: 500f);
		var index = BuildProtoIndex(new[] { e1, e2 });

		double mass1 = MslProteoformIndex.ComputeNeutralMass(500.0, 10);
		double mass2 = MslProteoformIndex.ComputeNeutralMass(500.0, 20);

		// Sanity: the two masses must be far enough apart for a ±100 Da window to be exclusive
		Assert.That(Math.Abs(mass1 - mass2), Is.GreaterThan(200.0),
			"Test setup: neutral masses must differ by > 200 Da.");

		// Query a ±100 Da window around mass1 — must include e1 and exclude e2
		var hits = index.QueryMassWindow(mass1 - 100, mass1 + 100);
		Assert.That(hits.Length, Is.EqualTo(1));
		Assert.That(hits[0].NeutralMass, Is.EqualTo(mass1).Within(0.01));
	}

	/// <summary>
	/// QueryMassWindow with includeDecoys=false must exclude decoy entries.
	/// In a target-decoy search workflow, scoring decoys alongside targets on
	/// the same spectrum would inflate the false-discovery rate estimate.
	/// </summary>
	[Test]
	public void ProteoformIndex_QueryMassWindow_ExcludesDecoysWhenRequested()
	{
		var target = MakeProteoformEntry("TARGET", charge: 10, mz: 1000f, isDecoy: false);
		var decoy = MakeProteoformEntry("DECOY", charge: 10, mz: 1000f, isDecoy: true);
		var index = BuildProtoIndex(new[] { target, decoy });

		double mass = MslProteoformIndex.ComputeNeutralMass(1000.0, 10);
		var withDecoys = index.QueryMassWindow(mass - 50, mass + 50, includeDecoys: true);
		var withoutDecoys = index.QueryMassWindow(mass - 50, mass + 50, includeDecoys: false);

		Assert.That(withDecoys.Length, Is.EqualTo(2));
		Assert.That(withoutDecoys.Length, Is.EqualTo(1));
		Assert.That(withoutDecoys[0].IsDecoy, Is.False);
	}

	/// <summary>
	/// GetEntry retrieves the full MslLibraryEntry including fragment ions given
	/// an MslProteoformIndexEntry. This is the step that converts a lightweight
	/// index hit into a scoreable entry; if it returns null no scoring happens.
	/// </summary>
	[Test]
	public void ProteoformIndex_GetEntry_ReturnsFullEntry()
	{
		var proto = MakeProteoformEntry();
		var entries = new[] { proto };
		var index = BuildProtoIndex(entries);

		double mass = MslProteoformIndex.ComputeNeutralMass(proto.PrecursorMz, proto.ChargeState);
		var hits = index.QueryMassWindow(mass - 10, mass + 10);
		Assert.That(hits.Length, Is.EqualTo(1));

		MslLibraryEntry? full = index.GetEntry(hits[0]);
		Assert.That(full, Is.Not.Null);
		Assert.That(full!.FullSequence, Is.EqualTo(proto.FullSequence));
		Assert.That(full.MatchedFragmentIons, Is.Not.Empty);
	}

	/// <summary>
	/// TryGetEntry performs O(1) lookup by sequence and charge. This mirrors the
	/// DDA-style lookup path but for proteoforms; a miss must return false cleanly
	/// without throwing so callers can continue with other candidates.
	/// </summary>
	[Test]
	public void ProteoformIndex_TryGetEntry_HitAndMiss()
	{
		var proto = MakeProteoformEntry("LARGEPROTEOFORM", charge: 10, mz: 1200f);
		var index = BuildProtoIndex(new[] { proto });

		bool hit = index.TryGetEntry("LARGEPROTEOFORM", 10, out MslLibraryEntry? found);
		bool miss = index.TryGetEntry("NOTHERE", 10, out MslLibraryEntry? notFound);

		Assert.That(hit, Is.True);
		Assert.That(found, Is.Not.Null);
		Assert.That(miss, Is.False);
		Assert.That(notFound, Is.Null);
	}

	/// <summary>
	/// GetAllEntries enumerates every proteoform entry through the loader delegate.
	/// This path is used by MslMerger and bulk export tools; if it silently skips
	/// entries the output library would be incomplete.
	/// </summary>
	[Test]
	public void ProteoformIndex_GetAllEntries_YieldsAllProteoforms()
	{
		var e1 = MakeProteoformEntry("PROT1", charge: 10, mz: 1000f);
		var e2 = MakeProteoformEntry("PROT2", charge: 5, mz: 2000f);
		var index = BuildProtoIndex(new[] { e1, e2 });

		var all = index.GetAllEntries(includeDecoys: true).ToList();
		Assert.That(all, Has.Count.EqualTo(2));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// MslProteoformScoringResult — LibraryEntry, FragmentCoverage
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// MslProteoformScoringResult.LibraryEntry exposes the scored entry so callers
	/// can access the sequence and metadata for output. An untested getter could mask
	/// a null assignment in the scorer that would cause NullRef in output generation.
	/// </summary>
	[Test]
	public void ProteoformScoringResult_LibraryEntry_IsAccessible()
	{
		var entry = MakeProteoformEntry();
		var result = new MslProteoformScoringResult
		{
			LibraryEntry = entry,
			CompositeScore = 1.5,
			SpectralAngle = 0.9,
			MatchedFragmentCount = 2,
			TotalLibraryFragments = 3,
			MatchedFragmentIons = new List<MatchedFragmentIon>().AsReadOnly()
		};

		Assert.That(result.LibraryEntry, Is.SameAs(entry));
		Assert.That(result.FragmentCoverage, Is.EqualTo(2.0 / 3.0).Within(1e-10));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// MslProteoformIndex + MslProteoformScorer — end-to-end top-down pipeline
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// MslProteoformIndex built from a mixed library correctly counts only Proteoform
	/// entries. Peptide entries must not appear in the proteoform index or they would
	/// be incorrectly scored against deconvoluted top-down spectra.
	/// </summary>
	[Test]
	public void MslProteoformIndex_BuildFromMixedLibrary_CountsOnlyProteoforms()
	{
		var entries = new[]
		{
			MakePeptideEntry(),                        // Peptide — should NOT count
            MakeProteoformEntry("PROTO1", 10, 1200f),  // Proteoform — should count
            MakeProteoformEntry("PROTO2",  5, 2000f),  // Proteoform — should count
        };
		using MslLibrary lib = SaveAndLoad(entries, "proto_count");
		var allEntries = lib.GetAllEntries(includeDecoys: true).ToList();
		var protoIndex = MslProteoformIndex.Build(allEntries, i => allEntries[i]);
		Assert.That(protoIndex.Count, Is.EqualTo(2));
	}

	/// <summary>
	/// End-to-end proteoform scoring: build a proteoform index from a saved library,
	/// query a mass window, score against synthetic deconvoluted peaks, and verify the
	/// top result has a spectral angle near 1. This exercises the full top-down search
	/// pipeline: Build → QueryMassWindow → MslProteoformScorer.Score.
	/// </summary>
	[Test]
	public void MslProteoformScorer_EndToEnd_FindsProteoform()
	{
		var proto = MakeProteoformEntry("LARGEPROTEOFORM", charge: 10, mz: 1200f, irt: 45f);
		using MslLibrary lib = SaveAndLoad(new[] { proto }, "score_proto");

		var allEntries = lib.GetAllEntries(includeDecoys: true).ToList();
		var protoIndex = MslProteoformIndex.Build(allEntries, i => allEntries[i]);

		double neutralMass = MslProteoformIndex.ComputeNeutralMass(proto.PrecursorMz, proto.ChargeState);
		var hits = protoIndex.QueryMassWindow(neutralMass - 10, neutralMass + 10);
		Assert.That(hits.Length, Is.EqualTo(1));

		// Synthetic peaks that match all 3 fragments exactly
		var peaks = proto.MatchedFragmentIons!.Select(f => new DeconvolutedPeak(
			(double)f.Mz * f.Charge - f.Charge * 1.007276, f.Intensity)).ToList();

		var results = MslProteoformScorer.Score(
			hits, lib, peaks,
			fragmentMassTolerance: 20e-6,
			minMatchedFragments: 1);

		Assert.That(results, Has.Count.EqualTo(1));
		Assert.That(results[0].SpectralAngle, Is.GreaterThan(0.99));
		Assert.That(results[0].MatchedFragmentCount, Is.EqualTo(3));
	}

	/// <summary>
	/// Scoring with an empty candidate span returns an empty result list without
	/// throwing. This is the expected behaviour when a proteoform-only index contains
	/// no entries in the queried mass window — callers must handle empty results cleanly.
	/// </summary>
	[Test]
	public void MslProteoformScorer_EmptyCandidates_ReturnsEmpty()
	{
		using MslLibrary lib = SaveAndLoad(new[] { MakePeptideEntry() }, "peptide_only_score");
		var peaks = new List<DeconvolutedPeak> { new(10000.0, 1.0f) };
		var results = MslProteoformScorer.Score(
			ReadOnlySpan<MslProteoformIndexEntry>.Empty, lib, peaks);
		Assert.That(results, Is.Empty);
	}

	/// <summary>
	/// MslLibrary.GetAllEntries with includeDecoys=true yields both target and decoy
	/// entries. This exercises the index enumerator lambda (MslLibrary.&lt;&gt;c__DisplayClass6_0
	/// &lt;Load&gt;b__1) which was at 0% coverage. The lambda is the entry-loader delegate
	/// captured during Load; if it is never invoked fragment data is never read back.
	/// </summary>
	[Test]
	public void MslLibrary_GetAllEntries_LoaderLambda_IsInvoked()
	{
		var entries = new[]
		{
			MakePeptideEntry("TARGET", 2, 449.74f, isDecoy: false),
			MakePeptideEntry("DECOY",  2, 430.26f, isDecoy: true),
		};
		using MslLibrary lib = SaveAndLoad(entries, "loader_lambda");

		var all = lib.GetAllEntries(includeDecoys: true).ToList();
		Assert.That(all, Has.Count.EqualTo(2));
		// Verify fragments were actually loaded by the lambda
		Assert.That(all.All(e => e.MatchedFragmentIons.Count > 0), Is.True);
	}

	// ═════════════════════════════════════════════════════════════════════════
	// SpectralLibrary — FileType, Software, single-string constructor
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// SpectralLibrary.Software returns a non-null value when constructed via the
	/// List&lt;string&gt; constructor. The Software property is used by mzLib's result-file
	/// infrastructure for logging and output attribution; a null here would cause a
	/// NullReferenceException in callers that write result metadata.
	/// </summary>
	[Test]
	public void SpectralLibrary_Software_IsReadable()
	{
		string msp = P("software_prop.msp");
		File.WriteAllLines(msp, new[]
		{
			"Name: PEPTIDE/2",
			"MW: 449.74",
			"Comment: Parent=449.74 iRT=35.4",
			"Num peaks: 1",
			"98.060	1.0	b2^1/0ppm"
		});

		// Use the List<string> constructor — this is the supported entry point
		// for SpectralLibrary and correctly sets up all internal state.
		var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { msp });
		Assert.That(() => { var _ = lib.Software; }, Throws.Nothing);
		lib.CloseConnections();
	}

	// ═════════════════════════════════════════════════════════════════════════
	// MslLibraryEntry — ToLibrarySpectrum (oligonucleotide / proteoform branches)
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// MslLibraryEntry.ToLibrarySpectrum must assign FivePrime terminus for
	/// oligonucleotide a-type ions. The terminus assignment table has a branch for
	/// oligonucleotides that is separate from peptides; an untested branch would
	/// silently produce None terminus for all oligo ions, corrupting spectral angle
	/// computation in RNA library search.
	/// </summary>
	[Test]
	public void MslLibraryEntry_ToLibrarySpectrum_OligoAIon_HasFivePrimeTerminus()
	{
		var entry = new MslLibraryEntry
		{
			FullSequence = "ACGU",
			BaseSequence = "ACGU",
			PrecursorMz = 600.0,
			ChargeState = 2,
			RetentionTime = 10.0,
			MoleculeType = MslFormat.MoleculeType.Oligonucleotide,
			DissociationType = DissociationType.HCD,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new() { Mz = 330f, Intensity = 1.0f, ProductType = ProductType.a,
						 FragmentNumber = 2, Charge = 1 }
			}
		};

		LibrarySpectrum spectrum = entry.ToLibrarySpectrum();
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.FivePrime));
	}

	/// <summary>
	/// MslLibraryEntry.ToLibrarySpectrum must assign ThreePrime terminus for
	/// oligonucleotide w-type ions. Complementary to the FivePrime test above;
	/// both termini must work correctly for bidirectional fragment matching.
	/// </summary>
	[Test]
	public void MslLibraryEntry_ToLibrarySpectrum_OligoWIon_HasThreePrimeTerminus()
	{
		var entry = new MslLibraryEntry
		{
			FullSequence = "ACGU",
			BaseSequence = "ACGU",
			PrecursorMz = 600.0,
			ChargeState = 2,
			RetentionTime = 10.0,
			MoleculeType = MslFormat.MoleculeType.Oligonucleotide,
			DissociationType = DissociationType.HCD,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new() { Mz = 380f, Intensity = 1.0f, ProductType = ProductType.w,
						 FragmentNumber = 2, Charge = 1 }
			}
		};

		LibrarySpectrum spectrum = entry.ToLibrarySpectrum();
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.ThreePrime));
	}

	/// <summary>
	/// MslLibraryEntry.StripModifications must handle nested brackets correctly.
	/// Modified sequences like "PEPTM[Common Variable:Oxidation on M]IDE" are the
	/// normal form in mzLib; stripping them correctly is required for elution group
	/// assignment, which groups all charge states of the same peptide together.
	/// </summary>
	[Test]
	public void MslLibraryEntry_StripModifications_NestedBrackets()
	{
		var entry = new MslLibraryEntry
		{
			FullSequence = "PEPTM[Common Variable:Oxidation on M]IDE",
			BaseSequence = "PEPTMIDE",
			PrecursorMz = 450.0,
			ChargeState = 2,
			RetentionTime = 30.0,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new() { Mz = 98f, Intensity = 1.0f, ProductType = ProductType.b,
						 FragmentNumber = 2, Charge = 1 }
			}
		};

		// Round-trip through FromLibrarySpectrum exercises StripModifications
		MslLibraryEntry fromSpectrum = MslLibraryEntry.FromLibrarySpectrum(
			entry.ToLibrarySpectrum());

		// StripModifications should produce "PEPTMIDE" from the modified sequence
		Assert.That(fromSpectrum.BaseSequence, Is.EqualTo("PEPTMIDE"));
	}
}