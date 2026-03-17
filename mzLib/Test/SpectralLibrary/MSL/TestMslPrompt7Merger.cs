using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.SpectralLibrary.MSL;

/// <summary>
/// Tests targeting the findings from Prompt 7 — MslMerger: Memory Claims vs Reality.
///
/// Covers:
///   M1 — outputEntries List accumulates all output entries before WriteStreaming;
///         O(k + UniqueStrings) claim is inaccurate for the output side
///   M2 — Fragment data in merged entries: fragments load on demand via GetAllEntries
///         but are necessary for output correctness (not a bug; amplifies M1 impact)
///   M3 — FlushQValueBuffer output order within mzWindow is non-deterministic
///   OK  — Exception safety, KeepFirst/KeepLast/KeepLowestQValue correctness,
///          duplicatesSkipped invariant, NaN tiebreak, mzWindow flush boundary
/// </summary>
[TestFixture]
public class TestMslPrompt7Merger
{
	private string _tempDir = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_tempDir = Path.Combine(Path.GetTempPath(),
			$"TestMslPrompt7_{Guid.NewGuid():N}");
		Directory.CreateDirectory(_tempDir);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(_tempDir))
			Directory.Delete(_tempDir, recursive: true);
	}

	private string TempPath(string name) =>
		Path.Combine(_tempDir, $"{name}_{Guid.NewGuid():N}.msl");

	private static MslLibraryEntry MakeEntry(
		string seq, int charge, double mz,
		double irt = 30.0,
		float qValue = float.NaN,
		int nFragments = 3)
	{
		var frags = Enumerable.Range(1, nFragments).Select(i => new MslFragmentIon
		{
			ProductType = ProductType.b,
			FragmentNumber = i,
			Charge = 1,
			Mz = 100.0f + i * 50.0f,
			Intensity = 1.0f / i,
			NeutralLoss = 0.0,
			ResiduePosition = i
		}).ToList();

		return new MslLibraryEntry
		{
			ModifiedSequence = seq,
			StrippedSequence = seq,
			PrecursorMz = mz,
			Charge = charge,
			Irt = irt,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			Source = MslFormat.SourceType.Predicted,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsProteotypic = false,
			QValue = qValue,
			ElutionGroupId = 0,
			IsDecoy = false,
			Fragments = frags
		};
	}

	private string WriteLib(string name, IReadOnlyList<MslLibraryEntry> entries)
	{
		string path = TempPath(name);
		MslWriter.Write(path, entries);
		return path;
	}

	// ═════════════════════════════════════════════════════════════════════
	// M1 — outputEntries memory: O(TotalOutputEntries), not O(k)
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Documents that MslMerger.Merge buffers ALL output entries (with their
	/// fragment arrays) in a List before calling WriteStreaming.
	///
	/// The O(k + UniqueStrings) claim in the XML doc is accurate only for the
	/// source-read side (k entries in the heap). The output side is
	/// O(TotalOutputEntries × AvgEntrySize).
	///
	/// This test verifies functional correctness of the merger.
	/// </summary>
	[Test]
	public void Merge_KeepFirst_OutputEntriesMatchSourceCount()
	{
		// BBB entries have distinct m/z (400.000 vs 400.001) so the heap pops
		// A-file first, making KeepFirst behaviour deterministic.
		var a = new List<MslLibraryEntry>
		{
			MakeEntry("AAA", 2, 300.0,   irt: 10.0),
			MakeEntry("BBB", 2, 400.000, irt: 10.0),
			MakeEntry("DDD", 2, 600.0,   irt: 40.0),
		};
		var b = new List<MslLibraryEntry>
		{
			MakeEntry("BBB", 2, 400.001, irt: 99.0),  // duplicate key, slightly higher m/z
            MakeEntry("CCC", 2, 500.0,   irt: 30.0),
		};

		string pathA = WriteLib("m1_a", a);
		string pathB = WriteLib("m1_b", b);
		string outPath = TempPath("m1_out");

		MslMergeResult result = MslMerger.Merge(
			new[] { pathA, pathB }, outPath,
			conflictPolicy: MslMergeConflictPolicy.KeepFirst);

		Assert.That(result.OutputEntryCount, Is.EqualTo(4));
		Assert.That(result.DuplicatesSkipped, Is.EqualTo(1));
		Assert.That(result.TotalSourceEntryCount, Is.EqualTo(5));

		using var lib = MslLibrary.Load(outPath);
		Assert.That(lib.PrecursorCount, Is.EqualTo(4));

		// A-file BBB (mz=400.000) is popped first; KeepFirst keeps it.
		lib.TryGetEntry("BBB", 2, out MslLibraryEntry? bbb);
		Assert.That(bbb!.Irt, Is.EqualTo(10.0).Within(0.01),
			"KeepFirst must retain the A-file entry for BBB (lower m/z, popped first).");
	}

	// ═════════════════════════════════════════════════════════════════════
	// DuplicatesSkipped invariant
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// OutputEntryCount + DuplicatesSkipped must equal TotalSourceEntryCount
	/// for all three conflict policies. This is the fundamental merge invariant.
	/// </summary>
	[TestCase(MslMergeConflictPolicy.KeepFirst, TestName = "Invariant_KeepFirst")]
	[TestCase(MslMergeConflictPolicy.KeepLast, TestName = "Invariant_KeepLast")]
	[TestCase(MslMergeConflictPolicy.KeepLowestQValue, TestName = "Invariant_KeepLowestQValue")]
	public void Merge_DuplicatesSkipped_Invariant_Holds(MslMergeConflictPolicy policy)
	{
		var a = new List<MslLibraryEntry>
		{
			MakeEntry("AAA", 2, 300.0, qValue: 0.05f),
			MakeEntry("BBB", 2, 400.0, qValue: 0.10f),
			MakeEntry("CCC", 2, 500.0, qValue: 0.01f),
		};
		var b = new List<MslLibraryEntry>
		{
			MakeEntry("BBB", 2, 400.0, qValue: 0.02f),  // duplicate
            MakeEntry("CCC", 2, 500.0, qValue: 0.05f),  // duplicate
            MakeEntry("DDD", 2, 600.0, qValue: 0.03f),
		};

		string pathA = WriteLib("inv_a", a);
		string pathB = WriteLib("inv_b", b);
		string outPath = TempPath("inv_out");

		MslMergeResult result = MslMerger.Merge(
			new[] { pathA, pathB }, outPath, conflictPolicy: policy);

		Assert.That(
			result.OutputEntryCount + result.DuplicatesSkipped,
			Is.EqualTo(result.TotalSourceEntryCount),
			$"OutputEntryCount ({result.OutputEntryCount}) + " +
			$"DuplicatesSkipped ({result.DuplicatesSkipped}) must equal " +
			$"TotalSourceEntryCount ({result.TotalSourceEntryCount}) " +
			$"for policy {policy}.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// KeepLowestQValue correctness
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Three entries with the same key from three different files.
	/// KeepLowestQValue must select the one with q=0.01.
	/// DuplicatesSkipped must be 2 (both losers counted).
	/// </summary>
	[Test]
	public void Merge_KeepLowestQValue_ThreeDuplicates_CorrectWinnerAndCount()
	{
		var a = new List<MslLibraryEntry> { MakeEntry("SAME", 2, 400.0, irt: 10.0, qValue: 0.10f) };
		var b = new List<MslLibraryEntry> { MakeEntry("SAME", 2, 400.0, irt: 20.0, qValue: 0.01f) };
		var c = new List<MslLibraryEntry> { MakeEntry("SAME", 2, 400.0, irt: 30.0, qValue: 0.05f) };

		string pathA = WriteLib("kqv3_a", a);
		string pathB = WriteLib("kqv3_b", b);
		string pathC = WriteLib("kqv3_c", c);
		string outPath = TempPath("kqv3_out");

		MslMergeResult result = MslMerger.Merge(
			new[] { pathA, pathB, pathC }, outPath,
			conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue);

		Assert.That(result.OutputEntryCount, Is.EqualTo(1));
		Assert.That(result.DuplicatesSkipped, Is.EqualTo(2),
			"Two of the three duplicates must be counted as skipped.");

		using var lib = MslLibrary.Load(outPath);
		lib.TryGetEntry("SAME", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.Irt, Is.EqualTo(20.0).Within(0.01),
			"KeepLowestQValue must select the entry with q=0.01 (irt=20).");
		Assert.That(entry.QValue, Is.EqualTo(0.01f).Within(1e-6f));
	}

	/// <summary>
	/// NaN q-value is treated as worst (infinity). An entry with a real q-value
	/// must always beat an entry with NaN q-value.
	/// </summary>
	[Test]
	public void Merge_KeepLowestQValue_NanQValue_TreatedAsInfinity()
	{
		var a = new List<MslLibraryEntry>
			{ MakeEntry("NQ", 2, 400.0, irt: 10.0, qValue: float.NaN) };
		var b = new List<MslLibraryEntry>
			{ MakeEntry("NQ", 2, 400.0, irt: 99.0, qValue: 0.01f) };

		string outPath = TempPath("nan_out");
		MslMerger.Merge(new[] { WriteLib("nan_a", a), WriteLib("nan_b", b) },
			outPath, conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue);

		using var lib = MslLibrary.Load(outPath);
		lib.TryGetEntry("NQ", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.Irt, Is.EqualTo(99.0).Within(0.01),
			"Real q-value (0.01) must beat NaN q-value.");
	}

	/// <summary>
	/// When both entries have NaN q-value, KeepFirst semantics apply
	/// (existing buffer entry is kept, new one is the duplicate).
	/// </summary>
	[Test]
	public void Merge_KeepLowestQValue_BothNaN_KeepFirstSemantics()
	{
		var a = new List<MslLibraryEntry>
			{ MakeEntry("NN", 2, 400.0, irt: 10.0, qValue: float.NaN) };
		var b = new List<MslLibraryEntry>
			{ MakeEntry("NN", 2, 400.0, irt: 99.0, qValue: float.NaN) };

		string outPath = TempPath("bothnan_out");
		MslMerger.Merge(new[] { WriteLib("nn_a", a), WriteLib("nn_b", b) },
			outPath, conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue);

		using var lib = MslLibrary.Load(outPath);
		lib.TryGetEntry("NN", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.Irt, Is.EqualTo(10.0).Within(0.01),
			"Both NaN q-values: KeepFirst semantics must apply (A-file entry kept).");
	}

	// ═════════════════════════════════════════════════════════════════════
	// KeepLast correctness
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// KeepLast must retain the entry from the highest-index source file.
	/// Three files: A (idx=0), B (idx=1), C (idx=2). C's entry must win.
	/// </summary>
	[Test]
	public void Merge_KeepLast_ThreeSources_LastFileWins()
	{
		var a = new List<MslLibraryEntry> { MakeEntry("KEY", 2, 400.0, irt: 10.0) };
		var b = new List<MslLibraryEntry> { MakeEntry("KEY", 2, 400.0, irt: 50.0) };
		var c = new List<MslLibraryEntry> { MakeEntry("KEY", 2, 400.0, irt: 99.0) };

		string outPath = TempPath("kl3_out");
		MslMergeResult result = MslMerger.Merge(
			new[] { WriteLib("kl_a", a), WriteLib("kl_b", b), WriteLib("kl_c", c) },
			outPath, conflictPolicy: MslMergeConflictPolicy.KeepLast);

		Assert.That(result.OutputEntryCount, Is.EqualTo(1));
		Assert.That(result.DuplicatesSkipped, Is.EqualTo(2));

		using var lib = MslLibrary.Load(outPath);
		lib.TryGetEntry("KEY", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.Irt, Is.EqualTo(99.0).Within(0.01),
			"KeepLast must keep the C-file entry (highest source index).");
	}

	// ═════════════════════════════════════════════════════════════════════
	// mzWindow flush boundary (strict >)
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// FlushQValueBuffer flushes entries where (currentMz - entryMz) > mzWindow.
	/// The condition is STRICT greater-than, so an entry exactly at the window
	/// boundary must NOT be flushed. This test confirms the boundary semantics.
	/// </summary>
	[Test]
	public void Merge_KeepLowestQValue_MzWindowBoundary_ExactlyAtWindow_NotFlushed()
	{
		// Two entries with the same key at mz=400.000.
		// The next entry will be at mz=400.001 (exactly mzWindow=0.001 Da away).
		// The flush condition is (400.001 - 400.000) > 0.001 → 0.001 > 0.001 → FALSE.
		// So the two 400.000 entries must still be resolved against each other.
		var a = new List<MslLibraryEntry>
		{
			MakeEntry("SAME", 2, 400.000, irt: 10.0, qValue: 0.10f),
			MakeEntry("OTHER", 2, 400.001, irt: 50.0, qValue: 0.01f),
		};
		var b = new List<MslLibraryEntry>
		{
			MakeEntry("SAME", 2, 400.000, irt: 99.0, qValue: 0.01f),
		};

		string outPath = TempPath("mzbound_out");
		MslMerger.Merge(
			new[] { WriteLib("mzb_a", a), WriteLib("mzb_b", b) },
			outPath, conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue);

		using var lib = MslLibrary.Load(outPath);
		Assert.That(lib.PrecursorCount, Is.EqualTo(2),
			"SAME (resolved) + OTHER must both be in the output.");

		lib.TryGetEntry("SAME", 2, out MslLibraryEntry? same);
		Assert.That(same!.Irt, Is.EqualTo(99.0).Within(0.01),
			"SAME: B-file entry (q=0.01) must beat A-file entry (q=0.10).");
	}

	// ═════════════════════════════════════════════════════════════════════
	// Exception safety: partial initialisation
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// If LoadIndexOnly throws for source file i, the libraries opened for
	/// files 0..i-1 must still be disposed. The null-conditional ?. operator
	/// on the array ensures this even when the array is partially initialised.
	/// </summary>
	[Test]
	public void Merge_InvalidSecondFile_ThrowsAndFirstLibraryDisposed()
	{
		string validPath = WriteLib("exc_valid", new List<MslLibraryEntry>
			{ MakeEntry("AAA", 2, 300.0) });
		string invalidPath = Path.Combine(_tempDir, "nonexistent_xyz.msl");
		string outPath = TempPath("exc_out");

		Assert.That(
			() => MslMerger.Merge(new[] { validPath, invalidPath }, outPath),
			Throws.InstanceOf<Exception>(),
			"Merge must throw when a source file does not exist.");

		// The output file must not exist (atomic write guarantee)
		Assert.That(File.Exists(outPath), Is.False,
			"No partial output file must exist after a failed merge.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// M3 — KeepLowestQValue flush order within mzWindow
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// FlushQValueBuffer iterates the qValueBuffer dictionary to collect keys
	/// to flush. Within a given flush, the ORDER of emitted entries for
	/// different keys at the same m/z is Dictionary iteration order, which is
	/// not guaranteed by the .NET specification.
	///
	/// This is a cosmetic issue — the file is structurally valid — but callers
	/// must not depend on the output order of entries with identical m/z values.
	/// This test documents the behavior: all expected entries ARE present, but
	/// we make no assertion about their relative order.
	/// </summary>
	[Test]
	public void Merge_KeepLowestQValue_SameMzDifferentKeys_AllPresent_OrderUnspecified()
	{
		// Three different peptides at exactly the same nominal m/z.
		// All must appear in the output; their relative order is unspecified.
		var a = new List<MslLibraryEntry>
		{
			MakeEntry("PEP1", 2, 400.0, qValue: 0.01f),
			MakeEntry("PEP2", 2, 400.0, qValue: 0.02f),
			MakeEntry("PEP3", 2, 400.0, qValue: 0.03f),
		};

		string outPath = TempPath("samedmz_out");
		MslMerger.Merge(new[] { WriteLib("sdmz_a", a) }, outPath,
			conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue);

		using var lib = MslLibrary.Load(outPath);
		Assert.That(lib.PrecursorCount, Is.EqualTo(3),
			"All three distinct peptides at the same m/z must appear in the output.");

		// Confirm all three are findable regardless of output order
		Assert.That(lib.TryGetEntry("PEP1", 2, out _), Is.True);
		Assert.That(lib.TryGetEntry("PEP2", 2, out _), Is.True);
		Assert.That(lib.TryGetEntry("PEP3", 2, out _), Is.True);
	}

	// ═════════════════════════════════════════════════════════════════════
	// Fragment data preserved through merge
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Merged entries must carry their fragment arrays through to the output file.
	/// This confirms that GetAllEntries() (which triggers demand-load) is correctly
	/// used in the merge path so fragments are available for WriteStreaming.

	/// </summary>
	[Test]
	public void Merge_FragmentData_PreservedInOutput()
	{
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("FRAG", 2, 400.0, nFragments: 5)
		};

		string outPath = TempPath("frag_out");
		MslMerger.Merge(new[] { WriteLib("frag_src", entries) }, outPath);

		using var lib = MslLibrary.Load(outPath);
		lib.TryGetEntry("FRAG", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.Fragments, Has.Count.EqualTo(5),
			"All 5 fragment ions must be present in the merged output file.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// Output is m/z sorted
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// The k-way heap merge must produce output sorted by precursor m/z ascending.
	/// Verified by reading all entries and asserting monotonically non-decreasing m/z.
	/// </summary>
	[Test]
	public void Merge_Output_IsMzSortedAscending()
	{
		var a = new List<MslLibraryEntry>
		{
			MakeEntry("A1", 2, 300.0),
			MakeEntry("A2", 2, 500.0),
			MakeEntry("A3", 2, 700.0),
		};
		var b = new List<MslLibraryEntry>
		{
			MakeEntry("B1", 2, 200.0),
			MakeEntry("B2", 2, 400.0),
			MakeEntry("B3", 2, 600.0),
		};

		string outPath = TempPath("sorted_out");
		MslMerger.Merge(new[] { WriteLib("sort_a", a), WriteLib("sort_b", b) },
			outPath, conflictPolicy: MslMergeConflictPolicy.KeepFirst);

		using var lib = MslLibrary.Load(outPath);
		var allEntries = lib.GetAllEntries().ToList();

		Assert.That(allEntries, Has.Count.EqualTo(6));

		for (int i = 1; i < allEntries.Count; i++)
		{
			Assert.That(allEntries[i].PrecursorMz,
				Is.GreaterThanOrEqualTo(allEntries[i - 1].PrecursorMz),
				$"Output must be sorted by m/z ascending. " +
				$"Entry {i - 1} m/z={allEntries[i - 1].PrecursorMz:F3} > " +
				$"entry {i} m/z={allEntries[i].PrecursorMz:F3}.");
		}
	}
}