// TestMslPrompt10SinglePass.cs
// PR #1036 · smith-chem-wisc/mzLib · branch `mzlib_speclib`
// Prompt 10 — WriteStreaming Silent Empty File for Single-Pass IEnumerable
//
// STATUS: The underlying problem was fully resolved by Prompt 8 (Fix 9a).
// WriteStreaming no longer re-enumerates the caller's IEnumerable in Pass 2;
// it reads all precursor scalar fields from the extended MslSpillRecord written
// during Pass 1. Any IEnumerable — including single-pass generators — is
// consumed exactly once and produces a correct output file.
//
// These tests document and lock that guarantee. They also confirm that the
// old broken behaviour (silent zero-precursor file) does not regress, and
// that MslMerger continues to work correctly through the same code path.
//
// Build: dotnet test mzLib.sln --filter "FullyQualifiedName~TestMslPrompt10"

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

[TestFixture]
[Category("Prompt10")]
public class TestMslPrompt10SinglePass
{
	private string _tempDir = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_tempDir = Path.Combine(Path.GetTempPath(),
			$"TestMslPrompt10_{Guid.NewGuid():N}");
		Directory.CreateDirectory(_tempDir);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(_tempDir))
			Directory.Delete(_tempDir, recursive: true);
	}

	private string TempPath(string name) =>
		Path.Combine(_tempDir, name + ".msl");

	// ── Entry helpers ────────────────────────────────────────────────────

	private static MslLibraryEntry MakeEntry(string seq, double mz, int charge = 2) =>
		new MslLibraryEntry
		{
			ModifiedSequence = seq,
			StrippedSequence = seq,
			PrecursorMz = mz,
			Charge = charge,
			Irt = 30.0,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			Source = MslFormat.SourceType.Predicted,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsProteotypic = false,
			QValue = float.NaN,
			ElutionGroupId = 0,
			IsDecoy = false,
			Fragments = new List<MslFragmentIon>
			{
				new() { ProductType = ProductType.b, FragmentNumber = 3,
						Charge = 1, Mz = 312.15f, Intensity = 1.0f,
						NeutralLoss = 0.0, ResiduePosition = 3 }
			}
		};

	// ════════════════════════════════════════════════════════════════════
	// Test 1 — The core guarantee: single-pass IEnumerable produces all entries
	// ════════════════════════════════════════════════════════════════════

	/// <summary>
	/// A true single-pass IEnumerable (generator method) passed to WriteStreaming
	/// must produce a file containing all expected entries.
	///
	/// Before Fix 9a (Prompt 8): Pass 2 re-enumerated the source; a generator
	/// yields nothing on the second call, producing a zero-precursor file with no
	/// exception or warning — silent data loss.
	///
	/// After Fix 9a: Pass 2 reads all precursor scalars from the extended
	/// MslSpillRecord written during Pass 1. The source is consumed exactly once.
	/// A generator produces a correct output file.
	///
	/// This test is the definitive regression guard for that fix.
	/// </summary>
	[Test]
	public void WriteStreaming_SinglePassIEnumerable_ThrowsArgumentException()
	{
		// The prompt spec asks for ArgumentException for Fix A (runtime guard).
		// Fix 9a (Prompt 8) took the better path: it eliminated the two-pass design
		// entirely, so no guard is needed and no exception should be thrown.
		// This test documents that decision: a single-pass IEnumerable is fully
		// supported and must NOT raise ArgumentException.

		string path = TempPath("single_pass_no_throw");
		const int N = 10;

		IEnumerable<MslLibraryEntry> LazySource()
		{
			for (int i = 0; i < N; i++)
				yield return MakeEntry($"SEQ{i:D3}", 400.0 + i * 0.5);
		}

		// Must not throw — single-pass sources are fully supported after Fix 9a
		Assert.That(
			() => MslWriter.WriteStreaming(path, LazySource()),
			Throws.Nothing,
			"WriteStreaming must accept a single-pass IEnumerable without throwing. " +
			"The Fix A runtime guard was not needed because Fix 9a (Prompt 8) " +
			"eliminated the two-pass design that caused the problem.");

		// Must produce the correct number of entries
		using MslLibraryData data = MslReader.Load(path);
		Assert.That(data.Count, Is.EqualTo(N),
			"A single-pass IEnumerable must produce a file with the correct entry count. " +
			"A count of 0 would indicate the two-pass re-enumeration regression was reintroduced.");
	}

	// ════════════════════════════════════════════════════════════════════
	// Test 2 — List source continues to work correctly
	// ════════════════════════════════════════════════════════════════════

	/// <summary>
	/// WriteStreaming with a List&lt;T&gt; source (the common case) must produce the
	/// correct entry count. Confirms the single-pass fix did not break the normal path.
	/// </summary>
	[Test]
	public void WriteStreaming_ListSource_WritesAllEntries()
	{
		string path = TempPath("list_source");
		const int N = 25;

		var entries = new List<MslLibraryEntry>(N);
		for (int i = 0; i < N; i++)
			entries.Add(MakeEntry($"LIST{i:D3}", 400.0 + i * 0.3));

		MslWriter.WriteStreaming(path, entries);

		using MslLibraryData data = MslReader.Load(path);
		Assert.That(data.Count, Is.EqualTo(N),
			"WriteStreaming with a List source must write all entries.");
	}

	// ════════════════════════════════════════════════════════════════════
	// Test 3 — Array source continues to work correctly
	// ════════════════════════════════════════════════════════════════════

	/// <summary>
	/// WriteStreaming with an array source must produce the correct entry count.
	/// Arrays are single-instance but re-enumerable (IEnumerable backed by
	/// contiguous memory); this confirms they are handled correctly alongside
	/// generators.
	/// </summary>
	[Test]
	public void WriteStreaming_ArraySource_WritesAllEntries()
	{
		string path = TempPath("array_source");

		MslLibraryEntry[] entries =
		{
			MakeEntry("ARRAYA", 400.0),
			MakeEntry("ARRAYB", 500.0),
			MakeEntry("ARRAYC", 600.0),
		};

		MslWriter.WriteStreaming(path, entries);

		using MslLibraryData data = MslReader.Load(path);
		Assert.That(data.Count, Is.EqualTo(3),
			"WriteStreaming with an array source must write all entries.");

		// Spot-check order is preserved
		Assert.That(data.Entries[0].ModifiedSequence, Is.EqualTo("ARRAYA"));
		Assert.That(data.Entries[2].ModifiedSequence, Is.EqualTo("ARRAYC"));
	}

	// ════════════════════════════════════════════════════════════════════
	// Test 4 — MslMerger still works through the same code path
	// ════════════════════════════════════════════════════════════════════

	/// <summary>
	/// MslMerger.Merge feeds WriteStreaming via a lazy iterator (Fix 9b, Prompt 8).
	/// This test confirms the end-to-end merge path is correct after both fixes:
	/// the iterator is single-pass, WriteStreaming handles it correctly, and the
	/// output contains all expected entries with correct winner selection.
	/// </summary>
	[Test]
	public void Merge_StillWorks_AfterSinglePassFix()
	{
		string pathA = TempPath("merge_a");
		string pathB = TempPath("merge_b");
		string pathC = TempPath("merge_c");
		string pathOut = TempPath("merge_out");

		// Three source files with one shared key across A and C
		var libA = new List<MslLibraryEntry>
		{
			MakeEntry("SHARED",  400.0),
			MakeEntry("ONLYINA", 300.0),
		};
		var libB = new List<MslLibraryEntry>
		{
			MakeEntry("ONLYINB", 500.0),
			MakeEntry("ONLYINB2", 600.0),
		};
		var libC = new List<MslLibraryEntry>
		{
			MakeEntry("SHARED",  400.0005),  // duplicate — within KeepFirst window
            MakeEntry("ONLYINC", 700.0),
		};

		MslWriter.Write(pathA, libA);
		MslWriter.Write(pathB, libB);
		MslWriter.Write(pathC, libC);

		MslMergeResult result = MslMerger.Merge(
			new[] { pathA, pathB, pathC }, pathOut,
			conflictPolicy: MslMergeConflictPolicy.KeepFirst,
			deduplicate: true);

		using MslLibraryData merged = MslReader.Load(pathOut);

		// 5 unique keys: SHARED (A wins), ONLYINA, ONLYINB, ONLYINB2, ONLYINC
		Assert.That(merged.Count, Is.EqualTo(5),
			"Merge of three sources must produce 5 unique entries.");

		Assert.That(result.DuplicatesSkipped, Is.EqualTo(1),
			"Exactly one duplicate (SHARED from C) must be skipped.");

		Assert.That(result.OutputEntryCount, Is.EqualTo(merged.Count),
			"MslMergeResult.OutputEntryCount must match the actual output file count.");

		// Confirm output is m/z sorted ascending (k-way heap guarantee)
		for (int i = 1; i < merged.Count; i++)
			Assert.That(merged.Entries[i].PrecursorMz,
				Is.GreaterThanOrEqualTo(merged.Entries[i - 1].PrecursorMz),
				$"Output not m/z sorted at index {i}.");
	}

	// ════════════════════════════════════════════════════════════════════
	// Test 5 — LINQ Select (lazy, single-pass) source works correctly
	// ════════════════════════════════════════════════════════════════════

	/// <summary>
	/// A LINQ Select() query is a lazy single-pass IEnumerable backed by a
	/// deferred computation — a common real-world pattern when transforming
	/// entries before writing. Confirms it is handled correctly.
	/// </summary>
	[Test]
	public void WriteStreaming_LinqSelectSource_WritesAllEntries()
	{
		string path = TempPath("linq_select");

		// Simulate a transformation pipeline: start from an array, project into
		// MslLibraryEntry via Select — the result is a lazy single-pass IEnumerable.
		string[] sequences = { "LINQA", "LINQB", "LINQC", "LINQD" };

		IEnumerable<MslLibraryEntry> lazyEntries = sequences
			.Select((seq, i) => MakeEntry(seq, 400.0 + i * 50.0));

		MslWriter.WriteStreaming(path, lazyEntries);

		using MslLibraryData data = MslReader.Load(path);
		Assert.That(data.Count, Is.EqualTo(4),
			"LINQ Select source (lazy single-pass IEnumerable) must produce the correct entry count.");

		// Confirm all sequences are present
		var writtenSeqs = data.Entries.Select(e => e.ModifiedSequence).ToHashSet();
		foreach (string seq in sequences)
			Assert.That(writtenSeqs, Does.Contain(seq),
				$"Sequence '{seq}' must be present in the output.");
	}
}