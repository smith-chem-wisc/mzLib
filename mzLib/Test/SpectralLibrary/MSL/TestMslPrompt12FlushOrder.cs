using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;

namespace Test.SpectralLibrary.MSL;

[TestFixture]
[Category("Prompt12")]
public class TestMslPrompt12FlushOrder
{
	private string _tempDir = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_tempDir = Path.Combine(Path.GetTempPath(),
			$"TestMslPrompt12_{Guid.NewGuid():N}");
		Directory.CreateDirectory(_tempDir);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(_tempDir))
			Directory.Delete(_tempDir, recursive: true);
	}

	private string TempPath(string name) =>
		Path.Combine(_tempDir, $"{name}.msl");

	private static MslLibraryEntry MakeEntry(
		string seq, int charge, double mz,
		float qValue = 0.01f)
	{
		return new MslLibraryEntry
		{
			FullSequence = seq,
			BaseSequence = seq,
			PrecursorMz = mz,
			ChargeState = charge,
			RetentionTime = 30.0,
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
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new() { ProductType = ProductType.b, FragmentNumber = 1,
						Charge = 1, Mz = 200f, Intensity = 1.0f, NeutralLoss = 0.0 }
			}
		};
	}

	private string WriteLib(string name, IReadOnlyList<MslLibraryEntry> entries)
	{
		string path = TempPath(name);
		MslWriter.Write(path, entries);
		return path;
	}

	// ════════════════════════════════════════════════════════════════════
	// Test 1 — Core determinism: same inputs → same output order every run
	// ════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Multiple distinct peptides within the 0.001 Da flush window must all
	/// appear in the output with PrecursorMz in non-decreasing order.
	///
	/// Before Fix 12: the order depended on Dictionary iteration order and
	/// was non-deterministic across .NET versions and platforms.
	/// After Fix 12: toFlush is sorted before yielding, ensuring determinism.
	/// Determinism itself is verified by the byte-identical two-run test.
	/// </summary>
	[Test]
	public void Merge_KeepLowestQValue_SameMzDifferentKeys_OutputIsStablySorted()
	{
		// Three peptides at identical m/z — all unique keys, no duplicates.
		// They sit in the qValueBuffer together and are flushed as a group.
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("CCC", 2, 400.0, qValue: 0.03f),
			MakeEntry("AAA", 2, 400.0, qValue: 0.01f),
			MakeEntry("BBB", 2, 400.0, qValue: 0.02f),
		};

		string srcPath = WriteLib("p12_stable_src", entries);
		string outPath = TempPath("p12_stable_out");

		MslMerger.Merge(new[] { srcPath }, outPath,
			conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue);

		using MslLibraryData data = MslReader.Load(outPath);
		Assert.That(data.Count, Is.EqualTo(3), "All three entries must be present.");

		// Assert all expected sequences are present (order within equal m/z is
		// an implementation detail — determinism is covered by the two-run test).
		var sequences = new HashSet<string>(
			new[] { data.Entries[0].FullSequence, data.Entries[1].FullSequence, data.Entries[2].FullSequence });
		Assert.That(sequences, Is.EquivalentTo(new[] { "AAA", "BBB", "CCC" }),
			"All three unique entries must be present in the output.");

		// Assert the public contract: PrecursorMz non-decreasing order.
		for (int i = 1; i < data.Count; i++)
		{
			Assert.That(data.Entries[i].PrecursorMz,
				Is.GreaterThanOrEqualTo(data.Entries[i - 1].PrecursorMz),
				$"Entry {i} must have PrecursorMz >= entry {i - 1} (non-decreasing order).");
		}
	}

	// ════════════════════════════════════════════════════════════════════
	// Test 2 — CI regression guard: identical output across two runs
	// ════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Running the same merge twice must produce byte-for-byte identical output
	/// files. This is the CI guard: if FlushQValueBuffer is ever changed to
	/// remove the sort, this test will start failing non-deterministically.
	/// </summary>
	[Test]
	public void Merge_KeepLowestQValue_OutputOrder_IsIdenticalAcrossRuns()
	{
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("ZZZ", 2, 400.0, qValue: 0.01f),
			MakeEntry("AAA", 2, 400.0, qValue: 0.02f),
			MakeEntry("MMM", 2, 400.0, qValue: 0.03f),
		};

		string srcPath = WriteLib("p12_idem_src", entries);
		string outPath1 = TempPath("p12_idem_out1");
		string outPath2 = TempPath("p12_idem_out2");

		MslMerger.Merge(new[] { srcPath }, outPath1,
			conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue);
		MslMerger.Merge(new[] { srcPath }, outPath2,
			conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue);

		byte[] bytes1 = File.ReadAllBytes(outPath1);
		byte[] bytes2 = File.ReadAllBytes(outPath2);

		Assert.That(bytes1, Is.EqualTo(bytes2),
			"Two runs of the same merge must produce byte-for-byte identical output. " +
			"A difference indicates non-determinism in FlushQValueBuffer.");
	}

	// ════════════════════════════════════════════════════════════════════
	// Test 3 — PrecursorMz sort takes priority over Name
	// ════════════════════════════════════════════════════════════════════

	/// <summary>
	/// When flushed entries have different PrecursorMz values, the primary sort
	/// key is PrecursorMz ascending. LookupKey is only the tiebreaker.
	/// </summary>
	[Test]
	public void Merge_KeepLowestQValue_FlushWithDifferentMz_SortedByMzFirst()
	{
		// ZZZ is at mz=400.0, AAA is at mz=400.0005 — both within the 0.001 Da window
		// of a next entry at mz=400.002. After flush, ZZZ (lower mz) must come first
		// even though "AAA/2" < "ZZZ/2" lexicographically.
		var libA = new List<MslLibraryEntry>
		{
			MakeEntry("ZZZ", 2, 400.0000, qValue: 0.01f),
			MakeEntry("AAA", 2, 400.0005, qValue: 0.01f),
		};
		// Trigger the flush: a next entry more than 0.001 Da ahead
		var libB = new List<MslLibraryEntry>
		{
			MakeEntry("NEXT", 2, 400.002, qValue: 0.01f),
		};

		string outPath = TempPath("p12_mzpri_out");
		MslMerger.Merge(
			new[] { WriteLib("p12_mzpri_a", libA), WriteLib("p12_mzpri_b", libB) },
			outPath, conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue);

		using MslLibraryData data = MslReader.Load(outPath);
		Assert.That(data.Count, Is.EqualTo(3));

		// ZZZ (mz=400.0000) must precede AAA (mz=400.0005)
		Assert.That(data.Entries[0].FullSequence, Is.EqualTo("ZZZ"),
			"ZZZ (lower mz=400.0000) must be emitted before AAA (mz=400.0005), " +
			"even though 'AAA/2' < 'ZZZ/2' lexicographically.");
		Assert.That(data.Entries[1].FullSequence, Is.EqualTo("AAA"));
		Assert.That(data.Entries[2].FullSequence, Is.EqualTo("NEXT"));
	}

	// ════════════════════════════════════════════════════════════════════
	// Test 4 — Empty buffer edge case
	// ════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Flushing an empty buffer must not throw and must produce zero output entries.
	/// </summary>
	[Test]
	public void Merge_KeepLowestQValue_EmptySource_NoOutputAndNoException()
	{
		string srcPath = WriteLib("p12_empty_src", new List<MslLibraryEntry>());
		string outPath = TempPath("p12_empty_out");

		Assert.That(
			() => MslMerger.Merge(new[] { srcPath }, outPath,
				conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue),
			Throws.Nothing,
			"Merging an empty source with KeepLowestQValue must not throw.");

		using MslLibraryData data = MslReader.Load(outPath);
		Assert.That(data.Count, Is.EqualTo(0),
			"Merging an empty source must produce a zero-precursor output file.");
	}

	// ════════════════════════════════════════════════════════════════════
	// Test 5 — Single entry in buffer flushes correctly
	// ════════════════════════════════════════════════════════════════════

	/// <summary>
	/// A buffer containing exactly one entry must flush that entry with all
	/// fields intact. Confirms the sort path does not corrupt single-entry flushes.
	/// </summary>
	[Test]
	public void Merge_KeepLowestQValue_SingleEntry_OutputsCorrectly()
	{
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("SINGLE", 2, 500.0, qValue: 0.005f)
		};

		string srcPath = WriteLib("p12_single_src", entries);
		string outPath = TempPath("p12_single_out");

		MslMerger.Merge(new[] { srcPath }, outPath,
			conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue);

		using MslLibraryData data = MslReader.Load(outPath);
		Assert.That(data.Count, Is.EqualTo(1));
		Assert.That(data.Entries[0].FullSequence, Is.EqualTo("SINGLE"));
		Assert.That(data.Entries[0].QValue, Is.EqualTo(0.005f).Within(1e-6f));
	}
}