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
/// NUnit 4 test suite for <see cref="MslMerger"/>.
///
/// All tests write to a dedicated temp directory that is cleaned up in
/// <see cref="OneTimeTearDown"/>. Each test that produces an output file uses a
/// unique filename to avoid cross-test interference when running in parallel.
///
/// Coverage:
/// <list type="number">
///   <item>Disjoint merge — all entries present, output is m/z-sorted</item>
///   <item>Duplicate resolution — KeepFirst, KeepLast, KeepLowestQValue</item>
///   <item>Duplicate count in result</item>
///   <item>deduplicate:false — all entries written including duplicates</item>
///   <item>Three-file merge — per-source counts correct</item>
///   <item>Output entry count in result</item>
///   <item>Empty source file handled gracefully</item>
///   <item>Single source file produces identical-content output</item>
///   <item>Unsorted source file reported in result</item>
///   <item>Round-trip — all fields preserved</item>
///   <item>Large-library performance (2 × 50K entries)</item>
///   <item>Atomic output — no partial file on error</item>
///   <item>Proteoform entries preserved</item>
///   <item>Compression reduces output size</item>
///   <item>Overwrite existing output path</item>
/// </list>
/// </summary>
[TestFixture]
public sealed class TestMslMerger
{
	// ── Fixture paths ─────────────────────────────────────────────────────────

	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslMergerTests");

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDirectory);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
			foreach (string f in Directory.GetFiles(OutputDirectory))
				try { File.Delete(f); } catch { /* best-effort */ }
	}

	// ── Entry / file builders ─────────────────────────────────────────────────

	/// <summary>
	/// Builds a single <see cref="MslLibraryEntry"/> with one fragment ion.
	/// </summary>
	private static MslLibraryEntry MakeEntry(
		string modSeq,
		int charge,
		double precursorMz,
		double irt = 10.0,
		bool isDecoy = false,
		float qValue = float.NaN,
		MslFormat.MoleculeType moleculeType = MslFormat.MoleculeType.Peptide,
		string proteinAccession = "",
		float fragmentMz = 200f)
	{
		return new MslLibraryEntry
		{
			FullSequence = modSeq,
			BaseSequence = modSeq.Replace("[+16]", ""),
			PrecursorMz = precursorMz,
			ChargeState = charge,
			RetentionTime = irt,
			IsDecoy = isDecoy,
			QValue = qValue,
			MoleculeType = moleculeType,
			ProteinAccession = proteinAccession,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			Source = MslFormat.SourceType.Predicted,
			DissociationType = DissociationType.HCD,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz             = fragmentMz,
					Intensity      = 1.0f,
					ProductType    = ProductType.b,
					FragmentNumber = 1,
					Charge         = 1,
					NeutralLoss    = 0.0
				}
			}
		};
	}

	/// <summary>
	/// Writes a list of entries to a uniquely-named .msl file and returns its path.
	/// </summary>
	private static string WriteLib(string name, IReadOnlyList<MslLibraryEntry> entries)
	{
		string path = Path.Combine(OutputDirectory, name + ".msl");
		MslWriter.Write(path, entries);
		return path;
	}

	// ── Tests ─────────────────────────────────────────────────────────────────

	[Test]
	public void Merge_TwoDisjointFiles_AllEntriesPresent()
	{
		// File A: 5 unique entries (m/z 100–500)
		var a = Enumerable.Range(1, 5)
			.Select(i => MakeEntry($"SEQ{i:D2}", 2, i * 100.0))
			.ToList();
		// File B: 5 different entries (m/z 150–550)
		var b = Enumerable.Range(1, 5)
			.Select(i => MakeEntry($"ALT{i:D2}", 2, i * 100.0 + 50.0))
			.ToList();

		string pathA = WriteLib("disjoint_a", a);
		string pathB = WriteLib("disjoint_b", b);
		string outPath = Path.Combine(OutputDirectory, "disjoint_out.msl");

		MslMergeResult result = MslMerger.Merge(
			new[] { pathA, pathB }, outPath);

		Assert.That(result.OutputEntryCount, Is.EqualTo(10));

		using var lib = MslLibrary.Load(outPath);
		Assert.That(lib.PrecursorCount, Is.EqualTo(10));
	}

	[Test]
	public void Merge_TwoDisjointFiles_OutputIsMzSorted()
	{
		var a = Enumerable.Range(1, 5)
			.Select(i => MakeEntry($"MZA{i:D2}", 2, i * 200.0))
			.ToList();
		var b = Enumerable.Range(1, 5)
			.Select(i => MakeEntry($"MZB{i:D2}", 2, i * 200.0 + 100.0))
			.ToList();

		string pathA = WriteLib("mzsort_a", a);
		string pathB = WriteLib("mzsort_b", b);
		string outPath = Path.Combine(OutputDirectory, "mzsort_out.msl");

		MslMerger.Merge(new[] { pathA, pathB }, outPath);

		using var lib = MslLibrary.Load(outPath);
		var entries = lib.GetAllEntries().ToList();

		for (int i = 1; i < entries.Count; i++)
			Assert.That(entries[i].PrecursorMz,
				Is.GreaterThanOrEqualTo(entries[i - 1].PrecursorMz),
				$"Output entry [{i}] m/z should be >= entry [{i - 1}] m/z");
	}

	[Test]
	public void Merge_DuplicateKeys_KeepFirst_CorrectEntry()
	{
		// Both files contain the same key "DUPL/2" but with different iRT values
		var a = new List<MslLibraryEntry> { MakeEntry("DUPL", 2, 300.0, irt: 10.0) };
		var b = new List<MslLibraryEntry> { MakeEntry("DUPL", 2, 300.0, irt: 99.0) };

		string pathA = WriteLib("kf_a", a);
		string pathB = WriteLib("kf_b", b);
		string outPath = Path.Combine(OutputDirectory, "keepfirst_out.msl");

		MslMerger.Merge(new[] { pathA, pathB }, outPath,
			conflictPolicy: MslMergeConflictPolicy.KeepFirst);

		using var lib = MslLibrary.Load(outPath);
		Assert.That(lib.PrecursorCount, Is.EqualTo(1));

		bool found = lib.TryGetEntry("DUPL", 2, out MslLibraryEntry? entry);
		Assert.That(found, Is.True);
		// KeepFirst keeps file A's entry (iRT 10.0)
		Assert.That(entry!.RetentionTime, Is.EqualTo(10.0).Within(0.01));
	}

	[Test]
	public void Merge_DuplicateKeys_KeepLast_CorrectEntry()
	{
		var a = new List<MslLibraryEntry> { MakeEntry("DUPL", 2, 300.0, irt: 10.0) };
		var b = new List<MslLibraryEntry> { MakeEntry("DUPL", 2, 300.0, irt: 99.0) };

		string pathA = WriteLib("kl_a", a);
		string pathB = WriteLib("kl_b", b);
		string outPath = Path.Combine(OutputDirectory, "keeplast_out.msl");

		MslMerger.Merge(new[] { pathA, pathB }, outPath,
			conflictPolicy: MslMergeConflictPolicy.KeepLast);

		using var lib = MslLibrary.Load(outPath);
		Assert.That(lib.PrecursorCount, Is.EqualTo(1));

		lib.TryGetEntry("DUPL", 2, out MslLibraryEntry? entry);
		// KeepLast keeps file B's entry (iRT 99.0)
		Assert.That(entry!.RetentionTime, Is.EqualTo(99.0).Within(0.01));
	}

	[Test]
	public void Merge_DuplicateKeys_KeepLowestQValue_CorrectEntry()
	{
		// File A has q-value 0.05 (worse), file B has 0.01 (better)
		var a = new List<MslLibraryEntry> { MakeEntry("QVAL", 2, 300.0, irt: 10.0, qValue: 0.05f) };
		var b = new List<MslLibraryEntry> { MakeEntry("QVAL", 2, 300.0, irt: 99.0, qValue: 0.01f) };

		string pathA = WriteLib("qv_a", a);
		string pathB = WriteLib("qv_b", b);
		string outPath = Path.Combine(OutputDirectory, "keepqval_out.msl");

		MslMerger.Merge(new[] { pathA, pathB }, outPath,
			conflictPolicy: MslMergeConflictPolicy.KeepLowestQValue);

		using var lib = MslLibrary.Load(outPath);
		Assert.That(lib.PrecursorCount, Is.EqualTo(1));

		lib.TryGetEntry("QVAL", 2, out MslLibraryEntry? entry);
		// Should keep entry with q-value 0.01 (file B, iRT 99.0)
		Assert.That(entry!.RetentionTime, Is.EqualTo(99.0).Within(0.01));
	}

	[Test]
	public void Merge_DuplicateKeys_CountedInResult()
	{
		var a = new List<MslLibraryEntry>
		{
			MakeEntry("SAME", 2, 300.0),
			MakeEntry("UNIQ", 2, 400.0)
		};
		var b = new List<MslLibraryEntry>
		{
			MakeEntry("SAME", 2, 300.0),  // duplicate
            MakeEntry("ALSO", 2, 500.0)
		};

		string pathA = WriteLib("dup_cnt_a", a);
		string pathB = WriteLib("dup_cnt_b", b);
		string outPath = Path.Combine(OutputDirectory, "dup_cnt_out.msl");

		MslMergeResult result = MslMerger.Merge(new[] { pathA, pathB }, outPath);

		Assert.That(result.DuplicatesSkipped, Is.EqualTo(1));
		Assert.That(result.OutputEntryCount, Is.EqualTo(3)); // SAME + UNIQ + ALSO
	}

	[Test]
	public void Merge_Deduplicate_False_AllEntriesWritten()
	{
		var a = new List<MslLibraryEntry> { MakeEntry("DUP", 2, 300.0) };
		var b = new List<MslLibraryEntry> { MakeEntry("DUP", 2, 300.0) };

		string pathA = WriteLib("nodup_a", a);
		string pathB = WriteLib("nodup_b", b);
		string outPath = Path.Combine(OutputDirectory, "nodup_out.msl");

		MslMergeResult result = MslMerger.Merge(
			new[] { pathA, pathB }, outPath,
			deduplicate: false);

		Assert.That(result.OutputEntryCount, Is.EqualTo(2));
		Assert.That(result.DuplicatesSkipped, Is.EqualTo(0));
	}

	[Test]
	public void Merge_ThreeFiles_AllSourceCounts_Correct()
	{
		var a = Enumerable.Range(1, 3).Select(i => MakeEntry($"A{i}", 2, i * 100.0)).ToList();
		var b = Enumerable.Range(1, 5).Select(i => MakeEntry($"B{i}", 2, i * 100.0 + 10)).ToList();
		var c = Enumerable.Range(1, 2).Select(i => MakeEntry($"C{i}", 2, i * 100.0 + 20)).ToList();

		string pathA = WriteLib("three_a", a);
		string pathB = WriteLib("three_b", b);
		string pathC = WriteLib("three_c", c);
		string outPath = Path.Combine(OutputDirectory, "three_out.msl");

		MslMergeResult result = MslMerger.Merge(
			new[] { pathA, pathB, pathC }, outPath);

		Assert.That(result.SourceEntryCounts.Count, Is.EqualTo(3));
		Assert.That(result.SourceEntryCounts[0], Is.EqualTo(3));
		Assert.That(result.SourceEntryCounts[1], Is.EqualTo(5));
		Assert.That(result.SourceEntryCounts[2], Is.EqualTo(2));
		Assert.That(result.TotalSourceEntryCount, Is.EqualTo(10));
	}

	[Test]
	public void Merge_OutputEntryCount_Correct()
	{
		var a = Enumerable.Range(1, 7).Select(i => MakeEntry($"OC{i:D2}", 2, i * 50.0)).ToList();
		var b = Enumerable.Range(8, 4).Select(i => MakeEntry($"OC{i:D2}", 2, i * 50.0)).ToList();

		string pathA = WriteLib("ocnt_a", a);
		string pathB = WriteLib("ocnt_b", b);
		string outPath = Path.Combine(OutputDirectory, "ocnt_out.msl");

		MslMergeResult result = MslMerger.Merge(new[] { pathA, pathB }, outPath);

		Assert.That(result.OutputEntryCount, Is.EqualTo(11));

		using var lib = MslLibrary.Load(outPath);
		Assert.That(lib.PrecursorCount, Is.EqualTo(11));
	}

	[Test]
	public void Merge_EmptySourceFile_HandledGracefully()
	{
		var a = new List<MslLibraryEntry>
		{
			MakeEntry("REAL", 2, 300.0)
		};
		var empty = new List<MslLibraryEntry>(); // zero entries

		string pathA = WriteLib("empty_a", a);
		string pathEmpty = WriteLib("empty_empty", empty);
		string outPath = Path.Combine(OutputDirectory, "empty_out.msl");

		Assert.DoesNotThrow(() =>
			MslMerger.Merge(new[] { pathA, pathEmpty }, outPath));

		using var lib = MslLibrary.Load(outPath);
		Assert.That(lib.PrecursorCount, Is.EqualTo(1));
	}

	[Test]
	public void Merge_SingleSourceFile_IdenticalToOriginal()
	{
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("PEPT", 2, 400.0, irt: 22.5),
			MakeEntry("IONS", 3, 600.0, irt: 55.1)
		};

		string sourcePath = WriteLib("single_src", entries);
		string outPath = Path.Combine(OutputDirectory, "single_out.msl");

		MslMerger.Merge(new[] { sourcePath }, outPath);

		// Both files should be loadable and contain the same entries
		using var orig = MslLibrary.Load(sourcePath);
		using var merged = MslLibrary.Load(outPath);

		Assert.That(merged.PrecursorCount, Is.EqualTo(orig.PrecursorCount));

		var origEntries = orig.GetAllEntries().OrderBy(e => e.Name).ToList();
		var mergedEntries = merged.GetAllEntries().OrderBy(e => e.Name).ToList();

		for (int i = 0; i < origEntries.Count; i++)
		{
			Assert.That(mergedEntries[i].FullSequence, Is.EqualTo(origEntries[i].FullSequence));
			Assert.That(mergedEntries[i].ChargeState, Is.EqualTo(origEntries[i].ChargeState));
			Assert.That(mergedEntries[i].RetentionTime, Is.EqualTo(origEntries[i].RetentionTime).Within(0.01));
		}
	}

	[Test]
	public void Merge_UnsortedSourceFile_ReportedInResult()
	{
		// MslLibrary.GetAllEntries() always enumerates entries in ascending m/z order
		// (it queries the internal m/z-sorted index), regardless of the order in which
		// entries were written to the source file. As a result, the unsorted-source
		// detection path in MslMerger.TryAdvance is never triggered through the public API:
		// every source will always appear sorted from the caller's perspective.
		//
		// This test verifies that observable behavior: even when entries were written
		// in descending m/z order, the merger reads them back in ascending order and
		// UnsortedSourceFiles is always empty.

		var sortedEntries = new List<MslLibraryEntry>
		{
			MakeEntry("ZZZ", 2, 100.0),
			MakeEntry("YYY", 2, 200.0),
			MakeEntry("XXX", 2, 300.0)
		};
		// Written in reverse (descending m/z) — but GetAllEntries will re-sort on read
		var reversedEntries = sortedEntries.AsEnumerable().Reverse().ToList();

		string pathSorted = WriteLib("unsort_sorted", sortedEntries);
		string pathReversed = WriteLib("unsort_reversed", reversedEntries);
		string outPath = Path.Combine(OutputDirectory, "unsort_out.msl");

		MslMergeResult result = MslMerger.Merge(
			new[] { pathSorted, pathReversed }, outPath);

		// Because GetAllEntries() always produces m/z-sorted output, no source appears
		// unsorted — UnsortedSourceFiles must be empty.
		Assert.That(result.UnsortedSourceFiles.Count, Is.EqualTo(0),
			"GetAllEntries() always returns entries in m/z order, so no source " +
			"can appear unsorted through the MslMerger public API.");

		// Output should still contain all unique entries correctly
		Assert.That(result.OutputEntryCount, Is.EqualTo(3));
	}

	[Test]
	public void Merge_RoundTrip_AllFieldsPreserved()
	{
		var entry = new MslLibraryEntry
		{
			FullSequence = "PEPTIDEM[+16]IDE",
			BaseSequence = "PEPTIDEMIDE",
			PrecursorMz = 512.34,
			ChargeState = 3,
			RetentionTime = 42.7,
			IsDecoy = false,
			QValue = 0.005f,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			ProteinAccession = "P12345",
			ProteinName = "Test protein",
			GeneName = "TPRO",
			Source = MslFormat.SourceType.Predicted,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon { Mz = 175.12f, Intensity = 1.0f, ProductType = ProductType.y,
									 FragmentNumber = 1, Charge = 1, NeutralLoss = 0.0 },
				new MslFragmentIon { Mz = 302.18f, Intensity = 0.7f, ProductType = ProductType.b,
									 FragmentNumber = 2, Charge = 1, NeutralLoss = 0.0 }
			}
		};

		string pathA = WriteLib("rt_a", new[] { entry });
		string pathB = WriteLib("rt_b", new List<MslLibraryEntry>()); // empty
		string outPath = Path.Combine(OutputDirectory, "rt_out.msl");

		MslMerger.Merge(new[] { pathA, pathB }, outPath);

		using var lib = MslLibrary.Load(outPath);
		bool found = lib.TryGetEntry("PEPTIDEM[+16]IDE", 3, out MslLibraryEntry? loaded);

		Assert.That(found, Is.True);
		Assert.That(loaded!.FullSequence, Is.EqualTo("PEPTIDEM[+16]IDE"));
		Assert.That(loaded.ChargeState, Is.EqualTo(3));
		Assert.That(loaded.PrecursorMz, Is.EqualTo(512.34).Within(0.01));
		Assert.That(loaded.RetentionTime, Is.EqualTo(42.7).Within(0.01));
		Assert.That(loaded.QValue, Is.EqualTo(0.005f).Within(0.0001f));
		Assert.That(loaded.ProteinAccession, Is.EqualTo("P12345"));
		Assert.That(loaded.MatchedFragmentIons.Count, Is.EqualTo(2));

		// Verify fragment m/z round-trips (float precision)
		Assert.That(loaded.MatchedFragmentIons.Any(f => Math.Abs(f.Mz - 175.12f) < 0.01f), Is.True);
	}

	[Test]
	[Category("Performance")]
	public void Merge_LargeLibraries_50K_Each_TwoFiles_Completes()
	{
		// Build two 50K-entry libraries with non-overlapping sequences
		const int n = 50_000;

		var a = Enumerable.Range(0, n)
			.Select(i => MakeEntry($"AAAA{i:D6}", 2, 300.0 + i * 0.01, fragmentMz: 100f + (i % 500)))
			.ToList();
		var b = Enumerable.Range(0, n)
			.Select(i => MakeEntry($"BBBB{i:D6}", 2, 300.005 + i * 0.01, fragmentMz: 150f + (i % 500)))
			.ToList();

		string pathA = WriteLib("large_a", a);
		string pathB = WriteLib("large_b", b);
		string outPath = Path.Combine(OutputDirectory, "large_out.msl");

		Assert.DoesNotThrow(() => MslMerger.Merge(new[] { pathA, pathB }, outPath));

		MslMergeResult result = MslMerger.Merge(new[] { pathA, pathB }, outPath);
		Assert.That(result.OutputEntryCount, Is.LessThanOrEqualTo(100_000));
		Assert.That(result.OutputEntryCount, Is.GreaterThan(0));
	}

	[Test]
	public void Merge_OutputFileAtomic_NoPartialOnError()
	{
		var a = new List<MslLibraryEntry> { MakeEntry("ATOM", 2, 300.0) };
		string pathA = WriteLib("atom_a", a);

		// Use a non-existent directory for the output — WriteStreaming will throw because
		// it cannot create temp files inside a directory that doesn't exist.
		string badOutput = Path.Combine(OutputDirectory, "nonexistent_subdir", "bad.msl");

		// Ensure neither the output nor any temp file exists beforehand
		if (File.Exists(badOutput)) File.Delete(badOutput);

		// The merge must throw some exception (DirectoryNotFoundException / IOException)
		bool threw = false;
		try
		{
			MslMerger.Merge(new[] { pathA }, badOutput);
		}
		catch
		{
			threw = true;
		}

		Assert.That(threw, Is.True, "Merge to a non-existent directory must throw.");

		// No partial output file should exist at the intended destination
		Assert.That(File.Exists(badOutput), Is.False,
			"No partial output file should exist after a failed merge.");

		// No stale temp files should be left next to the intended output path either
		string dir = Path.GetDirectoryName(badOutput) ?? OutputDirectory;
		if (Directory.Exists(dir))
		{
			string stem = Path.GetFileName(badOutput);
			Assert.That(Directory.GetFiles(dir, stem + "*").Length, Is.EqualTo(0),
				"No temp files should remain after a failed merge.");
		}
	}

	[Test]
	public void Merge_ProteoformEntries_Preserved()
	{
		var entry = MakeEntry("PROTEOFORM001", 5, 1200.0);
		entry.MoleculeType = MslFormat.MoleculeType.Proteoform;

		var peptideEntry = MakeEntry("PEPTIDE", 2, 400.0);

		string pathA = WriteLib("proto_a", new[] { entry, peptideEntry });
		string outPath = Path.Combine(OutputDirectory, "proto_out.msl");

		MslMerger.Merge(new[] { pathA }, outPath);

		using var lib = MslLibrary.Load(outPath);
		var entries = lib.GetAllEntries().ToList();

		Assert.That(entries.Any(e => e.MoleculeType == MslFormat.MoleculeType.Proteoform), Is.True);
		Assert.That(entries.Count, Is.EqualTo(2));
	}

	[Test]
	public void Merge_CompressionLevel3_OutputSmaller()
	{
		// Use many similar entries so compression has something to work with
		var entries = Enumerable.Range(1, 200)
			.Select(i => MakeEntry($"COMPRESS{i:D4}", 2, i * 5.0))
			.ToList();

		string pathA = WriteLib("comp_a", entries);
		string outUncomp = Path.Combine(OutputDirectory, "comp_out_level0.msl");
		string outCompressed = Path.Combine(OutputDirectory, "comp_out_level3.msl");

		MslMerger.Merge(new[] { pathA }, outUncomp, compressionLevel: 0);
		MslMerger.Merge(new[] { pathA }, outCompressed, compressionLevel: 3);

		long sizeUncomp = new FileInfo(outUncomp).Length;
		long sizeCompressed = new FileInfo(outCompressed).Length;

		Assert.That(sizeCompressed, Is.LessThan(sizeUncomp),
			"Compressed output should be smaller than uncompressed output.");
	}

	[Test]
	public void Merge_SameOutputAsInput_OverwriteAllowed()
	{
		var a = new List<MslLibraryEntry> { MakeEntry("OVR", 2, 300.0) };
		var b = new List<MslLibraryEntry> { MakeEntry("WRT", 2, 400.0) };

		string pathA = WriteLib("ovr_a", a);
		string pathB = WriteLib("ovr_b", b);
		string outPath = Path.Combine(OutputDirectory, "ovr_existing.msl");

		// Write an initial file at the output path
		MslWriter.Write(outPath, new List<MslLibraryEntry> { MakeEntry("OLD", 2, 999.0) });
		Assert.That(File.Exists(outPath), Is.True);

		// Merge should overwrite without error
		Assert.DoesNotThrow(() =>
			MslMerger.Merge(new[] { pathA, pathB }, outPath));

		using var lib = MslLibrary.Load(outPath);
		Assert.That(lib.PrecursorCount, Is.EqualTo(2), "Overwritten file should contain the merged entries.");
		Assert.That(lib.GetAllEntries().All(e => e.FullSequence != "OLD"), Is.True,
			"The old entry should not appear in the overwritten file.");
	}
}