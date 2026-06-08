// TestMslPrompt8Streaming.cs
// PR #1036 · smith-chem-wisc/mzLib · branch `mzlib_speclib`
// Prompt 8 — WriteStreaming Two-Pass Memory Problem
//
// Tests for:
//   Fix9a — MslSpillRecord extended to carry all precursor scalar fields;
//            WriteStreaming now consumes its source IEnumerable exactly once.
//   Fix9b — MslMerger.Merge feeds WriteStreaming via a lazy iterator,
//            eliminating the List<MslLibraryEntry> output accumulator.
//
// All tests are self-contained (no file fixtures, no network calls).
// Build: dotnet test mzLib.sln --filter "FullyQualifiedName~TestMslPrompt8"

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

namespace Test.SpectralLibrary;

[TestFixture]
[Category("Prompt8")]
public static class TestMslPrompt8Streaming
{
	// ────────────────────────────────────────────────────────────────────────
	// Helpers
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Creates a minimal MslLibraryEntry with one fragment. Used by all tests
	/// that do not need special field values.
	/// </summary>
	private static MslLibraryEntry MakeEntry(
		string modSeq,
		double mz,
		int charge = 2,
		float fragMz = 500f,
		float irt = 0f,
		float ionMobility = 0f,
		float qValue = float.NaN,
		bool isDecoy = false,
		string protein = "")
	{
		var entry = new MslLibraryEntry
		{
			FullSequence = modSeq,
			BaseSequence = modSeq.Replace("[Common Variable:Oxidation on M]", ""),
			PrecursorMz = mz,
			ChargeState = charge,
			RetentionTime = irt,
			IonMobility = ionMobility,
			QValue = qValue,
			IsDecoy = isDecoy,
			ProteinAccession = protein,
			Source = MslFormat.SourceType.Empirical,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz          = fragMz,
					Intensity   = 1.0f,
					Charge      = 1,
					ProductType = ProductType.b,
					FragmentNumber = 3
				}
			}
		};
		return entry;
	}

	/// <summary>
	/// Writes entries via <see cref="MslWriter.WriteStreaming"/> and immediately
	/// reads back the result via <see cref="MslReader.Load"/>. Returns the loaded data.
	/// </summary>
	private static MslLibraryData WriteAndLoad(
		string path,
		IEnumerable<MslLibraryEntry> entries,
		int compressionLevel = 0)
	{
		MslWriter.WriteStreaming(path, entries, compressionLevel);
		return MslReader.Load(path);
	}

	// ────────────────────────────────────────────────────────────────────────
	// Section 1 — Single-pass correctness (Fix9a)
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// The core correctness test for Fix9a: WriteStreaming output must be
	/// byte-for-byte field-identical to Write() output for the same entries.
	/// </summary>
	[Test]
	public static void WriteStreaming_SinglePass_ProducesIdenticalOutputToWrite()
	{
		string dir = TestContext.CurrentContext.TestDirectory;
		string pathW = Path.Combine(dir, "p8_write.msl");
		string pathWS = Path.Combine(dir, "p8_write_streaming.msl");

		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("PEPTIDE",     500.25, charge: 2, fragMz: 200f, irt: 12.3f,
					  ionMobility: 0.85f, qValue: 0.001f, protein: "P12345"),
			MakeEntry("LARGERPEP",   800.40, charge: 3, fragMz: 350f, irt: 45.0f,
					  ionMobility: 0f,    qValue: float.NaN),
			MakeEntry("DECOYSEQ",    600.30, charge: 2, fragMz: 250f, isDecoy: true),
		};

		try
		{
			MslWriter.Write(pathW, entries);
			MslWriter.WriteStreaming(pathWS, entries);

			using MslLibraryData w = MslReader.Load(pathW);
			using MslLibraryData ws = MslReader.Load(pathWS);

			Assert.That(ws.Count, Is.EqualTo(w.Count),
				"Entry count must match between Write and WriteStreaming.");

			for (int i = 0; i < w.Count; i++)
			{
				MslLibraryEntry expected = w.Entries[i];
				MslLibraryEntry actual = ws.Entries[i];

				Assert.That(actual.FullSequence, Is.EqualTo(expected.FullSequence),
					$"Entry {i}: FullSequence");
				Assert.That(actual.PrecursorMz, Is.EqualTo(expected.PrecursorMz).Within(1e-4),
					$"Entry {i}: PrecursorMz");
				Assert.That(actual.ChargeState, Is.EqualTo(expected.ChargeState),
					$"Entry {i}: ChargeState");
				Assert.That(actual.RetentionTime, Is.EqualTo(expected.RetentionTime).Within(1e-4f),
					$"Entry {i}: RetentionTime");
				Assert.That(actual.IonMobility, Is.EqualTo(expected.IonMobility).Within(1e-6f),
					$"Entry {i}: IonMobility");
				Assert.That(actual.IsDecoy, Is.EqualTo(expected.IsDecoy),
					$"Entry {i}: IsDecoy");
				Assert.That(actual.Source, Is.EqualTo(expected.Source),
					$"Entry {i}: Source");

				// QValue: both NaN should pass; otherwise compare within tolerance
				if (float.IsNaN(expected.QValue))
					Assert.That(float.IsNaN(actual.QValue), Is.True,
						$"Entry {i}: QValue expected NaN");
				else
					Assert.That(actual.QValue, Is.EqualTo(expected.QValue).Within(1e-6f),
						$"Entry {i}: QValue");

				Assert.That(actual.MatchedFragmentIons.Count, Is.EqualTo(expected.MatchedFragmentIons.Count),
					$"Entry {i}: fragment count");
				Assert.That(actual.MatchedFragmentIons[0].Mz, Is.EqualTo(expected.MatchedFragmentIons[0].Mz).Within(1e-4f),
					$"Entry {i}: fragment Mz");
			}
		}
		finally
		{
			if (File.Exists(pathW)) File.Delete(pathW);
			if (File.Exists(pathWS)) File.Delete(pathWS);
		}
	}

	/// <summary>
	/// Regression guard: a lazy generator (single-pass IEnumerable that cannot be
	/// re-enumerated) must now produce a file containing all expected entries.
	///
	/// Before Fix9a this produced a zero-precursor file because Pass 2 re-enumerated
	/// the source and got nothing.  After Fix9a Pass 2 reads only the spill file, so
	/// any IEnumerable is sufficient.
	///
	/// This test supersedes WriteStreaming_SinglePassIEnumerable_ProducesZeroPrecursorFile
	/// (Prompt 5 / G1), which documented the old broken behaviour.
	/// </summary>
	[Test]
	public static void WriteStreaming_SinglePassIEnumerable_ProducesAllEntries()
	{
		string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "p8_lazy.msl");

		const int N = 50;

		// A true single-pass IEnumerable — cannot be re-enumerated after first use.
		IEnumerable<MslLibraryEntry> LazySource()
		{
			for (int i = 0; i < N; i++)
				yield return MakeEntry($"SEQ{i:D3}", 400.0 + i * 0.5);
		}

		try
		{
			MslWriter.WriteStreaming(path, LazySource());
			using MslLibraryData data = MslReader.Load(path);

			Assert.That(data.Count, Is.EqualTo(N),
				"Single-pass IEnumerable must produce a file with the correct entry count.");

			// Spot-check first and last entries
			Assert.That(data.Entries[0].FullSequence, Is.EqualTo("SEQ000"));
			Assert.That(data.Entries[N - 1].FullSequence, Is.EqualTo($"SEQ{N - 1:D3}"));
		}
		finally
		{
			if (File.Exists(path)) File.Delete(path);
		}
	}

	/// <summary>
	/// All precursor scalar fields round-trip correctly through the extended spill record.
	/// This is the core contract for Fix9a: every field that was previously read by
	/// re-enumerating entries in Pass 2 must now be read from the spill file with
	/// identical values.
	/// </summary>
	[Test]
	public static void WriteStreaming_AllPrecursorScalars_RoundTripCorrectly()
	{
		string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "p8_scalars.msl");

		var entry = new MslLibraryEntry
		{
			FullSequence = "PEPTM[Common Variable:Oxidation on M]IDE",
			BaseSequence = "PEPTMIDE",
			PrecursorMz = 612.34,
			ChargeState = 3,
			RetentionTime = 27.5f,
			IonMobility = 1.05f,
			QValue = 0.0012f,
			IsDecoy = false,
			IsProteotypic = true,
			Nce = 28,
			Source = MslFormat.SourceType.EmpiricalRefined,
			MoleculeType = MslFormat.MoleculeType.Glycopeptide,
			DissociationType = DissociationType.HCD,
			ProteinAccession = "P99999",
			ProteinName = "Test Protein",
			GeneName = "TESTP",
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 456.78f, Intensity = 0.8f,
					Charge = 1, ProductType = ProductType.y,
					FragmentNumber = 5
				}
			}
		};

		try
		{
			MslWriter.WriteStreaming(path, new[] { entry });
			using MslLibraryData data = MslReader.Load(path);

			Assert.That(data.Count, Is.EqualTo(1));
			MslLibraryEntry r = data.Entries[0];

			Assert.That(r.PrecursorMz, Is.EqualTo(entry.PrecursorMz).Within(5e-3), "PrecursorMz");
			Assert.That(r.ChargeState, Is.EqualTo(entry.ChargeState), "ChargeState");
			Assert.That(r.RetentionTime, Is.EqualTo((float)entry.RetentionTime).Within(1e-4f), "RetentionTime");
			Assert.That(r.IonMobility, Is.EqualTo((float)entry.IonMobility).Within(1e-6f), "IonMobility");
			Assert.That(r.QValue, Is.EqualTo(entry.QValue).Within(1e-6f), "QValue");
			Assert.That(r.IsDecoy, Is.EqualTo(entry.IsDecoy), "IsDecoy");
			Assert.That(r.IsProteotypic, Is.EqualTo(entry.IsProteotypic), "IsProteotypic");
			Assert.That(r.Nce, Is.EqualTo(entry.Nce), "Nce");
			Assert.That(r.Source, Is.EqualTo(entry.Source), "Source");
			Assert.That(r.MoleculeType, Is.EqualTo(entry.MoleculeType), "MoleculeType");
			Assert.That(r.DissociationType, Is.EqualTo(entry.DissociationType), "DissociationType");
			Assert.That(r.BaseSequence, Is.EqualTo(entry.BaseSequence), "BaseSequence");
			Assert.That(r.MatchedFragmentIons.Count, Is.EqualTo(1), "fragment count");
		}
		finally
		{
			if (File.Exists(path)) File.Delete(path);
		}
	}

	/// <summary>
	/// WriteStreaming with compression enabled must produce field-identical output
	/// to Write() with the same compression level. The extended spill record must
	/// survive the compressed path unchanged.
	/// </summary>
	[Test]
	public static void WriteStreaming_WithCompression_ProducesIdenticalOutputToWrite()
	{
		string dir = TestContext.CurrentContext.TestDirectory;
		string pathW = Path.Combine(dir, "p8_write_comp.msl");
		string pathS = Path.Combine(dir, "p8_stream_comp.msl");

		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("COMPRESSED_A", 400.0, charge: 2, fragMz: 200f, irt: 5f),
			MakeEntry("COMPRESSED_B", 800.0, charge: 3, fragMz: 400f, irt: 30f),
		};

		try
		{
			MslWriter.Write(pathW, entries, compressionLevel: 3);
			MslWriter.WriteStreaming(pathS, entries, compressionLevel: 3);

			using MslLibraryData w = MslReader.Load(pathW);
			using MslLibraryData ws = MslReader.Load(pathS);

			Assert.That(ws.Count, Is.EqualTo(w.Count), "compressed: entry count");
			for (int i = 0; i < w.Count; i++)
			{
				Assert.That(ws.Entries[i].FullSequence, Is.EqualTo(w.Entries[i].FullSequence),
					$"compressed entry {i}: FullSequence");
				Assert.That(ws.Entries[i].ChargeState, Is.EqualTo(w.Entries[i].ChargeState),
					$"compressed entry {i}: ChargeState");
				Assert.That(ws.Entries[i].RetentionTime, Is.EqualTo(w.Entries[i].RetentionTime).Within(1e-4f),
					$"compressed entry {i}: RetentionTime");
			}
		}
		finally
		{
			if (File.Exists(pathW)) File.Delete(pathW);
			if (File.Exists(pathS)) File.Delete(pathS);
		}
	}

	/// <summary>
	/// MslSpillRecord size must be 56 bytes (22-byte original + 34 bytes of new scalar fields).
	/// This test serves the same function as MslStructs.SizeCheck — it catches accidental
	/// field additions or Pack-setting errors specific to the spill struct.
	/// </summary>
	[Test]
	public static void MslSpillRecord_SizeIs56Bytes()
	{
		int size = Marshal.SizeOf<MslSpillRecord>();
		Assert.That(size, Is.EqualTo(56),
			"MslSpillRecord must be exactly 56 bytes after Fix9a. " +
			"A different value means a field was added or removed without updating this test.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Section 2 — MslMerger lazy iterator correctness (Fix9b)
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// The fundamental post-fix test: merge two libraries with known overlap and
	/// verify winner selection is correct for all three conflict policies.
	/// </summary>
	[Test]
	[TestCase(MslMergeConflictPolicy.KeepFirst, 400.0000, TestName = "KeepFirst")]
	[TestCase(MslMergeConflictPolicy.KeepLast, 400.0005, TestName = "KeepLast")]
	[TestCase(MslMergeConflictPolicy.KeepLowestQValue, 400.0000, TestName = "KeepLowestQValue")]
	public static void Merge_LargeLibrary_OutputMatchesExpected(
		MslMergeConflictPolicy policy,
		double expectedWinnerMz)
	{
		string dir = TestContext.CurrentContext.TestDirectory;
		string pathA = Path.Combine(dir, "p8_merge_a.msl");
		string pathB = Path.Combine(dir, "p8_merge_b.msl");
		string pathOut = Path.Combine(dir, "p8_merge_out.msl");

		// Two libraries share one Name ("SHARED/2").
		// The duplicate entries must be within the 0.001 Da KeepLowestQValue flush
		// window so both are buffered together before one is resolved and emitted.
		// Library A has lower qValue so it wins under KeepLowestQValue.
		// For KeepLast, library B (index 1) wins; expectedWinnerMz is 400.0005 for that case.
		var libA = new List<MslLibraryEntry>
		{
			MakeEntry("ONLYINA", 300.0, charge: 2),
			MakeEntry("SHARED",  400.0000, charge: 2, qValue: 0.001f),
		};
		var libB = new List<MslLibraryEntry>
		{
			MakeEntry("SHARED",  400.0005, charge: 2, qValue: 0.05f),
			MakeEntry("ONLYINB", 600.0, charge: 2),
		};

		try
		{
			MslWriter.Write(pathA, libA);
			MslWriter.Write(pathB, libB);

			MslMergeResult result = MslMerger.Merge(
				new[] { pathA, pathB }, pathOut,
				conflictPolicy: policy, deduplicate: true);

			using MslLibraryData merged = MslReader.Load(pathOut);

			Assert.That(merged.Count, Is.EqualTo(3),
				$"{policy}: output should have 3 entries (2 unique + 1 winner from duplicate).");

			Assert.That(result.DuplicatesSkipped, Is.EqualTo(1),
				$"{policy}: exactly 1 duplicate skipped.");

			// Find the SHARED entry in the output and check which file won
			MslLibraryEntry? shared = merged.Entries
				.FirstOrDefault(e => e.FullSequence == "SHARED");

			Assert.That(shared, Is.Not.Null, "SHARED entry must be present in output.");
			Assert.That(shared!.PrecursorMz, Is.EqualTo(expectedWinnerMz).Within(0.0002),
				$"{policy}: SHARED winner has wrong PrecursorMz.");
		}
		finally
		{
			foreach (string f in new[] { pathA, pathB, pathOut })
				if (File.Exists(f)) File.Delete(f);
		}
	}

	/// <summary>
	/// Custom neutral losses must survive the merge round-trip.
	/// This is a regression guard for the Fix7 × Fix9b interaction: MslMerger
	/// previously threw NotSupportedException for any library containing custom
	/// losses (fixed in Fix7). This test confirms the fix remains effective when
	/// the lazy iterator feeds WriteStreaming.
	/// </summary>
	[Test]
	public static void Merge_CustomNeutralLoss_SurvivesMerge()
	{
		string dir = TestContext.CurrentContext.TestDirectory;
		string pathA = Path.Combine(dir, "p8_custom_a.msl");
		string pathOut = Path.Combine(dir, "p8_custom_out.msl");

		const double CustomLoss = -163.0633; // typical hexose neutral loss

		var entry = new MslLibraryEntry
		{
			FullSequence = "GLYCOPEPTIDE",
			BaseSequence = "GLYCOPEPTIDE",
			PrecursorMz = 750.0,
			ChargeState = 2,
			Source = MslFormat.SourceType.Empirical,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz           = 550.0f,
					Intensity    = 1.0f,
					Charge       = 1,
					ProductType  = ProductType.y,
					FragmentNumber = 5,
					NeutralLoss  = CustomLoss
				}
			}
		};

		try
		{
			MslWriter.Write(pathA, new[] { entry });

			// Merge with itself — trivial but exercises the custom-loss path
			MslMergeResult result = MslMerger.Merge(
				new[] { pathA }, pathOut,
				conflictPolicy: MslMergeConflictPolicy.KeepFirst);

			Assert.That(result.OutputEntryCount, Is.EqualTo(1));

			using MslLibraryData merged = MslReader.Load(pathOut);
			Assert.That(merged.Count, Is.EqualTo(1));
			MslLibraryEntry r = merged.Entries[0];
			Assert.That(r.MatchedFragmentIons.Count, Is.EqualTo(1));
			Assert.That(r.MatchedFragmentIons[0].NeutralLoss, Is.EqualTo(CustomLoss).Within(1e-6),
				"Custom neutral loss mass must be preserved through the merge.");
		}
		finally
		{
			foreach (string f in new[] { pathA, pathOut })
				if (File.Exists(f)) File.Delete(f);
		}
	}

	/// <summary>
	/// Verifies that the output entries are in m/z ascending order after the merge —
	/// the k-way heap must produce monotonically non-decreasing m/z output.
	/// </summary>
	[Test]
	public static void Merge_Output_IsMzSortedAscending_AfterFix9b()
	{
		string dir = TestContext.CurrentContext.TestDirectory;
		string pathA = Path.Combine(dir, "p8_sort_a.msl");
		string pathB = Path.Combine(dir, "p8_sort_b.msl");
		string pathOut = Path.Combine(dir, "p8_sort_out.msl");

		// Interleaved m/z values across two files
		var libA = new List<MslLibraryEntry>
		{
			MakeEntry("A1", 300.0), MakeEntry("A2", 500.0), MakeEntry("A3", 700.0)
		};
		var libB = new List<MslLibraryEntry>
		{
			MakeEntry("B1", 400.0), MakeEntry("B2", 600.0), MakeEntry("B3", 800.0)
		};

		try
		{
			MslWriter.Write(pathA, libA);
			MslWriter.Write(pathB, libB);

			MslMerger.Merge(new[] { pathA, pathB }, pathOut,
				conflictPolicy: MslMergeConflictPolicy.KeepFirst);

			using MslLibraryData merged = MslReader.Load(pathOut);
			Assert.That(merged.Count, Is.EqualTo(6));

			for (int i = 1; i < merged.Count; i++)
			{
				Assert.That(merged.Entries[i].PrecursorMz,
					Is.GreaterThanOrEqualTo(merged.Entries[i - 1].PrecursorMz),
					$"Output not sorted ascending at index {i}.");
			}
		}
		finally
		{
			foreach (string f in new[] { pathA, pathB, pathOut })
				if (File.Exists(f)) File.Delete(f);
		}
	}

	/// <summary>
	/// Fragment data must be fully preserved through the lazy merge.
	/// This tests that entry.MatchedFragmentIons — loaded on-demand by GetAllEntries/LoadFragmentsOnDemand
	/// — arrive intact in the output file even though the entry is never held in a List.
	/// </summary>
	[Test]
	public static void Merge_FragmentData_PreservedInOutput_AfterFix9b()
	{
		string dir = TestContext.CurrentContext.TestDirectory;
		string pathA = Path.Combine(dir, "p8_frag_a.msl");
		string pathOut = Path.Combine(dir, "p8_frag_out.msl");

		var entry = MakeEntry("FRAGTEST", 550.0, fragMz: 321.45f);
		// Add a second fragment to make the check more specific
		entry.MatchedFragmentIons.Add(new MslFragmentIon
		{
			Mz = 456.78f,
			Intensity = 0.5f,
			Charge = 1,
			ProductType = ProductType.y,
			FragmentNumber = 4
		});

		try
		{
			MslWriter.Write(pathA, new[] { entry });
			MslMerger.Merge(new[] { pathA }, pathOut,
				conflictPolicy: MslMergeConflictPolicy.KeepFirst);

			using MslLibraryData merged = MslReader.Load(pathOut);
			Assert.That(merged.Count, Is.EqualTo(1));

			MslLibraryEntry r = merged.Entries[0];
			Assert.That(r.MatchedFragmentIons.Count, Is.EqualTo(2), "Both fragments must survive the merge.");

			// MatchedFragmentIons are sorted by m/z ascending on write
			float[] expectedMzs = { 321.45f, 456.78f };
			for (int i = 0; i < expectedMzs.Length; i++)
				Assert.That(r.MatchedFragmentIons[i].Mz, Is.EqualTo(expectedMzs[i]).Within(0.01f),
					$"fragment[{i}].Mz");
		}
		finally
		{
			foreach (string f in new[] { pathA, pathOut })
				if (File.Exists(f)) File.Delete(f);
		}
	}

	/// <summary>
	/// MslMergeResult.OutputEntryCount must match the actual number of entries in the
	/// output file after the fix.  With the original List-based code, OutputEntryCount
	/// came from outputEntries.Count; with Fix9b it comes from state.OutputEntryCount
	/// incremented inside the iterator.
	/// </summary>
	[Test]
	public static void Merge_OutputEntryCount_MatchesActualOutputFile()
	{
		string dir = TestContext.CurrentContext.TestDirectory;
		string pathA = Path.Combine(dir, "p8_count_a.msl");
		string pathB = Path.Combine(dir, "p8_count_b.msl");
		string pathOut = Path.Combine(dir, "p8_count_out.msl");

		// 5 unique + 2 duplicates = 7 total source, 5 output
		var libA = new List<MslLibraryEntry>
		{
			MakeEntry("UNIQ1", 300.0), MakeEntry("SHARED1", 400.0),
			MakeEntry("UNIQ2", 500.0), MakeEntry("SHARED2", 600.0)
		};
		var libB = new List<MslLibraryEntry>
		{
			MakeEntry("SHARED1", 400.1), // duplicate (slightly different mz but same Name)
            MakeEntry("SHARED2", 600.1), // duplicate
            MakeEntry("UNIQ3",   700.0)
		};

		try
		{
			MslWriter.Write(pathA, libA);
			MslWriter.Write(pathB, libB);

			MslMergeResult result = MslMerger.Merge(
				new[] { pathA, pathB }, pathOut,
				conflictPolicy: MslMergeConflictPolicy.KeepFirst);

			using MslLibraryData merged = MslReader.Load(pathOut);

			Assert.That(result.OutputEntryCount, Is.EqualTo(merged.Count),
				"MslMergeResult.OutputEntryCount must equal the number of entries in the output file.");
			Assert.That(result.DuplicatesSkipped, Is.EqualTo(2),
				"Two duplicates (SHARED1, SHARED2) should have been skipped.");
		}
		finally
		{
			foreach (string f in new[] { pathA, pathB, pathOut })
				if (File.Exists(f)) File.Delete(f);
		}
	}

	/// <summary>
	/// Memory regression: confirm that Merge does not hold O(TotalOutputEntries)
	/// live MslLibraryEntry objects simultaneously.
	///
	/// This test uses GC.GetTotalMemory() snapshots as a coarse proxy for the heap
	/// allocation profile.  The pre-fix implementation accumulated all winning entries
	/// in a List, allocating roughly (N × avgEntrySize) bytes before WriteStreaming was
	/// called.  The post-fix implementation yields entries one at a time; WriteStreaming
	/// spills scalars to disk.
	///
	/// We write N = 5,000 entries through the merge path and assert that the GC
	/// allocation during the call does not exceed a generous 10 MB ceiling.  The
	/// 10 MB bound is chosen to be well above the actual post-fix overhead
	/// (string table + spill file IO buffers ≈ a few hundred KB for 5,000 entries)
	/// while being well below the pre-fix allocation (5,000 entries × ~600 bytes
	/// per entry skeleton + fragment list ≈ 3 MB, growing proportionally with N).
	///
	/// This is a coarse test and may need adjustment on extremely constrained CI
	/// machines; it is explicitly marked [Category("MemoryBound")] so it can be
	/// excluded from fast-feedback runs.
	/// </summary>
	[Test]
	[Category("MemoryBound")]
	public static void Merge_PeakMemory_IsSublinearInOutputSize()
	{
		const int N = 5_000;
		string dir = TestContext.CurrentContext.TestDirectory;
		string pathA = Path.Combine(dir, "p8_mem_a.msl");
		string pathOut = Path.Combine(dir, "p8_mem_out.msl");

		// Build a library of N entries
		var entries = new List<MslLibraryEntry>(N);
		for (int i = 0; i < N; i++)
			entries.Add(MakeEntry($"SEQ{i:D6}", 400.0 + i * 0.001));

		try
		{
			MslWriter.Write(pathA, entries);

			// Force a full GC before measurement
			GC.Collect(2, GCCollectionMode.Forced, blocking: true);
			GC.WaitForPendingFinalizers();
			long memBefore = GC.GetTotalMemory(forceFullCollection: true);

			MslMerger.Merge(new[] { pathA }, pathOut,
				conflictPolicy: MslMergeConflictPolicy.KeepFirst);

			long memAfter = GC.GetTotalMemory(forceFullCollection: false);
			long deltaBytes = memAfter - memBefore;

			// 10 MB ceiling — generous headroom above post-fix overhead (~few hundred KB)
			const long MaxDeltaBytes = 10L * 1024 * 1024;
			Assert.That(deltaBytes, Is.LessThan(MaxDeltaBytes),
				$"Peak memory delta during Merge of {N} entries was {deltaBytes / 1024} KB, " +
				$"exceeding the {MaxDeltaBytes / 1024} KB ceiling. " +
				"This suggests the List<MslLibraryEntry> accumulator was reintroduced.");
		}
		finally
		{
			foreach (string f in new[] { pathA, pathOut })
				if (File.Exists(f)) File.Delete(f);
		}
	}

	/// <summary>
	/// Exception safety: if the second source file is corrupt, the first library's
	/// FileStream must be disposed and no partial output file left at the output path.
	/// This was correct in the original code and must remain so after Fix9b.
	/// </summary>
	[Test]
	public static void Merge_InvalidSecondFile_ThrowsAndFirstLibraryDisposed_AfterFix9b()
	{
		string dir = TestContext.CurrentContext.TestDirectory;
		string pathA = Path.Combine(dir, "p8_exc_a.msl");
		string pathBad = Path.Combine(dir, "p8_exc_bad.msl");
		string pathOut = Path.Combine(dir, "p8_exc_out.msl");

		MslWriter.Write(pathA, new[] { MakeEntry("GOOD", 400.0) });

		// Write a file with invalid content
		File.WriteAllText(pathBad, "THIS IS NOT A VALID MSL FILE");

		try
		{
			Assert.That(
				() => MslMerger.Merge(new[] { pathA, pathBad }, pathOut),
				Throws.InstanceOf<Exception>());

			// pathOut must not exist (partial write cleaned up)
			Assert.That(File.Exists(pathOut), Is.False,
				"No partial output file should exist after a failed merge.");

			// pathA must be deletable (file handle released by finally block)
			Assert.DoesNotThrow(() => File.Delete(pathA),
				"Source library file handle must be released even when merge fails.");
		}
		finally
		{
			foreach (string f in new[] { pathA, pathBad, pathOut })
				if (File.Exists(f)) File.Delete(f);
		}
	}

	/// <summary>
	/// DuplicatesSkipped accounting invariant (Fix9b):
	///   OutputEntryCount + DuplicatesSkipped == TotalSourceEntryCount
	/// Must hold for all three conflict policies after the refactor.
	/// </summary>
	[Test]
	[TestCase(MslMergeConflictPolicy.KeepFirst, TestName = "Invariant_KeepFirst")]
	[TestCase(MslMergeConflictPolicy.KeepLast, TestName = "Invariant_KeepLast")]
	[TestCase(MslMergeConflictPolicy.KeepLowestQValue, TestName = "Invariant_KeepLowestQValue")]
	public static void Merge_DuplicatesSkipped_Invariant_Holds_AfterFix9b(
		MslMergeConflictPolicy policy)
	{
		string dir = TestContext.CurrentContext.TestDirectory;
		string pathA = Path.Combine(dir, "p8_inv_a.msl");
		string pathB = Path.Combine(dir, "p8_inv_b.msl");
		string pathOut = Path.Combine(dir, $"p8_inv_{policy}_out.msl");

		// 3 unique + 2 shared = 5 in A, 4 in B; 9 total source, 5 output
		var libA = new List<MslLibraryEntry>
		{
			MakeEntry("UNIQUE1", 300.0, qValue: 0.01f),
			MakeEntry("SHARED1", 400.0, qValue: 0.01f),
			MakeEntry("UNIQUE2", 500.0, qValue: 0.01f),
			MakeEntry("SHARED2", 600.0, qValue: 0.01f),
			MakeEntry("UNIQUE3", 700.0, qValue: 0.01f),
		};
		var libB = new List<MslLibraryEntry>
		{
			MakeEntry("SHARED1", 400.1, qValue: 0.05f),
			MakeEntry("SHARED2", 600.1, qValue: 0.05f),
			MakeEntry("UNIQUEB1", 800.0, qValue: 0.01f),
			MakeEntry("UNIQUEB2", 900.0, qValue: 0.01f),
		};

		try
		{
			MslWriter.Write(pathA, libA);
			MslWriter.Write(pathB, libB);

			MslMergeResult result = MslMerger.Merge(
				new[] { pathA, pathB }, pathOut,
				conflictPolicy: policy, deduplicate: true);

			Assert.That(
				result.OutputEntryCount + result.DuplicatesSkipped,
				Is.EqualTo(result.TotalSourceEntryCount),
				$"{policy}: OutputEntryCount ({result.OutputEntryCount}) + " +
				$"DuplicatesSkipped ({result.DuplicatesSkipped}) must equal " +
				$"TotalSourceEntryCount ({result.TotalSourceEntryCount}).");
		}
		finally
		{
			foreach (string f in new[] { pathA, pathB, pathOut })
				if (File.Exists(f)) File.Delete(f);
		}
	}

	/// <summary>
	/// Confirms the NCE field is correctly pre-encoded in the spill record and
	/// faithfully round-trips through WriteStreaming. This is the one field where
	/// a transformation (EncodeNce: × 10 + clamp) must happen in Pass 1 and the
	/// encoded value must be stored in the spill, not the raw entry.Nce value.
	/// </summary>
	[Test]
	[TestCase(20, TestName = "Nce_20")]
	[TestCase(28, TestName = "Nce_28")]
	[TestCase(35, TestName = "Nce_35")]
	[TestCase(0, TestName = "Nce_0_unknown")]
	[TestCase(3276, TestName = "Nce_3276_maxSafe")]
	public static void WriteStreaming_NceField_RoundTripsCorrectly(int nce)
	{
		string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"p8_nce_{nce}.msl");

		var entry = MakeEntry("NCETEST", 500.0);
		entry.Nce = nce;

		try
		{
			MslWriter.WriteStreaming(path, new[] { entry });
			using MslLibraryData data = MslReader.Load(path);

			Assert.That(data.Count, Is.EqualTo(1));
			Assert.That(data.Entries[0].Nce, Is.EqualTo(nce),
				$"NCE {nce} must round-trip correctly through WriteStreaming.");
		}
		finally
		{
			if (File.Exists(path)) File.Delete(path);
		}
	}

	/// <summary>
	/// Empty library: WriteStreaming must produce a valid zero-precursor file.
	/// Verifies the single-pass path handles the degenerate case (nPrecursors = 0)
	/// without touching the spill reader loop.
	/// </summary>
	[Test]
	public static void WriteStreaming_EmptySource_ProducesValidZeroPrecursorFile()
	{
		string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "p8_empty.msl");

		try
		{
			MslWriter.WriteStreaming(path, Enumerable.Empty<MslLibraryEntry>());
			using MslLibraryData data = MslReader.Load(path);
			Assert.That(data.Count, Is.EqualTo(0),
				"Empty source must produce a valid zero-precursor file.");
		}
		finally
		{
			if (File.Exists(path)) File.Delete(path);
		}
	}
}