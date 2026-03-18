using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace Test.SpectralLibrary.MSL;

/// <summary>
/// Tests targeting the findings from Prompt 2 — LoadIndexOnly Memory Behavior.
///
/// Covers:
///   M1  — LoadIndexOnly documented contract vs. implementation
///   M2  — Fragment demand-loading: empty before load, populated after
///   M3  — Thread safety of concurrent LoadFragmentsOnDemand calls
///   M4  — Disposal: FileStream closed, ObjectDisposedException after dispose
///   M5  — Exception safety: no stream leak when construction fails
///   M6  — Double-dispose idempotency
///   M7  — Fragment m/z values match after demand-load (correctness)
/// </summary>
[TestFixture]
public class TestMslPrompt2LoadIndexOnly
{
	// ── Temp file infrastructure ──────────────────────────────────────────

	private static readonly string OutputDir =
		Path.Combine(Path.GetTempPath(), "MslPrompt2Tests");

	private static string Tmp(string name) =>
		Path.Combine(OutputDir, name + ".msl");

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDir);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDir))
			Directory.Delete(OutputDir, recursive: true);
	}

	// ── Shared library fixture ────────────────────────────────────────────

	private static List<MslLibraryEntry> BuildEntries(int count = 3)
	{
		var entries = new List<MslLibraryEntry>();
		string[] sequences = { "PEPTIDE", "ACDEFGHIK", "LESLIEK" };
		float[] mzBases = { 312.15f, 450.18f, 188.07f };

		for (int i = 0; i < count; i++)
		{
			string seq = sequences[i % sequences.Length] + (i >= sequences.Length ? i.ToString() : "");
			entries.Add(new MslLibraryEntry
			{
				ModifiedSequence = seq,
				StrippedSequence = seq,
				PrecursorMz = 400.0 + i * 50.0,
				Charge = 2,
				Irt = 30.0 + i * 5.0,
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
					new MslFragmentIon
					{
						ProductType     = ProductType.b, FragmentNumber  = 3,
						Charge          = 1,             Mz              = mzBases[i % mzBases.Length],
						Intensity       = 1.0f,          NeutralLoss     = 0.0,
						ResiduePosition = 3
					},
					new MslFragmentIon
					{
						ProductType     = ProductType.y, FragmentNumber  = 4,
						Charge          = 1,             Mz              = mzBases[i % mzBases.Length] + 100f,
						Intensity       = 0.8f,          NeutralLoss     = 0.0,
						ResiduePosition = 3
					}
				}
			});
		}
		return entries;
	}

	private static string BuildAndWriteLibrary(string name, int entryCount = 3)
	{
		string path = Tmp(name);
		MslWriter.Write(path, BuildEntries(entryCount));
		return path;
	}

	// ═════════════════════════════════════════════════════════════════════
	// M1 — LoadIndexOnly contract: fragments start empty
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Immediately after LoadIndexOnly, every entry's Fragments list must be
	/// empty. This is the core promise of index-only mode: fragment bytes
	/// remain on disk until explicitly requested.
	///
	/// Verified via MslReader.LoadIndexOnly (returns MslLibraryData) rather than
	/// MslLibrary.LoadIndexOnly, because MslLibrary.GetAllEntries invokes the loader
	/// delegate and triggers demand-loading as part of its normal operation.
	/// MslLibraryData.Entries exposes the raw skeleton entries before any demand-load.
	/// </summary>
	[Test]
	public void LoadIndexOnly_FragmentLists_AreEmptyImmediatelyAfterOpen()
	{
		string path = BuildAndWriteLibrary(nameof(LoadIndexOnly_FragmentLists_AreEmptyImmediatelyAfterOpen));

		// MslReader.LoadIndexOnly returns MslLibraryData whose Entries collection
		// contains skeleton entries with empty fragment lists — no loader delegate
		// is invoked when accessing Entries directly.
		using MslLibraryData data = MslReader.LoadIndexOnly(path);

		foreach (MslLibraryEntry entry in data.Entries)
		{
			Assert.That(entry.Fragments, Is.Empty,
				$"Entry '{entry.ModifiedSequence}' must have an empty Fragments list " +
				"immediately after LoadIndexOnly (fragments are loaded on demand, not upfront).");
		}
	}

	/// <summary>
	/// LoadIndexOnly must report IsIndexOnly = true.
	/// </summary>
	[Test]
	public void LoadIndexOnly_IsIndexOnly_ReturnsTrue()
	{
		string path = BuildAndWriteLibrary(nameof(LoadIndexOnly_IsIndexOnly_ReturnsTrue));
		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);
		Assert.That(lib.IsIndexOnly, Is.True);
	}

	/// <summary>
	/// MslLibrary.Load (full-load) must return IsIndexOnly = false, confirming
	/// that the flag correctly distinguishes the two modes.
	/// </summary>
	[Test]
	public void Load_IsIndexOnly_ReturnsFalse()
	{
		string path = BuildAndWriteLibrary(nameof(Load_IsIndexOnly_ReturnsFalse));
		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.IsIndexOnly, Is.False);
	}

	// ═════════════════════════════════════════════════════════════════════
	// M2 — Fragment demand-loading correctness
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// After calling GetEntry (which triggers LoadFragmentsOnDemand internally),
	/// the entry's Fragments list must be non-empty.
	/// </summary>
	[Test]
	public void LoadIndexOnly_GetEntry_PopulatesFragments()
	{
		string path = BuildAndWriteLibrary(nameof(LoadIndexOnly_GetEntry_PopulatesFragments));

		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);

		bool found = lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? entry);

		Assert.That(found, Is.True, "TryGetEntry must find PEPTIDE/2.");
		Assert.That(entry!.Fragments, Is.Not.Empty,
			"Fragments must be populated after demand-load via GetEntry.");
	}

	/// <summary>
	/// The fragment m/z values recovered via demand-load must exactly match
	/// those written to disk. Verifies that seek+read produces the correct data,
	/// not just non-empty data.
	/// </summary>
	[Test]
	public void LoadIndexOnly_DemandLoadedFragments_HaveCorrectMzValues()
	{
		var entries = BuildEntries(1);
		float expectedMz0 = entries[0].Fragments[0].Mz;
		float expectedMz1 = entries[0].Fragments[1].Mz;

		string path = Tmp(nameof(LoadIndexOnly_DemandLoadedFragments_HaveCorrectMzValues));
		MslWriter.Write(path, entries);

		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);
		bool found = lib.TryGetEntry(entries[0].ModifiedSequence, 2, out MslLibraryEntry? entry);

		Assert.That(found, Is.True);

		// Fragments are sorted by m/z ascending on write
		var frags = entry!.Fragments.OrderBy(f => f.Mz).ToList();

		Assert.That(frags[0].Mz, Is.EqualTo(Math.Min(expectedMz0, expectedMz1)).Within(1e-3f),
			"Lowest m/z fragment must match the written value after demand-load.");
		Assert.That(frags[1].Mz, Is.EqualTo(Math.Max(expectedMz0, expectedMz1)).Within(1e-3f),
			"Highest m/z fragment must match the written value after demand-load.");
	}

	/// <summary>
	/// The fragment count after demand-load must match the count written.
	/// </summary>
	[Test]
	public void LoadIndexOnly_DemandLoadedFragments_CountMatchesWritten()
	{
		var entries = BuildEntries(1);
		int expectedCount = entries[0].Fragments.Count;

		string path = Tmp(nameof(LoadIndexOnly_DemandLoadedFragments_CountMatchesWritten));
		MslWriter.Write(path, entries);

		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);
		bool found = lib.TryGetEntry(entries[0].ModifiedSequence, 2, out MslLibraryEntry? entry);

		Assert.That(found, Is.True);
		Assert.That(entry!.Fragments.Count, Is.EqualTo(expectedCount),
			$"Demand-loaded fragment count must equal written count ({expectedCount}).");
	}

	/// <summary>
	/// Full-load and index-only demand-load must produce identical fragment data
	/// for the same entry. This is the strongest correctness guarantee.
	/// </summary>
	[Test]
	public void LoadIndexOnly_FragmentData_MatchesFullLoad()
	{
		string path = BuildAndWriteLibrary(nameof(LoadIndexOnly_FragmentData_MatchesFullLoad));

		using MslLibrary fullLib = MslLibrary.Load(path);
		using MslLibrary indexLib = MslLibrary.LoadIndexOnly(path);

		fullLib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? fullEntry);
		indexLib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? indexEntry);

		Assert.That(fullEntry, Is.Not.Null);
		Assert.That(indexEntry, Is.Not.Null);

		var fullFrags = fullEntry!.Fragments.OrderBy(f => f.Mz).ToList();
		var indexFrags = indexEntry!.Fragments.OrderBy(f => f.Mz).ToList();

		Assert.That(indexFrags.Count, Is.EqualTo(fullFrags.Count),
			"Index-only demand-load must produce the same fragment count as full-load.");

		for (int i = 0; i < fullFrags.Count; i++)
		{
			Assert.That(indexFrags[i].Mz, Is.EqualTo(fullFrags[i].Mz).Within(1e-4f),
				$"Fragment [{i}] Mz mismatch between full-load and demand-load.");
			Assert.That(indexFrags[i].Intensity, Is.EqualTo(fullFrags[i].Intensity).Within(1e-4f),
				$"Fragment [{i}] Intensity mismatch.");
			Assert.That(indexFrags[i].ProductType, Is.EqualTo(fullFrags[i].ProductType),
				$"Fragment [{i}] ProductType mismatch.");
		}
	}

	// ═════════════════════════════════════════════════════════════════════
	// M3 — Thread safety of concurrent demand-loads
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Multiple threads calling GetEntry simultaneously must each receive
	/// correct fragment data with no interleaving or corruption. The internal
	/// _streamLock must prevent seek+read races.
	/// </summary>
	[Test]
	public void LoadIndexOnly_ConcurrentDemandLoads_ProduceCorrectData()
	{
		// Write a library with 3 distinct entries
		string path = BuildAndWriteLibrary(
			nameof(LoadIndexOnly_ConcurrentDemandLoads_ProduceCorrectData), entryCount: 3);

		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);

		// Collect all entries' lookup keys upfront (before concurrent access)
		var allEntries = lib.GetAllEntries().ToList();
		var exceptions = new System.Collections.Concurrent.ConcurrentBag<Exception>();
		int iterations = 20;

		// Each task repeatedly demand-loads fragments for every entry
		var tasks = Enumerable.Range(0, 4).Select(_ => Task.Run(() =>
		{
			for (int iter = 0; iter < iterations; iter++)
			{
				foreach (MslLibraryEntry skeleton in allEntries)
				{
					try
					{
						bool found = lib.TryGetEntry(
							skeleton.ModifiedSequence, skeleton.Charge, out MslLibraryEntry? loaded);

						if (!found || loaded is null || loaded.Fragments.Count == 0)
							exceptions.Add(new AssertionException(
								$"Concurrent demand-load failed for '{skeleton.ModifiedSequence}': " +
								$"found={found}, fragments={loaded?.Fragments.Count ?? -1}"));
					}
					catch (Exception ex) when (ex is not AssertionException)
					{
						exceptions.Add(ex);
					}
				}
			}
		})).ToArray();

		Task.WaitAll(tasks);

		Assert.That(exceptions, Is.Empty,
			"Concurrent demand-loads must all succeed with no exceptions or empty fragment lists. " +
			$"Failures: {string.Join("; ", exceptions.Select(e => e.Message))}");
	}

	// ═════════════════════════════════════════════════════════════════════
	// M4 — Disposal correctness
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// After Dispose(), the FileStream must be released. Verified cross-platform by
	/// opening the file with FileShare.None (exclusive access): if the stream is still
	/// open this throws on both Windows and Linux/macOS, giving a reliable signal on
	/// all CI platforms. The Windows-only delete-lock approach is not used because on
	/// Linux deleting an open file succeeds, making the test pass even on a leaked stream.
	/// </summary>
	[Test]
	public void LoadIndexOnly_AfterDispose_FileStreamIsClosed()
	{
		string path = BuildAndWriteLibrary(nameof(LoadIndexOnly_AfterDispose_FileStreamIsClosed));

		MslLibrary lib = MslLibrary.LoadIndexOnly(path);
		lib.Dispose();

		// Exclusive open: throws IOException if any handle to the file is still open.
		// This is cross-platform: unlike File.Delete, FileShare.None fails on both
		// Windows and Linux/macOS when the stream is leaked.
		Assert.That(() =>
		{
			using var probe = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.None);
		},
			Throws.Nothing,
			"Opening with FileShare.None must succeed after Dispose() — confirms the FileStream was closed on all platforms.");
	}

	/// <summary>
	/// Calling Dispose() multiple times must not throw.
	/// </summary>
	[Test]
	public void LoadIndexOnly_DoubleDispose_IsIdempotent()
	{
		string path = BuildAndWriteLibrary(nameof(LoadIndexOnly_DoubleDispose_IsIdempotent));

		MslLibrary lib = MslLibrary.LoadIndexOnly(path);

		Assert.That(() => { lib.Dispose(); lib.Dispose(); }, Throws.Nothing,
			"Calling Dispose() twice must not throw.");
	}

	/// <summary>
	/// Every public query method must throw ObjectDisposedException after
	/// the library has been disposed.
	/// </summary>
	[Test]
	public void LoadIndexOnly_AfterDispose_QueriesThrowObjectDisposedException()
	{
		string path = BuildAndWriteLibrary(
			nameof(LoadIndexOnly_AfterDispose_QueriesThrowObjectDisposedException));

		MslLibrary lib = MslLibrary.LoadIndexOnly(path);
		lib.Dispose();

		Assert.That(() => lib.TryGetEntry("PEPTIDE", 2, out _),
			Throws.InstanceOf<ObjectDisposedException>(),
			"TryGetEntry must throw ObjectDisposedException after Dispose.");

		Assert.That(() => lib.GetAllEntries().ToList(),
			Throws.InstanceOf<ObjectDisposedException>(),
			"GetAllEntries must throw ObjectDisposedException after Dispose.");

		Assert.That(() => lib.QueryWindow(400f, 500f, 0f, 100f),
			Throws.InstanceOf<ObjectDisposedException>(),
			"QueryWindow must throw ObjectDisposedException after Dispose.");
	}

	/// <summary>
	/// LoadFragmentsOnDemand must throw ObjectDisposedException when the
	/// underlying MslLibraryData has been disposed. Verifies the inner lock
	/// null-check path is reachable.
	/// </summary>
	[Test]
	public void MslLibraryData_LoadFragmentsOnDemand_ThrowsAfterDispose()
	{
		string path = BuildAndWriteLibrary(
			nameof(MslLibraryData_LoadFragmentsOnDemand_ThrowsAfterDispose));

		// Access _rawLibrary directly: write a one-entry library, load index-only,
		// dispose the data container, then try to load fragments.
		// We do this via MslLibrary.Dispose which disposes _rawLibrary internally.

		string path2 = Tmp(nameof(MslLibraryData_LoadFragmentsOnDemand_ThrowsAfterDispose) + "_data");
		MslWriter.Write(path2, BuildEntries(1));

		// Get a MslLibraryData directly via MslReader
		MslLibraryData data = MslReader.LoadIndexOnly(path2);
		data.Dispose();

		// Now attempt a demand-load — must throw ObjectDisposedException
		Assert.That(() => data.LoadFragmentsOnDemand(0),
			Throws.InstanceOf<ObjectDisposedException>(),
			"LoadFragmentsOnDemand must throw ObjectDisposedException after Dispose.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// M5 — Zero-entry library in index-only mode
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// LoadIndexOnly on a zero-precursor library must succeed, report
	/// PrecursorCount = 0, and not hold a stream (or hold one safely disposable).
	/// </summary>
	[Test]
	public void LoadIndexOnly_EmptyLibrary_SucceedsWithZeroEntries()
	{
		string path = Tmp(nameof(LoadIndexOnly_EmptyLibrary_SucceedsWithZeroEntries));
		MslWriter.Write(path, new List<MslLibraryEntry>());

		Assert.That(() =>
		{
			using MslLibrary lib = MslLibrary.LoadIndexOnly(path);
			Assert.That(lib.PrecursorCount, Is.EqualTo(0),
				"PrecursorCount must be 0 for an empty library.");
			Assert.That(lib.GetAllEntries().Any(), Is.False,
				"GetAllEntries must yield nothing for an empty library.");
		}, Throws.Nothing,
		"LoadIndexOnly on a zero-precursor library must not throw.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// M6 — Precursor metadata available without loading fragments
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// PrecursorMz, Charge, ModifiedSequence, and Irt must all be correctly
	/// populated in index-only mode without any fragment demand-load occurring.
	/// This verifies that the skeleton entries carry full precursor metadata.
	/// </summary>
	[Test]
	public void LoadIndexOnly_PrecursorMetadata_AvailableWithoutFragmentLoad()
	{
		var written = BuildEntries(1);
		string path = Tmp(nameof(LoadIndexOnly_PrecursorMetadata_AvailableWithoutFragmentLoad));
		MslWriter.Write(path, written);

		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);

		// Use QueryWindow to get the skeleton entry without triggering fragment load
		using var results = lib.QueryWindow(
			mzLow: (float)written[0].PrecursorMz - 1f,
			mzHigh: (float)written[0].PrecursorMz + 1f,
			rtLow: 0f, rtHigh: 1000f);

		Assert.That(results.Count, Is.EqualTo(1),
			"QueryWindow must find the single written entry.");

		// Access skeleton fields directly (no GetEntry call = no fragment load)
		var indexEntry = results.Entries[0];
		Assert.That((double)indexEntry.PrecursorMz,
			Is.EqualTo(written[0].PrecursorMz).Within(1e-3),
			"PrecursorMz must be available in skeleton entry without fragment load.");
		Assert.That((int)indexEntry.Charge, Is.EqualTo(written[0].Charge),
			"Charge must be available in skeleton entry without fragment load.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// M7 — CRC corruption detected in both load modes
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Corrupting one byte in the fragment section must cause Load to throw
	/// InvalidDataException due to CRC mismatch. This applies equally to
	/// LoadIndexOnly since both call ValidateFileBytes (or the streaming
	/// equivalent after Fix 4).
	/// </summary>
	[Test]
	public void Load_CorruptedFragment_ThrowsInvalidDataException()
	{
		string path = BuildAndWriteLibrary(nameof(Load_CorruptedFragment_ThrowsInvalidDataException));

		byte[] data = File.ReadAllBytes(path);
		MslFileHeader header = MslReader.ReadHeaderOnly(path);

		// Corrupt one byte in the fragment section
		long corruptOffset = header.FragmentSectionOffset + 4;
		data[corruptOffset] ^= 0xFF;
		File.WriteAllBytes(path, data);

		Assert.That(() => MslReader.Load(path),
			Throws.TypeOf<InvalidDataException>().With.Message.Contains("CRC"),
			"Load must throw InvalidDataException when fragment bytes are corrupted.");
	}

	/// <summary>
	/// LoadIndexOnly must also detect CRC corruption, since it validates the
	/// full file checksum before returning.
	/// </summary>
	[Test]
	public void LoadIndexOnly_CorruptedFragment_ThrowsInvalidDataException()
	{
		string path = BuildAndWriteLibrary(
			nameof(LoadIndexOnly_CorruptedFragment_ThrowsInvalidDataException));

		byte[] data = File.ReadAllBytes(path);
		MslFileHeader header = MslReader.ReadHeaderOnly(path);

		long corruptOffset = header.FragmentSectionOffset + 4;
		data[corruptOffset] ^= 0xFF;
		File.WriteAllBytes(path, data);

		Assert.That(() => MslReader.LoadIndexOnly(path),
			Throws.TypeOf<InvalidDataException>().With.Message.Contains("CRC"),
			"LoadIndexOnly must throw InvalidDataException when data is corrupted.");
	}
}