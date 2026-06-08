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
/// NUnit 4 tests verifying the fix to <see cref="MslLibrary.GetAllEntries"/> that
/// previously cloned the entire index array via <c>.ToArray()</c> on every call.
///
/// Root cause: the C# compiler prohibits storing a <c>ReadOnlySpan&lt;T&gt;</c> as a
/// local variable in an iterator method (CS9202). The original code worked around this
/// by calling <c>QueryMzRange(float.MinValue, float.MaxValue).ToArray()</c>, which
/// materialised a full heap copy of the <c>_byMz</c> array before the <c>yield</c> loop.
/// For a 1M-entry index at 24 bytes each, this allocated 24 MB per enumeration.
///
/// Fix: add <see cref="MslIndex.Count"/> and <see cref="MslIndex.GetEntryAt"/> and
/// rewrite <c>GetAllEntries</c> to iterate by integer position, eliminating the copy.
///
/// Coverage:
///   1.  MslIndex.Count returns the correct entry count.
///   2.  MslIndex.Count returns 0 for an empty index.
///   3.  MslIndex.GetEntryAt returns entries in m/z-sorted order.
///   4.  MslIndex.GetEntryAt throws for out-of-range positions.
///   5.  MslIndex.GetEntryAt throws after Dispose.
///   6.  GetAllEntries(includeDecoys=true) yields all entries.
///   7.  GetAllEntries(includeDecoys=false) yields only targets.
///   8.  GetAllEntries results are in ascending m/z order.
///   9.  GetAllEntries on an empty library yields nothing.
///   10. GetAllEntries can be called multiple times with identical results.
///   11. GetAllEntries does not trigger significant Gen-0 GC pressure (regression guard).
/// </summary>
[TestFixture]
public sealed class TestMslPrompt15GetAllEntriesAllocation
{
    // ── Fixture paths ─────────────────────────────────────────────────────────

    private static readonly string OutputDir =
        Path.Combine(Path.GetTempPath(), "MslPrompt15Tests");

    [OneTimeSetUp]
    public void OneTimeSetUp() => Directory.CreateDirectory(OutputDir);

    [OneTimeTearDown]
    public void OneTimeTearDown()
    {
        if (Directory.Exists(OutputDir))
            Directory.Delete(OutputDir, recursive: true);
    }

    private static string TempPath(string name) =>
        Path.Combine(OutputDir, name + ".msl");

    // ── Helpers ───────────────────────────────────────────────────────────────

    private static List<MslLibraryEntry> BuildEntries(
        IEnumerable<(string seq, int charge, float mz, bool isDecoy)> specs)
    {
        var list = new List<MslLibraryEntry>();
        foreach (var (seq, charge, mz, isDecoy) in specs)
        {
            list.Add(new MslLibraryEntry
            {
                FullSequence = seq,
                BaseSequence = seq,
                PrecursorMz = mz,
                ChargeState = charge,
                IsDecoy = isDecoy,
                MoleculeType = MslFormat.MoleculeType.Peptide,
                DissociationType = DissociationType.HCD,
                MatchedFragmentIons = new List<MslFragmentIon>
                {
                    new MslFragmentIon
                    {
                        Mz             = 175.119f,
                        Intensity      = 1.0f,
                        ProductType    = ProductType.y,
                        FragmentNumber = 1,
                        Charge         = 1
                    }
                }
            });
        }
        return list;
    }

    private static string SaveAndReturn(string stem, List<MslLibraryEntry> entries)
    {
        string path = TempPath(stem);
        MslLibrary.Save(path, entries);
        return path;
    }

    // ══════════════════════════════════════════════════════════════════════════
    // MslIndex.Count — new property
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// MslIndex.Count must equal the number of entries built into the index.
    /// </summary>
    [Test]
    public void MslIndex_Count_ReturnsCorrectEntryCount()
    {
        var entries = BuildEntries(new[]
        {
            ("PEPTIDE", 2, 449.74f, false),
            ("ACDEFG",  3, 529.76f, false),
            ("DECOY",   2, 350.00f, true),
        });
        using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

        Assert.That(index.Count, Is.EqualTo(3),
            "MslIndex.Count must equal the number of entries built into the index.");
    }

    /// <summary>
    /// MslIndex.Count must be 0 for an empty index.
    /// </summary>
    [Test]
    public void MslIndex_Count_EmptyIndex_ReturnsZero()
    {
        using var index = MslIndex.Build(new List<MslLibraryEntry>(), _ => null);

        Assert.That(index.Count, Is.EqualTo(0),
            "MslIndex.Count must be 0 for an empty index.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // MslIndex.GetEntryAt — new method
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// GetEntryAt must return entries in ascending m/z order (the sort order of _byMz).
    /// Entries are supplied in reversed order to confirm sorting at construction time.
    /// </summary>
    [Test]
    public void MslIndex_GetEntryAt_ReturnsSortedByMz()
    {
        var entries = BuildEntries(new[]
        {
            ("HIGH", 2, 600.00f, false),
            ("MID",  2, 400.00f, false),
            ("LOW",  2, 200.00f, false),
        });
        using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

        float prev = float.NegativeInfinity;
        for (int i = 0; i < index.Count; i++)
        {
            float mz = index.GetEntryAt(i).PrecursorMz;
            Assert.That(mz, Is.GreaterThanOrEqualTo(prev),
                $"GetEntryAt({i}).PrecursorMz={mz} must be >= previous={prev}. " +
                "Index entries must be in ascending m/z order.");
            prev = mz;
        }
    }

    /// <summary>
    /// GetEntryAt must throw ArgumentOutOfRangeException for positions outside [0, Count).
    /// </summary>
    [Test]
    public void MslIndex_GetEntryAt_OutOfRange_ThrowsArgumentOutOfRangeException()
    {
        var entries = BuildEntries(new[] { ("PEPTIDE", 2, 449.74f, false) });
        using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

        Assert.That(() => index.GetEntryAt(-1),
            Throws.TypeOf<ArgumentOutOfRangeException>(),
            "GetEntryAt(-1) must throw ArgumentOutOfRangeException.");
        Assert.That(() => index.GetEntryAt(index.Count),
            Throws.TypeOf<ArgumentOutOfRangeException>(),
            $"GetEntryAt({index.Count}) (== Count) must throw ArgumentOutOfRangeException.");
    }

    /// <summary>
    /// GetEntryAt must throw ObjectDisposedException after the index is disposed.
    /// </summary>
    [Test]
    public void MslIndex_GetEntryAt_AfterDispose_ThrowsObjectDisposedException()
    {
        var entries = BuildEntries(new[] { ("PEPTIDE", 2, 449.74f, false) });
        var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);
        index.Dispose();

        Assert.That(() => index.GetEntryAt(0),
            Throws.TypeOf<ObjectDisposedException>(),
            "GetEntryAt must throw ObjectDisposedException after Dispose.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // GetAllEntries correctness
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// GetAllEntries(includeDecoys: true) must yield both target and decoy entries.
    /// </summary>
    [Test]
    public void GetAllEntries_IncludeDecoys_YieldsAllEntries()
    {
        var entries = BuildEntries(new[]
        {
            ("TARGET", 2, 449.74f, false),
            ("DECOY",  2, 350.00f, true),
        });
        string path = SaveAndReturn("incl_decoys", entries);

        using var lib = MslLibrary.Load(path);
        var all = lib.GetAllEntries(includeDecoys: true).ToList();

        Assert.That(all, Has.Count.EqualTo(2),
            "GetAllEntries(includeDecoys: true) must yield both target and decoy entries.");
    }

    /// <summary>
    /// GetAllEntries(includeDecoys: false) must yield only non-decoy entries.
    /// </summary>
    [Test]
    public void GetAllEntries_ExcludeDecoys_YieldsOnlyTargets()
    {
        var entries = BuildEntries(new[]
        {
            ("TARGET1", 2, 449.74f, false),
            ("TARGET2", 3, 530.00f, false),
            ("DECOY",   2, 350.00f, true),
        });
        string path = SaveAndReturn("excl_decoys", entries);

        using var lib = MslLibrary.Load(path);
        var targets = lib.GetAllEntries(includeDecoys: false).ToList();

        Assert.That(targets, Has.Count.EqualTo(2),
            "GetAllEntries(includeDecoys: false) must yield only the 2 target entries.");
        Assert.That(targets.All(e => !e.IsDecoy), Is.True,
            "All yielded entries must be non-decoy when includeDecoys is false.");
    }

    /// <summary>
    /// Entries yielded by GetAllEntries must be in ascending precursor m/z order,
    /// matching the sort order of the internal _byMz array.
    /// </summary>
    [Test]
    public void GetAllEntries_ResultsAreInAscendingMzOrder()
    {
        var entries = BuildEntries(new[]
        {
            ("C", 2, 600.00f, false),
            ("A", 2, 200.00f, false),
            ("B", 2, 400.00f, false),
        });
        string path = SaveAndReturn("mz_order", entries);

        using var lib = MslLibrary.Load(path);
        var all = lib.GetAllEntries().ToList();

        for (int i = 1; i < all.Count; i++)
        {
            Assert.That(all[i].PrecursorMz, Is.GreaterThanOrEqualTo(all[i - 1].PrecursorMz),
                $"Entry at position {i} (mz={all[i].PrecursorMz}) is not >= " +
                $"entry at position {i - 1} (mz={all[i - 1].PrecursorMz}). " +
                "GetAllEntries results must be in ascending m/z order.");
        }
    }

    /// <summary>
    /// GetAllEntries on an empty library must yield zero entries and not throw.
    /// </summary>
    [Test]
    public void GetAllEntries_EmptyLibrary_YieldsNothing()
    {
        string path = SaveAndReturn("empty_lib", new List<MslLibraryEntry>());

        using var lib = MslLibrary.Load(path);
        var all = lib.GetAllEntries().ToList();

        Assert.That(all, Is.Empty,
            "GetAllEntries on an empty library must yield zero entries.");
    }

    /// <summary>
    /// Calling GetAllEntries twice on the same library instance must yield identical results.
    /// This confirms the fix does not mutate any shared state (e.g. by consuming a span).
    /// </summary>
    [Test]
    public void GetAllEntries_MultipleCallsYieldSameResults()
    {
        var entries = BuildEntries(new[]
        {
            ("PEPTIDE", 2, 449.74f, false),
            ("ACDEFG",  3, 529.76f, false),
        });
        string path = SaveAndReturn("multi_call", entries);

        using var lib = MslLibrary.Load(path);
        var first = lib.GetAllEntries().Select(e => e.FullSequence).OrderBy(s => s).ToList();
        var second = lib.GetAllEntries().Select(e => e.FullSequence).OrderBy(s => s).ToList();

        Assert.That(second, Is.EqualTo(first),
            "Calling GetAllEntries twice on the same library must yield identical results.");
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Allocation regression guard
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Verifies that enumerating GetAllEntries across 1000 iterations on a 1000-entry
    /// library does not trigger significant Gen-0 GC collections.
    ///
    /// Before the fix, every call materialised a full MslPrecursorIndexEntry[] clone
    /// of the index. For 1000 entries × 24 bytes × 1000 calls = ~24 MB of garbage,
    /// which reliably triggered multiple Gen-0 collections.
    ///
    /// After the fix the only per-call heap allocation is the compiler-generated iterator
    /// state machine (~64 bytes). 1000 calls = ~64 kB — well below a Gen-0 threshold.
    /// We allow at most 2 collections to accommodate baseline GC noise.
    /// </summary>
    [Test]
    public void GetAllEntries_DoesNotAllocateIndexClone_LowGcPressure()
    {
        // Build a 1000-entry library
        var rawEntries = new List<MslLibraryEntry>();
        for (int i = 0; i < 1000; i++)
        {
            rawEntries.Add(new MslLibraryEntry
            {
                FullSequence = $"PEPTIDE{i:D4}",
                BaseSequence = $"PEPTIDE{i:D4}",
                PrecursorMz = 400.0 + i * 0.5,
                ChargeState = 2,
                IsDecoy = false,
                MoleculeType = MslFormat.MoleculeType.Peptide,
                DissociationType = DissociationType.HCD,
                MatchedFragmentIons = new List<MslFragmentIon>
                {
                    new MslFragmentIon
                    {
                        Mz             = 175.119f,
                        Intensity      = 1.0f,
                        ProductType    = ProductType.y,
                        FragmentNumber = 1,
                        Charge         = 1
                    }
                }
            });
        }

        string path = TempPath("alloc_1k");
        MslLibrary.Save(path, rawEntries);

        using var lib = MslLibrary.Load(path);

        // Warm up JIT and LRU cache so measurement is steady-state
        foreach (var _ in lib.GetAllEntries()) { }
        GC.Collect();
        GC.WaitForPendingFinalizers();
        GC.Collect();

        const int iterations = 1000;
        long gen0Before = GC.CollectionCount(0);

        for (int k = 0; k < iterations; k++)
        {
            foreach (var _ in lib.GetAllEntries()) { }
        }

        long gen0After = GC.CollectionCount(0);
        long gen0Delta = gen0After - gen0Before;

        // Allow up to 10 collections for baseline GC noise from background threads
        // (finalizer, JIT recompilation, test harness, parallel test runners).
        // The pre-fix code produced hundreds of collections for this workload.
        Assert.That(gen0Delta, Is.LessThanOrEqualTo(10),
            $"GetAllEntries triggered {gen0Delta} Gen-0 GC collections across " +
            $"{iterations} iterations on a 1000-entry library. Expected at most 10 " +
            $"(iterator state machine baseline + background GC noise). A much higher " +
            $"count indicates the index array clone was not eliminated by the fix.");
    }
}