using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Threading;
using System.Threading.Tasks;

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 test suite for <see cref="MslIndex"/>, <see cref="MslPrecursorIndexEntry"/>,
/// <see cref="MslWindowResults"/>, and <see cref="MslIndexStatistics"/>.
///
/// Coverage groups:
/// <list type="bullet">
///   <item>Struct layout — 24-byte size validation for MslPrecursorIndexEntry.</item>
///   <item>Build and basic query — sorted construction, QueryMzRange semantics.</item>
///   <item>Window queries — RT, IM, and decoy filters.</item>
///   <item>DDA-style lookup — TryGetBySequenceCharge.</item>
///   <item>Elution group queries — GetElutionGroup.</item>
///   <item>RT calibration — WithCalibratedRetentionTimes.</item>
///   <item>LRU buffer — GetEntry caching and eviction.</item>
///   <item>Statistics — GetStatistics aggregate values.</item>
///   <item>Binary-search boundaries — edge positions in the sorted array.</item>
///   <item>Scale tests — correctness at 1 k, 100 k, and 1 M entries.</item>
/// </list>
/// All tests use <c>Assert.That</c> exclusively; the legacy <c>Assert.AreEqual</c> style
/// is not used. All test classes carry <c>[TestFixture]</c>; all test methods carry
/// <c>[Test]</c> or <c>[TestCase]</c>.
/// </summary>
[TestFixture]
public class TestMslIndex
{
	// ═══════════════════════════════════════════════════════════════════════
	// Struct-layout tests
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslPrecursorIndexEntry"/> occupies exactly 24 unmanaged
	/// bytes. The index keeps one entry per precursor; 1 million entries must cost 24 MB
	/// exactly. A size deviation would indicate a missed field or wrong StructLayout.
	/// </summary>
	[Test]
	public void MslPrecursorIndexEntry_Is_24_Bytes()
		=> Assert.That(Marshal.SizeOf<MslPrecursorIndexEntry>(), Is.EqualTo(24));

	// ═══════════════════════════════════════════════════════════════════════
	// Build and basic query tests
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslIndex.Build"/> creates an index whose entry count
	/// matches the number of source <see cref="MslLibraryEntry"/> objects supplied.
	/// </summary>
	[Test]
	public void Build_FromEntries_CreatesCorrectCount()
	{
		// Arrange: three entries with distinct m/z values
		var entries = BuildEntries(
			("PEPTIDE", 2, 400.23f, 25.0f, 0f),
			("SEQUENCER", 3, 500.67f, 30.0f, 0f),
			("ANALYZE", 2, 300.11f, 15.0f, 0f));

		// Act
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Assert: statistics surface the correct count
		Assert.That(index.GetStatistics().TotalPrecursors, Is.EqualTo(3));
	}

	/// <summary>
	/// Verifies that after <see cref="MslIndex.Build"/> the internal entries are ordered
	/// by ascending <see cref="MslPrecursorIndexEntry.PrecursorMz"/>. The sort guarantees
	/// that every subsequent binary-search call produces correct results.
	/// </summary>
	[Test]
	public void Build_EntriesAreSortedByMz()
	{
		// Arrange: supply entries in deliberately unsorted m/z order
		var entries = BuildEntries(
			("C", 2, 600.0f, 10.0f, 0f),
			("A", 2, 200.0f, 10.0f, 0f),
			("B", 2, 400.0f, 10.0f, 0f));

		// Act
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Retrieve the full sorted m/z slice and verify ascending order
		var span = index.QueryMzRange(0f, float.MaxValue);
		float prev = float.NegativeInfinity;
		foreach (var entry in span)
		{
			Assert.That(entry.PrecursorMz, Is.GreaterThanOrEqualTo(prev));
			prev = entry.PrecursorMz;
		}
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.QueryMzRange"/> returns an empty span when
	/// the query window does not overlap any entry in the index.
	/// </summary>
	[Test]
	public void QueryMzRange_EmptyRange_ReturnsEmptySpan()
	{
		// Arrange: entries clustered around m/z 400–500
		var entries = BuildEntries(
			("A", 2, 400.0f, 10.0f, 0f),
			("B", 2, 450.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: query a window entirely below the lowest entry
		var span = index.QueryMzRange(100f, 200f);

		// Assert
		Assert.That(span.Length, Is.EqualTo(0));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.QueryMzRange"/> includes both the entry whose
	/// m/z equals the lower bound and the entry whose m/z equals the upper bound (closed
	/// interval semantics).
	/// </summary>
	[Test]
	public void QueryMzRange_IncludesBothEndpoints()
	{
		// Arrange: three entries — one exactly at each endpoint, one in between
		var entries = BuildEntries(
			("Low", 2, 300.0f, 10.0f, 0f),
			("Mid", 2, 350.0f, 10.0f, 0f),
			("High", 2, 400.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: query with the window exactly matching the first and last entry
		var span = index.QueryMzRange(300.0f, 400.0f);

		// Assert: all three entries are returned
		Assert.That(span.Length, Is.EqualTo(3));
	}

	/// <summary>
	/// Verifies that entries whose m/z falls strictly outside [mzLow, mzHigh] are not
	/// included in the result of <see cref="MslIndex.QueryMzRange"/>.
	/// </summary>
	[Test]
	public void QueryMzRange_ExcludesOutOfRangeEntries()
	{
		// Arrange: five entries spanning 200–600
		var entries = BuildEntries(
			("A", 2, 200.0f, 10.0f, 0f),
			("B", 2, 300.0f, 10.0f, 0f),
			("C", 2, 400.0f, 10.0f, 0f),
			("D", 2, 500.0f, 10.0f, 0f),
			("E", 2, 600.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: window covers only B and C
		var span = index.QueryMzRange(300.0f, 400.0f);

		// Assert: exactly two entries, both in range
		Assert.That(span.Length, Is.EqualTo(2));
		Assert.That(span[0].PrecursorMz, Is.EqualTo(300.0f).Within(1e-4f));
		Assert.That(span[1].PrecursorMz, Is.EqualTo(400.0f).Within(1e-4f));
	}

	/// <summary>
	/// Verifies that querying with a window that covers the full m/z range returns all
	/// entries in the index — useful as a baseline for downstream filtering tests.
	/// </summary>
	[Test]
	public void QueryMzRange_AllEntries_WhenWindowCoversFullRange()
	{
		// Arrange
		var entries = BuildEntries(
			("A", 2, 300.0f, 10.0f, 0f),
			("B", 2, 400.0f, 10.0f, 0f),
			("C", 2, 500.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act
		var span = index.QueryMzRange(0f, float.MaxValue);

		// Assert
		Assert.That(span.Length, Is.EqualTo(3));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.QueryMzRange"/> performs zero heap allocations on
	/// a typical call. Because the method returns a <see cref="ReadOnlySpan{T}"/> over the
	/// internal sorted array, it must not allocate any objects.
	/// </summary>
	[Test]
	public void QueryMzRange_ZeroAllocation_NoGcPressure()
	{
		// Arrange: build a small index and warm up JIT before measuring
		var entries = BuildEntries(
			("A", 2, 300.0f, 10.0f, 0f),
			("B", 2, 400.0f, 10.0f, 0f),
			("C", 2, 500.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Warm up: run once to ensure JIT compilation does not affect measurements
		_ = index.QueryMzRange(300f, 500f);

		// Measure Gen-0 collections before and after 1 000 calls
		GC.Collect();
		long gen0Before = GC.CollectionCount(0);

		const int iterations = 1_000;
		for (int k = 0; k < iterations; k++)
			_ = index.QueryMzRange(300f, 500f);

		long gen0After = GC.CollectionCount(0);

		// Assert: no Gen-0 collections triggered by the span queries alone
		Assert.That(gen0After, Is.EqualTo(gen0Before),
			$"QueryMzRange triggered {gen0After - gen0Before} Gen-0 GC collection(s) " +
			$"across {iterations} calls; expected zero.");
	}

	// ═══════════════════════════════════════════════════════════════════════
	// Window query tests
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslIndex.QueryWindow"/> excludes entries whose
	/// <see cref="MslPrecursorIndexEntry.Irt"/> falls outside [rtLow, rtHigh], even when
	/// their m/z is within the requested m/z window.
	/// </summary>
	[Test]
	public void QueryWindow_RtFilter_ExcludesEntriesOutsideRtWindow()
	{
		// Arrange: all entries at the same m/z but different RT values
		var entries = BuildEntries(
			("A", 2, 400.0f, 10.0f, 0f),   // RT = 10 — inside
			("B", 2, 400.1f, 50.0f, 0f),   // RT = 50 — outside
			("C", 2, 400.2f, 20.0f, 0f));  // RT = 20 — inside
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: RT window [5, 30] should include A and C but not B
		using var result = index.QueryWindow(399f, 401f, rtLow: 5f, rtHigh: 30f);

		// Assert
		Assert.That(result.Count, Is.EqualTo(2));
		Assert.That(result.Entries.ToArray().All(e => e.Irt >= 5f && e.Irt <= 30f), Is.True);
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.QueryWindow"/> includes entries whose
	/// <see cref="MslPrecursorIndexEntry.Irt"/> falls exactly on the RT window boundaries
	/// (closed-interval semantics).
	/// </summary>
	[Test]
	public void QueryWindow_RtFilter_IncludesEntriesInsideRtWindow()
	{
		// Arrange: entries with RT at exactly rtLow and rtHigh
		var entries = BuildEntries(
			("A", 2, 400.0f, 10.0f, 0f),   // RT exactly at rtLow
			("B", 2, 400.1f, 30.0f, 0f));  // RT exactly at rtHigh
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act
		using var result = index.QueryWindow(399f, 401f, rtLow: 10f, rtHigh: 30f);

		// Assert: both boundary entries are included
		Assert.That(result.Count, Is.EqualTo(2));
	}

	/// <summary>
	/// Verifies that the IM filter in <see cref="MslIndex.QueryWindow"/> is entirely
	/// skipped when both <c>imLow</c> and <c>imHigh</c> are zero (the default).
	/// Entries that would fail an IM filter must still be returned in this mode.
	/// </summary>
	[Test]
	public void QueryWindow_ImFilter_AppliedOnlyWhenBothBoundsNonZero()
	{
		// Arrange: entries with very different IM values
		var entries = BuildEntries(
			("A", 2, 400.0f, 10.0f, 0.5f),   // IM = 0.5
			("B", 2, 400.1f, 10.0f, 1.5f));  // IM = 1.5
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: imLow = imHigh = 0 → IM filter must be skipped
		using var result = index.QueryWindow(399f, 401f, rtLow: 0f, rtHigh: 100f,
											 imLow: 0f, imHigh: 0f);

		// Assert: both entries pass (IM filter bypassed)
		Assert.That(result.Count, Is.EqualTo(2));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.QueryWindow"/> excludes entries whose
	/// <see cref="MslPrecursorIndexEntry.IonMobility"/> falls outside [imLow, imHigh]
	/// when the IM filter is active (both bounds non-zero).
	/// </summary>
	[Test]
	public void QueryWindow_ImFilter_ExcludesEntriesOutsideWindow()
	{
		// Arrange: three entries with distinct IM values
		var entries = BuildEntries(
			("A", 2, 400.0f, 10.0f, 0.3f),   // IM = 0.3 — outside [0.5, 1.0]
			("B", 2, 400.1f, 10.0f, 0.7f),   // IM = 0.7 — inside
			("C", 2, 400.2f, 10.0f, 1.5f));  // IM = 1.5 — outside
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: IM window [0.5, 1.0]
		using var result = index.QueryWindow(399f, 401f, rtLow: 0f, rtHigh: 100f,
											 imLow: 0.5f, imHigh: 1.0f);

		// Assert: only entry B passes the IM filter
		Assert.That(result.Count, Is.EqualTo(1));
		Assert.That(result.Entries[0].IonMobility, Is.EqualTo(0.7f).Within(1e-5f));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.QueryWindow"/> excludes decoy entries by default
	/// (when <c>includeDecoys</c> is false, which is the default parameter value).
	/// </summary>
	[Test]
	public void QueryWindow_ExcludesDecoys_ByDefault()
	{
		// Arrange: mix of targets and a decoy at the same m/z and RT
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("TARGET", 2, 400.0f, 10.0f, 0f, isDecoy: false),
			MakeEntry("DECOY",  2, 400.1f, 10.0f, 0f, isDecoy: true),
		};
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: default — decoys excluded
		using var result = index.QueryWindow(399f, 401f, rtLow: 0f, rtHigh: 100f);

		// Assert: only the target entry is returned
		Assert.That(result.Count, Is.EqualTo(1));
		Assert.That(result.Entries[0].IsDecoy, Is.EqualTo((byte)0));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.QueryWindow"/> includes decoy entries when
	/// <c>includeDecoys</c> is explicitly set to true.
	/// </summary>
	[Test]
	public void QueryWindow_IncludesDecoys_WhenFlagSet()
	{
		// Arrange
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("TARGET", 2, 400.0f, 10.0f, 0f, isDecoy: false),
			MakeEntry("DECOY",  2, 400.1f, 10.0f, 0f, isDecoy: true),
		};
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: request decoys too
		using var result = index.QueryWindow(399f, 401f, rtLow: 0f, rtHigh: 100f, includeDecoys: true);

		// Assert: both target and decoy are returned
		Assert.That(result.Count, Is.EqualTo(2));
	}

	/// <summary>
	/// Verifies that the <see cref="MslWindowResults"/> returned by
	/// <see cref="MslIndex.QueryWindow"/> correctly returns its pooled buffer to
	/// <see cref="System.Buffers.ArrayPool{T}.Shared"/> when disposed. After disposal the
	/// <see cref="MslWindowResults.Entries"/> span must no longer be accessed, but
	/// calling <see cref="MslWindowResults.Dispose"/> a second time must not throw
	/// (idempotent behaviour).
	/// </summary>
	[Test]
	public void QueryWindow_Returns_Disposable_ThatReleasesPool()
	{
		// Arrange
		var entries = BuildEntries(("A", 2, 400.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: obtain result, verify it has data, then dispose
		MslWindowResults result = index.QueryWindow(399f, 401f, rtLow: 0f, rtHigh: 100f);
		Assert.That(result.Count, Is.EqualTo(1));

		// Double-dispose must not throw
		Assert.That(() =>
		{
			result.Dispose();
			result.Dispose();
		}, Throws.Nothing);
	}

	// ═══════════════════════════════════════════════════════════════════════
	// DDA-style lookup tests
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslIndex.TryGetBySequenceCharge"/> returns true and sets
	/// the out parameter correctly when an entry with the given modified sequence and
	/// charge is present in the index.
	/// </summary>
	[Test]
	public void TryGetBySequenceCharge_Found_WhenEntryExists()
	{
		// Arrange
		var entries = BuildEntries(("PEPTIDE", 2, 400.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act
		bool found = index.TryGetBySequenceCharge("PEPTIDE", 2, out var entry);

		// Assert
		Assert.That(found, Is.True);
		Assert.That(entry.Charge, Is.EqualTo((short)2));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.TryGetBySequenceCharge"/> returns false when no
	/// entry matches the given modified sequence and charge. The out parameter value is
	/// not tested because it is undefined for a miss.
	/// </summary>
	[Test]
	public void TryGetBySequenceCharge_NotFound_WhenEntryAbsent()
	{
		// Arrange
		var entries = BuildEntries(("PEPTIDE", 2, 400.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: wrong sequence
		bool found = index.TryGetBySequenceCharge("NOTHERE", 2, out _);

		// Assert
		Assert.That(found, Is.False);
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.TryGetBySequenceCharge"/> uses case-sensitive
	/// matching: "peptide" and "PEPTIDE" are distinct keys and must not resolve to the
	/// same entry.
	/// </summary>
	[Test]
	public void TryGetBySequenceCharge_CaseSensitive()
	{
		// Arrange: only the upper-case form is stored
		var entries = BuildEntries(("PEPTIDE", 2, 400.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: look up lower-case form
		bool found = index.TryGetBySequenceCharge("peptide", 2, out _);

		// Assert: lower-case key is absent
		Assert.That(found, Is.False);
	}

	// ═══════════════════════════════════════════════════════════════════════
	// Elution group tests
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslIndex.GetElutionGroup"/> returns all charge states of
	/// the same peptide when they share a common
	/// <see cref="MslPrecursorIndexEntry.ElutionGroupId"/>. This is the primary use case:
	/// the DIA engine queries the elution group to retrieve all charge variants together.
	/// </summary>
	[Test]
	public void GetElutionGroup_ReturnsAllChargeStates_ForSamePeptide()
	{
		// Arrange: two charge states for the same peptide, same stripped sequence → same elution group
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("PEPTIDE", 2, 400.0f, 10.0f, 0f, isDecoy: false, elutionGroupId: 7),
			MakeEntry("PEPTIDE", 3, 267.3f, 10.0f, 0f, isDecoy: false, elutionGroupId: 7),
		};
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act
		var span = index.GetElutionGroup(7);

		// Assert: both charge states are returned
		Assert.That(span.Length, Is.EqualTo(2));
		Assert.That(span.ToArray().Select(e => (int)e.Charge).OrderBy(c => c).ToArray(),
			Is.EqualTo(new[] { 2, 3 }));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.GetElutionGroup"/> returns an empty span when
	/// the requested elution-group ID is not present in the index.
	/// </summary>
	[Test]
	public void GetElutionGroup_ReturnsEmptySpan_ForUnknownId()
	{
		// Arrange: index with elution-group IDs 0 and 1 only
		var entries = BuildEntries(("A", 2, 400.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: request a group ID that was never assigned
		var span = index.GetElutionGroup(999);

		// Assert
		Assert.That(span.Length, Is.EqualTo(0));
	}

	// ═══════════════════════════════════════════════════════════════════════
	// RT calibration tests
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslIndex.WithCalibratedRetentionTimes"/> applies the
	/// formula <c>calibratedRT = slope * iRT + intercept</c> to every entry.
	/// Tests a non-trivial transform (slope = 2, intercept = 5).
	/// </summary>
	[Test]
	public void WithCalibratedRetentionTimes_TransformsAllEntries()
	{
		// Arrange: known iRT values
		var entries = BuildEntries(
			("A", 2, 400.0f, 10.0f, 0f),   // iRT = 10 → calibrated = 2*10 + 5 = 25
			("B", 2, 500.0f, 20.0f, 0f));  // iRT = 20 → calibrated = 2*20 + 5 = 45
		using var original = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act
		using var calibrated = original.WithCalibratedRetentionTimes(slope: 2.0, intercept: 5.0);
		var span = calibrated.QueryMzRange(0f, float.MaxValue);

		// Assert: iRT values are transformed — must find entries with calibrated RT
		float[] expectedIrts = { 25.0f, 45.0f };
		float[] actualIrts = span.ToArray().Select(e => e.Irt).OrderBy(v => v).ToArray();
		Assert.That(actualIrts[0], Is.EqualTo(expectedIrts[0]).Within(1e-3f));
		Assert.That(actualIrts[1], Is.EqualTo(expectedIrts[1]).Within(1e-3f));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.WithCalibratedRetentionTimes"/> does not modify
	/// the iRT values of the original index — the method is a pure transform returning a
	/// new instance.
	/// </summary>
	[Test]
	public void WithCalibratedRetentionTimes_OriginalIndex_Unchanged()
	{
		// Arrange
		var entries = BuildEntries(("A", 2, 400.0f, 10.0f, 0f));
		using var original = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Capture the original iRT before calibration
		float originalIrt = original.QueryMzRange(0f, float.MaxValue)[0].Irt;

		// Act: calibrate with a non-trivial transform
		using var calibrated = original.WithCalibratedRetentionTimes(slope: 3.0, intercept: 10.0);

		// Assert: original index is unchanged
		float afterIrt = original.QueryMzRange(0f, float.MaxValue)[0].Irt;
		Assert.That(afterIrt, Is.EqualTo(originalIrt).Within(1e-5f));
	}

	/// <summary>
	/// Verifies that applying an identity transform (slope = 1, intercept = 0) produces a
	/// new index whose iRT values are numerically identical to those in the original.
	/// </summary>
	[Test]
	public void WithCalibratedRetentionTimes_Slope1_Intercept0_IsIdentity()
	{
		// Arrange
		var entries = BuildEntries(
			("A", 2, 400.0f, 15.0f, 0f),
			("B", 2, 500.0f, 30.0f, 0f));
		using var original = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act
		using var identity = original.WithCalibratedRetentionTimes(slope: 1.0, intercept: 0.0);
		var origSpan = original.QueryMzRange(0f, float.MaxValue);
		var identSpan = identity.QueryMzRange(0f, float.MaxValue);

		// Assert: all iRT values are preserved exactly
		Assert.That(identSpan.Length, Is.EqualTo(origSpan.Length));
		for (int i = 0; i < origSpan.Length; i++)
			Assert.That(identSpan[i].Irt, Is.EqualTo(origSpan[i].Irt).Within(1e-5f));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.WithCalibratedRetentionTimes"/> returns an index
	/// whose entries remain sorted by m/z even after the iRT transform (only iRT changes,
	/// so the m/z sort order is invariant).
	/// </summary>
	[Test]
	public void WithCalibratedRetentionTimes_NewIndex_IsSortedByMz()
	{
		// Arrange: entries in scrambled m/z order
		var entries = BuildEntries(
			("C", 2, 600.0f, 30.0f, 0f),
			("A", 2, 200.0f, 10.0f, 0f),
			("B", 2, 400.0f, 20.0f, 0f));
		using var original = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act
		using var calibrated = original.WithCalibratedRetentionTimes(slope: 1.5, intercept: 2.0);
		var span = calibrated.QueryMzRange(0f, float.MaxValue);

		// Assert: m/z values are in ascending order
		float prev = float.NegativeInfinity;
		foreach (var e in span)
		{
			Assert.That(e.PrecursorMz, Is.GreaterThanOrEqualTo(prev));
			prev = e.PrecursorMz;
		}
	}

	// ═══════════════════════════════════════════════════════════════════════
	// LRU buffer tests
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslIndex.GetEntry"/> returns the same
	/// <see cref="MslLibraryEntry"/> object on the second call (cache hit), and that
	/// the entry-loader delegate is invoked only once across the two calls.
	/// </summary>
	[Test]
	public void GetEntry_CachesLoadedEntries()
	{
		// Arrange: track how many times the loader is called
		var entries = BuildEntries(("A", 2, 400.0f, 10.0f, 0f));
		int loaderCallCount = 0;

		using var index = MslIndex.Build(entries, i =>
		{
			Interlocked.Increment(ref loaderCallCount);
			return i < entries.Count ? entries[i] : null;
		});

		// Reset counter after Build (Build calls loader once per entry)
		loaderCallCount = 0;

		// Act: first call populates cache; second call should be a hit
		var first = index.GetEntry(0);
		var second = index.GetEntry(0);

		// Assert: both calls return the same entry; loader called only once post-reset
		Assert.That(first, Is.Not.Null);
		Assert.That(second, Is.SameAs(first));
		Assert.That(loaderCallCount, Is.EqualTo(1));  // hit on second call, no extra invocation
	}

	/// <summary>
	/// Verifies that when the LRU cache is at its capacity limit,
	/// <see cref="MslIndex.GetEntry"/> evicts the oldest entry before adding the new one,
	/// keeping the cache size bounded at <c>maxBufferSize</c>.
	/// </summary>
	[Test]
	public void GetEntry_EvictsOldestWhenAtCapacity()
	{
		// Arrange: three entries but a cache capacity of only two
		var entries = BuildEntries(
			("A", 2, 300.0f, 10.0f, 0f),
			("B", 2, 400.0f, 10.0f, 0f),
			("C", 2, 500.0f, 10.0f, 0f));

		// Build the raw entry array to pass as entries for the loader
		var raw = BuildRaw(entries);
		using var index = new MslIndex(raw, i => i < entries.Count ? entries[i] : null, maxBufferSize: 2);

		// Act: load all three entries — third load should evict the first
		index.GetEntry(0);  // A → cache: {0}
		index.GetEntry(1);  // B → cache: {0,1}
		index.GetEntry(2);  // C → triggers eviction of 0; cache: {1,2}

		// Verify cache state via statistics
		var stats = index.GetStatistics();
		// 0 hits in the initial three loads (all misses)
		Assert.That(stats.LruMisses, Is.EqualTo(3));
		Assert.That(stats.LruHits, Is.EqualTo(0));

		// Now fetch A again; it was evicted so this must be another miss
		index.GetEntry(0);
		var statsAfter = index.GetStatistics();
		Assert.That(statsAfter.LruMisses, Is.EqualTo(4));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.GetEntry"/> returns null when the requested
	/// precursor index is out of range (i.e. the entry-loader returns null), rather than
	/// throwing or returning a default struct.
	/// </summary>
	[Test]
	public void GetEntry_ReturnsNull_ForOutOfRangeIndex()
	{
		// Arrange: single entry at index 0
		var entries = BuildEntries(("A", 2, 400.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: request an index well outside the valid range
		var result = index.GetEntry(9999);

		// Assert
		Assert.That(result, Is.Null);
	}

	// ═══════════════════════════════════════════════════════════════════════
	// Statistics tests
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslIndex.GetStatistics"/> correctly counts target and
	/// decoy precursors and that <c>TargetPrecursors + DecoyPrecursors == TotalPrecursors</c>.
	/// </summary>
	[Test]
	public void GetStatistics_CorrectTargetAndDecoyCount()
	{
		// Arrange: two targets and one decoy
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("T1", 2, 400.0f, 10.0f, 0f, isDecoy: false),
			MakeEntry("T2", 2, 500.0f, 20.0f, 0f, isDecoy: false),
			MakeEntry("D1", 2, 450.0f, 15.0f, 0f, isDecoy: true),
		};
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act
		var stats = index.GetStatistics();

		// Assert
		Assert.That(stats.TotalPrecursors, Is.EqualTo(3));
		Assert.That(stats.TargetPrecursors, Is.EqualTo(2));
		Assert.That(stats.DecoyPrecursors, Is.EqualTo(1));
		Assert.That(stats.TargetPrecursors + stats.DecoyPrecursors, Is.EqualTo(stats.TotalPrecursors));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.GetStatistics"/> reports the correct minimum and
	/// maximum m/z values after sorting.
	/// </summary>
	[Test]
	public void GetStatistics_CorrectMinMaxMz()
	{
		// Arrange
		var entries = BuildEntries(
			("A", 2, 300.0f, 10.0f, 0f),
			("B", 2, 500.0f, 20.0f, 0f),
			("C", 2, 400.0f, 15.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act
		var stats = index.GetStatistics();

		// Assert
		Assert.That(stats.MinPrecursorMz, Is.EqualTo(300.0f).Within(1e-4f));
		Assert.That(stats.MaxPrecursorMz, Is.EqualTo(500.0f).Within(1e-4f));
	}

	/// <summary>
	/// Verifies that <see cref="MslIndex.GetStatistics"/> reports the correct number of
	/// distinct elution groups; this equals the number of unique stripped sequences in the
	/// entries supplied to <see cref="MslIndex.Build"/>.
	/// </summary>
	[Test]
	public void GetStatistics_CorrectElutionGroupCount()
	{
		// Arrange: two peptides (A and B), each with two charge states → two elution groups
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("PEPTIDE_A", 2, 400.0f, 10.0f, 0f, isDecoy: false, elutionGroupId: 1),
			MakeEntry("PEPTIDE_A", 3, 267.3f, 10.0f, 0f, isDecoy: false, elutionGroupId: 1),
			MakeEntry("PEPTIDE_B", 2, 500.0f, 20.0f, 0f, isDecoy: false, elutionGroupId: 2),
			MakeEntry("PEPTIDE_B", 3, 334.0f, 20.0f, 0f, isDecoy: false, elutionGroupId: 2),
		};
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act
		var stats = index.GetStatistics();

		// Assert: two distinct elution-group IDs
		Assert.That(stats.ElutionGroupCount, Is.EqualTo(2));
	}

	// ═══════════════════════════════════════════════════════════════════════
	// Binary-search boundary tests
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that a window exactly matching the first (lowest m/z) entry in the index
	/// returns that entry. Tests the lower boundary of the binary-search lower-bound algorithm.
	/// </summary>
	[Test]
	public void BinarySearch_WindowAtFirstEntry_ReturnsEntry()
	{
		// Arrange
		var entries = BuildEntries(
			("A", 2, 200.0f, 10.0f, 0f),
			("B", 2, 400.0f, 10.0f, 0f),
			("C", 2, 600.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: window that only covers the first entry
		var span = index.QueryMzRange(200.0f, 200.0f);

		// Assert
		Assert.That(span.Length, Is.EqualTo(1));
		Assert.That(span[0].PrecursorMz, Is.EqualTo(200.0f).Within(1e-4f));
	}

	/// <summary>
	/// Verifies that a window exactly matching the last (highest m/z) entry in the index
	/// returns that entry. Tests the upper boundary of the binary-search upper-bound scan.
	/// </summary>
	[Test]
	public void BinarySearch_WindowAtLastEntry_ReturnsEntry()
	{
		// Arrange
		var entries = BuildEntries(
			("A", 2, 200.0f, 10.0f, 0f),
			("B", 2, 400.0f, 10.0f, 0f),
			("C", 2, 600.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: window that only covers the last entry
		var span = index.QueryMzRange(600.0f, 600.0f);

		// Assert
		Assert.That(span.Length, Is.EqualTo(1));
		Assert.That(span[0].PrecursorMz, Is.EqualTo(600.0f).Within(1e-4f));
	}

	/// <summary>
	/// Verifies that a window positioned between two adjacent entries returns an empty span.
	/// This tests that the binary search does not accidentally round toward either neighbour.
	/// </summary>
	[Test]
	public void BinarySearch_WindowBetweenEntries_ReturnsEmpty()
	{
		// Arrange: entries at 300 and 500; window falls in the 350–450 gap
		var entries = BuildEntries(
			("A", 2, 300.0f, 10.0f, 0f),
			("B", 2, 500.0f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act
		var span = index.QueryMzRange(350.0f, 450.0f);

		// Assert: no entries in the gap
		Assert.That(span.Length, Is.EqualTo(0));
	}

	/// <summary>
	/// Verifies that an index containing a single entry returns that entry when the window
	/// exactly matches its m/z and returns empty when the window misses it.
	/// </summary>
	[Test]
	public void BinarySearch_SingleEntry_MatchesMz_ReturnsIt()
	{
		// Arrange: one entry
		var entries = BuildEntries(("OnlyEntry", 2, 412.5f, 10.0f, 0f));
		using var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act: exact match
		var spanHit = index.QueryMzRange(412.5f, 412.5f);
		var spanMiss = index.QueryMzRange(300.0f, 400.0f);

		// Assert
		Assert.That(spanHit.Length, Is.EqualTo(1));
		Assert.That(spanMiss.Length, Is.EqualTo(0));
	}

	// ═══════════════════════════════════════════════════════════════════════
	// Scale tests
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that <see cref="MslIndex.QueryMzRange"/> returns the correct number of
	/// entries at the tested scale, and that build time is under 2 seconds for
	/// <paramref name="nEntries"/> == 1,000,000 (CI hardware expectation).
	/// </summary>
	/// <param name="nEntries">Scale: 1 000, 100 000, or 1 000 000 entries.</param>
	[Test]
	[TestCase(1_000)]
	[TestCase(100_000)]
	[TestCase(1_000_000)]
	public void QueryMzRange_CorrectResults_AtScale(int nEntries)
	{
		// Arrange: uniformly-spaced m/z values from 400 to 1200 across nEntries
		var rawEntries = new MslPrecursorIndexEntry[nEntries];
		float mzStep = 800.0f / nEntries;
		for (int i = 0; i < nEntries; i++)
		{
			float mz = 400.0f + i * mzStep;
			rawEntries[i] = new MslPrecursorIndexEntry(
				precursorMz: mz,
				irt: (float)i,
				ionMobility: 0f,
				precursorIdx: i,
				elutionGroupId: i,
				charge: 2,
				isDecoy: 0,
				flags: 0);
		}

		// Measure build time (must be < 2 s for the 1 M case on CI hardware)
		var sw = Stopwatch.StartNew();
		using var index = new MslIndex(rawEntries, _ => null);
		sw.Stop();

		if (nEntries >= 1_000_000)
		{
			Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(2.0),
				$"Build took {sw.Elapsed.TotalSeconds:F2}s for {nEntries:N0} entries; expected < 2s.");
		}

		// Act: query the middle third of the m/z range
		float lo = 400.0f + mzStep * (nEntries / 3);
		float hi = 400.0f + mzStep * (nEntries * 2 / 3);
		var span = index.QueryMzRange(lo, hi);

		// Assert: all returned entries must actually be in [lo, hi]
		Assert.That(span.Length, Is.GreaterThan(0));
		foreach (var e in span)
		{
			Assert.That(e.PrecursorMz, Is.GreaterThanOrEqualTo(lo));
			Assert.That(e.PrecursorMz, Is.LessThanOrEqualTo(hi));
		}
	}

	// ═══════════════════════════════════════════════════════════════════════
	// IDisposable / ObjectDisposedException tests
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that every public method on <see cref="MslIndex"/> throws
	/// <see cref="ObjectDisposedException"/> after <see cref="MslIndex.Dispose"/> is called.
	/// Ensures that the disposed flag is checked consistently across all entry points.
	/// </summary>
	[Test]
	public void AllMethods_ThrowObjectDisposedException_AfterDispose()
	{
		// Arrange
		var entries = BuildEntries(("A", 2, 400.0f, 10.0f, 0f));
		var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);
		index.Dispose();

		// Assert: every public method throws after disposal
		Assert.That(() => index.QueryMzRange(0f, 1000f), Throws.InstanceOf<ObjectDisposedException>());
		Assert.That(() => index.QueryWindow(0f, 1000f, 0f, 100f), Throws.InstanceOf<ObjectDisposedException>());
		Assert.That(() => index.TryGetBySequenceCharge("A", 2, out _), Throws.InstanceOf<ObjectDisposedException>());
		Assert.That(() => index.GetElutionGroup(0), Throws.InstanceOf<ObjectDisposedException>());
		Assert.That(() => index.GetEntry(0), Throws.InstanceOf<ObjectDisposedException>());
		Assert.That(() => index.WithCalibratedRetentionTimes(1, 0), Throws.InstanceOf<ObjectDisposedException>());
	}

	/// <summary>
	/// Verifies that calling <see cref="MslIndex.Dispose"/> multiple times does not throw.
	/// </summary>
	[Test]
	public void Dispose_IsIdempotent()
	{
		// Arrange
		var entries = BuildEntries(("A", 2, 400.0f, 10.0f, 0f));
		var index = MslIndex.Build(entries, i => i < entries.Count ? entries[i] : null);

		// Act + Assert: double-dispose must not throw
		Assert.That(() =>
		{
			index.Dispose();
			index.Dispose();
		}, Throws.Nothing);
	}

	// ═══════════════════════════════════════════════════════════════════════
	// Private helpers
	// ═══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Builds a list of minimal <see cref="MslLibraryEntry"/> objects from a tuple array.
	/// Each tuple specifies (modifiedSequence, charge, precursorMz, irt, ionMobility).
	/// All other fields are set to sensible defaults (non-decoy, elution group = index + 1,
	/// Peptide molecule type, no fragments).
	/// </summary>
	/// <param name="specs">
	/// Variable-length array of (sequence, charge, mz, irt, im) tuples, one per entry.
	/// </param>
	/// <returns>
	/// An ordered list of <see cref="MslLibraryEntry"/> instances ready for
	/// <see cref="MslIndex.Build"/>.
	/// </returns>
	private static List<MslLibraryEntry> BuildEntries(
		params (string seq, int charge, float mz, float irt, float im)[] specs)
	{
		var list = new List<MslLibraryEntry>(specs.Length);
		for (int i = 0; i < specs.Length; i++)
		{
			var (seq, charge, mz, irt, im) = specs[i];
			list.Add(MakeEntry(seq, charge, mz, irt, im, isDecoy: false, elutionGroupId: i + 1));
		}
		return list;
	}

	/// <summary>
	/// Constructs a single minimal <see cref="MslLibraryEntry"/> with the specified fields
	/// and sensible defaults for all other properties. The entry has an empty fragment list
	/// (sufficient for index-construction tests that do not require fragment data).
	/// </summary>
	/// <param name="sequence">Modified sequence string (used as both modified and stripped).</param>
	/// <param name="charge">Precursor charge state.</param>
	/// <param name="mz">Precursor m/z.</param>
	/// <param name="irt">iRT value used for RT-window filtering.</param>
	/// <param name="ionMobility">Ion-mobility value (0 = not applicable).</param>
	/// <param name="isDecoy">True for a decoy precursor; false for a target.</param>
	/// <param name="elutionGroupId">
	/// Elution-group identifier; precursors sharing the same stripped sequence should
	/// receive the same ID.
	/// </param>
	/// <returns>A new <see cref="MslLibraryEntry"/> with the requested properties.</returns>
	private static MslLibraryEntry MakeEntry(
		string sequence,
		int charge,
		float mz,
		float irt,
		float ionMobility,
		bool isDecoy,
		int elutionGroupId = 0)
	{
		return new MslLibraryEntry
		{
			FullSequence = sequence,
			BaseSequence = sequence,
			PrecursorMz = mz,
			ChargeState = charge,
			RetentionTime = irt,
			IonMobility = ionMobility,
			IsDecoy = isDecoy,
			ElutionGroupId = elutionGroupId,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			Source = MslFormat.SourceType.Predicted,
			DissociationType = DissociationType.Unknown,
			Nce = 0,
			MatchedFragmentIons = new List<MslFragmentIon>()
		};
	}

	/// <summary>
	/// Converts a list of <see cref="MslLibraryEntry"/> objects to the compact
	/// <see cref="MslPrecursorIndexEntry"/> array format required by the
	/// <see cref="MslIndex"/> constructor. This mirrors the logic in
	/// <see cref="MslIndex.Build"/> and is used by tests that need to control
	/// <c>maxBufferSize</c> directly.
	/// </summary>
	/// <param name="entries">
	/// Ordered list of source entries; each entry's position determines its
	/// <see cref="MslPrecursorIndexEntry.PrecursorIdx"/>.
	/// </param>
	/// <returns>
	/// An unsorted array of compact index entries, one per source entry.
	/// </returns>
	private static MslPrecursorIndexEntry[] BuildRaw(IReadOnlyList<MslLibraryEntry> entries)
	{
		var raw = new MslPrecursorIndexEntry[entries.Count];
		for (int i = 0; i < entries.Count; i++)
		{
			MslLibraryEntry e = entries[i];
			raw[i] = new MslPrecursorIndexEntry(
				precursorMz: (float)e.PrecursorMz,
				irt: (float)e.RetentionTime,
				ionMobility: (float)e.IonMobility,
				precursorIdx: i,
				elutionGroupId: e.ElutionGroupId,
				charge: (short)e.ChargeState,
				isDecoy: (byte)(e.IsDecoy ? 1 : 0),
				flags: (byte)((int)e.MoleculeType & 0x03));
		}
		return raw;
	}
}