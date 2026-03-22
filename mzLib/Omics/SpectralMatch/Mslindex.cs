using System.Buffers;
using System.Collections.Concurrent;
using System.Runtime.InteropServices;
using System.Threading;

namespace Omics.SpectralMatch.MslSpectralLibrary;

// ─────────────────────────────────────────────────────────────────────────────
// MslPrecursorIndexEntry — 24 bytes, readonly struct
// ─────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Compact 24-byte summary of one library precursor, kept sorted by
/// <see cref="PrecursorMz"/> inside the index's hot-path array. The struct carries
/// only the fields needed for the two most common filter operations (m/z window and
/// RT/IM window); all other metadata lives in the full <see cref="MslLibraryEntry"/>
/// and is fetched on demand. At 24 bytes per entry, 1 million entries occupy 24 MB,
/// and the full human-proteome index (~5 million entries) occupies ~120 MB — easily
/// L3-cache-resident on modern hardware.
/// </summary>
[StructLayout(LayoutKind.Sequential)]
public readonly struct MslPrecursorIndexEntry : IComparable<MslPrecursorIndexEntry>
{
	// ── Filter fields (primary search keys) ─────────────────────────────

	/// <summary>
	/// Observed or predicted precursor m/z value (float32, 4 bytes).
	/// Primary sort key for the index array; all binary-search operations target this field.
	/// </summary>
	public readonly float PrecursorMz;

	/// <summary>
	/// Indexed Retention Time (iRT) value for this precursor (float32, 4 bytes).
	/// Used as the RT filter key in <see cref="MslIndex.QueryWindow"/>. When
	/// <see cref="MslIndex.WithCalibratedRetentionTimes"/> has been applied this field
	/// holds the run-specific calibrated RT rather than the raw iRT.
	/// </summary>
	public readonly float Irt;

	/// <summary>
	/// Ion-mobility value for this precursor (float32, 4 bytes).
	/// Used as the IM filter key when both IM bounds are non-zero in
	/// <see cref="MslIndex.QueryWindow"/>. Zero indicates ion mobility is not applicable
	/// or was not measured for this precursor.
	/// </summary>
	public readonly float IonMobility;

	// ── Index and grouping fields ────────────────────────────────────────

	/// <summary>
	/// Zero-based position of this precursor in the source entry array from which the
	/// index was built (int32, 4 bytes). Used by <see cref="MslIndex.GetEntry"/> to
	/// retrieve the full <see cref="MslLibraryEntry"/> via the entry-loader delegate.
	/// </summary>
	public readonly int PrecursorIdx;

	/// <summary>
	/// Elution group identifier shared by all charge states of the same peptide
	/// (int32, 4 bytes). Precursors with the same stripped sequence receive the same ID.
	/// Used by <see cref="MslIndex.GetElutionGroup"/> to retrieve all charge variants
	/// of a peptide in a single call.
	/// </summary>
	public readonly int ElutionGroupId;

	// ── Classification fields ────────────────────────────────────────────

	/// <summary>
	/// Precursor charge state (int16, 2 bytes).
	/// Stored here so callers can filter by charge without loading the full entry.
	/// </summary>
	public readonly short Charge;

	/// <summary>
	/// Decoy flag stored as a byte (0 = target, 1 = decoy) to avoid padding that a bool
	/// field would introduce in the sequential layout (1 byte).
	/// Evaluated in <see cref="MslIndex.QueryWindow"/> when <c>includeDecoys</c> is false.
	/// </summary>
	public readonly byte IsDecoy;

	/// <summary>
	/// Packed flags byte (1 byte). Bits 0–1 encode <see cref="MslFormat.MoleculeType"/>;
	/// bits 2–7 are reserved for future use and must be written as zero.
	/// </summary>
	public readonly byte Flags;

	// ── Constructor ──────────────────────────────────────────────────────

	/// <summary>
	/// Constructs an <see cref="MslPrecursorIndexEntry"/> from its individual components.
	/// All parameters map directly to the public readonly fields of the same name.
	/// </summary>
	/// <param name="precursorMz">
	/// Precursor m/z, used as the primary sort key. Must be a positive, finite float.
	/// </param>
	/// <param name="irt">
	/// iRT (or calibrated RT) value; used for RT-window filtering.
	/// </param>
	/// <param name="ionMobility">
	/// Ion-mobility value; pass 0 when IM is not available.
	/// </param>
	/// <param name="precursorIdx">
	/// Zero-based index into the source entry array for full-entry retrieval.
	/// </param>
	/// <param name="elutionGroupId">
	/// Elution group identifier; shared by all charge variants of the same peptide.
	/// </param>
	/// <param name="charge">Precursor charge state.</param>
	/// <param name="isDecoy">1 if this is a decoy precursor; 0 otherwise.</param>
	/// <param name="flags">Packed flags byte (MoleculeType in bits 0–1).</param>
	public MslPrecursorIndexEntry(
		float precursorMz,
		float irt,
		float ionMobility,
		int precursorIdx,
		int elutionGroupId,
		short charge,
		byte isDecoy,
		byte flags)
	{
		PrecursorMz = precursorMz;
		Irt = irt;
		IonMobility = ionMobility;
		PrecursorIdx = precursorIdx;
		ElutionGroupId = elutionGroupId;
		Charge = charge;
		IsDecoy = isDecoy;
		Flags = flags;
	}

	// ── IComparable<T> ───────────────────────────────────────────────────

	/// <summary>
	/// Compares this entry to <paramref name="other"/> by <see cref="PrecursorMz"/>.
	/// Required by <see cref="Array.Sort{T}(T[])"/> and by the binary-search lower-bound
	/// implementation inside <see cref="MslIndex"/>. NaN sorts to the end of the array
	/// per the float32 total ordering defined by <see cref="float.CompareTo(float)"/>.
	/// </summary>
	/// <param name="other">The entry to compare against.</param>
	/// <returns>
	/// A negative integer when this entry's m/z is smaller, zero when equal, and a
	/// positive integer when larger — matching the contract of <see cref="IComparable{T}"/>.
	/// </returns>
	public int CompareTo(MslPrecursorIndexEntry other)
		=> PrecursorMz.CompareTo(other.PrecursorMz);
}

// ─────────────────────────────────────────────────────────────────────────────
// MslWindowResults — disposable result wrapper
// ─────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Disposable, zero-extra-allocation result container returned by
/// <see cref="MslIndex.QueryWindow"/>. Internally backed by a pooled
/// <see cref="MslPrecursorIndexEntry"/> array rented from
/// <see cref="ArrayPool{T}.Shared"/>. Callers <b>must</b> call
/// <see cref="Dispose"/> (or use a <c>using</c> statement) after processing the
/// <see cref="Entries"/> span to return the buffer to the pool. Failing to dispose
/// will not cause data corruption but will increase GC pressure over time.
/// </summary>
public readonly struct MslWindowResults : IDisposable
{
	// ── Private state ────────────────────────────────────────────────────

	/// <summary>
	/// Pooled array rented from <see cref="ArrayPool{T}.Shared"/>. May be larger than
	/// <see cref="_count"/> because the pool returns the next-larger power-of-two size.
	/// Null only when <see cref="_count"/> is zero (empty result, no rental required).
	/// </summary>
	private readonly MslPrecursorIndexEntry[]? _buffer;

	/// <summary>
	/// Number of valid entries written into <see cref="_buffer"/>.
	/// May be less than <c>_buffer.Length</c> when the pool over-allocated.
	/// </summary>
	private readonly int _count;

	// ── Construction ─────────────────────────────────────────────────────

	/// <summary>
	/// Constructs an <see cref="MslWindowResults"/> wrapping a populated pooled buffer.
	/// Only <see cref="MslIndex"/> should call this constructor.
	/// </summary>
	/// <param name="buffer">
	/// Pooled array containing the matching entries at indices [0, count).
	/// May be null when count is 0.
	/// </param>
	/// <param name="count">
	/// Number of valid entries in <paramref name="buffer"/>. Must not exceed
	/// <c>buffer.Length</c> when buffer is non-null.
	/// </param>
	internal MslWindowResults(MslPrecursorIndexEntry[]? buffer, int count)
	{
		_buffer = buffer;
		_count = count;
	}

	// ── Public API ───────────────────────────────────────────────────────

	/// <summary>
	/// Read-only view over the matching index entries. The span is valid only while
	/// this <see cref="MslWindowResults"/> is alive — do not hold the span past
	/// <see cref="Dispose"/>. The entries are in ascending <see cref="MslPrecursorIndexEntry.PrecursorMz"/>
	/// order, reflecting the sort order of the backing index array.
	/// </summary>
	public ReadOnlySpan<MslPrecursorIndexEntry> Entries
		=> _buffer is null
			? ReadOnlySpan<MslPrecursorIndexEntry>.Empty
			: _buffer.AsSpan(0, _count);

	/// <summary>
	/// Number of matching index entries in this result set.
	/// Zero means no precursors passed all applied filters for the given query window.
	/// </summary>
	public int Count => _count;

	// ── IDisposable ──────────────────────────────────────────────────────

	/// <summary>
	/// Returns the rented backing array to <see cref="ArrayPool{T}.Shared"/>.
	/// Safe to call multiple times (idempotent due to null guard).
	/// Must be called once the caller has finished processing <see cref="Entries"/>;
	/// the span must not be used after this call.
	/// </summary>
	public void Dispose()
	{
		if (_buffer is not null)
			ArrayPool<MslPrecursorIndexEntry>.Shared.Return(_buffer);
	}
}

// ─────────────────────────────────────────────────────────────────────────────
// MslIndexStatistics
// ─────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Immutable snapshot of aggregate statistics computed from an <see cref="MslIndex"/>
/// instance. All fields are computed once at construction time from the sorted index
/// array and supporting dictionaries; the record is cheap to pass and log.
/// </summary>
public record MslIndexStatistics
{
	/// <summary>
	/// Total number of precursor entries in the index, including both targets and decoys.
	/// Equals <c>TargetPrecursors + DecoyPrecursors</c>.
	/// </summary>
	public int TotalPrecursors { get; init; }

	/// <summary>
	/// Number of non-decoy (target) precursor entries in the index.
	/// Computed by counting entries where <see cref="MslPrecursorIndexEntry.IsDecoy"/> == 0.
	/// </summary>
	public int TargetPrecursors { get; init; }

	/// <summary>
	/// Number of decoy precursor entries in the index.
	/// Computed by counting entries where <see cref="MslPrecursorIndexEntry.IsDecoy"/> != 0.
	/// </summary>
	public int DecoyPrecursors { get; init; }

	/// <summary>
	/// Smallest <see cref="MslPrecursorIndexEntry.PrecursorMz"/> value in the index.
	/// Equal to <c>float.PositiveInfinity</c> when the index is empty.
	/// </summary>
	public float MinPrecursorMz { get; init; }

	/// <summary>
	/// Largest <see cref="MslPrecursorIndexEntry.PrecursorMz"/> value in the index.
	/// Equal to <c>float.NegativeInfinity</c> when the index is empty.
	/// </summary>
	public float MaxPrecursorMz { get; init; }

	/// <summary>
	/// Smallest iRT value across all entries; useful for validating RT calibration bounds.
	/// Equal to <c>double.PositiveInfinity</c> when the index is empty.
	/// </summary>
	public double MinIrt { get; init; }

	/// <summary>
	/// Largest iRT value across all entries; useful for validating RT calibration bounds.
	/// Equal to <c>double.NegativeInfinity</c> when the index is empty.
	/// </summary>
	public double MaxIrt { get; init; }

	/// <summary>
	/// Number of distinct elution-group IDs in the index.
	/// Precursors sharing a stripped sequence share an elution-group ID; this value
	/// therefore represents the number of unique peptides (ignoring charge states).
	/// </summary>
	public int ElutionGroupCount { get; init; }

	/// <summary>
	/// Number of cache hits recorded by the LRU entry buffer since the index was created.
	/// A hit occurs when <see cref="MslIndex.GetEntry"/> finds the requested entry already
	/// cached and does not invoke the entry-loader delegate.
	/// </summary>
	public long LruHits { get; init; }

	/// <summary>
	/// Number of cache misses recorded by the LRU entry buffer since the index was created.
	/// A miss occurs when <see cref="MslIndex.GetEntry"/> must call the entry-loader delegate
	/// to fetch the full entry and add it to the cache.
	/// </summary>
	public long LruMisses { get; init; }
}

// ─────────────────────────────────────────────────────────────────────────────
// MslIndex
// ─────────────────────────────────────────────────────────────────────────────

/// <summary>
/// In-memory query engine for an .msl spectral library. Maintains a sorted array of
/// compact 24-byte <see cref="MslPrecursorIndexEntry"/> structs — one per precursor —
/// sorted by <see cref="MslPrecursorIndexEntry.PrecursorMz"/>. Provides:
/// <list type="bullet">
///   <item>O(log N) m/z-range queries via manual lower-bound binary search.</item>
///   <item>O(k) RT and IM filtering via linear scan over the (typically small) m/z slice.</item>
///   <item>O(1) DDA-style lookup by modified-sequence + charge via a dictionary.</item>
///   <item>O(k) elution-group queries via a pre-built group dictionary.</item>
///   <item>O(N) RT calibration producing a new immutable index instance.</item>
///   <item>LRU-buffered full-entry access to avoid repeated disk reads in index-only mode.</item>
/// </list>
/// <para>
/// This class is the hot path of DIA search — <see cref="QueryWindow"/> is called hundreds
/// of thousands of times per raw file. All code on that path is allocation-free or
/// pool-backed.
/// </para>
/// </summary>
public sealed class MslIndex : IDisposable
{
	// ── Private fields ───────────────────────────────────────────────────────

	/// <summary>
	/// The core sorted array of compact index entries, one per precursor.
	/// Sorted in ascending order by <see cref="MslPrecursorIndexEntry.PrecursorMz"/> so
	/// that binary search can locate the start of any m/z window in O(log N).
	/// </summary>
	private readonly MslPrecursorIndexEntry[] _byMz;

	/// <summary>
	/// Dictionary from "modifiedSequence/charge" key to the matching index entry,
	/// enabling O(1) DDA-style lookups via <see cref="TryGetBySequenceCharge"/>.
	/// Keys are case-sensitive; values reference entries in <see cref="_byMz"/>.
	/// </summary>
	private readonly Dictionary<string, MslPrecursorIndexEntry> _bySeqCharge;

	/// <summary>
	/// Dictionary from elution-group ID to the set of index entries sharing that ID.
	/// All charge states of the same peptide map to the same elution-group entry list,
	/// supporting rapid all-charge-variant queries via <see cref="GetElutionGroup"/>.
	/// </summary>
	private readonly Dictionary<int, List<MslPrecursorIndexEntry>> _byElutionGroup;

	/// <summary>
	/// Delegate that loads a full <see cref="MslLibraryEntry"/> by its zero-based
	/// precursor index. Supplied at construction time; typically backed by
	/// <see cref="MslLibrary.Entries"/> for full-load mode or
	/// <see cref="MslLibrary.LoadFragmentsOnDemand"/> for index-only mode.
	/// May return null for out-of-range indices.
	/// </summary>
	private readonly Func<int, MslLibraryEntry?> _entryLoader;

	/// <summary>
	/// Maximum number of <see cref="MslLibraryEntry"/> objects kept in the LRU cache.
	/// When the cache is full and a new entry must be added, the oldest entry (by insertion
	/// order) is evicted to keep memory bounded.
	/// </summary>
	private readonly int _maxBufferSize;

	// ── LRU buffer ───────────────────────────────────────────────────────────

	/// <summary>
	/// Thread-safe key-to-entry map that forms the "hot" half of the LRU cache.
	/// Populated on every <see cref="GetEntry"/> miss; checked on every call.
	/// </summary>
	private readonly ConcurrentDictionary<int, MslLibraryEntry> _lruCache;

	/// <summary>
	/// FIFO queue that tracks insertion order of keys in <see cref="_lruCache"/>.
	/// When eviction is required, the front key is dequeued and removed from the cache.
	/// </summary>
	private readonly ConcurrentQueue<int> _lruOrder;

	/// <summary>
	/// Running count of <see cref="GetEntry"/> calls that found the requested entry
	/// already present in <see cref="_lruCache"/>. Incremented with
	/// <see cref="Interlocked.Increment"/> to allow lock-free observation.
	/// </summary>
	private long _lruHits;

	/// <summary>
	/// Running count of <see cref="GetEntry"/> calls that did not find the requested
	/// entry in the cache and had to invoke <see cref="_entryLoader"/>. Incremented with
	/// <see cref="Interlocked.Increment"/>.
	/// </summary>
	private long _lruMisses;

	// ── Disposed flag ────────────────────────────────────────────────────────

	/// <summary>
	/// Set to 1 by <see cref="Dispose"/> to prevent use after disposal.
	/// Checked at the start of every public method.
	/// </summary>
	private int _disposed;

	// ── Constructor ──────────────────────────────────────────────────────────

	/// <summary>
	/// Builds an <see cref="MslIndex"/> from an array of compact precursor entries.
	/// <para>
	/// The constructor makes a defensive copy of <paramref name="entries"/> before sorting
	/// it by <see cref="MslPrecursorIndexEntry.PrecursorMz"/>, so the caller's array is
	/// never modified. After sorting, it builds the sequence/charge dictionary and the
	/// elution-group map in a single O(N) pass.
	/// </para>
	/// </summary>
	/// <param name="entries">
	/// Source array of compact precursor entries. A copy is made internally; the original
	/// array is not modified. Must not be null; may be empty (zero entries).
	/// </param>
	/// <param name="entryLoader">
	/// Delegate invoked by <see cref="GetEntry"/> when the requested full entry is not in
	/// the LRU cache. The integer parameter is the zero-based precursor index stored in
	/// <see cref="MslPrecursorIndexEntry.PrecursorIdx"/>. May return null for invalid indices.
	/// Must not be null.
	/// </param>
	/// <param name="maxBufferSize">
	/// Maximum number of full entries kept in the LRU cache. Default is 10,000.
	/// Larger values trade memory for reduced loader invocations in index-only mode.
	/// </param>
	/// <exception cref="ArgumentNullException">
	/// Thrown when <paramref name="entries"/> or <paramref name="entryLoader"/> is null.
	/// </exception>
	public MslIndex(
		MslPrecursorIndexEntry[] entries,
		Func<int, MslLibraryEntry?> entryLoader,
		int maxBufferSize = 10_000)
	{
		if (entries is null) throw new ArgumentNullException(nameof(entries));
		if (entryLoader is null) throw new ArgumentNullException(nameof(entryLoader));

		_entryLoader = entryLoader;
		_maxBufferSize = maxBufferSize;
		_lruCache = new ConcurrentDictionary<int, MslLibraryEntry>();
		_lruOrder = new ConcurrentQueue<int>();

		// Defensive copy then sort by m/z ascending
		_byMz = new MslPrecursorIndexEntry[entries.Length];
		entries.CopyTo(_byMz, 0);
		Array.Sort(_byMz);

		// Build the sequence/charge dictionary and elution-group map in one O(N) pass
		_bySeqCharge = new Dictionary<string, MslPrecursorIndexEntry>(_byMz.Length);
		_byElutionGroup = new Dictionary<int, List<MslPrecursorIndexEntry>>();

		for (int i = 0; i < _byMz.Length; i++)
		{
			// To populate the sequence/charge dictionary we need the modified sequence,
			// which is not stored in the compact struct. We obtain it via the loader.
			var entry = entryLoader(_byMz[i].PrecursorIdx);
			if (entry is not null)
			{
				string key = BuildSeqChargeKey(entry.FullSequence, _byMz[i].Charge);
				_bySeqCharge.TryAdd(key, _byMz[i]);
			}

			// Elution-group map: all charge states of the same peptide share one list
			int egId = _byMz[i].ElutionGroupId;
			if (!_byElutionGroup.TryGetValue(egId, out var list))
			{
				list = new List<MslPrecursorIndexEntry>(2);
				_byElutionGroup[egId] = list;
			}
			list.Add(_byMz[i]);
		}
	}

	// ── Static factory ───────────────────────────────────────────────────────

	/// <summary>
	/// Builds an <see cref="MslIndex"/> directly from a list of
	/// <see cref="MslLibraryEntry"/> objects. Typically called by
	/// <see cref="MslLibrary"/> after a full or index-only load.
	/// <para>
	/// Each entry is converted to an <see cref="MslPrecursorIndexEntry"/> struct using
	/// the entry's own fields; the loader delegate is wired to the list so that
	/// <see cref="GetEntry"/> can retrieve the full entry by index without re-reading disk.
	/// </para>
	/// </summary>
	/// <param name="entries">
	/// The ordered list of library entries as returned by <see cref="MslReader.Load"/>
	/// or <see cref="MslReader.LoadIndexOnly"/>. Must not be null; may be empty.
	/// </param>
	/// <param name="loader">
	/// Delegate used by <see cref="GetEntry"/> for LRU-cached full-entry access.
	/// Should return <c>entries[idx]</c> when <c>idx</c> is in range, null otherwise.
	/// Must not be null.
	/// </param>
	/// <returns>A fully-built <see cref="MslIndex"/> instance.</returns>
	/// <exception cref="ArgumentNullException">
	/// Thrown when <paramref name="entries"/> or <paramref name="loader"/> is null.
	/// </exception>
	public static MslIndex Build(
		IReadOnlyList<MslLibraryEntry> entries,
		Func<int, MslLibraryEntry?> loader)
	{
		if (entries is null) throw new ArgumentNullException(nameof(entries));
		if (loader is null) throw new ArgumentNullException(nameof(loader));

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

		return new MslIndex(raw, loader);
	}

	// ── Query: m/z range (zero-allocation) ───────────────────────────────────

	/// <summary>
	/// Returns all index entries whose <see cref="MslPrecursorIndexEntry.PrecursorMz"/>
	/// falls in the closed interval [<paramref name="mzLow"/>, <paramref name="mzHigh"/>].
	/// <para>
	/// This method is completely allocation-free: it returns a
	/// <see cref="ReadOnlySpan{T}"/> over the internal sorted array and performs no heap
	/// allocation on any call path. The span is valid until the next mutating operation on
	/// this index (which in practice never occurs — <see cref="MslIndex"/> is immutable
	/// after construction).
	/// </para>
	/// </summary>
	/// <param name="mzLow">
	/// Inclusive lower bound of the m/z window. Must be &lt;= <paramref name="mzHigh"/>.
	/// </param>
	/// <param name="mzHigh">
	/// Inclusive upper bound of the m/z window.
	/// </param>
	/// <returns>
	/// A <see cref="ReadOnlySpan{T}"/> slice of the sorted array containing only entries
	/// whose m/z is in [mzLow, mzHigh]. Returns an empty span when no entries match or
	/// when <paramref name="mzLow"/> &gt; <paramref name="mzHigh"/>.
	/// </returns>
	/// <exception cref="ObjectDisposedException">Thrown after <see cref="Dispose"/> is called.</exception>
	public ReadOnlySpan<MslPrecursorIndexEntry> QueryMzRange(float mzLow, float mzHigh)
	{
		ThrowIfDisposed();

		if (mzLow > mzHigh || _byMz.Length == 0)
			return ReadOnlySpan<MslPrecursorIndexEntry>.Empty;

		int lo = BinarySearchLowerBound(mzLow);
		if (lo >= _byMz.Length)
			return ReadOnlySpan<MslPrecursorIndexEntry>.Empty;

		// Find the exclusive upper bound: first index where PrecursorMz > mzHigh
		int hi = lo;
		while (hi < _byMz.Length && _byMz[hi].PrecursorMz <= mzHigh)
			hi++;

		int count = hi - lo;
		return count == 0
			? ReadOnlySpan<MslPrecursorIndexEntry>.Empty
			: _byMz.AsSpan(lo, count);
	}

	// ── Query: window (RT + IM filters, pooled result) ───────────────────────

	/// <summary>
	/// Returns all index entries that pass all of the following filters simultaneously:
	/// <list type="number">
	///   <item>
	///     <term>m/z window</term>
	///     <description>
	///       <see cref="MslPrecursorIndexEntry.PrecursorMz"/> in
	///       [<paramref name="mzLow"/>, <paramref name="mzHigh"/>].
	///     </description>
	///   </item>
	///   <item>
	///     <term>RT window</term>
	///     <description>
	///       <see cref="MslPrecursorIndexEntry.Irt"/> in
	///       [<paramref name="rtLow"/>, <paramref name="rtHigh"/>].
	///     </description>
	///   </item>
	///   <item>
	///     <term>IM window (optional)</term>
	///     <description>
	///       Applied only when both <paramref name="imLow"/> and <paramref name="imHigh"/>
	///       are non-zero. When either is zero the IM filter is skipped entirely.
	///     </description>
	///   </item>
	///   <item>
	///     <term>Decoy filter (optional)</term>
	///     <description>
	///       When <paramref name="includeDecoys"/> is false (the default), entries with
	///       <see cref="MslPrecursorIndexEntry.IsDecoy"/> != 0 are excluded.
	///     </description>
	///   </item>
	/// </list>
	/// <para>
	/// The result is backed by a pooled <see cref="MslPrecursorIndexEntry"/> array rented
	/// from <see cref="ArrayPool{T}.Shared"/>. The caller <b>must</b> dispose the returned
	/// <see cref="MslWindowResults"/> to return the buffer to the pool.
	/// </para>
	/// </summary>
	/// <param name="mzLow">Inclusive lower bound of the precursor m/z window.</param>
	/// <param name="mzHigh">Inclusive upper bound of the precursor m/z window.</param>
	/// <param name="rtLow">Inclusive lower bound of the RT (iRT) window.</param>
	/// <param name="rtHigh">Inclusive upper bound of the RT (iRT) window.</param>
	/// <param name="imLow">
	/// Inclusive lower bound of the ion-mobility window. Pass 0 (the default) to skip
	/// the IM filter entirely.
	/// </param>
	/// <param name="imHigh">
	/// Inclusive upper bound of the ion-mobility window. Pass 0 (the default) to skip
	/// the IM filter entirely.
	/// </param>
	/// <param name="includeDecoys">
	/// When true, decoy precursors are included in the results.
	/// When false (the default), decoy precursors are excluded.
	/// </param>
	/// <returns>
	/// An <see cref="MslWindowResults"/> instance whose <see cref="MslWindowResults.Entries"/>
	/// span contains all entries passing the applied filters in ascending m/z order.
	/// Must be disposed after use.
	/// </returns>
	/// <exception cref="ObjectDisposedException">Thrown after <see cref="Dispose"/> is called.</exception>
	public MslWindowResults QueryWindow(
		float mzLow, float mzHigh,
		float rtLow, float rtHigh,
		float imLow = 0f,
		float imHigh = 0f,
		bool includeDecoys = false)
	{
		ThrowIfDisposed();

		// Fast-path: empty index or inverted m/z range
		if (_byMz.Length == 0 || mzLow > mzHigh)
			return new MslWindowResults(null, 0);

		// Determine the m/z slice using binary search
		int start = BinarySearchLowerBound(mzLow);
		if (start >= _byMz.Length)
			return new MslWindowResults(null, 0);

		int end = start;
		while (end < _byMz.Length && _byMz[end].PrecursorMz <= mzHigh)
			end++;

		int sliceLen = end - start;
		if (sliceLen == 0)
			return new MslWindowResults(null, 0);

		// Determine whether the IM filter is active
		bool applyIm = imLow != 0f || imHigh != 0f;

		// Rent a buffer large enough to hold the entire m/z slice (worst-case output)
		MslPrecursorIndexEntry[] buffer = ArrayPool<MslPrecursorIndexEntry>.Shared.Rent(sliceLen);
		int writeIdx = 0;

		// Linear scan over the m/z slice, applying RT, IM, and decoy filters
		for (int i = start; i < end; i++)
		{
			ref readonly MslPrecursorIndexEntry e = ref _byMz[i];

			// Decoy filter: skip decoys when not requested
			if (!includeDecoys && e.IsDecoy != 0)
				continue;

			// RT filter
			if (e.Irt < rtLow || e.Irt > rtHigh)
				continue;

			// IM filter (only when both bounds are non-zero)
			if (applyIm && (e.IonMobility < imLow || e.IonMobility > imHigh))
				continue;

			buffer[writeIdx++] = e;
		}

		// If nothing survived the filters, return the buffer immediately
		if (writeIdx == 0)
		{
			ArrayPool<MslPrecursorIndexEntry>.Shared.Return(buffer);
			return new MslWindowResults(null, 0);
		}

		return new MslWindowResults(buffer, writeIdx);
	}

	// ── Query: DDA-style sequence/charge lookup ───────────────────────────────

	/// <summary>
	/// O(1) lookup of a precursor by its modified sequence and charge state.
	/// Uses the dictionary built during construction; key comparison is case-sensitive
	/// and uses ordinal string equality, matching the conventions of the mzLib notation
	/// (e.g. "PEPTM[Common Variable:Oxidation on M]IDE/2").
	/// </summary>
	/// <param name="modifiedSequence">
	/// Modified peptide sequence in mzLib bracket notation, e.g.
	/// "PEPTM[Common Variable:Oxidation on M]IDE". Case-sensitive.
	/// </param>
	/// <param name="charge">
	/// Precursor charge state as a positive integer. Combined with
	/// <paramref name="modifiedSequence"/> to form the lookup key.
	/// </param>
	/// <param name="entry">
	/// When the method returns true, contains the matching
	/// <see cref="MslPrecursorIndexEntry"/>; otherwise, contains a default value.
	/// </param>
	/// <returns>
	/// True when a matching entry was found; false when no entry with the given
	/// sequence and charge is present in the index.
	/// </returns>
	/// <exception cref="ObjectDisposedException">Thrown after <see cref="Dispose"/> is called.</exception>
	public bool TryGetBySequenceCharge(
		string modifiedSequence,
		int charge,
		out MslPrecursorIndexEntry entry)
	{
		ThrowIfDisposed();
		return _bySeqCharge.TryGetValue(
			BuildSeqChargeKey(modifiedSequence, (short)charge),
			out entry);
	}

	// ── Query: elution group ─────────────────────────────────────────────────

	/// <summary>
	/// Returns all index entries belonging to the given elution group. An elution group
	/// collects all charge states of the same peptide (same stripped sequence) so that the
	/// DIA search engine can consider them together during precursor scoring.
	/// </summary>
	/// <param name="elutionGroupId">
	/// The integer elution-group ID as stored in
	/// <see cref="MslPrecursorIndexEntry.ElutionGroupId"/>. Use 0 as the sentinel for
	/// ungrouped precursors when the library was built without elution-group assignments.
	/// </param>
	/// <returns>
	/// A read-only span over the entries in the group in ascending m/z order.
	/// Returns an empty span when no entries belong to the specified group.
	/// </returns>
	/// <exception cref="ObjectDisposedException">Thrown after <see cref="Dispose"/> is called.</exception>
	public ReadOnlySpan<MslPrecursorIndexEntry> GetElutionGroup(int elutionGroupId)
	{
		ThrowIfDisposed();

		return _byElutionGroup.TryGetValue(elutionGroupId, out var list)
			? list.ToArray().AsSpan()         // defensive copy so the caller can't mutate the list
			: ReadOnlySpan<MslPrecursorIndexEntry>.Empty;
	}

	// ── Full entry access (LRU-buffered) ─────────────────────────────────────

	/// <summary>
	/// Loads the full <see cref="MslLibraryEntry"/> for the precursor at the given
	/// <paramref name="precursorIdx"/>, using a bounded LRU cache to avoid repeated
	/// delegate invocations for recently accessed entries.
	/// <para>
	/// On a cache hit the cached entry is returned immediately with no delegate call.
	/// On a cache miss the entry-loader delegate supplied at construction time is invoked;
	/// if the result is non-null it is added to the cache. When the cache is full
	/// (at <c>maxBufferSize</c> entries) the oldest entry is evicted first.
	/// </para>
	/// </summary>
	/// <param name="precursorIdx">
	/// Zero-based precursor index as stored in
	/// <see cref="MslPrecursorIndexEntry.PrecursorIdx"/>. Values outside the valid range
	/// cause the delegate to return null, which is propagated to the caller.
	/// </param>
	/// <returns>
	/// The full <see cref="MslLibraryEntry"/> when the index is in range and the loader
	/// returns a non-null value; null otherwise.
	/// </returns>
	/// <exception cref="ObjectDisposedException">Thrown after <see cref="Dispose"/> is called.</exception>
	public MslLibraryEntry? GetEntry(int precursorIdx)
	{
		ThrowIfDisposed();

		// Cache hit: return the cached entry without calling the loader
		if (_lruCache.TryGetValue(precursorIdx, out MslLibraryEntry? cached))
		{
			Interlocked.Increment(ref _lruHits);
			return cached;
		}

		// Cache miss: invoke the loader
		Interlocked.Increment(ref _lruMisses);
		MslLibraryEntry? loaded = _entryLoader(precursorIdx);

		if (loaded is null)
			return null;

		// Evict the oldest entry when the cache is at capacity
		if (_lruCache.Count >= _maxBufferSize)
		{
			if (_lruOrder.TryDequeue(out int oldestKey))
				_lruCache.TryRemove(oldestKey, out _);
		}

		// Add the new entry to the cache and record its insertion order
		if (_lruCache.TryAdd(precursorIdx, loaded))
			_lruOrder.Enqueue(precursorIdx);

		return loaded;
	}

	// ── RT calibration ───────────────────────────────────────────────────────

	/// <summary>
	/// Applies a linear iRT → run-RT transform to every entry's
	/// <see cref="MslPrecursorIndexEntry.Irt"/> field and returns a brand-new
	/// <see cref="MslIndex"/> containing the transformed entries.
	/// The formula applied to each entry is:
	/// <code>calibratedRT = (float)(slope * entry.RetentionTime + intercept)</code>
	/// <para>
	/// The original index is not modified in any way. The new index shares the same
	/// entry-loader delegate and LRU configuration but starts with an empty LRU cache.
	/// </para>
	/// </summary>
	/// <param name="slope">
	/// Multiplicative scale factor for the iRT → run-RT conversion.
	/// A slope of 1.0 combined with an intercept of 0.0 produces an identity transform
	/// (output values are numerically identical to the input).
	/// </param>
	/// <param name="intercept">
	/// Additive offset applied after scaling. Represents the iRT value corresponding to
	/// a run-RT of zero.
	/// </param>
	/// <returns>
	/// A new <see cref="MslIndex"/> whose entries have calibrated RT values and whose
	/// m/z sort order is preserved.
	/// </returns>
	/// <exception cref="ObjectDisposedException">Thrown after <see cref="Dispose"/> is called.</exception>
	public MslIndex WithCalibratedRetentionTimes(double slope, double intercept)
	{
		ThrowIfDisposed();

		// Build a new entry array with transformed RetentionTime values; all other fields are unchanged
		var calibrated = new MslPrecursorIndexEntry[_byMz.Length];
		for (int i = 0; i < _byMz.Length; i++)
		{
			ref readonly MslPrecursorIndexEntry src = ref _byMz[i];
			float newIrt = (float)(slope * src.Irt + intercept);
			calibrated[i] = new MslPrecursorIndexEntry(
				precursorMz: src.PrecursorMz,
				irt: newIrt,
				ionMobility: src.IonMobility,
				precursorIdx: src.PrecursorIdx,
				elutionGroupId: src.ElutionGroupId,
				charge: src.Charge,
				isDecoy: src.IsDecoy,
				flags: src.Flags);
		}

		// The new index is already sorted (only RetentionTime changed, PrecursorMz is unchanged)
		// but we pass through the constructor to rebuild dictionaries correctly.
		return new MslIndex(calibrated, _entryLoader, _maxBufferSize);
	}

	// ── Statistics ───────────────────────────────────────────────────────────

	/// <summary>
	/// Computes and returns an <see cref="MslIndexStatistics"/> snapshot from the current
	/// state of the index. Performs a single O(N) pass over <see cref="_byMz"/> to gather
	/// aggregate values; LRU counters are read with <see cref="Interlocked.Read"/> to
	/// ensure a consistent snapshot in multi-threaded scenarios.
	/// </summary>
	/// <returns>
	/// An immutable <see cref="MslIndexStatistics"/> record reflecting the current
	/// state of the index and its LRU cache.
	/// </returns>
	public MslIndexStatistics GetStatistics()
	{
		int totalPrecursors = _byMz.Length;
		int targetPrecursors = 0;
		int decoyPrecursors = 0;
		float minMz = float.PositiveInfinity;
		float maxMz = float.NegativeInfinity;
		double minIrt = double.PositiveInfinity;
		double maxIrt = double.NegativeInfinity;

		for (int i = 0; i < _byMz.Length; i++)
		{
			ref readonly MslPrecursorIndexEntry e = ref _byMz[i];

			if (e.IsDecoy != 0) decoyPrecursors++;
			else targetPrecursors++;

			if (e.PrecursorMz < minMz) minMz = e.PrecursorMz;
			if (e.PrecursorMz > maxMz) maxMz = e.PrecursorMz;
			if (e.Irt < minIrt) minIrt = e.Irt;
			if (e.Irt > maxIrt) maxIrt = e.Irt;
		}

		return new MslIndexStatistics
		{
			TotalPrecursors = totalPrecursors,
			TargetPrecursors = targetPrecursors,
			DecoyPrecursors = decoyPrecursors,
			MinPrecursorMz = totalPrecursors > 0 ? minMz : float.PositiveInfinity,
			MaxPrecursorMz = totalPrecursors > 0 ? maxMz : float.NegativeInfinity,
			MinIrt = totalPrecursors > 0 ? minIrt : double.PositiveInfinity,
			MaxIrt = totalPrecursors > 0 ? maxIrt : double.NegativeInfinity,
			ElutionGroupCount = _byElutionGroup.Count,
			LruHits = Interlocked.Read(ref _lruHits),
			LruMisses = Interlocked.Read(ref _lruMisses),
		};
	}

	// ── IDisposable ───────────────────────────────────────────────────────────

	/// <summary>
	/// Clears the LRU cache and marks this instance as disposed. All subsequent calls
	/// to public methods will throw <see cref="ObjectDisposedException"/>.
	/// Safe to call multiple times (idempotent).
	/// </summary>
	public void Dispose()
	{
		if (Interlocked.Exchange(ref _disposed, 1) != 0)
			return;  // Already disposed

		_lruCache.Clear();

		// Drain the eviction queue (no finalizer needed; all resources are managed)
		while (_lruOrder.TryDequeue(out _)) { }
	}

	// ── Private helpers ───────────────────────────────────────────────────────

	/// <summary>
	/// Returns the index of the first element in <see cref="_byMz"/> whose
	/// <see cref="MslPrecursorIndexEntry.PrecursorMz"/> is greater than or equal to
	/// <paramref name="mzLow"/>. Returns <c>_byMz.Length</c> when all entries are below
	/// the lower bound.
	/// <para>
	/// This is a standard lower-bound binary search (analogous to <c>std::lower_bound</c>
	/// in C++). <see cref="Array.BinarySearch{T}"/> is intentionally avoided because its
	/// contract for duplicate keys is implementation-defined; the manual algorithm always
	/// returns the leftmost matching position.
	/// </para>
	/// </summary>
	/// <param name="mzLow">
	/// The inclusive lower bound to search for. All returned positions have
	/// <c>_byMz[pos].PrecursorMz >= mzLow</c>.
	/// </param>
	/// <returns>
	/// The smallest index <c>i</c> such that <c>_byMz[i].PrecursorMz >= mzLow</c>,
	/// or <c>_byMz.Length</c> when no such index exists.
	/// </returns>
	private int BinarySearchLowerBound(float mzLow)
	{
		int lo = 0;
		int hi = _byMz.Length;

		while (lo < hi)
		{
			// Use unsigned right-shift to avoid overflow for large mid values
			int mid = lo + ((hi - lo) >> 1);

			if (_byMz[mid].PrecursorMz < mzLow)
				lo = mid + 1;  // Answer is to the right of mid
			else
				hi = mid;      // Answer is at mid or to the left of mid
		}

		return lo;
	}

	/// <summary>
	/// Builds the dictionary key used by <see cref="_bySeqCharge"/> from a modified sequence
	/// string and a charge value. The key format is <c>"modifiedSequence/charge"</c>,
	/// e.g. <c>"PEPTM[Common Variable:Oxidation on M]IDE/2"</c>.
	/// This format is consistent with the DDA lookup convention used throughout mzLib.
	/// </summary>
	/// <param name="modifiedSequence">
	/// Modified sequence in mzLib bracket notation. Must not be null.
	/// </param>
	/// <param name="charge">
	/// Precursor charge state stored in the index entry (int16).
	/// </param>
	/// <returns>
	/// The string key for dictionary lookup, e.g. <c>"PEPTIDE/2"</c>.
	/// </returns>
	private static string BuildSeqChargeKey(string modifiedSequence, short charge)
		=> $"{modifiedSequence}/{charge}";

	/// <summary>
	/// Throws <see cref="ObjectDisposedException"/> when this instance has been disposed.
	/// Called at the start of every public method that accesses internal state.
	/// </summary>
	private void ThrowIfDisposed()
	{
		if (_disposed != 0)
			throw new ObjectDisposedException(nameof(MslIndex));
	}
}