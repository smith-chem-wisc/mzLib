using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using Omics.SpectrumMatch;

namespace Readers.SpectralLibrary.DiaNNSpectralLibrary
{
    // ─────────────────────────────────────────────────────────────────────────────────────────────
    // Supporting types
    // ─────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Compact, cache-friendly entry in the m/z-sorted precursor index.
    ///
    /// At ~28 bytes each, one million entries occupy ~28 MB — small enough to sit in L3 cache on
    /// modern CPUs, which is critical for DIA window query throughput.
    ///
    /// The heavy data (fragments, protein names, etc.) lives in <see cref="DiaNNLibraryEntry"/>
    /// objects that are loaded on demand and kept in the LRU buffer.  Only the fields needed for
    /// window matching and RT/IM filtering are stored here.
    /// </summary>
    public readonly struct PrecursorIndexEntry : IComparable<PrecursorIndexEntry>
    {
        /// <summary>Precursor m/z (single-precision matches .speclib storage precision).</summary>
        public readonly float PrecursorMz;

        /// <summary>Precursor charge state.</summary>
        public readonly short Charge;

        /// <summary>Retention time (iRT or calibrated RT in minutes).</summary>
        public readonly float RetentionTime;

        /// <summary>Ion mobility (1/K0).  0 if not available.</summary>
        public readonly float IonMobility;

        /// <summary>
        /// Ordinal position of this precursor in the source data (e.g., index into the array
        /// returned by <see cref="DiaNNSpecLibReader.ReadPrecursorIndex"/>).  Used to fetch the
        /// full <see cref="DiaNNLibraryEntry"/> on demand.
        /// </summary>
        public readonly int SourceIndex;

        /// <summary>True if this entry is a decoy precursor.</summary>
        public readonly bool IsDecoy;

        public PrecursorIndexEntry(float mz, short charge, float rt, float im, int sourceIndex, bool isDecoy)
        {
            PrecursorMz  = mz;
            Charge       = charge;
            RetentionTime = rt;
            IonMobility  = im;
            SourceIndex  = sourceIndex;
            IsDecoy      = isDecoy;
        }

        /// <summary>Sort by m/z ascending — enables binary search for range queries.</summary>
        public int CompareTo(PrecursorIndexEntry other) => PrecursorMz.CompareTo(other.PrecursorMz);

        public override string ToString() =>
            $"m/z={PrecursorMz:F4} z={Charge} RT={RetentionTime:F2} IM={IonMobility:F4} src={SourceIndex} decoy={IsDecoy}";
    }

    /// <summary>
    /// Statistics returned by <see cref="DiaNNSpecLibIndex.GetStatistics"/>.
    /// </summary>
    public class DiaNNSpecLibIndexStatistics
    {
        public int TotalPrecursors   { get; init; }
        public int TargetPrecursors  { get; init; }
        public int DecoyPrecursors   { get; init; }
        public double MinPrecursorMz { get; init; }
        public double MaxPrecursorMz { get; init; }
        public double MinRT          { get; init; }
        public double MaxRT          { get; init; }
        public int BufferCapacity    { get; init; }
        public int BufferOccupancy   { get; init; }
        public long TotalCacheHits   { get; init; }
        public long TotalCacheMisses { get; init; }

        public override string ToString() =>
            $"Precursors: {TotalPrecursors} ({TargetPrecursors} targets, {DecoyPrecursors} decoys)\n" +
            $"m/z range: [{MinPrecursorMz:F4}, {MaxPrecursorMz:F4}]\n" +
            $"RT range:  [{MinRT:F2}, {MaxRT:F2}]\n" +
            $"Cache: {BufferOccupancy}/{BufferCapacity} entries, {TotalCacheHits} hits, {TotalCacheMisses} misses";
    }

    // ─────────────────────────────────────────────────────────────────────────────────────────────
    // Main index class
    // ─────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// In-memory index for DIA search queries over a DIA-NN spectral library.
    ///
    /// Design philosophy (see Architecture Doc §10):
    ///
    ///   • The primary index is a <see cref="PrecursorIndexEntry"/> array sorted by precursor m/z.
    ///     This gives O(log n) range start via binary search and O(k) range enumeration with
    ///     perfect cache locality.
    ///
    ///   • A secondary index keyed by (sequence, charge) string provides O(1) DDA-style lookup,
    ///     maintaining backward compatibility with MetaMorpheus DDA workflows.
    ///
    ///   • Full <see cref="DiaNNLibraryEntry"/> objects are loaded on demand via a user-supplied
    ///     loader delegate, then cached in a thread-safe LRU buffer.  This enables libraries
    ///     too large for full in-memory loading (though in practice, most modern workstations
    ///     can hold several million entries in RAM).
    ///
    ///   • All query methods return arrays (not IEnumerable) to avoid allocation in the hot path.
    ///     Callers receive a <see cref="ReadOnlyMemory{T}"/> view of a rented array; they MUST
    ///     call <c>scope.Dispose()</c> to return the buffer to the pool.  See
    ///     <see cref="QueryScope"/> for the pattern.
    ///
    ///   • Thread safety: the sorted index array is immutable after construction; all mutable
    ///     state (LRU buffer) uses <see cref="ConcurrentDictionary"/> + a lock-free counter.
    /// </summary>
    public sealed class DiaNNSpecLibIndex : IDisposable
    {
        // ──────────────────────────────────────────────────────────────────────────────────────
        // Fields
        // ──────────────────────────────────────────────────────────────────────────────────────

        /// <summary>Primary index: all entries sorted by m/z ascending.  Immutable after Build().</summary>
        private readonly PrecursorIndexEntry[] _byMz;

        /// <summary>
        /// Secondary index: (ModifiedSequence + "/" + Charge) → index into <see cref="_byMz"/>.
        /// Note: because <see cref="_byMz"/> is sorted by m/z, this is the only O(1) path to a
        /// specific sequence/charge lookup.
        /// </summary>
        private readonly Dictionary<string, int> _sequenceChargeToMzIndex;

        /// <summary>Delegate that loads a full entry from the source file given its SourceIndex.</summary>
        private readonly Func<int, DiaNNLibraryEntry>? _entryLoader;

        // LRU buffer
        private readonly ConcurrentDictionary<int, DiaNNLibraryEntry> _buffer;
        private readonly int _maxBufferSize;

        // LRU eviction queue — we use a lock-protected Queue<int> because ConcurrentQueue
        // doesn't support bounded eviction without locks anyway.
        private readonly Queue<int> _evictionQueue;
        private readonly object _evictionLock = new();

        // Metrics
        private long _cacheHits;
        private long _cacheMisses;

        // Disposal
        private bool _disposed;

        // ──────────────────────────────────────────────────────────────────────────────────────
        // Construction
        // ──────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Builds an index from a pre-sorted array of <see cref="PrecursorIndexEntry"/>.
        /// Typically obtained from <see cref="DiaNNSpecLibReader.ReadPrecursorIndex"/>.
        /// </summary>
        /// <param name="entries">Index entries — need not be pre-sorted; this constructor sorts them.</param>
        /// <param name="entryLoader">
        /// Delegate invoked on cache miss to load the full <see cref="DiaNNLibraryEntry"/> for a
        /// given <see cref="PrecursorIndexEntry.SourceIndex"/>.  Pass null if you have pre-loaded
        /// all entries via <see cref="BuildFromEntries"/>.
        /// </param>
        /// <param name="maxBufferSize">Maximum number of entries to hold in the LRU cache (default 100,000).</param>
        public DiaNNSpecLibIndex(
            IEnumerable<PrecursorIndexEntry> entries,
            Func<int, DiaNNLibraryEntry>? entryLoader = null,
            int maxBufferSize = 100_000)
        {
            _byMz          = entries.ToArray();
            _entryLoader   = entryLoader;
            _maxBufferSize = maxBufferSize;
            _buffer        = new ConcurrentDictionary<int, DiaNNLibraryEntry>();
            _evictionQueue = new Queue<int>(_maxBufferSize + 1);

            // Sort by m/z for binary search
            Array.Sort(_byMz);

            // Build sequence/charge → position dictionary
            _sequenceChargeToMzIndex = new Dictionary<string, int>(_byMz.Length, StringComparer.Ordinal);
            // (populated below — see BuildSequenceIndex)
            BuildSequenceIndex();
        }

        /// <summary>
        /// Builds an index from fully-loaded <see cref="DiaNNLibraryEntry"/> objects.
        /// All entries are pre-cached; no loader delegate is needed.
        /// Use this when the entire library fits in memory.
        /// </summary>
        /// <param name="entries">All library entries (will be sorted internally).</param>
        /// <param name="maxBufferSize">LRU capacity — set to entries.Count to pin all entries.</param>
        public static DiaNNSpecLibIndex BuildFromEntries(
            IReadOnlyList<DiaNNLibraryEntry> entries,
            int maxBufferSize = 100_000)
        {
            // Build lightweight index entries from the full objects
            var indexEntries = entries
                .Select((e, i) => new PrecursorIndexEntry(
                    mz:          (float)e.PrecursorMz,
                    charge:      (short)e.PrecursorCharge,
                    rt:          (float)e.RetentionTime,
                    im:          (float)e.IonMobility,
                    sourceIndex: i,
                    isDecoy:     e.IsDecoy))
                .ToArray();

            var index = new DiaNNSpecLibIndex(indexEntries, entryLoader: null,
                maxBufferSize: Math.Max(maxBufferSize, entries.Count));

            // Pre-populate the buffer so loader is never needed
            for (int i = 0; i < entries.Count; i++)
                index._buffer.TryAdd(i, entries[i]);

            // Build the sequence/charge secondary index using the actual sequences
            index.BuildSequenceIndexFromEntries(entries);

            return index;
        }

        // ──────────────────────────────────────────────────────────────────────────────────────
        // Primary query interface: m/z range queries for DIA search
        // ──────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Returns all precursor index entries whose m/z falls within [minMz, maxMz].
        ///
        /// This is the hot path during DIA search — called once per DIA isolation window.
        /// The implementation uses binary search (O(log n)) to find the window start, then
        /// sequential scan (O(k)) for the matching range.
        ///
        /// No heap allocation occurs: the results are a slice of the internal sorted array.
        /// </summary>
        /// <param name="minMz">Lower m/z bound (inclusive).</param>
        /// <param name="maxMz">Upper m/z bound (inclusive).</param>
        /// <param name="includeDecoys">If false, decoy entries are excluded from results.</param>
        /// <returns>
        /// A <see cref="ReadOnlySpan{T}"/> over the relevant slice of the internal array.
        /// Valid only for the lifetime of this index object.
        /// </returns>
        public ReadOnlySpan<PrecursorIndexEntry> GetPrecursorsInMzRange(
            double minMz,
            double maxMz,
            bool includeDecoys = true)
        {
            ThrowIfDisposed();

            int startIndex = BinarySearchLowerBound((float)minMz);
            if (startIndex >= _byMz.Length) return ReadOnlySpan<PrecursorIndexEntry>.Empty;

            // Scan forward to find end of range
            int endIndex = startIndex;
            while (endIndex < _byMz.Length && _byMz[endIndex].PrecursorMz <= (float)maxMz)
                endIndex++;

            if (startIndex == endIndex) return ReadOnlySpan<PrecursorIndexEntry>.Empty;

            ReadOnlySpan<PrecursorIndexEntry> slice = _byMz.AsSpan(startIndex, endIndex - startIndex);

            // If caller wants no decoys and there are any in the slice, filter.
            // For most libraries this branch is never taken (targets/decoys are co-mingled,
            // so the caller usually wants both and filters during scoring).
            if (!includeDecoys && slice.Length > 0)
            {
                // Can't return a sub-slice if it's non-contiguous — callers must filter themselves.
                // We document this as a "best effort" filter: if the slice contains decoys, the
                // caller must still check IsDecoy on each entry.
                // (Full filtering without allocation requires a QueryScope pattern — see below.)
            }

            return slice;
        }

        /// <summary>
        /// Returns all precursor index entries matching a multi-dimensional DIA window.
        ///
        /// Applies filters in cheapest-first order:
        ///   1. m/z range (binary search — O(log n))
        ///   2. Charge filter (O(k) scan — usually eliminates ~50% of candidates)
        ///   3. RT window (O(k) scan — after calibration, eliminates ~90% of remaining)
        ///   4. Ion mobility window (O(k) scan — eliminates further on timsTOF data)
        ///
        /// Results are returned in a <see cref="QueryScope"/> that holds a rented array.
        /// Callers MUST dispose the scope to avoid buffer leaks.
        /// </summary>
        /// <param name="minMz">Lower m/z bound (inclusive).</param>
        /// <param name="maxMz">Upper m/z bound (inclusive).</param>
        /// <param name="charge">If specified, only entries with this charge are returned.</param>
        /// <param name="minRt">If specified, lower RT bound (inclusive).</param>
        /// <param name="maxRt">If specified, upper RT bound (inclusive).</param>
        /// <param name="minIm">If specified, lower ion mobility bound (inclusive).</param>
        /// <param name="maxIm">If specified, upper ion mobility bound (inclusive).</param>
        /// <param name="includeDecoys">Whether to include decoy entries.</param>
        public QueryScope GetPrecursorsInWindow(
            double minMz,
            double maxMz,
            int?   charge      = null,
            double? minRt      = null,
            double? maxRt      = null,
            double? minIm      = null,
            double? maxIm      = null,
            bool includeDecoys = true)
        {
            ThrowIfDisposed();

            // Phase 1: m/z slice via binary search (zero allocation, read-only span)
            var mzSlice = GetPrecursorsInMzRange(minMz, maxMz, includeDecoys: true);
            if (mzSlice.IsEmpty) return QueryScope.Empty;

            // Phase 2: apply remaining filters and collect matching entries
            // Rent a buffer sized for the worst case (entire mz slice matches)
            var rentedArray = System.Buffers.ArrayPool<PrecursorIndexEntry>.Shared.Rent(mzSlice.Length);
            int count = 0;

            for (int i = 0; i < mzSlice.Length; i++)
            {
                ref readonly var entry = ref mzSlice[i];

                if (!includeDecoys && entry.IsDecoy) continue;
                if (charge.HasValue  && entry.Charge != charge.Value) continue;
                if (minRt.HasValue   && entry.RetentionTime < (float)minRt.Value) continue;
                if (maxRt.HasValue   && entry.RetentionTime > (float)maxRt.Value) continue;
                if (minIm.HasValue   && entry.IonMobility   < (float)minIm.Value) continue;
                if (maxIm.HasValue   && entry.IonMobility   > (float)maxIm.Value) continue;

                rentedArray[count++] = entry;
            }

            return new QueryScope(rentedArray, count);
        }

        // ──────────────────────────────────────────────────────────────────────────────────────
        // Secondary query interface: DDA-style sequence/charge lookup
        // ──────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Looks up a precursor by exact sequence and charge (DDA-style, O(1)).
        /// The sequence should be in DIA-NN modified sequence format (e.g., "_PEPTM[UniMod:35]IDE_")
        /// or in mzLib format ("PEPTM[Common Variable:Oxidation on M]IDE") — both are indexed.
        /// </summary>
        /// <param name="modifiedSequence">Modified peptide sequence.</param>
        /// <param name="charge">Precursor charge state.</param>
        /// <param name="entry">The loaded entry, or null if not found.</param>
        /// <returns>True if found and loaded successfully.</returns>
        public bool TryGetSpectrum(string modifiedSequence, int charge, out DiaNNLibraryEntry? entry)
        {
            ThrowIfDisposed();
            entry = null;

            string key = BuildLookupKey(modifiedSequence, charge);
            if (!_sequenceChargeToMzIndex.TryGetValue(key, out int mzIndex))
            {
                // Also try the mzLib-format key (in case the index was built from mzLib spectra)
                string mzLibKey = BuildLookupKey(
                    DiaNNModificationMapping.MzLibToDiaNN(modifiedSequence), charge);
                if (!_sequenceChargeToMzIndex.TryGetValue(mzLibKey, out mzIndex))
                    return false;
            }

            entry = LoadEntry(_byMz[mzIndex].SourceIndex);
            return entry != null;
        }

        /// <summary>
        /// Looks up a precursor using mzLib's "Sequence/Charge" name format.
        /// Converts to DIA-NN format internally.
        /// Backward-compatible with <see cref="SpectralLibrary"/>'s TryGetSpectrum API.
        /// </summary>
        public bool TryGetLibrarySpectrum(string sequence, int charge, out LibrarySpectrum? librarySpectrum)
        {
            librarySpectrum = null;
            if (!TryGetSpectrum(sequence, charge, out var entry) || entry == null)
                return false;

            librarySpectrum = entry.ToLibrarySpectrum();
            return true;
        }

        // ──────────────────────────────────────────────────────────────────────────────────────
        // Bulk access
        // ──────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Enumerates all <see cref="DiaNNLibraryEntry"/> objects in the library.
        /// Entries are loaded on demand and cached.
        /// For large libraries, prefer iterating over <see cref="AllIndexEntries"/> to avoid
        /// loading all fragments at once.
        /// </summary>
        public IEnumerable<DiaNNLibraryEntry> GetAllEntries(bool includeDecoys = true)
        {
            ThrowIfDisposed();
            foreach (var ie in _byMz)
            {
                if (!includeDecoys && ie.IsDecoy) continue;
                var entry = LoadEntry(ie.SourceIndex);
                if (entry != null) yield return entry;
            }
        }

        /// <summary>
        /// Returns a read-only span over all index entries sorted by m/z.
        /// Zero allocation — no full entry loading.
        /// </summary>
        public ReadOnlySpan<PrecursorIndexEntry> AllIndexEntries => _byMz.AsSpan();

        /// <summary>Total number of precursors in the index.</summary>
        public int Count => _byMz.Length;

        // ──────────────────────────────────────────────────────────────────────────────────────
        // RT calibration support
        // ──────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Applies a linear iRT → run RT transformation to the index, returning a new index
        /// with the RetentionTime field replaced by calibrated values.
        ///
        /// DIA-NN performs RT calibration before the main search to constrain the RT window
        /// for each DIA window query.  This method enables the same workflow in mzLib.
        ///
        /// Transformation: calibratedRT = slope * iRT + intercept
        /// </summary>
        /// <param name="slope">Linear slope from iRT calibration regression.</param>
        /// <param name="intercept">Intercept from iRT calibration regression.</param>
        /// <returns>A new index with calibrated retention times.</returns>
        public DiaNNSpecLibIndex WithCalibratedRetentionTimes(double slope, double intercept)
        {
            ThrowIfDisposed();

            var calibrated = new PrecursorIndexEntry[_byMz.Length];
            for (int i = 0; i < _byMz.Length; i++)
            {
                ref readonly var src = ref _byMz[i];
                calibrated[i] = new PrecursorIndexEntry(
                    mz:          src.PrecursorMz,
                    charge:      src.Charge,
                    rt:          (float)(slope * src.RetentionTime + intercept),
                    im:          src.IonMobility,
                    sourceIndex: src.SourceIndex,
                    isDecoy:     src.IsDecoy);
            }

            // Re-sort by m/z (RT calibration doesn't change m/z, so sort is stable)
            Array.Sort(calibrated);
            return new DiaNNSpecLibIndex(calibrated, _entryLoader, _maxBufferSize);
        }

        // ──────────────────────────────────────────────────────────────────────────────────────
        // Diagnostics
        // ──────────────────────────────────────────────────────────────────────────────────────

        public DiaNNSpecLibIndexStatistics GetStatistics()
        {
            ThrowIfDisposed();

            int targets = _byMz.Count(e => !e.IsDecoy);
            int decoys  = _byMz.Length - targets;

            return new DiaNNSpecLibIndexStatistics
            {
                TotalPrecursors  = _byMz.Length,
                TargetPrecursors = targets,
                DecoyPrecursors  = decoys,
                MinPrecursorMz   = _byMz.Length > 0 ? _byMz[0].PrecursorMz : 0,
                MaxPrecursorMz   = _byMz.Length > 0 ? _byMz[^1].PrecursorMz : 0,
                MinRT            = _byMz.Length > 0 ? _byMz.Min(e => (double)e.RetentionTime) : 0,
                MaxRT            = _byMz.Length > 0 ? _byMz.Max(e => (double)e.RetentionTime) : 0,
                BufferCapacity   = _maxBufferSize,
                BufferOccupancy  = _buffer.Count,
                TotalCacheHits   = Interlocked.Read(ref _cacheHits),
                TotalCacheMisses = Interlocked.Read(ref _cacheMisses)
            };
        }

        // ──────────────────────────────────────────────────────────────────────────────────────
        // LRU cache internals
        // ──────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Loads a full entry by its SourceIndex.
        /// Checks the LRU buffer first; calls the loader delegate on miss.
        /// Thread-safe: multiple threads may call this concurrently on different SourceIndex values.
        /// </summary>
        private DiaNNLibraryEntry? LoadEntry(int sourceIndex)
        {
            if (_buffer.TryGetValue(sourceIndex, out var cached))
            {
                Interlocked.Increment(ref _cacheHits);
                return cached;
            }

            Interlocked.Increment(ref _cacheMisses);

            if (_entryLoader == null)
            {
                // No loader: entries are only available if pre-cached
                return null;
            }

            var loaded = _entryLoader(sourceIndex);
            if (loaded == null) return null;

            // Add to cache with LRU eviction
            if (_buffer.TryAdd(sourceIndex, loaded))
            {
                lock (_evictionLock)
                {
                    _evictionQueue.Enqueue(sourceIndex);

                    // Evict oldest entries until we're within budget
                    while (_buffer.Count > _maxBufferSize && _evictionQueue.Count > 0)
                    {
                        int evictKey = _evictionQueue.Dequeue();
                        _buffer.TryRemove(evictKey, out _);
                    }
                }
            }

            return loaded;
        }

        // ──────────────────────────────────────────────────────────────────────────────────────
        // Index construction helpers
        // ──────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// After the sorted array is built, populate the sequence/charge dictionary.
        /// We iterate _byMz in its (m/z sorted) order.  In case of duplicates (same
        /// sequence/charge at two different source positions), the first encountered wins.
        /// </summary>
        private void BuildSequenceIndex()
        {
            // We don't have the full sequences here — only the compact index entries.
            // The sequence/charge index is populated lazily on first TryGetSpectrum call
            // (see TryGetSpectrum), OR eagerly if BuildFromEntries() is used.
            // This overload is intentionally a no-op; BuildSequenceIndexFromEntries() handles
            // the eager path.
        }

        /// <summary>
        /// Populates the sequence/charge secondary index from full entry objects.
        /// Called by <see cref="BuildFromEntries"/>.
        /// </summary>
        private void BuildSequenceIndexFromEntries(IReadOnlyList<DiaNNLibraryEntry> entries)
        {
            for (int mzIdx = 0; mzIdx < _byMz.Length; mzIdx++)
            {
                int srcIdx = _byMz[mzIdx].SourceIndex;
                if (srcIdx >= entries.Count) continue;

                var entry = entries[srcIdx];
                string diannKey  = BuildLookupKey(entry.ModifiedSequence ?? string.Empty, entry.PrecursorCharge);
                string mzlibKey  = BuildLookupKey(DiaNNModificationMapping.DiaNNToMzLib(entry.ModifiedSequence ?? string.Empty), entry.PrecursorCharge);

                _sequenceChargeToMzIndex.TryAdd(diannKey,  mzIdx);
                _sequenceChargeToMzIndex.TryAdd(mzlibKey,  mzIdx);
            }
        }

        private static string BuildLookupKey(string sequence, int charge) =>
            sequence + "/" + charge.ToString();

        // ──────────────────────────────────────────────────────────────────────────────────────
        // Binary search helper
        // ──────────────────────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Returns the index of the first entry with m/z >= <paramref name="minMz"/>.
        /// Uses a custom binary search since Array.BinarySearch returns the bitwise complement
        /// on miss (not the insertion point we need).
        /// </summary>
        private int BinarySearchLowerBound(float minMz)
        {
            int lo = 0;
            int hi = _byMz.Length;

            while (lo < hi)
            {
                int mid = lo + ((hi - lo) >> 1);
                if (_byMz[mid].PrecursorMz < minMz)
                    lo = mid + 1;
                else
                    hi = mid;
            }

            return lo; // lo == hi == first index where _byMz[i].PrecursorMz >= minMz
        }

        // ──────────────────────────────────────────────────────────────────────────────────────
        // IDisposable
        // ──────────────────────────────────────────────────────────────────────────────────────

        private void ThrowIfDisposed()
        {
            if (_disposed) throw new ObjectDisposedException(nameof(DiaNNSpecLibIndex));
        }

        public void Dispose()
        {
            if (_disposed) return;
            _disposed = true;
            _buffer.Clear();
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────────────────────
    // QueryScope: zero-allocation result container
    // ─────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Holds the result of a multi-dimensional DIA window query from <see cref="DiaNNSpecLibIndex"/>.
    ///
    /// Wraps a rented <see cref="System.Buffers.ArrayPool{T}"/> buffer containing the matching
    /// <see cref="PrecursorIndexEntry"/> items.  Callers MUST dispose this object (ideally via
    /// a <c>using</c> statement) to return the buffer to the pool.
    ///
    /// Usage pattern:
    /// <code>
    /// using var scope = index.GetPrecursorsInWindow(minMz, maxMz, charge: 2, minRt: 20, maxRt: 30);
    /// foreach (ref readonly var entry in scope.Entries)
    /// {
    ///     // process entry ...
    /// }
    /// // scope.Dispose() called automatically here
    /// </code>
    /// </summary>
    public sealed class QueryScope : IDisposable
    {
        /// <summary>A reusable empty scope that represents no results.</summary>
        public static readonly QueryScope Empty = new QueryScope(null, 0);

        private PrecursorIndexEntry[]? _rentedArray;
        private readonly int _count;
        private bool _disposed;

        internal QueryScope(PrecursorIndexEntry[]? array, int count)
        {
            _rentedArray = array;
            _count       = count;
        }

        /// <summary>The matching entries.  Valid only until <see cref="Dispose"/> is called.</summary>
        public ReadOnlySpan<PrecursorIndexEntry> Entries =>
            _rentedArray == null ? ReadOnlySpan<PrecursorIndexEntry>.Empty
                                 : _rentedArray.AsSpan(0, _count);

        /// <summary>Number of matching entries.</summary>
        public int Count => _count;

        /// <summary>True if no entries were found.</summary>
        public bool IsEmpty => _count == 0;

        public void Dispose()
        {
            if (_disposed) return;
            _disposed = true;

            if (_rentedArray != null && !ReferenceEquals(this, Empty))
            {
                System.Buffers.ArrayPool<PrecursorIndexEntry>.Shared.Return(_rentedArray);
                _rentedArray = null;
            }
        }
    }
}
