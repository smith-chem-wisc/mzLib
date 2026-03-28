using Omics.SpectralMatch.MslSpectralLibrary;

namespace Readers.SpectralLibrary;

// ─────────────────────────────────────────────────────────────────────────────
// MslProteoformIndexEntry
// ─────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Compact index entry for a proteoform library precursor.
/// Sorted by <see cref="NeutralMass"/> ascending in the
/// <see cref="MslProteoformIndex"/> primary array.
/// At ~28 bytes per entry, 100,000 proteoform entries occupy ~2.8 MB —
/// entirely in L3 cache.
/// </summary>
public readonly struct MslProteoformIndexEntry : IComparable<MslProteoformIndexEntry>
{
	/// <summary>
	/// Neutral monoisotopic mass in daltons, computed as
	/// (PrecursorMz × ChargeState) − (ChargeState × 1.007276).
	/// Primary sort key; used for O(log N) binary search on mass window queries.
	/// </summary>
	public readonly double NeutralMass;

	/// <summary>Precursor m/z as stored in the binary record (float32).</summary>
	public readonly float PrecursorMz;

	/// <summary>Precursor charge state.</summary>
	public readonly short Charge;

	/// <summary>Indexed retention time (iRT or calibrated RT).</summary>
	public readonly float Irt;

	/// <summary>True for decoy entries.</summary>
	public readonly bool IsDecoy;

	/// <summary>
	/// Zero-based ordinal position of this entry in the source entry list.
	/// Used to retrieve the full <see cref="MslLibraryEntry"/> on demand.
	/// </summary>
	public readonly int OrdinalIndex;

	public MslProteoformIndexEntry(
		double neutralMass, float precursorMz, short charge,
		float irt, bool isDecoy, int ordinalIndex)
	{
		NeutralMass = neutralMass;
		PrecursorMz = precursorMz;
		Charge = charge;
		Irt = irt;
		IsDecoy = isDecoy;
		OrdinalIndex = ordinalIndex;
	}

	public int CompareTo(MslProteoformIndexEntry other)
		=> NeutralMass.CompareTo(other.NeutralMass);
}

// ─────────────────────────────────────────────────────────────────────────────
// MslProteoformIndex
// ─────────────────────────────────────────────────────────────────────────────

/// <summary>
/// In-memory neutral-mass index for top-down proteoform library search.
///
/// <para>
/// <b>Primary index:</b> a <see cref="MslProteoformIndexEntry"/> array sorted by
/// <see cref="MslProteoformIndexEntry.NeutralMass"/> ascending. Binary search on
/// neutral mass gives O(log N) range start; sequential scan gives O(k) for the
/// matched range.
/// </para>
///
/// <para>
/// <b>Secondary index:</b> a <c>Dictionary&lt;string, int&gt;</c> keyed by
/// <c>FullSequence + "/" + ChargeState</c> for O(1) DDA-style lookup by sequence
/// and charge. The value is the ordinal index into the source entry list.
/// </para>
///
/// <para>
/// The index stores only metadata; full <see cref="MslLibraryEntry"/> objects
/// (including fragment ions) are retrieved on demand via the loader delegate
/// supplied at <see cref="Build"/> time. In full-load mode the delegate is a
/// direct array lookup; in index-only mode it triggers a disk read through the
/// <c>MslLibraryData</c> LRU cache.
/// </para>
///
/// <para>
/// <b>Neutral mass computation:</b>
/// <c>NeutralMass = (PrecursorMz × ChargeState) − (ChargeState × ProtonMass)</c>
/// where <c>ProtonMass = 1.007276 Da</c>. Computation is performed in
/// <see cref="double"/> precision to preserve the mass accuracy required for
/// top-down search.
/// </para>
///
/// <para>
/// <b>Thread safety:</b> all query methods are safe for concurrent reads after
/// construction. The index is immutable; <see cref="WithCalibratedRetentionTimes"/>
/// returns a new instance.
/// </para>
/// </summary>
public sealed class MslProteoformIndex
{
	// ── Constants ─────────────────────────────────────────────────────────────

	/// <summary>
	/// Proton mass in daltons used for neutral mass computation.
	/// Value: 1.007276 Da.
	/// </summary>
	public const double ProtonMass = 1.007276;

	// ── Private state ─────────────────────────────────────────────────────────

	/// <summary>
	/// Primary array: all proteoform index entries sorted by
	/// <see cref="MslProteoformIndexEntry.NeutralMass"/> ascending.
	/// Backed by a plain array for zero-allocation span slicing.
	/// </summary>
	private readonly MslProteoformIndexEntry[] _entries;

	/// <summary>
	/// Secondary lookup: <c>"{FullSequence}/{ChargeState}"</c> → ordinal index in
	/// the source entry list. Ordinal is passed to <see cref="_loader"/> to retrieve
	/// the full <see cref="MslLibraryEntry"/>.
	/// </summary>
	private readonly Dictionary<string, int> _sequenceChargeToOrdinal;

	/// <summary>
	/// Delegate that loads a full <see cref="MslLibraryEntry"/> by ordinal index.
	/// In full-load mode this is a direct array access; in index-only mode it
	/// invokes the <c>MslLibraryData.LoadFragmentsOnDemand</c> path.
	/// May be null when the library was built without a loader (then
	/// <see cref="GetEntry(MslProteoformIndexEntry)"/> cannot return fragment data).
	/// </summary>
	private readonly Func<int, MslLibraryEntry?>? _loader;

	// ── Private constructor ───────────────────────────────────────────────────

	private MslProteoformIndex(
		MslProteoformIndexEntry[] entries,
		Dictionary<string, int> sequenceChargeToOrdinal,
		Func<int, MslLibraryEntry?>? loader)
	{
		_entries = entries;
		_sequenceChargeToOrdinal = sequenceChargeToOrdinal;
		_loader = loader;
	}

	// ── Construction ──────────────────────────────────────────────────────────

	/// <summary>
	/// Builds an <see cref="MslProteoformIndex"/> from a list of library entries.
	/// Only entries with <see cref="MslFormat.MoleculeType.Proteoform"/> are indexed;
	/// all other molecule types are silently skipped.
	/// </summary>
	/// <param name="entries">
	///   Source entry list. Must not be null. Non-proteoform entries are ignored.
	/// </param>
	/// <param name="loader">
	///   Delegate that loads a full <see cref="MslLibraryEntry"/> by ordinal index.
	///   In full-load mode pass <c>i => entries[i]</c>; in index-only mode pass the
	///   demand-loader closure. May be null when only index-level queries are needed.
	/// </param>
	/// <returns>
	///   A fully constructed <see cref="MslProteoformIndex"/> sorted by neutral mass.
	///   Returns an empty index when no proteoform entries are present.
	/// </returns>
	public static MslProteoformIndex Build(
		IReadOnlyList<MslLibraryEntry> entries,
		Func<int, MslLibraryEntry?>? loader)
	{
		ArgumentNullException.ThrowIfNull(entries);

		// Collect proteoform entries and build both indexes in a single pass.
		var rawList = new List<MslProteoformIndexEntry>(capacity: entries.Count / 4);
		var seqDict = new Dictionary<string, int>(StringComparer.Ordinal);

		for (int i = 0; i < entries.Count; i++)
		{
			MslLibraryEntry e = entries[i];
			if (e.MoleculeType != MslFormat.MoleculeType.Proteoform)
				continue;

			double neutralMass = ComputeNeutralMass((double)e.PrecursorMz, e.ChargeState);

			var indexEntry = new MslProteoformIndexEntry(
				neutralMass: neutralMass,
				precursorMz: (float)e.PrecursorMz,
				charge: (short)e.ChargeState,
				irt: (float)e.RetentionTime,
				isDecoy: e.IsDecoy,
				ordinalIndex: i);

			rawList.Add(indexEntry);

			// Secondary index key: "Sequence/ChargeState"
			string key = e.FullSequence + "/" + e.ChargeState;
			seqDict.TryAdd(key, i);
		}

		// Sort primary array by neutral mass ascending for binary search.
		MslProteoformIndexEntry[] sorted = rawList.ToArray();
		Array.Sort(sorted);

		return new MslProteoformIndex(sorted, seqDict, loader);
	}

	// ── Neutral mass helper ───────────────────────────────────────────────────

	/// <summary>
	/// Computes the neutral monoisotopic mass of a precursor from its observed m/z
	/// and charge state.
	/// <para>
	/// Formula: <c>NeutralMass = (mz × charge) − (charge × ProtonMass)</c>
	/// </para>
	/// Computation is performed entirely in <see cref="double"/> precision; the
	/// <paramref name="mz"/> input is widened from float32 storage to double
	/// before arithmetic to avoid accumulating float rounding error across the
	/// charge multiplication.
	/// </summary>
	public static double ComputeNeutralMass(double mz, int charge)
		=> (mz * charge) - (charge * ProtonMass);

	// ── Properties ────────────────────────────────────────────────────────────

	/// <summary>Number of proteoform entries in the index.</summary>
	public int Count => _entries.Length;

	/// <summary>
	/// Minimum neutral mass in the index (Da).
	/// Returns 0 when the index is empty.
	/// </summary>
	public double MinNeutralMass
		=> _entries.Length > 0 ? _entries[0].NeutralMass : 0.0;

	/// <summary>
	/// Maximum neutral mass in the index (Da).
	/// Returns 0 when the index is empty.
	/// </summary>
	public double MaxNeutralMass
		=> _entries.Length > 0 ? _entries[_entries.Length - 1].NeutralMass : 0.0;

	// ── Primary query: neutral mass window ────────────────────────────────────

	/// <summary>
	/// Returns all index entries whose neutral mass falls within
	/// [<paramref name="minMass"/>, <paramref name="maxMass"/>].
	///
	/// <para>
	/// <b>Zero-allocation.</b> The returned span is a slice of the internal sorted
	/// array; no heap allocation occurs. The span is valid for the lifetime of this
	/// <see cref="MslProteoformIndex"/> instance.
	/// </para>
	///
	/// <para>
	/// Implementation uses <see cref="Array.BinarySearch{T}"/> twice — once to find
	/// the inclusive lower bound and once to find the inclusive upper bound — giving
	/// O(log N) complexity. The sequential scan over the matched range is O(k) where
	/// k is the number of results.
	/// </para>
	/// </summary>
	/// <param name="minMass">Lower neutral mass bound in daltons (inclusive).</param>
	/// <param name="maxMass">Upper neutral mass bound in daltons (inclusive).</param>
	/// <param name="includeDecoys">
	///   When <see langword="true"/> (default), decoy entries are included.
	///   When <see langword="false"/>, decoy entries are excluded. Note that
	///   excluding decoys allocates a filtered copy; use <see langword="true"/>
	///   when allocation must be avoided on the hot path.
	/// </param>
	/// <returns>
	///   A <see cref="ReadOnlySpan{T}"/> of matching entries in ascending neutral
	///   mass order. Empty span when no entries fall in the requested window.
	/// </returns>
	public ReadOnlySpan<MslProteoformIndexEntry> QueryMassWindow(
		double minMass,
		double maxMass,
		bool includeDecoys = true)
	{
		if (_entries.Length == 0 || minMass > maxMass)
			return ReadOnlySpan<MslProteoformIndexEntry>.Empty;

		// Find the first index whose NeutralMass >= minMass using binary search.
		// We search with a sentinel entry that has NeutralMass == minMass and all
		// other fields at their minimum so the comparer places it at the correct
		// lower bound.
		int lo = LowerBound(minMass);
		if (lo >= _entries.Length || _entries[lo].NeutralMass > maxMass)
			return ReadOnlySpan<MslProteoformIndexEntry>.Empty;

		// Find the last index whose NeutralMass <= maxMass.
		int hi = UpperBound(maxMass);

		// hi is the index of the last entry ≤ maxMass; length = hi - lo + 1.
		int length = hi - lo + 1;
		if (length <= 0)
			return ReadOnlySpan<MslProteoformIndexEntry>.Empty;

		ReadOnlySpan<MslProteoformIndexEntry> slice =
			new ReadOnlySpan<MslProteoformIndexEntry>(_entries, lo, length);

		if (includeDecoys)
			return slice;

		// Decoy filtering requires a copy because we cannot return a non-contiguous
		// subset of the backing array as a zero-allocation span.
		MslProteoformIndexEntry[] filtered =
			slice.ToArray().Where(e => !e.IsDecoy).ToArray();
		return new ReadOnlySpan<MslProteoformIndexEntry>(filtered);
	}

	// ── Secondary query: sequence + charge lookup ─────────────────────────────

	/// <summary>
	/// O(1) lookup by modified sequence and charge state.
	/// Populates <paramref name="entry"/> with fragment ion data on hit.
	/// </summary>
	/// <param name="modifiedSequence">
	///   Modified sequence in mzLib bracket notation. Case-sensitive.
	/// </param>
	/// <param name="charge">Precursor charge state.</param>
	/// <param name="entry">
	///   On success, the full <see cref="MslLibraryEntry"/> with fragments populated;
	///   otherwise <see langword="null"/>.
	/// </param>
	/// <returns>
	///   <see langword="true"/> when found; <see langword="false"/> otherwise.
	/// </returns>
	public bool TryGetEntry(string modifiedSequence, int charge, out MslLibraryEntry? entry)
	{
		string key = modifiedSequence + "/" + charge;

		if (!_sequenceChargeToOrdinal.TryGetValue(key, out int ordinal))
		{
			entry = null;
			return false;
		}

		entry = _loader?.Invoke(ordinal);
		return entry is not null;
	}

	// ── Full entry retrieval ──────────────────────────────────────────────────

	/// <summary>
	/// Retrieves the full <see cref="MslLibraryEntry"/> for a given index entry,
	/// loading fragment data from disk if in index-only mode.
	/// Returns <see langword="null"/> when no loader was supplied at build time or
	/// when the ordinal is out of range.
	/// </summary>
	public MslLibraryEntry? GetEntry(MslProteoformIndexEntry indexEntry)
		=> _loader?.Invoke(indexEntry.OrdinalIndex);

	/// <summary>
	/// Enumerates all proteoform entries in ascending neutral mass order.
	/// Fragment data is loaded on demand for each entry.
	/// Decoy entries are included when <paramref name="includeDecoys"/> is
	/// <see langword="true"/> (default).
	/// </summary>
	public IEnumerable<MslLibraryEntry> GetAllEntries(bool includeDecoys = true)
	{
		foreach (MslProteoformIndexEntry e in _entries)
		{
			if (!includeDecoys && e.IsDecoy)
				continue;

			MslLibraryEntry? full = _loader?.Invoke(e.OrdinalIndex);
			if (full is not null)
				yield return full;
		}
	}

	// ── RT calibration ────────────────────────────────────────────────────────

	/// <summary>
	/// Returns a new <see cref="MslProteoformIndex"/> with all iRT values transformed
	/// by the linear function: <c>RT_calibrated = slope × RT_stored + intercept</c>.
	/// The original index is not modified.
	/// </summary>
	/// <param name="slope">
	///   Slope of the linear iRT → calibrated RT regression. Typical range: 0.1–10.
	/// </param>
	/// <param name="intercept">
	///   Intercept in the same units as the calibrated RT axis.
	/// </param>
	/// <returns>A new <see cref="MslProteoformIndex"/> with transformed RetentionTime values.</returns>
	public MslProteoformIndex WithCalibratedRetentionTimes(double slope, double intercept)
	{
		MslProteoformIndexEntry[] calibrated = new MslProteoformIndexEntry[_entries.Length];

		for (int i = 0; i < _entries.Length; i++)
		{
			MslProteoformIndexEntry src = _entries[i];
			float newIrt = (float)(slope * src.Irt + intercept);
			calibrated[i] = new MslProteoformIndexEntry(
				src.NeutralMass, src.PrecursorMz, src.Charge,
				newIrt, src.IsDecoy, src.OrdinalIndex);
		}

		// The calibrated array is already sorted by NeutralMass (RT calibration does
		// not affect neutral mass ordering), so no re-sort is needed.
		return new MslProteoformIndex(calibrated, _sequenceChargeToOrdinal, _loader);
	}

	// ── Binary search helpers ─────────────────────────────────────────────────

	/// <summary>
	/// Returns the index of the first element whose <c>NeutralMass >= minMass</c>,
	/// or <c>_entries.Length</c> if all elements are below <paramref name="minMass"/>.
	/// </summary>
	private int LowerBound(double minMass)
	{
		int lo = 0, hi = _entries.Length;

		while (lo < hi)
		{
			int mid = lo + ((hi - lo) >> 1);
			if (_entries[mid].NeutralMass < minMass)
				lo = mid + 1;
			else
				hi = mid;
		}

		return lo;
	}

	/// <summary>
	/// Returns the index of the last element whose <c>NeutralMass <= maxMass</c>,
	/// or -1 if all elements exceed <paramref name="maxMass"/>.
	/// </summary>
	private int UpperBound(double maxMass)
	{
		int lo = 0, hi = _entries.Length - 1, result = -1;

		while (lo <= hi)
		{
			int mid = lo + ((hi - lo) >> 1);
			if (_entries[mid].NeutralMass <= maxMass)
			{
				result = mid;
				lo = mid + 1;
			}
			else
			{
				hi = mid - 1;
			}
		}

		return result;
	}
}