using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;

namespace Readers.SpectralLibrary;

/// <summary>
/// Top-level user-facing façade for the .msl binary spectral library format.
///
/// <para>
/// <b>Two load modes are supported:</b>
/// <list type="bullet">
///   <item>
///     <term>Full load (<see cref="Load"/>)</term>
///     <description>
///       All precursor entries and fragment ions are read into memory. No file handle is held
///       open after loading. Best for libraries that fit comfortably in RAM. <see cref="Dispose"/>
///       is a no-op but safe to call.
///     </description>
///   </item>
///   <item>
///     <term>Index-only load (<see cref="LoadIndexOnly"/>)</term>
///     <description>
///       Only precursor metadata is loaded into memory; fragment blocks remain on disk and are
///       fetched on demand through an open <see cref="System.IO.FileStream"/>. Best for
///       multi-gigabyte libraries where loading all fragments would exhaust available RAM.
///       <b>Callers must call <see cref="Dispose"/> when finished</b> to release the file handle
///       (mandatory on Windows due to OS-level file locking).
///     </description>
///   </item>
/// </list>
/// </para>
///
/// <para>
/// All query operations produce identical results in both modes. The only observable difference
/// is that <see cref="IsIndexOnly"/> returns <see langword="true"/> in index-only mode and that
/// fragment ions for a given entry may be populated lazily rather than eagerly.
/// </para>
///
/// <para>
/// <b>Thread safety:</b> all read-only query methods (<see cref="QueryMzWindow"/>,
/// <see cref="TryGetEntry"/>, <see cref="GetAllEntries"/>, etc.) are safe to call from multiple
/// threads simultaneously. <see cref="Dispose"/> must not be called concurrently with any query
/// operation.
/// </para>
/// </summary>
public sealed class MslLibrary : IDisposable
{
	// ── Private state ─────────────────────────────────────────────────────────

	/// <summary>
	/// The in-memory query engine. Holds the m/z-sorted <see cref="MslPrecursorIndexEntry"/>
	/// array, the sequence/charge dictionary, the elution-group map, and the LRU entry cache.
	/// Non-null for the lifetime of this instance; nulled on <see cref="Dispose"/>.
	/// </summary>
	private MslIndex? _index;

	/// <summary>
	/// The Prompt 3 data container that owns the open <see cref="System.IO.FileStream"/> and
	/// the <see cref="MslLibraryData.LoadFragmentsOnDemand"/> implementation in index-only mode.
	/// Null in full-load mode and after <see cref="Dispose"/> has been called.
	/// </summary>
	private readonly MslLibraryData? _rawLibrary;

	/// <summary>
	/// Pre-computed statistics snapshot from the index. Cached once at construction so that
	/// property getters (<see cref="PrecursorCount"/>, <see cref="MinPrecursorMz"/>, etc.)
	/// execute in O(1) without re-querying the index.
	/// </summary>
	private readonly MslIndexStatistics _stats;

	/// <summary>
	/// Guards the disposed state so that concurrent callers see a consistent value and do not
	/// race with <see cref="Dispose"/>.
	/// </summary>
	private volatile bool _disposed;

	// ── Private constructor ───────────────────────────────────────────────────

	/// <summary>
	/// Internal constructor used by <see cref="Load"/>, <see cref="LoadIndexOnly"/>, and
	/// <see cref="WithCalibratedRetentionTimes"/>. All public entry points go through one of
	/// the static factory methods; direct construction is not exposed.
	/// </summary>
	/// <param name="index">
	///   Fully constructed <see cref="MslIndex"/> for this library. Ownership is transferred;
	///   will be disposed by <see cref="Dispose"/>.
	/// </param>
	/// <param name="header">
	///   Deserialized file header containing format version, flags, and section offsets.
	/// </param>
	/// <param name="isIndexOnly">
	///   <see langword="true"/> when the library was opened in index-only mode;
	///   <see langword="false"/> for full-load mode.
	/// </param>
	/// <param name="rawLibrary">
	///   The <see cref="MslLibraryData"/> data container holding the
	///   <see cref="System.IO.FileStream"/> and on-demand fragment loader in index-only mode.
	///   <see langword="null"/> in full-load mode.
	/// </param>
	private MslLibrary(
		MslIndex index,
		MslFileHeader header,
		bool isIndexOnly,
		MslLibraryData? rawLibrary)
	{
		_index = index;
		Header = header;
		IsIndexOnly = isIndexOnly;
		_rawLibrary = rawLibrary;
		_stats = index.GetStatistics();
	}

	// ── Static factory methods ────────────────────────────────────────────────

	/// <summary>
	/// Opens a .msl file and reads all precursor entries and fragment ions into memory.
	///
	/// <para>
	/// No file handle is retained after this method returns. The returned library is safe to use
	/// without calling <see cref="Dispose"/>, though disposing is always harmless.
	/// </para>
	///
	/// <para>
	/// Best for libraries that fit comfortably in RAM (typically &lt;~2 GB). For larger
	/// libraries use <see cref="LoadIndexOnly"/> instead.
	/// </para>
	/// </summary>
	/// <param name="filePath">
	///   Absolute or relative path to the .msl file. Must not be null or empty. The file must
	///   pass all five validation checks (magic, version, footer magic, precursor count, CRC-32).
	/// </param>
	/// <returns>
	///   A fully populated <see cref="MslLibrary"/> in full-load mode where
	///   <see cref="IsIndexOnly"/> is <see langword="false"/> and every entry's
	///   <c>Fragments</c> list is populated.
	/// </returns>
	/// <exception cref="ArgumentNullException"><paramref name="filePath"/> is null.</exception>
	/// <exception cref="FileNotFoundException">No file exists at <paramref name="filePath"/>.</exception>
	/// <exception cref="FormatException">Magic, version, or structural validation failed.</exception>
	/// <exception cref="InvalidDataException">CRC-32 checksum mismatch (data corruption).</exception>
	public static MslLibrary Load(string filePath)
	{
		ArgumentNullException.ThrowIfNull(filePath);

		// Delegate full deserialization to MslReader; receives an MslLibraryData container
		MslLibraryData rawLib = MslReader.Load(filePath);

		// Build the query index from the fully-loaded entries; the loader delegate is a direct
		// array index lookup because all fragments are already in memory.
		MslIndex index = MslIndex.Build(
			rawLib.Entries,
			i => i >= 0 && i < rawLib.Count ? rawLib.Entries[i] : null);

		return new MslLibrary(index, rawLib.Header, isIndexOnly: false, rawLibrary: null);
	}

	/// <summary>
	/// Opens a .msl file in index-only mode: precursor metadata is loaded into memory but
	/// fragment blocks remain on disk and are read lazily on demand.
	///
	/// <para>
	/// The returned library holds an open <see cref="System.IO.FileStream"/>.
	/// <b>Callers must call <see cref="Dispose"/> when finished</b> to release the file handle.
	/// On Windows, failing to dispose prevents the file from being moved or deleted while the
	/// library is in scope.
	/// </para>
	///
	/// <para>
	/// All query APIs produce identical results to those returned by <see cref="Load"/>. The
	/// difference is that <see cref="GetEntry"/> and <see cref="GetAllEntries"/> trigger disk
	/// reads for fragment data rather than returning pre-populated objects.
	/// </para>
	///
	/// <para>
	/// <b>Compressed files:</b> when <see cref="MslFormat.FileFlagIsCompressed"/> is set in the
	/// file header, index-only mode is not available because fragment block offsets are relative
	/// to the decompressed buffer rather than absolute file positions. This method transparently
	/// falls back to full decompression in that case and returns a library with
	/// <see cref="IsIndexOnly"/> == <see langword="false"/> and
	/// <see cref="IsCompressed"/> == <see langword="true"/>. No exception is thrown.
	/// </para>
	/// </summary>
	/// <param name="filePath">
	///   Absolute or relative path to the .msl file. Must not be null or empty.
	/// </param>
	/// <returns>
	///   An <see cref="MslLibrary"/> in index-only mode for uncompressed files, or in full-load
	///   mode for compressed files (see compressed-file note above).
	/// </returns>
	/// <exception cref="ArgumentNullException"><paramref name="filePath"/> is null.</exception>
	/// <exception cref="FileNotFoundException">No file exists at <paramref name="filePath"/>.</exception>
	/// <exception cref="FormatException">Magic, version, or structural validation failed.</exception>
	/// <exception cref="InvalidDataException">CRC-32 checksum mismatch (data corruption).</exception>
	public static MslLibrary LoadIndexOnly(string filePath)
	{
		ArgumentNullException.ThrowIfNull(filePath);

		// Delegate to MslReader. For uncompressed files this returns a true index-only
		// MslLibraryData (IsIndexOnly = true) with an open FileStream for on-demand reads.
		// For compressed files, MslReader.LoadIndexOnly falls back to full-load and returns
		// an MslLibraryData with IsIndexOnly = false and all fragments already in memory.
		MslLibraryData rawLib = MslReader.LoadIndexOnly(filePath);

		MslIndex index;

		if (!rawLib.IsIndexOnly)
		{
			// Compressed-file fallback: rawLib is fully loaded — use direct array access,
			// same as the MslLibrary.Load path. Do not call LoadFragmentsOnDemand.
			index = MslIndex.Build(
				rawLib.Entries,
				i => i >= 0 && i < rawLib.Count ? rawLib.Entries[i] : null);

			// rawLib has no open stream for compressed files; no rawLibrary to retain.
			return new MslLibrary(index, rawLib.Header, isIndexOnly: false, rawLibrary: null);
		}

		// Uncompressed index-only path: capture rawLib in the loader closure.
		// On index cache miss the delegate seeks the open FileStream and reads the fragment block.
		index = MslIndex.Build(
			rawLib.Entries,
			i =>
			{
				if (i < 0 || i >= rawLib.Count)
					return null;
				MslLibraryEntry skeleton = rawLib.Entries[i];
				skeleton.Fragments = rawLib.LoadFragmentsOnDemand(i);
				return skeleton;
			});

		return new MslLibrary(index, rawLib.Header, isIndexOnly: true, rawLibrary: rawLib);
	}

	/// <summary>
	/// Serializes a list of <see cref="MslLibraryEntry"/> objects to a .msl binary file.
	///
	/// <para>
	/// The write is atomic: bytes are first written to a temporary file and then the temporary
	/// file is renamed to <paramref name="filePath"/>, so a failed write never leaves a
	/// partially-written file at the destination path.
	/// </para>
	///
	/// <para>
	/// Fragment intensities are normalized per precursor (max = 1.0) and fragment records within
	/// each precursor are sorted by m/z ascending before writing.
	/// </para>
	/// </summary>
	/// <param name="filePath">
	///   Destination file path. Created or overwritten. The parent directory must already exist.
	/// </param>
	/// <param name="entries">
	///   Ordered list of library entries to write. Must not be null. May be empty, which produces
	///   a valid zero-precursor .msl file.
	/// </param>
	/// <exception cref="ArgumentNullException">
	///   <paramref name="filePath"/> or <paramref name="entries"/> is null.
	/// </exception>
	/// <exception cref="IOException">An I/O error occurred while writing the file.</exception>
	public static void Save(string filePath, IReadOnlyList<MslLibraryEntry> entries)
	{
		ArgumentNullException.ThrowIfNull(filePath);
		ArgumentNullException.ThrowIfNull(entries);

		MslWriter.Write(filePath, entries);
	}

	/// <summary>
	/// Converts a list of <see cref="LibrarySpectrum"/> objects to <see cref="MslLibraryEntry"/>
	/// instances and serializes them to a .msl binary file.
	///
	/// <para>
	/// Extended metadata fields that have no representation in <see cref="LibrarySpectrum"/>
	/// (ion mobility, protein accession, q-value, etc.) are initialized to sensible defaults
	/// by <see cref="MslLibraryEntry.FromLibrarySpectrum"/>.
	/// </para>
	///
	/// <para>The write is atomic; see <see cref="Save"/> for details.</para>
	/// </summary>
	/// <param name="filePath">
	///   Destination file path. Created or overwritten. The parent directory must already exist.
	/// </param>
	/// <param name="spectra">
	///   Library spectra to convert and write. Must not be null. May be empty.
	/// </param>
	/// <exception cref="ArgumentNullException">
	///   <paramref name="filePath"/> or <paramref name="spectra"/> is null.
	/// </exception>
	/// <exception cref="IOException">An I/O error occurred while writing the file.</exception>
	public static void SaveFromLibrarySpectra(string filePath, IReadOnlyList<LibrarySpectrum> spectra)
	{
		ArgumentNullException.ThrowIfNull(filePath);
		ArgumentNullException.ThrowIfNull(spectra);

		// Convert each LibrarySpectrum to an MslLibraryEntry, then delegate to the primary Save path
		var entries = new List<MslLibraryEntry>(spectra.Count);
		foreach (LibrarySpectrum spectrum in spectra)
			entries.Add(MslLibraryEntry.FromLibrarySpectrum(spectrum));

		MslWriter.Write(filePath, entries);
	}

	// ── Public properties ─────────────────────────────────────────────────────

	/// <summary>
	/// Total number of precursor entries in the library, including both targets and decoys.
	/// Equivalent to <c>TargetCount + DecoyCount</c>. Available in O(1) in both load modes
	/// because it is derived from the pre-computed statistics snapshot.
	/// </summary>
	public int PrecursorCount => _stats.TotalPrecursors;

	/// <summary>
	/// Number of non-decoy (target) precursor entries in the library.
	/// Satisfies the invariant <c>TargetCount + DecoyCount == PrecursorCount</c>.
	/// </summary>
	public int TargetCount => _stats.TargetPrecursors;

	/// <summary>
	/// Number of decoy precursor entries in the library.
	/// Satisfies the invariant <c>TargetCount + DecoyCount == PrecursorCount</c>.
	/// </summary>
	public int DecoyCount => _stats.DecoyPrecursors;

	/// <summary>
	/// Smallest precursor m/z value among all entries, in Thomson (m/z) units.
	/// Returns 0 when the library contains no entries.
	/// </summary>
	public float MinPrecursorMz => _stats.MinPrecursorMz;

	/// <summary>
	/// Largest precursor m/z value among all entries, in Thomson (m/z) units.
	/// Returns 0 when the library contains no entries.
	/// </summary>
	public float MaxPrecursorMz => _stats.MaxPrecursorMz;

	/// <summary>
	/// Deserialized file header providing format version, file-level flags, and byte offsets
	/// of each binary section within the .msl file. Immutable after construction.
	/// </summary>
	public MslFileHeader Header { get; }

	/// <summary>
	/// <see langword="true"/> when this instance was opened via <see cref="LoadIndexOnly"/>;
	/// <see langword="false"/> when opened via <see cref="Load"/>.
	///
	/// In index-only mode, fragment ions for each precursor are loaded from disk on demand
	/// rather than pre-populated in memory. Callers must <see cref="Dispose"/> index-only
	/// libraries to release the underlying file handle.
	/// </summary>
	public bool IsIndexOnly { get; }

	/// <summary>
	/// True when the library file was written with zstd block compression
	/// (<see cref="MslFormat.FileFlagIsCompressed"/> is set in <see cref="Header"/>).
	/// Compressed libraries always load in full-load mode regardless of which
	/// Load method was called; <see cref="IsIndexOnly"/> is always <see langword="false"/>
	/// for compressed files.
	/// </summary>
	public bool IsCompressed => (Header.FileFlags & MslFormat.FileFlagIsCompressed) != 0;

	// ── Disposal guard ────────────────────────────────────────────────────────

	/// <summary>
	/// Throws <see cref="ObjectDisposedException"/> when called after <see cref="Dispose"/>.
	/// Every public method that touches <see cref="_index"/> must call this first.
	/// </summary>
	/// <exception cref="ObjectDisposedException">This instance has already been disposed.</exception>
	private void ThrowIfDisposed()
	{
		if (_disposed)
			throw new ObjectDisposedException(nameof(MslLibrary),
				"This MslLibrary instance has been disposed and can no longer be used.");
	}

	// ── DDA-style lookup ──────────────────────────────────────────────────────

	/// <summary>
	/// Attempts to find the precursor whose modified sequence and charge state match the
	/// given arguments and returns the corresponding <see cref="MslLibraryEntry"/>.
	///
	/// <para>
	/// In full-load mode the returned entry already contains fully populated fragment ions.
	/// In index-only mode fragment ions are loaded from disk on demand via the LRU cache inside
	/// <see cref="MslIndex"/>; the returned entry is therefore fully populated in both modes.
	/// </para>
	///
	/// <para>
	/// The lookup is O(1): a pre-built dictionary keyed on
	/// <c>"{modifiedSequence}/{charge}"</c> is consulted inside <see cref="MslIndex"/>.
	/// Key comparison is case-sensitive ordinal.
	/// </para>
	/// </summary>
	/// <param name="modifiedSequence">
	///   The modified sequence in mzLib bracket notation, e.g.
	///   <c>"PEPTM[Common Variable:Oxidation on M]IDE"</c>. Must not be null.
	/// </param>
	/// <param name="charge">Precursor charge state (typically 1–6).</param>
	/// <param name="entry">
	///   When this method returns <see langword="true"/>, contains the matching
	///   <see cref="MslLibraryEntry"/> with fragment ions populated; otherwise
	///   <see langword="null"/>.
	/// </param>
	/// <returns>
	///   <see langword="true"/> if a matching precursor was found; <see langword="false"/>
	///   otherwise.
	/// </returns>
	/// <exception cref="ObjectDisposedException">This instance has been disposed.</exception>
	public bool TryGetEntry(string modifiedSequence, int charge, out MslLibraryEntry? entry)
	{
		ThrowIfDisposed();

		if (!_index!.TryGetBySequenceCharge(modifiedSequence, charge, out MslPrecursorIndexEntry indexEntry))
		{
			entry = null;
			return false;
		}

		// Delegate to the index's LRU-buffered loader, which handles both full-load
		// (direct array access) and index-only (on-demand disk read) transparently.
		entry = _index.GetEntry(indexEntry.PrecursorIdx);
		return entry is not null;
	}

	/// <summary>
	/// Attempts to find a precursor by modified sequence and charge state and returns it as a
	/// <see cref="LibrarySpectrum"/> compatible with the MetaMorpheus spectral-angle scoring
	/// pipeline.
	///
	/// <para>
	/// This method has the same signature as
	/// <see cref="SpectralLibrary.TryGetSpectrum(string, int, out LibrarySpectrum)"/> and is
	/// the primary interop point between <see cref="MslLibrary"/> and MetaMorpheus. Callers
	/// that currently accept a <see cref="SpectralLibrary"/> can be adapted to accept an
	/// <see cref="MslLibrary"/> by replacing the <c>TryGetSpectrum</c> call with this method.
	/// </para>
	///
	/// <para>
	/// The conversion from <see cref="MslLibraryEntry"/> to <see cref="LibrarySpectrum"/> is
	/// performed by <see cref="MslLibraryEntry.ToLibrarySpectrum"/>. Internal-fragment metadata
	/// (<c>SecondaryProductType</c>, <c>SecondaryFragmentNumber</c>) is not round-trippable
	/// through <see cref="LibrarySpectrum"/> and is therefore dropped during conversion.
	/// </para>
	/// </summary>
	/// <param name="sequence">
	///   Modified sequence in mzLib bracket notation. Must not be null. Case-sensitive.
	/// </param>
	/// <param name="charge">Precursor charge state.</param>
	/// <param name="librarySpectrum">
	///   When this method returns <see langword="true"/>, contains the matching spectrum;
	///   otherwise <see langword="null"/>.
	/// </param>
	/// <returns>
	///   <see langword="true"/> if the precursor was found and successfully converted;
	///   <see langword="false"/> otherwise.
	/// </returns>
	/// <exception cref="ObjectDisposedException">This instance has been disposed.</exception>
	public bool TryGetLibrarySpectrum(string sequence, int charge, out LibrarySpectrum? librarySpectrum)
	{
		ThrowIfDisposed();

		if (!TryGetEntry(sequence, charge, out MslLibraryEntry? entry) || entry is null)
		{
			librarySpectrum = null;
			return false;
		}

		librarySpectrum = entry.ToLibrarySpectrum();
		return true;
	}

	// ── DIA window queries ────────────────────────────────────────────────────

	/// <summary>
	/// Returns a span of all <see cref="MslPrecursorIndexEntry"/> values whose
	/// <c>PrecursorMz</c> falls within the inclusive range [<paramref name="mzLow"/>,
	/// <paramref name="mzHigh"/>].
	///
	/// <para>
	/// <b>Zero-allocation.</b> The span is a slice of the internal sorted array inside
	/// <see cref="MslIndex"/>; no heap allocation occurs on any call path.
	/// The span is valid until this library instance is disposed.
	/// </para>
	///
	/// <para>
	/// This is the hot-path entry point for DIA isolation window queries. To load the
	/// full entry for a candidate, call <see cref="GetEntry"/> with its
	/// <see cref="MslPrecursorIndexEntry.PrecursorIdx"/>.
	/// </para>
	/// </summary>
	/// <param name="mzLow">Inclusive lower bound of the m/z window (Thomson).</param>
	/// <param name="mzHigh">Inclusive upper bound of the m/z window (Thomson).</param>
	/// <returns>
	///   A <see cref="ReadOnlySpan{T}"/> of matching index entries in ascending
	///   <c>PrecursorMz</c> order. Empty span when no entries fall in the window.
	/// </returns>
	/// <exception cref="ObjectDisposedException">This instance has been disposed.</exception>
	public ReadOnlySpan<MslPrecursorIndexEntry> QueryMzWindow(float mzLow, float mzHigh)
	{
		ThrowIfDisposed();
		return _index!.QueryMzRange(mzLow, mzHigh);
	}

	/// <summary>
	/// Returns all precursor index entries that pass m/z, retention-time, and optionally
	/// ion-mobility window filters.
	///
	/// <para>
	/// The result buffer is rented from <see cref="System.Buffers.ArrayPool{T}"/>.
	/// <b>The caller must call <see cref="MslWindowResults.Dispose"/> on the returned value</b>
	/// (typically via a <c>using</c> statement) to return the buffer to the pool. Failure to
	/// dispose does not corrupt data but will defeat the purpose of pooling and generate GC
	/// pressure.
	/// </para>
	///
	/// <para>
	/// The ion-mobility filter is skipped entirely when both <paramref name="imLow"/> and
	/// <paramref name="imHigh"/> are 0.
	/// </para>
	/// </summary>
	/// <param name="mzLow">Inclusive lower bound of the precursor m/z window.</param>
	/// <param name="mzHigh">Inclusive upper bound of the precursor m/z window.</param>
	/// <param name="rtLow">Inclusive lower bound of the retention-time (iRT) window.</param>
	/// <param name="rtHigh">Inclusive upper bound of the retention-time (iRT) window.</param>
	/// <param name="imLow">
	///   Inclusive lower bound of the ion-mobility window. Pass 0 (along with
	///   <paramref name="imHigh"/> == 0) to skip the ion-mobility filter entirely.
	/// </param>
	/// <param name="imHigh">
	///   Inclusive upper bound of the ion-mobility window. Pass 0 (along with
	///   <paramref name="imLow"/> == 0) to skip the ion-mobility filter entirely.
	/// </param>
	/// <param name="includeDecoys">
	///   When <see langword="true"/>, decoy precursors are included in the result.
	///   When <see langword="false"/> (default), only target precursors are returned.
	/// </param>
	/// <returns>
	///   An <see cref="MslWindowResults"/> containing all matching entries in ascending
	///   <c>PrecursorMz</c> order. The caller must dispose this value.
	/// </returns>
	/// <exception cref="ObjectDisposedException">This instance has been disposed.</exception>
	public MslWindowResults QueryWindow(
		float mzLow, float mzHigh,
		float rtLow, float rtHigh,
		float imLow = 0f, float imHigh = 0f,
		bool includeDecoys = false)
	{
		ThrowIfDisposed();
		return _index!.QueryWindow(mzLow, mzHigh, rtLow, rtHigh, imLow, imHigh, includeDecoys);
	}

	/// <summary>
	/// Returns the full <see cref="MslLibraryEntry"/> for the precursor at
	/// zero-based position <paramref name="precursorIdx"/> in the m/z-sorted index.
	///
	/// <para>
	/// In full-load mode the entry is retrieved from the in-memory array.
	/// In index-only mode fragment ions are loaded from disk on the first access to each
	/// entry; subsequent accesses for the same entry are served from the LRU cache inside
	/// <see cref="MslIndex"/> without additional disk I/O (until the entry is evicted).
	/// </para>
	///
	/// <para>
	/// The <see cref="MslPrecursorIndexEntry.PrecursorIdx"/> field returned by
	/// <see cref="QueryMzWindow"/> and <see cref="QueryWindow"/> is the correct value
	/// to pass here.
	/// </para>
	/// </summary>
	/// <param name="precursorIdx">
	///   Zero-based precursor index as stored in
	///   <see cref="MslPrecursorIndexEntry.PrecursorIdx"/>. Out-of-range values return
	///   <see langword="null"/>.
	/// </param>
	/// <returns>
	///   The fully populated <see cref="MslLibraryEntry"/>, or <see langword="null"/> when
	///   the index is out of range or the entry cannot be loaded.
	/// </returns>
	/// <exception cref="ObjectDisposedException">This instance has been disposed.</exception>
	public MslLibraryEntry? GetEntry(int precursorIdx)
	{
		ThrowIfDisposed();
		return _index!.GetEntry(precursorIdx);
	}

	// ── Bulk enumeration ──────────────────────────────────────────────────────

	/// <summary>
	/// Lazily enumerates all library entries in ascending precursor m/z order.
	///
	/// <para>
	/// In full-load mode each entry is already fully populated in memory.
	/// In index-only mode fragment ions are loaded from disk as each entry is yielded.
	/// Because the LRU cache has a finite capacity, do not retain references to previously
	/// yielded entries if memory is constrained — the fragments backing an evicted entry will
	/// not be reloaded automatically.
	/// </para>
	///
	/// <para>
	/// Iterating the sequence to completion while the library is disposed is safe — the
	/// enumerator checks for disposal before yielding each entry and stops cleanly.
	/// </para>
	/// </summary>
	/// <param name="includeDecoys">
	///   When <see langword="true"/> (default), decoy entries are included in the sequence.
	///   When <see langword="false"/>, only target entries are yielded.
	/// </param>
	/// <returns>
	///   An <see cref="IEnumerable{T}"/> of <see cref="MslLibraryEntry"/> in ascending m/z
	///   order. Each entry has its <c>Fragments</c> list populated before being yielded.
	/// </returns>
	/// <exception cref="ObjectDisposedException">
	///   Thrown when the first call to <see cref="IEnumerator{T}.MoveNext"/> is made after
	///   this instance has been disposed.
	/// </exception>
	public IEnumerable<MslLibraryEntry> GetAllEntries(bool includeDecoys = true)
	{
		ThrowIfDisposed();

		// Snapshot the m/z-sorted index entries into an array before entering the iterator.
		// ReadOnlySpan<T> cannot be stored as a local in an iterator method (yield return)
		// under C# 12 (CS9202); copying to an array is the idiomatic workaround.
		MslPrecursorIndexEntry[] all = _index!.QueryMzRange(float.MinValue, float.MaxValue).ToArray();

		foreach (MslPrecursorIndexEntry indexEntry in all)
		{
			// Apply the decoy filter before triggering a potentially expensive disk read
			if (!includeDecoys && indexEntry.IsDecoy == 1)
				continue;

			MslLibraryEntry? entry = _index.GetEntry(indexEntry.PrecursorIdx);
			if (entry is not null)
				yield return entry;
		}
	}

	// ── RT calibration ────────────────────────────────────────────────────────

	/// <summary>
	/// Returns a new <see cref="MslLibrary"/> that shares the same fragment data as this
	/// instance but uses a new <see cref="MslIndex"/> with calibrated retention times.
	///
	/// <para>
	/// The calibration formula applied to every precursor's stored iRT value is:
	/// <code>calibratedRT = (float)(slope × iRT + intercept)</code>
	/// </para>
	///
	/// <para>
	/// <b>This instance is not modified.</b> The new library is independent and holds its
	/// own index. Dispose the new library when finished; it does not share a file stream with
	/// the original.
	/// </para>
	///
	/// <para>
	/// <b>Important:</b> in index-only mode, applying calibration internally rebuilds the
	/// sequence/charge dictionary, which would trigger on-demand fragment reads for every
	/// entry. It is therefore strongly recommended to call this method only on full-load
	/// libraries. See §9.3 of the Prompt 4 handoff for the full rationale.
	/// </para>
	/// </summary>
	/// <param name="slope">
	///   Slope of the linear iRT → calibrated RT regression. Typical values: 0.1–10.0.
	///   No guard against zero or NaN; callers are responsible for providing a well-conditioned
	///   regression fit.
	/// </param>
	/// <param name="intercept">
	///   Intercept of the linear iRT → calibrated RT regression (in the same units as the
	///   calibrated RT axis, typically minutes or seconds depending on the experiment).
	/// </param>
	/// <returns>
	///   A new <see cref="MslLibrary"/> whose precursor index entries have transformed
	///   <c>Irt</c> values. The new library has the same <see cref="Header"/> and
	///   <see cref="IsIndexOnly"/> flag as the original.
	/// </returns>
	/// <exception cref="ObjectDisposedException">This instance has been disposed.</exception>
	public MslLibrary WithCalibratedRetentionTimes(double slope, double intercept)
	{
		ThrowIfDisposed();

		// MslIndex.WithCalibratedRetentionTimes creates a brand-new index with transformed
		// Irt values; the original index (_index) is not modified.
		MslIndex calibratedIndex = _index!.WithCalibratedRetentionTimes(slope, intercept);

		return new MslLibrary(calibratedIndex, Header, IsIndexOnly, _rawLibrary);
	}

	// ── Elution group access ──────────────────────────────────────────────────

	/// <summary>
	/// Returns all <see cref="MslPrecursorIndexEntry"/> values that belong to the given
	/// elution group, i.e. all charge states and modification forms of a single peptide that
	/// were assigned to the same <see cref="MslPrecursorIndexEntry.ElutionGroupId"/>.
	///
	/// <para>
	/// The span is backed by a defensive copy allocated per call; it is therefore safe to
	/// retain across subsequent queries. The allocation is one
	/// <see cref="MslPrecursorIndexEntry"/> array per call; this is acceptable because elution
	/// group lookups are not on the microsecond hot path.
	/// </para>
	/// </summary>
	/// <param name="elutionGroupId">
	///   The elution group identifier assigned during library construction. All precursors
	///   sharing a stripped sequence receive the same identifier.
	/// </param>
	/// <returns>
	///   A <see cref="ReadOnlySpan{T}"/> of all matching index entries. Empty span when no
	///   precursors carry the given <paramref name="elutionGroupId"/>.
	/// </returns>
	/// <exception cref="ObjectDisposedException">This instance has been disposed.</exception>
	public ReadOnlySpan<MslPrecursorIndexEntry> GetElutionGroup(int elutionGroupId)
	{
		ThrowIfDisposed();
		return _index!.GetElutionGroup(elutionGroupId);
	}

	// ── IDisposable ───────────────────────────────────────────────────────────

	/// <summary>
	/// Releases all resources held by this <see cref="MslLibrary"/> instance.
	///
	/// <para>
	/// In <b>full-load mode</b>: clears the LRU entry cache inside <see cref="MslIndex"/> and
	/// prevents further queries. No file handle is involved.
	/// </para>
	///
	/// <para>
	/// In <b>index-only mode</b>: additionally closes and releases the
	/// <see cref="System.IO.FileStream"/> that backs on-demand fragment reads. On Windows this
	/// is necessary to allow the file to be moved, deleted, or written by other processes.
	/// </para>
	///
	/// <para>
	/// <b>Idempotent:</b> safe to call multiple times. Subsequent calls after the first are
	/// silently ignored.
	/// </para>
	///
	/// <para>
	/// After disposal, all public methods except <see cref="Dispose"/> throw
	/// <see cref="ObjectDisposedException"/>.
	/// </para>
	/// </summary>
	public void Dispose()
	{
		if (_disposed)
			return;

		_disposed = true;

		// Dispose the index first to clear the LRU cache and mark it unusable.
		// The index holds a reference to the entry-loader delegate which may close over
		// _rawLibrary; we therefore dispose _rawLibrary after the index.
		_index?.Dispose();
		_index = null;

		// In index-only mode, dispose the raw data container, which closes the FileStream.
		// This is safe to call even if _rawLibrary was never used after index construction
		// because MslLibraryData.Dispose is itself idempotent.
		_rawLibrary?.Dispose();
	}
}