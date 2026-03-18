using Omics.SpectralMatch.MslSpectralLibrary;

namespace Readers.SpectralLibrary;

/// <summary>
/// In-memory container for a loaded .msl spectral library. Returned by both
/// <see cref="MslReader.Load"/> (full-load mode) and <see cref="MslReader.LoadIndexOnly"/>
/// (index-only mode). The public API is identical in both modes; the difference is whether
/// fragment ions are pre-populated in each entry's <c>Fragments</c> list or must be fetched
/// on demand via <see cref="LoadFragmentsOnDemand"/>.
///
/// <para>
/// In index-only mode the object holds an open <see cref="System.IO.FileStream"/> to support
/// on-demand fragment reads. Callers <b>must</b> call <see cref="Dispose"/> when finished to
/// release the file handle. On Windows this is mandatory because the OS uses file locking.
/// </para>
///
/// <para>
/// In full-load mode <see cref="Dispose"/> is a no-op but is still safe to call.
/// </para>
/// </summary>
public sealed class MslLibraryData : IDisposable
{
	// ── Private state ─────────────────────────────────────────────────────────

	/// <summary>
	/// Open file stream held in index-only mode for on-demand fragment reads.
	/// Set to null after <see cref="Dispose"/> is called, or in full-load mode.
	/// </summary>
	private FileStream? _onDemandStream;

	/// <summary>
	/// Synchronization lock that serialises seek+read sequences against
	/// <see cref="_onDemandStream"/> so concurrent callers of
	/// <see cref="LoadFragmentsOnDemand"/> do not interleave their I/O.
	/// </summary>
	private readonly object _streamLock = new();

	/// <summary>
	/// Snapshot of all raw precursor records kept for index-only on-demand reads.
	/// Provides <c>FragmentBlockOffset</c> and <c>FragmentCount</c> per precursor
	/// without re-parsing the file. Null in full-load mode.
	/// </summary>
	private readonly MslPrecursorRecord[]? _precursorRecords;

	/// <summary>
	/// Snapshot of the complete string table, kept so <see cref="LoadFragmentsOnDemand"/>
	/// can resolve string indices for any newly-constructed entries if needed in future.
	/// Null in full-load mode (strings are already baked into the entry objects).
	/// </summary>
	private readonly string[]? _strings;

	/// <summary>
	/// Snapshot of the complete protein table, kept for symmetry with
	/// <see cref="_strings"/> and potential future use. Null in full-load mode.
	/// </summary>
	private readonly MslProteinRecord[]? _proteins;

	/// <summary>
	/// Extended annotation table masses read from the file's optional extended annotation
	/// section (format version 2+, <see cref="MslFormat.FileFlagHasExtAnnotations"/>).
	/// Index 0 is the reserved sentinel (0.0 = no loss); valid custom entries start at 1.
	/// Empty array when the file contains no custom neutral losses (version 1 or flag absent).
	/// Retained so on-demand fragment reads can decode custom-loss indices without re-opening
	/// the file.
	/// </summary>
	private readonly double[] _customLossMasses;

	/// <summary>
	/// True when this instance was created in index-only mode; false for full-load mode.
	/// Controls the behaviour of <see cref="LoadFragmentsOnDemand"/>.
	/// </summary>
	private readonly bool _isIndexOnly;

	// ── Construction ──────────────────────────────────────────────────────────

	/// <summary>
	/// Constructs a full-load library instance. All fragments are already populated in
	/// <paramref name="entries"/>. No file stream is retained.
	/// </summary>
	/// <param name="entries">All precursor entries with their fragment ions fully loaded.</param>
	/// <param name="header">Deserialized file header providing version and flag metadata.</param>
	/// <exception cref="ArgumentNullException"><paramref name="entries"/> is null.</exception>
	internal MslLibraryData(List<MslLibraryEntry> entries, MslFileHeader header)
	{
		Entries = entries ?? throw new ArgumentNullException(nameof(entries));
		Header = header;
		_customLossMasses = Array.Empty<double>();
		_isIndexOnly = false;
	}

	/// <summary>
	/// Constructs an index-only library instance. Entries have empty <c>Fragments</c> lists;
	/// the retained <paramref name="onDemandStream"/> is used to fetch fragments on demand.
	/// Ownership of the stream is transferred to this instance and will be closed by
	/// <see cref="Dispose"/>.
	/// </summary>
	/// <param name="entries">Skeleton precursor entries (metadata populated, fragment lists empty).</param>
	/// <param name="header">Deserialized file header.</param>
	/// <param name="precursorRecords">
	/// Raw precursor records retained for on-demand fragment reads; provides
	/// <c>FragmentBlockOffset</c> and <c>FragmentCount</c> per entry.
	/// </param>
	/// <param name="strings">Complete string table snapshot.</param>
	/// <param name="proteins">Complete protein table snapshot.</param>
	/// <param name="onDemandStream">
	/// An open <see cref="System.IO.FileStream"/> used to seek and read fragment blocks.
	/// Ownership is transferred; will be disposed by <see cref="Dispose"/>.
	/// </param>
	/// <param name="customLossMasses">
	/// Extended annotation table masses from the file (format version 2+). Index 0 is the
	/// reserved sentinel (0.0). Pass <see cref="Array.Empty{T}"/> for version-1 files or
	/// files that contain no custom neutral losses.
	/// </param>
	/// <exception cref="ArgumentNullException">Any required parameter is null.</exception>
	internal MslLibraryData(
		List<MslLibraryEntry> entries,
		MslFileHeader header,
		MslPrecursorRecord[] precursorRecords,
		string[] strings,
		MslProteinRecord[] proteins,
		FileStream onDemandStream,
		double[] customLossMasses)
	{
		Entries = entries ?? throw new ArgumentNullException(nameof(entries));
		Header = header;
		_precursorRecords = precursorRecords ?? throw new ArgumentNullException(nameof(precursorRecords));
		_strings = strings ?? throw new ArgumentNullException(nameof(strings));
		_proteins = proteins ?? throw new ArgumentNullException(nameof(proteins));
		_onDemandStream = onDemandStream ?? throw new ArgumentNullException(nameof(onDemandStream));
		_customLossMasses = customLossMasses ?? Array.Empty<double>();
		_isIndexOnly = true;
	}

	// ── Public properties ─────────────────────────────────────────────────────

	/// <summary>
	/// All precursor entries in the library, in file order. In full-load mode each entry's
	/// <c>Fragments</c> list is fully populated. In index-only mode the list is empty until
	/// <see cref="LoadFragmentsOnDemand"/> is called for that entry's index.
	/// </summary>
	public IReadOnlyList<MslLibraryEntry> Entries { get; }

	/// <summary>
	/// The deserialized file header, exposing format version, precursor count, file-level
	/// flags, and the byte offsets of each section within the file.
	/// </summary>
	public MslFileHeader Header { get; }

	/// <summary>
	/// Total number of precursor entries. Equivalent to <c>Entries.Count</c>.
	/// Available in both full-load and index-only mode without loading any fragment data.
	/// </summary>
	public int Count => Entries.Count;

	/// <summary>
	/// True when this instance was loaded via <see cref="MslReader.LoadIndexOnly"/>.
	/// In this mode the fragment lists in <see cref="Entries"/> are empty until fragments are
	/// retrieved explicitly via <see cref="LoadFragmentsOnDemand"/>.
	/// False when loaded via <see cref="MslReader.Load"/> (full-load mode).
	/// </summary>
	public bool IsIndexOnly => _isIndexOnly;

	/// <summary>
	/// True when the library file was written with zstd block compression
	/// (<see cref="MslFormat.FileFlagIsCompressed"/> is set in <see cref="Header"/>).
	/// Compressed libraries always load in full-load mode regardless of which
	/// <see cref="MslReader"/> method was called; <see cref="IsIndexOnly"/> is always
	/// <c>false</c> for compressed files.
	/// </summary>
	public bool IsCompressed => (Header.FileFlags & MslFormat.FileFlagIsCompressed) != 0;

	// ── On-demand fragment loading ────────────────────────────────────────────

	/// <summary>
	/// Fetches the fragment ions for the precursor at zero-based index
	/// <paramref name="precursorIndex"/> by seeking to its stored fragment block offset
	/// and reading exactly <c>FragmentCount × 20</c> bytes.
	///
	/// <para>
	/// This method is only available in index-only mode. In full-load mode the fragment list
	/// is already populated in <see cref="Entries"/>[<paramref name="precursorIndex"/>].Fragments;
	/// calling this method on a full-load instance throws <see cref="InvalidOperationException"/>.
	/// </para>
	///
	/// <para>Thread-safe: multiple callers may request different precursors concurrently.
	/// The seek+read sequence is protected by an internal lock on the shared stream.</para>
	/// </summary>
	/// <param name="precursorIndex">
	/// Zero-based index into <see cref="Entries"/>. Must be in range [0, Count).
	/// </param>
	/// <returns>
	/// A new list of <see cref="MslFragmentIon"/> objects for the requested precursor,
	/// in m/z ascending order as written by the writer.
	/// </returns>
	/// <exception cref="InvalidOperationException">
	/// Called on a full-load instance (<see cref="IsIndexOnly"/> == false).
	/// </exception>
	/// <exception cref="ArgumentOutOfRangeException">
	/// <paramref name="precursorIndex"/> is outside [0, Count).
	/// </exception>
	/// <exception cref="ObjectDisposedException">
	/// <see cref="Dispose"/> has already been called.
	/// </exception>
	public List<MslFragmentIon> LoadFragmentsOnDemand(int precursorIndex)
	{
		if (!_isIndexOnly)
			throw new InvalidOperationException(
				"LoadFragmentsOnDemand is only available in index-only mode. " +
				"Use MslReader.LoadIndexOnly() to open a library in this mode.");

		if (precursorIndex < 0 || precursorIndex >= Count)
			throw new ArgumentOutOfRangeException(nameof(precursorIndex),
				$"Index {precursorIndex} is out of range [0, {Count}).");

		// Serialise the seek+read so concurrent callers on different threads do not interleave
		lock (_streamLock)
		{
			if (_onDemandStream is null)
				throw new ObjectDisposedException(nameof(MslLibraryData),
					"Cannot load fragments: the library has been disposed.");

			return MslReader.ReadFragmentBlockFromStream(
				_onDemandStream,
				_precursorRecords![precursorIndex],
				_customLossMasses);
		}
	}

	// ── IDisposable ───────────────────────────────────────────────────────────

	/// <summary>
	/// Releases the open file stream held in index-only mode. After disposal,
	/// <see cref="LoadFragmentsOnDemand"/> will throw <see cref="ObjectDisposedException"/>.
	///
	/// Safe to call on a full-load instance (no-op in that case).
	/// Safe to call multiple times (idempotent).
	/// </summary>
	public void Dispose()
	{
		// Null the stream inside the lock so in-flight LoadFragmentsOnDemand callers
		// that are waiting will see null and throw ObjectDisposedException rather than
		// attempting to use an already-closed stream.
		lock (_streamLock)
		{
			_onDemandStream?.Dispose();
			_onDemandStream = null;
		}
	}
}