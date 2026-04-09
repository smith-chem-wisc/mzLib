using Omics.SpectralMatch.MslSpectralLibrary;

namespace Readers.SpectralLibrary;

/// <summary>
/// Static helper that encapsulates all .msl file-type detection and opening logic.
/// Keeps the routing decision in one place so both <see cref="SpectralLibrary"/> and any
/// future callers can share identical behaviour without duplicating extension checks.
/// </summary>
public static class MslFileTypeHandler
{
	// ── Size threshold constant ───────────────────────────────────────────────

	/// <summary>
	/// Default file-size boundary (in bytes) above which <see cref="Open"/> switches from a
	/// full in-memory load to an index-only load.  1 GiB = 1 073 741 824 bytes.
	/// Callers may override this via the <c>indexOnlyThresholdBytes</c> parameter of
	/// <see cref="Open"/>.
	/// </summary>
	public const long DefaultIndexOnlyThresholdBytes = 1_073_741_824L;

	// ── File-type detection ───────────────────────────────────────────────────

	/// <summary>
	/// Returns true when <paramref name="filePath"/> has the <c>.msl</c> extension,
	/// regardless of case.  Returns false for all other extensions (including <c>.msp</c>,
	/// <c>.tsv</c>, <c>.parquet</c>, etc.) and for null or empty input.
	/// </summary>
	/// <param name="filePath">
	///   Absolute or relative path whose extension is inspected.  May be null or empty,
	///   in which case the method returns false rather than throwing.
	/// </param>
	/// <returns>
	///   <c>true</c> if the extension is <c>.msl</c> (case-insensitive); <c>false</c> otherwise.
	/// </returns>
	public static bool IsMslFile(string filePath)
	{
		// Guard against null / empty paths — return false rather than propagating an exception
		if (string.IsNullOrEmpty(filePath))
			return false;

		// Path.GetExtension returns the dot-prefixed extension (e.g. ".msl") or "" if none.
		// OrdinalIgnoreCase handles Windows paths like "Library.MSL" from Explorer drag-and-drop.
		string extension = Path.GetExtension(filePath);
		return string.Equals(extension, ".msl", StringComparison.OrdinalIgnoreCase);
	}

	// ── Library opening ───────────────────────────────────────────────────────

	/// <summary>
	/// Opens the .msl file at <paramref name="filePath"/> and returns the resulting
	/// <see cref="MslLibrary"/> instance.
	///
	/// <para>
	/// Load-mode selection is automatic based on file size:
	/// <list type="bullet">
	///   <item>
	///     Files <b>smaller than</b> <paramref name="indexOnlyThresholdBytes"/> are opened with
	///     <see cref="MslReader.Load"/>, which reads every precursor and all fragment blocks into
	///     memory.  No file handle is retained after the call returns.
	///   </item>
	///   <item>
	///     Files <b>at or above</b> the threshold are opened with
	///     <see cref="MslReader.LoadIndexOnly"/>, which reads only precursor metadata and fetches
	///     fragment ions from disk on demand.  The returned <see cref="MslLibrary"/> holds an open
	///     <see cref="FileStream"/> and <b>must</b> be disposed by the caller when no longer needed.
	///   </item>
	/// </list>
	/// </para>
	/// </summary>
	/// <param name="filePath">
	///   Absolute or relative path to a valid .msl binary spectral library file.
	///   The file must exist and be readable.  This method does not validate the extension —
	///   call <see cref="IsMslFile"/> first if extension validation is required.
	/// </param>
	/// <param name="indexOnlyThresholdBytes">
	///   File-size boundary in bytes.  Defaults to <see cref="DefaultIndexOnlyThresholdBytes"/>
	///   (1 GiB).  Pass a smaller value in unit tests to exercise index-only mode on small files.
	/// </param>
	/// <returns>
	///   A fully initialised <see cref="MslLibrary"/> in either full-load or index-only mode.
	///   In index-only mode the caller is responsible for disposing the returned
	///   <see cref="MslLibrary"/>, which will release the underlying file handle.
	/// </returns>
	/// <exception cref="FileNotFoundException">
	///   Thrown when <paramref name="filePath"/> does not exist on the file system.
	/// </exception>
	/// <exception cref="InvalidDataException">
	///   Thrown when the file does not begin with the MZLB magic bytes or has an unsupported
	///   format version.
	/// </exception>
	public static MslLibrary Open(string filePath,
		long indexOnlyThresholdBytes = DefaultIndexOnlyThresholdBytes)
	{
		// Verify existence up front so we surface a clear FileNotFoundException rather than
		// a cryptic IO error deep inside the reader
		if (!File.Exists(filePath))
			throw new FileNotFoundException($"MSL library file not found: {filePath}", filePath);

		// FileInfo.Length is a single cheap OS stat call — no need to open the file
		long fileSize = new FileInfo(filePath).Length;

		// Select load mode: index-only for large files, full load for everything else.
		// MslLibrary.Load / MslLibrary.LoadIndexOnly are the public facade entry points;
		// they wrap MslReader internally and return MslLibrary rather than MslLibraryData.
		return fileSize >= indexOnlyThresholdBytes
			? MslLibrary.LoadIndexOnly(filePath)
			: MslLibrary.Load(filePath);
	}
}