using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;

namespace Readers.SpectralLibrary;

/// <summary>
/// Converts text-based spectral library formats (MSP, pDeep, ms2pip) into fully in-memory
/// <see cref="MslLibrary"/> objects without writing any intermediate binary file to disk.
///
/// <para>
/// The returned <see cref="MslLibrary"/> is identical in behavior to one produced by
/// <see cref="MslLibrary.Load"/>: all precursor entries and fragment ions are held in memory,
/// the m/z-sorted <see cref="MslIndex"/> is fully constructed, and all query APIs
/// (<see cref="MslLibrary.QueryWindow"/>, <see cref="MslLibrary.TryGetEntry"/>,
/// <see cref="MslLibrary.QueryMzWindow"/>, <see cref="MslLibrary.GetAllEntries"/>, etc.)
/// are immediately available. No file handle is held and <see cref="MslLibrary.IsIndexOnly"/>
/// is always <see langword="false"/> for libraries produced by this class.
/// </para>
///
/// <para>
/// <b>Typical usage (DDA/DIA search startup):</b>
/// <code>
/// MslLibrary lib = MslConverter.FromMspFile("myLibrary.msp");
/// // lib is ready for QueryWindow, TryGetEntry, etc. — no .msl file written
/// </code>
/// </para>
///
/// <para>
/// <b>Thread safety:</b> both public methods are stateless and safe to call from multiple
/// threads. The returned <see cref="MslLibrary"/> inherits the thread-safety guarantees
/// of that class (concurrent reads are safe; disposal must not race with queries).
/// </para>
/// </summary>
public static class MslConverter
{
	// ── Public API ────────────────────────────────────────────────────────────

	/// <summary>
	/// Parses a text-based spectral library file (MSP, pDeep, or ms2pip) and returns a
	/// fully in-memory <see cref="MslLibrary"/>. No binary .msl file is written to disk.
	///
	/// <para>
	/// All spectra present in the file are loaded eagerly. For very large text libraries
	/// the parse cost is paid once here; all subsequent queries run against the in-memory
	/// index with the same performance as a binary-loaded library.
	/// </para>
	///
	/// <para>
	/// Format detection follows the same convention used by <see cref="SpectralLibrary"/>:
	/// files whose path contains "pdeep" are parsed as pDeep format; files whose path
	/// contains "ms2pip" are parsed as ms2pip format; all other files are parsed as MSP.
	/// </para>
	/// </summary>
	/// <param name="textLibraryPath">
	///   Absolute or relative path to the text spectral library file (.msp or equivalent).
	///   Must not be null or empty. The file must exist and be readable.
	/// </param>
	/// <returns>
	///   A fully populated <see cref="MslLibrary"/> in full-load mode containing all
	///   spectra parsed from the file. <see cref="MslLibrary.IsIndexOnly"/> is
	///   <see langword="false"/>. No file handle is retained after this method returns.
	/// </returns>
	/// <exception cref="ArgumentNullException">
	///   <paramref name="textLibraryPath"/> is null or empty.
	/// </exception>
	/// <exception cref="FileNotFoundException">
	///   No file exists at <paramref name="textLibraryPath"/>.
	/// </exception>
	/// <exception cref="InvalidDataException">
	///   The file could not be parsed as a valid spectral library (e.g. empty file,
	///   malformed entries, or no spectra found).
	/// </exception>
	public static MslLibrary FromMspFile(string textLibraryPath)
	{
		if (string.IsNullOrEmpty(textLibraryPath))
			throw new ArgumentNullException(nameof(textLibraryPath),
				"Text library path must not be null or empty.");

		if (!File.Exists(textLibraryPath))
			throw new FileNotFoundException(
				$"Text spectral library file not found: '{textLibraryPath}'.", textLibraryPath);

		// Use the existing SpectralLibrary text parser to enumerate all spectra.
		// SpectralLibrary does not implement IDisposable; call CloseConnections() manually
		// in a finally block to ensure StreamReaders are released even on parse failure.
		var textLib = new SpectralLibrary(new List<string> { textLibraryPath });
		List<LibrarySpectrum> spectra;
		try
		{
			spectra = textLib.GetAllLibrarySpectra().ToList();
		}
		finally
		{
			textLib.CloseConnections();
		}

		if (spectra.Count == 0)
			throw new InvalidDataException(
				$"No spectra could be parsed from '{textLibraryPath}'. " +
				$"Verify the file is a valid MSP, pDeep, or ms2pip spectral library.");

		return FromMspLibrarySpectra(spectra);
	}

	/// <summary>
	/// Converts an already-parsed collection of <see cref="LibrarySpectrum"/> objects into
	/// a fully in-memory <see cref="MslLibrary"/>. No file I/O of any kind is performed.
	///
	/// <para>
	/// This overload is useful when the caller has already obtained spectra from Koina
	/// predictions, MetaMorpheus in-silico generation, or any other source that produces
	/// <see cref="LibrarySpectrum"/> objects directly, and wants the full
	/// <see cref="MslLibrary"/> query API without persisting a binary file.
	/// </para>
	///
	/// <para>
	/// Each <see cref="LibrarySpectrum"/> is converted via
	/// <see cref="MslLibraryEntry.FromLibrarySpectrum"/>. Fields not present in
	/// <see cref="LibrarySpectrum"/> (ion mobility, protein accession, molecule type, NCE,
	/// dissociation type, q-value) receive the same sensible defaults as
	/// <see cref="MslLibraryEntry.FromLibrarySpectrum"/> applies.
	/// </para>
	///
	/// <para>
	/// Duplicate sequence/charge keys within <paramref name="spectra"/> are handled by the
	/// <see cref="MslIndex"/> LRU cache: the last entry with a given key wins in the
	/// sequence/charge dictionary, consistent with how the binary reader handles duplicates
	/// in a file.
	/// </para>
	/// </summary>
	/// <param name="spectra">
	///   The spectra to convert. Must not be null. May be empty, in which case an empty
	///   <see cref="MslLibrary"/> with a precursor count of zero is returned.
	/// </param>
	/// <returns>
	///   A fully populated <see cref="MslLibrary"/> in full-load mode.
	///   <see cref="MslLibrary.IsIndexOnly"/> is <see langword="false"/>.
	///   The library is ready for all query operations immediately.
	/// </returns>
	/// <exception cref="ArgumentNullException">
	///   <paramref name="spectra"/> is null.
	/// </exception>
	public static MslLibrary FromMspLibrarySpectra(IReadOnlyList<LibrarySpectrum> spectra)
	{
		ArgumentNullException.ThrowIfNull(spectra);

		// Convert each LibrarySpectrum to the richer MslLibraryEntry model.
		// FromLibrarySpectrum handles modification notation, fragment ion conversion,
		// and neutral-loss classification.
		var entries = new List<MslLibraryEntry>(spectra.Count);
		foreach (LibrarySpectrum spectrum in spectra)
			entries.Add(MslLibraryEntry.FromLibrarySpectrum(spectrum));

		// Build a synthetic header that correctly describes the in-memory content.
		// Byte-offset fields are meaningless for a file-less library; they are zeroed.
		// The magic and version are set correctly so that any code path that inspects
		// the header for format identification behaves consistently.
		var header = new MslFileHeader
		{
			Magic = MslFormat.MagicAsUInt32,
			FormatVersion = MslFormat.CurrentVersion,
			FileFlags = 0,
			NPrecursors = entries.Count,
			NProteins = 0,
			NElutionGroups = 0,   // MslIndex.Build computes this automatically
			NStrings = 0,   // no string table; data lives in entry objects
			ProteinTableOffset = 0,
			StringTableOffset = 0,
			PrecursorSectionOffset = 0,
			FragmentSectionOffset = 0,
		};

		// Build the in-memory query index from the entry list.
		// The loader delegate is a direct array lookup — all fragments are already in memory.
		MslIndex index = MslIndex.Build(
			entries,
			i => i >= 0 && i < entries.Count ? entries[i] : null);

		// Invoke the internal MslLibrary constructor via the same pattern used by Load().
		// isIndexOnly = false  (everything is in memory)
		// rawLibrary  = null   (no FileStream to hold)
		return MslLibrary.CreateInMemory(index, header);
	}
}