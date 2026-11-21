using Omics.SpectralMatch.MslSpectralLibrary;

namespace Readers.SpectralLibrary;

// ────────────────────────────────────────────────────────────────────────────
// MslMergeConflictPolicy
// ────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Determines which entry is kept when the same Name appears in multiple source files.
/// </summary>
public enum MslMergeConflictPolicy
{
	/// <summary>
	/// Keep the entry from the source file with the lowest index in the input list
	/// (i.e. the first file that provided this key wins).
	/// Default policy — predictable and fast.
	/// </summary>
	KeepFirst,

	/// <summary>
	/// Keep the entry from the source file with the highest index in the input list
	/// (i.e. the last file that provided this key wins — later files override earlier ones).
	/// Useful for incremental library updates.
	/// </summary>
	KeepLast,

	/// <summary>
	/// Keep the entry with the lowest q-value (highest confidence).
	/// When q-values are NaN or equal, falls back to KeepFirst.
	/// Useful when merging empirical and predicted libraries.
	/// </summary>
	KeepLowestQValue
}

// ────────────────────────────────────────────────────────────────────────────
// MslMergeResult
// ────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Describes the outcome of an <see cref="MslMerger.Merge"/> call.
/// </summary>
public record MslMergeResult
{
	/// <summary>Path of the output file that was written.</summary>
	public required string OutputPath { get; init; }

	/// <summary>Total number of entries written to the output file.</summary>
	public required int OutputEntryCount { get; init; }

	/// <summary>Total number of entries read across all source files.</summary>
	public required int TotalSourceEntryCount { get; init; }

	/// <summary>
	/// Number of entries skipped due to duplicate Name resolution.
	/// Zero when <c>deduplicate: false</c>.
	/// </summary>
	public required int DuplicatesSkipped { get; init; }

	/// <summary>
	/// Source file paths that were found to contain entries not in m/z ascending order.
	/// The k-way merge still produces correct output; only the global m/z sort is not
	/// guaranteed when this list is non-empty.
	/// </summary>
	public required IReadOnlyList<string> UnsortedSourceFiles { get; init; }

	/// <summary>
	/// Per-source-file entry counts, in the same order as the input paths.
	/// </summary>
	public required IReadOnlyList<int> SourceEntryCounts { get; init; }
}

// ────────────────────────────────────────────────────────────────────────────
// MslMerger
// ────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Merges multiple .msl spectral library files into a single output file using a
/// k-way merge algorithm.
///
/// <para>
/// <b>Source-read memory</b>: O(k) — only k entries live in the heap at any one time.
/// Source files are read via <see cref="MslLibrary.LoadIndexOnly"/> and enumerated with
/// <see cref="MslLibrary.GetAllEntries"/>. Fragment bytes are loaded on demand per entry.
/// </para>
///
/// <para>
/// <b>Output memory</b>: O(spill file) — winning entries are yielded lazily to
/// <see cref="MslWriter.WriteStreaming"/> via an iterator method and are never
/// accumulated in a list.  WriteStreaming spills precursor scalar fields to a temp
/// file during Pass 1 (56 bytes × N entries on disk) and reads the source enumerable
/// exactly once.  Peak RAM during the merge is determined by the source heap (k entries)
/// plus the string and protein deduplication tables in WriteStreaming.
/// </para>
///
/// <para>
/// <b>Winner pre-scan</b>: when <c>deduplicate</c> is true, a first pass reads each
/// source's precursor metadata (no fragment bytes) and records, for every Name, which
/// copy the conflict policy keeps. The merge pass then writes exactly that copy. The
/// pre-scan holds one small record per unique Name. The winner never depends on
/// precursor m/z, so two copies of one Name whose m/z values differ (or tie) still
/// resolve by file order and q-value as the policy states.
/// </para>
///
/// <para>
/// Output entries are sorted by precursor m/z ascending (k-way merge order) unless any
/// source file is internally unsorted, in which case a warning is recorded in
/// <see cref="MslMergeResult.UnsortedSourceFiles"/>.
/// </para>
/// </summary>
public static class MslMerger
{
	/// <summary>
	/// Merges multiple .msl files into a single output .msl file.
	/// </summary>
	/// <param name="sourcePaths">
	///   Paths to the source .msl files to merge. Must not be null or empty.
	///   Files are merged in the order supplied; this order is significant for
	///   <see cref="MslMergeConflictPolicy.KeepFirst"/> and
	///   <see cref="MslMergeConflictPolicy.KeepLast"/>.
	/// </param>
	/// <param name="outputPath">
	///   Destination .msl file path. Created or overwritten. Must not be null.
	/// </param>
	/// <param name="conflictPolicy">
	///   Policy for handling duplicate LookupKeys across source files.
	///   Default: <see cref="MslMergeConflictPolicy.KeepFirst"/>.
	/// </param>
	/// <param name="deduplicate">
	///   When <see langword="true"/> (default), each unique Name is written only once.
	///   When <see langword="false"/>, all entries are written including duplicates.
	/// </param>
	/// <param name="compressionLevel">
	///   zstd compression level for the output (0 = no compression, 1–22 = zstd).
	/// </param>
	/// <returns>
	///   An <see cref="MslMergeResult"/> containing counts and diagnostic information.
	/// </returns>
	/// <exception cref="ArgumentNullException">
	///   Thrown when <paramref name="sourcePaths"/> or <paramref name="outputPath"/> is null.
	/// </exception>
	/// <exception cref="ArgumentException">
	///   Thrown when <paramref name="sourcePaths"/> is empty.
	/// </exception>
	public static MslMergeResult Merge(
	IReadOnlyList<string> sourcePaths,
	string outputPath,
	MslMergeConflictPolicy conflictPolicy = MslMergeConflictPolicy.KeepFirst,
	bool deduplicate = true,
	int compressionLevel = 0)
	{
		if (sourcePaths is null) throw new ArgumentNullException(nameof(sourcePaths));
		if (outputPath is null) throw new ArgumentNullException(nameof(outputPath));
		if (sourcePaths.Count == 0)
			throw new ArgumentException("At least one source path is required.", nameof(sourcePaths));
		MslWriter.ThrowIfInvalidCompressionLevel(compressionLevel, nameof(compressionLevel));
		ThrowIfOutputIsASource(sourcePaths, outputPath);

		if (deduplicate && !Enum.IsDefined(conflictPolicy))
			throw new ArgumentOutOfRangeException(nameof(conflictPolicy), conflictPolicy, "Unknown conflict policy.");

		int k = sourcePaths.Count;

		// ── Open source libraries ─────────────────────────────────────────────────
		MslLibrary[] libraries = new MslLibrary[k];
		IEnumerator<MslLibraryEntry>[] enumerators = new IEnumerator<MslLibraryEntry>[k];
		int[] sourceCounts = new int[k];
		float[] lastMz = new float[k];
		var unsortedSourceFiles = new List<string>();

		try
		{
			for (int i = 0; i < k; i++)
				libraries[i] = MslLibrary.LoadIndexOnly(sourcePaths[i]);

			// ── Winner pre-scan (every policy), over the libraries already open ────
			Dictionary<string, MergeWinner>? winners = deduplicate
				? BuildWinnerMap(libraries, conflictPolicy)
				: null;

			for (int i = 0; i < k; i++)
			{
				enumerators[i] = libraries[i].GetAllEntries().GetEnumerator();
				lastMz[i] = float.MinValue;
			}

			// ── Shared mutable state passed into the iterator ─────────────────────
			// We use a tiny reference-type wrapper so the iterator can write back
			// the counters that Merge() reports in MslMergeResult.
			var state = new MslMergeState();

			// ── Build the lazy merge sequence and feed it directly to WriteStreaming ─
			// WriteStreaming is now single-pass (Fix9a), so this IEnumerable is
			// consumed exactly once — no List materialisation needed.
			IEnumerable<MslLibraryEntry> mergedSequence = MergeEntries(
				sourcePaths, k, enumerators, lastMz, unsortedSourceFiles,
				sourceCounts, winners, state);

			MslWriter.WriteStreaming(outputPath, mergedSequence, compressionLevel);

			// sourceCounts are incremented inside TryAdvance (called from MergeEntries);
			// state.OutputEntryCount and state.DuplicatesSkipped are incremented inside
			// MergeEntries as entries are yielded / skipped.
			int totalSourceCount = 0;
			for (int i = 0; i < k; i++) totalSourceCount += sourceCounts[i];

			return new MslMergeResult
			{
				OutputPath = outputPath,
				OutputEntryCount = state.OutputEntryCount,
				TotalSourceEntryCount = totalSourceCount,
				DuplicatesSkipped = state.DuplicatesSkipped,
				UnsortedSourceFiles = unsortedSourceFiles.Distinct().ToList(),
				SourceEntryCounts = sourceCounts
			};
		}
		finally
		{
			for (int i = 0; i < k; i++)
			{
				try { enumerators[i]?.Dispose(); } catch { }
				try { libraries[i]?.Dispose(); } catch { }
			}
		}
	}

	// ────────────────────────────────────────────────────────────────────────
	// Private helpers
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>The copy of one Name that the conflict policy keeps.</summary>
	/// <param name="SourceIdx">Index of the source file in the input list.</param>
	/// <param name="Ordinal">
	///   Zero-based position of the copy in that source's <see cref="MslLibrary.GetAllEntries"/> order.
	/// </param>
	/// <param name="QValue">The copy's q-value (NaN when it has none).</param>
	private readonly record struct MergeWinner(int SourceIdx, int Ordinal, float QValue);

	/// <summary>
	/// First pass: for every Name, records the copy that <paramref name="conflictPolicy"/> keeps.
	/// Sources are read in input order and, within a source, in <see cref="MslLibrary.GetAllEntries"/>
	/// order, from precursor metadata only, so no fragment bytes are read. It reads the libraries
	/// the merge has already opened, so no source is opened twice.
	/// <list type="bullet">
	///   <item>KeepFirst: the first copy seen.</item>
	///   <item>KeepLast: the last copy seen.</item>
	///   <item>KeepLowestQValue: the copy with the lowest q-value. A NaN q-value never replaces
	///   a number, and a tie keeps the earlier copy (KeepFirst).</item>
	/// </list>
	/// Precursor m/z plays no part, so copies of one Name whose m/z values differ still resolve
	/// by the policy.
	/// </summary>
	private static Dictionary<string, MergeWinner> BuildWinnerMap(
		MslLibrary[] libraries, MslMergeConflictPolicy conflictPolicy)
	{
		var map = new Dictionary<string, MergeWinner>(StringComparer.Ordinal);

		for (int s = 0; s < libraries.Length; s++)
		{
			MslLibrary lib = libraries[s];

			// The index walk GetAllEntries makes, so each ordinal matches the merge pass.
			ReadOnlySpan<MslPrecursorIndexEntry> all = lib.QueryMzWindow(float.MinValue, float.MaxValue);
			int ordinal = 0;
			foreach (MslPrecursorIndexEntry indexEntry in all)
			{
				MslLibraryEntry? entry = lib.GetEntry(indexEntry.PrecursorIdx, includeFragments: false);
				if (entry is null)
					continue;

				var candidate = new MergeWinner(s, ordinal++, entry.QValue);
				if (!map.TryGetValue(entry.Name, out MergeWinner current))
				{
					map[entry.Name] = candidate;
					continue;
				}

				bool replace = conflictPolicy switch
				{
					MslMergeConflictPolicy.KeepLast => true,
					MslMergeConflictPolicy.KeepLowestQValue =>
						!float.IsNaN(candidate.QValue)
						&& (float.IsNaN(current.QValue) || candidate.QValue < current.QValue),
					_ => false // KeepFirst
				};
				if (replace)
					map[entry.Name] = candidate;
			}
		}

		return map;
	}

	/// <summary>
	/// Refuses an output path that names one of the sources. Every source stays open for the
	/// whole merge, so writing over one fails on Windows (the file is locked) and elsewhere
	/// replaces a file that is still being read.
	/// </summary>
	private static void ThrowIfOutputIsASource(IReadOnlyList<string> sourcePaths, string outputPath)
	{
		StringComparison comparison = OperatingSystem.IsWindows() || OperatingSystem.IsMacOS()
			? StringComparison.OrdinalIgnoreCase
			: StringComparison.Ordinal;
		string output = Path.GetFullPath(outputPath);
		foreach (string source in sourcePaths)
		{
			if (source is not null && string.Equals(Path.GetFullPath(source), output, comparison))
				throw new ArgumentException(
					$"The output path is also one of the sources ({outputPath}). Write the merge to a new file.",
					nameof(outputPath));
		}
	}

	/// <summary>
	/// Core k-way merge iterator. Yields entries in ascending precursor m/z; entries of equal m/z
	/// come out in source order, so the output is deterministic. When <paramref name="winners"/>
	/// is not null, only the copy of each Name that <see cref="BuildWinnerMap"/> chose is yielded.
	/// Increments counters in <paramref name="state"/> as it runs so that
	/// <see cref="Merge"/> can report accurate totals in <see cref="MslMergeResult"/>
	/// after <see cref="MslWriter.WriteStreaming"/> has consumed the sequence.
	/// </summary>
	private static IEnumerable<MslLibraryEntry> MergeEntries(
		IReadOnlyList<string> sourcePaths,
		int k,
		IEnumerator<MslLibraryEntry>[] enumerators,
		float[] lastMz,
		List<string> unsortedSourceFiles,
		int[] sourceCounts,
		Dictionary<string, MergeWinner>? winners,
		MslMergeState state)
	{
		// k-way min-heap keyed on (PrecursorMz, source index). At most one entry per source is
		// queued, so the pair orders every queued entry.
		var heap = new PriorityQueue<(MslLibraryEntry Entry, int SourceIdx, int Ordinal), (float Mz, int SourceIdx)>();

		// Seed the heap with one entry from each source
		for (int i = 0; i < k; i++)
		{
			if (TryAdvance(enumerators[i], sourcePaths[i], ref lastMz[i],
						   unsortedSourceFiles, ref sourceCounts[i], out MslLibraryEntry? first))
				heap.Enqueue((first!, i, sourceCounts[i] - 1), ((float)first!.PrecursorMz, i));
		}

		while (heap.TryDequeue(out var item, out _))
		{
			var (entry, sourceIdx, ordinal) = item;

			// Advance this source and push its next entry onto the heap
			if (TryAdvance(enumerators[sourceIdx], sourcePaths[sourceIdx], ref lastMz[sourceIdx],
						   unsortedSourceFiles, ref sourceCounts[sourceIdx], out MslLibraryEntry? next))
				heap.Enqueue((next!, sourceIdx, sourceCounts[sourceIdx] - 1), ((float)next!.PrecursorMz, sourceIdx));

			if (winners is not null
				&& !(winners.TryGetValue(entry.Name, out MergeWinner winner)
					 && winner.SourceIdx == sourceIdx
					 && winner.Ordinal == ordinal))
			{
				state.DuplicatesSkipped++;
				continue;
			}

			state.OutputEntryCount++;
			yield return entry;
		}
	}

	/// <summary>
	/// Advances the given enumerator by one position, updating per-source state.
	/// Detects out-of-order m/z and records the source path in <paramref name="unsortedList"/>
	/// if this is the first unsorted violation for that source.
	/// Returns <see langword="true"/> when a new entry was available; <see langword="false"/>
	/// when the source is exhausted.
	/// </summary>
	private static bool TryAdvance(
		IEnumerator<MslLibraryEntry> enumerator,
		string sourcePath,
		ref float lastMzForSource,
		List<string> unsortedList,
		ref int sourceCount,
		out MslLibraryEntry? entry)
	{
		if (!enumerator.MoveNext())
		{
			entry = null;
			return false;
		}

		entry = enumerator.Current;
		sourceCount++;

		float mz = (float)entry.PrecursorMz;

		// Detect unsorted source (m/z decreased from the previous entry in this source)
		if (mz < lastMzForSource && !unsortedList.Contains(sourcePath))
			unsortedList.Add(sourcePath);

		lastMzForSource = mz;
		return true;
	}

	/// <summary>
	/// Mutable counter bag shared between <see cref="MslMerger.Merge"/> and the
	/// <see cref="MslMerger.MergeEntries"/> iterator so that post-write result
	/// construction can report accurate totals.
	///
	/// Using a reference type avoids the "iterator cannot use ref parameters"
	/// restriction that would prevent passing <c>ref int</c> counters into an
	/// iterator method.
	/// </summary>
	internal sealed class MslMergeState
	{
		/// <summary>Total number of entries yielded to WriteStreaming.</summary>
		public int OutputEntryCount { get; set; }

		/// <summary>Total number of entries discarded by the conflict policy.</summary>
		public int DuplicatesSkipped { get; set; }
	}
}