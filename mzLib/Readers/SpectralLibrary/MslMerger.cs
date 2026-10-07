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
/// For <see cref="MslMergeConflictPolicy.KeepLowestQValue"/>, an additional buffer
/// bounded by the number of entries within a 0.001 Da m/z window is held while winners
/// are resolved. This buffer is typically very small (at most a few dozen entries).
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

		int k = sourcePaths.Count;

		// ── KeepLast two-pass pre-scan ────────────────────────────────────────────
		Dictionary<string, int>? keepLastMap = null;
		if (deduplicate && conflictPolicy == MslMergeConflictPolicy.KeepLast)
			keepLastMap = BuildKeepLastMap(sourcePaths);

		// ── Open source libraries ─────────────────────────────────────────────────
		MslLibrary[] libraries = new MslLibrary[k];
		IEnumerator<MslLibraryEntry>[] enumerators = new IEnumerator<MslLibraryEntry>[k];
		int[] sourceCounts = new int[k];
		float[] lastMz = new float[k];
		var unsortedSourceFiles = new List<string>();

		try
		{
			for (int i = 0; i < k; i++)
			{
				libraries[i] = MslLibrary.LoadIndexOnly(sourcePaths[i]);
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
				sourceCounts, keepLastMap, deduplicate, conflictPolicy, state);

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

	/// <summary>
	/// Builds a dictionary mapping each Name to the index of the last source file
	/// that contains an entry with that key. Used for the KeepLast conflict policy.
	/// This is the first pass of the two-pass KeepLast algorithm; it requires O(UniqueKeys)
	/// RAM which is acceptable (~30 MB for 1M unique keys).
	/// </summary>
	private static Dictionary<string, int> BuildKeepLastMap(IReadOnlyList<string> sourcePaths)
	{
		var map = new Dictionary<string, int>(StringComparer.Ordinal);

		for (int i = 0; i < sourcePaths.Count; i++)
		{
			using var lib = MslLibrary.LoadIndexOnly(sourcePaths[i]);
			foreach (MslLibraryEntry entry in lib.GetAllEntries())
				map[entry.Name] = i;   // overwrite: last writer wins
		}

		return map;
	}
	/// <summary>
	/// Core k-way merge iterator.  Yields each winning <see cref="MslLibraryEntry"/>
	/// in precursor-m/z ascending order according to <paramref name="conflictPolicy"/>.
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
		Dictionary<string, int>? keepLastMap,
		bool deduplicate,
		MslMergeConflictPolicy conflictPolicy,
		MslMergeState state)
	{
		// k-way min-heap keyed on PrecursorMz
		var heap = new PriorityQueue<(MslLibraryEntry Entry, IEnumerator<MslLibraryEntry> Enum, int SourceIdx), float>();

		// Seed the heap with one entry from each source
		for (int i = 0; i < k; i++)
		{
			if (TryAdvance(enumerators[i], sourcePaths[i], ref lastMz[i],
						   unsortedSourceFiles, ref sourceCounts[i], out MslLibraryEntry? first))
				heap.Enqueue((first!, enumerators[i], i), (float)first!.PrecursorMz);
		}

		// KeepFirst: HashSet of already-yielded keys
		var emittedKeys = (deduplicate && conflictPolicy == MslMergeConflictPolicy.KeepFirst)
			? new HashSet<string>(StringComparer.Ordinal)
			: null;

		// KeepLowestQValue: pending buffer for entries within a narrow m/z window
		var qValueBuffer = (deduplicate && conflictPolicy == MslMergeConflictPolicy.KeepLowestQValue)
			? new Dictionary<string, (float QValue, MslLibraryEntry Entry)>(StringComparer.Ordinal)
			: null;

		while (heap.Count > 0)
		{
			heap.TryDequeue(out var item, out float _);
			var (entry, enumerator, sourceIdx) = item;

			// Advance this source and push its next entry onto the heap
			if (TryAdvance(enumerator, sourcePaths[sourceIdx], ref lastMz[sourceIdx],
						   unsortedSourceFiles, ref sourceCounts[sourceIdx], out MslLibraryEntry? next))
				heap.Enqueue((next!, enumerator, sourceIdx), (float)next!.PrecursorMz);

			if (!deduplicate)
			{
				state.OutputEntryCount++;
				yield return entry;
				continue;
			}

			string key = entry.Name;

			switch (conflictPolicy)
			{
				case MslMergeConflictPolicy.KeepFirst:
					if (emittedKeys!.Add(key))
					{
						state.OutputEntryCount++;
						yield return entry;
					}
					else
					{
						state.DuplicatesSkipped++;
					}
					break;

				case MslMergeConflictPolicy.KeepLast:
					// Only yield if this source is the last source that has this key
					if (keepLastMap!.TryGetValue(key, out int lastSource) && lastSource == sourceIdx)
					{
						state.OutputEntryCount++;
						yield return entry;
					}
					else
					{
						state.DuplicatesSkipped++;
					}
					break;

				case MslMergeConflictPolicy.KeepLowestQValue:
					// Flush buffer entries that are more than 0.001 Da behind current m/z
					foreach (MslLibraryEntry flushed in FlushQValueBuffer(
						qValueBuffer!, (float)entry.PrecursorMz, 0.001f))
					{
						state.OutputEntryCount++;
						yield return flushed;
					}

					// Update buffer with the better (lower) q-value entry for this key
					float q = entry.QValue;
					if (qValueBuffer!.TryGetValue(key, out var existing))
					{
						float existingQ = existing.QValue;
						bool replaceWithNew = (!float.IsNaN(q) && (float.IsNaN(existingQ) || q < existingQ));
						if (replaceWithNew)
							qValueBuffer[key] = (q, entry);
						// Either way, one duplicate was resolved
						state.DuplicatesSkipped++;
					}
					else
					{
						qValueBuffer[key] = (q, entry);
					}
					break;
			}
		}

		// Flush any remaining KeepLowestQValue buffer entries
		if (qValueBuffer is { Count: > 0 })
		{
			foreach (MslLibraryEntry flushed in FlushQValueBuffer(qValueBuffer, float.MaxValue, 0f))
			{
				state.OutputEntryCount++;
				yield return flushed;
			}
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
	/// Yields all KeepLowestQValue buffer entries whose precursor m/z is more than
	/// <paramref name="mzWindow"/> Da below <paramref name="currentMz"/>.  When
	/// <paramref name="currentMz"/> is <see cref="float.MaxValue"/> all remaining
	/// entries are yielded regardless of m/z.
	///
	/// Entries are removed from <paramref name="buffer"/> as they are yielded.
	///
	/// Output order within a single flush is deterministic: entries are sorted by
	/// PrecursorMz ascending, then Name lexicographic (ordinal), so the
	/// output file is stable across .NET versions, platforms, and insertion orders.
	/// </summary>
	private static IEnumerable<MslLibraryEntry> FlushQValueBuffer(
		Dictionary<string, (float QValue, MslLibraryEntry Entry)> buffer,
		float currentMz,
		float mzWindow)
	{
		// Collect keys eligible for flushing this pass.
		var toFlush = new List<string>();

		foreach (var (key, (_, entry)) in buffer)
		{
			if (currentMz == float.MaxValue || (currentMz - (float)entry.PrecursorMz) > mzWindow)
				toFlush.Add(key);
		}

		// Sort by PrecursorMz ascending, then Name lexicographic, to produce
		// a stable deterministic output order regardless of Dictionary iteration order.
		// This is the only change from the original: without this sort the relative
		// ordering of entries flushed in the same pass was non-deterministic.
		toFlush.Sort((a, b) =>
		{
			float mzA = (float)buffer[a].Entry.PrecursorMz;
			float mzB = (float)buffer[b].Entry.PrecursorMz;
			int mzCmp = ((float)mzA).CompareTo((float)mzB);
			return mzCmp != 0 ? mzCmp : string.Compare(a, b, StringComparison.Ordinal);
		});

		foreach (string key in toFlush)
		{
			MslLibraryEntry winner = buffer[key].Entry;
			buffer.Remove(key);
			yield return winner;
		}
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