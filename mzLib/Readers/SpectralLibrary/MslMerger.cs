using Omics.SpectralMatch.MslSpectralLibrary;

namespace Readers.SpectralLibrary;

// ────────────────────────────────────────────────────────────────────────────
// MslMergeConflictPolicy
// ────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Determines which entry is kept when the same LookupKey appears in multiple source files.
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
	/// Number of entries skipped due to duplicate LookupKey resolution.
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
/// streaming k-way merge algorithm. Memory usage is O(k + UniqueStrings), not
/// O(TotalEntries), making this suitable for merging very large libraries.
///
/// <para>
/// Source files are read via <see cref="MslLibrary.LoadIndexOnly"/> and enumerated with
/// <see cref="MslLibrary.GetAllEntries"/>. The output is written via
/// <see cref="MslWriter.WriteStreaming"/>. All source file handles are released when the
/// method returns.
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
	///   When <see langword="true"/> (default), each unique LookupKey is written only once.
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

		// ── KeepLast two-pass pre-scan ────────────────────────────────────────
		// For KeepLast we need to know which source index last contains each key.
		Dictionary<string, int>? keepLastMap = null;
		if (deduplicate && conflictPolicy == MslMergeConflictPolicy.KeepLast)
		{
			keepLastMap = BuildKeepLastMap(sourcePaths);
		}

		// ── Open source libraries ─────────────────────────────────────────────
		MslLibrary[] libraries = new MslLibrary[k];
		IEnumerator<MslLibraryEntry>[] enumerators = new IEnumerator<MslLibraryEntry>[k];
		int[] sourceCounts = new int[k];
		float[] lastMz = new float[k]; // for unsorted-source detection
		var unsortedSourceFiles = new List<string>();

		try
		{
			for (int i = 0; i < k; i++)
			{
				libraries[i] = MslLibrary.LoadIndexOnly(sourcePaths[i]);
				enumerators[i] = libraries[i].GetAllEntries().GetEnumerator();
				lastMz[i] = float.MinValue;
			}

			// ── k-way min-heap: element = (entry, enumerator, sourceIdx), priority = PrecursorMz
			var heap = new PriorityQueue<(MslLibraryEntry Entry, IEnumerator<MslLibraryEntry> Enum, int SourceIdx), float>();

			// Seed the heap
			for (int i = 0; i < k; i++)
			{
				if (TryAdvance(enumerators[i], sourcePaths[i], ref lastMz[i], unsortedSourceFiles, ref sourceCounts[i], out MslLibraryEntry? first))
					heap.Enqueue((first!, enumerators[i], i), (float)first!.PrecursorMz);
			}

			// ── Deduplication state ───────────────────────────────────────────
			// KeepFirst: HashSet of already-emitted keys
			var emittedKeys = (deduplicate && conflictPolicy == MslMergeConflictPolicy.KeepFirst)
				? new HashSet<string>(StringComparer.Ordinal)
				: null;

			// KeepLowestQValue: buffer pending entries within a narrow m/z window
			// Key = LookupKey, Value = (bestQValue, entry)
			var qValueBuffer = (deduplicate && conflictPolicy == MslMergeConflictPolicy.KeepLowestQValue)
				? new Dictionary<string, (float QValue, MslLibraryEntry Entry)>(StringComparer.Ordinal)
				: null;

			int totalRead = 0;
			int duplicatesSkipped = 0;
			var outputEntries = new List<MslLibraryEntry>();

			while (heap.Count > 0)
			{
				heap.TryDequeue(out var item, out float _);
				var (entry, enumerator, sourceIdx) = item;

				totalRead++;

				// Advance this source and push its next entry onto the heap
				if (TryAdvance(enumerator, sourcePaths[sourceIdx], ref lastMz[sourceIdx], unsortedSourceFiles, ref sourceCounts[sourceIdx], out MslLibraryEntry? next))
					heap.Enqueue((next!, enumerator, sourceIdx), (float)next!.PrecursorMz);

				if (!deduplicate)
				{
					outputEntries.Add(entry);
					continue;
				}

				string key = entry.LookupKey;

				switch (conflictPolicy)
				{
					case MslMergeConflictPolicy.KeepFirst:
						if (emittedKeys!.Add(key))
							outputEntries.Add(entry);
						else
							duplicatesSkipped++;
						break;

					case MslMergeConflictPolicy.KeepLast:
						// Only emit if this source is the last source that has this key
						if (keepLastMap!.TryGetValue(key, out int lastSource) && lastSource == sourceIdx)
							outputEntries.Add(entry);
						else
							duplicatesSkipped++;
						break;

					case MslMergeConflictPolicy.KeepLowestQValue:
						// Flush buffer entries that are far enough behind the current m/z
						FlushQValueBuffer(qValueBuffer!, outputEntries, (float)entry.PrecursorMz, 0.001f, ref duplicatesSkipped);

						// Update buffer with the better (lower) q-value entry
						float q = entry.QValue;
						if (qValueBuffer!.TryGetValue(key, out var existing))
						{
							// Keep lower q-value; NaN is treated as infinity (worst)
							float existingQ = existing.QValue;
							bool replaceWithNew = (!float.IsNaN(q) && (float.IsNaN(existingQ) || q < existingQ));
							if (replaceWithNew)
								qValueBuffer[key] = (q, entry);
							// Either way, one entry for this key is a duplicate
							duplicatesSkipped++;
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
				FlushQValueBuffer(qValueBuffer, outputEntries, float.MaxValue, 0f, ref duplicatesSkipped);

			// Recalculate totalRead from per-source counts (loop above counts heap pops which
			// match entries read, but sourceCounts were incremented in TryAdvance on MoveNext)
			int totalSourceCount = 0;
			for (int i = 0; i < k; i++) totalSourceCount += sourceCounts[i];

			// Write the merged output atomically via WriteStreaming
			MslWriter.WriteStreaming(outputPath, outputEntries, compressionLevel);

			return new MslMergeResult
			{
				OutputPath = outputPath,
				OutputEntryCount = outputEntries.Count,
				TotalSourceEntryCount = totalSourceCount,
				DuplicatesSkipped = duplicatesSkipped,
				UnsortedSourceFiles = unsortedSourceFiles.Distinct().ToList(),
				SourceEntryCounts = sourceCounts
			};
		}
		finally
		{
			for (int i = 0; i < k; i++)
			{
				try { enumerators[i]?.Dispose(); } catch { /* best-effort */ }
				try { libraries[i]?.Dispose(); } catch { /* best-effort */ }
			}
		}
	}

	// ────────────────────────────────────────────────────────────────────────
	// Private helpers
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Builds a dictionary mapping each LookupKey to the index of the last source file
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
				map[entry.LookupKey] = i;   // overwrite: last writer wins
		}

		return map;
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
	/// Flushes all KeepLowestQValue buffer entries whose m/z is more than
	/// <paramref name="mzWindow"/> Da below <paramref name="currentMz"/> into
	/// <paramref name="output"/>. When <paramref name="currentMz"/> is
	/// <see cref="float.MaxValue"/> all remaining entries are flushed regardless of m/z.
	/// </summary>
	private static void FlushQValueBuffer(
		Dictionary<string, (float QValue, MslLibraryEntry Entry)> buffer,
		List<MslLibraryEntry> output,
		float currentMz,
		float mzWindow,
		ref int duplicatesSkipped)
	{
		// Collect keys to flush in this pass
		var toFlush = new List<string>();

		foreach (var (key, (_, entry)) in buffer)
		{
			if (currentMz == float.MaxValue || (currentMz - (float)entry.PrecursorMz) > mzWindow)
				toFlush.Add(key);
		}

		foreach (string key in toFlush)
		{
			output.Add(buffer[key].Entry);
			buffer.Remove(key);
		}
	}
}