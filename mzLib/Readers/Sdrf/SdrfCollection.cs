using System.Collections;

namespace Readers
{
    /// <summary>
    /// Several SDRF documents treated as one body of experiments.
    ///
    /// This is the type the whole component exists for: annotating each experiment consistently is
    /// only worth doing if the annotations can then be pooled and queried together. Merging is the
    /// easy half; <see cref="SdrfDriftLint"/> is the half that decides whether the pooled table
    /// means anything, because per-document validity says nothing about whether two documents used
    /// the same word for the same thing.
    /// </summary>
    public sealed class SdrfCollection : IReadOnlyList<SdrfDocument>
    {
        /// <summary>
        /// Column added by <see cref="Merge"/> naming the document each row came from. A comment[...]
        /// column so the merged table is still a well-formed SDRF document.
        /// </summary>
        public const string SourceDocumentColumn = "comment[source document]";

        private readonly List<SdrfDocument> _documents;
        private readonly List<string> _labels;

        /// <param name="documents">The documents to pool.</param>
        /// <param name="labels">
        /// One label per document, used as the provenance value in <see cref="Merge"/>. When null,
        /// each document's file name without extension is used.
        /// </param>
        public SdrfCollection(IEnumerable<SdrfDocument> documents, IEnumerable<string> labels = null)
        {
            _documents = documents?.ToList() ?? throw new ArgumentNullException(nameof(documents));
            _labels = labels?.ToList() ?? _documents.Select(DefaultLabel).ToList();

            if (_labels.Count != _documents.Count)
                throw new ArgumentException(
                    $"Got {_labels.Count} labels for {_documents.Count} documents.", nameof(labels));
        }

        public static SdrfCollection FromFiles(IEnumerable<string> filePaths)
        {
            if (filePaths is null) throw new ArgumentNullException(nameof(filePaths));
            var paths = filePaths.ToList();
            return new SdrfCollection(paths.Select(p => new SdrfDocument(p)));
        }

        public IReadOnlyList<string> Labels => _labels;

        public SdrfDocument this[int index] => _documents[index];

        public int Count => _documents.Count;

        /// <summary>
        /// The union of every document's columns.
        ///
        /// Ordered by SDRF's own block structure -- source name, characteristics, assay name and
        /// technology type, comments, then factor values -- because the union of several documents
        /// has no natural order of its own and an arbitrary one would make the merged table differ
        /// between runs. Within a block, columns keep the order in which they were first seen, so a
        /// single-document collection round-trips its own column order.
        ///
        /// A column that repeats within one document (comment[modification parameters] repeats up
        /// to 8 times in the corpus) is carried over at the highest multiplicity any single document
        /// uses, so no values are dropped.
        /// </summary>
        public SdrfHeader MergedHeader
        {
            get
            {
                var maxMultiplicity = new Dictionary<string, int>(StringComparer.Ordinal);
                var firstSeen = new List<string>();

                foreach (var document in _documents)
                {
                    var counts = new Dictionary<string, int>(StringComparer.Ordinal);
                    foreach (var name in document.Header)
                    {
                        counts.TryGetValue(name, out int c);
                        counts[name] = c + 1;
                        if (!maxMultiplicity.ContainsKey(name))
                        {
                            maxMultiplicity[name] = 0;
                            firstSeen.Add(name);
                        }
                    }
                    foreach (var kv in counts)
                        if (kv.Value > maxMultiplicity[kv.Key])
                            maxMultiplicity[kv.Key] = kv.Value;
                }

                var ordered = firstSeen
                    .Select((name, index) => (name, index))
                    .OrderBy(x => BlockRank(x.name))
                    .ThenBy(x => x.index)
                    .SelectMany(x => Enumerable.Repeat(x.name, maxMultiplicity[x.name]))
                    .ToList();

                // Only append the provenance column if the inputs do not already carry one.
                // Appending unconditionally duplicated the name when a MERGED document was fed back
                // into a collection -- the obvious incremental workflow, and one Merge_IsWritableAsSdrf
                // shows is supported. IndexOf then bound the label to the inherited column,
                // overwriting the original document's provenance, while the appended column filled
                // with "not available" on every row and the header grew by one per merge generation.
                if (!ordered.Contains(SourceDocumentColumn, StringComparer.Ordinal))
                    ordered.Add(SourceDocumentColumn);

                return new SdrfHeader(ordered);
            }
        }

        /// <summary>
        /// Pools every row into one document over <see cref="MergedHeader"/>, with a
        /// comment[source document] column recording provenance. A column a document does not have
        /// is filled with the SDRF reserved word "not available" rather than left blank, so the
        /// merged table says "this experiment did not record it" instead of being silently empty.
        ///
        /// The result is an ANALYSIS table, not something to deposit. source name + assay name +
        /// comment[label] is unique within a document but nothing stops two experiments both having
        /// a "Sample 1"/"run 1", so the merged table will usually fail
        /// <see cref="SdrfValidator"/>'s uniqueness rule. Use comment[source document] as part of
        /// the key when querying.
        /// </summary>
        public SdrfDocument Merge()
        {
            var header = MergedHeader;
            int sourceColumn = header.IndexOf(SourceDocumentColumn);
            var rows = new List<SdrfRow>();

            for (int d = 0; d < _documents.Count; d++)
            {
                var document = _documents[d];
                var sourceHeader = document.Header;

                // Which source indexes feed each merged column, respecting multiplicity: the Nth
                // occurrence of a name in the merged header takes the Nth occurrence in the source.
                var used = new Dictionary<string, int>(StringComparer.Ordinal);
                var map = new int[header.Count];
                for (int c = 0; c < header.Count; c++)
                {
                    string name = header[c];
                    used.TryGetValue(name, out int taken);
                    var candidates = sourceHeader.IndexesOf(name);
                    map[c] = taken < candidates.Count ? candidates[taken] : -1;
                    used[name] = taken + 1;
                }

                foreach (var row in document.Results)
                {
                    var cells = new string[header.Count];
                    for (int c = 0; c < header.Count; c++)
                    {
                        int source = map[c];
                        cells[c] = source >= 0 && source < row.Cells.Count
                            ? row.Cells[source]
                            : "not available";
                    }
                    // Only stamp our label where the source document had none. Overwriting an
                    // inherited value destroyed the per-experiment provenance when a merged
                    // document was merged again -- every row became the outer label, and which
                    // original experiment a row came from was no longer recoverable. Fixing the
                    // duplicate column without fixing this left the data loss in place.
                    if (string.IsNullOrEmpty(cells[sourceColumn])
                        || cells[sourceColumn] == "not available")
                        cells[sourceColumn] = _labels[d];
                    rows.Add(new SdrfRow(header, cells));
                }
            }

            return new SdrfDocument(header, rows);
        }

        public IEnumerator<SdrfDocument> GetEnumerator() => _documents.GetEnumerator();

        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

        /// <summary>
        /// A label that stays unique across documents.
        ///
        /// Deliberately NOT the bare file name. Every search writes its SDRF into its own output
        /// folder, so N re-searches of one experiment produce N files with the SAME base name in
        /// different directories -- labels collided, and Merge then attributed every row to
        /// whichever document was seen last. That is exactly the path-implied identity the
        /// experiment-in-a-column design exists to avoid, so the fallback keeps enough of the path
        /// to stay distinguishable.
        /// </summary>
        private static string DefaultLabel(SdrfDocument document)
        {
            if (string.IsNullOrEmpty(document.FilePath))
                return "(in memory)";

            string stem = Path.GetFileName(document.FilePath);
            string extension = SupportedFileType.Sdrf.GetFileExtension();
            if (stem.EndsWith(extension, StringComparison.OrdinalIgnoreCase))
                stem = stem.Substring(0, stem.Length - extension.Length);

            // Qualify with the containing folder, which is what distinguishes one search's output
            // from another's. Callers that have a better identity should pass explicit labels.
            string? folder = Path.GetFileName(Path.GetDirectoryName(document.FilePath));
            return string.IsNullOrEmpty(folder) ? stem : $"{folder}/{stem}";
        }

        /// <summary>
        /// SDRF's block order: sample metadata, then data-file metadata, then factor values.
        /// </summary>
        private static int BlockRank(string columnName)
        {
            if (string.Equals(columnName, "source name", StringComparison.OrdinalIgnoreCase)) return 0;
            if (columnName.StartsWith("characteristics[", StringComparison.OrdinalIgnoreCase)) return 1;
            if (string.Equals(columnName, "assay name", StringComparison.OrdinalIgnoreCase)) return 2;
            if (string.Equals(columnName, "technology type", StringComparison.OrdinalIgnoreCase)) return 3;
            if (columnName.StartsWith("comment[", StringComparison.OrdinalIgnoreCase)) return 4;
            // Case-INSENSITIVE, matching SdrfValidator.BlockRank. An ordinal test put
            // "Factor Value[organism part]" (which the corpus contains) in the unrecognised
            // bucket, so a merged header ordered it before the comments and the validator then
            // reported ColumnOrdering on this tool's own output.
            if (columnName.StartsWith("factor value[", StringComparison.OrdinalIgnoreCase)) return 6;
            return 5; // unrecognised columns sit with the comments, before factor values
        }
    }
}
