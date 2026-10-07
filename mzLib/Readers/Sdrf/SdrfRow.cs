namespace Readers
{
    /// <summary>
    /// One row of an SDRF-Proteomics document: one sample linked to one data file.
    ///
    /// Cells are held as RAW STRINGS and are never interpreted at read time. That is the whole
    /// reason a byte-identical round trip is achievable, and it is also a correctness requirement
    /// rather than laziness: the SDRF key=value grammar ("NT=...;AC=...") is not safely detectable
    /// from a cell's shape. In the curated corpus, comment[file uri] cells carry pre-signed
    /// download URLs whose query strings contain "Signature=", "Expires=", and "Id=" -- 5,750
    /// occurrences each. A parser that treated any cell containing '=' and ';' as a controlled
    /// vocabulary term would happily decode those as CV descriptors.
    ///
    /// Interpretation is therefore an opt-in projection layered on top (see the CV work in
    /// PR-4), not something baked into the document model.
    /// </summary>
    public sealed class SdrfRow
    {
        private readonly List<string> _cells;

        public SdrfRow(SdrfHeader header, IEnumerable<string> cells)
        {
            Header = header ?? throw new ArgumentNullException(nameof(header));
            _cells = cells?.ToList() ?? throw new ArgumentNullException(nameof(cells));
        }

        /// <summary>The document header these cells are positioned against. Shared, not copied.</summary>
        public SdrfHeader Header { get; }

        /// <summary>
        /// The cells exactly as they appeared, in order.
        ///
        /// The count is NOT guaranteed to equal <see cref="SdrfHeader.Count"/>. Real files are
        /// ragged: PXD059974 in the curated corpus has a 46-column header with 17 of its 23 rows
        /// carrying only 42 cells. The reader preserves that rather than padding, so the file
        /// round-trips; reporting it is the validator's job.
        /// </summary>
        public IReadOnlyList<string> Cells => _cells;

        /// <summary>
        /// The value under the first column with this name, or null when the column is absent or
        /// the row is too short to reach it. Null means "not present in this document"; the SDRF
        /// reserved word "not available" is a real, distinct value and is returned as such.
        /// </summary>
        public string? this[string columnName]
        {
            get
            {
                int i = Header.IndexOf(columnName);
                return i >= 0 && i < _cells.Count ? _cells[i] : null;
            }
        }

        /// <summary>
        /// Every value under this name, in document order -- the accessor for multi-cardinality
        /// columns such as comment[modification parameters], which legitimately repeats.
        /// Indexes past the end of a short row are skipped.
        /// </summary>
        public IReadOnlyList<string> All(string columnName)
        {
            var values = new List<string>();
            foreach (int i in Header.IndexesOf(columnName))
            {
                if (i < _cells.Count)
                    values.Add(_cells[i]);
            }
            return values;
        }

        public override string ToString() => string.Join('\t', _cells);
    }
}
