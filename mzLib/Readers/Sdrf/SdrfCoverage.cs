namespace Readers
{
    /// <summary>
    /// How completely one column is actually filled in across a collection.
    /// </summary>
    /// <param name="Column">The column name.</param>
    /// <param name="Documents">Documents that declare the column at all.</param>
    /// <param name="Rows">Rows in those documents.</param>
    /// <param name="Filled">Rows whose value is a real answer.</param>
    /// <param name="Absent">
    /// Rows that say nothing: empty, or one of the SDRF reserved words. Reserved words are the
    /// SPEC-CORRECT way to write "no value", which is exactly why they have to be counted as
    /// unfilled here -- a validator sees them as right, and they are, but a corpus made of them
    /// answers no question.
    /// </param>
    /// <param name="DistinctValues">
    /// How many different real answers appear. A column filled to 100% with one repeated value
    /// cannot group anything, so fill rate alone is not enough.
    /// </param>
    public sealed record SdrfColumnCoverage(
        string Column,
        int Documents,
        int Rows,
        int Filled,
        int Absent,
        int DistinctValues)
    {
        /// <summary>Fraction of rows carrying a real answer, 0 to 1. Zero when the column has no rows.</summary>
        public double FillRate => Rows == 0 ? 0d : (double)Filled / Rows;

        /// <summary>
        /// True when the column is present but says nothing useful -- either almost never filled, or
        /// filled with a single value everywhere. Both make it useless for grouping across
        /// experiments, which is the only reason to have annotated it.
        /// </summary>
        public bool IsUninformative => FillRate < 0.5d || DistinctValues <= 1;

        public override string ToString() =>
            $"{Column}: {Filled}/{Rows} filled ({FillRate:P0}), {DistinctValues} distinct value(s)" +
            (IsUninformative ? "  <-- uninformative" : "");
    }

    /// <summary>
    /// Measures how much a pooled set of SDRF documents actually SAYS, as opposed to how consistent
    /// or how valid it is.
    ///
    /// This exists because the other two instruments are structurally blind to an empty corpus.
    /// <see cref="SdrfValidator"/> treats the reserved words -- "not available", "not applicable" --
    /// as the correct way to state an absence, and it is right to: they are the spec's own
    /// mechanism. <see cref="SdrfDriftLint"/> skips them outright, since a reserved word cannot
    /// drift. So a corpus whose sample columns are 90% "not available" produces zero validator
    /// errors and zero drift findings, and reads as flawless.
    ///
    /// That is a real risk rather than a hypothetical one. The assay columns -- enzyme,
    /// modifications, tolerances, instrument -- come from a search's own parameters and will always
    /// be complete. The sample columns -- organism part, disease, cell type, replicate structure --
    /// come from a human, and nothing in a search knows them. An implementation that quietly writes
    /// a reserved word when it cannot find one produces a large, uniform, perfectly consistent
    /// corpus that cannot answer a single biological question, because grouping across experiments
    /// means grouping by sample properties, not by our own search settings.
    ///
    /// Fill rate is the metric that makes that visible, and distinct-value count is what stops a
    /// column filled with one repeated value from looking healthy.
    /// </summary>
    public static class SdrfCoverage
    {
        /// <summary>
        /// Per-column coverage across the collection, ordered worst-first so the empty columns are
        /// the ones you see.
        /// </summary>
        public static IReadOnlyList<SdrfColumnCoverage> Measure(SdrfCollection collection)
        {
            if (collection is null) throw new ArgumentNullException(nameof(collection));

            var documents = new Dictionary<string, HashSet<int>>(StringComparer.Ordinal);
            var rows = new Dictionary<string, int>(StringComparer.Ordinal);
            var filled = new Dictionary<string, int>(StringComparer.Ordinal);
            var values = new Dictionary<string, HashSet<string>>(StringComparer.Ordinal);

            for (int d = 0; d < collection.Count; d++)
            {
                var document = collection[d];
                var header = document.Header;

                foreach (var row in document.Results)
                {
                    for (int c = 0; c < header.Count; c++)
                    {
                        string column = header[c];
                        if (string.IsNullOrWhiteSpace(column)) continue;

                        if (!documents.TryGetValue(column, out var seenIn))
                        {
                            documents[column] = seenIn = new HashSet<int>();
                            rows[column] = 0;
                            filled[column] = 0;
                            values[column] = new HashSet<string>(StringComparer.Ordinal);
                        }
                        seenIn.Add(d);
                        rows[column]++;

                        string cell = c < row.Cells.Count ? row.Cells[c] : "";
                        if (IsAnswer(cell))
                        {
                            filled[column]++;
                            values[column].Add(cell);
                        }
                    }
                }
            }

            return rows.Keys
                .Select(column => new SdrfColumnCoverage(
                    column,
                    documents[column].Count,
                    rows[column],
                    filled[column],
                    rows[column] - filled[column],
                    values[column].Count))
                .OrderBy(c => c.FillRate)
                .ThenBy(c => c.DistinctValues)
                .ThenBy(c => c.Column, StringComparer.Ordinal)
                .ToList();
        }

        /// <summary>
        /// The columns that are present but say nothing -- the ones a validator and a drift lint
        /// both report as fine. This is the list worth failing a build on.
        /// </summary>
        public static IReadOnlyList<SdrfColumnCoverage> Uninformative(SdrfCollection collection) =>
            Measure(collection).Where(c => c.IsUninformative).ToList();

        /// <summary>
        /// A cell counts as an answer when it is neither empty nor a reserved word. Comparison is
        /// case-insensitive on purpose: "Not Available" says just as little as "not available",
        /// even though the validator separately -- and correctly -- reports the casing.
        /// </summary>
        private static bool IsAnswer(string cell) =>
            !string.IsNullOrWhiteSpace(cell)
            && !SdrfValidator.ReservedWords.Any(w => string.Equals(cell.Trim(), w, StringComparison.OrdinalIgnoreCase));
    }
}
