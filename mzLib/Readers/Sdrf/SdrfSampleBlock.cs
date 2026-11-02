using System;
using System.Collections.Generic;
using System.Linq;

namespace Readers
{
    /// <summary>
    /// The sample half of one <c>source name</c>, lifted out of an SDRF somebody else wrote, so that
    /// a search can copy it into the SDRF it writes rather than re-stating facts it does not know.
    ///
    /// This is the read-side half of REQ-2 (sample carry-through). A MetaMorpheus search knows its
    /// assay -- enzyme, tolerances, instrument, modifications -- and knows nothing at all about the
    /// biology: organism part, disease, age, sex, individual. Those live in the deposited SDRF, or in
    /// one aging built from a PRIDE project record, and today they are dropped on the floor. A block
    /// is what survives that journey: every sample cell for one sample, keyed and merged, ready to be
    /// handed back to <see cref="SdrfBuilder"/> through <c>SdrfSample.RawCharacteristics</c> and
    /// <c>FactorValues</c>.
    ///
    /// <para><b>Named cross-assembly caller (D19):</b> MetaMorpheus's
    /// <c>PostSearchAnalysisTaskSdrf</c>, which reads an input SDRF beside the spectra files and
    /// copies the matched block into each row it writes.</para>
    ///
    /// <para><b>What a block is NOT.</b> It does not interpret a single cell. Values come out exactly
    /// as they went in -- reserved words included, because an SDRF built from a config design is
    /// allowed to carry <c>not available</c> (D27) and rewriting that would erase the one honest
    /// statement in the file. Nothing here resolves a CV term, normalises an age, or decides whether
    /// a block is worth having; <see cref="SdrfAge"/> and <see cref="SdrfSampleInformativeness"/>
    /// exist for those and take their input from the same document.</para>
    /// </summary>
    public sealed class SdrfSampleBlock
    {
        private const string SourceNameColumn = "source name";
        private const string CharacteristicsPrefix = "characteristics[";
        private const string FactorValuePrefix = "factor value[";

        /// <summary>
        /// The keys of <see cref="Cells"/> in header order. Kept separately because a dictionary's
        /// enumeration order is not part of its contract.
        /// </summary>
        private readonly IReadOnlyList<string> _columnOrder;

        internal SdrfSampleBlock(string sourceName,
            IReadOnlyDictionary<string, IReadOnlyList<string>> cells,
            IReadOnlyList<string> columnOrder,
            IReadOnlyList<string> conflictingColumns,
            int rowCount)
        {
            SourceName = sourceName;
            Cells = cells;
            _columnOrder = columnOrder;
            ConflictingColumns = conflictingColumns;
            RowCount = rowCount;
        }

        /// <summary>
        /// The sample's name as the FIRST row that named it spelled it, trimmed of surrounding
        /// whitespace because this is the join key. <c>Cells["source name"]</c> keeps the cell
        /// verbatim, padding and all, like every other cell in the block.
        /// </summary>
        public string SourceName { get; }

        /// <summary>
        /// Every sample column this block agrees on, keyed by the column name AS THE HEADER SPELLS IT.
        ///
        /// The value is a list because SDRF columns may repeat -- nine corpus files repeat
        /// <c>characteristics[organism part]</c> -- and joining repeats into one string would invent a
        /// value that no row contains.
        ///
        /// Keys are compared ORDINALLY and case-sensitively, so a file writing both
        /// <c>characteristics[age]</c> and <c>characteristics[Age]</c> yields two entries. That is the
        /// specification's rule for column names and it is deliberate here: unifying spellings is a
        /// curation decision, and a reader that quietly merged them would make it silently and
        /// irreversibly.
        /// </summary>
        public IReadOnlyDictionary<string, IReadOnlyList<string>> Cells { get; }

        /// <summary>
        /// Sample columns whose rows DISAGREE about this sample, excluded from <see cref="Cells"/>.
        ///
        /// Copying one of the disagreeing values through would be the "first row wins" rule, which is
        /// how a sample silently acquires the wrong organism part when a document is inconsistent.
        /// Nothing about the conflict is guessed at: the column is withheld and named, and a caller
        /// that wants to resolve it has the document.
        /// </summary>
        public IReadOnlyList<string> ConflictingColumns { get; }

        /// <summary>How many rows of the document carry this source name.</summary>
        public int RowCount { get; }

        /// <summary>
        /// The first value of a column, or null when the block has no agreed value for it -- which
        /// covers both "the document has no such column" and "its rows disagree". Mirrors
        /// <see cref="SdrfRow.this[string]"/>.
        /// </summary>
        public string? this[string columnName] =>
            columnName is not null && Cells.TryGetValue(columnName, out var values) && values.Count > 0
                ? values[0]
                : null;

        /// <summary>Every value of a repeating column, or empty. Mirrors <see cref="SdrfRow.All"/>.</summary>
        public IReadOnlyList<string> All(string columnName) =>
            columnName is not null && Cells.TryGetValue(columnName, out var values)
                ? values
                : Array.Empty<string>();

        /// <summary>The <c>characteristics[...]</c> columns this block agrees on, in header order.</summary>
        public IEnumerable<string> CharacteristicColumns =>
            _columnOrder.Where(c => c.StartsWith(CharacteristicsPrefix, StringComparison.OrdinalIgnoreCase));

        /// <summary>The <c>factor value[...]</c> columns this block agrees on, in header order.</summary>
        public IEnumerable<string> FactorValueColumns =>
            _columnOrder.Where(c => c.StartsWith(FactorValuePrefix, StringComparison.OrdinalIgnoreCase));

        public override string ToString() =>
            $"{SourceName}: {Cells.Count} column(s) over {RowCount} row(s)" +
            (ConflictingColumns.Count == 0 ? "" : $", {ConflictingColumns.Count} in conflict");

        /// <summary>
        /// Every sample block in a document, keyed by <c>source name</c>.
        ///
        /// <para><b>Keyed on source name ONLY</b> (D27). Not on source name plus data file, not on the
        /// assay: one sample measured in three runs is one sample, and it is the join that lets a
        /// search match a row it is writing to a sample somebody else described.</para>
        ///
        /// <para><b>The two matching rules are different rules, and both are applied here.</b> Column
        /// NAMES are ordinal and case-sensitive, per the specification and
        /// <see cref="SdrfHeader.IndexOf"/>; cell VALUES are compared case-insensitively, because
        /// <c>PXD023158</c> writes <c>tmt126</c> and means the channel every other file spells
        /// <c>TMT126</c>. So <c>Sample1</c> and <c>sample1</c> are ONE block, while
        /// <c>characteristics[age]</c> and <c>characteristics[Age]</c> are TWO columns.</para>
        ///
        /// <para><b>The sample half is chosen by column NAME, not by position.</b> It is
        /// <c>source name</c>, every <c>characteristics[...]</c> and every <c>factor value[...]</c>.
        /// Position would have been the specification's answer -- the sample columns run before
        /// <c>assay name</c> -- but an extension column smuggled in through the wrong door lands on the
        /// wrong side of it, and this way a ragged or reordered document still yields the right cells.
        /// The prefix match ignores case because two corpus files write <c>Factor Value[</c>, and
        /// missing a factor column loses the one cell that says what the study varied.</para>
        /// </summary>
        /// <param name="document">The document to read. Its results are loaded if they are not yet.</param>
        /// <param name="problems">
        /// Rows that could not be placed: a document with no <c>source name</c> column, or a row whose
        /// source name is blank. Never null, and empty on a well-formed document.
        /// </param>
        public static IReadOnlyDictionary<string, SdrfSampleBlock> BySourceName(
            SdrfDocument document, out IReadOnlyList<string> problems)
        {
            if (document is null) throw new ArgumentNullException(nameof(document));

            var found = new List<string>();
            problems = found;
            var blocks = new Dictionary<string, SdrfSampleBlock>(StringComparer.OrdinalIgnoreCase);

            SdrfHeader header = document.Header;
            string? sourceColumn = header
                .FirstOrDefault(n => string.Equals(n, SourceNameColumn, StringComparison.OrdinalIgnoreCase));
            if (sourceColumn is null)
            {
                found.Add($"The document has no {SourceNameColumn} column, so it describes no samples.");
                return blocks;
            }

            // The header's OWN spelling of each sample column, de-duplicated but order-preserving.
            // Using the header's spelling rather than a canonical one is what makes SdrfRow.All --
            // which is ordinal by design -- return the cells of a mis-cased column instead of nothing.
            var sampleColumns = header
                .Where(IsSampleColumn)
                .Distinct(StringComparer.Ordinal)
                .ToList();

            var rowsBySample = new Dictionary<string, List<SdrfRow>>(StringComparer.OrdinalIgnoreCase);
            var firstSpelling = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);

            for (int i = 0; i < document.Results.Count; i++)
            {
                SdrfRow row = document.Results[i];
                string name = row[sourceColumn]?.Trim() ?? "";
                if (name.Length == 0)
                {
                    // Row 1 is the first DATA row: the header is not a result.
                    found.Add($"Row {i + 1} has a blank {SourceNameColumn} and cannot be keyed to a sample.");
                    continue;
                }

                if (!rowsBySample.TryGetValue(name, out var rows))
                {
                    rowsBySample[name] = rows = new List<SdrfRow>();
                    firstSpelling[name] = name;
                }
                rows.Add(row);
            }

            foreach (var (name, rows) in rowsBySample)
            {
                var agreed = new Dictionary<string, IReadOnlyList<string>>(StringComparer.Ordinal);
                var order = new List<string>();
                var conflicts = new List<string>();

                foreach (string column in sampleColumns)
                {
                    IReadOnlyList<string> first = rows[0].All(column);
                    bool agrees = rows.Skip(1).All(r => SameValues(first, r.All(column)));

                    if (agrees)
                    {
                        agreed[column] = first;
                        order.Add(column);
                    }
                    else conflicts.Add(column);
                }

                blocks[name] = new SdrfSampleBlock(firstSpelling[name], agreed, order, conflicts, rows.Count);
            }

            return blocks;
        }

        private static bool IsSampleColumn(string column) =>
            string.Equals(column, SourceNameColumn, StringComparison.OrdinalIgnoreCase)
            || column.StartsWith(CharacteristicsPrefix, StringComparison.OrdinalIgnoreCase)
            || column.StartsWith(FactorValuePrefix, StringComparison.OrdinalIgnoreCase);

        /// <summary>
        /// Whether two rows say the same thing in one column. Values are compared case-insensitively
        /// and without surrounding whitespace -- the cell rule -- but the values KEPT are whichever
        /// the first row wrote, verbatim, because a block is a copy and not a normalisation.
        /// </summary>
        private static bool SameValues(IReadOnlyList<string> left, IReadOnlyList<string> right)
        {
            if (left.Count != right.Count) return false;
            for (int i = 0; i < left.Count; i++)
                if (!string.Equals(left[i]?.Trim(), right[i]?.Trim(), StringComparison.OrdinalIgnoreCase))
                    return false;
            return true;
        }
    }
}
