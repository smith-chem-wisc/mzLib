using System.Text.RegularExpressions;

namespace Readers
{
    /// <summary>
    /// Structural validation of an SDRF-Proteomics document against specification v1.1.0.
    ///
    /// Deliberately separate from <see cref="SdrfDocument.LoadResults"/>: a malformed-but-readable
    /// file must still parse, so that a caller can read it, be told what is wrong with it, and fix
    /// it. A reader that refused to load an imperfect file could not round-trip one, and a reader
    /// that repaired one would silently rewrite a curated annotation.
    ///
    /// Scope is STRUCTURE ONLY -- column presence, ordering, uniqueness, cell shape. It does not
    /// resolve controlled-vocabulary terms; "is MS:1001911 a real accession, and is it the term it
    /// claims to be" needs the ontology snapshot and arrives with the CV resolver.
    ///
    /// Rule severities were calibrated against the 1,236-file curated corpus rather than read off
    /// the specification. Where the two disagreed the corpus won, because a rule that fires on most
    /// curated files is a wrong rule. See SdrfValidatorCorpusCalibration in the test project, which
    /// pins the observed rates so that a future rule change cannot quietly start condemning the
    /// reference set.
    /// </summary>
    public static class SdrfValidator
    {
        /// <summary>
        /// Columns present in EVERY one of the 1,236 curated corpus files, so demanding them
        /// condemns nothing the community considers valid — which is the whole calibration rule.
        ///
        /// There are thirteen, not the three an earlier revision listed. That mistake came from
        /// reading the survey's OCCURRENCE counts as FILE counts: a column that repeats inflates its
        /// occurrences, so comment[modification parameters] shows 2,564 occurrences while sitting in
        /// only 1,176 files, and characteristics[organism] shows 1,242 while being genuinely
        /// universal. Ten columns that mining actually depends on were being reported as warnings.
        /// </summary>
        public static readonly IReadOnlyList<string> RequiredColumns = new[]
        {
            "source name",
            "assay name",
            "technology type",
            "characteristics[organism]",
            "characteristics[organism part]",
            "characteristics[biological replicate]",
            "comment[data file]",
            "comment[instrument]",
            "comment[label]",
            "comment[cleavage agent details]",
            "comment[technical replicate]",
            "comment[fraction identifier]",
            "comment[proteomics data acquisition method]"
        };

        /// <summary>
        /// Columns the specification lists as mandatory but which are genuinely absent from part of
        /// the curated corpus, so requiring them would invalidate curated data. Warnings.
        /// characteristics[disease] misses 4 files, characteristics[cell type] 15,
        /// comment[modification parameters] 60.
        /// </summary>
        public static readonly IReadOnlyList<string> RecommendedColumns = new[]
        {
            "characteristics[disease]",
            "characteristics[cell type]",
            "comment[modification parameters]"
        };

        /// <summary>
        /// Reserved words. The specification requires them lowercase; a cased variant is the single
        /// most common way a hand-edited file drifts, and it breaks equality-based mining.
        /// </summary>
        public static readonly IReadOnlyList<string> ReservedWords = new[]
        {
            "not available",
            "not applicable",
            "anonymized",
            "pooled"
        };

        /// <summary>
        /// Columns whose value must be a positive integer or a reserved word. Getting these wrong
        /// silently corrupts any grouping done downstream.
        /// </summary>
        private static readonly string[] IntegerColumns =
        {
            "characteristics[biological replicate]",
            "comment[technical replicate]",
            "comment[fraction identifier]"
        };

        // "characteristics [organism]" (space before the bracket) and "Characteristics[organism]"
        // both read as a different column under the spec's case- and space-sensitive comparison,
        // so they silently become an unrecognized extra column rather than an error.
        private static readonly Regex StructuredColumn =
            new(@"^(characteristics|comment|factor value)\s*\[", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        private static readonly Regex WellFormedStructuredColumn =
            new(@"^(characteristics|comment|factor value)\[[^\]]+\]$", RegexOptions.Compiled);

        public static SdrfValidationResult Validate(SdrfDocument document)
        {
            if (document is null) throw new ArgumentNullException(nameof(document));

            var messages = new List<SdrfValidationMessage>();
            var header = document.Header;
            var rows = document.Results;

            ValidateHeader(header, messages);
            ValidateRows(header, rows, messages);
            ValidateRowKeyUniqueness(header, rows, messages);

            return new SdrfValidationResult(messages);
        }

        private static void ValidateHeader(SdrfHeader header, List<SdrfValidationMessage> messages)
        {
            if (header.Count == 0)
            {
                messages.Add(new SdrfValidationMessage(
                    SdrfValidationSeverity.Error, "EmptyHeader", "The document has no columns."));
                return;
            }

            foreach (var required in RequiredColumns)
            {
                if (!header.Contains(required))
                    messages.Add(new SdrfValidationMessage(
                        SdrfValidationSeverity.Error, "RequiredColumn",
                        $"Required column '{required}' is missing.", null, required));
            }

            foreach (var recommended in RecommendedColumns)
            {
                if (!header.Contains(recommended))
                    messages.Add(new SdrfValidationMessage(
                        SdrfValidationSeverity.Warning, "RecommendedColumn",
                        $"Column '{recommended}' is listed as mandatory for MS proteomics but is absent.",
                        null, recommended));
            }

            for (int i = 0; i < header.Count; i++)
            {
                string name = header[i];

                if (string.IsNullOrWhiteSpace(name))
                {
                    // Almost always trailing tabs on the header line. One corpus file has 23.
                    messages.Add(new SdrfValidationMessage(
                        SdrfValidationSeverity.Warning, "EmptyColumnName",
                        $"Column {i} has an empty name; this is usually a trailing tab.", null, name));
                    continue;
                }

                if (name != name.ToLowerInvariant())
                    messages.Add(new SdrfValidationMessage(
                        SdrfValidationSeverity.Warning, "ColumnNameCase",
                        $"Column name '{name}' is not lowercase. Column names are case-sensitive, so " +
                        "this will not match the same column in another document.", null, name));

                // A structured name that is not well-formed is worse than a cased one: it does not
                // parse as characteristics/comment/factor value at all, so consumers see an
                // unrecognized free-text column instead of the property it was meant to be.
                if (StructuredColumn.IsMatch(name) && !WellFormedStructuredColumn.IsMatch(name))
                    messages.Add(new SdrfValidationMessage(
                        SdrfValidationSeverity.Warning, "MalformedColumnName",
                        $"Column name '{name}' looks like a structured column but is not of the form " +
                        "'characteristics[...]', 'comment[...]' or 'factor value[...]' -- check for a " +
                        "space before the bracket.", null, name));
            }

            ValidateColumnOrdering(header, messages);
        }

        /// <summary>
        /// The specification fixes the block order: sample metadata (source name, characteristics)
        /// first, then data-file metadata (assay name, comment), then factor values. Order is
        /// checked by block, not by exact position, because the set of columns within a block is
        /// open-ended.
        /// </summary>
        private static void ValidateColumnOrdering(SdrfHeader header, List<SdrfValidationMessage> messages)
        {
            // Rewritten: the original compared only two pairs of extremes (last characteristic vs
            // first comment, last non-factor vs first factor) and consequently fired ZERO times
            // across all 1,236 curated files while missing the violation its own summary describes.
            // ["source name","assay name","technology type","characteristics[organism]",
            //  "comment[data file]"] has sample metadata after data-file metadata and produced no
            // message, because lastCharacteristic(3) > firstComment(4) is false. It also never
            // checked the position of source name, assay name or technology type at all.
            //
            // Instead: give every column its block rank and require the sequence to be
            // non-decreasing. That detects every out-of-place column, including a leading
            // characteristics[...] before source name.
            int previousRank = int.MinValue;
            int previousIndex = -1;

            for (int i = 0; i < header.Count; i++)
            {
                string name = header[i];
                if (string.IsNullOrWhiteSpace(name)) continue;

                int? maybeRank = BlockRank(name);
                if (maybeRank is null)
                    continue; // unrecognised column: the specification gives it no position

                int rank = maybeRank.Value;
                if (rank < previousRank)
                {
                    // Report the column that is out of place, then ADOPT its rank. Holding the old
                    // high-water mark produced a cascade: one misplaced factor value[...] made every
                    // following column look wrong, so a single defect emitted three messages, none
                    // of which named the actual offender.
                    messages.Add(new SdrfValidationMessage(
                        SdrfValidationSeverity.Warning, "ColumnOrdering",
                        $"'{header[previousIndex]}' (column {previousIndex}) appears before " +
                        $"'{name}' (column {i}) but belongs after it. The specification orders " +
                        "columns as source name, then characteristics[...], then assay name and " +
                        "technology type, then comment[...], then factor value[...].",
                        null, header[previousIndex]));
                }

                previousRank = rank;
                previousIndex = i;
            }
        }

        /// <summary>
        /// The specification's block order. Compared case-INSENSITIVELY on purpose: the corpus
        /// contains "Factor Value[organism part]", and an ordinal test silently disarmed the whole
        /// rule on those files while a mixed-case pair elsewhere produced a false positive. Casing
        /// is reported separately by ColumnNameCase, which is where it belongs.
        /// </summary>
        private static int? BlockRank(string columnName)
        {
            if (string.Equals(columnName, "source name", StringComparison.OrdinalIgnoreCase)) return 0;
            if (columnName.StartsWith("characteristics[", StringComparison.OrdinalIgnoreCase)) return 1;
            if (string.Equals(columnName, "assay name", StringComparison.OrdinalIgnoreCase)) return 2;
            if (string.Equals(columnName, "technology type", StringComparison.OrdinalIgnoreCase)) return 3;
            if (columnName.StartsWith("comment[", StringComparison.OrdinalIgnoreCase)) return 4;
            if (columnName.StartsWith("factor value[", StringComparison.OrdinalIgnoreCase)) return 5;

            // Null, NOT a rank after the comment block. Ranking unknowns forced legitimate
            // sample-section columns to the end and produced a false positive on 28 curated files:
            // "material type" is a real MAGE-TAB sample column, so every file writing it before
            // "assay name" was reported out of order. Skipping unknowns takes the rule from 31
            // flagged files (28 of them wrong) to 3, all genuine.
            return null;
        }

        private static void ValidateRows(SdrfHeader header, IReadOnlyList<SdrfRow> rows,
            List<SdrfValidationMessage> messages)
        {
            if (rows.Count == 0)
            {
                messages.Add(new SdrfValidationMessage(
                    SdrfValidationSeverity.Warning, "NoRows", "The document has a header but no rows."));
                return;
            }

            var integerColumnIndexes = IntegerColumns
                .Select(c => (Column: c, Indexes: header.IndexesOf(c)))
                .Where(x => x.Indexes.Count > 0)
                .ToList();

            for (int r = 0; r < rows.Count; r++)
            {
                var row = rows[r];

                // Ragged rows are an Error: the cells no longer line up with the header, so every
                // value past the first missing column is being read under the wrong name.
                if (row.Cells.Count != header.Count)
                    messages.Add(new SdrfValidationMessage(
                        SdrfValidationSeverity.Error, "RowWidth",
                        $"Row has {row.Cells.Count} cells but the header has {header.Count}.", r));

                for (int c = 0; c < row.Cells.Count; c++)
                {
                    string cell = row.Cells[c];
                    string column = c < header.Count ? header[c] : $"(column {c})";

                    if (string.IsNullOrWhiteSpace(cell))
                    {
                        messages.Add(new SdrfValidationMessage(
                            SdrfValidationSeverity.Warning, "EmptyCell",
                            "Cell is empty; use a reserved word such as 'not available' or " +
                            "'not applicable' to say so explicitly.", r, column));
                        continue;
                    }

                    foreach (var reserved in ReservedWords)
                    {
                        if (!string.Equals(cell, reserved, StringComparison.OrdinalIgnoreCase))
                            continue;
                        if (!string.Equals(cell, reserved, StringComparison.Ordinal))
                            messages.Add(new SdrfValidationMessage(
                                SdrfValidationSeverity.Warning, "ReservedWordCase",
                                $"Reserved word '{cell}' should be lowercase '{reserved}'.", r, column));
                        break;
                    }
                }

                foreach (var (column, indexes) in integerColumnIndexes)
                {
                    foreach (int i in indexes)
                    {
                        if (i >= row.Cells.Count) continue;
                        string value = row.Cells[i];
                        if (string.IsNullOrWhiteSpace(value)) continue;
                        if (ReservedWords.Any(w => string.Equals(value, w, StringComparison.OrdinalIgnoreCase)))
                            continue;
                        if (int.TryParse(value, out int parsed) && parsed > 0) continue;

                        messages.Add(new SdrfValidationMessage(
                            SdrfValidationSeverity.Warning, "NonIntegerValue",
                            $"'{value}' is not a positive integer or a reserved word. Replicate and " +
                            "fraction identifiers are 1-based integers.", r, column));
                    }
                }
            }
        }

        /// <summary>
        /// The specification's one hard row constraint: source name + assay name + comment[label]
        /// must be unique. A duplicate means two rows are indistinguishable, so anything keyed on
        /// them silently collapses or double-counts.
        /// </summary>
        private static void ValidateRowKeyUniqueness(SdrfHeader header, IReadOnlyList<SdrfRow> rows,
            List<SdrfValidationMessage> messages)
        {
            if (!header.Contains("source name") || !header.Contains("assay name"))
                return; // already reported as a missing required column

            // A short row cannot be keyed. SdrfRow's indexer returns null both when a column is
            // absent and when the row is too short to reach it, so two rows that both fall short of
            // assay name collapsed onto the same "" key and were reported as duplicates of each
            // other -- a second, wrong diagnosis on top of the RowWidth error they already had.
            // Reach only as far as the columns a row MUST have to be keyable. comment[label] is
            // optional -- a missing one reads as "", which is a legitimate key component -- so
            // including it exempted perfectly keyable rows from the check whenever it happened to
            // be the rightmost of the three, turning a false positive into a false negative.
            int keyReach = new[] { "source name", "assay name" }
                .Select(header.IndexOf)
                .Where(i => i >= 0)
                .DefaultIfEmpty(-1)
                .Max();

            var seen = new Dictionary<string, int>(StringComparer.Ordinal);
            for (int r = 0; r < rows.Count; r++)
            {
                var row = rows[r];
                if (row.Cells.Count <= keyReach)
                    continue;
                // Tab is a safe key separator: cells come from splitting on tab, so no cell can
                // contain one and the composite key cannot collide by accident.
                string key = string.Join('\t', new[]
                {
                    row["source name"] ?? "",
                    row["assay name"] ?? "",
                    row["comment[label]"] ?? ""
                });

                if (seen.TryGetValue(key, out int first))
                    messages.Add(new SdrfValidationMessage(
                        SdrfValidationSeverity.Error, "RowKeyUniqueness",
                        "source name + assay name + comment[label] duplicates the combination first " +
                        $"seen on line {first + 2}; the two rows cannot be told apart.", r));
                else
                    seen[key] = r;
            }
        }
    }
}
