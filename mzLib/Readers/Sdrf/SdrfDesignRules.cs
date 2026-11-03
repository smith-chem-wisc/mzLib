using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Readers
{
    /// <summary>
    /// The read rules every SDRF → design projection shares (QuantProject M7): which column names the
    /// file (MAP-12), how the condition is built (MAP-07), how replicate and fraction numbers are read
    /// (MAP-21, MAP-22, MAP-34), and how the searched files are matched. Kept in one place so the
    /// label-free and isobaric projections cannot drift apart on any of them.
    /// </summary>
    internal static class SdrfDesignRules
    {
        /// <summary>
        /// Reserved words that state a value is unknown. A condition built from one would pair
        /// samples nobody said were alike (PXD049018: +cGAMP and −cGAMP both written
        /// <c>not available</c>).
        /// </summary>
        internal static readonly string[] UnknownWords = { "not available", "not applicable" };

        internal static bool IsUnknownWord(string value) =>
            UnknownWords.Any(w => string.Equals(value.Trim(), w, StringComparison.OrdinalIgnoreCase));

        /// <summary>
        /// MAP-12: the file the search read when the document says which, else the acquired one.
        /// Null (and a refusal) when neither column is present.
        /// </summary>
        internal static string? ChooseFileKeyColumn(SdrfHeader header, List<string> refusals)
        {
            string? keyColumn = header.Contains(SdrfLabelFreeDesign.SearchedDataFileColumn) ? SdrfLabelFreeDesign.SearchedDataFileColumn
                : header.Contains(SdrfLabelFreeDesign.DataFileColumn) ? SdrfLabelFreeDesign.DataFileColumn
                : null;
            if (keyColumn == null)
                refusals.Add($"The SDRF has neither '{SdrfLabelFreeDesign.SearchedDataFileColumn}' nor '{SdrfLabelFreeDesign.DataFileColumn}', so no row names a file.");
            return keyColumn;
        }

        /// <summary>The bare file name a row's key cell names, or empty when the cell is blank.</summary>
        internal static string FileNameOf(SdrfRow row, string keyColumn)
        {
            string? cell = row[keyColumn];
            return string.IsNullOrWhiteSpace(cell) ? string.Empty : Path.GetFileName(cell.Trim());
        }

        // MAP-07 as amended by QP-S14: the declared columns, in order; else the only factor column.
        internal static List<string> ResolveConditionColumns(SdrfHeader header, IReadOnlyList<string>? declared,
            List<string> refusals)
        {
            var factorColumns = header
                .Where(name => name.StartsWith(SdrfLabelFreeDesign.FactorValuePrefix, StringComparison.Ordinal))
                .Distinct(StringComparer.Ordinal)
                .ToList();
            string available = factorColumns.Count == 0
                ? "the SDRF has no factor value columns"
                : "its factor value columns are " + string.Join(", ", factorColumns.Select(c => $"'{c}'"));

            if (declared != null && declared.Count > 0)
            {
                var missing = declared.Where(c => !header.Contains(c)).ToList();
                foreach (var column in missing)
                    refusals.Add($"The declared condition column '{column}' is not in the SDRF; {available}.");

                var repeated = declared.GroupBy(c => c, StringComparer.Ordinal).Where(g => g.Count() > 1).Select(g => g.Key);
                foreach (var column in repeated)
                    refusals.Add($"The condition column '{column}' is declared more than once.");

                return declared.ToList();
            }

            if (factorColumns.Count == 1)
                return factorColumns;

            refusals.Add(factorColumns.Count == 0
                ? "No condition: the SDRF has no factor value columns and none was declared. Declare the column(s) the condition is built from."
                : $"Several factor value columns and none declared ({string.Join(", ", factorColumns.Select(c => $"'{c}'"))}). " +
                  "Declare which of them make up the condition; using all of them could split conditions on a nuisance factor.");
            return new List<string>();
        }

        /// <summary>
        /// Reads one row's factor values for the condition (MAP-07). An empty or unknown value is a
        /// refusal (SDRF-P1), and the refusal names its remedy, which depends on how the column was chosen.
        /// </summary>
        internal static List<string> ReadFactorValues(SdrfRow row, int line, string fileName, List<string> conditionColumns,
            bool conditionDeclared, List<string> refusals)
        {
            var factorValues = new List<string>();
            foreach (var column in conditionColumns)
            {
                string value = (row[column] ?? string.Empty).Trim();
                if (value.Length == 0)
                    refusals.Add($"Line {line}{Describe(fileName)}: '{column}' is empty.");
                else if (IsUnknownWord(value))
                    refusals.Add($"Line {line}{Describe(fileName)}: '{column}' is '{value}'. A condition cannot be built from an unknown " +
                                 (conditionDeclared
                                     ? "factor; fill it in, or leave the column out of the declared condition columns to pool these rows."
                                     : "factor. It was used because it is the only factor value column and none was declared; " +
                                       "fill it in, or declare the column(s) the condition should be built from instead."));
                factorValues.Add(value);
            }

            return factorValues;
        }

        /// <summary>The condition a row's factor values make: joined with <c>_</c> in declared order.</summary>
        internal static string JoinCondition(IReadOnlyList<string> factorValues) => string.Join("_", factorValues);

        // MAP-34 (fractions) and MAP-21 (technical replicates): copied as they are, never renumbered;
        // anything that is not an integer >= 1 is refused. An absent optional column means 1.
        internal static int ReadPositiveInteger(SdrfRow row, string column, int line, string fileName, bool required,
            List<string> refusals)
        {
            string? cell = row[column];
            if (cell == null && !row.Header.Contains(column) && !required)
                return 1;

            string text = (cell ?? string.Empty).Trim();
            if (int.TryParse(text, NumberStyles.None, CultureInfo.InvariantCulture, out int value) && value >= 1)
                return value;

            refusals.Add($"Line {line}{Describe(fileName)}: '{column}' is '{text}', not an integer of 1 or more.");
            return 0;
        }

        internal static string Describe(string fileName) => fileName.Length == 0 ? string.Empty : $" ({fileName})";

        // QP-S14: two different combinations of factor values that join to the same text are refused,
        // never merged. So is one condition spelled two ways, because cell values match case-insensitively.
        internal static void RefuseConditionCollisions(IEnumerable<(string Condition, IReadOnlyList<string> FactorValues)> rows,
            List<string> refusals)
        {
            foreach (var group in rows.GroupBy(r => r.Condition, StringComparer.OrdinalIgnoreCase))
            {
                var spellings = group
                    .Select(r => r.FactorValues)
                    .Distinct(new SequenceComparer())
                    .ToList();
                if (spellings.Count < 2)
                    continue;

                refusals.Add($"Condition '{group.Key}' comes from {spellings.Count} different sets of factor values: " +
                             string.Join("; ", spellings.Select(s => "(" + string.Join(", ", s.Select(v => $"'{v}'")) + ")")) +
                             ". They are refused rather than merged.");
            }
        }

        // Exact names, as MetaMorpheus's reader compares them (ordinal, extension included). Two searched
        // paths with one name are refused: the design names a file by its name alone, and MetaMorpheus
        // matches each row to the first path of that name, so the other is never defined.
        internal static Dictionary<string, string> IndexSearchedFiles(IReadOnlyCollection<string> searchedFiles,
            List<string> refusals)
        {
            var searchedByName = new Dictionary<string, string>(StringComparer.Ordinal);
            var usable = searchedFiles.Where(p => !string.IsNullOrWhiteSpace(p)).ToList();
            if (usable.Count == 0)
            {
                refusals.Add("SearchedFiles was given but names no file, so no row could be kept.");
                return searchedByName;
            }

            foreach (var group in usable.GroupBy(p => Path.GetFileName(p.Trim()), StringComparer.Ordinal))
            {
                var paths = group.Distinct(StringComparer.Ordinal).ToList();
                if (paths.Count > 1)
                {
                    refusals.Add($"Searched files {string.Join(", ", paths.Select(p => $"'{p}'"))} share the name '{group.Key}'. " +
                                 "A design names a file by its name alone, so MetaMorpheus cannot tell them apart; rename them or search them separately.");
                }
                searchedByName[group.Key] = paths[0];
            }

            return searchedByName;
        }

        /// <summary>
        /// Refuses every searched file no kept row names. A file whose rows were refused for their
        /// contents is not reported again: saying it has no row would be wrong.
        /// </summary>
        internal static void RefuseSearchedFilesWithoutRows(IEnumerable<string> namedFiles, HashSet<string> unreadFiles,
            List<(int Line, string FileName)> dropped, Dictionary<string, string> searchedByName, List<string> refusals)
        {
            var named = new HashSet<string>(namedFiles, StringComparer.Ordinal);

            foreach (var name in searchedByName.Keys.Where(n => !named.Contains(n) && !unreadFiles.Contains(n)))
            {
                string stem = Path.GetFileNameWithoutExtension(name);
                var nearly = dropped.FirstOrDefault(d => string.Equals(
                    Path.GetFileNameWithoutExtension(d.FileName), stem, StringComparison.OrdinalIgnoreCase));
                refusals.Add(nearly.FileName == null
                    ? $"Searched file '{name}' has no SDRF row. MetaMorpheus skips quantification when any searched file is missing from the design."
                    : $"Searched file '{name}' has no SDRF row naming it exactly; line {nearly.Line} names '{nearly.FileName}'. " +
                      "Restrict the SDRF to the searched files first (SdrfSearchScope), which records the searched name.");
            }
        }

        /// <summary>The note recorded for a row dropped because the search does not read its file.</summary>
        internal static string DroppedRowNote(int line, string fileName, string keyColumn) =>
            fileName.Length == 0
                ? $"Line {line} dropped: '{keyColumn}' is empty, so it names no searched file."
                : $"Line {line} ('{fileName}') dropped: the search does not read that file.";

        private sealed class SequenceComparer : IEqualityComparer<IReadOnlyList<string>>
        {
            public bool Equals(IReadOnlyList<string>? x, IReadOnlyList<string>? y) =>
                x != null && y != null && x.SequenceEqual(y, StringComparer.Ordinal);

            public int GetHashCode(IReadOnlyList<string> values) =>
                values.Aggregate(17, (hash, v) => unchecked(hash * 31 + StringComparer.Ordinal.GetHashCode(v)));
        }
    }
}
