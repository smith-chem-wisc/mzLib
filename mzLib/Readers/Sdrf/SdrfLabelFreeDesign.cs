using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using MassSpectrometry;

namespace Readers
{
    /// <summary>
    /// What a caller tells <see cref="SdrfLabelFreeDesign.Read(SdrfDocument, SdrfLabelFreeDesignOptions?)"/>
    /// that the SDRF alone cannot say.
    /// </summary>
    public sealed class SdrfLabelFreeDesignOptions
    {
        /// <summary>
        /// The <c>factor value[...]</c> columns the condition is built from, in order, by their exact
        /// column names (column names match case-sensitively). The condition is their values joined
        /// with <c>_</c>. Leave null or empty to use the document's only factor column; a document
        /// with several and none declared is refused (MAP-07, N5).
        /// </summary>
        public IReadOnlyList<string>? ConditionColumns { get; init; }

        /// <summary>
        /// The files the search will read, as paths or bare names. When given, the design is checked
        /// the way MetaMorpheus reads it: each searched file must be named EXACTLY (case and
        /// extension) by one row, and rows naming a file that is not searched are dropped and
        /// reported before any row is checked. Two searched paths sharing a file name, or a list naming
        /// no file, are refused. No stem matching happens here: joining searched files to rows is
        /// <c>SdrfSearchScope</c>'s job, and it writes the <c>comment[searched data file]</c> this
        /// reader keys on.
        /// </summary>
        public IReadOnlyCollection<string>? SearchedFiles { get; init; }
    }

    /// <summary>
    /// A label-free experimental design read from an SDRF (QuantProject M7): either a design that
    /// MetaMorpheus will accept, or a refusal that lists every reason it would not.
    ///
    /// <para><b>Why refuse rather than repair.</b> An invalid <c>ExperimentalDesign.tsv</c> is worse
    /// than none: MetaMorpheus skips quantification with one warning and no error exit. So this reader
    /// checks everything MetaMorpheus's own validator checks (MAP-32), reports all of it at once rather
    /// than the first failure, and writes nothing when anything fails. It relabels in only two ways,
    /// both reported in <see cref="Notes"/>: biological replicates are ranked within each condition
    /// (MAP-33), and rows naming a file the search does not read are dropped.</para>
    ///
    /// <para><b>Numbering.</b> The model (<see cref="SpectraFileInfo"/>) is 0-based; SDRF and
    /// <c>ExperimentalDesign.tsv</c> are 1-based. The conversion happens once on the way in and once
    /// in <see cref="WriteExperimentalDesignTsv"/>, and nowhere else.</para>
    ///
    /// <para>The rules are the design ↔ SDRF mapping table's read rules: MAP-07 (condition), MAP-12
    /// (file key), MAP-21/22 (technical replicate, fraction), MAP-32..34.</para>
    /// </summary>
    public sealed class SdrfLabelFreeDesign
    {
        public const string SearchedDataFileColumn = "comment[searched data file]";
        public const string DataFileColumn = "comment[data file]";
        public const string BiologicalReplicateColumn = "characteristics[biological replicate]";
        public const string TechnicalReplicateColumn = "comment[technical replicate]";
        public const string FractionColumn = "comment[fraction identifier]";
        public const string LabelColumn = "comment[label]";
        public const string SourceNameColumn = "source name";
        public const string FactorValuePrefix = "factor value[";

        /// <summary>The header MetaMorpheus writes and reads, <c>EngineLayer.ExperimentalDesign</c>.</summary>
        public const string ExperimentalDesignHeader = "FileName\tCondition\tBiorep\tFraction\tTechrep";

        /// <summary>
        /// Reserved words that state a factor is unknown. A condition built from one would pair
        /// samples nobody said were alike (PXD049018: +cGAMP and −cGAMP both written
        /// <c>not available</c>), so a factor column holding one is refused.
        /// </summary>
        private static readonly string[] UnknownFactorWords = { "not available", "not applicable" };

        private readonly List<SpectraFileInfo> _files;

        private SdrfLabelFreeDesign(List<SpectraFileInfo> files, List<string> refusals, List<string> notes,
            string? fileKeyColumn, List<string> conditionColumns)
        {
            _files = files;
            Refusals = refusals;
            Notes = notes;
            FileKeyColumn = fileKeyColumn;
            ConditionColumns = conditionColumns;
        }

        /// <summary>True when there is nothing in <see cref="Refusals"/> and the design can be used.</summary>
        public bool IsValid => Refusals.Count == 0;

        /// <summary>
        /// One entry per file, 0-based, in SDRF row order. Empty when the design is refused, so a
        /// caller that forgets to check <see cref="IsValid"/> gets no design rather than a wrong one.
        /// </summary>
        public IReadOnlyList<SpectraFileInfo> Files => IsValid ? _files : Array.Empty<SpectraFileInfo>();

        /// <summary>Every reason the design was refused. Empty when it is valid.</summary>
        public IReadOnlyList<string> Refusals { get; }

        /// <summary>
        /// What was relabelled or dropped on the way to a valid design: biological replicates
        /// renumbered within a condition (with the mapping), and rows whose file is not searched.
        /// </summary>
        public IReadOnlyList<string> Notes { get; }

        /// <summary>The column the files were keyed on (MAP-12), or null when neither was present.</summary>
        public string? FileKeyColumn { get; }

        /// <summary>The factor columns the condition was built from, in join order.</summary>
        public IReadOnlyList<string> ConditionColumns { get; }

        /// <summary>
        /// The design as mzLib quantification takes it. Throws when the design was refused.
        /// </summary>
        public SampleExperimentalDesign ToExperimentalDesign()
        {
            ThrowIfRefused();
            return SampleExperimentalDesign.LabelFree(_files);
        }

        /// <summary>
        /// Writes <c>ExperimentalDesign.tsv</c> in MetaMorpheus's format: 1-based, one row per file,
        /// the file named with its extension. Throws when the design was refused, so a refusal can
        /// never leave a file behind for MetaMorpheus to find.
        /// </summary>
        public void WriteExperimentalDesignTsv(string path)
        {
            ThrowIfRefused();

            var text = new StringBuilder();
            text.Append(ExperimentalDesignHeader).Append('\n');
            foreach (var file in _files)
            {
                text.Append(Path.GetFileName(file.FullFilePathWithExtension)).Append('\t')
                    .Append(file.Condition).Append('\t')
                    .Append((file.BiologicalReplicate + 1).ToString(CultureInfo.InvariantCulture)).Append('\t')
                    .Append((file.Fraction + 1).ToString(CultureInfo.InvariantCulture)).Append('\t')
                    .Append((file.TechnicalReplicate + 1).ToString(CultureInfo.InvariantCulture)).Append('\n');
            }

            File.WriteAllText(path, text.ToString(), new UTF8Encoding(false));
        }

        /// <summary>
        /// A human-readable account of the projection: the columns used, what was relabelled or
        /// dropped, and every refusal. Print it wherever the design is built; it is the only record
        /// of a renumbering (MAP-33).
        /// </summary>
        public string Report()
        {
            var report = new StringBuilder();
            report.AppendLine(IsValid
                ? $"Label-free design read from SDRF: {_files.Count} file(s), " +
                  $"{_files.Select(f => f.Condition).Distinct(StringComparer.Ordinal).Count()} condition(s)."
                : $"Label-free design REFUSED: {Refusals.Count} reason(s). No design was produced.");
            report.AppendLine($"Files keyed on: {FileKeyColumn ?? "(none)"}");
            report.AppendLine($"Condition from: {(ConditionColumns.Count == 0 ? "(none)" : string.Join(" + ", ConditionColumns))}");
            foreach (var note in Notes)
                report.AppendLine("  note: " + note);
            foreach (var refusal in Refusals)
                report.AppendLine("  refused: " + refusal);
            return report.ToString();
        }

        /// <summary>Reads an SDRF file. See <see cref="Read(SdrfDocument, SdrfLabelFreeDesignOptions?)"/>.</summary>
        public static SdrfLabelFreeDesign Read(string sdrfPath, SdrfLabelFreeDesignOptions? options = null)
        {
            if (sdrfPath == null)
                throw new ArgumentNullException(nameof(sdrfPath));
            return Read(new SdrfDocument(sdrfPath), options);
        }

        /// <summary>
        /// Projects an SDRF onto a label-free design. Never throws for a bad document: every problem
        /// is a refusal, listed in <see cref="Refusals"/>.
        /// </summary>
        public static SdrfLabelFreeDesign Read(SdrfDocument sdrf, SdrfLabelFreeDesignOptions? options = null)
        {
            if (sdrf == null)
                throw new ArgumentNullException(nameof(sdrf));
            options ??= new SdrfLabelFreeDesignOptions();

            var refusals = new List<string>();
            var notes = new List<string>();
            var header = sdrf.Header;
            var rows = sdrf.Results.ToList();

            // MAP-12: key on the file the search read when the document says which, else the acquired one.
            string? keyColumn = header.Contains(SearchedDataFileColumn) ? SearchedDataFileColumn
                : header.Contains(DataFileColumn) ? DataFileColumn
                : null;
            if (keyColumn == null)
                refusals.Add($"The SDRF has neither '{SearchedDataFileColumn}' nor '{DataFileColumn}', so no row names a file.");

            if (!header.Contains(BiologicalReplicateColumn))
                refusals.Add($"The SDRF has no '{BiologicalReplicateColumn}' column. MetaMorpheus needs a biological replicate for every file.");

            var conditionColumns = ResolveConditionColumns(header, options.ConditionColumns, refusals);
            bool conditionDeclared = options.ConditionColumns is { Count: > 0 };

            if (rows.Count == 0)
                refusals.Add("The SDRF has no rows.");

            if (refusals.Count > 0)
                return new SdrfLabelFreeDesign(new List<SpectraFileInfo>(), refusals, notes, keyColumn, conditionColumns);

            Dictionary<string, string>? searchedByName = null;
            if (options.SearchedFiles != null)
            {
                searchedByName = IndexSearchedFiles(options.SearchedFiles, refusals);
                if (refusals.Count > 0)
                    return new SdrfLabelFreeDesign(new List<SpectraFileInfo>(), refusals, notes, keyColumn, conditionColumns);
            }

            // Read every row, collecting every problem rather than stopping at the first. A row for a
            // file the search does not read is dropped BEFORE it is checked, as MetaMorpheus's reader
            // skips it: a problem in a row nothing reads must not refuse the design.
            var parsed = new List<ParsedRow>();
            var unreadFiles = new HashSet<string>(StringComparer.Ordinal);
            var dropped = new List<(int Line, string FileName)>();
            for (int i = 0; i < rows.Count; i++)
            {
                int line = i + 2;
                string? fileCell = rows[i][keyColumn!];
                string fileName = string.IsNullOrWhiteSpace(fileCell) ? string.Empty : Path.GetFileName(fileCell.Trim());
                if (searchedByName != null && !searchedByName.ContainsKey(fileName))
                {
                    dropped.Add((line, fileName));
                    notes.Add(fileName.Length == 0
                        ? $"Line {line} dropped: '{keyColumn}' is empty, so it names no searched file."
                        : $"Line {line} ('{fileName}') dropped: the search does not read that file.");
                    continue;
                }

                var row = ParseRow(rows[i], line, keyColumn!, conditionColumns, conditionDeclared, refusals);
                if (row != null)
                    parsed.Add(row);
                else if (fileName.Length > 0)
                    unreadFiles.Add(fileName);
            }

            RefuseRepeatedFiles(parsed, refusals);
            RefuseConditionCollisions(parsed, refusals);

            if (searchedByName != null)
                UseSearchedPaths(parsed, unreadFiles, dropped, searchedByName, refusals);

            if (refusals.Count > 0)
                return new SdrfLabelFreeDesign(new List<SpectraFileInfo>(), refusals, notes, keyColumn, conditionColumns);

            var bioreps = RankBiologicalReplicates(parsed, notes);

            var files = parsed
                .Select(r => new SpectraFileInfo(r.FilePath, r.Condition, bioreps[r] - 1, r.TechnicalReplicate - 1, r.Fraction - 1))
                .ToList();

            RefuseWhatMetaMorpheusWouldReject(files, refusals);

            return new SdrfLabelFreeDesign(files, refusals, notes, keyColumn, conditionColumns);
        }

        private sealed class ParsedRow
        {
            public required int Line { get; init; }
            public required string FileName { get; init; }
            public required string FilePath { get; set; }
            public required string Condition { get; init; }
            public required IReadOnlyList<string> FactorValues { get; init; }
            public required int BiologicalReplicate { get; init; }
            public required int Fraction { get; init; }
            public required int TechnicalReplicate { get; init; }
        }

        // MAP-07 as amended by QP-S14: the declared columns, in order; else the only factor column.
        private static List<string> ResolveConditionColumns(SdrfHeader header, IReadOnlyList<string>? declared,
            List<string> refusals)
        {
            var factorColumns = header
                .Where(name => name.StartsWith(FactorValuePrefix, StringComparison.Ordinal))
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

        private static ParsedRow? ParseRow(SdrfRow row, int line, string keyColumn, List<string> conditionColumns,
            bool conditionDeclared, List<string> refusals)
        {
            int before = refusals.Count;

            string? label = row[LabelColumn];
            if (label != null && label.IndexOf("label free sample", StringComparison.OrdinalIgnoreCase) < 0)
            {
                refusals.Add($"Line {line}: '{LabelColumn}' is '{label}', not label free. " +
                             "Only label-free designs are read from SDRF so far; an isobaric design has no projection yet.");
            }

            string? fileCell = row[keyColumn];
            string fileName = string.IsNullOrWhiteSpace(fileCell) ? string.Empty : Path.GetFileName(fileCell.Trim());
            if (fileName.Length == 0)
                refusals.Add($"Line {line}: '{keyColumn}' is empty.");

            var factorValues = new List<string>();
            foreach (var column in conditionColumns)
            {
                string value = (row[column] ?? string.Empty).Trim();
                if (value.Length == 0)
                    refusals.Add($"Line {line}{Describe(fileName)}: '{column}' is empty.");
                else if (UnknownFactorWords.Any(w => string.Equals(value, w, StringComparison.OrdinalIgnoreCase)))
                    // SDRF-P1 (sdrf 012): the refusal names its remedy, and the remedy depends on how the column was chosen.
                    refusals.Add($"Line {line}{Describe(fileName)}: '{column}' is '{value}'. A condition cannot be built from an unknown " +
                                 (conditionDeclared
                                     ? "factor; fill it in, or leave the column out of the declared condition columns to pool these rows."
                                     : "factor. It was used because it is the only factor value column and none was declared; " +
                                       "fill it in, or declare the column(s) the condition should be built from instead."));
                factorValues.Add(value);
            }

            int biorep = ReadPositiveInteger(row, BiologicalReplicateColumn, line, fileName, required: true, refusals);
            int fraction = ReadPositiveInteger(row, FractionColumn, line, fileName, required: false, refusals);
            int techrep = ReadPositiveInteger(row, TechnicalReplicateColumn, line, fileName, required: false, refusals);

            if (refusals.Count > before)
                return null;

            return new ParsedRow
            {
                Line = line,
                FileName = fileName,
                FilePath = fileName,
                Condition = string.Join("_", factorValues),
                FactorValues = factorValues,
                BiologicalReplicate = biorep,
                Fraction = fraction,
                TechnicalReplicate = techrep,
            };
        }

        // MAP-34 (fractions) and MAP-21 (technical replicates): copied as they are, never renumbered;
        // anything that is not an integer >= 1 is refused. An absent optional column means 1.
        private static int ReadPositiveInteger(SdrfRow row, string column, int line, string fileName, bool required,
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

        private static string Describe(string fileName) => fileName.Length == 0 ? string.Empty : $" ({fileName})";

        // Label-free measures a file once. Case-insensitive, like SampleExperimentalDesign's keys.
        private static void RefuseRepeatedFiles(List<ParsedRow> rows, List<string> refusals)
        {
            foreach (var group in rows.GroupBy(r => r.FileName, StringComparer.OrdinalIgnoreCase).Where(g => g.Count() > 1))
            {
                refusals.Add($"'{group.Key}' is named by {group.Count()} rows (lines {string.Join(", ", group.Select(r => r.Line))}). " +
                             "A label-free design has one row per file.");
            }
        }

        // QP-S14: two different combinations of factor values that join to the same text are refused,
        // never merged. So is one condition spelled two ways, because cell values match case-insensitively.
        private static void RefuseConditionCollisions(List<ParsedRow> rows, List<string> refusals)
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
        private static Dictionary<string, string> IndexSearchedFiles(IReadOnlyCollection<string> searchedFiles,
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

        private static void UseSearchedPaths(List<ParsedRow> rows, HashSet<string> unreadFiles,
            List<(int Line, string FileName)> dropped, Dictionary<string, string> searchedByName, List<string> refusals)
        {
            var named = new HashSet<string>(rows.Select(r => r.FileName), StringComparer.Ordinal);

            // A row already refused for its contents names its file; saying it has no row would be wrong.
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

            foreach (var row in rows)
                row.FilePath = searchedByName[row.FileName];
        }

        // MAP-33: rank the biological replicates within each condition, keeping their order.
        // Study-wide numbering (control 1-3, treated 4-6) is common and MetaMorpheus rejects it.
        private static Dictionary<ParsedRow, int> RankBiologicalReplicates(List<ParsedRow> rows, List<string> notes)
        {
            var ranks = new Dictionary<ParsedRow, int>();
            foreach (var condition in rows.GroupBy(r => r.Condition, StringComparer.Ordinal))
            {
                var values = condition.Select(r => r.BiologicalReplicate).Distinct().OrderBy(v => v).ToList();
                for (int i = 0; i < values.Count; i++)
                {
                    foreach (var row in condition.Where(r => r.BiologicalReplicate == values[i]))
                        ranks[row] = i + 1;
                }

                if (values.Select((v, i) => v != i + 1).Any(changed => changed))
                {
                    notes.Add($"Condition '{condition.Key}': biological replicates renumbered " +
                              string.Join(", ", values.Select((v, i) => $"{v} -> {i + 1}")) + ".");
                }
            }

            return ranks;
        }

        // The checks of MetaMorpheus's ExperimentalDesign.GetErrorsInExperimentalDesign, reporting
        // every failure instead of the first: within each condition and biological replicate the
        // fractions run 1..N with no gap (a missing LAST fraction is fine), within each fraction the
        // technical replicates run 1..N, and no (condition, biorep, fraction, techrep) repeats.
        private static void RefuseWhatMetaMorpheusWouldReject(List<SpectraFileInfo> files, List<string> refusals)
        {
            foreach (var condition in files.GroupBy(f => f.Condition, StringComparer.Ordinal))
            {
                foreach (var biorep in condition.GroupBy(f => f.BiologicalReplicate).OrderBy(g => g.Key))
                {
                    int fractions = biorep.Max(f => f.Fraction) + 1;
                    for (int fraction = 0; fraction < fractions; fraction++)
                    {
                        var inFraction = biorep.Where(f => f.Fraction == fraction).ToList();
                        string where = $"Condition '{condition.Key}' biorep {biorep.Key + 1} fraction {fraction + 1}";
                        if (inFraction.Count == 0)
                        {
                            refusals.Add($"{where} is missing. Fractions are copied as they are, never renumbered, " +
                                         "because they are matched across samples by number.");
                            continue;
                        }

                        int techreps = inFraction.Max(f => f.TechnicalReplicate) + 1;
                        for (int techrep = 0; techrep < techreps; techrep++)
                        {
                            var inTechrep = inFraction.Where(f => f.TechnicalReplicate == techrep).ToList();
                            if (inTechrep.Count == 0)
                                refusals.Add($"{where} techrep {techrep + 1} is missing.");
                            else if (inTechrep.Count > 1)
                                refusals.Add($"{where} techrep {techrep + 1} is named by {inTechrep.Count} files: " +
                                             string.Join(", ", inTechrep.Select(f => $"'{Path.GetFileName(f.FullFilePathWithExtension)}'")) +
                                             ". Different samples cannot share one replicate; declare the factor that tells them apart.");
                        }
                    }
                }
            }
        }

        private void ThrowIfRefused()
        {
            if (!IsValid)
                throw new InvalidOperationException("The design was refused and cannot be used:" + Environment.NewLine + Report());
        }

        private sealed class SequenceComparer : IEqualityComparer<IReadOnlyList<string>>
        {
            public bool Equals(IReadOnlyList<string>? x, IReadOnlyList<string>? y) =>
                x != null && y != null && x.SequenceEqual(y, StringComparer.Ordinal);

            public int GetHashCode(IReadOnlyList<string> values) =>
                values.Aggregate(17, (hash, v) => unchecked(hash * 31 + StringComparer.Ordinal.GetHashCode(v)));
        }
    }
}
