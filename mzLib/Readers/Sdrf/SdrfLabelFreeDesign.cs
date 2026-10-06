using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
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
    /// than the first failure, and writes nothing when anything fails.</para>
    ///
    /// <para><b>Numbers are the SDRF's.</b> A biological replicate number the SDRF gives is written exactly
    /// as given, never ranked or closed up: replicate 4 of one condition can be the same subject as replicate 4
    /// of another, and a gap can be a lost sample. Only where the SDRF gives no number (no column, or
    /// <c>not available</c>) is one added, per sample: the one number the sample's other rows give, or else one after
    /// the highest the condition uses, so an added number never fills a gap. That,
    /// and dropping rows that name a file the search does not read, are the only changes, and both are reported
    /// in <see cref="Notes"/>, as are gaps in the numbering and a sample name whose one number disagrees with
    /// its replicate.</para>
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
        /// What was added, dropped or worth checking on the way to a valid design: biological replicates added
        /// where the SDRF gives none (with the numbers), rows whose file is not searched, conditions whose
        /// replicates are not numbered 1..N, and sample names whose one number disagrees with their replicate.
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
        /// A human-readable account of the projection: the columns used, what was added or dropped, and
        /// every refusal. Print it wherever the design is built; it is the only record of a biological
        /// replicate number this reader added.
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

            string? keyColumn = SdrfDesignRules.ChooseFileKeyColumn(header, refusals);

            var conditionColumns = SdrfDesignRules.ResolveConditionColumns(header, options.ConditionColumns, refusals);
            bool conditionDeclared = options.ConditionColumns is { Count: > 0 };

            if (rows.Count == 0)
                refusals.Add("The SDRF has no rows.");

            if (refusals.Count > 0)
                return new SdrfLabelFreeDesign(new List<SpectraFileInfo>(), refusals, notes, keyColumn, conditionColumns);

            Dictionary<string, string>? searchedByName = null;
            if (options.SearchedFiles != null)
            {
                searchedByName = SdrfDesignRules.IndexSearchedFiles(options.SearchedFiles, refusals);
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
                string fileName = SdrfDesignRules.FileNameOf(rows[i], keyColumn!);
                if (searchedByName != null && !searchedByName.ContainsKey(fileName))
                {
                    dropped.Add((line, fileName));
                    notes.Add(SdrfDesignRules.DroppedRowNote(line, fileName, keyColumn!));
                    continue;
                }

                var row = ParseRow(rows[i], line, keyColumn!, conditionColumns, conditionDeclared, refusals);
                if (row != null)
                    parsed.Add(row);
                else if (fileName.Length > 0)
                    unreadFiles.Add(fileName);
            }

            RefuseRepeatedFiles(parsed, refusals);
            SdrfDesignRules.RefuseConditionCollisions(parsed.Select(r => (r.Condition, r.FactorValues)), refusals);

            if (searchedByName != null)
            {
                SdrfDesignRules.RefuseSearchedFilesWithoutRows(parsed.Select(r => r.FileName), unreadFiles, dropped, searchedByName, refusals);
                foreach (var row in parsed)
                    row.FilePath = searchedByName[row.FileName];
            }

            if (refusals.Count > 0)
                return new SdrfLabelFreeDesign(new List<SpectraFileInfo>(), refusals, notes, keyColumn, conditionColumns);

            var bioreps = NumberBiologicalReplicates(parsed, rowsWereDropped: dropped.Count > 0, notes, refusals);
            if (refusals.Count > 0)
                return new SdrfLabelFreeDesign(new List<SpectraFileInfo>(), refusals, notes, keyColumn, conditionColumns);

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
            /// <summary>The SDRF's number, or null where it gives none.</summary>
            public required int? BiologicalReplicate { get; init; }
            public required string SourceName { get; init; }
            public required int Fraction { get; init; }
            public required int TechnicalReplicate { get; init; }
        }

        private static ParsedRow? ParseRow(SdrfRow row, int line, string keyColumn, List<string> conditionColumns,
            bool conditionDeclared, List<string> refusals)
        {
            int before = refusals.Count;

            string? label = row[LabelColumn];
            if (label != null && label.IndexOf("label free sample", StringComparison.OrdinalIgnoreCase) < 0)
            {
                refusals.Add($"Line {line}: '{LabelColumn}' is '{label}', not label free. " +
                             $"An isobaric design is read with {nameof(SdrfIsobaricDesign)}.");
            }

            string fileName = SdrfDesignRules.FileNameOf(row, keyColumn);
            if (fileName.Length == 0)
                refusals.Add($"Line {line}: '{keyColumn}' is empty.");

            var factorValues = SdrfDesignRules.ReadFactorValues(row, line, fileName, conditionColumns, conditionDeclared, refusals);

            int? biorep = ReadBiologicalReplicate(row, line, fileName, refusals);
            int fraction = SdrfDesignRules.ReadPositiveInteger(row, FractionColumn, line, fileName, required: false, refusals);
            int techrep = SdrfDesignRules.ReadPositiveInteger(row, TechnicalReplicateColumn, line, fileName, required: false, refusals);

            if (refusals.Count > before)
                return null;

            return new ParsedRow
            {
                Line = line,
                FileName = fileName,
                FilePath = fileName,
                Condition = SdrfDesignRules.JoinCondition(factorValues),
                FactorValues = factorValues,
                BiologicalReplicate = biorep,
                SourceName = (row[SourceNameColumn] ?? string.Empty).Trim(),
                Fraction = fraction,
                TechnicalReplicate = techrep,
            };
        }

        // Label-free measures a file once. Case-insensitive, like SampleExperimentalDesign's keys.
        private static void RefuseRepeatedFiles(List<ParsedRow> rows, List<string> refusals)
        {
            foreach (var group in rows.GroupBy(r => r.FileName, StringComparer.OrdinalIgnoreCase).Where(g => g.Count() > 1))
            {
                refusals.Add($"'{group.Key}' is named by {group.Count()} rows (lines {string.Join(", ", group.Select(r => r.Line))}). " +
                             "A label-free design has one row per file.");
            }
        }

        /// <summary>
        /// The SDRF's biological replicate number, or null where it gives none: no column, an empty cell, or
        /// <c>not available</c> / <c>not applicable</c>. Anything else that is not an integer of 1 or more is
        /// refused, because it is a number the SDRF gives that cannot be kept.
        /// </summary>
        private static int? ReadBiologicalReplicate(SdrfRow row, int line, string fileName, List<string> refusals)
        {
            string text = (row[BiologicalReplicateColumn] ?? string.Empty).Trim();
            if (text.Length == 0 || SdrfDesignRules.IsUnknownWord(text))
                return null;

            if (int.TryParse(text, NumberStyles.None, CultureInfo.InvariantCulture, out int value) && value >= 1)
                return value;

            refusals.Add($"Line {line}{SdrfDesignRules.Describe(fileName)}: '{BiologicalReplicateColumn}' is '{text}', not an integer of 1 or more.");
            return null;
        }

        /// <summary>
        /// Every row's biological replicate, 1-based. A number the SDRF gives is kept exactly. A row the SDRF gives no
        /// number whose <c>source name</c> carries exactly one number on its other rows in the condition is that sample
        /// and takes its number; several numbers there is refused, since which one is meant cannot be told. Any other
        /// sample without a number (its rows share a <c>source name</c>; a row without one is a sample of its own) is
        /// numbered after the highest number the condition uses, in row order, so an added number never fills a gap.
        /// Reports what it added, any condition whose given numbers are not 1..N, and any sample name whose only number
        /// is not its replicate.
        /// </summary>
        private static Dictionary<ParsedRow, int> NumberBiologicalReplicates(List<ParsedRow> rows, bool rowsWereDropped,
            List<string> notes, List<string> refusals)
        {
            var numbers = new Dictionary<ParsedRow, int>();
            foreach (var condition in rows.GroupBy(r => r.Condition, StringComparer.Ordinal))
            {
                var given = condition.Where(r => r.BiologicalReplicate.HasValue).Select(r => r.BiologicalReplicate!.Value).ToHashSet();
                foreach (var row in condition.Where(r => r.BiologicalReplicate.HasValue))
                    numbers[row] = row.BiologicalReplicate!.Value;

                var givenBySource = condition.Where(r => r.BiologicalReplicate.HasValue && r.SourceName.Length > 0)
                    .GroupBy(r => r.SourceName, StringComparer.Ordinal)
                    .ToDictionary(g => g.Key, g => g.Select(r => r.BiologicalReplicate!.Value).Distinct().OrderBy(v => v).ToList(), StringComparer.Ordinal);

                int highest = given.Count > 0 ? given.Max() : 0;
                var added = new List<string>();
                foreach (var sample in condition.Where(r => !r.BiologicalReplicate.HasValue)
                             .GroupBy(r => r.SourceName.Length > 0 ? r.SourceName : "line " + r.Line, StringComparer.Ordinal))
                {
                    int number;
                    if (givenBySource.TryGetValue(sample.Key, out var ofSource) && ofSource.Count == 1)
                    {
                        number = ofSource[0];
                        added.Add($"'{sample.Key}' -> {number} (the number its other rows give)");
                    }
                    else if (ofSource != null)
                    {
                        foreach (var row in sample)
                        {
                            refusals.Add($"Line {row.Line}{SdrfDesignRules.Describe(row.FileName)}: the SDRF gives no biological replicate, " +
                                         $"and the other rows of '{sample.Key}' in condition '{condition.Key}' give several " +
                                         $"({string.Join(", ", ofSource)}), so which one it is cannot be told.");
                        }
                        number = ofSource[0];
                    }
                    else
                    {
                        number = ++highest;
                        added.Add($"'{sample.Key}' -> {number}");
                    }

                    foreach (var row in sample)
                        numbers[row] = number;
                }

                if (added.Count > 0)
                {
                    notes.Add($"Condition '{condition.Key}': the SDRF gives no biological replicate for {added.Count} sample(s), " +
                              "so these were numbered here: " + string.Join(", ", added) + ".");
                }

                // Over the numbers the SDRF gives, before any are added: an added number must not hide a gap.
                var present = given.OrderBy(v => v).ToList();
                if (present.Select((v, i) => v != i + 1).Any(gap => gap))
                {
                    notes.Add($"Condition '{condition.Key}': biological replicates {string.Join(", ", present)}, kept as the SDRF numbers them. " +
                              (rowsWereDropped
                                  ? "Rows were dropped from the SDRF, so a missing number may be a sample whose file is not searched."
                                  : "Every row of the SDRF is searched, so the gaps are how the SDRF numbers its samples, not lost samples."));
                }

                foreach (var sample in condition.Where(r => r.BiologicalReplicate.HasValue && r.SourceName.Length > 0)
                             .GroupBy(r => r.SourceName, StringComparer.Ordinal))
                {
                    var numbersInName = Regex.Matches(sample.Key, @"\d+");
                    if (numbersInName.Count != 1
                        || !int.TryParse(numbersInName[0].Value, NumberStyles.None, CultureInfo.InvariantCulture, out int inName))
                        continue;

                    var disagreeing = sample.Select(r => r.BiologicalReplicate!.Value).Distinct().Where(v => v != inName).ToList();
                    if (disagreeing.Count > 0)
                    {
                        notes.Add($"Line {sample.First().Line}: '{SourceNameColumn}' is '{sample.Key}', whose only number ({inName}) is not " +
                                  $"its biological replicate ({string.Join(", ", disagreeing)}). Kept as the SDRF gives it; check which is meant.");
                    }
                }
            }

            return numbers;
        }

        // The checks of MetaMorpheus's ExperimentalDesign.GetErrorsInExperimentalDesign, reporting
        // every failure instead of the first: within each condition and biological replicate the
        // fractions run 1..N with no gap (a missing LAST fraction is fine), within each fraction the
        // technical replicates run 1..N, and no (condition, biorep, fraction, techrep) repeats.
        // Its biological replicate rule (each condition numbered 1..N with no gap) is deliberately NOT
        // mirrored: numbers the SDRF gives are kept, so a gap is reported in Notes, not refused. Current
        // MetaMorpheus still rejects such a design; it ships with the MetaMorpheus change that turns that
        // refusal into a warning, together with #1422, which quantifies gapped designs.
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
    }
}
