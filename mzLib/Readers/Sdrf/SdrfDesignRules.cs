using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

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

        // GR-20 (sdrf QP-S28 and QP-S31, ruled 2026-10-08): the sample-type words that mark a pooled reference, and the
        // pooled-sample word that marks none.
        private static readonly HashSet<string> ReferenceSampleTypes = new(StringComparer.OrdinalIgnoreCase) { "reference", SdrfReserved.Pooled };
        private const string NotPooled = "not pooled";

        /// <summary>True for the biological replicate the SDRF specification prescribes for a pooled sample.</summary>
        internal static bool IsPooledWord(string value) =>
            string.Equals(value.Trim(), SdrfReserved.Pooled, StringComparison.OrdinalIgnoreCase);

        /// <summary>
        /// The pooled-reference signals a row carries (GR-20), each read on its own; any one makes the row's sample a
        /// reference. They come in two kinds. <c>Role</c> states the sample's role: <c>characteristics[sample type]</c>
        /// <c>reference</c> or <c>pooled</c>, and, for an isobaric channel, <c>comment[reference channel]</c> naming it
        /// (<paramref name="namesThisChannel"/>; a label-free file has no channel). <c>Preparation</c> can also say only how
        /// the sample was made: <c>characteristics[pooled sample]</c> other than <c>not pooled</c>, and a biological replicate
        /// of <c>pooled</c>. A kind every sample carries tells no sample apart and is cleared (<see cref="EverySampleIsMarked"/>).
        /// Beside them, what the row says that reads as a study sample (a sample type naming something else,
        /// <c>not pooled</c>), so a disagreement can be reported. Empty cells and unknown words say neither. A role is never
        /// read from a name.
        /// </summary>
        internal static (List<string> Role, List<string> Preparation, List<string> StudySample) ReadReferenceSignals(SdrfRow row,
            Func<string, bool>? namesThisChannel)
        {
            var role = new List<string>();
            var preparation = new List<string>();
            var studySample = new List<string>();

            string sampleTypeCell = (row[SdrfLabelFreeDesign.SampleTypeColumn] ?? string.Empty).Trim();
            string sampleType = sampleTypeCell.Length == 0 ? string.Empty : SdrfQuantAuditor.ReadLabel(sampleTypeCell).Name;
            if (ReferenceSampleTypes.Contains(sampleType))
                role.Add($"sample type '{sampleType}'");
            else if (sampleType.Length > 0 && !IsUnknownWord(sampleType))
                studySample.Add($"sample type '{sampleType}'");

            // The deposits write the plex's reference channel on every row (TMT127N), so it marks only the row it names.
            string referenceChannel = (row[SdrfIsobaricDesign.ReferenceChannelColumn] ?? string.Empty).Trim();
            if (namesThisChannel != null && referenceChannel.Length > 0 && !IsUnknownWord(referenceChannel)
                && namesThisChannel(referenceChannel))
            {
                role.Add($"reference channel '{referenceChannel}'");
            }

            string pooledSample = (row[SdrfLabelFreeDesign.PooledSampleColumn] ?? string.Empty).Trim();
            if (string.Equals(pooledSample, NotPooled, StringComparison.OrdinalIgnoreCase))
                studySample.Add($"pooled sample '{pooledSample}'");
            else if (pooledSample.Length > 0 && !IsUnknownWord(pooledSample))
                preparation.Add("pooled sample");

            string biorep = (row[SdrfLabelFreeDesign.BiologicalReplicateColumn] ?? string.Empty).Trim();
            if (IsPooledWord(biorep))
                preparation.Add($"biological replicate '{biorep}'");

            return (role, preparation, studySample);
        }

        /// <summary>
        /// sdrf 041: two signals that disagree are reported, and the sample is read as a reference, because a pool read as a
        /// study sample pulls its group towards the mean of all groups.
        /// </summary>
        internal static string DisagreementNote(string what, IEnumerable<string> reference, IEnumerable<string> studySample) =>
            $"{what} is marked a reference by {string.Join(", ", reference.Distinct())}, but also says " +
            $"{string.Join(", ", studySample.Distinct())}. It is read as a reference: a pool read as a study sample would pull " +
            "its group towards the mean of all groups.";

        /// <summary>
        /// Which kinds of a row's reference signals a whole design clears (<see cref="EverySampleIsMarked"/>).
        /// </summary>
        internal readonly record struct ClearedReferenceSignals(bool Role, bool Preparation)
        {
            /// <summary>The row's reference signals the design still reads, roles first; any one makes its sample a reference.</summary>
            internal List<string> Kept(IEnumerable<string> role, IEnumerable<string> preparation) =>
                (Role ? Enumerable.Empty<string>() : role).Concat(Preparation ? Enumerable.Empty<string>() : preparation).ToList();
        }

        /// <summary>
        /// GR-20: a reference signal every sample of a whole design carries tells no sample apart, so it is cleared, and a note
        /// says so. A role (sample type <c>reference</c> or <c>pooled</c>, a named reference channel) is cleared only when every
        /// sample has one (ruled 2026-11-03). A preparation mark (<c>pooled sample</c>, a <c>pooled</c> replicate) says how the
        /// samples were made, each itself a pool, not that they are references (ruled 2026-10-26). It is cleared when every
        /// sample without a role that stands carries one, so a role on only some samples still marks them references.
        /// Judged over the whole design, never per plex: a plex holding only reference channels beside plexes of study samples
        /// keeps its references. A signal on only some samples marks them references.
        /// </summary>
        /// <param name="signalsByRow">The role and preparation signals of every row the design reads (dropped rows excluded).</param>
        /// <param name="notes">The design's notes; one sentence is added for each kind cleared, naming the signal(s).</param>
        internal static ClearedReferenceSignals EverySampleIsMarked(
            IReadOnlyList<(IReadOnlyList<string> Role, IReadOnlyList<string> Preparation)> signalsByRow, List<string> notes)
        {
            if (signalsByRow.Count == 0)
                return default;

            bool role = signalsByRow.All(s => s.Role.Count > 0);
            var judged = role ? signalsByRow : signalsByRow.Where(s => s.Role.Count == 0).ToList();
            bool preparation = judged.Count > 0 && judged.All(s => s.Preparation.Count > 0);

            if (role)
            {
                var roles = signalsByRow.SelectMany(s => s.Role).Distinct(StringComparer.OrdinalIgnoreCase).ToList();
                bool preparationMarksSome = !preparation && judged.Any(s => s.Preparation.Count > 0);
                notes.Add($"Every sample in this design is marked a reference by {string.Join(" or ", roles)}, so the role tells no " +
                          "sample apart: all are read as study samples" +
                          (preparationMarksSome ? ", except those that pooled sample or a pooled replicate marks as a reference." : "."));
            }

            if (preparation)
            {
                var signals = judged.SelectMany(s => s.Preparation).Distinct(StringComparer.OrdinalIgnoreCase).ToList();
                notes.Add($"Every sample in this design is marked by {string.Join(" or ", signals)}, so the mark says how the samples " +
                          "were made, not that they are pooled references: all are read as study samples" +
                          (judged.Count == signalsByRow.Count ? "." : ", except those that sample type or reference channel marks as a reference."));
            }

            return new ClearedReferenceSignals(role, preparation);
        }

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

        /// <summary>
        /// The bare file name a row's key cell names, or empty when the cell is blank. A path is split at either separator,
        /// whatever the platform (<see cref="SdrfFileNamePattern.FileName"/>): an SDRF written on Windows and read on Linux
        /// names <c>D:\search\X-calib.mzML</c>.
        /// </summary>
        internal static string FileNameOf(SdrfRow row, string keyColumn)
        {
            string? cell = row[keyColumn];
            return string.IsNullOrWhiteSpace(cell) ? string.Empty : SdrfFileNamePattern.FileName(cell);
        }

        /// <summary>How a design's condition was chosen, which decides what a refusal offers as the way on.</summary>
        internal enum ConditionSource
        {
            /// <summary>The caller declared the columns.</summary>
            Declared,
            /// <summary>None was declared and the SDRF has one factor value column.</summary>
            OnlyColumn,
            /// <summary>None was declared and the SDRF has several factor value columns: all of them (QP-S42 c).</summary>
            EveryColumn,
            /// <summary>The caller stated that every sample is one condition (QP-S42 b).</summary>
            Stated,
            /// <summary>No condition could be chosen; the design is refused.</summary>
            None,
        }

        /// <summary>
        /// How a design's condition is built (MAP-07, as amended by QP-S14 and by the rulings of 2026-11-11, QP-S42 to
        /// QP-S44): the columns whose values are joined, or the one condition the caller stated, and whether a factor
        /// value nobody stated is kept as a level of its own.
        /// </summary>
        internal sealed class ConditionRule
        {
            internal ConditionRule(List<string> columns, ConditionSource source, string? singleCondition, bool unknownIsOwnCondition)
            {
                Columns = columns;
                Source = source;
                SingleCondition = singleCondition;
                UnknownIsOwnCondition = unknownIsOwnCondition;
            }

            /// <summary>The columns read, in join order. Empty when the condition was stated or none could be chosen.</summary>
            internal List<string> Columns { get; }

            internal ConditionSource Source { get; }

            /// <summary>The one condition every sample belongs to, when the caller stated it; else null.</summary>
            internal string? SingleCondition { get; }

            /// <summary>The caller stated that a factor value nobody stated is a level of its own (QP-S44).</summary>
            internal bool UnknownIsOwnCondition { get; }

            /// <summary>The report's "Condition from" line.</summary>
            internal string Description => SingleCondition != null
                ? $"one condition, '{SingleCondition}', as stated"
                : Columns.Count == 0 ? "(none)" : string.Join(" + ", Columns);
        }

        /// <summary>The <c>not available</c> a blank factor value is kept as, when the caller states such files are their own condition.</summary>
        internal const string BlankFactorValue = "not available";

        // sdrf's census (2026-11-07) leaves these out of "condition-like": they differ between the replicates of one group
        // by nature. Used only to point a refusal at what might be the groups, never to build one.
        private static readonly Regex ReplicateLevelCharacteristic = new(
            @"replicate|pooled sample|sample type|individual|\[sex\]|\[age\]|ancestry|cell identifier|sample id|biosample|mass\]|" +
            @"weight|bmi|height|ethnicity|race|sample (name|number)|patient|donor|subject|animal",
            RegexOptions.IgnoreCase | RegexOptions.CultureInvariant);

        /// <summary>
        /// Chooses how the condition is built (MAP-07): the declared columns, in order; else the only factor value column;
        /// else every factor value column, in the SDRF's order, with a note (QP-S42 c). With no factor value column and
        /// none declared, the design is refused, and the refusal names any condition-like characteristic that differs
        /// between <paramref name="rows"/> (QP-S42 b). A stated single condition makes every sample that condition and
        /// sets every factor value column aside, with a note naming them (ruled 2026-11-11); it cannot be given with
        /// declared columns.
        /// </summary>
        internal static ConditionRule ResolveConditionColumns(SdrfHeader header, IReadOnlyList<SdrfRow> rows,
            IReadOnlyList<string>? declared, string? singleCondition, bool unknownIsOwnCondition,
            List<string> refusals, List<string> notes)
        {
            var factorColumns = header
                .Where(name => name.StartsWith(SdrfLabelFreeDesign.FactorValuePrefix, StringComparison.Ordinal))
                .Distinct(StringComparer.Ordinal)
                .ToList();
            string factorList = string.Join(", ", factorColumns.Select(c => $"'{c}'"));
            string available = factorColumns.Count == 0
                ? "the SDRF has no factor value columns"
                : "its factor value columns are " + factorList;
            bool declaredAny = declared != null && declared.Count > 0;

            if (singleCondition != null)
            {
                string name = singleCondition.Trim();
                if (name.Length == 0)
                {
                    refusals.Add("The single condition is blank. Name the one condition every sample belongs to.");
                    return new ConditionRule(new List<string>(), ConditionSource.None, null, unknownIsOwnCondition);
                }

                if (declaredAny)
                {
                    refusals.Add($"Condition columns ({string.Join(", ", declared!.Select(c => $"'{c}'"))}) and a single condition " +
                                 $"('{name}') were both given. Give one: the columns build the condition from the SDRF, and a " +
                                 "single condition makes every sample one condition.");
                    return new ConditionRule(new List<string>(), ConditionSource.None, null, unknownIsOwnCondition);
                }

                if (factorColumns.Count > 0)
                    notes.Add($"Every sample is condition '{name}', as stated, so the SDRF's factor value column(s) {factorList} " +
                              "are set aside.");
                return new ConditionRule(new List<string>(), ConditionSource.Stated, name, unknownIsOwnCondition);
            }

            if (declaredAny)
            {
                var missing = declared!.Where(c => !header.Contains(c)).ToList();
                foreach (var column in missing)
                    refusals.Add($"The declared condition column '{column}' is not in the SDRF; {available}.");

                var repeated = declared!.GroupBy(c => c, StringComparer.Ordinal).Where(g => g.Count() > 1).Select(g => g.Key);
                foreach (var column in repeated)
                    refusals.Add($"The condition column '{column}' is declared more than once.");

                return new ConditionRule(declared!.ToList(), ConditionSource.Declared, null, unknownIsOwnCondition);
            }

            if (factorColumns.Count == 1)
                return new ConditionRule(factorColumns, ConditionSource.OnlyColumn, null, unknownIsOwnCondition);

            if (factorColumns.Count > 1)
            {
                notes.Add($"No condition column was declared, so the condition is built from every factor value column, in the " +
                          $"SDRF's order: {string.Join(" + ", factorColumns.Select(c => $"'{c}'"))}. Declare the column(s) to " +
                          "build it from fewer.");
                return new ConditionRule(factorColumns, ConditionSource.EveryColumn, null, unknownIsOwnCondition);
            }

            refusals.Add("No condition: the SDRF has no factor value columns and none was declared." +
                         DifferingCharacteristics(header, rows) +
                         " Declare the column(s) the condition is built from, or state that every sample is one condition " +
                         "(SingleCondition).");
            return new ConditionRule(new List<string>(), ConditionSource.None, null, unknownIsOwnCondition);
        }

        /// <summary>
        /// The condition-like characteristics whose stated values differ between the rows, for the no-condition refusal
        /// (QP-S42 b): " These characteristics differ between the samples, so they may be what tells the groups apart:
        /// ...", or empty when none does. Unknown words and blanks are not values. A guide for the reader of the refusal
        /// only: the groups are never built from it.
        /// </summary>
        private static string DifferingCharacteristics(SdrfHeader header, IReadOnlyList<SdrfRow> rows)
        {
            var differing = new List<string>();
            foreach (var column in header
                         .Where(name => name.StartsWith("characteristics[", StringComparison.OrdinalIgnoreCase)
                                        && !ReplicateLevelCharacteristic.IsMatch(name))
                         .Distinct(StringComparer.Ordinal))
            {
                var values = rows
                    .Select(r => (r[column] ?? string.Empty).Trim())
                    .Where(v => v.Length > 0 && !IsUnknownWord(v))
                    .Distinct(StringComparer.OrdinalIgnoreCase)
                    .ToList();
                if (values.Count < 2)
                    continue;
                differing.Add($"'{column}' ({string.Join(", ", values.Take(3))}" +
                              (values.Count > 3 ? $", and {values.Count - 3} more" : string.Empty) + ")");
            }

            return differing.Count == 0
                ? string.Empty
                : " These characteristics differ between the samples, so they may be what tells the groups apart: " +
                  string.Join("; ", differing) + ".";
        }

        /// <summary>
        /// Reads one row's factor values for the condition (MAP-07); with a stated single condition, that condition is
        /// the row's only value. An empty or unknown value is a refusal (SDRF-P1), and the refusal names its remedy, which
        /// depends on how the column was chosen, unless the caller stated that such files are a condition of their own
        /// (QP-S44): then the value is kept as a level (a blank as <c>not available</c>) and recorded in
        /// <paramref name="keptUnknown"/> for one note. A pooled reference's unknown value (a pool of every condition says
        /// <c>not applicable</c>) is kept as written: a reference is in no comparison, so its condition pairs it with
        /// nothing, and a pool is never refused (GR-20).
        /// </summary>
        internal static List<string> ReadFactorValues(SdrfRow row, int line, string fileName, ConditionRule condition,
            List<string> refusals, List<(string Column, string Value, string FileName)> keptUnknown, bool isReference = false)
        {
            if (condition.SingleCondition != null)
                return new List<string> { condition.SingleCondition };

            var factorValues = new List<string>();
            foreach (var column in condition.Columns)
            {
                string value = (row[column] ?? string.Empty).Trim();
                if (value.Length == 0)
                {
                    if (condition.UnknownIsOwnCondition)
                    {
                        keptUnknown.Add((column, string.Empty, fileName));
                        value = BlankFactorValue;
                    }
                    else
                        refusals.Add($"Line {line}{Describe(fileName)}: '{column}' is empty. Fill it in, or state that files whose " +
                                     "factor nobody stated are a condition of their own (UnknownFactorIsOwnCondition).");
                }
                else if (IsUnknownWord(value) && !isReference)
                {
                    if (condition.UnknownIsOwnCondition)
                        keptUnknown.Add((column, value, fileName));
                    else
                        refusals.Add($"Line {line}{Describe(fileName)}: '{column}' is '{value}'. A condition cannot be built from an unknown " +
                                     condition.Source switch
                                     {
                                         ConditionSource.Declared =>
                                             "factor; fill it in, leave the column out of the declared condition columns to pool these rows, ",
                                         ConditionSource.EveryColumn =>
                                             "factor. It was used because none was declared, so every factor value column makes the " +
                                             "condition; fill it in, declare the column(s) the condition should be built from instead, ",
                                         _ =>
                                             "factor. It was used because it is the only factor value column and none was declared; " +
                                             "fill it in, declare the column(s) the condition should be built from instead, ",
                                     } +
                                     "or state that files whose factor nobody stated are a condition of their own (UnknownFactorIsOwnCondition).");
                }
                factorValues.Add(value);
            }

            return factorValues;
        }

        /// <summary>
        /// One note per column and value for the factor values kept as a level of their own, as the caller stated (QP-S44),
        /// naming the files, so a condition called <c>not available</c> is never a surprise in the output.
        /// </summary>
        internal static void NoteKeptUnknownFactorValues(List<(string Column, string Value, string FileName)> keptUnknown,
            List<string> notes)
        {
            foreach (var group in keptUnknown
                         .GroupBy(k => (k.Column, Value: k.Value.ToLowerInvariant()))
                         .OrderBy(g => g.Key.Column, StringComparer.Ordinal).ThenBy(g => g.Key.Value, StringComparer.Ordinal))
            {
                var files = group.Select(k => k.FileName).Distinct(StringComparer.Ordinal).ToList();
                string what = group.Key.Value.Length == 0
                    ? $"is empty for {files.Count} file(s), written '{BlankFactorValue}' and"
                    : $"is '{group.First().Value}' for {files.Count} file(s),";
                notes.Add($"'{group.Key.Column}' {what} kept as a level of its own, as stated, so those files make a condition of " +
                          $"their own: {string.Join(", ", files.Take(10).Select(f => $"'{f}'"))}" +
                          (files.Count > 10 ? $", and {files.Count - 10} more." : "."));
            }
        }

        // QP-S45 (ruled 2026-11-11): the words before a number that make a file name read like a fraction, or like a repeat
        // injection, of one sample. sdrf's census (2026-11-07) used the same families; 'b' alone is left out, being as often
        // a batch or a biological replicate as a band.
        private static readonly HashSet<string> FractionWords = new(StringComparer.OrdinalIgnoreCase)
            { "f", "fr", "fx", "frac", "fraction", "band", "slice", "gel", "scx", "sax", "hprp", "hph", "hilic", "sec", "ief" };
        private static readonly HashSet<string> RepeatInjectionWords = new(StringComparer.OrdinalIgnoreCase)
            { "r", "rep", "replicate", "run", "inj", "injection", "tr", "tech", "techrep" };
        private static readonly Regex NameToken = new(@"\d+|[A-Za-z]+|[^A-Za-z\d]+", RegexOptions.CultureInvariant);

        /// <summary>
        /// A hint for files refused for sharing one condition, replicate, fraction and technical replicate (QP-S45, ruled
        /// 2026-11-11): when their names differ only by one number, and the word before it reads like a fraction (F02,
        /// F03) or a repeat injection (rep1, rep2), or the number ends the name after a hyphen (-01, -02), the hint says
        /// so and names the column that tells them apart. Null otherwise. A guess from names, so it is only ever a hint,
        /// never a renumbering.
        /// </summary>
        internal static string? SharedReplicateHint(IReadOnlyList<string> fileNames)
        {
            var tokens = fileNames
                .Select(f => NameToken.Matches(Path.GetFileNameWithoutExtension(f)).Select(m => m.Value).ToList())
                .ToList();
            if (tokens.Count < 2 || tokens.Any(t => t.Count != tokens[0].Count))
                return null;

            int differing = -1;
            for (int k = 0; k < tokens[0].Count; k++)
            {
                if (tokens.All(t => string.Equals(t[k], tokens[0][k], StringComparison.OrdinalIgnoreCase)))
                    continue;
                if (differing >= 0 || !tokens.All(t => char.IsDigit(t[k][0])))
                    return null;
                differing = k;
            }
            if (differing < 0)
                return null;
            var numbers = tokens.Select(t => t[differing].TrimStart('0')).ToList();
            if (numbers.Distinct(StringComparer.Ordinal).Count() != numbers.Count)
                return null;

            // The word right before the number, directly (F02) or across one '_' or '-' (rep_1, X-01).
            var first = tokens[0];
            string separator = differing >= 1 && !char.IsLetterOrDigit(first[differing - 1][0]) ? first[differing - 1] : string.Empty;
            int wordAt = separator.Length == 0 ? differing - 1 : differing - 2;
            string word = wordAt >= 0 && char.IsLetter(first[wordAt][0]) && (separator.Length == 0 || separator is "_" or "-")
                ? first[wordAt]
                : string.Empty;
            int digits = tokens.Max(t => t[differing].Length);

            string column;
            string kind;
            string shownWord;
            if (FractionWords.Contains(word) && digits <= 3)
                (column, kind, shownWord) = (SdrfLabelFreeDesign.FractionColumn, "a fraction", word);
            else if (RepeatInjectionWords.Contains(word) && digits <= 2)
                (column, kind, shownWord) = (SdrfLabelFreeDesign.TechnicalReplicateColumn, "a repeat injection", word);
            else if (separator == "-" && differing == first.Count - 1 && digits <= 2)
                (column, kind, shownWord) = (SdrfLabelFreeDesign.TechnicalReplicateColumn, "a repeat injection", string.Empty);
            else
                return null;

            string examples = string.Join(", ", tokens.Take(3).Select(t =>
                $"'{(shownWord.Length > 0 ? t[wordAt] : string.Empty)}{separator}{t[differing]}'"));
            string plural = kind == "a fraction" ? "fractions" : "repeat injections";
            return $"Their names differ only by a number that reads like {kind} ({examples}): if they are {plural} of one sample, " +
                   $"give each its '{column}'.";
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

            // Split at either separator on any platform: a search run on Windows can be read on Linux.
            foreach (var group in usable.GroupBy(SdrfFileNamePattern.FileName, StringComparer.Ordinal))
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
