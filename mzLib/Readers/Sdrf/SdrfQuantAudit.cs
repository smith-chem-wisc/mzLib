using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;

namespace Readers
{
    /// <summary>What kind of quantification design, if any, an SDRF document actually describes.</summary>
    public enum SdrfQuantKind
    {
        /// <summary>No isobaric label anywhere: label-free, SILAC, dimethyl, and so on.</summary>
        NotIsobaric,

        /// <summary>
        /// An isobaric experiment whose <c>comment[label]</c> names the REAGENT KIT rather than the
        /// channel -- "TMT10plex" on every row. The channels are never enumerated, so no channel-level
        /// design exists in the file to recover. Detecting this is the point: read as if it were
        /// channel-level, every row of such a file claims the same channel.
        /// </summary>
        KitOnly,

        /// <summary>At least one row names an individual reporter channel, e.g. "TMT127N".</summary>
        ChannelLevel
    }

    /// <summary>How one <c>comment[label]</c> cell was written.</summary>
    public enum SdrfLabelForm
    {
        /// <summary>Free text, e.g. <c>TMT126</c>. About six sevenths of the corpus.</summary>
        Bare,

        /// <summary>
        /// The key=value grammar, e.g. <c>NT=TMT126;AC=PRIDE:0000516</c>. The key ORDER varies in the
        /// wild -- <c>AC=</c> first is commoner than <c>NT=</c> first for label-free -- so this is
        /// decided by parsing, never by matching a prefix.
        /// </summary>
        Accessioned
    }

    /// <summary>Where the plex a row belongs to came from, if anywhere.</summary>
    public enum SdrfPlexSource
    {
        /// <summary>Nothing in the document says. This is the usual answer for bulk TMT.</summary>
        None,

        /// <summary>A column named it. The column is reported, because its spelling varies.</summary>
        Column,

        /// <summary>Recovered from a pattern in the data file names, which is reported for inspection.</summary>
        FileName
    }

    /// <summary>Whether one of the eleven facts survived into this document.</summary>
    public enum SdrfFactPresence
    {
        /// <summary>At least one row carries a real answer.</summary>
        Present,

        /// <summary>
        /// The column is missing, or every cell is empty or an SDRF reserved word. Reserved words
        /// count as absent deliberately: "not available" is the spec-correct way to say nothing, and
        /// a validator is right to accept it, but a consumer still gets no answer.
        /// </summary>
        Absent,

        /// <summary>
        /// Answered, but not in a form this fact can be read in -- a fraction identifier of "19b"
        /// where an integer is required. Distinguished from Absent because the two need different
        /// fixes: one is a curation gap, the other a curation error.
        /// </summary>
        Unparseable
    }

    /// <summary>One of the eleven per-(data file x channel) facts, and how this document fared on it.</summary>
    /// <param name="Number">Its number in the project's fact set, so a report can be read beside the plan.</param>
    /// <param name="Name">Short name of the fact.</param>
    /// <param name="Column">The SDRF column it lives in, or null when SDRF has no column for it.</param>
    /// <param name="Presence">Present, Absent or Unparseable.</param>
    /// <param name="Detail">Why, when that is not obvious -- the offending value, or what is missing.</param>
    public sealed record SdrfQuantFact(
        int Number,
        string Name,
        string? Column,
        SdrfFactPresence Presence,
        string? Detail = null);

    /// <summary>
    /// A read-only audit of what one SDRF document can and cannot tell a quantification engine.
    ///
    /// This exists because the failure it describes is silent. A file whose <c>comment[label]</c>
    /// names a kit rather than a channel parses cleanly, validates cleanly, and then yields one
    /// channel for a ten-channel experiment. Likewise a file with two rows per (data file, channel),
    /// or with no isobaric modification declared at all. None of that is a structural defect, so
    /// <see cref="SdrfValidator"/> is right not to report it; it is a defect in what the annotation
    /// MEANS, and it needs its own reader.
    ///
    /// Deliberately scoped to every SDRF rather than to TMT. A user points this at whatever they just
    /// downloaded, and "this is not an isobaric experiment" is a useful answer.
    /// </summary>
    public sealed record SdrfQuantAudit
    {
        /// <summary>The document audited, for reports that pool many.</summary>
        public string? FilePath { get; init; }

        /// <summary>Channel-level, kit-only, or not isobaric at all.</summary>
        public SdrfQuantKind Kind { get; init; }

        /// <summary>
        /// The distinct reporter channels named, upper-cased so that <c>itraq114</c> and
        /// <c>ITRAQ114</c> are one channel. The spellings actually used are in
        /// <see cref="NonCanonicalChannelSpellings"/>.
        /// </summary>
        public IReadOnlyList<string> Channels { get; init; } = [];

        /// <summary>
        /// How many channels the document enumerates, which is the plex size it implies. Note this is
        /// what the FILE says, not what the kit provides: a 10-plex run that used 8 channels implies 8.
        /// </summary>
        public int ImpliedPlex => Channels.Count;

        /// <summary>Kit names found, e.g. <c>TMT10PLEX</c>. Populated for <see cref="SdrfQuantKind.KitOnly"/>.</summary>
        public IReadOnlyList<string> KitLabels { get; init; } = [];

        /// <summary>
        /// The reagent families named, upper-cased: <c>TMT</c> or <c>ITRAQ</c>. TMTpro folds into TMT,
        /// being the same reporter series extended.
        ///
        /// Worth reporting separately from the channels because it is what makes counts of this corpus
        /// comparable. Three earlier surveys of it disagreed -- 49, 51 and 52 channel-level files --
        /// and the spread turns out to be scope, not arithmetic: 51 channel-level and 31 kit-only is
        /// right for the TMT family alone, and iTRAQ adds 10 and 6 on top.
        /// </summary>
        public IReadOnlyList<string> Families { get; init; } = [];

        /// <summary>Where the plex assignment came from.</summary>
        public SdrfPlexSource PlexSource { get; init; }

        /// <summary>
        /// The evidence for <see cref="PlexSource"/> -- the column name, or the filename pattern --
        /// shown rather than summarised, because every inference in this corpus has needed checking.
        /// </summary>
        public string? PlexEvidence { get; init; }

        /// <summary>The distinct plex values recovered, in document order.</summary>
        public IReadOnlyList<string> Plexes { get; init; } = [];

        /// <summary>
        /// Whether any <c>comment[modification parameters]</c> cell declares a TMT or iTRAQ
        /// modification. A channel-level file that declares none cannot be searched from its own
        /// annotation, and this is common enough to be worth its own line in the report.
        /// </summary>
        public bool DeclaresIsobaricModification { get; init; }

        /// <summary>The eleven facts, in order.</summary>
        public IReadOnlyList<SdrfQuantFact> Facts { get; init; } = [];

        /// <summary>How many label cells were written as free text.</summary>
        public int BareLabelCells { get; init; }

        /// <summary>How many label cells were written in the key=value grammar.</summary>
        public int AccessionedLabelCells { get; init; }

        /// <summary>
        /// True when one document uses both forms. Exactly one file in the curated corpus does, which
        /// is the whole reason the form is decided per cell and never sniffed once per document.
        /// </summary>
        public bool MixesLabelForms => BareLabelCells > 0 && AccessionedLabelCells > 0;

        /// <summary>
        /// Channel spellings that are not the canonical upper-case form -- <c>itraq114</c>,
        /// <c>TMT10Plex</c>, <c>tmt126</c>. Harmless to a case-insensitive reader and a real finding
        /// for anyone joining on the raw string.
        /// </summary>
        public IReadOnlyList<string> NonCanonicalChannelSpellings { get; init; } = [];

        /// <summary>
        /// (data file, label) pairs occurring on more than one row, worst first. A channel is one
        /// measurement of one file, so a repeat means the pair does not identify a row and any
        /// dictionary keyed on it silently keeps one. For a kit-only document this degenerates into
        /// "this data file appears more than once", which is the same defect one level up.
        /// </summary>
        public IReadOnlyList<SdrfDuplicateRows> DuplicateFileLabelPairs { get; init; } = [];

        /// <summary>The largest number of rows any single data file occupies.</summary>
        public int MostRowsForOneDataFile { get; init; }

        /// <summary>Rows whose cell count differs from the header's.</summary>
        public int RaggedRows { get; init; }

        /// <summary>Data files named by the document that were found on disk. Empty when no directory was given.</summary>
        public IReadOnlyList<string> MatchedDataFiles { get; init; } = [];

        /// <summary>Data files named by the document that were not found. Empty when no directory was given.</summary>
        public IReadOnlyList<string> UnmatchedDataFiles { get; init; } = [];

        /// <summary>True when a data directory was searched at all, so "0 matched" can be told from "not looked".</summary>
        public bool SearchedForDataFiles { get; init; }

        /// <summary>A human-readable report, which is what the CMD verb prints.</summary>
        public string ToReport()
        {
            var text = new StringBuilder();
            text.AppendLine($"SDRF quantification audit: {FilePath ?? "(document)"}");
            text.AppendLine($"  design            : {Kind}");

            if (Kind == SdrfQuantKind.ChannelLevel)
                text.AppendLine($"  channels          : {Channels.Count} ({string.Join(", ", Channels)})");
            if (KitLabels.Count > 0)
                text.AppendLine($"  kit               : {string.Join(", ", KitLabels)} -- channels are never enumerated");

            text.AppendLine(PlexSource == SdrfPlexSource.None
                ? "  plex              : none stated"
                : $"  plex              : {PlexSource} -- {PlexEvidence} -> {string.Join(", ", Plexes)}");

            if (Kind != SdrfQuantKind.NotIsobaric)
                text.AppendLine($"  isobaric mod      : {(DeclaresIsobaricModification ? "declared" : "NOT DECLARED")}");

            text.AppendLine($"  label form        : {BareLabelCells} bare, {AccessionedLabelCells} accessioned"
                            + (MixesLabelForms ? " -- MIXED within this document" : string.Empty));

            if (NonCanonicalChannelSpellings.Count > 0)
                text.AppendLine($"  spelling variants : {string.Join(", ", NonCanonicalChannelSpellings)}");
            if (DuplicateFileLabelPairs.Count > 0)
                text.AppendLine($"  duplicate rows    : {DuplicateFileLabelPairs.Count} (data file, label) pairs repeat; "
                                + $"worst {DuplicateFileLabelPairs[0].Rows} rows for {DuplicateFileLabelPairs[0].DataFile} / {DuplicateFileLabelPairs[0].Label}");
            if (RaggedRows > 0)
                text.AppendLine($"  ragged rows       : {RaggedRows}");
            if (SearchedForDataFiles)
                text.AppendLine($"  data files        : {MatchedDataFiles.Count} found, {UnmatchedDataFiles.Count} missing");

            text.AppendLine("  facts:");
            foreach (var fact in Facts)
            {
                string mark = fact.Presence switch
                {
                    SdrfFactPresence.Present => "  ok  ",
                    SdrfFactPresence.Absent => " ABSENT",
                    _ => " UNPARSEABLE"
                };
                text.AppendLine($"    {fact.Number,2}. {fact.Name,-22}{mark}"
                                + (fact.Detail == null ? string.Empty : $"  {fact.Detail}"));
            }
            return text.ToString();
        }
    }

    /// <summary>One repeated (data file, label) pair and how many rows carry it.</summary>
    public sealed record SdrfDuplicateRows(string DataFile, string Label, int Rows);

    /// <summary>
    /// Reads an <see cref="SdrfDocument"/> for what it says about quantification. Read-only, and it
    /// never throws on a bad document: the point is to describe one.
    /// </summary>
    public static class SdrfQuantAuditor
    {
        // A reporter channel: family, nominal mass, and an optional N/C position with an optional
        // deuterium marker ("TMT131ND"). Case-insensitive because the corpus writes itraq114,
        // ITRAQ114 and TMT10Plex, and those are the same reagents as their canonical spellings.
        private static readonly Regex ChannelLabel =
            new(@"^(?<family>tmtpro|tmt|itraq)(?<mass>\d{3})(?<position>[nc]?d?)$",
                RegexOptions.IgnoreCase | RegexOptions.Compiled);

        // A kit: "TMT10plex", "TMTpro-16plex", "iTRAQ 4plex".
        private static readonly Regex KitLabel =
            new(@"^(?<family>tmtpro|tmt|itraq)[-_ ]?(?<size>\d+)\s*plex$",
                RegexOptions.IgnoreCase | RegexOptions.Compiled);

        // The family name alone. Enumerates nothing, so it is a kit-level statement, not a channel.
        private static readonly Regex FamilyOnlyLabel =
            new(@"^(tmtpro|tmt|itraq)$", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        // Column spellings observed to carry a plex or batch. comment[batch] is a misspelling of the
        // characteristics form and is accepted because a real corpus file uses it.
        private static readonly string[] PlexColumns =
        [
            "characteristics[biological replicate batch]",
            "characteristics[batch]",
            "comment[batch]",
            "characteristics[plex]",
            "comment[plex]",
            "characteristics[tmt plex]"
        ];

        private const string LabelColumn = "comment[label]";
        private const string DataFileColumn = "comment[data file]";
        private const string ModificationColumn = "comment[modification parameters]";

        private static readonly HashSet<string> ReservedWords = new(StringComparer.OrdinalIgnoreCase)
        {
            "not available", "not applicable", "anonymized", "pooled"
        };

        /// <summary>Audits the document at <paramref name="sdrfPath"/>.</summary>
        /// <param name="sdrfPath">The SDRF file to read.</param>
        /// <param name="dataDirectory">
        /// Optional folder of downloaded data files. When given, every <c>comment[data file]</c> is
        /// looked for in it and reported found or missing.
        /// </param>
        public static SdrfQuantAudit Audit(string sdrfPath, string? dataDirectory = null)
        {
            var document = new SdrfDocument(sdrfPath);
            document.LoadResults();
            return Audit(document, dataDirectory);
        }

        /// <summary>Audits a document already in memory.</summary>
        public static SdrfQuantAudit Audit(SdrfDocument document, string? dataDirectory = null)
        {
            ArgumentNullException.ThrowIfNull(document);

            var rows = document.Results ?? [];
            var header = document.Header;

            var channels = new List<string>();
            var kits = new List<string>();
            var families = new List<string>();
            var spellings = new List<string>();
            var pairCounts = new Dictionary<(string File, string Label), int>();
            var rowsPerFile = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);
            int bare = 0, accessioned = 0;

            foreach (var row in rows)
            {
                string? dataFile = Value(row, DataFileColumn);
                if (dataFile != null)
                    rowsPerFile[dataFile] = rowsPerFile.GetValueOrDefault(dataFile) + 1;

                foreach (string cell in row.All(LabelColumn))
                {
                    if (string.IsNullOrWhiteSpace(cell)) continue;

                    (string name, SdrfLabelForm form) = ReadLabel(cell);
                    if (form == SdrfLabelForm.Accessioned) accessioned++; else bare++;

                    var channel = ChannelLabel.Match(name);
                    if (channel.Success)
                    {
                        AddFamily(families, channel.Groups["family"].Value);
                        string canonical = name.ToUpperInvariant();
                        if (!channels.Contains(canonical)) channels.Add(canonical);
                        if (!string.Equals(name, canonical, StringComparison.Ordinal)
                            && !spellings.Contains(name))
                        {
                            spellings.Add(name);
                        }
                    }
                    else if (KitLabel.IsMatch(name) || FamilyOnlyLabel.IsMatch(name))
                    {
                        var kit = KitLabel.Match(name);
                        AddFamily(families, kit.Success ? kit.Groups["family"].Value : name);
                        string canonical = name.ToUpperInvariant();
                        if (!kits.Contains(canonical)) kits.Add(canonical);
                        if (!string.Equals(name, canonical, StringComparison.Ordinal)
                            && !spellings.Contains(name))
                        {
                            spellings.Add(name);
                        }
                    }
                    else
                    {
                        continue; // label free, SILAC, dimethyl: not an isobaric statement
                    }

                    if (dataFile != null)
                    {
                        var key = (dataFile, name.ToUpperInvariant());
                        pairCounts[key] = pairCounts.GetValueOrDefault(key) + 1;
                    }
                }
            }

            var kind = channels.Count > 0 ? SdrfQuantKind.ChannelLevel
                : kits.Count > 0 ? SdrfQuantKind.KitOnly
                : SdrfQuantKind.NotIsobaric;

            (SdrfPlexSource plexSource, string? evidence, var plexes) = ReadPlex(header, rows);

            var duplicates = pairCounts
                .Where(kv => kv.Value > 1)
                .OrderByDescending(kv => kv.Value)
                .Select(kv => new SdrfDuplicateRows(kv.Key.File, kv.Key.Label, kv.Value))
                .ToList();

            (var matched, var unmatched) = MatchDataFiles(rowsPerFile.Keys, dataDirectory);

            return new SdrfQuantAudit
            {
                FilePath = document.FilePath,
                Kind = kind,
                Channels = channels,
                KitLabels = kits,
                Families = families,
                PlexSource = plexSource,
                PlexEvidence = evidence,
                Plexes = plexes,
                DeclaresIsobaricModification = DeclaresIsobaricModification(rows),
                Facts = ReadFacts(header, rows, kind),
                BareLabelCells = bare,
                AccessionedLabelCells = accessioned,
                NonCanonicalChannelSpellings = spellings,
                DuplicateFileLabelPairs = duplicates,
                MostRowsForOneDataFile = rowsPerFile.Count == 0 ? 0 : rowsPerFile.Values.Max(),
                RaggedRows = rows.Count(r => r.Cells.Count != header.Count),
                MatchedDataFiles = matched,
                UnmatchedDataFiles = unmatched,
                SearchedForDataFiles = dataDirectory != null
            };
        }

        /// <summary>TMTpro folds into TMT: the same reporter series, extended.</summary>
        private static void AddFamily(List<string> families, string family)
        {
            string canonical = family.StartsWith("tmt", StringComparison.OrdinalIgnoreCase) ? "TMT" : "ITRAQ";
            if (!families.Contains(canonical)) families.Add(canonical);
        }

        /// <summary>
        /// The name a label cell carries, and the form it was written in. Both key orders occur --
        /// <c>NT=..;AC=..</c> and <c>AC=..;NT=..</c> -- so the cell is parsed rather than prefix-matched.
        /// </summary>
        private static (string Name, SdrfLabelForm Form) ReadLabel(string cell)
        {
            var pairs = SdrfCell.ParseKeyValues(cell);
            if (pairs.Count == 0) return (cell.Trim(), SdrfLabelForm.Bare);
            return pairs.TryGetValue("NT", out string? name) && !string.IsNullOrWhiteSpace(name)
                ? (name.Trim(), SdrfLabelForm.Accessioned)
                : (cell.Trim(), SdrfLabelForm.Accessioned);
        }

        private static (SdrfPlexSource, string?, IReadOnlyList<string>) ReadPlex(
            SdrfHeader header, IReadOnlyList<SdrfRow> rows)
        {
            foreach (string column in PlexColumns)
            {
                if (!header.Contains(column)) continue;
                var values = rows.Select(r => Value(r, column)).Where(v => v != null).Distinct().ToList();
                if (values.Count > 0)
                    return (SdrfPlexSource.Column, column, values!);
            }

            // No filename inference is attempted. The project measured the source-name partition at
            // 3 correct out of 9 where ground truth exists, failing by OVER-splitting, and no bulk
            // dataset in the curated corpus carries a usable batch column. Reporting "none stated" is
            // the honest answer; guessing here would manufacture plexes that are not in the file.
            return (SdrfPlexSource.None, null, []);
        }

        private static bool DeclaresIsobaricModification(IReadOnlyList<SdrfRow> rows)
        {
            foreach (var row in rows)
            {
                foreach (string cell in row.All(ModificationColumn))
                {
                    if (string.IsNullOrWhiteSpace(cell)) continue;
                    string name = ReadLabel(cell).Name;
                    if (name.Contains("tmt", StringComparison.OrdinalIgnoreCase)
                        || name.Contains("itraq", StringComparison.OrdinalIgnoreCase))
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        private static IReadOnlyList<SdrfQuantFact> ReadFacts(
            SdrfHeader header, IReadOnlyList<SdrfRow> rows, SdrfQuantKind kind)
        {
            var facts = new List<SdrfQuantFact>
            {
                Simple(1, "data file", DataFileColumn),
                LabelFact(2, kind),
                new(3, "reporter m/z", null, SdrfFactPresence.Absent,
                    "SDRF has no column for it; it is derived from the channel"),
                Simple(4, "sample", "source name"),
                PlexFact(5),
                FactorValue(6),
                Integer(7, "biological replicate", "characteristics[biological replicate]"),
                Integer(8, "technical replicate", "comment[technical replicate]"),
                Integer(9, "fraction", "comment[fraction identifier]"),
                Simple(10, "sample type", "characteristics[sample type]"),
                Simple(11, "assay / run", "assay name")
            };
            return facts;

            SdrfQuantFact Simple(int number, string name, string column) =>
                new(number, name, column,
                    header.Contains(column)
                        ? (rows.Any(r => Value(r, column) != null) ? SdrfFactPresence.Present : SdrfFactPresence.Absent)
                        : SdrfFactPresence.Absent,
                    header.Contains(column) ? null : "column absent");

            SdrfQuantFact LabelFact(int number, SdrfQuantKind k) => k switch
            {
                SdrfQuantKind.ChannelLevel => new(number, "channel", LabelColumn, SdrfFactPresence.Present),
                SdrfQuantKind.KitOnly => new(number, "channel", LabelColumn, SdrfFactPresence.Unparseable,
                    "names the kit, not the channel"),
                _ => new(number, "channel", LabelColumn, SdrfFactPresence.Absent, "not an isobaric experiment")
            };

            SdrfQuantFact PlexFact(int number)
            {
                string? found = PlexColumns.FirstOrDefault(header.Contains);
                return found == null
                    ? new(number, "plex", null, SdrfFactPresence.Absent, "no column names it")
                    : new(number, "plex", found,
                        rows.Any(r => Value(r, found) != null) ? SdrfFactPresence.Present : SdrfFactPresence.Absent);
            }

            SdrfQuantFact FactorValue(int number)
            {
                string? column = header.FirstOrDefault(h =>
                    h.StartsWith("factor value[", StringComparison.OrdinalIgnoreCase));
                return column == null
                    ? new(number, "condition", null, SdrfFactPresence.Absent, "no factor value column")
                    : new(number, "condition", column,
                        rows.Any(r => Value(r, column) != null) ? SdrfFactPresence.Present : SdrfFactPresence.Absent);
            }

            SdrfQuantFact Integer(int number, string name, string column)
            {
                if (!header.Contains(column))
                    return new(number, name, column, SdrfFactPresence.Absent, "column absent");

                var answers = rows.Select(r => Value(r, column)).Where(v => v != null).ToList();
                if (answers.Count == 0)
                    return new(number, name, column, SdrfFactPresence.Absent);

                string? bad = answers.FirstOrDefault(v =>
                    !int.TryParse(v, NumberStyles.Integer, CultureInfo.InvariantCulture, out _));
                return bad == null
                    ? new(number, name, column, SdrfFactPresence.Present)
                    : new(number, name, column, SdrfFactPresence.Unparseable, $"not an integer: '{bad}'");
            }
        }

        /// <summary>
        /// The value under a column, or null when the column is absent, the row is short, the cell is
        /// blank, or the cell is a reserved word. Reserved words are spec-correct ways to say nothing,
        /// so they are absence, not answers -- see <see cref="SdrfFactPresence.Absent"/>.
        /// </summary>
        private static string? Value(SdrfRow row, string column)
        {
            string? raw = row[column];
            if (string.IsNullOrWhiteSpace(raw)) return null;
            raw = raw.Trim();
            return ReservedWords.Contains(raw) ? null : raw;
        }

        private static (IReadOnlyList<string> Matched, IReadOnlyList<string> Unmatched) MatchDataFiles(
            IEnumerable<string> dataFiles, string? dataDirectory)
        {
            if (dataDirectory == null || !Directory.Exists(dataDirectory))
                return ([], []);

            var onDisk = Directory
                .EnumerateFiles(dataDirectory, "*", SearchOption.AllDirectories)
                .Select(Path.GetFileName)
                .Where(n => n != null)
                .ToHashSet(StringComparer.OrdinalIgnoreCase);

            var matched = new List<string>();
            var unmatched = new List<string>();
            foreach (string named in dataFiles.OrderBy(f => f, StringComparer.OrdinalIgnoreCase))
            {
                // The column may carry a path or a bare name; the file on disk is matched by name,
                // because a downloaded copy rarely sits where the submitter's did.
                string name = Path.GetFileName(named.Replace('\\', '/'));
                (onDisk.Contains(name) ? matched : unmatched).Add(named);
            }
            return (matched, unmatched);
        }
    }
}
