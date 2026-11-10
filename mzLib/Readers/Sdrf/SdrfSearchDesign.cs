using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;

namespace Readers
{
    /// <summary>
    /// What a caller tells <see cref="SdrfSearchDesign.Read(SdrfDocument, IEnumerable{string}, SdrfSearchDesignOptions?)"/>
    /// that the SDRF alone cannot say: whether the search is isobaric, with which kit, where each file's plex is
    /// written, and which columns make up the condition.
    /// </summary>
    public sealed class SdrfSearchDesignOptions
    {
        /// <summary>
        /// The kit an isobaric search used, or null for a label-free search. It decides which reader is used:
        /// <see cref="SdrfIsobaricDesign"/> with a kit, <see cref="SdrfLabelFreeDesign"/> without. MetaMorpheus
        /// passes the tag of <c>SearchParameters.MultiplexModId</c>.
        /// </summary>
        public IsobaricMassTag? Tag { get; init; }

        /// <summary>Isobaric only: as <see cref="SdrfIsobaricDesignOptions.PlexColumn"/>.</summary>
        public string? PlexColumn { get; init; }

        /// <summary>Isobaric only: as <see cref="SdrfIsobaricDesignOptions.PlexFileNamePattern"/>.</summary>
        public string? PlexFileNamePattern { get; init; }

        /// <summary>Isobaric only: as <see cref="SdrfIsobaricDesignOptions.SinglePlex"/>.</summary>
        public string? SinglePlex { get; init; }

        /// <summary>
        /// The columns the condition is built from, in order. A bare name is its factor value column:
        /// <c>genotype</c> means <c>factor value[genotype]</c>, as MetaMorpheus's <c>--sdrfCondition</c> reads it.
        /// The <c>factor value[</c> prefix is matched in any case and written as the specification names it. Any
        /// other name holding a <c>[</c>, such as <c>characteristics[cell type]</c>, is used as given. Otherwise
        /// as <see cref="SdrfLabelFreeDesignOptions.ConditionColumns"/>: leave null or empty to use the document's
        /// only factor column.
        /// </summary>
        public IReadOnlyList<string>? ConditionColumns { get; init; }

        /// <summary>
        /// Optional: searched file (as given) -> the acquired file it came from, when the caller knows it. Passed to
        /// <see cref="SdrfSearchScope.Restrict"/>, where it overrides the suffix rule for those files.
        /// </summary>
        public IReadOnlyDictionary<string, string>? AcquiredNameOf { get; init; }
    }

    /// <summary>
    /// The experimental design ONE SEARCH will quantify with, made from an SDRF: either the design MetaMorpheus will
    /// accept for exactly these files, or every reason there is none.
    ///
    /// <para><b>Why one entry point.</b> Two callers need the same answer. A pipeline (PXReprise) asks it before a
    /// search starts, so that a search whose SDRF gives no usable design is never run: an invalid design makes
    /// MetaMorpheus skip quantification after the whole search. MetaMorpheus asks it when it writes the design from
    /// the SDRF. If the two asked different code they could disagree, and a search that passed the check would still
    /// quantify nothing.</para>
    ///
    /// <para><b>What it does.</b> Restricts the SDRF to the searched files with <see cref="SdrfSearchScope.Restrict"/>,
    /// so a row naming <c>X.raw</c> stands for a search of <c>X.mzML</c> or <c>X-calib.mzML</c>, then reads the kept
    /// rows with <see cref="SdrfLabelFreeDesign"/> or <see cref="SdrfIsobaricDesign"/>. Every rule is theirs: this
    /// type adds none, so a design it accepts is one those readers accept. An SDRF with no <c>comment[data file]</c>
    /// column cannot be joined, so it is read as given and the reader keys on its <c>comment[searched data file]</c>
    /// or refuses it.</para>
    ///
    /// <para><b>Never throws for a bad SDRF or a bad combination of options.</b> A missing, empty or unreadable file,
    /// a join that keeps no row, and a plex source given without a kit are refusals, like the readers' own. Only a
    /// null argument throws. Nothing is written unless <see cref="WriteDesign"/> is called, and it throws on a refused
    /// design, so a refusal can never leave a design file behind.</para>
    /// </summary>
    public sealed class SdrfSearchDesign
    {
        private const string FactorValuePrefix = SdrfLabelFreeDesign.FactorValuePrefix;

        private readonly List<string> _refusals;
        private readonly List<string> _notes;
        private readonly List<string> _joinReport;

        private SdrfSearchDesign(SdrfLabelFreeDesign? labelFree, SdrfIsobaricDesign? isobaric, SdrfSearchScoping? scoping,
            List<string> refusals, List<string> notes, List<string> joinReport)
        {
            LabelFree = labelFree;
            Isobaric = isobaric;
            Scoping = scoping;
            _refusals = refusals;
            _notes = notes;
            _joinReport = joinReport;
        }

        /// <summary>True when there is nothing in <see cref="Refusals"/>: MetaMorpheus will accept this design.</summary>
        public bool IsValid => Refusals.Count == 0;

        /// <summary>The label-free reader's result, or null for an isobaric search or when no row was left to read.</summary>
        public SdrfLabelFreeDesign? LabelFree { get; }

        /// <summary>The isobaric reader's result, or null for a label-free search or when no row was left to read.</summary>
        public SdrfIsobaricDesign? Isobaric { get; }

        /// <summary>
        /// The SDRF restricted to the searched files, with what the restriction could not place, or null when the
        /// SDRF could not be read or has no <c>comment[data file]</c> column.
        /// </summary>
        public SdrfSearchScoping? Scoping { get; }

        /// <summary>Every reason there is no design: this type's own, then the reader's. Empty when it is valid.</summary>
        public IReadOnlyList<string> Refusals =>
            _refusals.Concat(LabelFree?.Refusals ?? Isobaric?.Refusals ?? Array.Empty<string>()).ToList();

        /// <summary>
        /// What is worth checking on the way to the design: rows set aside because the search does not read their
        /// files, then the reader's own notes (replicate numbers added, gaps in the numbering, and so on).
        /// </summary>
        public IReadOnlyList<string> Notes =>
            _notes.Concat(LabelFree?.Notes ?? Isobaric?.Notes ?? Array.Empty<string>()).ToList();

        /// <summary>
        /// A human-readable account: how the SDRF was joined to the searched files, then the reader's report, or this
        /// type's own refusals when no reader ran. Print it wherever the design is checked or built. Its wording is
        /// for people; test <see cref="IsValid"/> and <see cref="Refusals"/>, never the text.
        /// </summary>
        public string Report()
        {
            var report = new StringBuilder();
            if (LabelFree is null && Isobaric is null)
                report.AppendLine($"Design REFUSED: {Refusals.Count} reason(s). No design was produced.");
            foreach (var line in _joinReport)
                report.AppendLine(line);
            foreach (var note in _notes)
                report.AppendLine("  note: " + note);
            if (LabelFree is not null)
                report.Append(LabelFree.Report());
            else if (Isobaric is not null)
                report.Append(Isobaric.Report());
            else
                foreach (var refusal in _refusals)
                    report.AppendLine("  refused: " + refusal);
            return report.ToString();
        }

        /// <summary>
        /// The design as mzLib quantification takes it. Throws <see cref="InvalidOperationException"/> when the
        /// design was refused.
        /// </summary>
        public SampleExperimentalDesign ToExperimentalDesign()
        {
            ThrowIfRefused();
            return Isobaric?.ToExperimentalDesign() ?? LabelFree!.ToExperimentalDesign();
        }

        /// <summary>
        /// Writes the design file MetaMorpheus reads: <c>ExperimentalDesign.tsv</c>'s format for a label-free search,
        /// <c>TmtDesign.txt</c>'s for an isobaric one, to <paramref name="path"/>. Throws
        /// <see cref="InvalidOperationException"/>, and writes nothing, when the design was refused.
        /// </summary>
        public void WriteDesign(string path)
        {
            if (path == null)
                throw new ArgumentNullException(nameof(path));
            ThrowIfRefused();
            if (Isobaric is not null)
                Isobaric.WriteTmtDesign(path);
            else
                LabelFree!.WriteExperimentalDesignTsv(path);
        }

        /// <summary>
        /// Reads an SDRF file for one search. A file that is missing, empty or cannot be read is refused, not thrown.
        /// See <see cref="Read(SdrfDocument, IEnumerable{string}, SdrfSearchDesignOptions?)"/>.
        /// </summary>
        public static SdrfSearchDesign Read(string sdrfPath, IEnumerable<string> searchedFiles, SdrfSearchDesignOptions? options = null)
        {
            if (sdrfPath == null)
                throw new ArgumentNullException(nameof(sdrfPath));
            if (searchedFiles == null)
                throw new ArgumentNullException(nameof(searchedFiles));

            if (!File.Exists(sdrfPath))
                return Refused($"The SDRF file does not exist: {sdrfPath}");

            SdrfDocument document;
            try
            {
                // Loaded once, here, so that a file that cannot be read is a refusal rather than an exception thrown
                // later from inside a reader. The copy in memory never reloads from disk.
                var onDisk = new SdrfDocument(sdrfPath);
                document = new SdrfDocument(onDisk.Header, onDisk.Results.ToList());
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException or MzLibException)
            {
                return Refused($"The SDRF file could not be read: {sdrfPath}. {e.Message}");
            }

            return Read(document, searchedFiles, options);
        }

        /// <summary>
        /// Restricts the SDRF to the searched files and reads the design those rows describe. Never throws for a bad
        /// document or a bad combination of options: every problem is in <see cref="Refusals"/>.
        /// </summary>
        /// <param name="sdrf">The SDRF.</param>
        /// <param name="searchedFiles">The files the search will read, as paths or bare names, exactly as the search
        /// names them (a calibrated search names its <c>-calib.mzML</c> files).</param>
        /// <param name="options">The kit, plex source and condition columns; null for a label-free search with the
        /// document's only factor column.</param>
        public static SdrfSearchDesign Read(SdrfDocument sdrf, IEnumerable<string> searchedFiles, SdrfSearchDesignOptions? options = null)
        {
            if (sdrf == null)
                throw new ArgumentNullException(nameof(sdrf));
            if (searchedFiles == null)
                throw new ArgumentNullException(nameof(searchedFiles));
            options ??= new SdrfSearchDesignOptions();

            var searched = searchedFiles.Where(f => !string.IsNullOrWhiteSpace(f)).ToList();
            if (searched.Count == 0)
                return Refused("No searched file was given, so there is no search for a design to describe.");

            if (options.Tag is null && (options.PlexColumn ?? options.PlexFileNamePattern ?? options.SinglePlex) is not null)
                return Refused("A plex source was given but no isobaric kit, so it is not known whether the search is isobaric. " +
                    "Give the kit for an isobaric search, or no plex source for a label-free one.");

            var conditionColumns = options.ConditionColumns?.Select(ConditionColumnName).ToList();

            var notes = new List<string>();
            var joinReport = new List<string>();
            SdrfSearchScoping? scoping = null;
            var scoped = sdrf;
            if (sdrf.Header.Contains(SdrfLabelFreeDesign.DataFileColumn))
            {
                scoping = SdrfSearchScope.Restrict(sdrf, searched, options.AcquiredNameOf);
                int rowsRead = sdrf.Results.Count;
                int rowsKept = scoping.Document.Results.Count;
                ReportJoin(scoping, rowsRead, rowsKept, joinReport);

                // When no row is kept the reader is not run: its "The SDRF has no rows" would point at the SDRF
                // rather than at the searched files.
                if (rowsKept == 0 && rowsRead > 0)
                    return new SdrfSearchDesign(null, null, scoping,
                        new List<string> { "No SDRF row names a searched file, so there is no design to read. " +
                            "Check that the searched files are the ones this SDRF describes." },
                        notes, joinReport);

                // The reader sees only the kept rows, so two of its notes describe them, not the file.
                if (rowsKept < rowsRead)
                    notes.Add($"{rowsRead - rowsKept} of {rowsRead} SDRF row(s) were set aside because the search does not " +
                        "read their files. Line numbers in the design report count the kept rows, not the lines of the SDRF " +
                        "file, and where it says every row of the SDRF is searched it means the kept rows, so a gap in " +
                        "replicate or fraction numbers may be a sample whose file is not searched.");
                scoped = scoping.Document;
            }

            if (options.Tag is not null)
            {
                var isobaric = SdrfIsobaricDesign.Read(scoped, new SdrfIsobaricDesignOptions
                {
                    Tag = options.Tag,
                    PlexColumn = options.PlexColumn,
                    PlexFileNamePattern = options.PlexFileNamePattern,
                    SinglePlex = options.SinglePlex,
                    ConditionColumns = conditionColumns,
                    SearchedFiles = searched,
                });
                return new SdrfSearchDesign(null, isobaric, scoping, new List<string>(), notes, joinReport);
            }

            var labelFree = SdrfLabelFreeDesign.Read(scoped, new SdrfLabelFreeDesignOptions
            {
                ConditionColumns = conditionColumns,
                SearchedFiles = searched,
            });
            return new SdrfSearchDesign(labelFree, null, scoping, new List<string>(), notes, joinReport);
        }

        /// <summary>
        /// The join's own lists, only those with entries, SharedStem first: it means one searched file whose rows
        /// name two different acquired files, which only the depositor can resolve.
        /// </summary>
        private static void ReportJoin(SdrfSearchScoping scoping, int rowsRead, int rowsKept, List<string> joinReport)
        {
            var lists = new (string Name, IReadOnlyList<string> Entries)[]
            {
                ("SharedStem (its rows name different acquired files and every one is kept, so the design is refused; " +
                 "if they list one run under two names, such as X.raw and X.mzML, keep one row, and otherwise only the " +
                 "depositor can say which file is right)", scoping.SharedStem),
                ("SearchedWithoutRow", scoping.SearchedWithoutRow),
                ("DroppedDataFiles (the search does not read them)", scoping.DroppedDataFiles),
                ("Ambiguous (claimed by two or more searched files; the first, in name order, keeps the rows, so the " +
                 "others have none and the design is refused)", scoping.Ambiguous),
            };
            if (rowsKept == rowsRead && lists.All(l => l.Entries.Count == 0))
                return;
            joinReport.Add($"SDRF joined to the searched files: {rowsKept} of {rowsRead} SDRF rows kept.");
            foreach (var (name, entries) in lists)
                foreach (string entry in entries)
                    joinReport.Add($"  {name}: {entry}");
        }

        /// <summary>A bare name is its factor value column; a <c>factor value[</c> prefix in any case is written as specified.</summary>
        private static string ConditionColumnName(string column)
        {
            string name = column.Trim();
            if (name.StartsWith(FactorValuePrefix, StringComparison.OrdinalIgnoreCase))
                return FactorValuePrefix + name.Substring(FactorValuePrefix.Length);
            return name.Contains('[') ? name : FactorValuePrefix + name + "]";
        }

        private static SdrfSearchDesign Refused(string reason) =>
            new(null, null, null, new List<string> { reason }, new List<string>(), new List<string>());

        private void ThrowIfRefused()
        {
            if (!IsValid)
                throw new InvalidOperationException("The design was refused and cannot be used:" + Environment.NewLine + Report());
        }
    }
}
