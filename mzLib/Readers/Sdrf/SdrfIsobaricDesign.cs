using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using MassSpectrometry;
using Omics.Modifications;

namespace Readers
{
    /// <summary>
    /// What a caller tells <see cref="SdrfIsobaricDesign.Read(SdrfDocument, SdrfIsobaricDesignOptions)"/>
    /// that the SDRF alone cannot say: the kit, and where each file's plex is written.
    /// </summary>
    public sealed class SdrfIsobaricDesignOptions
    {
        /// <summary>
        /// The kit the search used. Required. It supplies the channel list and reporter m/z order;
        /// labels are never used to guess it, because one kit's labels can be a subset of another's
        /// (a TMT11 plex that left 131C unannotated reads exactly like TMT10). MetaMorpheus passes the
        /// tag of <c>SearchParameters.MultiplexModId</c>.
        /// </summary>
        public IsobaricMassTag? Tag { get; init; }

        /// <summary>
        /// The column each row's plex is read from, by its exact name. One of <see cref="PlexColumn"/>,
        /// <see cref="PlexFileNamePattern"/> and <see cref="SinglePlex"/> must be given. No column is
        /// read as the plex unless named here, <c>comment[sample preparation batch]</c> included: in
        /// PXD061609 that column holds three SAX fractionation arms of one plex (MAP-09, N7).
        /// </summary>
        public string? PlexColumn { get; init; }

        /// <summary>
        /// A regular expression matched against each row's file name (the name the design is keyed
        /// on, MAP-12). The plex is its group named <c>plex</c>, else its first group. A file name it
        /// does not match is refused. Matching ignores case.
        /// </summary>
        public string? PlexFileNamePattern { get; init; }

        /// <summary>The name of the one plex every file belongs to.</summary>
        public string? SinglePlex { get; init; }

        /// <summary>
        /// The <c>factor value[...]</c> columns the condition is built from, exactly as
        /// <see cref="SdrfLabelFreeDesignOptions.ConditionColumns"/> (MAP-07).
        /// </summary>
        public IReadOnlyList<string>? ConditionColumns { get; init; }

        /// <summary>
        /// The files the search will read, exactly as <see cref="SdrfLabelFreeDesignOptions.SearchedFiles"/>:
        /// rows naming any other file are dropped and reported before any row is checked, and every
        /// searched file must be named exactly by some row.
        /// </summary>
        public IReadOnlyCollection<string>? SearchedFiles { get; init; }
    }

    /// <summary>One annotated channel of a plex.</summary>
    public sealed class SdrfIsobaricChannel
    {
        /// <summary>The bare channel label as the kit spells it (<c>127N</c>).</summary>
        public required string Label { get; init; }

        /// <summary>The SDRF <c>source name</c> (MAP-01).</summary>
        public required string SampleName { get; init; }

        public required string Condition { get; init; }

        /// <summary>
        /// 1-based and passed through as written: isobaric numbers are never renumbered (N3). Only a reference
        /// channel whose SDRF gives no number (<c>pooled</c>) is numbered, and the <see cref="SdrfIsobaricDesign.Notes"/>
        /// say so.
        /// </summary>
        public required int BiologicalReplicate { get; init; }

        /// <summary>
        /// True for a pooled reference channel (GR-20): kept and carried, and never a sample in a comparison. The SDRF
        /// marks one in any of four ways: <c>characteristics[sample type]</c> <c>reference</c> or <c>pooled</c>,
        /// <c>comment[reference channel]</c> naming the channel, <c>characteristics[pooled sample]</c> other than
        /// <c>not pooled</c>, or a biological replicate of <c>pooled</c>. A role is never inferred from a name.
        /// </summary>
        public bool IsReference { get; init; }
    }

    /// <summary>One file of an isobaric design and the channels its plex annotates.</summary>
    public sealed class SdrfIsobaricFile
    {
        /// <summary>The searched path when <see cref="SdrfIsobaricDesignOptions.SearchedFiles"/> was given, else the file name.</summary>
        public required string FilePath { get; init; }

        public required string Plex { get; init; }

        /// <summary>1-based, as written (MAP-34).</summary>
        public required int Fraction { get; init; }

        /// <summary>1-based, as written (MAP-21).</summary>
        public required int TechnicalReplicate { get; init; }

        /// <summary>The plex's annotated channels, in reporter m/z order. Unannotated channels are absent.</summary>
        public required IReadOnlyList<SdrfIsobaricChannel> Channels { get; init; }
    }

    /// <summary>
    /// An isobaric (TMT, iTRAQ, DiLeu) experimental design read from an SDRF (QuantProject M7): either a
    /// design MetaMorpheus will accept, or a refusal that lists every reason it would not.
    ///
    /// <para>Same discipline as <see cref="SdrfLabelFreeDesign"/>: every problem is reported at once and
    /// nothing is produced when anything fails, because MetaMorpheus quantifies nothing, with only a
    /// warning, from a design it rejects. The checks are MetaMorpheus's own
    /// (<c>TmtExperimentalDesign.Read</c>): one plex, fraction and technical replicate per file; a
    /// channel described once per plex; no (sample, biological replicate, fraction, technical
    /// replicate) twice within a plex.</para>
    ///
    /// <para><b>What is never guessed.</b> The plex comes only from what the caller declares (MAP-09);
    /// the kit only from <see cref="SdrfIsobaricDesignOptions.Tag"/>. A channel a drafted SDRF could
    /// not name (<c>source name</c> unknown, or <c>comment[biological replicate source]</c> =
    /// <c>default</c>) is refused by name rather than read (QP-S23). Biological replicates are passed
    /// through, not ranked or shifted (N3).</para>
    ///
    /// <para><b>Pooled reference channels</b> (GR-20, sdrf QP-S28 and QP-S31) are never refused. A channel the SDRF marks
    /// as one, in any of the four ways <see cref="SdrfIsobaricChannel.IsReference"/> lists, is kept and carried as a
    /// reference channel (<see cref="IsobaricQuantSampleInfo.IsReferenceChannel"/>, and <c>Sample Type</c> in
    /// <c>TmtDesign.txt</c>), so statistics leave it out of every comparison. Any one signal is enough; one that
    /// disagrees is reported, and the channel is still a reference. Its biological replicate, when the SDRF gives none,
    /// is numbered here. Every other sample type (carrier, empty, single cell, ...) is read as a study sample, as
    /// before.</para>
    ///
    /// <para><b>A plex is one channel-to-sample map.</b> As in MetaMorpheus, channels are collected per
    /// plex and every file of the plex carries all of them, so fractions share their samples.</para>
    /// </summary>
    public sealed class SdrfIsobaricDesign
    {
        public const string BiologicalReplicateSourceColumn = "comment[biological replicate source]";

        /// <summary>The value a drafter writes when it assigned a biological replicate nobody stated (QP-S23).</summary>
        public const string DefaultReplicateSource = "default";

        public const string SampleTypeColumn = "characteristics[sample type]";
        public const string PooledSampleColumn = "characteristics[pooled sample]";
        public const string ReferenceChannelColumn = "comment[reference channel]";

        /// <summary>
        /// The header written by <see cref="WriteTmtDesign"/>: MetaMorpheus's <c>TmtDesign.txt</c> without
        /// the optional <c>Sample Type</c> column, which MetaMorpheus reads as study sample when absent. The column is
        /// added, as <see cref="TmtDesignSampleTypeColumn"/>, only when a reference channel is marked.
        /// </summary>
        public const string TmtDesignHeader = "File\tPlex\tSample Name\tTMT Channel\tCondition\tBiological Replicate\tFraction\tTechnical Replicate";

        /// <summary>MetaMorpheus's optional <c>TmtDesign.txt</c> column, read since MetaMorpheus 1.1.10.</summary>
        public const string TmtDesignSampleTypeColumn = "Sample Type";

        // The sample-type words that mark a reference (sdrf 041), and the pooled-sample word that marks none.
        private static readonly HashSet<string> ReferenceSampleTypes = new(StringComparer.OrdinalIgnoreCase) { "reference", SdrfReserved.Pooled };
        private const string NotPooled = "not pooled";

        private readonly List<SdrfIsobaricFile> _files;

        private SdrfIsobaricDesign(List<SdrfIsobaricFile> files, List<string> refusals, List<string> notes,
            string? fileKeyColumn, List<string> conditionColumns, string plexSource, IsobaricMassTag? tag)
        {
            _files = files;
            Refusals = refusals;
            Notes = notes;
            FileKeyColumn = fileKeyColumn;
            ConditionColumns = conditionColumns;
            PlexSource = plexSource;
            Tag = tag;
        }

        /// <summary>True when there is nothing in <see cref="Refusals"/> and the design can be used.</summary>
        public bool IsValid => Refusals.Count == 0;

        /// <summary>One entry per file, in SDRF row order. Empty when the design is refused.</summary>
        public IReadOnlyList<SdrfIsobaricFile> Files => IsValid ? _files : Array.Empty<SdrfIsobaricFile>();

        /// <summary>Every reason the design was refused. Empty when it is valid.</summary>
        public IReadOnlyList<string> Refusals { get; }

        /// <summary>What was dropped or filled on the way to a valid design.</summary>
        public IReadOnlyList<string> Notes { get; }

        /// <summary>The column the files were keyed on (MAP-12), or null when neither was present.</summary>
        public string? FileKeyColumn { get; }

        /// <summary>The factor columns the condition was built from, in join order.</summary>
        public IReadOnlyList<string> ConditionColumns { get; }

        /// <summary>Where the plex was read from, as the report prints it.</summary>
        public string PlexSource { get; }

        /// <summary>The kit the design was read against.</summary>
        public IsobaricMassTag? Tag { get; }

        /// <summary>
        /// The plexes by the <see cref="IsobaricQuantSampleInfo.PlexId"/> each receives: 1..N in
        /// case-insensitive name order, the rule MetaMorpheus's <c>ToMzLibDesign</c> uses, so both paths
        /// number the same plexes the same way.
        /// </summary>
        public IReadOnlyDictionary<string, int> PlexIds => _files
            .Select(f => f.Plex)
            .Distinct(StringComparer.OrdinalIgnoreCase)
            .OrderBy(p => p, StringComparer.OrdinalIgnoreCase)
            .Select((plex, index) => (plex, id: index + 1))
            .ToDictionary(x => x.plex, x => x.id, StringComparer.OrdinalIgnoreCase);

        /// <summary>
        /// The design as mzLib quantification takes it: every channel of the kit for every file, in
        /// reporter m/z order, because quantification aligns samples with reporter intensities by
        /// position. A channel the plex does not annotate is empty (no condition, biological replicate
        /// 0, no sample name), as in MetaMorpheus. Throws when the design was refused.
        /// </summary>
        public SampleExperimentalDesign ToExperimentalDesign()
        {
            ThrowIfRefused();

            var plexIds = PlexIds;
            var samples = new List<ISampleInfo>();
            foreach (var file in _files)
            {
                var byLabel = file.Channels.ToDictionary(c => c.Label, StringComparer.OrdinalIgnoreCase);
                for (int i = 0; i < Tag!.ChannelLabels.Count; i++)
                {
                    string label = Tag.ChannelLabels[i];
                    byLabel.TryGetValue(label, out var channel);
                    samples.Add(new IsobaricQuantSampleInfo(
                        fullFilePathWithExtension: file.FilePath,
                        condition: channel?.Condition ?? string.Empty,
                        biologicalReplicate: channel?.BiologicalReplicate ?? 0,
                        technicalReplicate: file.TechnicalReplicate,
                        fraction: file.Fraction,
                        plexId: plexIds[file.Plex],
                        channelLabel: label,
                        reporterIonMz: Tag.ReporterIonMzs[i],
                        isReferenceChannel: channel?.IsReference ?? false)
                    {
                        SampleName = channel?.SampleName,
                    });
                }
            }

            return SampleExperimentalDesign.FromSamples(samples);
        }

        /// <summary>
        /// Writes MetaMorpheus's <c>TmtDesign.txt</c>: one row per file and annotated channel, 1-based,
        /// channels in reporter m/z order (MAP-29). When a reference channel is marked, MetaMorpheus's
        /// <c>Sample Type</c> column is added (<c>reference</c> or <c>study sample</c>); otherwise the file is written
        /// without it, exactly as before. Throws when the design was refused, so a refusal can never leave a file
        /// behind for MetaMorpheus to find.
        /// </summary>
        public void WriteTmtDesign(string path)
        {
            ThrowIfRefused();

            bool withSampleType = _files.Any(f => f.Channels.Any(c => c.IsReference));
            var text = new StringBuilder();
            text.Append(TmtDesignHeader);
            if (withSampleType)
                text.Append('\t').Append(TmtDesignSampleTypeColumn);
            text.Append('\n');
            foreach (var file in _files)
            {
                foreach (var channel in file.Channels)
                {
                    text.Append(file.FilePath).Append('\t')
                        .Append(file.Plex).Append('\t')
                        .Append(channel.SampleName).Append('\t')
                        .Append(channel.Label).Append('\t')
                        .Append(channel.Condition).Append('\t')
                        .Append(channel.BiologicalReplicate.ToString(CultureInfo.InvariantCulture)).Append('\t')
                        .Append(file.Fraction.ToString(CultureInfo.InvariantCulture)).Append('\t')
                        .Append(file.TechnicalReplicate.ToString(CultureInfo.InvariantCulture));
                    if (withSampleType)
                        text.Append('\t').Append(channel.IsReference ? "reference" : "study sample");
                    text.Append('\n');
                }
            }

            File.WriteAllText(path, text.ToString(), new UTF8Encoding(false));
        }

        /// <summary>A human-readable account of the projection: the columns used, notes, and every refusal.</summary>
        public string Report()
        {
            var report = new StringBuilder();
            report.AppendLine(IsValid
                ? $"Isobaric design read from SDRF: {_files.Count} file(s), {PlexIds.Count} plex(es), kit {Tag!.TagType}."
                : $"Isobaric design REFUSED: {Refusals.Count} reason(s). No design was produced.");
            report.AppendLine($"Files keyed on: {FileKeyColumn ?? "(none)"}");
            report.AppendLine($"Plex from: {PlexSource}");
            report.AppendLine($"Condition from: {(ConditionColumns.Count == 0 ? "(none)" : string.Join(" + ", ConditionColumns))}");
            if (IsValid)
            {
                var plexIds = PlexIds;
                var references = _files
                    .GroupBy(f => f.Plex, StringComparer.OrdinalIgnoreCase)
                    .OrderBy(p => plexIds[p.Key])
                    .Select(p => (Plex: p.Key, Labels: p.First().Channels.Where(c => c.IsReference).Select(c => c.Label).ToList()))
                    .Where(p => p.Labels.Count > 0)
                    .ToList();
                report.AppendLine(references.Count == 0
                    ? "Reference channels: none marked; every annotated channel is a study sample."
                    : $"Reference channels: {references.Sum(p => p.Labels.Count)} (" +
                      string.Join("; ", references.Select(p => $"plex '{p.Plex}': {string.Join(", ", p.Labels)}")) +
                      "), kept and carried, and left out of every comparison.");
            }
            foreach (var note in Notes)
                report.AppendLine("  note: " + note);
            foreach (var refusal in Refusals)
                report.AppendLine("  refused: " + refusal);
            return report.ToString();
        }

        /// <summary>Reads an SDRF file. See <see cref="Read(SdrfDocument, SdrfIsobaricDesignOptions)"/>.</summary>
        public static SdrfIsobaricDesign Read(string sdrfPath, SdrfIsobaricDesignOptions options)
        {
            if (sdrfPath == null)
                throw new ArgumentNullException(nameof(sdrfPath));
            return Read(new SdrfDocument(sdrfPath), options);
        }

        /// <summary>
        /// Projects an SDRF onto an isobaric design. Never throws for a bad document or bad options:
        /// every problem is a refusal, listed in <see cref="Refusals"/>.
        /// </summary>
        public static SdrfIsobaricDesign Read(SdrfDocument sdrf, SdrfIsobaricDesignOptions options)
        {
            if (sdrf == null)
                throw new ArgumentNullException(nameof(sdrf));
            if (options == null)
                throw new ArgumentNullException(nameof(options));

            var refusals = new List<string>();
            var notes = new List<string>();
            var header = sdrf.Header;
            var rows = sdrf.Results.ToList();
            var tag = options.Tag;

            if (tag == null)
                refusals.Add("No isobaric tag was given, so the channels are unknown. Pass the kit the search used; it is never guessed from the labels.");

            var plexReader = PlexReader.Create(options, header, refusals);

            string? keyColumn = SdrfDesignRules.ChooseFileKeyColumn(header, refusals);
            foreach (var column in new[] { SdrfLabelFreeDesign.LabelColumn, SdrfLabelFreeDesign.SourceNameColumn, SdrfLabelFreeDesign.BiologicalReplicateColumn })
            {
                if (!header.Contains(column))
                    refusals.Add($"The SDRF has no '{column}' column, which every channel needs.");
            }

            var conditionColumns = SdrfDesignRules.ResolveConditionColumns(header, options.ConditionColumns, refusals);
            bool conditionDeclared = options.ConditionColumns is { Count: > 0 };

            if (rows.Count == 0)
                refusals.Add("The SDRF has no rows.");

            if (refusals.Count > 0)
                return Refused(refusals, notes, keyColumn, conditionColumns, plexReader, tag);

            Dictionary<string, string>? searchedByName = null;
            if (options.SearchedFiles != null)
            {
                searchedByName = SdrfDesignRules.IndexSearchedFiles(options.SearchedFiles, refusals);
                if (refusals.Count > 0)
                    return Refused(refusals, notes, keyColumn, conditionColumns, plexReader, tag);
            }

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

                var row = ParseRow(rows[i], line, fileName, keyColumn!, tag!, plexReader!, conditionColumns, conditionDeclared, refusals);
                if (row != null)
                    parsed.Add(row);
                else if (fileName.Length > 0)
                    unreadFiles.Add(fileName);
            }

            NumberPooledReferences(parsed, tag!, notes);
            NoteDisagreeingReferenceSignals(parsed, notes);

            SdrfDesignRules.RefuseConditionCollisions(parsed.Select(r => (r.Condition, r.FactorValues)), refusals);
            RefuseInconsistentFiles(parsed, refusals);
            var channelsByPlex = CollectChannelsByPlex(parsed, tag!, refusals);
            RefuseRepeatedSamples(parsed, refusals);

            if (searchedByName != null)
                SdrfDesignRules.RefuseSearchedFilesWithoutRows(parsed.Select(r => r.FileName), unreadFiles, dropped, searchedByName, refusals);

            if (refusals.Count > 0)
                return Refused(refusals, notes, keyColumn, conditionColumns, plexReader, tag);

            var files = new List<SdrfIsobaricFile>();
            foreach (var file in parsed.GroupBy(r => r.FileName, StringComparer.Ordinal))
            {
                var first = file.First();
                var channels = channelsByPlex[first.Plex];
                var missing = channels.Select(c => c.Label)
                    .Except(file.Select(r => r.Label), StringComparer.OrdinalIgnoreCase)
                    .ToList();
                if (missing.Count > 0)
                {
                    notes.Add($"'{file.Key}' has no row for channel(s) {string.Join(", ", missing)}; " +
                              $"they take plex '{first.Plex}''s annotation, as every file of a plex does.");
                }

                files.Add(new SdrfIsobaricFile
                {
                    FilePath = searchedByName?[file.Key] ?? file.Key,
                    Plex = first.Plex,
                    Fraction = first.Fraction,
                    TechnicalReplicate = first.TechnicalReplicate,
                    Channels = channels,
                });
            }

            return new SdrfIsobaricDesign(files, refusals, notes, keyColumn, conditionColumns, plexReader!.Description, tag);
        }

        private static SdrfIsobaricDesign Refused(List<string> refusals, List<string> notes, string? keyColumn,
            List<string> conditionColumns, PlexReader? plexReader, IsobaricMassTag? tag) =>
            new(new List<SdrfIsobaricFile>(), refusals, notes, keyColumn, conditionColumns,
                plexReader?.Description ?? "(not declared)", tag);

        private sealed class ParsedRow
        {
            public required int Line { get; init; }
            public required string FileName { get; init; }
            public required string Plex { get; init; }
            public required string Label { get; init; }
            public required string SampleName { get; init; }
            public required string Condition { get; init; }
            public required IReadOnlyList<string> FactorValues { get; init; }
            /// <summary>The SDRF's number, or null for a reference channel that gives none (<c>pooled</c>).</summary>
            public required int? GivenBiologicalReplicate { get; init; }
            /// <summary>The biological replicate cell as written, for the note that reports an added number.</summary>
            public required string BiologicalReplicateText { get; init; }
            /// <summary>The given number, or the one <see cref="NumberPooledReferences"/> adds.</summary>
            public int BiologicalReplicate { get; set; }
            public required int Fraction { get; init; }
            public required int TechnicalReplicate { get; init; }
            /// <summary>The reference signals this row carries (GR-20); any one makes the channel a reference.</summary>
            public required IReadOnlyList<string> ReferenceSignals { get; init; }
            /// <summary>What this row says that reads as a study sample, reported when a reference signal is also present.</summary>
            public required IReadOnlyList<string> StudySampleSignals { get; init; }
            public bool IsReference => ReferenceSignals.Count > 0;
        }

        private static ParsedRow? ParseRow(SdrfRow row, int line, string fileName, string keyColumn, IsobaricMassTag tag,
            PlexReader plexReader, List<string> conditionColumns, bool conditionDeclared, List<string> refusals)
        {
            int before = refusals.Count;
            string where = $"Line {line}{SdrfDesignRules.Describe(fileName)}";

            if (fileName.Length == 0)
                refusals.Add($"Line {line}: '{keyColumn}' is empty.");

            string? label = ReadChannel(row, where, tag, refusals);

            // QP-S23: a drafted SDRF writes every channel of the plex, and marks the ones whose sample
            // nobody stated. Those are refused by name, never read as samples.
            string sampleName = (row[SdrfLabelFreeDesign.SourceNameColumn] ?? string.Empty).Trim();
            if (sampleName.Length == 0)
                refusals.Add($"{where}: '{SdrfLabelFreeDesign.SourceNameColumn}' is empty.");
            else if (SdrfDesignRules.IsUnknownWord(sampleName))
                refusals.Add($"{where}: channel {label ?? "?"} has '{SdrfLabelFreeDesign.SourceNameColumn}' '{sampleName}'. " +
                             "A channel whose sample is unknown cannot be quantified as a sample; name it, or remove the row.");

            string? replicateSource = row[BiologicalReplicateSourceColumn]?.Trim();
            if (string.Equals(replicateSource, DefaultReplicateSource, StringComparison.OrdinalIgnoreCase))
                refusals.Add($"{where}: channel {label ?? "?"} ('{sampleName}') has '{BiologicalReplicateSourceColumn}' '{replicateSource}', " +
                             "so its biological replicate was filled in rather than stated. State it, or remove the row.");

            var factorValues = SdrfDesignRules.ReadFactorValues(row, line, fileName, conditionColumns, conditionDeclared, refusals);
            string? plex = plexReader.Read(row, fileName, where, refusals);

            // GR-20: a reference channel is never refused for giving no number (`pooled`, or an unknown word);
            // one is added later, per plex and sample. A study sample's number is required as before.
            string biorepText = (row[SdrfLabelFreeDesign.BiologicalReplicateColumn] ?? string.Empty).Trim();
            var (referenceSignals, studySampleSignals) = ReadReferenceSignals(row, label, tag, biorepText);
            bool givesNoNumber = string.Equals(biorepText, SdrfReserved.Pooled, StringComparison.OrdinalIgnoreCase)
                                 || biorepText.Length == 0 || SdrfDesignRules.IsUnknownWord(biorepText);
            int? biorep = referenceSignals.Count > 0 && givesNoNumber
                ? null
                : SdrfDesignRules.ReadPositiveInteger(row, SdrfLabelFreeDesign.BiologicalReplicateColumn, line, fileName, required: true, refusals);
            int fraction = SdrfDesignRules.ReadPositiveInteger(row, SdrfLabelFreeDesign.FractionColumn, line, fileName, required: false, refusals);
            int techrep = SdrfDesignRules.ReadPositiveInteger(row, SdrfLabelFreeDesign.TechnicalReplicateColumn, line, fileName, required: false, refusals);

            if (refusals.Count > before)
                return null;

            return new ParsedRow
            {
                Line = line,
                FileName = fileName,
                Plex = plex!,
                Label = label!,
                SampleName = sampleName,
                Condition = SdrfDesignRules.JoinCondition(factorValues),
                FactorValues = factorValues,
                GivenBiologicalReplicate = biorep,
                BiologicalReplicateText = biorepText,
                BiologicalReplicate = biorep ?? 0,
                Fraction = fraction,
                TechnicalReplicate = techrep,
                ReferenceSignals = referenceSignals,
                StudySampleSignals = studySampleSignals,
            };
        }

        /// <summary>
        /// The four reference signals (sdrf 041, ruled 2026-10-08), each read on its own; any one makes the channel a
        /// reference. Beside them, what the row says that reads as a study sample (a sample type that names something
        /// else, a pooled sample of <c>not pooled</c>), so a disagreement can be reported. Empty cells and
        /// <c>not available</c> / <c>not applicable</c> say neither.
        /// </summary>
        private static (List<string> Reference, List<string> StudySample) ReadReferenceSignals(SdrfRow row, string? channel,
            IsobaricMassTag tag, string biorepText)
        {
            var reference = new List<string>();
            var studySample = new List<string>();

            string sampleTypeCell = (row[SampleTypeColumn] ?? string.Empty).Trim();
            string sampleType = sampleTypeCell.Length == 0 ? string.Empty : SdrfQuantAuditor.ReadLabel(sampleTypeCell).Name;
            if (ReferenceSampleTypes.Contains(sampleType))
                reference.Add($"sample type '{sampleType}'");
            else if (sampleType.Length > 0 && !SdrfDesignRules.IsUnknownWord(sampleType))
                studySample.Add($"sample type '{sampleType}'");

            // The deposits write the plex's reference channel on every row (TMT127N), so it marks only the row it names.
            string referenceChannel = (row[ReferenceChannelColumn] ?? string.Empty).Trim();
            if (channel != null && referenceChannel.Length > 0 && !SdrfDesignRules.IsUnknownWord(referenceChannel)
                && NamesChannel(referenceChannel, channel, tag))
            {
                reference.Add($"reference channel '{referenceChannel}'");
            }

            string pooledSample = (row[PooledSampleColumn] ?? string.Empty).Trim();
            if (string.Equals(pooledSample, NotPooled, StringComparison.OrdinalIgnoreCase))
                studySample.Add($"pooled sample '{pooledSample}'");
            else if (pooledSample.Length > 0 && !SdrfDesignRules.IsUnknownWord(pooledSample))
                reference.Add("pooled sample");

            if (string.Equals(biorepText, SdrfReserved.Pooled, StringComparison.OrdinalIgnoreCase))
                reference.Add($"biological replicate '{biorepText}'");

            return (reference, studySample);
        }

        // A reference-channel cell names one channel, or several separated by commas, in any of the label forms the
        // label column takes (TMT127N, 127N, PRIDE's TMT10 TMT131).
        private static bool NamesChannel(string cell, string channel, IsobaricMassTag tag) =>
            new[] { cell }.Concat(cell.Split(',')).Any(item =>
                SdrfIsobaricLabels.TryGetChannel(tag.TagType, item.Trim(), out string? named)
                && string.Equals(named, channel, StringComparison.OrdinalIgnoreCase));

        /// <summary>
        /// Numbers each reference channel whose SDRF gives no biological replicate. MetaMorpheus's TmtDesign.txt needs a
        /// number of 1 or more, and two channels of one sample in a plex need different ones. So, per plex and sample, in
        /// reporter m/z order, after the highest number that sample is given in its plex: an added number never repeats
        /// one. Every number added is reported.
        /// </summary>
        private static void NumberPooledReferences(List<ParsedRow> rows, IsobaricMassTag tag, List<string> notes)
        {
            foreach (var sample in rows.GroupBy(r => (Plex: r.Plex.ToUpperInvariant(), Sample: r.SampleName.ToUpperInvariant())))
            {
                var unnumbered = sample.Where(r => r.GivenBiologicalReplicate == null).ToList();
                if (unnumbered.Count == 0)
                    continue;

                int highest = sample.Where(r => r.GivenBiologicalReplicate != null)
                    .Select(r => r.GivenBiologicalReplicate!.Value).DefaultIfEmpty(0).Max();
                var labels = tag.ChannelLabels
                    .Where(l => unnumbered.Any(r => string.Equals(r.Label, l, StringComparison.OrdinalIgnoreCase)))
                    .ToList();
                var added = new List<string>();
                foreach (var label in labels)
                {
                    int number = ++highest;
                    foreach (var row in unnumbered.Where(r => string.Equals(r.Label, label, StringComparison.OrdinalIgnoreCase)))
                        row.BiologicalReplicate = number;
                    added.Add($"{label} -> {number}");
                }

                var first = sample.First();
                string written = string.Join("', '", unnumbered.Select(r => r.BiologicalReplicateText).Distinct(StringComparer.OrdinalIgnoreCase));
                notes.Add($"Plex '{first.Plex}' sample '{first.SampleName}': channel(s) {string.Join(", ", labels)} give '{written}' " +
                          "as their biological replicate, a pooled reference with no number, so these were numbered here: " +
                          string.Join(", ", added) + ".");
            }
        }

        // sdrf 041: two signals that disagree are reported, and the channel is read as a reference, because a pool read as
        // a study sample pulls its group towards the mean of all groups. Once per plex and channel, not per fraction.
        private static void NoteDisagreeingReferenceSignals(List<ParsedRow> rows, List<string> notes)
        {
            foreach (var channel in rows.Where(r => r.IsReference && r.StudySampleSignals.Count > 0)
                         .GroupBy(r => (Plex: r.Plex.ToUpperInvariant(), Label: r.Label.ToUpperInvariant())))
            {
                var first = channel.First();
                notes.Add($"Plex '{first.Plex}' channel {first.Label} ('{first.SampleName}') is marked a reference by " +
                          $"{string.Join(", ", channel.SelectMany(r => r.ReferenceSignals).Distinct())}, but also says " +
                          $"{string.Join(", ", channel.SelectMany(r => r.StudySampleSignals).Distinct())}. It is read as a reference: " +
                          "a pool read as a study sample would pull its group towards the mean of all groups.");
            }
        }

        // MAP-06 / N4: read per cell, in either form, case-insensitively, and return the kit's own
        // spelling of the channel. A family prefix that names a different reagent is refused, because
        // bare channel numbers overlap between kits (iTRAQ4 114-117, DiLeu4 115-118). The label rules,
        // TMT10's TMT131 included, are SdrfIsobaricLabels', which writers of SDRFs share (QP-S26).
        private static string? ReadChannel(SdrfRow row, string where, IsobaricMassTag tag, List<string> refusals)
        {
            string cell = (row[SdrfLabelFreeDesign.LabelColumn] ?? string.Empty).Trim();
            if (cell.Length == 0)
            {
                refusals.Add($"{where}: '{SdrfLabelFreeDesign.LabelColumn}' is empty.");
                return null;
            }

            string name = SdrfQuantAuditor.ReadLabel(cell).Name;
            if (name.IndexOf("label free", StringComparison.OrdinalIgnoreCase) >= 0)
            {
                refusals.Add($"{where}: '{SdrfLabelFreeDesign.LabelColumn}' is '{cell}', a label-free row in an isobaric design. " +
                             $"A label-free design is read with {nameof(SdrfLabelFreeDesign)}.");
                return null;
            }

            if (!SdrfIsobaricLabels.TryGetChannel(tag.TagType, cell, out string? channel))
            {
                refusals.Add($"{where}: '{SdrfLabelFreeDesign.LabelColumn}' is '{cell}', which is not a channel of {tag.TagType}. " +
                             $"Its channels are {string.Join(", ", tag.ChannelLabels)}.");
                return null;
            }

            return channel;
        }

        // MetaMorpheus: each file has one plex, one fraction and one technical replicate, and names each
        // channel at most once.
        private static void RefuseInconsistentFiles(List<ParsedRow> rows, List<string> refusals)
        {
            foreach (var file in rows.GroupBy(r => r.FileName, StringComparer.Ordinal))
            {
                var states = file.Select(r => (r.Plex, r.Fraction, r.TechnicalReplicate))
                    .Distinct(new PlexStateComparer())
                    .ToList();
                if (states.Count > 1)
                {
                    refusals.Add($"'{file.Key}' is given {states.Count} different plex/fraction/technical replicate combinations: " +
                                 string.Join("; ", states.Select(s => $"plex '{s.Plex}' fraction {s.Fraction} techrep {s.TechnicalReplicate}")) +
                                 ". A file is one plex, one fraction and one technical replicate.");
                }

                foreach (var channel in file.GroupBy(r => r.Label, StringComparer.OrdinalIgnoreCase).Where(g => g.Count() > 1))
                {
                    refusals.Add($"'{file.Key}' names channel {channel.Key} on {channel.Count()} rows (lines " +
                                 $"{string.Join(", ", channel.Select(r => r.Line))}). A channel holds one sample per file.");
                }
            }
        }

        // MetaMorpheus keys channel annotations per plex, so every row of a plex must describe a channel
        // the same way. A plex that disagrees with itself is refused, never resolved by row order.
        private static Dictionary<string, List<SdrfIsobaricChannel>> CollectChannelsByPlex(List<ParsedRow> rows,
            IsobaricMassTag tag, List<string> refusals)
        {
            var result = new Dictionary<string, List<SdrfIsobaricChannel>>(StringComparer.OrdinalIgnoreCase);
            foreach (var plex in rows.GroupBy(r => r.Plex, StringComparer.OrdinalIgnoreCase))
            {
                var channels = new List<SdrfIsobaricChannel>();
                foreach (var label in tag.ChannelLabels)
                {
                    var described = plex.Where(r => string.Equals(r.Label, label, StringComparison.OrdinalIgnoreCase)).ToList();
                    if (described.Count == 0)
                        continue;

                    // The role is part of the description, as MetaMorpheus's Sample Type is.
                    var descriptions = described
                        .GroupBy(r => (r.SampleName.ToUpperInvariant(), r.Condition.ToUpperInvariant(), r.BiologicalReplicate, r.IsReference))
                        .ToList();
                    if (descriptions.Count > 1)
                    {
                        refusals.Add($"Plex '{plex.Key}' channel {label} is described {descriptions.Count} ways: " +
                                     string.Join("; ", descriptions.Select(d =>
                                     {
                                         var r = d.First();
                                         return $"'{r.SampleName}' / '{r.Condition}' / biorep {r.BiologicalReplicate} / " +
                                                $"{(r.IsReference ? "reference" : "study sample")} (line(s) {string.Join(", ", d.Select(x => x.Line))})";
                                     })) +
                                     ". Every file of a plex shares one channel-to-sample map.");
                        continue;
                    }

                    var first = described[0];
                    channels.Add(new SdrfIsobaricChannel
                    {
                        Label = label,
                        SampleName = first.SampleName,
                        Condition = first.Condition,
                        BiologicalReplicate = first.BiologicalReplicate,
                        IsReference = first.IsReference,
                    });
                }

                result[plex.Key] = channels;
            }

            return result;
        }

        // MetaMorpheus's ValidateUniqueSampleBioFracTech: within a plex, no (sample, biological replicate,
        // fraction, technical replicate) twice. It is what refuses one plex split across fractionation arms
        // that reuse fraction numbers (PXD061609 read as a single plex).
        private static void RefuseRepeatedSamples(List<ParsedRow> rows, List<string> refusals)
        {
            var repeated = rows
                .GroupBy(r => (Plex: r.Plex.ToUpperInvariant(), Sample: r.SampleName.ToUpperInvariant(),
                    r.BiologicalReplicate, r.Fraction, r.TechnicalReplicate))
                .Where(g => g.Select(r => r.FileName).Distinct(StringComparer.Ordinal).Count() > 1
                            || g.Select(r => r.Label).Distinct(StringComparer.OrdinalIgnoreCase).Count() > 1);
            foreach (var group in repeated)
            {
                var first = group.First();
                refusals.Add($"Plex '{first.Plex}' has sample '{first.SampleName}' biorep {first.BiologicalReplicate} " +
                             $"fraction {first.Fraction} techrep {first.TechnicalReplicate} more than once: " +
                             string.Join(", ", group.Select(r => $"'{r.FileName}' {r.Label}").Distinct()) +
                             ". Files that share a fraction number need different technical replicates, or different plexes.");
            }
        }

        private void ThrowIfRefused()
        {
            if (!IsValid)
                throw new InvalidOperationException("The design was refused and cannot be used:" + Environment.NewLine + Report());
        }

        private sealed class PlexStateComparer : IEqualityComparer<(string Plex, int Fraction, int TechnicalReplicate)>
        {
            public bool Equals((string Plex, int Fraction, int TechnicalReplicate) x, (string Plex, int Fraction, int TechnicalReplicate) y) =>
                string.Equals(x.Plex, y.Plex, StringComparison.OrdinalIgnoreCase)
                && x.Fraction == y.Fraction && x.TechnicalReplicate == y.TechnicalReplicate;

            public int GetHashCode((string Plex, int Fraction, int TechnicalReplicate) state) =>
                HashCode.Combine(StringComparer.OrdinalIgnoreCase.GetHashCode(state.Plex), state.Fraction, state.TechnicalReplicate);
        }

        /// <summary>The one plex source the caller declared (MAP-09). Nothing else is ever read as the plex.</summary>
        private sealed class PlexReader
        {
            private readonly string? _column;
            private readonly Regex? _pattern;
            private readonly string? _single;

            private PlexReader(string? column, Regex? pattern, string? single, string description)
            {
                _column = column;
                _pattern = pattern;
                _single = single;
                Description = description;
            }

            public string Description { get; }

            public static PlexReader? Create(SdrfIsobaricDesignOptions options, SdrfHeader header, List<string> refusals)
            {
                int declared = new[] { options.PlexColumn, options.PlexFileNamePattern, options.SinglePlex }
                    .Count(v => !string.IsNullOrWhiteSpace(v));
                if (declared != 1)
                {
                    refusals.Add(declared == 0
                        ? "No plex source was declared. Name the column that holds the plex, give a file-name pattern " +
                          "whose group captures it, or name the single plex. The plex is never inferred (MAP-09): " +
                          "'comment[sample preparation batch]' is often not a plex."
                        : $"{declared} plex sources were declared; give exactly one of PlexColumn, PlexFileNamePattern and SinglePlex.");
                    return null;
                }

                if (!string.IsNullOrWhiteSpace(options.SinglePlex))
                    return new PlexReader(null, null, options.SinglePlex.Trim(), $"one plex, '{options.SinglePlex.Trim()}' (declared)");

                if (!string.IsNullOrWhiteSpace(options.PlexColumn))
                {
                    if (!header.Contains(options.PlexColumn))
                    {
                        refusals.Add($"The declared plex column '{options.PlexColumn}' is not in the SDRF.");
                        return null;
                    }
                    return new PlexReader(options.PlexColumn, null, null, $"column '{options.PlexColumn}' (declared)");
                }

                Regex pattern;
                try
                {
                    pattern = new Regex(options.PlexFileNamePattern!, RegexOptions.IgnoreCase | RegexOptions.CultureInvariant);
                }
                catch (ArgumentException e)
                {
                    refusals.Add($"The plex file-name pattern '{options.PlexFileNamePattern}' is not a valid regular expression: {e.Message}");
                    return null;
                }
                if (pattern.GetGroupNumbers().Length < 2)
                {
                    refusals.Add($"The plex file-name pattern '{options.PlexFileNamePattern}' has no group, so it captures no plex. " +
                                 "Put the plex in a group, e.g. 'TMT_?pool(?<plex>\\d+)'.");
                    return null;
                }
                return new PlexReader(null, pattern, null, $"file-name pattern '{options.PlexFileNamePattern}' (declared)");
            }

            public string? Read(SdrfRow row, string fileName, string where, List<string> refusals)
            {
                if (_single != null)
                    return _single;

                if (_column != null)
                {
                    string value = (row[_column] ?? string.Empty).Trim();
                    if (value.Length == 0 || SdrfDesignRules.IsUnknownWord(value))
                    {
                        refusals.Add($"{where}: the plex column '{_column}' is '{value}', which names no plex.");
                        return null;
                    }
                    return value;
                }

                var match = _pattern!.Match(fileName);
                Group group = match.Groups["plex"].Success ? match.Groups["plex"] : match.Groups[1];
                if (!match.Success || !group.Success || group.Value.Length == 0)
                {
                    refusals.Add($"{where}: the plex file-name pattern '{_pattern}' does not capture a plex from '{fileName}'.");
                    return null;
                }
                return group.Value;
            }
        }
    }
}
