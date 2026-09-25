using System.Text.RegularExpressions;
using MzLibUtil;
using UsefulProteomicsDatabases;

namespace Readers
{
    /// <summary>Where a drafted cell's value came from -- the provenance a drafted SDRF records (D31).</summary>
    internal enum SdrfDraftSource
    {
        /// <summary>Read off the file names and the record's text; the cell's evidence says how.</summary>
        Inferred,

        /// <summary>Taken from the PRIDE project record, where PRIDE lists exactly one value (D27).</summary>
        PrideProjectRecord,

        /// <summary>Nothing states it; written as <c>not available</c>.</summary>
        NotAvailable
    }

    /// <summary>
    /// One drafted cell: its value, where the value came from, and why. <see cref="Term"/> is the
    /// controlled-vocabulary term when the value is one (organism, organism part, disease, instrument).
    /// </summary>
    internal sealed record SdrfDraftCell(string Value, SdrfDraftSource Source, string Evidence, CvParam? Term = null)
    {
        internal static SdrfDraftCell NotAvailable(string why) => new(SdrfReserved.NotAvailable, SdrfDraftSource.NotAvailable, why);
    }

    /// <summary>One drafted row: one raw file of the deposit.</summary>
    internal sealed record SdrfDraftRow(
        string DataFile,
        SdrfDraftCell SourceName,
        SdrfDraftCell Organism,
        SdrfDraftCell OrganismPart,
        SdrfDraftCell Disease,
        SdrfDraftCell Instrument,
        SdrfDraftCell BiologicalReplicate,
        SdrfDraftCell TechnicalReplicate,
        SdrfDraftCell Fraction,
        IReadOnlyList<SdrfDraftCell> Factors);

    /// <summary>
    /// A drafted SDRF: one row per raw file, every cell with its provenance and evidence, and the factor
    /// columns the rows' <see cref="SdrfDraftRow.Factors"/> fill, in order.
    /// </summary>
    internal sealed record SdrfDraft(IReadOnlyList<SdrfDraftRow> Rows, IReadOnlyList<string> FactorColumns);

    /// <summary>
    /// Drafts an SDRF for a PRIDE deposit from what PRIDE gives: the project record and the raw-file
    /// names. It combines <see cref="SdrfFileNamePattern"/> (structure from how names vary),
    /// <see cref="SdrfRecordAnchors"/> (conditions the record's own text names) and
    /// <see cref="SdrfReplicateResolver"/> (replicates and what they count).
    ///
    /// <para><b>Every rule here was measured before it was written.</b> The sdrf project graded drafts blind
    /// against curated SDRFs for 240 deposits, two graders each, in rounds that each fixed what the last
    /// found (results/benchmark, 2026-09-23). The decisions it cites: D27 (a project fact only when PRIDE
    /// lists one value), D34 (disease), D35 (a deposited SDRF comes first -- this type is for when there is
    /// none, or to fill one's gaps), D36 (the SDRF governs quantification).</para>
    ///
    /// <para><b>What it will not do.</b> It never invents a value nothing states: a cell with no source is
    /// <c>not available</c>, and every inferred cell carries the evidence it was read from, so a caller can
    /// write it with provenance <c>inferred</c> and a reviewer can check it. With no structure in the names,
    /// every file is its own sample and biological replicate 1 -- never a count across the deposit.</para>
    ///
    /// <para>Pure: an already-fetched project and file list in, a draft out. Internal (D19) until a
    /// cross-assembly caller is named.</para>
    /// </summary>
    internal static class SdrfDrafter
    {
        /// <summary>A factor level naming the non-diseased arm of a case/control study (D34). Deliberately not
        /// "control", "DMSO" or "vehicle": untreated cancer cells are still cancer.</summary>
        private static readonly Regex CaseControlArm = new(
            @"^(neg|negative|sham|healthy|hc|normal|nondiseased|non-diseased)$", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        private static readonly Regex Parenthetical = new(@"\s*\(.*\)\s*$", RegexOptions.Compiled);

        /// <summary>
        /// Drafts one row per raw file. Throws only on a null argument; an empty file list drafts no rows.
        /// </summary>
        public static SdrfDraft Draft(PrideProject project, IEnumerable<string> rawFileNames)
        {
            if (project == null) throw new ArgumentNullException(nameof(project));
            if (rawFileNames == null) throw new ArgumentNullException(nameof(rawFileNames));
            var names = rawFileNames.Distinct(StringComparer.OrdinalIgnoreCase).OrderBy(n => n, StringComparer.Ordinal).ToList();
            if (names.Count == 0) return new SdrfDraft(Array.Empty<SdrfDraftRow>(), Array.Empty<string>());

            var text = new[] { project.Title, project.ProjectDescription, project.SampleProcessingProtocol, project.DataProcessingProtocol }
                .Concat(project.Keywords).Concat(project.ProjectTags).ToList();
            var structure = SdrfFileNamePattern.Read(names);
            var anchors = SdrfRecordAnchors.Read(names, text);
            var replicates = SdrfReplicateResolver.Read(names, text);
            var byName = structure.Files.ToDictionary(f => f.FileName, StringComparer.Ordinal);
            var markerOf = replicates.Files.ToDictionary(r => r.FileName, StringComparer.Ordinal);

            // ---- conditions: anchored in the record first, else what the names alone showed ----
            bool anchored = anchors.Factors.Count > 0;
            IReadOnlyList<string> LevelsOf(string file) => anchored
                ? anchors.LevelsByFile[file].Select(l => l.Length == 0 ? SdrfReserved.NotAvailable : l).ToList()
                : byName[file].FactorLevels;
            var factorEvidence = anchored
                ? anchors.Factors.Select(f => f.Evidence).ToList()
                : structure.Slots.Where(s => s.Role == SdrfFileNameRole.Factor).Select(s => s.Evidence).ToList();
            int factorCount = factorEvidence.Count;
            var factorColumns = Enumerable.Range(0, factorCount)
                .Select(i => i == 0 ? "factor value[condition]" : $"factor value[condition {i + 1}]").ToList();
            bool conditionKnown = factorCount > 0;
            string ConditionOf(string file) => byName[file].Family + "\u001f" + string.Join("\u001f", LevelsOf(file));

            // ---- samples and replicates (measured order: marker, named count, one stated sample, rank, 1) ----
            // A count the names do not start at 1 in every condition was ranked, and its evidence says so (MAP-33).
            var rankedFamilies = structure.Slots
                .Where(s => s.Renumbered && s.Role is SdrfFileNameRole.BiologicalReplicate or SdrfFileNameRole.Replicate)
                .Select(s => s.Family).ToHashSet();
            var reading = names.ToDictionary(n => n, n => ReadReplicates(n, byName[n], markerOf[n], replicates, structure.Found,
                rankedFamilies.Contains(byName[n].Family)), StringComparer.Ordinal);

            // A re-injection the names mark (NEG1rep beside NEG1) is its twin's sample: same key, same count.
            if (structure.Found)
                foreach (var twins in names.GroupBy(n => byName[n].SampleKey).Where(g => g.Any(n => byName[n].TechnicalReplicate > 1)))
                {
                    var first = twins.OrderBy(n => byName[n].TechnicalReplicate ?? 1).ThenBy(n => n, StringComparer.Ordinal).First();
                    foreach (var n in twins)
                        reading[n] = reading[n] with { Key = reading[first].Key, Bio = reading[first].Bio ?? reading[n].Bio };
                }

            var samplesInCondition = names.GroupBy(ConditionOf).ToDictionary(g => g.Key, g => g.Select(n => reading[n].Key).Distinct().Count());
            var rank = names.GroupBy(ConditionOf)
                .SelectMany(g => g.Select(n => reading[n].Key).Distinct().Select((key, i) => (g.Key, key, i + 1)))
                .ToDictionary(x => (x.Item1, x.Item2), x => x.Item3);
            var sourceName = names.GroupBy(n => reading[n].Key)
                .ToDictionary(g => g.Key, g => SdrfFileNamePattern.Stem(g.OrderBy(n => n, StringComparer.Ordinal).First()));

            // ---- project facts (D27), disease (D34) ----
            var organism = One(project.Organisms, "organism", Organism);
            var part = One(project.OrganismParts, "organism part", t => t);
            var instrument = One(project.Instruments, "instrument", t => t);
            var disease = One(project.Diseases, "disease", t => t);
            if (disease.Term != null && disease.Value.Equals("disease free", StringComparison.OrdinalIgnoreCase))
                disease = disease with { Value = "normal", Evidence = "PRIDE's project record lists 'Disease free'" };
            bool split = names.Any(n => LevelsOf(n).Any(CaseControlArm.IsMatch)) && names.Any(n => !LevelsOf(n).Any(CaseControlArm.IsMatch));

            var rows = new List<SdrfDraftRow>(names.Count);
            foreach (var n in names)
            {
                var r = reading[n];
                string cond = ConditionOf(n);
                // A file none of whose levels is known has no condition to rank within: ranking it with
                // the other unknowns would count across the deposit, which a draft never does.
                bool hasCondition = conditionKnown && LevelsOf(n).Any(l => l != SdrfReserved.NotAvailable);
                (int bio, string bioWhy) = r.Bio is int b ? (b, r.BioWhy)
                    : replicates.SingleBiologicalSample ? (1, replicates.SingleSampleEvidence)
                    : hasCondition && samplesInCondition[cond] <= SdrfReplicateResolver.MaxReplicates
                        ? (rank[(cond, r.Key)], "the sample's rank within its condition")
                        : (1, "no replicate structure found, so no count is claimed");
                var levels = LevelsOf(n);
                var factors = levels.Select((l, i) => l == SdrfReserved.NotAvailable
                    ? SdrfDraftCell.NotAvailable("this file carries none of the condition's levels")
                    : new SdrfDraftCell(l, SdrfDraftSource.Inferred, factorEvidence[i])).ToList();
                var rowDisease = split && disease.Term != null
                    ? levels.Any(CaseControlArm.IsMatch)
                        ? new SdrfDraftCell("normal", SdrfDraftSource.Inferred, "the control arm of a case/control split in the file names (D34)")
                        : disease with { Source = SdrfDraftSource.Inferred, Evidence = "a case arm of a case/control split; the project's disease (D34)" }
                    : disease;
                rows.Add(new SdrfDraftRow(
                    n,
                    new SdrfDraftCell(sourceName[r.Key], SdrfDraftSource.Inferred, r.KeyWhy),
                    organism, part, rowDisease, instrument,
                    new SdrfDraftCell(bio.ToString(System.Globalization.CultureInfo.InvariantCulture), SdrfDraftSource.Inferred, bioWhy),
                    new SdrfDraftCell((r.Tech ?? 1).ToString(System.Globalization.CultureInfo.InvariantCulture), SdrfDraftSource.Inferred, r.TechWhy),
                    new SdrfDraftCell((r.Frac ?? 1).ToString(System.Globalization.CultureInfo.InvariantCulture), SdrfDraftSource.Inferred, r.FracWhy),
                    factors));
            }
            return new SdrfDraft(rows, factorColumns);
        }

        /// <summary>
        /// Writes a draft as an SDRF through <see cref="SdrfBuilder"/>, with its provenance.
        ///
        /// <para><b>Provenance, option B (D31).</b> <c>comment[characteristics source]</c> is each row's
        /// default: the source most of its stated characteristics share. <c>comment[&lt;characteristic&gt;
        /// source]</c> is written only where one cell's source differs -- the control arm's inferred
        /// <c>normal</c> beside a project-record organism part. A cell nothing states is
        /// <c>not available</c> and has no source.</para>
        ///
        /// <para>Assay facts a draft cannot know (cleavage agent, modifications, tolerances) are
        /// <c>not available</c>: the search that uses this SDRF knows them and writes them in the SDRF it
        /// emits (D37). <c>assay name</c> is <c>run &lt;file stem&gt;</c> (MAP-13).</para>
        /// </summary>
        public static SdrfDocument ToDocument(SdrfDraft draft, string? proteomeXchangeAccession = null)
        {
            if (draft == null) throw new ArgumentNullException(nameof(draft));
            if (draft.Rows.Count == 0) throw new ArgumentException("A draft with no rows has no SDRF.", nameof(draft));
            ControlledVocabulary.Pride.TryGetByAccession("MS:1002038", out var labelFree);

            var rows = draft.Rows.Select(r =>
            {
                var stated = new List<(string Name, SdrfDraftCell Cell)>
                    { ("organism", r.Organism), ("organism part", r.OrganismPart), ("disease", r.Disease) }
                    .Where(x => x.Cell.Source != SdrfDraftSource.NotAvailable).ToList();
                var comments = new Dictionary<string, string>(StringComparer.Ordinal);
                if (stated.Count > 0)
                {
                    var byRow = stated.GroupBy(x => x.Cell.Source)
                        .OrderByDescending(g => g.Count()).ThenBy(g => g.Key == SdrfDraftSource.PrideProjectRecord ? 0 : 1)
                        .First().Key;
                    comments["comment[characteristics source]"] = SourceWord(byRow);
                    foreach (var (name, cell) in stated.Where(x => x.Cell.Source != byRow))
                        comments[$"comment[{name} source]"] = SourceWord(cell.Source);
                }

                var characteristics = new Dictionary<string, CvParam>(StringComparer.Ordinal);
                if (r.OrganismPart.Term != null) characteristics["characteristics[organism part]"] = r.OrganismPart.Term;
                if (r.Disease.Source != SdrfDraftSource.NotAvailable)
                    characteristics["characteristics[disease]"] = r.Disease.Value == "normal" ? Normal : r.Disease.Term!;

                var factors = draft.FactorColumns.Zip(r.Factors)
                    .Where(x => x.Second.Source != SdrfDraftSource.NotAvailable)
                    .ToDictionary(x => x.First, x => x.Second.Value, StringComparer.Ordinal);

                return new SdrfRowInput(
                    new SdrfSample
                    {
                        SourceName = r.SourceName.Value,
                        Organism = r.Organism.Term,
                        Characteristics = characteristics,
                        BiologicalReplicate = int.Parse(r.BiologicalReplicate.Value, System.Globalization.CultureInfo.InvariantCulture),
                        Label = labelFree,
                        FactorValues = factors,
                        Comments = comments
                    },
                    new SdrfAssay
                    {
                        DataFileName = r.DataFile,
                        AssayName = "run " + SdrfFileNamePattern.Stem(r.DataFile),
                        Instrument = r.Instrument.Term,
                        Fraction = int.Parse(r.Fraction.Value, System.Globalization.CultureInfo.InvariantCulture),
                        TechnicalReplicate = int.Parse(r.TechnicalReplicate.Value, System.Globalization.CultureInfo.InvariantCulture)
                    });
            }).ToList();

            // A draft states only what it read; everything else is "not available", never a refusal (D27's
            // exception to D17, carried over to drafted SDRFs by D30).
            return SdrfBuilder.Build(rows, new SdrfBuilderOptions
            {
                RequireSampleMetadata = false,
                ProteomeXchangeAccession = proteomeXchangeAccession
            });
        }

        // PATO's "normal": the control arm's disease cell, and PRIDE's "Disease free", as a term, so the
        // disease column never mixes terms with free text.
        private static readonly CvParam Normal = new("PATO", "PATO:0000461", "normal", "");

        private static string SourceWord(SdrfDraftSource source) => source switch
        {
            SdrfDraftSource.PrideProjectRecord => "pride project record",
            _ => "inferred"
        };

        private sealed record Replicates(string Key, string KeyWhy, int? Bio, string BioWhy, int? Tech, string TechWhy, int? Frac, string FracWhy);

        private static Replicates ReadReplicates(string file, SdrfFileNameReading t, SdrfReplicateReading m, SdrfReplicates all, bool found,
            bool ranked)
        {
            string key = found ? "T:" + t.SampleKey : "F:" + SdrfFileNamePattern.Stem(file);
            string keyWhy = found ? "files the names read as one sample" : "no structure in the names, so each file is its own sample";
            int? bio = t.BiologicalReplicate ?? t.Replicate, tech = t.TechnicalReplicate, frac = t.Fraction;
            string bioWhy = bio == null ? ""
                : ranked ? "a replicate count in the file names, ranked from 1 within each condition because the names do not count from 1 in every condition"
                : "a replicate count in the file names";
            string techWhy = tech != null ? "a re-injection or technical count in the file names" : "no re-injection marked";
            string fracWhy = frac != null ? "a fraction index in the file names" : "no fraction marked";
            if (m.Number != null)
            {
                var kind = m.Kind;
                if (all.SingleBiologicalSample && kind is SdrfReplicateKind.Unstated or SdrfReplicateKind.Biological)
                    kind = SdrfReplicateKind.Technical;
                if (m.Outer != null) { bio = m.Outer; bioWhy = "the outer of two replicate markers in the names"; }
                switch (kind)
                {
                    case SdrfReplicateKind.Technical:
                        tech = m.Number; techWhy = $"a replicate marker the record or the name says is technical ({all.MarkerEvidence})";
                        key = $"M:{m.Base}#{m.Outer}"; keyWhy = "re-injections of one base name";
                        break;
                    case SdrfReplicateKind.Fraction:
                        frac ??= m.Number; fracWhy = "a fraction marker in the names";
                        key = $"M:{m.Base}#{m.Outer}"; keyWhy = "fractions of one base name";
                        break;
                    default:
                        if (m.Outer != null) { tech = m.Number; techWhy = "the inner of two replicate markers in the names"; }
                        else { bio = m.Number; bioWhy = $"a replicate marker within its base name ({all.MarkerEvidence})"; }
                        key = $"M:{m.Base}#{m.Outer ?? m.Number}"; keyWhy = "a replicate of its base name";
                        break;
                }
            }
            return new Replicates(key, keyWhy, bio, bioWhy, tech, techWhy, frac, fracWhy);
        }

        /// <summary>A project fact, only when PRIDE lists exactly one value for it (D27).</summary>
        private static SdrfDraftCell One(List<CvParam> terms, string what, Func<CvParam, CvParam> normalise)
        {
            var distinct = terms.GroupBy(t => t.Accession, StringComparer.OrdinalIgnoreCase).Select(g => g.First()).ToList();
            if (distinct.Count != 1)
                return SdrfDraftCell.NotAvailable(distinct.Count == 0
                    ? $"PRIDE's project record lists no {what}"
                    : $"PRIDE's project record lists {distinct.Count} values for {what}, so none is written per sample");
            var term = normalise(distinct[0]);
            return new SdrfDraftCell(term.Name, SdrfDraftSource.PrideProjectRecord, $"the one {what} PRIDE's project record lists", term);
        }

        /// <summary>
        /// PRIDE labels organisms NEWT with a bare taxon (<c>NEWT:9606</c>, "Homo sapiens (human)"); SDRF writes
        /// <c>NCBITaxon:9606</c> and the name alone, lower case.
        /// </summary>
        private static CvParam Organism(CvParam t)
        {
            string accession = t.Accession.StartsWith("NEWT:", StringComparison.OrdinalIgnoreCase)
                ? "NCBITaxon:" + t.Accession["NEWT:".Length..]
                : t.Accession;
            string name = Parenthetical.Replace(t.Name, "").Trim().ToLowerInvariant();
            return t with { CvLabel = "NCBITaxon", Accession = accession, Name = name };
        }
    }
}
