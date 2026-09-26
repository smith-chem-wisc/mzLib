using System.Globalization;
using MzLibUtil;

namespace Readers
{
    /// <summary>A place where the draft read something different from what the depositor wrote. Never applied.</summary>
    /// <param name="DataFile">The row's <c>comment[data file]</c>.</param>
    /// <param name="Column">The column the two disagree in.</param>
    /// <param name="Deposited">What the deposited SDRF says, kept.</param>
    /// <param name="Drafted">What the draft read.</param>
    /// <param name="Evidence">Why the draft read it, so a reviewer can decide.</param>
    internal sealed record SdrfDisagreement(string DataFile, string Column, string Deposited, string Drafted, string Evidence);

    /// <summary>The improved document, what changed, and every disagreement left for a person.</summary>
    internal sealed record SdrfImprovement(
        SdrfDocument Document,
        IReadOnlyList<SdrfDisagreement> Disagreements,
        int FilledCells,
        int AddedColumns,
        int AddedRows);

    /// <summary>
    /// Improves a DEPOSITED SDRF with a draft of the same deposit (sdrf D35): the deposited SDRF is the
    /// base, and the draft only fills what it leaves open.
    ///
    /// <para><b>The rules.</b> A stated value is never changed: where the draft reads something different,
    /// that is a <see cref="SdrfDisagreement"/> for a person, because on good curated SDRFs the draft was
    /// judged right 7% of the time against the curated 93% (sdrf benchmark, fresh set 1). A GAP -- an empty
    /// cell or <c>not available</c> -- is filled from the draft and marked in
    /// <c>comment[&lt;characteristic&gt; source]</c>; <c>not applicable</c>, <c>anonymized</c> and
    /// <c>pooled</c> are statements and stay. A column the deposit lacks is added in its block and filled.
    /// A raw file the deposit does not list gets a drafted row. Every deposited row is marked
    /// <c>comment[characteristics source] = deposited</c> unless it already says otherwise.</para>
    ///
    /// <para><b>Channels.</b> A file with several rows is multiplexed: one row per channel. The draft reads
    /// files, not channels, so on such a file it fills only the per-file facts (fraction, technical
    /// replicate, instrument) and never a sample-level one.</para>
    ///
    /// <para><b>Factors.</b> A deposit that states its own factor columns keeps them untouched -- their
    /// names say what was compared, and the draft's cannot be matched to them. A deposit with none gets
    /// the draft's.</para>
    ///
    /// <para>The deposited table is edited in place, not rebuilt, so every column the depositor wrote keeps
    /// its position and spelling. Internal (D19).</para>
    /// </summary>
    internal static class SdrfImprover
    {
        private const string DataFile = "comment[data file]";
        private const string SourceName = "source name";
        private const string RowSource = "comment[characteristics source]";
        private const string BiologicalReplicate = "characteristics[biological replicate]";
        private const string AssayName = "assay name";
        private const string TechnologyType = "technology type";
        private const string TechnologyTypeValue = "proteomic profiling by mass spectrometry";

        private sealed record Target(string Column, string Name, bool PerFile, bool Characteristic, Func<SdrfDraftRow, SdrfDraftCell> Cell);

        private static readonly Target[] Targets =
        {
            new("characteristics[organism]", "organism", false, true, r => r.Organism),
            new("characteristics[organism part]", "organism part", false, true, r => r.OrganismPart),
            new("characteristics[disease]", "disease", false, true, r => r.Disease),
            new(BiologicalReplicate, "biological replicate", false, true, r => r.BiologicalReplicate),
            new("comment[instrument]", "instrument", true, false, r => r.Instrument),
            new("comment[fraction identifier]", "fraction identifier", true, false, r => r.Fraction),
            new("comment[technical replicate]", "technical replicate", true, false, r => r.TechnicalReplicate),
        };

        private static readonly CvParam Normal = new("PATO", "PATO:0000461", "normal", "");

        /// <summary>Columns that describe the whole assay, not one file or one sample: the only ones a drafted row
        /// may copy from the deposit. The builder writes each of them once per document's settings.</summary>
        private static readonly HashSet<string> AssayWide = new(StringComparer.Ordinal)
        {
            "technology type", "comment[label]", "comment[cleavage agent details]", "comment[modification parameters]",
            "comment[precursor mass tolerance]", "comment[fragment mass tolerance]", "comment[proteomics data acquisition method]",
            "comment[dissociation method]", "comment[proteomexchange accession number]"
        };

        /// <summary>
        /// Improves <paramref name="deposited"/> with <paramref name="draft"/>. Throws on a null argument, or
        /// on a deposited SDRF with no <c>comment[data file]</c> column, which leaves nothing to join on.
        /// </summary>
        public static SdrfImprovement Improve(SdrfDocument deposited, SdrfDraft draft)
        {
            if (deposited == null) throw new ArgumentNullException(nameof(deposited));
            if (draft == null) throw new ArgumentNullException(nameof(draft));
            var header = deposited.Header.ToList();
            if (!header.Contains(DataFile))
                throw new ArgumentException("The deposited SDRF has no comment[data file] column, so no row can be joined to a raw file.", nameof(deposited));

            var rows = deposited.Results.Select(r => r.Cells.Concat(Enumerable.Repeat("", Math.Max(0, header.Count - r.Cells.Count))).Take(header.Count).ToList()).ToList();
            var added = new List<List<string>>();
            int addedColumns = 0, filled = 0;
            var disagreements = new List<SdrfDisagreement>();

            int Col(string name) => header.IndexOf(name);
            int Ensure(string name)
            {
                int at = Col(name);
                if (at >= 0) return at;
                at = InsertAt(header, name);
                header.Insert(at, name);
                foreach (var r in rows.Concat(added)) r.Insert(at, "");
                addedColumns++;
                return at;
            }

            var draftByStem = draft.Rows.GroupBy(r => SdrfFileNamePattern.Stem(r.DataFile), StringComparer.OrdinalIgnoreCase)
                .ToDictionary(g => g.Key, g => g.First(), StringComparer.OrdinalIgnoreCase);
            string StemOf(List<string> row) => SdrfFileNamePattern.Stem(row[Col(DataFile)]);
            // Multiplexed = several rows naming ONE file. By the full name, not the stem: a run listed under two
            // extensions (S1.raw, S1.mzML) is two label-free rows, not two channels.
            var rowsPerFile = rows.GroupBy(r => r[Col(DataFile)], StringComparer.OrdinalIgnoreCase).ToDictionary(g => g.Key, g => g.Count(), StringComparer.OrdinalIgnoreCase);
            bool Multiplexed(List<string> row) => rowsPerFile[row[Col(DataFile)]] > 1;

            // ---- deposited rows: mark them, fill gaps, report disagreements ----
            int rowSource = Ensure(RowSource);
            foreach (var r in rows) if (IsGap(r[rowSource])) r[rowSource] = "deposited";
            foreach (var target in Targets)
            {
                bool any = rows.Any(r => draftByStem.TryGetValue(StemOf(r), out var d)
                    && target.Cell(d).Source is not (SdrfDraftSource.NotAvailable or SdrfDraftSource.Default)
                    && (target.PerFile || !Multiplexed(r)));
                if (any) Ensure(target.Column);
            }
            // A column the deposit writes only as free text gets its fills as free text too, not a term beside
            // free text (SdrfDriftLint's MixedTermAndFreeText). A new or term-written column gets terms.
            var freeText = Targets.ToDictionary(t => t.Column, t =>
            {
                var stated = Col(t.Column) < 0 ? new List<string>() : rows.Select(r => r[Col(t.Column)]).Where(v => !IsGap(v)).ToList();
                return stated.Count > 0 && stated.All(v => !SdrfCell.IsTerm(v));
            });
            string Write(Target target, SdrfDraftCell cell) => Written(target, cell, freeText[target.Column]);
            foreach (var r in rows)
            {
                if (!draftByStem.TryGetValue(StemOf(r), out var d)) continue;
                bool multiplexed = Multiplexed(r);
                foreach (var target in Targets)
                {
                    if (!target.PerFile && multiplexed) continue;
                    int at = Col(target.Column);
                    var cell = target.Cell(d);
                    // A default is not a reading: it neither fills nor disputes.
                    if (at < 0 || cell.Source is SdrfDraftSource.NotAvailable or SdrfDraftSource.Default) continue;
                    if (IsGap(r[at]))
                    {
                        r[at] = Write(target, cell);
                        int mark = Ensure($"comment[{target.Name} source]");
                        r[mark] = Word(cell.Source);
                        filled++;
                    }
                    // A project-level summary is weaker evidence than a per-file statement: it may fill a gap,
                    // never dispute (PRIDE's "LTQ Orbitrap" against a curated per-file "Q Exactive").
                    else if (cell.Source != SdrfDraftSource.PrideProjectRecord && !Same(r[at], cell))
                        disagreements.Add(new SdrfDisagreement(r[Col(DataFile)], target.Column, r[at], Write(target, cell), cell.Evidence));
                }
            }

            // ---- factors: the draft's only where the deposit states none ----
            bool draftFactors = !header.Any(h => h.StartsWith("factor value[", StringComparison.OrdinalIgnoreCase)) && draft.FactorColumns.Count > 0;
            if (draftFactors)
            {
                foreach (var f in draft.FactorColumns) { header.Add(f); foreach (var r in rows) r.Add(""); addedColumns++; }
                foreach (var r in rows)
                {
                    if (!draftByStem.TryGetValue(StemOf(r), out var d) || Multiplexed(r)) continue;
                    for (int k = 0; k < draft.FactorColumns.Count; k++)
                        if (d.Factors[k].Source != SdrfDraftSource.NotAvailable) { r[Col(draft.FactorColumns[k])] = d.Factors[k].Value; filled++; }
                }
            }

            // ---- every column the specification requires (G29) ----
            // A deposit missing one is not improved where it most visibly could be. Added in its block,
            // "not available" on every row -- except technology type, whose one specified value every
            // mass-spectrometry deposit has. Added before the drafted rows, so they carry it too.
            foreach (var required in SdrfValidator.RequiredColumns.Where(c => Col(c) < 0))
            {
                int at = Ensure(required);
                foreach (var r in rows) r[at] = required == TechnologyType ? TechnologyTypeValue : SdrfReserved.NotAvailable;
            }

            // ---- raw files the deposit does not list ----
            var listed = new HashSet<string>(rows.Select(StemOf), StringComparer.OrdinalIgnoreCase);
            var depositedNames = new HashSet<string>(rows.Select(r => r[Col(SourceName) < 0 ? 0 : Col(SourceName)]), StringComparer.Ordinal);
            // An ASSAY-WIDE column every deposited row fills with one value is carried to drafted rows (a label,
            // a cleavage agent). Nothing else: a one-row deposit makes every column "constant", and a file URI
            // or a sample's individual is not the drafted file's. By POSITION: SDRF repeats columns.
            var constant = Enumerable.Range(0, header.Count).Select(i =>
            {
                var values = rows.Select(r => r[i]).ToList();
                return AssayWide.Contains(header[i]) && values.Count > 0 && values.All(v => !IsGap(v)) && values.Distinct(StringComparer.Ordinal).Count() == 1
                    ? values[0] : SdrfReserved.NotAvailable;
            }).ToList();
            var usedAssays = new HashSet<string>(Col("assay name") >= 0 ? rows.Select(r => r[Col("assay name")]) : Enumerable.Empty<string>(), StringComparer.Ordinal);
            foreach (var d in draft.Rows.Where(d => !listed.Contains(SdrfFileNamePattern.Stem(d.DataFile))))
            {
                var row = constant.ToList();
                if (Col(SourceName) >= 0)
                    row[Col(SourceName)] = depositedNames.Contains(d.SourceName.Value) ? d.SourceName.Value + " (drafted)" : d.SourceName.Value;
                if (Col("assay name") >= 0)
                {
                    // Unique within the document: one run listed under two extensions (X.raw, X.mzXML) would
                    // otherwise give two rows one key (PXD001587). The full file name is the fallback.
                    string assay = "run " + SdrfFileNamePattern.Stem(d.DataFile);
                    if (!usedAssays.Add(assay)) { assay = "run " + d.DataFile; usedAssays.Add(assay); }
                    row[Col("assay name")] = assay;
                }
                row[Col(DataFile)] = d.DataFile;
                row[Col(RowSource)] = "inferred";
                foreach (var target in Targets)
                    if (Col(target.Column) >= 0)
                        row[Col(target.Column)] = target.Cell(d).Source == SdrfDraftSource.NotAvailable ? SdrfReserved.NotAvailable : Write(target, target.Cell(d));
                foreach (var h in header.Where(h => h.EndsWith(" source]", StringComparison.Ordinal) && h != RowSource))
                    row[Col(h)] = SdrfReserved.NotApplicable;
                for (int k = 0; k < draft.FactorColumns.Count && draftFactors; k++)
                    row[Col(draft.FactorColumns[k])] = d.Factors[k].Value;
                if (!draftFactors)
                    foreach (var h in header.Where(h => h.StartsWith("factor value[", StringComparison.OrdinalIgnoreCase)))
                        row[Col(h)] = SdrfReserved.NotAvailable;
                added.Add(row);
            }

            // An override column's empty cell means "no override on this row" -- not applicable, as the
            // builder writes it. Only columns this call added are touched.
            foreach (var h in header.Where(h => h.EndsWith(" source]", StringComparison.Ordinal) && h != RowSource && !deposited.Header.Contains(h)))
                foreach (var r in rows) if (r[Col(h)].Length == 0) r[Col(h)] = SdrfReserved.NotApplicable;
            foreach (var r in rows)
                for (int i = 0; i < header.Count; i++)
                    if (r[i].Length == 0 && !deposited.Header.Contains(header[i])) r[i] = SdrfReserved.NotAvailable;

            var newHeader = new SdrfHeader(header);
            var document = new SdrfDocument(newHeader, rows.Concat(added).Select(r => new SdrfRow(newHeader, r)));
            return new SdrfImprovement(document, disagreements, filled, addedColumns, added.Count);
        }

        private static bool IsGap(string? cell) =>
            string.IsNullOrWhiteSpace(cell) || string.Equals(cell.Trim(), SdrfReserved.NotAvailable, StringComparison.OrdinalIgnoreCase);

        private static string Word(SdrfDraftSource source) =>
            source == SdrfDraftSource.PrideProjectRecord ? "pride project record" : "inferred";

        private static string Written(Target target, SdrfDraftCell cell, bool freeText)
        {
            if (freeText) return cell.Value;
            if (target.Name == "disease" && cell.Value == "normal") return SdrfCell.ToCell(Normal);
            return cell.Term != null ? SdrfCell.ToCell(cell.Term) : cell.Value;
        }

        /// <summary>The same statement, however it is spelled: a number by value, a term by accession or name.</summary>
        private static bool Same(string deposited, SdrfDraftCell drafted)
        {
            string d = deposited.Trim();
            if (int.TryParse(d, NumberStyles.None, CultureInfo.InvariantCulture, out int a)
                && int.TryParse(drafted.Value, NumberStyles.None, CultureInfo.InvariantCulture, out int b))
                return a == b;
            string name = d, accession = "";
            if (SdrfCell.TryParseTerm(d, out var term)) { name = term.Name; accession = term.Accession; }
            if (accession.Length > 0 && drafted.Term != null && string.Equals(accession, drafted.Term.Accession, StringComparison.OrdinalIgnoreCase))
                return true;
            return string.Equals(name, drafted.Value, StringComparison.OrdinalIgnoreCase);
        }

        /// <summary>Where a new column belongs: a characteristic before the biological replicate, a comment before the factors.</summary>
        private static int InsertAt(List<string> header, string name)
        {
            // The order SdrfBuilder writes: source name, characteristics..., biological replicate, assay name,
            // technology type, the comment block, the factors.
            if (name == SourceName) return 0;
            if (name == AssayName)
            {
                int lastCharacteristic = header.FindLastIndex(h => h.StartsWith("characteristics[", StringComparison.OrdinalIgnoreCase));
                return lastCharacteristic >= 0 ? lastCharacteristic + 1 : Math.Min(1, header.Count);
            }
            if (name == TechnologyType)
            {
                int assay = header.IndexOf(AssayName);
                if (assay >= 0) return assay + 1;
            }
            if (name.StartsWith("characteristics[", StringComparison.OrdinalIgnoreCase))
            {
                int bio = header.IndexOf(BiologicalReplicate);
                if (bio >= 0 && name != BiologicalReplicate) return bio;
                int last = header.FindLastIndex(h => h.StartsWith("characteristics[", StringComparison.OrdinalIgnoreCase));
                return last >= 0 ? last + 1 : Math.Min(1, header.Count);
            }
            int factor = header.FindIndex(h => h.StartsWith("factor value[", StringComparison.OrdinalIgnoreCase));
            return factor >= 0 ? factor : header.Count;
        }
    }
}
