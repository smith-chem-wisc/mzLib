using System.Text;
using System.Text.RegularExpressions;
using System.Xml.Linq;
using Readers;

namespace SampleEvidenceModel
{
    /// <summary>
    /// The text a reader is given about a deposit: its paper's methods and tables (from JATS full text) and its
    /// supplement tables, each row labelled so a claim can cite it. Measured choices (sdrf results/benchmark,
    /// 2026-09-25/26): methods-like sections and tables only -- half the false hits of a whole-paper reading came
    /// from the bibliography and the author affiliations.
    /// </summary>
    internal static class PublicationText
    {
        private static readonly Regex MethodsTitle = new(
            @"method|experimental|materials|procedure|sample prep|patients|subjects|animals|cohort|study design|participants",
            RegexOptions.IgnoreCase | RegexOptions.Compiled);

        /// <summary>The methods-like sections and every table (caption and cells) of a JATS article, as plain text.</summary>
        internal static string MethodsAndTables(string jatsXml)
        {
            if (jatsXml == null) throw new ArgumentNullException(nameof(jatsXml));
            var doc = XDocument.Parse(jatsXml);
            var body = doc.Descendants().FirstOrDefault(e => e.Name.LocalName == "body");
            if (body == null) return string.Empty;
            var sb = new StringBuilder();
            foreach (var sec in body.Elements().Where(e => e.Name.LocalName == "sec"))
            {
                string title = sec.Elements().FirstOrDefault(e => e.Name.LocalName == "title")?.Value ?? "";
                if (MethodsTitle.IsMatch(title) || ((string?)sec.Attribute("sec-type") ?? "").Contains("method", StringComparison.OrdinalIgnoreCase))
                    sb.Append("## ").Append(Clean(title)).Append('\n').Append(Clean(sec.Value)).Append("\n\n");
            }
            foreach (var table in doc.Descendants().Where(e => e.Name.LocalName == "table-wrap"))
                sb.Append("## Table\n").Append(Clean(table.Value)).Append("\n\n");
            return sb.ToString();
        }

        /// <summary>
        /// Supplement tables as text, one line per row prefixed with the row's locator (<c>[mmc1.xlsx!S1!R14]</c>).
        /// Tables that describe samples come first and are given whole (up to <paramref name="maxSampleRows"/>); other
        /// tables follow, capped at <paramref name="maxRowsPerTable"/>; result tables (proteins, statistics) come last,
        /// as a preview of <paramref name="resultRows"/> rows. Measured (sdrf G36 batch 2, 2026-09-26): in document
        /// order, protein tables filled the budget before a cohort table was reached, and an ISA assay table of 160
        /// runs lost its last 10.
        /// </summary>
        internal static string Tables(IEnumerable<SupplementTable> tables, int maxRowsPerTable = 150, int maxChars = 200_000,
            int maxSampleRows = 2_000, int resultRows = 5)
        {
            var sb = new StringBuilder();
            var ranked = tables.Select(t => (t, Rank(t))).OrderBy(x => x.Item2).ToList();
            for (int i = 0; i < ranked.Count; i++)
            {
                var (t, rank) = ranked[i];
                if (sb.Length >= maxChars) { sb.Append("[further tables omitted for length]\n"); break; }
                string name = t.Sheet.Length > 0 ? $"{t.File}!{t.Sheet}" : t.File;
                sb.Append("### ").Append(name);
                if (t.Title.Length > 0) sb.Append(" -- ").Append(t.Title);
                if (rank == ResultRank) sb.Append(" [result table: preview only]");
                sb.Append('\n').Append("[header] ").Append(string.Join(" | ", t.Header)).Append('\n');
                int cap = rank switch { SampleRank => maxSampleRows, ResultRank => resultRows, _ => maxRowsPerTable };
                // Each table still to come keeps a small reserve, so a long table cannot hide the ones after it.
                int stop = maxChars - Math.Min(maxChars / 2, ReservePerTable * (ranked.Count - i - 1));
                int k = 0;
                for (; k < Math.Min(t.Rows.Count, cap) && sb.Length < stop; k++)
                    sb.Append('[').Append(name).Append("!R").Append(t.RowNumbers[k]).Append("] ").Append(string.Join(" | ", t.Rows[k])).Append('\n');
                if (t.Rows.Count > k) sb.Append($"[{t.Rows.Count - k} more rows]\n");
                sb.Append('\n');
            }
            return sb.ToString();
        }

        private const int ReservePerTable = 2_000;

        private const int SampleRank = 0, OtherRank = 1, ResultRank = 2;

        private static readonly Regex SampleTitle = new(@"sample|patient|subject|donor|cohort|participant|clinical|demograph|characteristic|channel|\btmt\b|itraq|label|design|run ?order|raw ?file|s_[^ ]*\.txt|a_[^ ]*\.txt",
            RegexOptions.IgnoreCase | RegexOptions.Compiled);

        /// <summary>
        /// Sample tables first: an ISA study or assay file, or a table with two or more sample-describing headers, or a
        /// sample-describing name with at least one. Result tables last (<see cref="SampleEvidenceExtractor.IsResultTable"/>).
        /// </summary>
        private static int Rank(SupplementTable t)
        {
            if (SampleEvidenceExtractor.IsResultTable(t)) return ResultRank;
            if (FigureSheet.IsMatch($"{t.Sheet} {t.Title}")) return OtherRank;
            // Distinct columns, not headers: "fraction_helix | fraction_sheet | fraction_coil" or "Cluster 1..5" are one
            // measurement repeated, and they made a figure's data outrank the cohort table (G36 batch 2, PXD017291).
            var described = t.Header.Select(HeaderMap.ColumnFor).Where(c => c != null).Distinct().Count();
            int result = t.Header.Count(h => ResultWord.IsMatch(h));
            if (result >= described) return OtherRank;
            return described >= 2 || (described >= 1 && SampleTitle.IsMatch($"{t.File} {t.Sheet} {t.Title}")) ? SampleRank : OtherRank;
        }

        private static readonly Regex FigureSheet = new(@"\bfig", RegexOptions.IgnoreCase | RegexOptions.Compiled);
        private static readonly Regex ResultWord = new(@"p[ ._-]?val|q[ ._-]?val|p\.adj|fdr|fold|log2|ratio|intensity|abundance|protein|peptide|gene|uniprot|score|enrich|pathway",
            RegexOptions.IgnoreCase | RegexOptions.Compiled);

        private static readonly Regex Digits = new(@"\d+", RegexOptions.Compiled);

        /// <summary>
        /// The raw file names, one per line when they fit in <paramref name="maxChars"/>. Otherwise each group of names
        /// that differ only in their numbers is written once as a pattern with its numbers listed, so no file is lost
        /// from a 650-file deposit (sdrf G36 batch 2: a 600-name cap cut PXD011967's last 50). A claim still names the
        /// file in full, and <see cref="ModelEvidenceReader.Interpret"/> still requires it to be one of the deposit's.
        /// </summary>
        internal static string Files(IReadOnlyList<string> files, int maxChars = 60_000)
        {
            var sb = new StringBuilder();
            if (files.Sum(f => f.Length + 1) <= maxChars)
            {
                foreach (var f in files) sb.Append(f).Append('\n');
                return sb.ToString();
            }
            sb.Append("Listed by pattern: each {n} stands for the numbers shown, in order, one file per line of numbers.\n");
            foreach (var g in files.GroupBy(f => Digits.Replace(f, "{n}")))
            {
                if (g.Count() == 1) { sb.Append(g.First()).Append('\n'); continue; }
                sb.Append(g.Key).Append($"  ({g.Count()} files)\n  ");
                sb.Append(string.Join(", ", g.Select(f => string.Join("/", Digits.Matches(f).Select(m => m.Value))))).Append('\n');
            }
            return sb.ToString();
        }

        internal static string Clean(string s) => Regex.Replace(s, @"\s+", " ").Trim();
    }
}
