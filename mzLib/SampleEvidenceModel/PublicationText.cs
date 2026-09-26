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
        /// Supplement tables as text, one line per row prefixed with the row's locator
        /// (<c>[mmc1.xlsx!S1!R14]</c>), capped per table and overall so one result table cannot crowd out the rest.
        /// Result tables (proteins, statistics) are left out entirely.
        /// </summary>
        internal static string Tables(IEnumerable<SupplementTable> tables, int maxRowsPerTable = 150, int maxChars = 120_000)
        {
            var sb = new StringBuilder();
            foreach (var t in tables)
            {
                if (sb.Length >= maxChars) { sb.Append("[further tables omitted for length]\n"); break; }
                string name = t.Sheet.Length > 0 ? $"{t.File}!{t.Sheet}" : t.File;
                sb.Append("### ").Append(name);
                if (t.Title.Length > 0) sb.Append(" -- ").Append(t.Title);
                sb.Append('\n').Append("[header] ").Append(string.Join(" | ", t.Header)).Append('\n');
                for (int k = 0; k < Math.Min(t.Rows.Count, maxRowsPerTable); k++)
                    sb.Append('[').Append(name).Append("!R").Append(t.RowNumbers[k]).Append("] ").Append(string.Join(" | ", t.Rows[k])).Append('\n');
                if (t.Rows.Count > maxRowsPerTable) sb.Append($"[{t.Rows.Count - maxRowsPerTable} more rows]\n");
                sb.Append('\n');
            }
            return sb.ToString();
        }

        internal static string Clean(string s) => Regex.Replace(s, @"\s+", " ").Trim();
    }
}
