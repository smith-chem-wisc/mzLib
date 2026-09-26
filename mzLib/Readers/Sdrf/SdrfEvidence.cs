using System.Text;
using System.Text.RegularExpressions;
using MzLibUtil;

namespace Readers
{
    /// <summary>How far a piece of publication evidence can be trusted. The drafter applies Certain and Likely; a
    /// Guess is kept for review only.</summary>
    internal enum SdrfEvidenceConfidence
    {
        /// <summary>Read by a deterministic rule and validated (an ISA-Tab or SDRF supplement).</summary>
        Certain,

        /// <summary>Read by a rule that can be wrong on an odd table (a table keyed by file name, a channel map).</summary>
        Likely,

        /// <summary>Plausible but unconfirmed. Never written into a draft.</summary>
        Guess
    }

    /// <summary>
    /// One claim about a deposit's samples, read from its paper or supplementary files (sdrf design
    /// SAMPLE-EVIDENCE.md; D41: read in mzLib, D42: rules first, a model reader only when the caller opts in).
    /// </summary>
    /// <param name="DataFile">The raw file the claim is about, as PRIDE names it; empty for the whole deposit.</param>
    /// <param name="Label">The isobaric channel (<c>TMT127N</c>, bare), or empty.</param>
    /// <param name="Column">The SDRF column exactly as an SDRF writes it, e.g. <c>characteristics[age]</c>.</param>
    /// <param name="Value">The value as it will be written, units and codes already resolved (<c>77Y</c>, <c>female</c>).</param>
    /// <param name="Source">Where it was read: <c>supplement</c>, <c>paper</c> or <c>pride record</c>.</param>
    /// <param name="Locator">Where exactly, so a reviewer can check it: <c>mmc2.xlsx!Sheet1!R14C3</c>, a section, a quote.</param>
    /// <param name="Method">How it was read: <c>isa-tab</c>, <c>sdrf</c>, <c>file-key</c>, <c>channel-map</c>, <c>text</c>, <c>model</c>.</param>
    /// <param name="Confidence">See <see cref="SdrfEvidenceConfidence"/>.</param>
    /// <param name="DataFilePattern">
    /// A set of the deposit's raw files the claim holds for, as a glob over their names (<c>*TMT6*</c>, <c>*_SetA_*</c>;
    /// <c>*</c> any run of characters, <c>?</c> one, case-insensitive), when it is neither one file nor all of them.
    /// Empty otherwise; <see cref="DataFile"/> is then empty too.
    /// </param>
    internal sealed record SdrfEvidence(
        string DataFile,
        string Label,
        string Column,
        string Value,
        string Source,
        string Locator,
        string Method,
        SdrfEvidenceConfidence Confidence,
        string DataFilePattern = "")
    {
        /// <summary>Whether a raw file name matches a <see cref="DataFilePattern"/> glob (the whole name, ignoring case).</summary>
        internal static bool GlobMatches(string pattern, string fileName) =>
            Regex.IsMatch(fileName, "^" + Regex.Escape(pattern).Replace(@"\*", ".*").Replace(@"\?", ".") + "$",
                RegexOptions.IgnoreCase | RegexOptions.CultureInvariant);
    }

    /// <summary>
    /// A claim the drafter did not apply, and why: it disagrees with a reading the drafter made, conflicts with another
    /// claim, or has no place in a draft yet.
    /// </summary>
    internal sealed record SdrfEvidenceNote(string DataFile, string Column, string DraftValue, string EvidenceValue, string Why);

    /// <summary>
    /// The serialised form of <see cref="SdrfEvidence"/>: tab-separated, UTF-8, a header naming every column. It is how
    /// evidence is reviewed, cached and used as a test fixture.
    /// </summary>
    internal static class SdrfEvidenceFile
    {
        internal static readonly string[] Columns = { "data file", "label", "column", "value", "source", "locator", "method", "confidence" };

        /// <summary>Written by every writer, optional to a reader: files written before it existed still read.</summary>
        internal const string PatternColumn = "data file pattern";

        /// <summary>Reads an evidence file. Throws <see cref="MzLibException"/> when a column is missing or a confidence is unknown.</summary>
        internal static IReadOnlyList<SdrfEvidence> Read(string path)
        {
            if (path == null) throw new ArgumentNullException(nameof(path));
            var lines = File.ReadAllLines(path, Encoding.UTF8);
            if (lines.Length == 0) throw new MzLibException($"Evidence file is empty: '{path}'");
            var header = lines[0].Split('\t').Select(h => h.Trim().ToLowerInvariant()).ToList();
            var missing = Columns.Where(c => !header.Contains(c)).ToList();
            if (missing.Count > 0)
                throw new MzLibException($"Evidence file '{path}' lacks the column(s) {string.Join(", ", missing)}.");
            int At(string c) => header.IndexOf(c);
            var claims = new List<SdrfEvidence>();
            for (int i = 1; i < lines.Length; i++)
            {
                if (lines[i].Length == 0) continue;
                var cells = lines[i].Split('\t');
                string Get(string c) => At(c) < cells.Length ? cells[At(c)] : "";
                if (!Enum.TryParse(Get("confidence"), ignoreCase: true, out SdrfEvidenceConfidence confidence))
                    throw new MzLibException($"Evidence file '{path}' line {i + 1}: unknown confidence '{Get("confidence")}'.");
                claims.Add(new SdrfEvidence(Get("data file"), Get("label"), Get("column"), Get("value"), Get("source"),
                    Get("locator"), Get("method"), confidence, At(PatternColumn) >= 0 ? Get(PatternColumn) : ""));
            }
            return claims;
        }

        /// <summary>Writes an evidence file. Tabs and line breaks inside a value become spaces.</summary>
        internal static void Write(string path, IEnumerable<SdrfEvidence> claims)
        {
            if (path == null) throw new ArgumentNullException(nameof(path));
            if (claims == null) throw new ArgumentNullException(nameof(claims));
            static string Clean(string s) => (s ?? "").Replace('\t', ' ').Replace('\r', ' ').Replace('\n', ' ');
            var sb = new StringBuilder(string.Join('\t', Columns)).Append('\t').Append(PatternColumn).Append('\n');
            foreach (var c in claims)
                sb.Append(string.Join('\t', new[] { c.DataFile, c.Label, c.Column, c.Value, c.Source, c.Locator, c.Method }.Select(Clean)))
                  .Append('\t').Append(c.Confidence.ToString().ToLowerInvariant()).Append('\t').Append(Clean(c.DataFilePattern)).Append('\n');
            File.WriteAllText(path, sb.ToString(), new UTF8Encoding(false));
        }
    }
}
