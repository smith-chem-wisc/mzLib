using System.Globalization;
using System.Text.RegularExpressions;

namespace Readers
{
    /// <summary>
    /// Rule-based sample evidence from a deposit's supplement tables (sdrf design SAMPLE-EVIDENCE.md, E3; D42: rules
    /// first). Ported from the prototype measured on 510 PRIDE-linked papers (sdrf results/benchmark, 2026-09-26).
    ///
    /// <para>Four rules, tried in this order, and a claim a stronger rule made is never repeated by a weaker one:</para>
    /// <list type="bullet">
    ///   <item><b>isa-tab</b> (Certain): an ISA-Tab study (<c>Sample Name</c>, <c>Characteristics[...]</c>) joined to its
    ///   assay's <c>Raw Data File</c>.</item>
    ///   <item><b>sdrf</b> (Certain): a supplement that is itself an SDRF; only rows naming this deposit's files.</item>
    ///   <item><b>channel-map</b> (Likely): an isobaric table giving each plex channel its sample, the plex named as a token
    ///   of the raw-file names -- in its own column, or in one cell with the channel (<c>TMT05_TMT-129</c>).</item>
    ///   <item><b>file-key</b> (Likely): a table with a column equal to the raw-file stems; its other columns mapped to
    ///   SDRF characteristics by header and decoded by it (<c>Gender (1F,0M)</c>, <c>Age (months)</c>).</item>
    /// </list>
    ///
    /// <para><b>What it will not read.</b> A result table (proteins, peptides, statistics), figure source data, a value
    /// that is a code it cannot decode, and anything it would have to guess: those are for a reader that understands
    /// the table, which is opt-in (D42).</para>
    /// </summary>
    internal static class SampleEvidenceExtractor
    {
        private const string Supplement = "supplement";

        /// <summary>The claims the tables support about <paramref name="rawFiles"/>. Never null; empty when nothing links.</summary>
        internal static IReadOnlyList<SdrfEvidence> Extract(IEnumerable<SupplementTable> tables, IReadOnlyList<string> rawFiles)
        {
            if (tables == null) throw new ArgumentNullException(nameof(tables));
            if (rawFiles == null) throw new ArgumentNullException(nameof(rawFiles));
            var list = tables.Where(t => t != null).ToList();
            if (list.Count == 0 || rawFiles.Count == 0) return Array.Empty<SdrfEvidence>();
            var index = new FileIndex(rawFiles);
            var claims = new List<SdrfEvidence>();
            claims.AddRange(IsaTab(list, index));
            claims.AddRange(SdrfTables(list, index));
            foreach (var t in list.Where(t => !IsResultTable(t)))
                claims.AddRange(ChannelMap(t, index));
            foreach (var t in list.Where(t => !IsResultTable(t)))
                claims.AddRange(FileKey(t, index, rawFiles.Count));

            // One claim per file, channel and column: the first rule to make it wins.
            var seen = new HashSet<(string, string, string)>();
            return claims.Where(c => seen.Add((c.DataFile, c.Label, c.Column))).ToList();
        }

        // ---------------- the rules ----------------

        private static readonly Regex IsaStudy = new(@"(^|__|/)s_[^/]*\.txt$", RegexOptions.IgnoreCase | RegexOptions.Compiled);
        private static readonly Regex IsaAssay = new(@"(^|__|/)a_[^/]*\.txt$", RegexOptions.IgnoreCase | RegexOptions.Compiled);
        private static readonly Regex SdrfColumn = new(@"^(characteristics|comment|factor value)\s*\[\s*([^\]]*?)\s*\]$", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        private static IEnumerable<SdrfEvidence> IsaTab(List<SupplementTable> tables, FileIndex index)
        {
            foreach (var study in tables.Where(t => IsaStudy.IsMatch(t.File)))
            {
                int sample = Find(study.Header, "Sample Name");
                if (sample < 0) continue;
                int source = Find(study.Header, "Source Name");
                var characteristics = study.Header.Select((h, j) => (j, m: SdrfColumn.Match(h.Trim())))
                    .Where(x => x.m.Success && x.m.Groups[1].Value.Equals("characteristics", StringComparison.OrdinalIgnoreCase) && x.m.Groups[2].Value.Length > 0)
                    .Select(x => (x.j, Column: $"characteristics[{x.m.Groups[2].Value.ToLowerInvariant()}]")).ToList();
                var rowOf = new Dictionary<string, int>(StringComparer.Ordinal);
                for (int k = 0; k < study.Rows.Count; k++) rowOf.TryAdd(study.Rows[k][sample].Trim(), k);
                foreach (var assay in tables.Where(t => IsaAssay.IsMatch(t.File)))
                {
                    int aSample = Find(assay.Header, "Sample Name"), raw = Find(assay.Header, "Raw Data File");
                    if (aSample < 0 || raw < 0) continue;
                    for (int k = 0; k < assay.Rows.Count; k++)
                    {
                        var files = index.Match(FileStem(assay.Rows[k][raw]));
                        if (files.Count != 1 || !rowOf.TryGetValue(assay.Rows[k][aSample].Trim(), out int s)) continue;
                        string file = files.First();
                        var srow = study.Rows[s];
                        string loc = assay.Locator(k, raw);
                        if (source >= 0 && srow[source].Trim().Length > 0)
                            yield return Claim(file, "", "source name", srow[source].Trim(), loc, "isa-tab", SdrfEvidenceConfidence.Certain);
                        foreach (var (j, column) in characteristics)
                            if (Usable(srow[j]))
                                yield return Claim(file, "", column, srow[j].Trim(), study.Locator(s, j), "isa-tab", SdrfEvidenceConfidence.Certain);
                    }
                }
            }
        }

        private static IEnumerable<SdrfEvidence> SdrfTables(List<SupplementTable> tables, FileIndex index)
        {
            foreach (var t in tables)
            {
                var header = t.Header.Select(h => h.Trim().ToLowerInvariant()).ToList();
                int data = header.FindIndex(h => h.StartsWith("comment[data file", StringComparison.Ordinal));
                if (data < 0 || !header.Contains("source name")) continue;
                for (int k = 0; k < t.Rows.Count; k++)
                {
                    var files = index.Match(FileStem(t.Rows[k][data]));
                    if (files.Count != 1) continue;
                    for (int j = 0; j < header.Count; j++)
                    {
                        if (j == data || !Usable(t.Rows[k][j])) continue;
                        var m = SdrfColumn.Match(header[j]);
                        string column = header[j] == "source name" ? "source name"
                            : m.Success && m.Groups[2].Value.Length > 0 ? $"{m.Groups[1].Value}[{m.Groups[2].Value}]" : "";
                        if (column.Length > 0)
                            yield return Claim(files.First(), "", column, t.Rows[k][j].Trim(), t.Locator(k, j), "sdrf", SdrfEvidenceConfidence.Certain);
                    }
                }
            }
        }

        private static readonly Regex Tag = new(@"^(?:tmt(?:pro)?[ -]?|itraq[ -]?)?(1[23]\d[NC]?|11[3-9]|121)$", RegexOptions.IgnoreCase | RegexOptions.Compiled);
        private static readonly Regex Combined = new(@"^(.+?)[ _-]+(?:tmt(?:pro)?[ -]?|itraq[ -]?)?(1[23]\d[NC]?|11[3-9]|121)$", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        private static IEnumerable<SdrfEvidence> ChannelMap(SupplementTable t, FileIndex index)
        {
            if (FigureData(t)) yield break;
            int sample = SampleColumn(t, -1, -1);
            // (a) plex and channel in one cell
            for (int j = 0; j < t.Header.Count; j++)
            {
                var values = Column(t, j);
                var parsed = values.Select(v => (v.Row, M: Combined.Match(v.Value.Replace(" ", "")))).Where(x => x.M.Success).ToList();
                if (values.Count < 4 || parsed.Count < 0.8 * values.Count) continue;
                int s = SampleColumn(t, j, -1);
                if (s < 0 || !IdentifierColumn(t, s)) continue;
                bool any = false;
                foreach (var (row, m) in parsed)
                    foreach (var c in ChannelClaims(t, index.Match(m.Groups[1].Value), row, m.Groups[2].Value, s, j, -1))
                    { any = true; yield return c; }
                if (any) yield break;
            }
            // (b) a channel column beside a plex column
            int tag = Enumerable.Range(0, t.Header.Count).FirstOrDefault(j =>
            {
                var v = Column(t, j);
                return v.Count >= 4 && v.Count(x => Tag.IsMatch(x.Value.Replace(" ", ""))) >= 0.8 * v.Count;
            }, -1);
            if (tag < 0) yield break;
            int plex = -1, best = 0;
            for (int j = 0; j < t.Header.Count; j++)
            {
                if (j == tag) continue;
                int n = Column(t, j).SelectMany(v => index.Match(v.Value)).Distinct().Count();
                if (n > best) (plex, best) = (j, n);
            }
            sample = SampleColumn(t, tag, plex);
            if (plex < 0 || sample < 0 || !IdentifierColumn(t, sample)) yield break;
            for (int k = 0; k < t.Rows.Count; k++)
            {
                var m = Tag.Match(t.Rows[k][tag].Trim().Replace(" ", ""));
                if (!m.Success) continue;
                foreach (var c in ChannelClaims(t, index.Match(t.Rows[k][plex].Trim()), k, m.Groups[1].Value, sample, tag, plex))
                    yield return c;
            }
        }

        private static IEnumerable<SdrfEvidence> ChannelClaims(SupplementTable t, IReadOnlyCollection<string> files, int row, string tag, int sample, int skipA, int skipB)
        {
            string value = t.Rows[row][sample].Trim();
            if (files.Count == 0 || value.Length == 0) yield break;
            string label = (Regex.IsMatch(tag, "^(11[3-9]|121)$") ? "iTRAQ" : "TMT") + tag.ToUpperInvariant();
            foreach (var file in files.OrderBy(f => f, StringComparer.Ordinal))
            {
                yield return Claim(file, label, "source name", value, t.Locator(row, sample), "channel-map", SdrfEvidenceConfidence.Likely);
                for (int j = 0; j < t.Header.Count; j++)
                {
                    if (j == sample || j == skipA || j == skipB) continue;
                    string? column = HeaderMap.ColumnFor(t.Header[j]);
                    if (column == null || !column.StartsWith("characteristics[", StringComparison.Ordinal)) continue;
                    if (HeaderMap.Decode(column, t.Header[j], t.Rows[row][j]) is { } v)
                        yield return Claim(file, label, column, v, t.Locator(row, j), "channel-map", SdrfEvidenceConfidence.Likely);
                }
            }
        }

        private static IEnumerable<SdrfEvidence> FileKey(SupplementTable t, FileIndex index, int fileCount)
        {
            int key = -1;
            Dictionary<string, int>? rowOf = null;
            for (int j = 0; j < t.Header.Count; j++)
            {
                var hits = new Dictionary<string, List<int>>(StringComparer.Ordinal);
                for (int k = 0; k < t.Rows.Count; k++)
                    foreach (var f in index.ByExactStem(FileStem(t.Rows[k][j])))
                        (hits.TryGetValue(f, out var l) ? l : hits[f] = new List<int>()).Add(k);
                var unique = hits.Where(h => h.Value.Count == 1).ToDictionary(h => h.Key, h => h.Value[0], StringComparer.Ordinal);
                if (unique.Count >= Math.Max(2, 0.3 * fileCount) && (rowOf == null || unique.Count > rowOf.Count))
                    (key, rowOf) = (j, unique);
            }
            if (rowOf == null) yield break;
            foreach (var (file, k) in rowOf.OrderBy(x => x.Key, StringComparer.Ordinal))
                for (int j = 0; j < t.Header.Count; j++)
                {
                    if (j == key) continue;
                    string? column = HeaderMap.ColumnFor(t.Header[j]);
                    if (column == null || !column.StartsWith("characteristics[", StringComparison.Ordinal)) continue;
                    if (HeaderMap.Decode(column, t.Header[j], t.Rows[k][j]) is { } v)
                        yield return Claim(file, "", column, v, t.Locator(k, j), "file-key", SdrfEvidenceConfidence.Likely);
                }
        }

        // ---------------- guards ----------------

        private static readonly Regex UniProt = new(@"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$", RegexOptions.Compiled);
        private static readonly Regex ResultWord = new(@"\b(protein ?ids?|proteins?|peptides?|sequence|gene ?names?|genes?|uniprot|accession|intensity|lfq|ibaq|abundance|ratio|fold|log2\w*|p[ -]?value|q[ -]?value|fdr|pep|score|m/?z|charge|modifications?|sites?|localization|spectral ?counts?|go ?term|pathway|enrichment|unique|razor|coverage|retention|scan|precursor|psm|t-?test|anova|regulated|annotation|symbol)\b", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        /// <summary>Proteins, peptides, statistics: rows are not samples.</summary>
        internal static bool IsResultTable(SupplementTable t)
        {
            var cells = t.Rows.Take(200).SelectMany(r => r).Where(c => c.Length > 0).ToList();
            if (cells.Count > 0 && cells.Count(c => UniProt.IsMatch(c.Split(';')[0])) > 0.05 * cells.Count) return true;
            int result = t.Header.Sum(h => ResultWord.Matches(h).Count > 0 ? 1 : 0);
            int design = t.Header.Count(h => HeaderMap.ColumnFor(h) != null);
            return result >= Math.Max(2, design + 1);
        }

        private static bool FigureData(SupplementTable t) => Regex.IsMatch($"{t.Sheet} {t.Title}", @"\bfig", RegexOptions.IgnoreCase);

        /// <summary>Sample identifiers are words, not small numbers (a figure's "sample 1" column is a measurement key).</summary>
        private static bool IdentifierColumn(SupplementTable t, int j)
        {
            var v = Column(t, j);
            return v.Count > 0 && v.Count(x => Regex.IsMatch(x.Value, @"^[-+.\deE]+$")) <= 0.5 * v.Count;
        }

        private static int SampleColumn(SupplementTable t, int skipA, int skipB) =>
            Enumerable.Range(0, t.Header.Count).FirstOrDefault(j => j != skipA && j != skipB
                && HeaderMap.ColumnFor(t.Header[j]) is "source name" or "characteristics[individual]", -1);

        // ---------------- helpers ----------------

        private static SdrfEvidence Claim(string file, string label, string column, string value, string locator, string method, SdrfEvidenceConfidence confidence) =>
            new(file, label, column, value, Supplement, locator, method, confidence);

        private static List<(int Row, string Value)> Column(SupplementTable t, int j) =>
            t.Rows.Select((r, k) => (k, r[j].Trim())).Where(x => x.Item2.Length > 0).ToList();

        private static int Find(IReadOnlyList<string> header, string name)
        {
            for (int j = 0; j < header.Count; j++)
                if (header[j].Trim().Equals(name, StringComparison.OrdinalIgnoreCase)) return j;
            return -1;
        }

        private static bool Usable(string v)
        {
            string s = v.Trim();
            return s.Length > 0 && !s.Equals(SdrfReserved.NotAvailable, StringComparison.OrdinalIgnoreCase)
                && !s.Equals(SdrfReserved.NotApplicable, StringComparison.OrdinalIgnoreCase);
        }

        private static string FileStem(string name) => SdrfFileNamePattern.Stem(name.Trim().Replace('\\', '/').Split('/')[^1]);

        /// <summary>Raw files by exact stem and by name token, so a table value finds its files without a scan.</summary>
        private sealed class FileIndex
        {
            private readonly Dictionary<string, List<string>> _byStem = new(StringComparer.OrdinalIgnoreCase);
            private readonly Dictionary<string, HashSet<string>> _byToken = new(StringComparer.OrdinalIgnoreCase);

            internal FileIndex(IEnumerable<string> files)
            {
                foreach (var f in files.Distinct(StringComparer.OrdinalIgnoreCase))
                {
                    string stem = FileStem(f);
                    (_byStem.TryGetValue(stem, out var l) ? l : _byStem[stem] = new List<string>()).Add(f);
                    foreach (var t in Tokens(stem))
                        (_byToken.TryGetValue(t, out var s) ? s : _byToken[t] = new HashSet<string>(StringComparer.Ordinal)).Add(f);
                }
            }

            internal IReadOnlyList<string> ByExactStem(string stem) =>
                stem.Length > 0 && _byStem.TryGetValue(stem, out var l) ? l : Array.Empty<string>();

            /// <summary>Files whose stem is the value, or whose name tokens include every token of it (not a lone 1- or 2-digit number).</summary>
            internal IReadOnlyCollection<string> Match(string value)
            {
                string v = value.Trim();
                var hit = new HashSet<string>(ByExactStem(v), StringComparer.Ordinal);
                if (v.Length >= 2 && !Regex.IsMatch(v, @"^\d{1,2}$"))
                {
                    var tokens = Tokens(v).ToList();
                    if (tokens.Count > 0 && tokens.All(_byToken.ContainsKey))
                        hit.UnionWith(tokens.Skip(1).Aggregate(new HashSet<string>(_byToken[tokens[0]]), (a, t) => { a.IntersectWith(_byToken[t]); return a; }));
                }
                return hit;
            }

            private static IEnumerable<string> Tokens(string s) =>
                Regex.Split(s.ToLowerInvariant(), @"[_\-. ()\[\]]+").Where(t => t.Length > 0);
        }
    }

    /// <summary>
    /// Supplement headers mapped to SDRF columns, and values decoded by their header. Order matters: the first rule
    /// that matches wins, most specific first (measured against curated SDRFs, sdrf results/benchmark 2026-09-26).
    /// </summary>
    internal static class HeaderMap
    {
        private static readonly (string Column, Regex Header)[] Rules = new (string, string)[]
        {
            ("comment[data file]", @"raw ?file|file ?name|data ?file|ms ?run|raw ?data|^run\b|run ?(name|id)"),
            ("comment[label]", @"\btmt\b|itraq|channel|reporter|\blabel\b|plex"),
            ("comment[fraction identifier]", @"\bfraction"),
            ("characteristics[organism part]", @"organism ?part"),
            ("characteristics[developmental stage]", @"age ?group|developmental|life ?stage"),
            ("characteristics[age]", @"\bage\b|age ?\(|age at|^age"),
            ("characteristics[sex]", @"\bsex\b|gender"),
            ("characteristics[cell line]", @"cell ?line(?! ?group)(?! ?class)"),
            ("characteristics[cell type]", @"cell ?type"),
            ("characteristics[phenotype]", @"subtype|phenotype|classification|cluster"),
            ("characteristics[disease]", @"disease|diagnos|cancer ?type|tumou?r ?type|histolog|pathology"),
            ("characteristics[organism part]", @"tissue(?! ?enrich)(?! ?specific)|organ\b|anatom|body ?site|specimen|biopsy ?site"),
            ("characteristics[individual]", @"patient|donor|subject|individual|case ?(id|no|number|#)|\bcase\b|mouse ?(id|no)|animal|participant"),
            ("characteristics[strain]", @"strain|genotype|mutant|mutation|knock"),
            ("characteristics[treatment]", @"treatment|drug|compound|inhibitor|stimul|agent"),
            ("characteristics[dose]", @"dose|concentration"),
            ("characteristics[time]", @"time ?point|timepoint|\btime\b|\bday\b|\bweek\b|\bhours?\b"),
            ("characteristics[biological replicate]", @"biological ?rep|bio ?rep|replicate"),
            ("comment[batch]", @"\bbatch"),
            ("characteristics[organism]", @"organism|species"),
            ("source name", @"^sample|sample ?(id|name|no|number|#)|^id$|^name$|sample$"),
        }.Select(r => (r.Item1, new Regex(r.Item2, RegexOptions.IgnoreCase | RegexOptions.Compiled))).ToArray();

        /// <summary>The SDRF column a header describes, or null.</summary>
        internal static string? ColumnFor(string header)
        {
            string h = header.Trim();
            if (h.Length == 0) return null;
            foreach (var (column, rx) in Rules) if (rx.IsMatch(h)) return column;
            return null;
        }

        private static readonly HashSet<string> Missing = new(StringComparer.OrdinalIgnoreCase)
            { "na", "n/a", "nd", "-", "--", "none", "unknown", "not available", "not applicable", "?" };

        /// <summary>The value as an SDRF writes it, or null when it cannot be decided without a guess.</summary>
        internal static string? Decode(string column, string header, string value)
        {
            string v = value.Trim();
            if (v.Length == 0 || Missing.Contains(v)) return null;
            switch (column)
            {
                case "characteristics[sex]":
                    if (v is "0" or "1" && SexCode(header) is { } code) v = v == "1" ? code.One : code.Zero;
                    return v.ToLowerInvariant() switch
                    {
                        "m" or "male" or "man" or "men" => "male",
                        "f" or "female" or "woman" or "women" => "female",
                        _ => null
                    };
                case "characteristics[age]":
                    return Age(header, v);
                case "characteristics[biological replicate]":
                    return Regex.IsMatch(v, @"^[1-9]\d{0,2}$") ? v : null;
                case "characteristics[organism]":
                    if (Regex.IsMatch(v, @"^(not|no|none|unknown|n/?a)\b", RegexOptions.IgnoreCase)) return null;
                    return Regex.IsMatch(v, @"^[A-Z][a-z]+ [a-z]+( [a-z.]+)*$") || Regex.IsMatch(v, @"^[A-Z]\. ?[a-z]+$") ? v : null;
                default:
                    return v;
            }
        }

        /// <summary><c>Gender (1F,0M)</c> / <c>Sex (0=male, 1=female)</c>: what 1 and 0 stand for.</summary>
        private static (string One, string Zero)? SexCode(string header)
        {
            var m = Regex.Match(header, @"([01])\s*=?\s*([fm])\w*\s*[,;/ ]+\s*([01])\s*=?\s*([fm])", RegexOptions.IgnoreCase);
            if (!m.Success || m.Groups[1].Value == m.Groups[3].Value) return null;
            string first = m.Groups[2].Value.ToLowerInvariant() == "f" ? "female" : "male";
            string second = m.Groups[4].Value.ToLowerInvariant() == "f" ? "female" : "male";
            return m.Groups[1].Value == "1" ? (first, second) : (second, first);
        }

        /// <summary><c>77</c> -> <c>77Y</c>; the unit from the value, else the header (months, weeks, days), else years.</summary>
        private static string? Age(string header, string v)
        {
            var m = Regex.Match(v, @"^(\d+(?:\.\d+)?)\s*(y|yr|yrs|years?|years old|m|mo|months?|w|wk|weeks?|d|days?)?$", RegexOptions.IgnoreCase);
            if (!m.Success) return null;
            string unit = m.Groups[2].Success && m.Groups[2].Value.Length > 0 ? m.Groups[2].Value[..1].ToUpperInvariant()
                : Regex.IsMatch(header, @"month", RegexOptions.IgnoreCase) ? "M"
                : Regex.IsMatch(header, @"week", RegexOptions.IgnoreCase) ? "W"
                : Regex.IsMatch(header, @"\bdays?\b", RegexOptions.IgnoreCase) ? "D" : "Y";
            string n = m.Groups[1].Value;
            if (n.Contains('.')) n = n.TrimEnd('0').TrimEnd('.');
            return double.Parse(n, CultureInfo.InvariantCulture) > 0 ? n + unit : null;
        }
    }
}
