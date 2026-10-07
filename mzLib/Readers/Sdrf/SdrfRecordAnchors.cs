using System.Text.RegularExpressions;

namespace Readers
{
    /// <summary>
    /// One condition found as file-name tokens the deposit's own record also uses.
    /// </summary>
    /// <param name="Levels">The values, as the file names write them, commonest first.</param>
    /// <param name="Evidence">Why these tokens were read as one condition, in words a curator can check.</param>
    internal sealed record SdrfAnchoredFactor(IReadOnlyList<string> Levels, string Evidence);

    /// <summary>
    /// What <see cref="SdrfRecordAnchors.Read"/> found. <see cref="LevelsByFile"/> gives, for every file,
    /// its level of each factor in <see cref="Factors"/> order, or an empty string when the file carries
    /// none of that factor's levels. Never null; empty lists when nothing was anchored.
    /// </summary>
    internal sealed record SdrfAnchoredFactors(
        IReadOnlyList<SdrfAnchoredFactor> Factors,
        IReadOnlyDictionary<string, IReadOnlyList<string>> LevelsByFile);

    /// <summary>
    /// Finds a deposit's CONDITIONS as the file-name tokens that its own PRIDE record also uses:
    /// <c>LPS</c>, <c>TNF</c>, <c>IL1BETA</c> in the names and "exposed to TNF, IL-1β, and LPS" in the
    /// description.
    ///
    /// <para><b>Why a second reader.</b> <see cref="SdrfFileNamePattern"/> reads structure from how names
    /// vary, one family of like-shaped names at a time. The blind benchmark's fresh set (sdrf project,
    /// 2026-09-23) found that it lost on factors in two ways this fixes. Conditions spanned several name
    /// shapes (<c>IL1BETA_R1</c> beside <c>PHA_CTRL_R1</c>), and no single family saw them all. And it
    /// read animal IDs (<c>0316</c>) as categories. A token the depositor ALSO wrote in the record is a
    /// condition with independent evidence, whatever shape its name has. A token the record never
    /// mentions is not anchored, however it varies.</para>
    ///
    /// <para><b>Deliberately narrow.</b> An anchor is an exact word of the record, or two or three of its
    /// words joined (<c>IL-1β</c> -> <c>il1beta</c>), of three letters or more, and neither an English
    /// function word nor an instrument or acquisition word. A short list of control words (<c>ctrl</c>,
    /// <c>wt</c>, <c>sham</c>, <c>mock</c>...) needs no anchor. Tokens that never share a file form one
    /// factor. A factor needs two or more levels, each on two or more files, covering at least half the
    /// deposit. Abbreviations the record spells out (<c>as_</c> for asthenozoospermic) are NOT guessed.</para>
    ///
    /// <para>Pure: names and text in, a reading out. No network, no PRIDE types.</para>
    /// </summary>
    internal static class SdrfRecordAnchors
    {
        private static readonly Regex Parts = new(@"[_\-\.\s]+", RegexOptions.Compiled);
        private static readonly Regex Letters = new(@"[A-Za-z]{2,}", RegexOptions.Compiled);
        private static readonly Regex Words = new(@"[a-z0-9]+", RegexOptions.Compiled);

        // Conditions depositors abbreviate without spelling out; no anchor needed.
        private static readonly HashSet<string> ControlWords = new(StringComparer.OrdinalIgnoreCase)
        {
            "ctrl", "ctl", "control", "mock", "sham", "untreated", "vehicle", "wt", "wildtype", "ko", "neg", "pos"
        };

        // Never a condition, however often the record uses it: function words, and the words of how a
        // file was acquired rather than what it contains.
        private static readonly HashSet<string> NeverAnchors = new(StringComparer.OrdinalIgnoreCase)
        {
            "and", "the", "for", "with", "from", "were", "was", "are", "not", "all", "each", "per", "into", "using",
            "then", "than", "that", "this", "these", "those", "which", "after", "before", "between", "both",
            "run", "runs", "rep", "replicate", "replicates", "sample", "samples", "data", "file", "raw", "mode",
            "method", "blank", "std", "standard", "test", "mix", "fraction", "fractions", "band", "frac", "inj",
            "injection", "dia", "dda", "prm", "srm", "mrm", "swath", "tmt", "itraq", "silac", "lfq", "msms",
            "orbitrap", "lumos", "exploris", "fusion", "eclipse", "velos", "elite", "timstof", "qexactive",
            "exactive", "astral", "ltq", "tof", "qtof", "nano", "hplc", "uplc", "lcms", "ms", "lc", "qc"
        };

        /// <summary>
        /// Reads one deposit's file names against its record text (title, description, protocols,
        /// keywords -- any strings; nulls are skipped). Throws only on a null argument.
        /// </summary>
        public static SdrfAnchoredFactors Read(IEnumerable<string> fileNames, IEnumerable<string?> recordText)
        {
            if (fileNames == null) throw new ArgumentNullException(nameof(fileNames));
            if (recordText == null) throw new ArgumentNullException(nameof(recordText));
            var names = fileNames.ToList();
            var vocabulary = Vocabulary(recordText);

            // file -> anchored tokens (lower case), and the first spelling of each as the names write it
            var written = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            var tokensOf = names.ToDictionary(n => n, n => AnchoredTokens(n, vocabulary, written), StringComparer.Ordinal);

            int n = names.Count;
            var filesOf = tokensOf.SelectMany(kv => kv.Value.Select(t => (t, file: kv.Key)))
                .GroupBy(x => x.t, StringComparer.OrdinalIgnoreCase)
                .ToDictionary(g => g.Key, g => g.Select(x => x.file).ToHashSet(StringComparer.Ordinal), StringComparer.OrdinalIgnoreCase);
            var levels = filesOf.Where(kv => kv.Value.Count >= 2 && kv.Value.Count < n)
                .OrderByDescending(kv => kv.Value.Count).ThenBy(kv => kv.Key, StringComparer.Ordinal)
                .Select(kv => kv.Key).ToList();

            // Tokens that never share a file are levels of one factor.
            var groups = new List<List<string>>();
            foreach (var level in levels)
            {
                var home = groups.FirstOrDefault(g => g.All(other => !filesOf[other].Overlaps(filesOf[level])));
                if (home != null) home.Add(level);
                else groups.Add(new List<string> { level });
            }

            var factors = new List<(List<string> Levels, int Covered)>();
            foreach (var g in groups)
            {
                int covered = g.Sum(l => filesOf[l].Count);
                if (g.Count >= 2 && covered * 2 >= n) factors.Add((g, covered));
            }
            factors.Sort((a, b) => b.Covered.CompareTo(a.Covered));

            var result = factors.Select(f => new SdrfAnchoredFactor(
                f.Levels.Select(l => written[l]).ToList(),
                $"'{string.Join("', '", f.Levels.Select(l => written[l]))}' each appear in the PRIDE record text and " +
                $"never together in one file name; together they cover {f.Covered} of {n} files")).ToList();
            var byFile = names.ToDictionary(
                name => name,
                name => (IReadOnlyList<string>)factors.Select(f =>
                    f.Levels.Where(l => filesOf[l].Contains(name)).Select(l => written[l]).FirstOrDefault() ?? "").ToList(),
                StringComparer.Ordinal);
            return new SdrfAnchoredFactors(result, byFile);
        }

        /// <summary>
        /// Every word of the record, lower-cased, plus every two and three consecutive words joined, so a
        /// name that writes <c>IL1BETA</c> finds a record that writes <c>IL-1β</c>.
        /// </summary>
        private static HashSet<string> Vocabulary(IEnumerable<string?> text)
        {
            var v = new HashSet<string>(StringComparer.Ordinal);
            foreach (var t in text)
            {
                if (string.IsNullOrWhiteSpace(t)) continue;
                var words = Words.Matches(Greek(t.ToLowerInvariant())).Select(m => m.Value).ToList();
                for (int i = 0; i < words.Count; i++)
                {
                    v.Add(words[i]);
                    if (i + 1 < words.Count) v.Add(words[i] + words[i + 1]);
                    if (i + 2 < words.Count) v.Add(words[i] + words[i + 1] + words[i + 2]);
                }
            }
            return v;
        }

        private static string Greek(string s) => s
            .Replace("α", "alpha").Replace("β", "beta").Replace("γ", "gamma").Replace("δ", "delta")
            .Replace("κ", "kappa").Replace("μ", "u");

        /// <summary>
        /// The anchored tokens of one file name: a whole part (<c>IL1BETA</c>) when the record uses it, else
        /// each run of letters inside the part (<c>PHA</c> in <c>PHACTRL</c>) that it uses.
        /// </summary>
        private static HashSet<string> AnchoredTokens(string fileName, HashSet<string> vocabulary, Dictionary<string, string> written)
        {
            var found = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
            foreach (var part in Parts.Split(SdrfFileNamePattern.Stem(fileName)).Where(p => p.Length > 0))
            {
                if (IsAnchor(part, vocabulary)) { Add(part); continue; }
                foreach (Match m in Letters.Matches(part))
                    if (m.Value.Length < part.Length && IsAnchor(m.Value, vocabulary)) Add(m.Value);
            }
            return found;

            void Add(string token)
            {
                found.Add(token);
                written.TryAdd(token, token);
            }
        }

        private static bool IsAnchor(string token, HashSet<string> vocabulary)
        {
            string t = token.ToLowerInvariant();
            if (!t.Any(char.IsLetter) || NeverAnchors.Contains(t)) return false;
            if (ControlWords.Contains(t)) return true;
            return t.Length >= 3 && vocabulary.Contains(t);
        }
    }
}
