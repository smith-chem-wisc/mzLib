using System.Globalization;
using System.Text.RegularExpressions;

namespace Readers
{
    /// <summary>What a replicate marker counts.</summary>
    internal enum SdrfReplicateKind
    {
        /// <summary>The record does not settle it: no stated number matches the marker's count, or two do.</summary>
        Unstated,

        /// <summary>Independent biological samples (cultures, animals, patients).</summary>
        Biological,

        /// <summary>Re-injections or re-runs of one sample.</summary>
        Technical,

        /// <summary>Fractions of one sample.</summary>
        Fraction
    }

    /// <summary>
    /// One file's replicate marker, if it has one.
    /// </summary>
    /// <param name="FileName">The file name as given.</param>
    /// <param name="Base">The name without its marker(s): files sharing it are replicates of one another.</param>
    /// <param name="Number">The inner marker's rank within its base (1-based), or null when the file has none.</param>
    /// <param name="Outer">The outer marker's rank, when the name carries two levels (<c>X_2_1</c>), else null.</param>
    internal sealed record SdrfReplicateReading(string FileName, string Base, int? Number, int? Outer);

    /// <summary>
    /// What <see cref="SdrfReplicateResolver.Read"/> found: every file's markers, what the inner and outer
    /// markers count, and whether the record says there is only one biological sample.
    /// </summary>
    internal sealed record SdrfReplicates(
        IReadOnlyList<SdrfReplicateReading> Files,
        SdrfReplicateKind MarkerKind,
        string MarkerEvidence,
        SdrfReplicateKind OuterKind,
        bool SingleBiologicalSample,
        string SingleSampleEvidence);

    /// <summary>
    /// Reads a deposit's REPLICATES the way a curator does: siblings whose names differ only in a final
    /// marker (<c>CT10A/B/C</c>, <c>_R1.._R3</c>, <c>_Rat1..3</c>, <c>_repeat1..3</c>, <c>-1..-3</c>) are
    /// replicates of one base, numbered within that base; and the record's own words decide what they
    /// count.
    ///
    /// <para><b>Why.</b> In the blind benchmark (sdrf project, fresh sets 1 and 2, 2026-09-23) the drafted
    /// biological replicate was judged right 25-29% of the time against the curated 55-67%. The graders'
    /// quotes gave five causes: numbering across the whole deposit when no condition was found; re-injections,
    /// fractions or a pooled sample counted as biological replicates; a replicate marker in the name ignored;
    /// individuals not ranked within their condition; and one stated sample numbered anyway. The first
    /// three are this type's job.</para>
    ///
    /// <para><b>Local, not global.</b> A marker is read per BASE NAME, so a deposit that mixes name shapes
    /// still numbers each base 1..n, and numbering can never run across conditions. A marker needs two or
    /// more siblings whose values run contiguously; a four-digit animal ID is not a marker.</para>
    ///
    /// <para><b>The record decides the kind, only on a matching number.</b> "three biological replicates",
    /// "grown in triplicate", "n = 3" make a 1..3 marker biological; "analyzed in triplicate", "injected three
    /// times" make it technical; "3 fractions" makes it a fraction. A stated number that does not match the
    /// marker's count decides nothing, and neither do two kinds stating the same number. With two marker
    /// levels (<c>X_2_1</c>), the outer defaults to biological and the inner to technical.</para>
    ///
    /// <para>Pure: names and text in, a reading out.</para>
    /// </summary>
    internal static class SdrfReplicateResolver
    {
        private const string Num = @"(\d+|one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|duplicates?|triplicates?|quadruplicates?|twice)";

        private static readonly Regex NumberMarker = new(
            @"^(?<base>.*?)[_\-\.\s]*(?:biorep|techrep|replicate|repeat|rep|rat|mouse|br|tr|r)?(?<![0-9])(?<n>\d{1,2})$",
            RegexOptions.IgnoreCase | RegexOptions.Compiled);
        private static readonly Regex LetterMarker = new(
            @"^(?<base>.*?[0-9A-Za-z])[_\-\.\s]?(?<l>[A-Da-d])$", RegexOptions.Compiled);

        private static readonly Regex[] Biological =
        {
            new(@"\b" + Num + @"\s+(?:independent\s+)?biological(?:ly)?(?:\s+independent)?\s+(?:replicates?|samples?|triplicates?|duplicates?|repeats?)", RegexOptions.IgnoreCase | RegexOptions.Compiled),
            new(@"\bbiological\s+" + Num, RegexOptions.IgnoreCase | RegexOptions.Compiled),
            new(@"\b(?:grown|cultured|cultivated|incubated|prepared|harvested|collected|treated|infected)\s+in\s+" + Num, RegexOptions.IgnoreCase | RegexOptions.Compiled),
            new(@"\bn\s*=\s*(\d+)", RegexOptions.IgnoreCase | RegexOptions.Compiled),
            new(@"\b" + Num + @"\s+(?:\S+\s+)?(?:mice|animals|patients|donors|rats|individuals|subjects|men|women|volunteers|plants|cultures)\b", RegexOptions.IgnoreCase | RegexOptions.Compiled),
            new(@"\b" + Num + @"\s+replicates\b", RegexOptions.IgnoreCase | RegexOptions.Compiled),
        };

        private static readonly Regex[] Technical =
        {
            new(@"\b" + Num + @"\s+technical\s+(?:replicates?|triplicates?|duplicates?|repeats?|injections?)", RegexOptions.IgnoreCase | RegexOptions.Compiled),
            new(@"\btechnical\s+" + Num, RegexOptions.IgnoreCase | RegexOptions.Compiled),
            new(@"\b(?:analy[sz]ed|injected|measured|run|acquired|performed)\s+(?:\w+\s+){0,3}?in\s+" + Num, RegexOptions.IgnoreCase | RegexOptions.Compiled),
            new(@"\b(?:injected|analy[sz]ed|measured|run|acquired)\s+(?:" + Num + @"\s+times|(twice))", RegexOptions.IgnoreCase | RegexOptions.Compiled),
            new(@"\b" + Num + @"\s+(?:repeated\s+|replicate\s+)?injections\b", RegexOptions.IgnoreCase | RegexOptions.Compiled),
        };

        private static readonly Regex[] Fractions =
        {
            new(@"\b" + Num + @"\s+(?:\S+\s+){0,2}?fractions\b", RegexOptions.IgnoreCase | RegexOptions.Compiled),
        };

        private static readonly Regex[] OneSample =
        {
            new(@"\bpooled\b", RegexOptions.IgnoreCase | RegexOptions.Compiled),
            new(@"\ba\s+single\s+(?:\S+\s+)?(?:lysate|culture|sample|cell\s+line|biosample|extract)", RegexOptions.IgnoreCase | RegexOptions.Compiled),
            new(@"\bone\s+(?:lysate|sample|culture|biosample|extract)\b", RegexOptions.IgnoreCase | RegexOptions.Compiled),
        };

        /// <summary>
        /// Reads one deposit's file names against its record text. Throws only on a null argument.
        /// </summary>
        public static SdrfReplicates Read(IEnumerable<string> fileNames, IEnumerable<string?> recordText)
        {
            if (fileNames == null) throw new ArgumentNullException(nameof(fileNames));
            if (recordText == null) throw new ArgumentNullException(nameof(recordText));
            var names = fileNames.ToList();
            string text = string.Join(" \n ", recordText.Where(t => !string.IsNullOrWhiteSpace(t)));

            // Level 1: the last marker of each name.
            var inner = Markers(names.ToDictionary(n => n, n => SdrfFileNamePattern.Stem(n), StringComparer.Ordinal));
            // Level 2: a marker on the BASES that level 1 left (X_2_1 -> base X_2 -> base X, outer 2).
            var bases = inner.Values.Where(m => m.Rank != null).Select(m => m.Base).Distinct(StringComparer.OrdinalIgnoreCase)
                .ToDictionary(b => b, b => b, StringComparer.OrdinalIgnoreCase);
            var outer = Markers(bases);

            var files = names.Select(n =>
            {
                var m = inner[n];
                if (m.Rank == null) return new SdrfReplicateReading(n, SdrfFileNamePattern.Stem(n), null, null);
                var o = outer.TryGetValue(m.Base, out var om) && om.Rank != null ? om : null;
                return new SdrfReplicateReading(n, o?.Base ?? m.Base, m.Rank, o?.Rank);
            }).ToList();

            int innerCount = files.Where(f => f.Number != null).Select(f => f.Number!.Value).DefaultIfEmpty(0).Max();
            int outerCount = files.Where(f => f.Outer != null).Select(f => f.Outer!.Value).DefaultIfEmpty(0).Max();
            bool twoLevels = outerCount > 0;

            var bio = Counts(text, Biological);
            var tech = Counts(text, Technical);
            var frac = Counts(text, Fractions);

            var (kind, evidence) = innerCount == 0
                ? (SdrfReplicateKind.Unstated, "no file name carries a replicate marker")
                : Decide(innerCount, bio, tech, frac, twoLevels ? SdrfReplicateKind.Technical : SdrfReplicateKind.Unstated);
            var outerKind = twoLevels ? Decide(outerCount, bio, tech, frac, SdrfReplicateKind.Biological).Kind : SdrfReplicateKind.Unstated;

            bool statesBiological = bio.Count > 0 || Regex.IsMatch(text, @"\bbiological(?:ly)?\b", RegexOptions.IgnoreCase);
            string? single = OneSample.Select(r => r.Match(text)).FirstOrDefault(m => m.Success)?.Value;
            if (single == null && tech.Count > 0 && !statesBiological) single = "technical replicates only";
            bool one = single != null && !statesBiological;

            return new SdrfReplicates(files, kind, evidence, outerKind, one,
                one ? $"the record says '{single}' and states no biological replicates" : "");
        }

        private sealed record Marker(string Base, int? Rank);

        /// <summary>
        /// Stem -> its final marker, ranked within its base. A base qualifies only with two or more members
        /// whose marker values are distinct and run contiguously.
        /// </summary>
        private static Dictionary<string, Marker> Markers(IReadOnlyDictionary<string, string> stems)
        {
            var parsed = stems.ToDictionary(kv => kv.Key, kv => Parse(kv.Value), StringComparer.Ordinal);
            var result = parsed.ToDictionary(kv => kv.Key, kv => new Marker(kv.Value.Base ?? stems[kv.Key], null), StringComparer.Ordinal);
            foreach (var group in parsed.Where(kv => kv.Value.Base != null)
                         .GroupBy(kv => kv.Value.Base!.ToLowerInvariant() + "|" + kv.Value.Letter))
            {
                var members = group.OrderBy(kv => kv.Value.Value).ToList();
                if (members.Count < 2) continue;
                var values = members.Select(kv => kv.Value.Value).ToList();
                if (values.Distinct().Count() != values.Count) continue;
                if (values.Where((v, i) => v != values[0] + i).Any()) continue;
                for (int i = 0; i < members.Count; i++)
                    result[members[i].Key] = new Marker(members[i].Value.Base!, i + 1);
            }
            return result;
        }

        private static (string? Base, int Value, bool Letter) Parse(string stem)
        {
            var m = NumberMarker.Match(stem);
            if (m.Success && m.Groups["base"].Value.Length > 0)
                return (m.Groups["base"].Value.TrimEnd('_', '-', '.', ' '), int.Parse(m.Groups["n"].Value, CultureInfo.InvariantCulture), false);
            var l = LetterMarker.Match(stem);
            if (l.Success)
                return (l.Groups["base"].Value.TrimEnd('_', '-', '.', ' '), char.ToUpperInvariant(l.Groups["l"].Value[0]) - 'A' + 1, true);
            return (null, 0, false);
        }

        private static (SdrfReplicateKind Kind, string Evidence) Decide(int count, HashSet<int> bio, HashSet<int> tech, HashSet<int> frac, SdrfReplicateKind fallback)
        {
            var matches = new List<SdrfReplicateKind>();
            if (bio.Contains(count)) matches.Add(SdrfReplicateKind.Biological);
            if (tech.Contains(count)) matches.Add(SdrfReplicateKind.Technical);
            if (frac.Contains(count)) matches.Add(SdrfReplicateKind.Fraction);
            if (matches.Count == 1)
                return (matches[0], $"the record states {count} {matches[0].ToString().ToLowerInvariant()} replicate(s)/fraction(s), matching the 1..{count} marker");
            return (fallback, matches.Count > 1
                ? $"the record states {count} of more than one kind, so the 1..{count} marker's kind is not decided"
                : $"no number in the record matches the 1..{count} marker");
        }

        private static HashSet<int> Counts(string text, Regex[] patterns)
        {
            var found = new HashSet<int>();
            foreach (var r in patterns)
                foreach (Match m in r.Matches(text))
                    foreach (Group g in m.Groups.Cast<Group>().Skip(1))
                        if (g.Success && ToNumber(g.Value) is int n) found.Add(n);
            return found;
        }

        private static int? ToNumber(string word)
        {
            string w = word.ToLowerInvariant().TrimEnd('s');
            if (int.TryParse(w, NumberStyles.None, CultureInfo.InvariantCulture, out int n)) return n;
            return w switch
            {
                "one" => 1, "two" => 2, "twice" => 2, "duplicate" => 2, "three" => 3, "triplicate" => 3,
                "four" => 4, "quadruplicate" => 4, "five" => 5, "six" => 6, "seven" => 7, "eight" => 8,
                "nine" => 9, "ten" => 10, "eleven" => 11, "twelve" => 12, _ => null
            };
        }
    }
}
