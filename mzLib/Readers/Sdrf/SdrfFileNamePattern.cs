using System.Globalization;
using System.Text.RegularExpressions;

namespace Readers
{
    /// <summary>
    /// What a varying part of a deposit's file names was read as.
    /// </summary>
    internal enum SdrfFileNameRole
    {
        /// <summary>
        /// An index named as a fraction (<c>Band_01..10</c>). Kept as written, not ranked per sample,
        /// so fractions still line up across samples and a missing file stays a visible gap.
        /// </summary>
        Fraction,

        /// <summary>An index named as a biological replicate (<c>BR1</c>, <c>bio2</c>).</summary>
        BiologicalReplicate,

        /// <summary>An index named as a technical replicate or re-injection (<c>TR1</c>, <c>inj2</c>).</summary>
        TechnicalReplicate,

        /// <summary>
        /// A replicate index whose kind the name does not state (<c>WT_DMSO1..3</c>). Biological and
        /// technical cannot be told apart from a trailing number, so the reading says so rather than
        /// choosing; the caller decides, and marks the choice as inferred.
        /// </summary>
        Replicate,

        /// <summary>An index that names the sample itself (<c>Lysate1</c>, <c>Patient2</c>).</summary>
        Sample,

        /// <summary>A word that takes two or more values, each shared by two or more files (<c>WT</c>/<c>KO</c>).</summary>
        Factor,

        /// <summary>A date that splits the files into groups of two or more (<c>20210304_...</c>).</summary>
        Batch
    }

    /// <summary>
    /// One varying part of the file names, the role it was read as, and why.
    /// </summary>
    /// <param name="Position">The token position in the aligned names, 0-based.</param>
    /// <param name="Role">What the part was read as.</param>
    /// <param name="Levels">The distinct values it takes, as written, in sorted order.</param>
    /// <param name="Evidence">The pattern that produced the role, in words a curator can check.</param>
    internal sealed record SdrfFileNameSlot(int Position, SdrfFileNameRole Role, IReadOnlyList<string> Levels, string Evidence);

    /// <summary>
    /// One file's reading. Numbers are 1-based, the convention <see cref="SdrfBuilder"/> takes. A replicate
    /// or sample index is renumbered within its group when the names count from elsewhere; a fraction
    /// never is, and is only shifted, for the whole deposit at once, when its count starts at 0. The
    /// slot's evidence says which happened.
    /// </summary>
    /// <param name="FileName">The file name as given.</param>
    /// <param name="SampleKey">
    /// Equal for files read as the same sample: everything that varies except the fraction, the
    /// technical replicate and the batch. Empty on every file when nothing but a fraction or a
    /// technical replicate varies: the whole deposit is then one sample.
    /// </param>
    /// <param name="FactorLevels">This file's value for each <see cref="SdrfFileNameRole.Factor"/> slot, in slot order.</param>
    internal sealed record SdrfFileNameReading(
        string FileName,
        string SampleKey,
        IReadOnlyList<string> FactorLevels,
        int? Fraction,
        int? BiologicalReplicate,
        int? TechnicalReplicate,
        int? Replicate,
        string? Batch,
        bool IsControl);

    /// <summary>
    /// What <see cref="SdrfFileNamePattern.Read"/> found. When <see cref="Found"/> is false,
    /// <see cref="NoStructureReason"/> says why, <see cref="Slots"/> is empty, and every reading
    /// carries only its file name as its sample key.
    /// </summary>
    internal sealed record SdrfFileNameStructure(
        IReadOnlyList<SdrfFileNameSlot> Slots,
        IReadOnlyList<SdrfFileNameReading> Files,
        string? NoStructureReason)
    {
        public bool Found => NoStructureReason is null;
    }

    /// <summary>
    /// Reads a design out of ONE deposit's file names, or says that it cannot.
    ///
    /// <para><b>How far to trust it.</b> Inferring a design from file names has not been measured
    /// yet. The nearest measurement -- partitioning channel-level files by their set of
    /// <c>source name</c> values, right in 3 of 9 files with a plex column -- failed by OVER-SPLITTING,
    /// so the rules here lean against that failure. (1) A token earns a role only from how it varies
    /// across the SET, never from one name: an unnamed index must run contiguously inside every group,
    /// and a word is a factor only if every value it takes is shared by two or more files.
    /// (2) Anything that names each file differently is an identifier, and an identifier is not
    /// structure. (3) When the evidence does not decide, the answer is "no structure found", with the
    /// reason, and every cell a caller writes from a found structure is marked inferred with its
    /// <see cref="SdrfFileNameSlot.Evidence"/>.</para>
    ///
    /// <para>Pure: file names in, a reading out. The PRIDE listing, the provenance marking and the
    /// SDRF itself belong to the caller.</para>
    /// </summary>
    internal static class SdrfFileNamePattern
    {
        // Longest first, so ".mzml.gz" is stripped whole before ".gz" alone could be.
        private static readonly string[] DataExtensions =
        {
            ".mzml.gz", ".raw.gz", ".mzxml.gz", ".mgf.gz", ".raw.zip", ".d.zip",
            ".raw", ".mzml", ".mzxml", ".mgf", ".wiff2", ".wiff", ".d", ".baf", ".tdf", ".gz", ".zip"
        };

        private static readonly Regex Separators = new(@"[_\-\.\s]+", RegexOptions.Compiled);
        private static readonly Regex WordThenNumber = new(@"^(?<w>[A-Za-z]+)(?<n>\d+)$", RegexOptions.Compiled);

        // The word just before an index, when it states what the index counts. Compared ignoring case.
        private static readonly Dictionary<string, SdrfFileNameRole> IndexWords = new(StringComparer.OrdinalIgnoreCase)
        {
            ["band"] = SdrfFileNameRole.Fraction, ["fraction"] = SdrfFileNameRole.Fraction,
            ["frac"] = SdrfFileNameRole.Fraction, ["fr"] = SdrfFileNameRole.Fraction,
            ["f"] = SdrfFileNameRole.Fraction, ["fx"] = SdrfFileNameRole.Fraction,
            ["slice"] = SdrfFileNameRole.Fraction, ["gel"] = SdrfFileNameRole.Fraction,
            ["br"] = SdrfFileNameRole.BiologicalReplicate, ["bio"] = SdrfFileNameRole.BiologicalReplicate,
            ["biorep"] = SdrfFileNameRole.BiologicalReplicate,
            ["tr"] = SdrfFileNameRole.TechnicalReplicate, ["tech"] = SdrfFileNameRole.TechnicalReplicate,
            ["techrep"] = SdrfFileNameRole.TechnicalReplicate, ["inj"] = SdrfFileNameRole.TechnicalReplicate,
            ["injection"] = SdrfFileNameRole.TechnicalReplicate,
            ["rep"] = SdrfFileNameRole.Replicate, ["replicate"] = SdrfFileNameRole.Replicate,
            ["r"] = SdrfFileNameRole.Replicate,
            ["sample"] = SdrfFileNameRole.Sample, ["s"] = SdrfFileNameRole.Sample,
            ["lysate"] = SdrfFileNameRole.Sample, ["patient"] = SdrfFileNameRole.Sample,
            ["donor"] = SdrfFileNameRole.Sample, ["subject"] = SdrfFileNameRole.Sample,
            ["mouse"] = SdrfFileNameRole.Sample, ["animal"] = SdrfFileNameRole.Sample,
        };

        // A factor level that names a control arm. An IP's IgG must never merge into its bait.
        private static readonly Regex ControlLevel = new(
            @"^(igg\w*|isotype|control|ctrl|ctl|mock|blank|neg|negative)$",
            RegexOptions.IgnoreCase | RegexOptions.Compiled);

        private static readonly string[] DateFormats = { "yyyyMMdd", "yyMMdd" };

        /// <summary>
        /// Reads the file names of one deposit. Throws only on a malformed argument; weak evidence is
        /// described, never thrown.
        /// </summary>
        /// <exception cref="ArgumentNullException"><paramref name="fileNames"/> is null.</exception>
        /// <exception cref="ArgumentException">A name is blank, or two names are the same.</exception>
        public static SdrfFileNameStructure Read(IEnumerable<string> fileNames)
        {
            if (fileNames == null) throw new ArgumentNullException(nameof(fileNames));
            var names = fileNames.ToList();
            for (int i = 0; i < names.Count; i++)
                if (string.IsNullOrWhiteSpace(names[i]))
                    throw new ArgumentException($"File name {i} is blank.", nameof(fileNames));
            var duplicate = names.GroupBy(n => n, StringComparer.OrdinalIgnoreCase).FirstOrDefault(g => g.Count() > 1);
            if (duplicate != null)
                throw new ArgumentException(
                    $"'{duplicate.Key}' is listed {duplicate.Count()} times; pass each file of the deposit once.",
                    nameof(fileNames));

            if (names.Count < 2)
                return None(names, "one file has no siblings to compare it with");

            var stems = names.Select(Stem).ToList();
            var tokens = Align(stems);
            if (tokens == null)
                return None(names, "the file names do not share one shape, so no part of them can be compared across the set");

            int width = tokens[0].Length;
            var varying = Enumerable.Range(0, width)
                .Where(p => tokens.Select(t => t[p]).Distinct(StringComparer.OrdinalIgnoreCase).Count() > 1)
                .ToList();
            if (varying.Count == 0)
                return None(names, "the file names differ only in their extension");

            var slots = new List<SdrfFileNameSlot>();
            var batchPositions = new List<int>();
            var designPositions = new List<int>();
            foreach (int p in varying)
            {
                var values = tokens.Select(t => t[p]).ToList();
                if (values.All(IsDate))
                {
                    var byDate = values.GroupBy(v => v).ToList();
                    if (byDate.Count < names.Count && byDate.All(g => g.Count() >= 2))
                    {
                        batchPositions.Add(p);
                        slots.Add(new SdrfFileNameSlot(p, SdrfFileNameRole.Batch, Sorted(values),
                            $"a date at part {p + 1} splits the {names.Count} files into {byDate.Count} groups of two or more"));
                    }
                    // A date on every file is when it was acquired, not a design: ignored, not refused.
                    continue;
                }
                designPositions.Add(p);
            }
            if (designPositions.Count == 0)
                return None(names, "the only part that varies is an acquisition date");

            // Two files the design cannot tell apart would be one sample written twice.
            var keys = tokens.Select(t => string.Join("\u001f", designPositions.Select(p => t[p].ToUpperInvariant()))).ToList();
            if (keys.Distinct().Count() < keys.Count)
                return None(names, "two files read as the same design cell once the dates are set aside");

            var factorPositions = designPositions.Where(p => !IsNumber(tokens[0][p])).ToList();
            var indexPositions = designPositions.Where(p => IsNumber(tokens[0][p])).ToList();
            if (designPositions.Any(p => tokens.Any(t => IsNumber(t[p]) != IsNumber(tokens[0][p]))))
                return None(names, "a part is a number in some names and a word in others");

            foreach (int p in factorPositions)
            {
                var values = tokens.Select(t => t[p]).ToList();
                var levels = values.GroupBy(v => v, StringComparer.OrdinalIgnoreCase).ToList();
                if (levels.Any(g => g.Count() < 2))
                    return None(names, $"part {p + 1} ('{levels.First(g => g.Count() < 2).Key}') names a single file, " +
                                       "so it reads as an identifier, not as a condition shared by replicates");
                string controls = string.Join(", ", levels.Select(g => g.Key).Where(k => ControlLevel.IsMatch(k)));
                slots.Add(new SdrfFileNameSlot(p, SdrfFileNameRole.Factor, Sorted(values),
                    $"part {p + 1} takes {levels.Count} values, each shared by two or more files" +
                    (controls.Length > 0 ? $"; '{controls}' read as a control" : "")));
            }

            var numbers = new Dictionary<int, int[]>();
            int unnamed = 0;
            foreach (int p in indexPositions)
            {
                var role = IndexRole(tokens[0], p, factorPositions, out string named);
                if (role == null)
                {
                    unnamed++;
                    role = SdrfFileNameRole.Replicate;
                }
                if (unnamed > 1)
                    return None(names, "two numbered parts carry no word saying what they count, so neither can be read");

                var others = designPositions.Where(q => q != p).ToList();
                if (role == SdrfFileNameRole.Fraction)
                {
                    // Fractions are matched ACROSS samples by index (FlashLFQ transfers only between
                    // fractions at most one apart), so they are never ranked per sample: a lysate whose
                    // band 1 was never uploaded keeps band 2 as fraction 2, and a gap stays a gap for the
                    // design reader to report (MAP-34). One offset for the whole deposit, and only to
                    // move a count that starts at 0 onto the 1-based SDRF scale.
                    var bands = tokens.Select(t => long.Parse(t[p], CultureInfo.InvariantCulture)).ToArray();
                    long shift = bands.Min() == 0 ? 1 : 0;
                    numbers[p] = bands.Select(v => (int)(v + shift)).ToArray();
                    slots.Add(new SdrfFileNameSlot(p, role.Value, Sorted(tokens.Select(t => t[p])),
                        $"part {p + 1} numbers the files{named}" +
                        (shift > 0 ? "; shifted by 1 for the whole deposit, because the count starts at 0" : "")));
                    continue;
                }
                if (!TryRank(tokens, p, others, out var rank, out bool renumbered))
                    return None(names, $"the numbers at part {p + 1} do not run contiguously inside each group, " +
                                       "so they read as run numbers, not as fractions or replicates");
                if (others.Count == 0 && role != SdrfFileNameRole.Fraction)
                    return None(names, $"the only part that varies is a number (part {p + 1}), which cannot tell " +
                                       "a replicate from a run number");

                numbers[p] = rank;
                slots.Add(new SdrfFileNameSlot(p, role.Value, Sorted(tokens.Select(t => t[p])),
                    $"part {p + 1} counts {Distinct(rank)} within each group" + named +
                    (renumbered ? "; renumbered from 1 within each group" : "")));
            }

            slots.Sort((a, b) => a.Position.CompareTo(b.Position));
            var readings = new List<SdrfFileNameReading>(names.Count);
            for (int i = 0; i < names.Count; i++)
            {
                var t = tokens[i];
                int? Number(SdrfFileNameRole r) =>
                    slots.Where(s => s.Role == r && numbers.ContainsKey(s.Position)).Select(s => (int?)numbers[s.Position][i]).FirstOrDefault();
                var sampleParts = slots
                    .Where(s => s.Role is SdrfFileNameRole.Factor or SdrfFileNameRole.Sample
                        or SdrfFileNameRole.BiologicalReplicate or SdrfFileNameRole.Replicate)
                    .Select(s => s.Role == SdrfFileNameRole.Factor ? t[s.Position] : numbers[s.Position][i].ToString(CultureInfo.InvariantCulture));
                var factorLevels = slots.Where(s => s.Role == SdrfFileNameRole.Factor).Select(s => t[s.Position]).ToList();
                readings.Add(new SdrfFileNameReading(
                    names[i],
                    string.Join(" ", sampleParts),
                    factorLevels,
                    Number(SdrfFileNameRole.Fraction),
                    Number(SdrfFileNameRole.BiologicalReplicate),
                    Number(SdrfFileNameRole.TechnicalReplicate),
                    Number(SdrfFileNameRole.Replicate),
                    batchPositions.Count > 0 ? string.Join(" ", batchPositions.Select(p => t[p])) : null,
                    factorLevels.Any(l => ControlLevel.IsMatch(l))));
            }
            return new SdrfFileNameStructure(slots, readings, null);
        }

        private static SdrfFileNameStructure None(IReadOnlyList<string> names, string reason) =>
            new(Array.Empty<SdrfFileNameSlot>(),
                names.Select(n => new SdrfFileNameReading(n, n, Array.Empty<string>(), null, null, null, null, null, false)).ToList(),
                reason);

        /// <summary>The file name without its folder or its data-file extensions.</summary>
        internal static string Stem(string fileName)
        {
            string stem = Path.GetFileName(fileName.Trim().TrimEnd('/', '\\'));
            foreach (string ext in DataExtensions)
                if (stem.Length > ext.Length && stem.EndsWith(ext, StringComparison.OrdinalIgnoreCase))
                {
                    stem = stem[..^ext.Length];
                    if (!ext.EndsWith(".gz") && !ext.EndsWith(".zip")) break;
                }
            return stem;
        }

        /// <summary>
        /// Splits every stem into tokens of one shape, or returns null. A word with a number stuck to it
        /// (<c>DMSO1</c>) is split everywhere first; when that leaves the names ragged -- a gene such as
        /// <c>CHD6</c> beside <c>IgG</c> -- only the last part is split; then nothing is.
        /// </summary>
        private static string[][]? Align(IReadOnlyList<string> stems)
        {
            foreach (var mode in new[] { SplitMode.Every, SplitMode.Last, SplitMode.None })
            {
                var split = stems.Select(s => Tokenize(s, mode)).ToArray();
                if (split.All(t => t.Length == split[0].Length && t.Length > 0)) return split;
            }
            return null;
        }

        private enum SplitMode { Every, Last, None }

        private static string[] Tokenize(string stem, SplitMode mode)
        {
            var parts = Separators.Split(stem).Where(p => p.Length > 0).ToArray();
            var tokens = new List<string>();
            for (int i = 0; i < parts.Length; i++)
            {
                var m = WordThenNumber.Match(parts[i]);
                bool split = m.Success && (mode == SplitMode.Every || (mode == SplitMode.Last && i == parts.Length - 1));
                if (split)
                {
                    tokens.Add(m.Groups["w"].Value);
                    tokens.Add(m.Groups["n"].Value);
                }
                else tokens.Add(parts[i]);
            }
            return tokens.ToArray();
        }

        /// <summary>
        /// The role a word just before the index states, if it is a constant word the vocabulary knows.
        /// A word that itself varies (<c>DMSO</c> in <c>WT_DMSO1</c>) states nothing: the index is then an
        /// unnamed replicate count, and the caller gets <see cref="SdrfFileNameRole.Replicate"/>.
        /// </summary>
        private static SdrfFileNameRole? IndexRole(string[] first, int p, IReadOnlyList<int> factorPositions, out string named)
        {
            named = "";
            if (p == 0 || factorPositions.Contains(p - 1)) return null;
            if (!IndexWords.TryGetValue(first[p - 1], out var role)) return null;
            named = $", after the word '{first[p - 1]}'";
            return role;
        }

        /// <summary>
        /// Ranks the numbers at <paramref name="p"/> within each group of files that agree on every
        /// other design part. False unless each group's numbers are distinct and contiguous.
        /// </summary>
        private static bool TryRank(string[][] tokens, int p, IReadOnlyList<int> others, out int[] rank, out bool renumbered)
        {
            rank = new int[tokens.Length];
            renumbered = false;
            var groups = Enumerable.Range(0, tokens.Length)
                .GroupBy(i => string.Join("\u001f", others.Select(q => tokens[i][q].ToUpperInvariant())));
            foreach (var group in groups)
            {
                var values = group.Select(i => (i, v: long.Parse(tokens[i][p], CultureInfo.InvariantCulture))).OrderBy(x => x.v).ToList();
                for (int k = 0; k < values.Count; k++)
                {
                    if (values[k].v != values[0].v + k) return false;
                    rank[values[k].i] = k + 1;
                }
                if (values[0].v != 1) renumbered = true;
            }
            return true;
        }

        private static bool IsNumber(string token) => token.Length is > 0 and <= 9 && token.All(char.IsAsciiDigit);

        private static bool IsDate(string token) =>
            token.Length is 6 or 8 && token.All(char.IsAsciiDigit) &&
            DateTime.TryParseExact(token, DateFormats, CultureInfo.InvariantCulture, DateTimeStyles.None, out var d) &&
            d.Year is >= 2000 and <= 2099;

        private static IReadOnlyList<string> Sorted(IEnumerable<string> values) =>
            values.Distinct(StringComparer.OrdinalIgnoreCase).OrderBy(v => v, StringComparer.OrdinalIgnoreCase).ToList();

        private static string Distinct(int[] rank) => rank.Length == 0 ? "0" : $"1..{rank.Max()}";
    }
}
