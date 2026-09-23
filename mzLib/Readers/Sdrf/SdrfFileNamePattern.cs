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
    /// <param name="Levels">The distinct values it takes, as written; numbers in numeric order, then words.</param>
    /// <param name="Evidence">The pattern that produced the role, in words a curator can check.</param>
    /// <remarks>
    /// <see cref="Family"/> says which family of like-shaped names the slot was read in (0 = the largest).
    /// A <see cref="Position"/> of -1 marks a re-injection slot, read across the deposit rather than at a part.
    /// </remarks>
    internal sealed record SdrfFileNameSlot(int Position, SdrfFileNameRole Role, IReadOnlyList<string> Levels, string Evidence)
    {
        public int Family { get; init; }
    }

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
        bool IsControl)
    {
        /// <summary>
        /// The family of like-shaped names this file was read in (0 = the largest), or -1 for a run read
        /// alone. A caller numbering replicates within a condition numbers within a family too: two
        /// families are two sub-experiments, and their samples do not share one count.
        /// </summary>
        public int Family { get; init; } = -1;
    }

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
    /// <para><b>What a deposit looks like, as measured.</b> A blind grading of 60 PRIDE deposits against
    /// their curated SDRFs (sdrf project, pilot 1, 2026-09-23) found the first version reading structure
    /// in 6 and refusing the rest, half because the whole deposit did not share one shape. So a deposit is
    /// read as FAMILIES of like-shaped names (design runs beside blanks, standards and libraries), each on
    /// its own. A re-injection marker (<c>NEG1rep</c>, <c>C10_rr</c>, <c>_Replicate_</c>, a trailing
    /// <c>_2</c>) is read only when the unmarked twin exists. A sidecar (<c>.wiff.scan</c>) folds into its
    /// run, and a dash-separated date is one date. A number that does not count is a category when each of
    /// its values is shared, and an identifier when not -- neither refuses the family.</para>
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
            ".wiff.scan", ".raw", ".mzml", ".mzxml", ".mgf", ".wiff2", ".wiff", ".d", ".baf", ".tdf", ".gz", ".zip"
        };

        private static readonly Regex Separators = new(@"[_\-\.\s]+", RegexOptions.Compiled);
        private static readonly Regex WordThenNumber = new(@"^(?<w>[A-Za-z]+)(?<n>\d+)$", RegexOptions.Compiled);
        private static readonly Regex LettersOrDigits = new(@"\d+|[^\d]+", RegexOptions.Compiled);

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

            // A RUN is a distinct stem: a sidecar (X.wiff.scan beside X.wiff) is the same run, and a date
            // written 2018-02-28 is one token, not three numbers.
            var runOf = names.ToDictionary(n => n, n => Normalize(Stem(n)), StringComparer.OrdinalIgnoreCase);
            var runs = runOf.Values.Distinct(StringComparer.OrdinalIgnoreCase).OrderBy(r => r, StringComparer.Ordinal).ToList();
            if (runs.Count < 2)
                return None(names, "the file names differ only in their extension");

            // A re-injection is read before anything else: a run whose name is another run's name plus a
            // marker (NEG1rep beside NEG1, C10_rr beside C10, AF_17_Replicate_Control beside AF_17_Control)
            // is that run again. Only when the unmarked twin exists: rep1..rep3 with no bare twin is a count.
            var again = Reinjections(runs);
            var design = runs.Where(r => !again.ContainsKey(r)).ToList();

            // A deposit is rarely one family of names: design runs sit beside blanks, standards and library
            // runs. Each family of like-shaped names is read on its own, and a lone run is its own sample.
            List<List<string>> families = Align(design) != null
                ? new List<List<string>> { design }
                : design.GroupBy(r => Tokenize(r, SplitMode.None).Length).Select(g => g.ToList())
                    .OrderByDescending(f => f.Count).ThenBy(f => f[0], StringComparer.Ordinal).ToList();

            var slots = new List<SdrfFileNameSlot>();
            var byRun = new Dictionary<string, SdrfFileNameReading>(StringComparer.OrdinalIgnoreCase);
            var reasons = new List<string>();
            bool several = families.Count > 1 || again.Count > 0;
            for (int k = 0; k < families.Count; k++)
            {
                var family = families[k];
                string? reason = null;
                var read = family.Count < 2 ? null : ReadFamily(family, out reason);
                if (read == null)
                {
                    if (reason != null) reasons.Add(reason);
                    foreach (var r in family) byRun[r] = Alone(r);
                    continue;
                }
                slots.AddRange(read.Value.Slots.Select(x => x with { Family = k }));
                foreach (var (r, reading) in read.Value.Readings)
                    byRun[r] = (several ? reading with { SampleKey = $"F{k + 1} {reading.SampleKey}".TrimEnd() } : reading) with { Family = k };
            }

            foreach (var (run, (twin, count)) in again)
            {
                var first = byRun[twin];
                byRun[twin] = first with { TechnicalReplicate = first.TechnicalReplicate ?? 1 };
                byRun[run] = byRun[twin] with { FileName = run, TechnicalReplicate = count };
            }
            if (again.Count > 0)
                slots.Add(new SdrfFileNameSlot(-1, SdrfFileNameRole.TechnicalReplicate,
                    Sorted(again.Values.Select(v => v.Count.ToString(CultureInfo.InvariantCulture)).Prepend("1")),
                    $"{again.Count} run(s) repeat another run's name with a re-injection marker, e.g. '{again.Keys.First()}'"));

            if (slots.Count == 0)
            {
                string why = reasons.Count == 0
                    ? "no family of two or more like-shaped names was found"
                    : families.Count == 1 ? reasons[0] : $"no family of like-shaped names showed structure; the largest: {reasons[0]}";
                return None(names, why);
            }
            var files = names.Select(n => byRun[runOf[n]] with { FileName = n }).ToList();
            return new SdrfFileNameStructure(slots, files, null);
        }

        private static SdrfFileNameReading Alone(string run) =>
            new(run, run, Array.Empty<string>(), null, null, null, null, null, false);

        // Words that mark a re-injection when the name without them is also in the deposit.
        private static readonly HashSet<string> ReinjectionWords = new(StringComparer.OrdinalIgnoreCase)
        {
            "rep", "re", "repeat", "rerun", "reinj", "reinjection", "replicate", "rr", "inj", "injection"
        };

        /// <summary>Run -> (its unmarked twin, which injection this one is).</summary>
        private static Dictionary<string, (string Twin, int Count)> Reinjections(IReadOnlyList<string> runs)
        {
            const string Sep = "\u001f";
            var byTokens = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            foreach (var r in runs) byTokens.TryAdd(string.Join(Sep, Tokenize(r, SplitMode.Every)), r);
            var found = new Dictionary<string, (string Twin, int Count)>(StringComparer.OrdinalIgnoreCase);
            foreach (var r in runs)
            {
                var t = Tokenize(r, SplitMode.Every);
                for (int i = 0; i < t.Length && !found.ContainsKey(r); i++)
                {
                    int take, count;
                    if (ReinjectionWords.Contains(t[i]))
                    {
                        bool numbered = i + 1 < t.Length && IsNumber(t[i + 1]);
                        take = numbered ? 2 : 1;
                        count = numbered ? int.Parse(t[i + 1], CultureInfo.InvariantCulture) : 2;
                    }
                    else if (i == t.Length - 1 && i > 0 && IsNumber(t[i]) && t[i].Length == 1 && t[i][0] is >= '2' and <= '9')
                    {
                        // A bare trailing 2..9 on a name that also exists without it: the run again.
                        take = 1;
                        count = t[i][0] - '0';
                    }
                    else continue;
                    if (count < 2) continue;
                    var without = t.Take(i).Concat(t.Skip(i + take));
                    if (byTokens.TryGetValue(string.Join(Sep, without), out var twin) &&
                        !string.Equals(twin, r, StringComparison.OrdinalIgnoreCase))
                        found[r] = (twin, count);
                }
            }
            // A twin that is itself a re-injection is not a first injection: keep only chains of one.
            foreach (var r in found.Keys.ToList())
                if (found.ContainsKey(found[r].Twin)) found.Remove(r);
            return found;
        }

        private static readonly Regex DashedDate = new(@"(?<!\d)(20\d{2})[-_.](\d{2})[-_.](\d{2})(?!\d)", RegexOptions.Compiled);

        private static string Normalize(string stem) => DashedDate.Replace(stem, "$1$2$3");

        /// <summary>
        /// Reads ONE family of like-shaped runs. Null, with the reason, when the evidence does not decide.
        /// </summary>
        private static (List<SdrfFileNameSlot> Slots, List<(string Run, SdrfFileNameReading Reading)> Readings)? ReadFamily(
            IReadOnlyList<string> runs, out string? reason)
        {
            reason = null;
            var tokens = Align(runs);
            if (tokens == null)
            {
                reason = "the file names do not share one shape, so no part of them can be compared across the set";
                return null;
            }

            int width = tokens[0].Length;
            var varying = Enumerable.Range(0, width)
                .Where(p => tokens.Select(t => t[p]).Distinct(StringComparer.OrdinalIgnoreCase).Count() > 1)
                .ToList();
            if (varying.Count == 0) { reason = "the file names differ only in their extension"; return null; }

            var slots = new List<SdrfFileNameSlot>();
            var batchPositions = new List<int>();
            var designPositions = new List<int>();
            foreach (int p in varying)
            {
                var values = tokens.Select(t => t[p]).ToList();
                if (values.All(IsDate))
                {
                    var byDate = values.GroupBy(v => v).ToList();
                    if (byDate.Count < runs.Count && byDate.All(g => g.Count() >= 2))
                    {
                        batchPositions.Add(p);
                        slots.Add(new SdrfFileNameSlot(p, SdrfFileNameRole.Batch, Sorted(values),
                            $"a date at part {p + 1}{(values[0].Length == 6 ? " (read as yyMMdd)" : "")} splits the {runs.Count} files into {byDate.Count} groups of two or more"));
                    }
                    // A date on every file is when it was acquired, not a design: ignored, not refused.
                    continue;
                }
                designPositions.Add(p);
            }
            if (designPositions.Count == 0) { reason = "the only part that varies is an acquisition date"; return null; }

            // Two files the design cannot tell apart would be one sample written twice.
            var keys = tokens.Select(t => string.Join("\u001f", designPositions.Select(p => t[p].ToUpperInvariant()))).ToList();
            if (keys.Distinct().Count() < keys.Count) { reason = "two files read as the same design cell once the dates are set aside"; return null; }

            if (designPositions.Any(p => tokens.Any(t => IsNumber(t[p]) != IsNumber(tokens[0][p]))))
            {
                reason = "a part is a number in some names and a word in others";
                return null;
            }

            bool Shared(int p) => tokens.GroupBy(t => t[p], StringComparer.OrdinalIgnoreCase).All(g => g.Count() >= 2);
            var identity = new List<int>();   // parts that name the sample and say nothing more
            var identityWord = new Dictionary<int, string>();
            var factors = new List<int>();

            foreach (int p in designPositions.Where(p => !IsNumber(tokens[0][p])))
                (Shared(p) ? factors : identity).Add(p);

            // Numbers. A word before one may say what it counts. Of the unnamed ones only the LAST can be a
            // replicate count, and only when no named replicate exists; an earlier unnamed number is a
            // category when each of its values is shared by two or more files, and an identifier otherwise.
            var indexPositions = designPositions.Where(p => IsNumber(tokens[0][p])).ToList();
            var named = indexPositions.ToDictionary(p => p, p => IndexRole(tokens[0], p, factors, out _));
            bool namedReplicate = named.Values.Any(r => r is SdrfFileNameRole.Replicate
                or SdrfFileNameRole.BiologicalReplicate or SdrfFileNameRole.TechnicalReplicate);
            var unnamed = indexPositions.Where(p => named[p] == null).ToList();
            int? countCandidate = namedReplicate || unnamed.Count == 0 ? null : unnamed[^1];
            foreach (int p in unnamed.Where(p => p != countCandidate))
                (Shared(p) ? factors : identity).Add(p);

            var numbers = new Dictionary<int, int[]>();
            foreach (int p in indexPositions)
            {
                if (named[p] == null && p != countCandidate) continue;
                var role = named[p] ?? SdrfFileNameRole.Replicate;
                IndexRole(tokens[0], p, factors, out string word);
                if (role == SdrfFileNameRole.Sample) { identity.Add(p); identityWord[p] = word; continue; }
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
                    slots.Add(new SdrfFileNameSlot(p, role, Sorted(tokens.Select(t => t[p])),
                        $"part {p + 1} numbers the files{word}" +
                        (shift > 0 ? "; shifted by 1 for the whole deposit, because the count starts at 0" : "")));
                    continue;
                }
                var others = designPositions.Where(q => q != p).ToList();
                if (!TryRank(tokens, p, others, out var rank, out bool renumbered) || rank.Max() == 1)
                {
                    // 1..1 in every group passes "contiguous" trivially and counts nothing.
                    // Not a count: a category if shared, an identifier if not. Never a refusal of the family.
                    (Shared(p) ? factors : identity).Add(p);
                    continue;
                }
                // Alone, an unnamed number could be a run number; a word that names it settles the kind.
                if (others.Count == 0 && named[p] == null)
                {
                    reason = $"the only part that varies is a number (part {p + 1}), which cannot tell a replicate from a run number";
                    return null;
                }
                numbers[p] = rank;
                slots.Add(new SdrfFileNameSlot(p, role, Sorted(tokens.Select(t => t[p])),
                    $"part {p + 1} counts {Distinct(rank)} within each group" + word +
                    (renumbered ? "; renumbered from 1 within each group" : "")));
            }

            foreach (int p in factors.OrderBy(p => p))
            {
                var levels = tokens.Select(t => t[p]).GroupBy(v => v, StringComparer.OrdinalIgnoreCase).ToList();
                string controls = string.Join(", ", levels.Select(g => g.Key).Where(k => ControlLevel.IsMatch(k)));
                slots.Add(new SdrfFileNameSlot(p, SdrfFileNameRole.Factor, Sorted(tokens.Select(t => t[p])),
                    $"part {p + 1} takes {levels.Count} values, each shared by two or more files" +
                    (IsNumber(tokens[0][p]) ? ", read as a category because the numbers do not count" : "") +
                    (controls.Length > 0 ? $"; '{controls}' read as a control" : "")));
            }
            foreach (int p in identity.OrderBy(p => p))
                slots.Add(new SdrfFileNameSlot(p, SdrfFileNameRole.Sample, Sorted(tokens.Select(t => t[p])),
                    $"part {p + 1} names the sample" + identityWord.GetValueOrDefault(p, "")));

            if (!slots.Any(x => x.Role is not (SdrfFileNameRole.Sample or SdrfFileNameRole.Batch)))
            {
                reason = "nothing but an identifier varies, so every file is its own sample and no structure is claimed";
                return null;
            }

            slots.Sort((a, b) => a.Position.CompareTo(b.Position));
            var readings = new List<(string, SdrfFileNameReading)>(runs.Count);
            for (int i = 0; i < runs.Count; i++)
            {
                var t = tokens[i];
                int? Number(SdrfFileNameRole r) =>
                    slots.Where(s => s.Role == r && numbers.ContainsKey(s.Position)).Select(s => (int?)numbers[s.Position][i]).FirstOrDefault();
                var sampleParts = slots
                    .Where(s => s.Role is SdrfFileNameRole.Factor or SdrfFileNameRole.Sample
                        or SdrfFileNameRole.BiologicalReplicate or SdrfFileNameRole.Replicate)
                    .Select(s => numbers.TryGetValue(s.Position, out var n) ? n[i].ToString(CultureInfo.InvariantCulture) : t[s.Position]);
                var factorLevels = slots.Where(s => s.Role == SdrfFileNameRole.Factor).Select(s => t[s.Position]).ToList();
                readings.Add((runs[i], new SdrfFileNameReading(
                    runs[i],
                    string.Join(" ", sampleParts),
                    factorLevels,
                    Number(SdrfFileNameRole.Fraction),
                    Number(SdrfFileNameRole.BiologicalReplicate),
                    Number(SdrfFileNameRole.TechnicalReplicate),
                    Number(SdrfFileNameRole.Replicate),
                    batchPositions.Count > 0 ? string.Join(" ", batchPositions.Select(p => t[p])) : null,
                    factorLevels.Any(l => ControlLevel.IsMatch(l)))));
            }
            return (slots, readings);
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
                if (mode == SplitMode.Every)
                {
                    // Every run of letters and every run of digits is its own token: NEG1rep -> NEG 1 rep.
                    tokens.AddRange(LettersOrDigits.Matches(parts[i]).Select(m => m.Value));
                    continue;
                }
                var m = WordThenNumber.Match(parts[i]);
                if (m.Success && mode == SplitMode.Last && i == parts.Length - 1)
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

        // Numbers by value (1, 2, 10), then everything else by text.
        private static IReadOnlyList<string> Sorted(IEnumerable<string> values) =>
            values.Distinct(StringComparer.OrdinalIgnoreCase)
                .OrderBy(v => IsNumber(v) ? long.Parse(v, CultureInfo.InvariantCulture) : long.MaxValue)
                .ThenBy(v => v, StringComparer.OrdinalIgnoreCase).ToList();

        private static string Distinct(int[] rank) => rank.Length == 0 ? "0" : $"1..{rank.Max()}";
    }
}
