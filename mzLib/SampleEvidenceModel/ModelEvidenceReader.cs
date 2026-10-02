using System.Text;
using System.Text.Json;
using System.Text.RegularExpressions;
using Anthropic;
using Anthropic.Models.Messages;
using Readers;

namespace SampleEvidenceModel
{
    /// <summary>What a reader is given about one deposit. Every string may be empty; lists are never null.</summary>
    /// <param name="Accession">The ProteomeXchange accession.</param>
    /// <param name="RecordText">The PRIDE record's title, description and protocols.</param>
    /// <param name="PaperText">The paper's methods and tables (<see cref="PublicationText.MethodsAndTables"/>), or empty.</param>
    /// <param name="Tables">The supplement tables (<see cref="SupplementTableReader"/>).</param>
    /// <param name="RawFiles">The deposit's raw file names, as PRIDE lists them.</param>
    /// <param name="RuleEvidence">What the rule-based readers already found, so the model looks for what is missing.</param>
    internal sealed record PublicationInput(
        string Accession,
        string RecordText,
        string PaperText,
        IReadOnlyList<SupplementTable> Tables,
        IReadOnlyList<string> RawFiles,
        IReadOnlyList<SdrfEvidence> RuleEvidence);

    /// <summary>
    /// The design a publication states (the user's counting check, sdrf G36 check 3); 0 means not stated.
    /// <see cref="Plexes"/> is the number of multiplexed (TMT/iTRAQ/SILAC) plexes or mixes: there several samples share one run.
    /// Every number here was quoted from the given text by a quote that states it (<see cref="ModelEvidenceReader.Interpret"/>):
    /// a count read off the file names would make the check agree with itself (G36 rerun: PXD006430, PXD011967).
    /// </summary>
    internal sealed record StatedDesign(IReadOnlyList<(string Name, int Samples)> Groups, int Plexes, int TechnicalReplicates,
        int FractionsPerSample, int RunsStated, IReadOnlyList<string> Quotes)
    {
        /// <summary>
        /// The raw files the design implies: label-free, samples x fractions x injections; multiplexed, PLEXES x fractions x
        /// injections (PXD010429: 29 x 12 x 2 = 696). 0 when it cannot be counted: no plexes and no groups, or a label-free
        /// design with a group whose size is not stated (a partial sum would read as a mismatch).
        /// </summary>
        public int PredictedRuns => (Plexes > 0 ? Plexes : Groups.Any(g => g.Samples == 0) ? 0 : Groups.Sum(g => g.Samples)) is var units and > 0
            ? units * Math.Max(1, TechnicalReplicates) * Math.Max(1, FractionsPerSample)
            : 0;
    }

    /// <summary>A model reading: the claims that passed validation, the ones that did not and why, and what it cost.</summary>
    internal sealed record ModelEvidenceResult(
        IReadOnlyList<SdrfEvidence> Evidence,
        IReadOnlyList<string> Rejected,
        StatedDesign? Design,
        string StopReason,
        long InputTokens,
        long OutputTokens,
        long CacheReadTokens);

    /// <summary>
    /// The OPT-IN model reader for publication evidence (sdrf design SAMPLE-EVIDENCE.md, E5; D42 rules first, model
    /// opt-in; D43 its own project on the Anthropic SDK, <see cref="DefaultModel"/>). It reads what the rules cannot --
    /// channel assignments written in a Methods paragraph, a column headed "group" that holds sex, the design a paper
    /// states -- and returns <see cref="SdrfEvidence"/> with method <c>model</c>.
    ///
    /// <para><b>Nothing it says is taken on trust.</b> Every claim must quote the text that states it, and the quote
    /// must occur in the text the model was given; the file must be one of the deposit's raw files (or none, for the
    /// whole deposit); the column must be an SDRF column. A claim that fails any of these is rejected and reported in
    /// <see cref="ModelEvidenceResult.Rejected"/>, never written. What passes still carries <c>source method = model</c>
    /// (aging 038), so a consumer can filter it with a column test.</para>
    ///
    /// <para><b>What the quote check does not stop.</b> It stops invented text, not text the paper itself plants: a
    /// sentence in a deposit's own paper or supplement can steer the model (prompt injection), and a claim quoting that
    /// sentence passes. Treat a model reading as untrusted input like the paper it came from: it is evidence for review,
    /// written with <c>source method = model</c>, never a fact the drafter acts on unseen.</para>
    ///
    /// <para>Costs money and needs a key (<c>ANTHROPIC_API_KEY</c> or <c>ant auth login</c>), so only a caller that opts
    /// in constructs one; its answers are meant to be saved as evidence and read once per deposit.</para>
    /// </summary>
    internal sealed class ModelEvidenceReader
    {
        /// <summary>The model, chosen for extraction accuracy (D43).</summary>
        internal const string DefaultModel = "claude-opus-5";

        private readonly AnthropicClient _client;
        private readonly string _model;

        internal ModelEvidenceReader(AnthropicClient client, string model = DefaultModel)
        {
            _client = client ?? throw new ArgumentNullException(nameof(client));
            _model = string.IsNullOrWhiteSpace(model) ? throw new ArgumentException("A model is required.", nameof(model)) : model;
        }

        /// <summary>Reads one deposit's publication. Throws the SDK's exceptions on a transport or API failure.</summary>
        internal async Task<ModelEvidenceResult> ReadAsync(PublicationInput input, CancellationToken cancellationToken = default)
        {
            if (input == null) throw new ArgumentNullException(nameof(input));
            // The beta endpoint, for the server-side refusal fallback: should a safety classifier decline (a paper on a
            // pathogen, say), the request is re-served by the fallback model inside the same call instead of stopping.
            // Streamed, so the model's whole output limit can be used without the request timing out: a deposit with a
            // thousand channel claims needs far more than a non-streamed call's 16k (Alexander-Sol on #1377, PXD007160).
            var response = await _client.Beta.Messages.CreateStreaming(new Anthropic.Models.Beta.Messages.MessageCreateParams
            {
                Model = _model,
                MaxTokens = MaxOutputTokens,
                System = new List<Anthropic.Models.Beta.Messages.BetaTextBlockParam>
                    { new() { Text = SystemPrompt, CacheControl = new Anthropic.Models.Beta.Messages.BetaCacheControlEphemeral() } },
                Messages = [new() { Role = Anthropic.Models.Beta.Messages.Role.User, Content = UserPrompt(input) }],
                Thinking = new Anthropic.Models.Beta.Messages.BetaThinkingConfigAdaptive(),
                OutputConfig = new Anthropic.Models.Beta.Messages.BetaOutputConfig
                {
                    Effort = Anthropic.Models.Beta.Messages.Effort.High,
                    Format = new Anthropic.Models.Beta.Messages.BetaJsonOutputFormat { Schema = Schema }
                },
                // "default": the server picks the fallback by refusal category, so no model list is maintained here.
                Betas = ["server-side-fallback-2026-07-01"],
                Fallbacks = new Anthropic.Models.Beta.Messages.Default(),
            }, cancellationToken).Aggregate().ConfigureAwait(false);

            // Raw(), the wire string: the SDK's enum wrapper does not ToString() to "refusal", so a refusal's text was read.
            string stop = response.StopReason?.Raw() ?? "";
            var text = new StringBuilder();
            foreach (var block in response.Content)
                if (block.TryPickText(out var t)) text.Append(t.Text);
            // A cut-off answer is JSON that ends mid-claim: say so, rather than "not JSON", so the caller knows to split it.
            var result = stop == "max_tokens"
                ? new ModelEvidenceResult(Array.Empty<SdrfEvidence>(), new[] { $"no reading: the answer was cut off at the output limit ({MaxOutputTokens} tokens), so nothing in it was kept; read the deposit in parts (per table or per plex)" }, null, stop, 0, 0, 0)
                : stop == "refusal" || text.Length == 0
                ? new ModelEvidenceResult(Array.Empty<SdrfEvidence>(), new[] { $"no reading: stop reason '{stop}'" }, null, stop, 0, 0, 0)
                : Interpret(text.ToString(), input, stop);
            return result with
            {
                InputTokens = response.Usage.InputTokens,
                OutputTokens = response.Usage.OutputTokens,
                CacheReadTokens = response.Usage.CacheReadInputTokens ?? 0
            };
        }

        // ---------------- the prompt ----------------

        internal const string SystemPrompt = """
            You read the publication behind a proteomics data deposit (ProteomeXchange/PRIDE) to fill its SDRF-Proteomics
            sample sheet: one row per raw mass-spectrometry file (and, for TMT/iTRAQ, per file and channel), describing the
            biological sample in it. You are given the deposit's record text, its paper's methods and tables, its
            supplementary tables (each row labelled like [mmc2.xlsx!Sheet1!R14]), its raw file names, and what a rule-based
            reader already extracted.

            Return claims about the samples that were measured by mass spectrometry in THIS deposit, and the design the
            publication states.

            Claims:
            - Claim only what the given text states about these samples. Background, other experiments, other datasets,
              references and reagents are not claims.
            - quote: copy, verbatim, the shortest passage (at most about 300 characters) of the given text that states the
              claim. For a table, quote the row's cells as they appear and put the row label, e.g. [mmc2.xlsx!Sheet1!R14],
              at the start of the quote.
            - data_file: one raw file name exactly as listed, or "" when the value holds for every file of the deposit.
              Assign per-file values only when the text links a sample to a file (a sample ID that appears in the file
              name, an explicit table, a stated run order). Never assign by guessing an order.
            - data_file_pattern: when the value holds for a SET of the files (the TMT6 files, Set A, one tissue), a glob
              over the file names (* any characters, ? one), with data_file "". It must match at least one listed file.
              Otherwise "".
            - label: for isobaric labelling, the channel as TMT126, TMT127N, TMTpro134C or iTRAQ114; for SILAC, the
              state as SILAC light, SILAC medium or SILAC heavy (one claim per state, each for its own sample);
              otherwise "".
            - column: an SDRF column: "source name", "characteristics[<name>]" with name one of organism part, disease,
              cell type, cell line, sex, age, developmental stage, individual, strain, genotype, treatment, compound,
              dose, time, biological replicate, phenotype, or "comment[fraction identifier]",
              "comment[technical replicate]", or "factor value[<name>]" for the variable the study compares.
            - value: age as 45Y, 9M, 8W or 3D; sex as male or female; replicate and fraction numbers as whole numbers
              from 1; anything else as the text writes it.
            - Do not repeat what the rule-based reader already found; add what is missing and correct nothing silently.
            - confidence: "likely" when the text states it for these samples; "guess" when you inferred it.

            Design: the groups the study compares with the number of biological samples in each; for multiplexed labelling
            (TMT, iTRAQ, SILAC) the number of plexes or mixes (several samples share one run); technical replicates (injections) per sample or plex;
            fractions per sample or plex; and the number of runs if stated. Give each number its own quote, copied
            verbatim, that contains that number (as digits or a word such as "twelve", "twice" or "triplicate"). A number
            the text does not state is 0 with an empty quote: never count files, file names or table rows to get one.
            Use an empty list of groups when the text states none.

            If nothing can be claimed, return an empty list of claims.
            """;

        internal static string UserPrompt(PublicationInput input)
        {
            var sb = new StringBuilder();
            sb.Append("# Deposit ").Append(input.Accession).Append("\n\n## PRIDE record\n").Append(input.RecordText).Append("\n\n");
            sb.Append("## Paper: methods and tables\n").Append(input.PaperText.Length == 0 ? "(no full text available)"
                : input.PaperText.Length <= MaxPaperChars ? input.PaperText
                : input.PaperText[..MaxPaperChars] + $"\n[paper cut at {MaxPaperChars:N0} of {input.PaperText.Length:N0} characters]").Append("\n\n");
            sb.Append("## Supplementary tables\n").Append(input.Tables.Count > 0 ? PublicationText.Tables(input.Tables) : "(none)").Append("\n\n");
            sb.Append($"## Raw files ({input.RawFiles.Count})\n");
            sb.Append(PublicationText.Files(input.RawFiles));
            sb.Append("\n## Already found by the rule-based reader\n");
            if (input.RuleEvidence.Count == 0) sb.Append("(nothing)\n");
            foreach (var e in input.RuleEvidence.Take(300))
                sb.Append($"{(e.DataFile.Length > 0 ? e.DataFile : "(all files)")}\t{e.Label}\t{e.Column}\t{e.Value}\n");
            return sb.ToString();
        }

        private static readonly Dictionary<string, JsonElement> Schema = BuildSchema();

        private static Dictionary<string, JsonElement> BuildSchema()
        {
            var str = new { type = "string" };
            var integer = new { type = "integer" };
            var claim = new
            {
                type = "object",
                additionalProperties = false,
                required = new[] { "data_file", "data_file_pattern", "label", "column", "value", "source", "quote", "confidence" },
                properties = new
                {
                    data_file = str, data_file_pattern = str, label = str, column = str, value = str,
                    source = new { type = "string", @enum = new[] { "paper", "supplement", "pride record" } },
                    quote = str,
                    confidence = new { type = "string", @enum = new[] { "likely", "guess" } }
                }
            };
            var group = new { type = "object", additionalProperties = false, required = new[] { "name", "samples", "quote" }, properties = new { name = str, samples = integer, quote = str } };
            var design = new
            {
                type = "object",
                additionalProperties = false,
                required = new[] { "groups", "plexes", "plexes_quote", "technical_replicates", "technical_replicates_quote",
                    "fractions_per_sample", "fractions_quote", "runs_stated", "runs_quote" },
                properties = new
                {
                    groups = new { type = "array", items = group },
                    plexes = integer, plexes_quote = str,
                    technical_replicates = integer, technical_replicates_quote = str,
                    fractions_per_sample = integer, fractions_quote = str,
                    runs_stated = integer, runs_quote = str
                }
            };
            return new Dictionary<string, JsonElement>
            {
                ["type"] = JsonSerializer.SerializeToElement("object"),
                ["additionalProperties"] = JsonSerializer.SerializeToElement(false),
                ["required"] = JsonSerializer.SerializeToElement(new[] { "claims", "design" }),
                ["properties"] = JsonSerializer.SerializeToElement(new { claims = new { type = "array", items = claim }, design }),
            };
        }

        // ---------------- validation ----------------

        private static readonly Regex Column = new(@"^(source name|(characteristics|factor value)\[[a-z][a-z0-9 ]*\]|comment\[(fraction identifier|technical replicate)\])$", RegexOptions.Compiled);
        private static readonly Regex Label = new(@"^(TMT(pro)?1[23]\d[NC]?|iTRAQ(11[3-9]|121)|SILAC (light|medium|heavy))$", RegexOptions.Compiled);
        private static readonly Regex RowRef = new(@"\[([^\[\]]+![^\[\]]*R\d+)\]", RegexOptions.Compiled);

        /// <summary>
        /// Turns the model's JSON into evidence, keeping only claims whose quote occurs in the given text, whose file is
        /// the deposit's, and whose column and label are well formed. Pure: no I/O.
        /// </summary>
        internal static ModelEvidenceResult Interpret(string json, PublicationInput input, string stopReason = "end_turn")
        {
            var rejected = new List<string>();
            JsonElement root;
            try { root = JsonDocument.Parse(json).RootElement; }
            catch (JsonException e) { return new(Array.Empty<SdrfEvidence>(), new[] { "not JSON: " + e.Message }, null, stopReason, 0, 0, 0); }

            string corpus = Norm(string.Join("\n", input.RecordText, input.PaperText, PublicationText.Tables(input.Tables, int.MaxValue, int.MaxValue, int.MaxValue, int.MaxValue)));
            var files = input.RawFiles.ToDictionary(f => f, f => f, StringComparer.OrdinalIgnoreCase);
            // Each row label, with that row's own cells: a quote citing a row must be whole cells of THAT row.
            var refs = new Dictionary<string, string[]>(StringComparer.Ordinal);
            foreach (var t in input.Tables)
                for (int k = 0; k < t.Rows.Count; k++)
                    refs.TryAdd($"{(t.Sheet.Length > 0 ? $"{t.File}!{t.Sheet}" : t.File)}!R{t.RowNumbers[k]}", t.Rows[k].Select(Norm).ToArray());
            var evidence = new List<SdrfEvidence>();
            var seen = new HashSet<(string, string, string, string)>();

            if (root.TryGetProperty("claims", out var claims) && claims.ValueKind == JsonValueKind.Array)
                foreach (var c in claims.EnumerateArray())
                {
                    string Get(string name) => c.TryGetProperty(name, out var v) && v.ValueKind == JsonValueKind.String ? v.GetString()!.Trim() : "";
                    string file = Get("data_file"), pattern = Get("data_file_pattern"), label = Get("label"), column = Get("column").ToLowerInvariant(),
                        value = Get("value"), source = Get("source"), quote = Get("quote"), confidence = Get("confidence");
                    string what = $"{(file.Length > 0 ? file : pattern.Length > 0 ? pattern : "(all files)")} {label} {column} = '{value}'";

                    if (value.Length == 0) { rejected.Add($"{what}: no value"); continue; }
                    if (!Column.IsMatch(column)) { rejected.Add($"{what}: not an SDRF column"); continue; }
                    if (label.Length > 0 && !Label.IsMatch(label)) { rejected.Add($"{what}: not a channel label"); continue; }
                    if (file.Length > 0 && !files.TryGetValue(file, out file!)) { rejected.Add($"{what}: not one of the deposit's raw files"); continue; }
                    if (pattern.Length > 0 && file.Length > 0) { rejected.Add($"{what}: names both a file and a file pattern"); continue; }
                    // A model-written glob becomes a backtracking regex; many wildcards that fail to match take minutes.
                    if (pattern.Count(ch => ch == '*') > MaxPatternWildcards) { rejected.Add($"{what}: the file pattern has more than {MaxPatternWildcards} '*' wildcards"); continue; }
                    if (pattern.Length > 0 && !input.RawFiles.Any(f => SdrfEvidence.GlobMatches(pattern, f))) { rejected.Add($"{what}: the file pattern matches none of the deposit's raw files"); continue; }
                    string rowRef = RowRef.Match(quote) is { Success: true } m && refs.ContainsKey(m.Groups[1].Value) ? m.Groups[1].Value : "";
                    string bare = Norm(RowRef.Replace(quote, " "));
                    // A row's cells can be short ("P1 | F"), so a cited row anchors a short quote; free text needs a phrase.
                    if (bare.Length < (rowRef.Length > 0 ? 1 : MinQuoteLength)) { rejected.Add($"{what}: the quote is too short to check"); continue; }
                    // A row quote is WHOLE cells of the row it cites: "male" is not a cell of "female", nor "7" of "77".
                    if (rowRef.Length > 0 && !WholeCells(bare, refs[rowRef])) { rejected.Add($"{what}: the quote is not whole cells of the row it cites"); continue; }
                    if (rowRef.Length == 0 && !corpus.Contains(bare, StringComparison.Ordinal)) { rejected.Add($"{what}: the quote is not in the given text"); continue; }
                    // The value is in the key: two claims that disagree both reach the drafter, which fills nothing from them.
                    if (!seen.Add((file + "|" + pattern, label, column, value.ToLowerInvariant()))) continue;

                    // Real text is not enough: a quote that does not state the value is kept only for review.
                    bool stated = QuoteStates(bare, column, value);
                    // A row describes one sample: it supports a per-file claim only when it names that file (or, for a
                    // pattern, one of its files). Another sample's row is kept for review, as a guess.
                    if (rowRef.Length > 0 && (file.Length > 0 || pattern.Length > 0)
                        && !(file.Length > 0 ? new[] { file } : input.RawFiles.Where(f => SdrfEvidence.GlobMatches(pattern, f))).Any(f => RowNamesFile(refs[rowRef], f)))
                        stated = false;
                    string locator = rowRef.Length > 0 ? rowRef : $"{source}: \"{Truncate(quote, 160)}\"";
                    evidence.Add(new SdrfEvidence(file, label, column, value,
                        source is "paper" or "supplement" or "pride record" ? source : "paper", locator, "model",
                        confidence == "guess" || !stated ? SdrfEvidenceConfidence.Guess : SdrfEvidenceConfidence.Likely, pattern));
                }

            StatedDesign? design = null;
            if (root.TryGetProperty("design", out var d) && d.ValueKind == JsonValueKind.Object)
            {
                var quotes = new List<string>();
                // A number is kept only with a quote that is in the given text, states that number, and (for what it
                // counts) says what it counts near it: "twelve" alone, or "twelve donors" for fractions, proves nothing.
                int Stated(string what, int n, string quote, DesignKeyword? keyword = null)
                {
                    if (n <= 0) return 0;
                    string q = Norm(quote);
                    string? why = q.Length < MinQuoteLength ? "the quote is too short to check"
                        : !corpus.Contains(q, StringComparison.Ordinal) ? "the quote is not in the given text"
                        : !States(q, n) ? "the quote does not state it"
                        : keyword != null && !NumberSpans(q, n, idioms: true).Any(span => keyword.Pattern.IsMatch(Window(q, span)))
                            ? $"the quote does not say '{keyword.Name}' near the number"
                        : null;
                    if (why == null)
                    {
                        if (!quotes.Contains(quote)) quotes.Add(quote);
                        return n;
                    }
                    rejected.Add($"design {what} = {n}: {why}");
                    return 0;
                }
                static string Str(JsonElement e, string n) => e.TryGetProperty(n, out var v) && v.ValueKind == JsonValueKind.String ? v.GetString()!.Trim() : "";
                static int Int(JsonElement e, string n) => e.TryGetProperty(n, out var v) && v.ValueKind == JsonValueKind.Number && v.TryGetInt32(out int i) ? Math.Max(0, i) : 0;

                var groups = new List<(string Name, int Samples)>();
                if (d.TryGetProperty("groups", out var g) && g.ValueKind == JsonValueKind.Array)
                    foreach (var x in g.EnumerateArray().Where(x => x.ValueKind == JsonValueKind.Object && Str(x, "name").Length > 0))
                        groups.Add((Str(x, "name"), Stated($"group '{Str(x, "name")}'", Int(x, "samples"), Str(x, "quote"))));
                int plexes = Stated("plexes", Int(d, "plexes"), Str(d, "plexes_quote"), PlexWords);
                int tech = Stated("technical replicates", Int(d, "technical_replicates"), Str(d, "technical_replicates_quote"), InjectionWords);
                int fractions = Stated("fractions", Int(d, "fractions_per_sample"), Str(d, "fractions_quote"), FractionWords);
                int runs = Stated("runs", Int(d, "runs_stated"), Str(d, "runs_quote"), RunWords);
                if (groups.Any(x => x.Samples > 0) || plexes > 0 || runs > 0)
                    design = new StatedDesign(groups, plexes, tech, fractions, runs, quotes);
            }
            return new ModelEvidenceResult(evidence, rejected, design, stopReason, 0, 0, 0);
        }

        /// <summary>The shortest free-text quote checked: shorter ("and", "Human") occurs in any text and proves nothing.</summary>
        internal const int MinQuoteLength = 10;

        /// <summary>The most '*' wildcards a model's file pattern may carry (see <see cref="SdrfEvidence.GlobMatches"/>).</summary>
        internal const int MaxPatternWildcards = 4;

        /// <summary>The most characters of the paper put in the prompt, like the tables' 200k cap.</summary>
        internal const int MaxPaperChars = 200_000;

        /// <summary>The output limit of <see cref="DefaultModel"/>, covering both its thinking and its answer.</summary>
        internal const int MaxOutputTokens = 128_000;

        private static readonly Regex AgeOrNumber = new(@"^(\d+)(\.\d+)?\s*([ymwd])?$", RegexOptions.Compiled | RegexOptions.IgnoreCase);
        private static readonly Regex Word = new(@"[a-z0-9]+", RegexOptions.Compiled);
        private static readonly HashSet<string> StopWords = new(StringComparer.Ordinal) { "and", "the", "with", "from", "for", "not", "of", "in" };

        /// <summary>Words that say nothing about which value is meant ("cell type", "tissue sample"), so they state nothing.</summary>
        private static readonly HashSet<string> GenericWords = new(StringComparer.Ordinal) { "type", "cell", "cells", "tissue", "sample", "samples" };

        /// <summary>The endings a quoted word may add to a value's stem: "ventricle" is stated by "ventricular", "normal" not by "normalized".</summary>
        private static readonly HashSet<string> Endings = new(StringComparer.Ordinal) { "", "s", "es", "al", "ar", "ic", "ial", "ous", "ular" };

        /// <summary>
        /// Whether a (normalised) quote states the claim's value: sex by its words or letter code; an age by its number
        /// with an age unit or "aged" beside it, or as a table cell of its own; a count by its number (<see cref="States"/>);
        /// anything else by at least half of the value's telling words, each as a whole word or an inflection of it.
        /// </summary>
        internal static bool QuoteStates(string quote, string column, string value)
        {
            string q = Norm(quote), v = Norm(value);
            if (column == "characteristics[sex]" && v is "male" or "female")
                return Regex.IsMatch(q, v == "male" ? @"\b(m|males?|man|men|boys?)\b" : @"\b(f|females?|wom[ae]n|girls?)\b");
            if (AgeOrNumber.Match(v) is { Success: true } n)
            {
                string number = n.Groups[1].Value + n.Groups[2].Value;
                bool isAge = n.Groups[3].Success || column.EndsWith("[age]", StringComparison.Ordinal);
                if (!isAge)
                    return n.Groups[2].Success ? Regex.IsMatch(q, $@"(?<![\d.,]){Regex.Escape(number)}(?!\d|[.,]\d)")
                        : int.TryParse(number, out int k) && States(q, k);
                return StatesAge(q, number, n.Groups[3].Success ? char.ToLowerInvariant(n.Groups[3].Value[0]) : 'y');
            }
            var words = Word.Matches(v).Select(w => w.Value).Where(w => w.Length >= 3 && !StopWords.Contains(w) && !GenericWords.Contains(w)).ToList();
            if (words.Count == 0) return Regex.IsMatch(q, $@"(?<![a-z0-9]){Regex.Escape(v)}(?![a-z0-9])");
            var quoted = Word.Matches(q).Select(w => w.Value).ToHashSet(StringComparer.Ordinal);
            return words.Count(w => quoted.Any(x => Inflects(w, x))) * 2 >= words.Count;
        }

        /// <summary>Whether quoted word <paramref name="x"/> is value word <paramref name="w"/> or an inflection of it.</summary>
        private static bool Inflects(string w, string x)
        {
            if (x == w) return true;
            int common = 0;
            while (common < w.Length && common < x.Length && w[common] == x[common]) common++;
            return common >= 5 && common >= w.Length - 2 && Endings.Contains(x[common..]);
        }

        private static readonly Dictionary<char, string> AgeUnits = new()
        {
            ['y'] = @"(y|yrs?|years?|y/o|yo)",
            ['m'] = @"(mo|mos|months?)",
            ['w'] = @"(wks?|weeks?|w)",
            ['d'] = @"(d|days?)",
        };

        /// <summary>
        /// Whether a quote states an age: the number with its unit after it ("35 years", "a 35-year-old"), after "aged"
        /// or "age" for years, or as a whole table cell. A bare number states nothing: "a 35 min gradient" is not 35Y.
        /// </summary>
        private static bool StatesAge(string q, string number, char unit)
        {
            if (q.Split(" | ").Any(cell => cell.Trim() == number)) return true;
            var unitAfter = new Regex($@"^\s*-?\s*{AgeUnits[unit]}(?![a-z])");
            var ageBefore = new Regex(@"\bage[ds]?\s*(of|at|:|=)?\s*$");
            var spans = int.TryParse(number, out int k) ? NumberSpans(q, k, idioms: false)
                : Regex.Matches(q, $@"(?<![\d.,]){Regex.Escape(number)}(?!\d|[.,]\d)").Select(m => (m.Index, m.Length));
            return spans.Any(span => unitAfter.IsMatch(q[(span.Index + span.Length)..])
                || unit == 'y' && ageBefore.IsMatch(q[..span.Index]));
        }

        /// <summary>
        /// Whether a quote names, as one of its cells, the file a claim is about: the file's stem, or a sample id with a
        /// letter and a digit that stands as a whole token in the stem ("S001" in "Plasma_S001_F01"). A row about
        /// another sample does not support a claim about this one.
        /// </summary>
        private static bool RowNamesFile(IEnumerable<string> cells, string file)
        {
            string stem = Norm(Path.GetFileNameWithoutExtension(file));
            static bool Token(string inner, string outer) => Regex.IsMatch(outer, $@"(?<![a-z0-9]){Regex.Escape(inner)}(?![a-z0-9])");
            return cells.Select(c => c.Trim()).Where(c => c.Length > 0).Any(c => c == stem
                || c.Any(char.IsLetter) && c.Any(char.IsDigit) && (Token(c, stem) || Token(stem, c)));
        }

        /// <summary>Whether a row quote is one or more of the cited row's cells, each whole.</summary>
        private static bool WholeCells(string bare, IReadOnlyCollection<string> cells) =>
            bare.Split(" | ").Select(p => p.Trim()).ToList() is var parts && parts.All(p => p.Length > 0 && cells.Contains(p));

        private static readonly string[] NumberWords =
            { "zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve",
              "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eighteen", "nineteen", "twenty" };

        /// <summary>Whether a quote states the number: as digits standing alone, as a word, or as once/twice/triplicate.</summary>
        internal static bool States(string quote, int n) => NumberSpans(Norm(quote), n, idioms: true).Any();

        /// <summary>Where a (normalised) quote states the number: digits standing alone, a number word, and with
        /// <paramref name="idioms"/> once/twice/triplicate and the like.</summary>
        private static IEnumerable<(int Index, int Length)> NumberSpans(string q, int n, bool idioms)
        {
            var patterns = new List<string> { $@"(?<![\d.,]){n}(?!\d|[.,]\d)" };
            if (n < NumberWords.Length) patterns.Add($@"\b{NumberWords[n]}\b");
            if (idioms)
                patterns.Add(n switch
                {
                    1 => @"\b(once|single|singly)\b",
                    2 => @"\b(twice|duplicates?|pairs?)\b",
                    3 => @"\b(thrice|triplicates?)\b",
                    4 => @"\bquadruplicates?\b",
                    _ => @"(?!)"
                });
            return patterns.SelectMany(p => Regex.Matches(q, p).Select(m => (m.Index, m.Length)));
        }

        /// <summary>What a design number counts, and the words that must stand near it to say so.</summary>
        private sealed record DesignKeyword(string Name, Regex Pattern);

        private static readonly DesignKeyword FractionWords = new("fraction", new(@"fraction", RegexOptions.Compiled));
        private static readonly DesignKeyword InjectionWords = new("injection or replicate", new(@"inject|replicat|plicate|twice|thrice|once|\brun|analy[sz]|measur|acquir", RegexOptions.Compiled));
        private static readonly DesignKeyword PlexWords = new("plex or set", new(@"plex|\bsets?\b|\bmix|batch|tmt|itraq|silac|label", RegexOptions.Compiled));
        private static readonly DesignKeyword RunWords = new("run", new(@"\bruns?\b|file|acquisition|injection|measurement|analys", RegexOptions.Compiled));

        /// <summary>The quote around a number, the words that can say what it counts ("separated into twelve fractions").</summary>
        private static string Window(string q, (int Index, int Length) span, int reach = 40)
        {
            int from = Math.Max(0, span.Index - reach), to = Math.Min(q.Length, span.Index + span.Length + reach);
            return q[from..to];
        }

        /// <summary>
        /// Compatibility-normalised (NFKC: a PDF's "ﬁ" ligature is "fi"), soft hyphens and zero-width characters removed,
        /// lower case, one space for any run of whitespace, typographic dashes and quotes made plain.
        /// </summary>
        private static string Norm(string s) =>
            Regex.Replace(Regex.Replace(s.Normalize(NormalizationForm.FormKC), "[\u00AD\u200B\u200C\u200D\uFEFF]", "")
                .Replace('‐', '-').Replace('‑', '-').Replace('–', '-').Replace('—', '-')
                .Replace('‘', '\'').Replace('’', '\'').Replace('“', '"').Replace('”', '"'),
                @"\s+", " ").Trim().ToLowerInvariant();

        private static string Truncate(string s, int n) => s.Length <= n ? s : s[..n] + "...";
    }
}
