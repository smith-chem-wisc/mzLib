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
    /// </summary>
    internal sealed record StatedDesign(IReadOnlyList<(string Name, int Samples)> Groups, int Plexes, int TechnicalReplicates,
        int FractionsPerSample, int RunsStated, string Quote)
    {
        /// <summary>
        /// The raw files the design implies: label-free, samples x fractions x injections; multiplexed, PLEXES x fractions x
        /// injections (PXD010429: 29 x 12 x 2 = 696). 0 when neither groups nor plexes are stated.
        /// </summary>
        public int PredictedRuns => (Plexes > 0 ? Plexes : Groups.Sum(g => g.Samples)) is var units and > 0
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
            var response = await _client.Beta.Messages.Create(new Anthropic.Models.Beta.Messages.MessageCreateParams
            {
                Model = _model,
                MaxTokens = 16000,
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
            }, cancellationToken).ConfigureAwait(false);

            string stop = response.StopReason?.ToString() ?? "";
            var text = new StringBuilder();
            foreach (var block in response.Content)
                if (block.TryPickText(out var t)) text.Append(t.Text);
            var result = stop == "refusal" || text.Length == 0
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
            fractions per sample or plex; and the number of runs if stated, with the quote that states it. Use 0 for
            anything not stated, and an empty list of groups when the text states none.

            If nothing can be claimed, return an empty list of claims.
            """;

        internal static string UserPrompt(PublicationInput input)
        {
            var sb = new StringBuilder();
            sb.Append("# Deposit ").Append(input.Accession).Append("\n\n## PRIDE record\n").Append(input.RecordText).Append("\n\n");
            sb.Append("## Paper: methods and tables\n").Append(input.PaperText.Length > 0 ? input.PaperText : "(no full text available)").Append("\n\n");
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
                required = new[] { "data_file", "label", "column", "value", "source", "quote", "confidence" },
                properties = new
                {
                    data_file = str, label = str, column = str, value = str,
                    source = new { type = "string", @enum = new[] { "paper", "supplement", "pride record" } },
                    quote = str,
                    confidence = new { type = "string", @enum = new[] { "likely", "guess" } }
                }
            };
            var group = new { type = "object", additionalProperties = false, required = new[] { "name", "samples" }, properties = new { name = str, samples = integer } };
            var design = new
            {
                type = "object",
                additionalProperties = false,
                required = new[] { "groups", "plexes", "technical_replicates", "fractions_per_sample", "runs_stated", "quote" },
                properties = new
                {
                    groups = new { type = "array", items = group },
                    plexes = integer,
                    technical_replicates = integer, fractions_per_sample = integer, runs_stated = integer, quote = str
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
            var refs = new HashSet<string>(input.Tables.SelectMany(t => Enumerable.Range(0, t.Rows.Count)
                .Select(k => $"{(t.Sheet.Length > 0 ? $"{t.File}!{t.Sheet}" : t.File)}!R{t.RowNumbers[k]}")), StringComparer.Ordinal);
            var evidence = new List<SdrfEvidence>();
            var seen = new HashSet<(string, string, string)>();

            if (root.TryGetProperty("claims", out var claims) && claims.ValueKind == JsonValueKind.Array)
                foreach (var c in claims.EnumerateArray())
                {
                    string Get(string name) => c.TryGetProperty(name, out var v) && v.ValueKind == JsonValueKind.String ? v.GetString()!.Trim() : "";
                    string file = Get("data_file"), label = Get("label"), column = Get("column").ToLowerInvariant(), value = Get("value"),
                        source = Get("source"), quote = Get("quote"), confidence = Get("confidence");
                    string what = $"{(file.Length > 0 ? file : "(all files)")} {label} {column} = '{value}'";

                    if (value.Length == 0) { rejected.Add($"{what}: no value"); continue; }
                    if (!Column.IsMatch(column)) { rejected.Add($"{what}: not an SDRF column"); continue; }
                    if (label.Length > 0 && !Label.IsMatch(label)) { rejected.Add($"{what}: not a channel label"); continue; }
                    if (file.Length > 0 && !files.TryGetValue(file, out file!)) { rejected.Add($"{what}: not one of the deposit's raw files"); continue; }
                    string rowRef = RowRef.Match(quote) is { Success: true } m && refs.Contains(m.Groups[1].Value) ? m.Groups[1].Value : "";
                    string bare = Norm(RowRef.Replace(quote, " "));
                    if (bare.Length < 3 || !corpus.Contains(bare, StringComparison.Ordinal)) { rejected.Add($"{what}: the quote is not in the given text"); continue; }
                    if (!seen.Add((file, label, column))) continue;

                    string locator = rowRef.Length > 0 ? rowRef : $"{source}: \"{Truncate(quote, 160)}\"";
                    evidence.Add(new SdrfEvidence(file, label, column, value,
                        source is "paper" or "supplement" or "pride record" ? source : "paper", locator, "model",
                        confidence == "guess" ? SdrfEvidenceConfidence.Guess : SdrfEvidenceConfidence.Likely));
                }

            StatedDesign? design = null;
            if (root.TryGetProperty("design", out var d) && d.ValueKind == JsonValueKind.Object)
            {
                int Int(string n) => d.TryGetProperty(n, out var v) && v.ValueKind == JsonValueKind.Number && v.TryGetInt32(out int i) ? Math.Max(0, i) : 0;
                var groups = d.TryGetProperty("groups", out var g) && g.ValueKind == JsonValueKind.Array
                    ? g.EnumerateArray().Select(x => (Name: x.TryGetProperty("name", out var n) ? n.GetString() ?? "" : "",
                        Samples: x.TryGetProperty("samples", out var s) && s.TryGetInt32(out int k) ? Math.Max(0, k) : 0)).ToList()
                    : new List<(string, int)>();
                string quote = d.TryGetProperty("quote", out var q) ? q.GetString() ?? "" : "";
                // A design is only as good as its quote.
                if (groups.Count > 0 && Norm(quote).Length >= 3 && corpus.Contains(Norm(quote), StringComparison.Ordinal))
                    design = new StatedDesign(groups, Int("plexes"), Int("technical_replicates"), Int("fractions_per_sample"), Int("runs_stated"), quote);
                else if (groups.Count > 0)
                    rejected.Add("design: the quote is not in the given text");
            }
            return new ModelEvidenceResult(evidence, rejected, design, stopReason, 0, 0, 0);
        }

        /// <summary>Lower case, one space for any run of whitespace, typographic dashes and quotes made plain.</summary>
        private static string Norm(string s) =>
            Regex.Replace(s.Replace('‐', '-').Replace('‑', '-').Replace('–', '-').Replace('—', '-')
                .Replace('‘', '\'').Replace('’', '\'').Replace('“', '"').Replace('”', '"').Replace(' ', ' '),
                @"\s+", " ").Trim().ToLowerInvariant();

        private static string Truncate(string s, int n) => s.Length <= n ? s : s[..n] + "...";
    }
}
