using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using Anthropic;
using NUnit.Framework;
using Readers;
using SampleEvidenceModel;

namespace Test.FileReadingTests
{
    /// <summary>
    /// The opt-in model reader's safeguards, offline: nothing the model says is written unless its quote is in the text
    /// it was given, its file is the deposit's and its column is an SDRF column (sdrf D42/D43).
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestModelEvidenceReader
    {
        private static PublicationInput Input() => new(
            "PXD060431",
            "Human heart proteomes from HFpEF patients and non-failing donors.",
            "## Methods Left ventricular biopsies were collected from 10 HFpEF patients and 10 non-failing donors. Each sample was injected once.",
            new[]
            {
                new SupplementTable("JAH3-s001.xlsx", "Data", "", new[] { "ID", "group", "Age" },
                    new IReadOnlyList<string>[] { new[] { "HumanHFpEF_1", "female", "77" }, new[] { "HumanControl_1", "male", "72" } },
                    new[] { 2, 3 }, false)
            },
            new[] { "HumanHFpEF_1.raw", "HumanControl_1.raw" },
            Array.Empty<SdrfEvidence>());

        private static string Claim(string file, string column, string value, string quote, string source = "supplement", string confidence = "likely", string label = "", string pattern = "") =>
            $$"""{"data_file":"{{file}}","data_file_pattern":"{{pattern}}","label":"{{label}}","column":"{{column}}","value":"{{value}}","source":"{{source}}","quote":"{{quote}}","confidence":"{{confidence}}"}""";

        private static string Answer(string design, params string[] claims) => $$"""{"claims":[{{string.Join(",", claims)}}],"design":{{design}}}""";

        private const string NoDesign = """{"groups":[],"plexes":0,"plexes_quote":"","technical_replicates":0,"technical_replicates_quote":"","fractions_per_sample":0,"fractions_quote":"","runs_stated":0,"runs_quote":""}""";

        private static string Design(string groups = "", int plexes = 0, string plexesQuote = "", int tech = 0, string techQuote = "",
            int fractions = 0, string fractionsQuote = "") =>
            $$"""{"groups":[{{groups}}],"plexes":{{plexes}},"plexes_quote":"{{plexesQuote}}","technical_replicates":{{tech}},"technical_replicates_quote":"{{techQuote}}","fractions_per_sample":{{fractions}},"fractions_quote":"{{fractionsQuote}}","runs_stated":0,"runs_quote":""}""";

        private static string Group(string name, int samples, string quote) => $$"""{"name":"{{name}}","samples":{{samples}},"quote":"{{quote}}"}""";

        [Test]
        public void AClaimQuotingATableRowIsKeptWithTheRowAsItsLocator()
        {
            var json = Answer(NoDesign, Claim("HumanHFpEF_1.raw", "characteristics[sex]", "female", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77"));

            var r = ModelEvidenceReader.Interpret(json, Input());

            var e = r.Evidence.Single();
            Assert.That((e.DataFile, e.Column, e.Value, e.Method, e.Locator), Is.EqualTo(
                ("HumanHFpEF_1.raw", "characteristics[sex]", "female", "model", "JAH3-s001.xlsx!Data!R2")));
            Assert.That(r.Rejected, Is.Empty);
        }

        [Test]
        public void AClaimFromThePaperIsLocatedByItsQuote()
        {
            var json = Answer(NoDesign, Claim("", "characteristics[organism part]", "heart left ventricle",
                "Left ventricular biopsies were collected from 10 HFpEF patients", source: "paper"));

            var e = ModelEvidenceReader.Interpret(json, Input()).Evidence.Single();

            Assert.That(e.DataFile, Is.Empty, "a deposit-wide claim");
            Assert.That(e.Locator, Does.StartWith("paper: \"Left ventricular biopsies"));
        }

        [TestCase("HumanHFpEF_1.raw", "characteristics[age]", "77Y", "HumanHFpEF_1 was 81 years old", "the quote is not in the given text")]
        [TestCase("HumanHFpEF_9.raw", "characteristics[age]", "77Y", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77", "not one of the deposit's raw files")]
        [TestCase("HumanHFpEF_1.raw", "age", "77Y", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77", "not an SDRF column")]
        [TestCase("HumanHFpEF_1.raw", "characteristics[age]", "", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77", "no value")]
        public void AnUnsupportedClaimIsRejectedAndSaysWhy(string file, string column, string value, string quote, string why)
        {
            var r = ModelEvidenceReader.Interpret(Answer(NoDesign, Claim(file, column, value, quote)), Input());

            Assert.That(r.Evidence, Is.Empty);
            Assert.That(r.Rejected.Single(), Does.EndWith(why));
        }

        [Test]
        public void ADesignIsKeptOnlyWithARealQuoteAndItsRunsAreCounted()
        {
            string design = Design(
                Group("HFpEF", 10, "collected from 10 HFpEF patients") + "," + Group("non-failing", 10, "10 non-failing donors"),
                tech: 1, techQuote: "Each sample was injected once.");

            var r = ModelEvidenceReader.Interpret(Answer(design), Input());

            Assert.That(r.Design!.Groups.Select(g => g.Samples), Is.EqualTo(new[] { 10, 10 }));
            Assert.That(r.Design.PredictedRuns, Is.EqualTo(20), "10 + 10 samples, one injection, no fractions");
            Assert.That(r.Design.Quotes, Has.Count.EqualTo(3));

            var invented = design.Replace("collected from 10 HFpEF", "collected from 12 HFpEF");
            var r2 = ModelEvidenceReader.Interpret(Answer(invented), Input());
            Assert.That(r2.Rejected.Single(), Is.EqualTo("design group 'HFpEF' = 10: the quote is not in the given text"));
            Assert.That(r2.Design!.Groups[0].Samples, Is.Zero);
            Assert.That(r2.Design.PredictedRuns, Is.Zero, "one group's size unknown: no count, rather than a false mismatch");
        }

        [Test]
        public void AnIsobaricDesignCountsRunsByPlexNotBySample()
        {
            // PXD010429: 174 samples in 29 TMT 6-plexes, 12 fractions, injected twice = 696 runs, the deposit's count.
            var input = Input() with { PaperText = "The 174 samples were distributed across 29 TMT 6-plexes and separated into twelve concatenated fractions, each injected twice." };
            string design = Design(plexes: 29, plexesQuote: "distributed across 29 TMT 6-plexes",
                tech: 2, techQuote: "each injected twice", fractions: 12, fractionsQuote: "twelve concatenated fractions");

            var d = ModelEvidenceReader.Interpret(Answer(design), input).Design!;

            Assert.That(d.Plexes, Is.EqualTo(29));
            Assert.That(d.PredictedRuns, Is.EqualTo(696), "plexes x fractions x injections, not samples");
        }

        [Test]
        public void ACountTheQuoteDoesNotStateIsDropped()
        {
            // G36 rerun, PXD011967: 20 plexes counted from Set1..Set20 in the file names, quoted from a sentence about one set.
            var input = Input() with { PaperText = "Each TMT6plex set contained one donor from each of five age groups and a reference." };
            string design = Design(plexes: 20, plexesQuote: "Each TMT6plex set contained one donor", fractions: 5, fractionsQuote: "each of five age groups");

            var r = ModelEvidenceReader.Interpret(Answer(design), input);

            Assert.That(r.Design, Is.Null, "nothing countable is left");
            Assert.That(r.Rejected, Does.Contain("design plexes = 20: the quote does not state it"));
            Assert.That(r.Rejected, Does.Contain("design fractions = 5: the quote does not say 'fraction' near the number"),
                "'five' states 5, but of age groups, not of fractions");
        }

        [TestCase("each injected twice", 2, true)]
        [TestCase("analysed in triplicate", 3, true)]
        [TestCase("twelve concatenated fractions", 12, true)]
        [TestCase("29 TMT 6-plexes", 6, true)]
        [TestCase("29 TMT 6-plexes", 9, false)]
        [TestCase("1.5 mg of protein", 5, false)]
        [TestCase("a 120 min gradient", 12, false)]
        [TestCase("n = 10 per group", 10, true)]
        public void AQuoteStatesANumberAsDigitsOrWords(string quote, int n, bool states) =>
            Assert.That(ModelEvidenceReader.States(quote, n), Is.EqualTo(states));

        [Test]
        public void AClaimAboutASetOfFilesKeepsItsPatternAndAPatternMatchingNothingIsRejected()
        {
            const string q = "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77";
            var r = ModelEvidenceReader.Interpret(Answer(NoDesign,
                Claim("", "characteristics[disease]", "HFpEF", q, pattern: "HumanHFpEF_*"),
                Claim("", "characteristics[disease]", "x", q, pattern: "*_TMT6_*"),
                Claim("HumanHFpEF_1.raw", "characteristics[sex]", "female", q, pattern: "Human*")), Input());

            Assert.That(r.Evidence.Single().DataFilePattern, Is.EqualTo("HumanHFpEF_*"));
            Assert.That(r.Rejected, Has.Some.EndsWith("the file pattern matches none of the deposit's raw files"));
            Assert.That(r.Rejected, Has.Some.EndsWith("names both a file and a file pattern"));
        }

        [Test]
        public void AGuessStaysAGuessAndDuplicatesCollapse()
        {
            var q = "[JAH3-s001.xlsx!Data!R3] HumanControl_1 | male | 72";
            var r = ModelEvidenceReader.Interpret(Answer(NoDesign,
                Claim("HumanControl_1.raw", "characteristics[sex]", "male", q, confidence: "guess"),
                Claim("HumanControl_1.raw", "characteristics[sex]", "male", q)), Input());

            Assert.That(r.Evidence.Single().Confidence, Is.EqualTo(SdrfEvidenceConfidence.Guess));
        }

        [Test]
        public void TwoClaimsThatDisagreeBothReachTheDrafter()
        {
            // pcruzparri on #1377: the key had no value, so the second, disagreeing claim was dropped and nothing said so.
            var r = ModelEvidenceReader.Interpret(Answer(NoDesign,
                Claim("HumanControl_1.raw", "characteristics[sex]", "male", "[JAH3-s001.xlsx!Data!R3] HumanControl_1 | male | 72", confidence: "guess"),
                Claim("HumanControl_1.raw", "characteristics[sex]", "female", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77"),
                Claim("HumanControl_1.raw", "characteristics[sex]", "MALE", "[JAH3-s001.xlsx!Data!R3] HumanControl_1 | male | 72")), Input());

            Assert.That(r.Evidence.Select(e => e.Value), Is.EqualTo(new[] { "male", "female" }), "a case-insensitive repeat still collapses");
            Assert.That(r.Rejected, Is.Empty);
        }

        [TestCase("and", "the quote is too short to check")]
        [TestCase("Human", "the quote is too short to check")]
        public void ATrivialQuoteIsRejected(string quote, string why)
        {
            var r = ModelEvidenceReader.Interpret(Answer(NoDesign, Claim("HumanHFpEF_1.raw", "characteristics[disease]", "type 2 diabetes mellitus", quote)), Input());

            Assert.That(r.Evidence, Is.Empty);
            Assert.That(r.Rejected.Single(), Does.EndWith(why));
        }

        [Test]
        public void AQuoteMustBeInTheRowItCites()
        {
            // R3 says male; the quoted text is R2's. The locator would have sent a checker to the wrong row.
            var r = ModelEvidenceReader.Interpret(Answer(NoDesign, Claim("HumanControl_1.raw", "characteristics[sex]", "female",
                "[JAH3-s001.xlsx!Data!R3] HumanHFpEF_1 | female")), Input());

            Assert.That(r.Evidence, Is.Empty);
            Assert.That(r.Rejected.Single(), Does.EndWith("the quote is not whole cells of the row it cites"));
        }

        [TestCase("characteristics[sex]", "female", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77", true)]
        [TestCase("characteristics[sex]", "male", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77", false)]
        [TestCase("characteristics[age]", "77Y", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77", true)]
        [TestCase("characteristics[age]", "35Y", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77", false)]
        [TestCase("characteristics[organism part]", "heart left ventricle", "Left ventricular biopsies were collected", true)]
        [TestCase("characteristics[disease]", "type 2 diabetes mellitus", "Left ventricular biopsies were collected", false)]
        public void AClaimWhoseValueItsQuoteDoesNotStateIsOnlyAGuess(string column, string value, string quote, bool stated)
        {
            var r = ModelEvidenceReader.Interpret(Answer(NoDesign, Claim("HumanHFpEF_1.raw", column, value, quote)), Input());

            Assert.That(r.Evidence.Single().Confidence, Is.EqualTo(stated ? SdrfEvidenceConfidence.Likely : SdrfEvidenceConfidence.Guess));
        }

        [Test]
        public void AFilePatternWithManyWildcardsIsRejectedWithoutBeingTried()
        {
            // pcruzparri on #1377: '*0' x 12 + '*ZZZ.raw' took 54 s against one 48-character name in a backtracking regex.
            string pattern = string.Concat(Enumerable.Repeat("*0", 12)) + "*ZZZ.raw";
            var input = Input() with { RawFiles = new[] { new string('0', 44) + ".raw", "HumanHFpEF_1.raw" } };
            var clock = System.Diagnostics.Stopwatch.StartNew();

            var r = ModelEvidenceReader.Interpret(Answer(NoDesign, Claim("", "characteristics[disease]", "HFpEF",
                "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77", pattern: pattern)), input);

            Assert.That(clock.Elapsed, Is.LessThan(TimeSpan.FromSeconds(1)));
            Assert.That(r.Rejected.Single(), Does.EndWith("the file pattern has more than 4 '*' wildcards"));
        }

        [Test]
        public void AnAnswerThatIsNotJsonIsRejectedNotThrown()
        {
            var r = ModelEvidenceReader.Interpret("I could not find anything.", Input());
            Assert.That(r.Evidence, Is.Empty);
            Assert.That(r.Rejected.Single(), Does.StartWith("not JSON"));
        }

        [Test]
        public void ThePromptCarriesTheFilesTheTablesAndWhatTheRulesFound()
        {
            var input = Input() with { RuleEvidence = new[] { new SdrfEvidence("HumanHFpEF_1.raw", "", "characteristics[age]", "77Y", "supplement", "x", "file-key", SdrfEvidenceConfidence.Likely) } };

            string prompt = ModelEvidenceReader.UserPrompt(input);

            Assert.That(prompt, Does.Contain("HumanControl_1.raw"));
            Assert.That(prompt, Does.Contain("[JAH3-s001.xlsx!Data!R3] HumanControl_1 | male | 72"));
            Assert.That(prompt, Does.Contain("HumanHFpEF_1.raw\t\tcharacteristics[age]\t77Y"));
        }

        private static SupplementTable Table(string file, string[] header, int rows, Func<int, string[]> row) => new(
            file, "S1", "", header, Enumerable.Range(1, rows).Select(k => (IReadOnlyList<string>)row(k)).ToArray(),
            Enumerable.Range(2, rows).ToArray(), false);

        private static SupplementTable ProteinTable() => Table("proteins.xlsx", new[] { "Protein Accession", "Gene", "log2 fold change", "p-value" }, 300,
            k => new[] { $"P{10000 + k}", $"GENE{k}", "1.5", "0.01" });

        private static SupplementTable CohortTable() => Table("cohort.xlsx", new[] { "Sample ID", "Age", "Sex", "Diagnosis" }, 200,
            k => new[] { $"S{k:000}", $"{40 + k % 30}", k % 2 == 0 ? "F" : "M", "HFpEF" });

        [Test]
        public void ASampleTableComesBeforeAResultTableAndIsGivenWhole()
        {
            // G36 batch 2: in document order, protein tables filled the budget before the cohort table was reached.
            string text = PublicationText.Tables(new[] { ProteinTable(), CohortTable() }, maxRowsPerTable: 150, maxChars: 20_000);

            Assert.That(text.IndexOf("### cohort.xlsx"), Is.LessThan(text.IndexOf("### proteins.xlsx")), "sample table first");
            Assert.That(text, Does.Contain("[cohort.xlsx!S1!R201] S200 | 60 | F | HFpEF"), "all 200 rows, past the 150-row cap");
            Assert.That(text, Does.Contain("### proteins.xlsx!S1 [result table: preview only]"));
            Assert.That(text, Does.Contain("[proteins.xlsx!S1!R6] P10005"));
            Assert.That(text, Does.Not.Contain("[proteins.xlsx!S1!R7]"), "a result table is previewed, not given");
        }

        [Test]
        public void RepeatedHeadersAndFigureDataDoNotOutrankTheCohort()
        {
            // G36 batch 2, PXD017291: three "fraction_*" headers made a figure's 3,800 rows the first table in the prompt.
            var figure = Table("s003.xlsx", new[] { "gene_name", "set", "fraction_helix", "fraction_sheet", "fraction_coil" }, 3000,
                k => new[] { $"G{k}", "aggregator", "0.1", "0.2", "0.7" }) with { Sheet = "Figure2G" };
            var clusters = Table("s007.xlsx", new[] { "", "Cluster 1", "Cluster 2", "Cluster 3" }, 3000, k => new[] { $"c{k}", "1", "2", "3" });
            string text = PublicationText.Tables(new[] { figure, clusters, CohortTable() }, maxChars: 30_000);

            Assert.That(text, Does.StartWith("### cohort.xlsx"));
            Assert.That(text, Does.Contain("### s003.xlsx!Figure2G"), "a long table cannot hide the ones after it");
            Assert.That(text, Does.Contain("### s007.xlsx"));

            var big = Table("big.xlsx", new[] { "Sample ID", "Age", "Sex" }, 2000, k => new[] { $"B{k:0000}", "50", "F" });
            string two = PublicationText.Tables(new[] { big, CohortTable() }, maxChars: 30_000);
            Assert.That(two, Does.Contain("### cohort.xlsx"), "a long table cannot take the whole budget");
        }

        [Test]
        public void ARowLeftOutOfThePromptIsStillValidText()
        {
            // The quote check reads every row: shortening the prompt must not turn a true quote into a rejection.
            var input = Input() with { Tables = new[] { ProteinTable() } };
            Assert.That(ModelEvidenceReader.UserPrompt(input), Does.Not.Contain("P10200"));

            var r = ModelEvidenceReader.Interpret(Answer(NoDesign, Claim("", "characteristics[organism]", "Homo sapiens",
                "[proteins.xlsx!S1!R201] P10200 | GENE200 | 1.5 | 0.01", confidence: "guess")), input);

            Assert.That(r.Rejected, Is.Empty);
            Assert.That(r.Evidence.Single().Confidence, Is.EqualTo(SdrfEvidenceConfidence.Guess), "the row is real, but it does not state the organism");
        }

        [Test]
        public void ManyFilesAreListedByPatternAndNoneIsLost()
        {
            // G36 batch 2: a 600-name cap cut PXD011967's last 50 files.
            var files = Enumerable.Range(1, 70).SelectMany(s => Enumerable.Range(1, 10).Select(f => $"Plasma_Sample{s:000}_F{f:00}.raw")).ToList();
            files.Add("QC_blank.raw");

            string text = PublicationText.Files(files, maxChars: 5_000);

            Assert.That(text, Does.Contain("Plasma_Sample{n}_F{n}.raw  (700 files)"));
            Assert.That(text, Does.Contain("001/01, 001/02"));
            Assert.That(text, Does.Contain("070/10"), "the last file is there");
            Assert.That(text, Does.Contain("QC_blank.raw"));
            Assert.That(text.Length, Is.LessThan(files.Sum(f => f.Length + 1)));

            Assert.That(PublicationText.Files(new[] { "a1.raw", "a2.raw" }), Is.EqualTo("a1.raw\na2.raw\n"), "a short list is given whole");
        }

        [Test]
        public void ASilacStateIsAChannelLabel()
        {
            // G36 batch 2, PXD006430: with no way to say which sample was heavy, both went into one source name.
            var input = Input() with { PaperText = "Proteins from ctrl cells (Light SILAC labeled) and from TRAP3high cells (Heavy SILAC labeled) were mixed." };
            var r = ModelEvidenceReader.Interpret(Answer(NoDesign,
                Claim("HumanHFpEF_1.raw", "source name", "TRAP3high", "from TRAP3high cells (Heavy SILAC labeled)", source: "paper", label: "SILAC heavy"),
                Claim("HumanHFpEF_1.raw", "source name", "ctrl", "from ctrl cells (Light SILAC labeled)", source: "paper", label: "SILAC light"),
                Claim("HumanHFpEF_1.raw", "source name", "x", "from ctrl cells (Light SILAC labeled)", source: "paper", label: "heavy")), input);

            Assert.That(r.Evidence.Select(e => e.Label), Is.EqualTo(new[] { "SILAC heavy", "SILAC light" }));
            Assert.That(r.Rejected.Single(), Does.EndWith("not a channel label"));
        }

        private sealed class RecordingHandler : HttpMessageHandler
        {
            private readonly string _answer;
            private readonly string _stopReason;
            public List<string> Bodies { get; } = new();
            public RecordingHandler(string answer, string stopReason = "end_turn") => (_answer, _stopReason) = (answer, stopReason);

            protected override async Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
            {
                string body = await request.Content!.ReadAsStringAsync(cancellationToken);
                Bodies.Add(body);
                var usage = new { input_tokens = 1200, output_tokens = 300, cache_read_input_tokens = 900, cache_creation_input_tokens = 0 };
                if (!System.Text.Json.JsonDocument.Parse(body).RootElement.TryGetProperty("stream", out var s) || !s.GetBoolean())
                {
                    string message = System.Text.Json.JsonSerializer.Serialize(new
                    {
                        id = "msg_test", type = "message", role = "assistant", model = "claude-opus-5",
                        content = new[] { new { type = "text", text = _answer } },
                        stop_reason = _stopReason, stop_sequence = (string?)null, usage
                    });
                    return new HttpResponseMessage(HttpStatusCode.OK) { Content = new StringContent(message, Encoding.UTF8, "application/json") };
                }
                // The streamed form: the same message as server-sent events.
                static string Event(string name, object data) => $"event: {name}\ndata: {System.Text.Json.JsonSerializer.Serialize(data)}\n\n";
                var sse = new StringBuilder()
                    .Append(Event("message_start", new { type = "message_start", message = new { id = "msg_test", type = "message", role = "assistant", model = "claude-opus-5", content = Array.Empty<object>(), stop_reason = (string?)null, stop_sequence = (string?)null, usage = new { input_tokens = 1200, output_tokens = 1, cache_read_input_tokens = 900, cache_creation_input_tokens = 0 } } }))
                    .Append(Event("content_block_start", new { type = "content_block_start", index = 0, content_block = new { type = "text", text = "" } }))
                    .Append(Event("content_block_delta", new { type = "content_block_delta", index = 0, delta = new { type = "text_delta", text = _answer } }))
                    .Append(Event("content_block_stop", new { type = "content_block_stop", index = 0 }))
                    .Append(Event("message_delta", new { type = "message_delta", delta = new { stop_reason = _stopReason, stop_sequence = (string?)null }, usage }))
                    .Append(Event("message_stop", new { type = "message_stop" }));
                return new HttpResponseMessage(HttpStatusCode.OK) { Content = new StringContent(sse.ToString(), Encoding.UTF8, "text/event-stream") };
            }
        }

        [Test]
        public async Task TheRequestAsksForTheSchemaCachesTheInstructionsAndFallsBackOnARefusal()
        {
            string answer = Answer(NoDesign, Claim("HumanHFpEF_1.raw", "characteristics[sex]", "female", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77"));
            var handler = new RecordingHandler(answer);
            var client = new AnthropicClient { ApiKey = "test-key", HttpClient = new HttpClient(handler), MaxRetries = 0 };

            var result = await new ModelEvidenceReader(client).ReadAsync(Input());

            using var request = System.Text.Json.JsonDocument.Parse(handler.Bodies.Single());
            var root = request.RootElement;
            Assert.That(root.GetProperty("model").GetString(), Is.EqualTo("claude-opus-5"));
            Assert.That(root.GetProperty("fallbacks").GetString(), Is.EqualTo("default"), "the server-side refusal fallback");
            Assert.That(root.GetProperty("output_config").GetProperty("format").GetProperty("type").GetString(), Is.EqualTo("json_schema"));
            Assert.That(root.GetProperty("system")[0].GetProperty("cache_control").GetProperty("type").GetString(), Is.EqualTo("ephemeral"));
            Assert.That(root.GetProperty("thinking").GetProperty("type").GetString(), Is.EqualTo("adaptive"));

            Assert.That(result.Evidence.Single().Value, Is.EqualTo("female"));
            Assert.That((result.InputTokens, result.OutputTokens, result.CacheReadTokens), Is.EqualTo((1200L, 300L, 900L)));
        }

        [Test]
        public async Task ARefusalYieldsNoClaimsEvenWithTextThatWouldPass()
        {
            // A refusal can still carry text; a claim that would pass validation must not be kept from it.
            string answer = Answer(NoDesign, Claim("HumanHFpEF_1.raw", "characteristics[sex]", "female", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77"));
            var client = new AnthropicClient { ApiKey = "test-key", HttpClient = new HttpClient(new RecordingHandler(answer, "refusal")), MaxRetries = 0 };

            var result = await new ModelEvidenceReader(client).ReadAsync(Input());

            Assert.That(result.Evidence, Is.Empty);
            Assert.That(result.StopReason, Is.EqualTo("refusal"));
            Assert.That(result.Rejected.Single(), Is.EqualTo("no reading: stop reason 'refusal'"));
            Assert.That(result.InputTokens, Is.EqualTo(1200L), "a refused call still reports what it cost");
        }

        // ---------------- Alexander-Sol's review of #1377 (2026-10-01) ----------------

        [TestCase("characteristics[sex]", "male", "[JAH3-s001.xlsx!Data!R2] male")]
        [TestCase("characteristics[age]", "7Y", "[JAH3-s001.xlsx!Data!R2] 7")]
        public void ARowQuoteMustBeWholeCellsOfThatRow(string column, string value, string quote)
        {
            // R2 is HumanHFpEF_1 | female | 77: "male" is inside "female", "7" inside "77", and neither is a cell.
            var r = ModelEvidenceReader.Interpret(Answer(NoDesign, Claim("HumanHFpEF_1.raw", column, value, quote)), Input());

            Assert.That(r.Evidence, Is.Empty);
            Assert.That(r.Rejected.Single(), Does.EndWith("the quote is not whole cells of the row it cites"));
        }

        [Test]
        public void ARowAboutAnotherSampleOnlyMakesAGuess()
        {
            // R2 is HumanHFpEF_1's row; nothing in it names HumanControl_1.
            var r = ModelEvidenceReader.Interpret(Answer(NoDesign,
                Claim("HumanControl_1.raw", "characteristics[sex]", "female", "[JAH3-s001.xlsx!Data!R2] female"),
                Claim("HumanHFpEF_1.raw", "characteristics[sex]", "female", "[JAH3-s001.xlsx!Data!R2] female")), Input());

            Assert.That(r.Evidence.Select(e => (e.DataFile, e.Confidence)), Is.EqualTo(new[]
            {
                ("HumanControl_1.raw", SdrfEvidenceConfidence.Guess),
                ("HumanHFpEF_1.raw", SdrfEvidenceConfidence.Likely)
            }), "a row names its own file, not another's");
        }

        [Test]
        public void ADesignNumberNeedsAFullQuoteThatSaysWhatItCounts()
        {
            var input = Input() with { PaperText = "Peptides from twelve donors were pooled and separated into twelve fractions." };

            var r = ModelEvidenceReader.Interpret(Answer(Design(fractions: 12, fractionsQuote: "twelve")), input);
            Assert.That(r.Rejected, Does.Contain("design fractions = 12: the quote is too short to check"));

            var r2 = ModelEvidenceReader.Interpret(Answer(Design(fractions: 12, fractionsQuote: "Peptides from twelve donors were pooled")), input);
            Assert.That(r2.Rejected, Does.Contain("design fractions = 12: the quote does not say 'fraction' near the number"));

            var r3 = ModelEvidenceReader.Interpret(Answer(Design(fractions: 12, fractionsQuote: "separated into twelve fractions")), input);
            Assert.That(r3.Rejected, Is.Empty);
        }

        [TestCase("characteristics[disease]", "normal", "intensities were normalized to the median", false)]
        [TestCase("characteristics[disease]", "normal", "normal tissue adjacent to the tumour", true)]
        [TestCase("factor value[treatment]", "control", "kept in a temperature controlled room", false)]
        [TestCase("characteristics[disease]", "type 2 diabetes mellitus", "the cell type was not recorded", false)]
        [TestCase("characteristics[disease]", "type 2 diabetes mellitus", "patients with type 2 diabetes", true)]
        [TestCase("characteristics[organism part]", "heart left ventricle", "left ventricular biopsies were collected", true)]
        [TestCase("characteristics[age]", "35Y", "separated on a 35 min gradient", false)]
        [TestCase("characteristics[age]", "35Y", "a 35-year-old donor", true)]
        [TestCase("characteristics[age]", "35Y", "patients aged 35 at biopsy", true)]
        [TestCase("characteristics[age]", "9M", "mice were 9 years old", false)]
        [TestCase("characteristics[age]", "9M", "mice at 9 months of age", true)]
        [TestCase("characteristics[age]", "77Y", "humanhfpef_1 | female | 77", true)]
        [TestCase("characteristics[biological replicate]", "2", "the second of 2 replicates", true)]
        public void AQuoteStatesAValueAsAWordOrAnAgeNotAsAFragment(string column, string value, string quote, bool states) =>
            Assert.That(ModelEvidenceReader.QuoteStates(quote, column, value), Is.EqualTo(states));

        [Test]
        public void LigaturesAndSoftHyphensInAPdfDoNotRejectATrueQuote()
        {
            var input = Input() with { PaperText = "Left ven­tricular biopsies were ﬁxed in formalin." };

            var r = ModelEvidenceReader.Interpret(Answer(NoDesign, Claim("", "characteristics[organism part]", "heart left ventricle",
                "Left ventricular biopsies were fixed in formalin", source: "paper")), input);

            Assert.That(r.Rejected, Is.Empty);
            Assert.That(r.Evidence, Has.Count.EqualTo(1));
        }

        [Test]
        public void APaperLongerThanTheCapIsCutInThePromptButStillChecked()
        {
            const string tail = "Each sample was injected in triplicate at the very end.";
            var input = Input() with { PaperText = new string('x', ModelEvidenceReader.MaxPaperChars) + " collected from 10 HFpEF patients. " + tail };

            string prompt = ModelEvidenceReader.UserPrompt(input);

            Assert.That(prompt.Length, Is.LessThan(ModelEvidenceReader.MaxPaperChars + 10_000));
            Assert.That(prompt, Does.Contain("[paper cut at"));
            Assert.That(prompt, Does.Not.Contain(tail));
            var d = ModelEvidenceReader.Interpret(Answer(Design(Group("HFpEF", 10, "collected from 10 HFpEF patients"), tech: 3, techQuote: tail)), input).Design!;
            Assert.That(d.TechnicalReplicates, Is.EqualTo(3), "the quote check reads the whole paper");
        }

        [Test]
        public async Task AnAnswerCutOffAtTheOutputLimitSaysSoAndKeepsNothing()
        {
            // A large deposit (PXD007160: 1,688 channel claims) can outrun max_tokens; the JSON then ends mid-claim.
            string full = Answer(NoDesign, Claim("HumanHFpEF_1.raw", "characteristics[sex]", "female", "[JAH3-s001.xlsx!Data!R2] HumanHFpEF_1 | female | 77"));
            var handler = new RecordingHandler(full[..(full.Length / 2)], "max_tokens");
            var client = new AnthropicClient { ApiKey = "test-key", HttpClient = new HttpClient(handler), MaxRetries = 0 };

            var result = await new ModelEvidenceReader(client).ReadAsync(Input());

            Assert.That(result.StopReason, Is.EqualTo("max_tokens"));
            Assert.That(result.Evidence, Is.Empty);
            Assert.That(result.Rejected.Single(), Does.StartWith("no reading: the answer was cut off at the output limit"));
            Assert.That(result.OutputTokens, Is.EqualTo(300L), "a cut-off call still reports what it cost");
            using var request = System.Text.Json.JsonDocument.Parse(handler.Bodies.Single());
            Assert.That(request.RootElement.GetProperty("stream").GetBoolean(), Is.True, "streamed, so a long answer cannot time out");
            Assert.That(request.RootElement.GetProperty("max_tokens").GetInt32(), Is.EqualTo(ModelEvidenceReader.MaxOutputTokens));
        }

        /// <summary>
        /// The one live call: catches drift in what the recorded stub cannot (the model id, the beta header, the fallback
        /// and adaptive thinking). Runs in the external-service job when ANTHROPIC_API_KEY is set; skipped otherwise.
        /// </summary>
        [Test]
        [Category("ExternalService")]
        public async Task LiveTheApiAcceptsTheRequestThisReaderSends()
        {
            string? key = Environment.GetEnvironmentVariable("ANTHROPIC_API_KEY");
            if (string.IsNullOrWhiteSpace(key)) Assert.Ignore("ANTHROPIC_API_KEY is not set; the live model check is skipped.");

            var result = await new ModelEvidenceReader(new AnthropicClient { ApiKey = key }).ReadAsync(Input());

            Assert.That(result.StopReason, Is.EqualTo("end_turn"), string.Join("; ", result.Rejected));
            Assert.That(result.Rejected, Has.None.StartsWith("not JSON"));
            Assert.That(result.InputTokens, Is.GreaterThan(0));
        }


        [Test]
        public void TheMethodsAndTablesOfAJatsArticleAreKeptAndTheReferencesAreNot()
        {
            const string jats = """
                <article><front><aff>Leipzig, Germany</aff></front><body>
                <sec><title>Introduction</title><p>Heart failure is common.</p></sec>
                <sec sec-type="methods"><title>Materials and methods</title><p>Biopsies from 10 patients.</p></sec>
                <table-wrap><caption>Table 1 Cohort</caption><table><tr><td>Age</td><td>77</td></tr></table></table-wrap>
                </body><back><ref-list><ref>Smith et al. Germany</ref></ref-list></back></article>
                """;

            string text = PublicationText.MethodsAndTables(jats);

            Assert.That(text, Does.Contain("Biopsies from 10 patients"));
            Assert.That(text, Does.Contain("Table 1 Cohort"));
            Assert.That(text, Does.Not.Contain("Heart failure is common"));
            Assert.That(text, Does.Not.Contain("Germany"));
        }
    }
}
