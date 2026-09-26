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

        private static string Claim(string file, string column, string value, string quote, string source = "supplement", string confidence = "likely", string label = "") =>
            $$"""{"data_file":"{{file}}","label":"{{label}}","column":"{{column}}","value":"{{value}}","source":"{{source}}","quote":"{{quote}}","confidence":"{{confidence}}"}""";

        private static string Answer(string design, params string[] claims) => $$"""{"claims":[{{string.Join(",", claims)}}],"design":{{design}}}""";

        private const string NoDesign = """{"groups":[],"plexes":0,"technical_replicates":0,"fractions_per_sample":0,"runs_stated":0,"quote":""}""";

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
            const string design = """{"groups":[{"name":"HFpEF","samples":10},{"name":"non-failing","samples":10}],"plexes":0,"technical_replicates":1,"fractions_per_sample":0,"runs_stated":0,"quote":"collected from 10 HFpEF patients and 10 non-failing donors"}""";

            var r = ModelEvidenceReader.Interpret(Answer(design), Input());

            Assert.That(r.Design!.Groups.Select(g => g.Samples), Is.EqualTo(new[] { 10, 10 }));
            Assert.That(r.Design.PredictedRuns, Is.EqualTo(20), "10 + 10 samples, one injection, no fractions");

            var invented = design.Replace("collected from 10 HFpEF", "collected from 12 HFpEF");
            var r2 = ModelEvidenceReader.Interpret(Answer(invented), Input());
            Assert.That(r2.Design, Is.Null);
            Assert.That(r2.Rejected.Single(), Does.StartWith("design"));
        }

        [Test]
        public void AnIsobaricDesignCountsRunsByPlexNotBySample()
        {
            // PXD010429: 174 samples in 29 TMT 6-plexes, 12 fractions, injected twice = 696 runs, the deposit's count.
            var input = Input() with { PaperText = "The 174 samples were distributed across 29 TMT 6-plexes and separated into twelve concatenated fractions, each injected twice." };
            const string design = """{"groups":[{"name":"tumor","samples":116},{"name":"pool","samples":58}],"plexes":29,"technical_replicates":2,"fractions_per_sample":12,"runs_stated":0,"quote":"distributed across 29 TMT 6-plexes"}""";

            var d = ModelEvidenceReader.Interpret(Answer(design), input).Design!;

            Assert.That(d.Plexes, Is.EqualTo(29));
            Assert.That(d.PredictedRuns, Is.EqualTo(696), "plexes x fractions x injections, not samples");
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
            public List<string> Bodies { get; } = new();
            public RecordingHandler(string answer) => _answer = answer;

            protected override async Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
            {
                Bodies.Add(await request.Content!.ReadAsStringAsync(cancellationToken));
                string message = System.Text.Json.JsonSerializer.Serialize(new
                {
                    id = "msg_test", type = "message", role = "assistant", model = "claude-opus-5",
                    content = new[] { new { type = "text", text = _answer } },
                    stop_reason = "end_turn", stop_sequence = (string?)null,
                    usage = new { input_tokens = 1200, output_tokens = 300, cache_read_input_tokens = 900, cache_creation_input_tokens = 0 }
                });
                return new HttpResponseMessage(HttpStatusCode.OK) { Content = new StringContent(message, Encoding.UTF8, "application/json") };
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
