using MzLibUtil;
using Newtonsoft.Json.Linq;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

/// <summary>
/// M1b: MassIVE's own search (<see cref="MassiveClient.SearchDatasetsAsync"/>) and ProteomeCentral's search across
/// every ProteomeXchange repository (<see cref="ProteomeCentralClient.SearchDatasetsAsync"/>), offline over pages
/// both services served on 2026-11-10 (MassiveFixtures/PROVENANCE.md). Backoffs are counted, not slept through.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class MassiveSearchTests
{
    /// <summary>Answers each request with <c>respond</c> and records the URIs asked for.</summary>
    private sealed class Handler(Func<HttpRequestMessage, HttpResponseMessage> respond) : HttpMessageHandler
    {
        public List<Uri> Requests { get; } = new();

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            Requests.Add(request.RequestUri!);
            return Task.FromResult(respond(request));
        }
    }

    private int _backoffs;

    [SetUp]
    public void SetUp() => _backoffs = 0;

    private static byte[] Fixture(string name) =>
        File.ReadAllBytes(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "MassiveFixtures", name));

    private static HttpResponseMessage Body(byte[] body, HttpStatusCode status = HttpStatusCode.OK) =>
        new(status) { Content = new ByteArrayContent(body) };

    private static HttpResponseMessage Body(string body, HttpStatusCode status = HttpStatusCode.OK) =>
        Body(Encoding.UTF8.GetBytes(body), status);

    private static int PageOf(HttpRequestMessage request) =>
        int.Parse(Regex.Match(request.RequestUri!.Query, @"pageNumber=(\d+)").Groups[1].Value);

    private static readonly byte[] NotFound = Encoding.ASCII.GetBytes("""{"code":404,"title":"Not Found","message":"No data found for the specified parameters."}""");

    /// <summary>MassIVE's four captured pages for "zebrafish" at pageSize 30: 30, 30, 7 rows, then 404.</summary>
    private static HttpResponseMessage MassivePage(int page) => page switch
    {
        <= 3 => Body(Fixture($"massive_search_zebrafish_p{page}.200.json")),
        _ => Body(NotFound, HttpStatusCode.NotFound),
    };

    /// <summary>ProteomeCentral's four captured pages for "zebrafish" in MassIVE at pageSize 5: 5, 5, 1, 0 rows.</summary>
    private static HttpResponseMessage ProteomeCentralPage(int page) =>
        Body(Fixture($"proteomecentral_search_zebrafish_massive_p{Math.Min(page, 4)}.200.json"));

    private MassiveClient Massive(Handler handler, int maxRetries = 3, int maxPages = 10000) =>
        new(new HttpClient(handler)) { MaxRetries = maxRetries, MaxPages = maxPages, RetryDelay = (_, _) => { _backoffs++; return Task.CompletedTask; } };

    private ProteomeCentralClient Central(Handler handler, int maxRetries = 3, int maxPages = 10000) =>
        new(new HttpClient(handler)) { MaxRetries = maxRetries, MaxPages = maxPages, RetryDelay = (_, _) => { _backoffs++; return Task.CompletedTask; } };

    // ---- MassIVE -----------------------------------------------------------------------------------------------

    [Test]
    public async Task Massive_EveryPageIsRead_UntilAShortOne()
    {
        var handler = new Handler(r => MassivePage(PageOf(r)));
        using MassiveClient client = Massive(handler);

        IReadOnlyList<MassiveDatasetSearchResult> hits = await client.SearchDatasetsAsync("zebrafish", pageSize: 30);

        Assert.Multiple(() =>
        {
            Assert.That(hits, Has.Count.EqualTo(67));
            Assert.That(hits.Select(h => h.MassiveAccession).Distinct().Count(), Is.EqualTo(67));
            Assert.That(hits[0].MassiveAccession, Is.EqualTo("MSV000078304"));
            Assert.That(hits[^1].MassiveAccession, Is.EqualTo("MSV000101937"));
            Assert.That(hits.Count(h => h.PxdAccession != null), Is.EqualTo(19));
            Assert.That(hits.Count(h => h.PxdAccession == null), Is.EqualTo(48), "MassIVE-only datasets ProteomeCentral cannot see");
            Assert.That(hits[0].Title, Does.Contain("zebrafish"));
            Assert.That(handler.Requests, Has.Count.EqualTo(3), "page 3 was short, so page 4 is never asked for");
            Assert.That(handler.Requests[0].AbsoluteUri, Is.EqualTo(
                MassiveClient.DefaultProxiBaseAddress + "datasets?resultType=compact&pageSize=30&pageNumber=1&search=zebrafish"));
        });
    }

    [Test]
    public async Task Massive_404AfterAFullPage_EndsTheSearch()
    {
        var handler = new Handler(r => PageOf(r) <= 2 ? MassivePage(PageOf(r)) : Body(NotFound, HttpStatusCode.NotFound));
        using MassiveClient client = Massive(handler);

        IReadOnlyList<MassiveDatasetSearchResult> hits = await client.SearchDatasetsAsync("zebrafish", pageSize: 30);

        Assert.That(hits, Has.Count.EqualTo(60));
        Assert.That(handler.Requests, Has.Count.EqualTo(3));
    }

    [Test]
    public async Task Massive_NoMatch_IsAnEmptyList_InOneRequest()
    {
        var handler = new Handler(_ => Body(NotFound, HttpStatusCode.NotFound));
        using MassiveClient client = Massive(handler);

        IReadOnlyList<MassiveDatasetSearchResult> hits = await client.SearchDatasetsAsync("qqzzxxvvnomatch");

        Assert.That(hits, Is.Not.Null.And.Empty);
        Assert.That(handler.Requests, Has.Count.EqualTo(1));
    }

    [Test]
    public async Task Massive_ADatasetServedTwiceAsTheIndexMoves_IsKeptOnce()
    {
        // Page 2 starts with page 1's last row, as a row inserted before it would push it along.
        JArray p1 = JArray.Parse(Encoding.UTF8.GetString(Fixture("massive_search_zebrafish_p1.200.json")));
        JArray p2 = JArray.Parse(Encoding.UTF8.GetString(Fixture("massive_search_zebrafish_p2.200.json")));
        JArray p3 = JArray.Parse(Encoding.UTF8.GetString(Fixture("massive_search_zebrafish_p3.200.json")));
        var shifted2 = new JArray(new[] { p1[^1] }.Concat(p2.Take(29)));
        var shifted3 = new JArray(new[] { p2[^1] }.Concat(p3));
        var handler = new Handler(r => PageOf(r) switch
        {
            1 => Body(p1.ToString()),
            2 => Body(shifted2.ToString()),
            _ => Body(shifted3.ToString()),
        });
        using MassiveClient client = Massive(handler);

        IReadOnlyList<MassiveDatasetSearchResult> hits = await client.SearchDatasetsAsync("zebrafish", pageSize: 30);

        Assert.That(hits, Has.Count.EqualTo(67));
        Assert.That(hits.Select(h => h.MassiveAccession).Distinct().Count(), Is.EqualTo(67));
    }

    [Test]
    public void Massive_ServerIgnoringThePaging_IsCaught()
    {
        var handler = new Handler(_ => MassivePage(1));
        using MassiveClient client = Massive(handler);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.SearchDatasetsAsync("zebrafish", pageSize: 30));

        Assert.That(e!.Message, Does.Contain("ignoring the paging"));
        Assert.That(handler.Requests, Has.Count.EqualTo(2));
    }

    [Test]
    public void Massive_PastMaxPages_Throws()
    {
        // Every page full and new, as a server that never stops would send.
        var handler = new Handler(r => Body($$"""[{"accession":[{"accession":"MS:1002487","value":"MSV{{PageOf(r):D9}}"}]}]"""));
        using MassiveClient client = Massive(handler, maxPages: 3);

        Assert.ThrowsAsync<HttpRequestException>(async () => await client.SearchDatasetsAsync("zebrafish", pageSize: 1));
        Assert.That(handler.Requests, Has.Count.EqualTo(3));
    }

    [Test]
    public async Task Massive_SpeciesAndEscaping_GoIntoTheQuery()
    {
        var handler = new Handler(_ => Body(NotFound, HttpStatusCode.NotFound));
        using MassiveClient client = Massive(handler);

        await client.SearchDatasetsAsync("a&pageSize=1#b", speciesTaxId: "7955", pageSize: 10);

        Assert.That(handler.Requests[0].Query, Is.EqualTo("?resultType=compact&pageSize=10&pageNumber=1&search=a%26pageSize%3D1%23b&species=7955"));
    }

    [Test]
    public void Massive_Outage_IsRetried_ThenTransient()
    {
        var handler = new Handler(_ => Body("busy", HttpStatusCode.ServiceUnavailable));
        using MassiveClient client = Massive(handler);

        var e = Assert.ThrowsAsync<HttpRequestException>(async () => await client.SearchDatasetsAsync("zebrafish"));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(e), Is.True);
            Assert.That(e!.Message, Does.Contain("failed with status 503"));
            Assert.That(handler.Requests, Has.Count.EqualTo(4));
        });
    }

    [TestCase("""{"datasets":[]}""")]
    [TestCase("not json")]
    [TestCase("""[{"title":"no accession"}]""")]
    [TestCase("""[42]""")]
    public void Massive_PageThatIsNotAListOfDatasets_IsAContractBreak(string body)
    {
        using MassiveClient client = Massive(new Handler(_ => Body(body)));

        Assert.ThrowsAsync<MzLibException>(async () => await client.SearchDatasetsAsync("zebrafish"));
    }

    [Test]
    public void Massive_ArgumentsAreCheckedBeforeAnyRequest()
    {
        var handler = new Handler(_ => Body(NotFound, HttpStatusCode.NotFound));
        using MassiveClient client = Massive(handler);

        Assert.Multiple(() =>
        {
            Assert.ThrowsAsync<ArgumentException>(async () => await client.SearchDatasetsAsync(" "));
            Assert.ThrowsAsync<ArgumentException>(async () => await client.SearchDatasetsAsync(new string('a', MassiveClient.MaxKeywordLength + 1)));
            Assert.ThrowsAsync<ArgumentException>(async () => await client.SearchDatasetsAsync("x", speciesTaxId: "Danio rerio"));
            Assert.ThrowsAsync<ArgumentException>(async () => await client.SearchDatasetsAsync("x", speciesTaxId: ""));
            Assert.ThrowsAsync<ArgumentOutOfRangeException>(async () => await client.SearchDatasetsAsync("x", pageSize: 0));
            Assert.That(handler.Requests, Is.Empty);
        });
        Assert.DoesNotThrowAsync(async () => await client.SearchDatasetsAsync(new string('a', MassiveClient.MaxKeywordLength)));
    }

    [Test]
    public void Massive_Cancelled_IsACancellation()
    {
        using MassiveClient client = Massive(new Handler(r => MassivePage(PageOf(r))));
        using var cancel = new CancellationTokenSource();
        cancel.Cancel();

        Assert.CatchAsync<OperationCanceledException>(async () => await client.SearchDatasetsAsync("zebrafish", cancellationToken: cancel.Token));
    }

    // ---- ProteomeCentral ---------------------------------------------------------------------------------------

    [Test]
    public async Task Central_EveryPageIsRead_ColumnsByName_EntitiesDecoded_NoLabHead()
    {
        var handler = new Handler(r => ProteomeCentralPage(PageOf(r)));
        using ProteomeCentralClient client = Central(handler);

        IReadOnlyList<ProteomeCentralDataset> hits = await client.SearchDatasetsAsync("zebrafish",
            new ProteomeCentralSearchFilter { Repository = "MassIVE" }, pageSize: 5);

        ProteomeCentralDataset first = hits[0];
        Assert.Multiple(() =>
        {
            Assert.That(hits, Has.Count.EqualTo(11));
            Assert.That(hits.All(h => h.Repository == "MassIVE"), Is.True);
            Assert.That(first.Accession, Is.EqualTo("PXD062174"));
            Assert.That(first.Title, Does.StartWith("CF-MS data of zebrafish brain lysate for 'VerteBrain"), "&apos; is decoded");
            Assert.That(first.Species, Is.EqualTo("Danio rerio"));
            Assert.That(first.Sdrf, Is.Null);
            Assert.That(first.AnnounceDate, Is.EqualTo("2026-04-28"));
            Assert.That(typeof(ProteomeCentralDataset).GetProperties().Any(p => p.Name.Contains("Lab")), Is.False);
            Assert.That(hits.SelectMany(h => typeof(ProteomeCentralDataset).GetProperties().Select(p => p.GetValue(h) as string))
                .Any(v => v != null && v.Contains("Redacted Lab Head")), Is.False);
            Assert.That(handler.Requests, Has.Count.EqualTo(3), "page 3 was short");
            Assert.That(handler.Requests[0].AbsoluteUri, Is.EqualTo(
                ProteomeCentralClient.DefaultBaseAddress + "datasets?resultType=compact&pageSize=5&pageNumber=1&search=zebrafish&repository=MassIVE"));
        });
    }

    [Test]
    public async Task Central_EmptyPageAfterAFullOne_EndsTheSearch()
    {
        var handler = new Handler(r => ProteomeCentralPage(PageOf(r) <= 2 ? PageOf(r) : 4));
        using ProteomeCentralClient client = Central(handler);

        IReadOnlyList<ProteomeCentralDataset> hits = await client.SearchDatasetsAsync("zebrafish", pageSize: 5);

        Assert.That(hits, Has.Count.EqualTo(10));
        Assert.That(handler.Requests, Has.Count.EqualTo(3));
    }

    [Test]
    public async Task Central_ColumnsThatMove_AreStillReadByName()
    {
        const string answer = """
            {"result_set":{"datasets_title_list":["title","repository","dataset identifier"]},
             "datasets":[["A title","PRIDE","PXD000001"]]}
            """;
        using ProteomeCentralClient client = Central(new Handler(_ => Body(answer)));

        ProteomeCentralDataset hit = (await client.SearchDatasetsAsync("x", pageSize: 5)).Single();

        Assert.Multiple(() =>
        {
            Assert.That(hit.Accession, Is.EqualTo("PXD000001"));
            Assert.That(hit.Title, Is.EqualTo("A title"));
            Assert.That(hit.Repository, Is.EqualTo("PRIDE"));
            Assert.That(hit.Species, Is.Null, "a column that is not there is null, never another column's value");
        });
    }

    [Test]
    public async Task Central_FilterOnly_SendsNoSearch_AndEveryFilterItKnows()
    {
        var handler = new Handler(_ => ProteomeCentralPage(4));
        using ProteomeCentralClient client = Central(handler);

        await client.SearchDatasetsAsync(null, new ProteomeCentralSearchFilter
        {
            Repository = "MassIVE", Species = "Danio rerio", Instrument = "Q Exactive", Year = "2026",
            Keywords = "DatasetType:Proteomics", Sdrf = "curated valid SDRF", Files = "1-5",
        });

        Assert.That(handler.Requests[0].Query, Is.EqualTo(
            "?resultType=compact&pageSize=100&pageNumber=1&repository=MassIVE&species=Danio%20rerio&instrument=Q%20Exactive" +
            "&year=2026&keywords=DatasetType%3AProteomics&sdrf=curated%20valid%20SDRF&files=1-5"));
    }

    [Test]
    public void Central_ServerIgnoringThePaging_IsCaught()
    {
        var handler = new Handler(_ => ProteomeCentralPage(1));
        using ProteomeCentralClient client = Central(handler);

        Assert.ThrowsAsync<MzLibException>(async () => await client.SearchDatasetsAsync("zebrafish", pageSize: 5));
        Assert.That(handler.Requests, Has.Count.EqualTo(2));
    }

    [Test]
    public void Central_PastMaxPages_Throws()
    {
        var handler = new Handler(r => Body($$"""{"result_set":{"datasets_title_list":["dataset identifier"]},"datasets":[["PXD{{PageOf(r):D6}}"]]}"""));
        using ProteomeCentralClient client = Central(handler, maxPages: 3);

        Assert.ThrowsAsync<HttpRequestException>(async () => await client.SearchDatasetsAsync("x", pageSize: 1));
        Assert.That(handler.Requests, Has.Count.EqualTo(3));
    }

    [Test]
    public void Central_Outage_IsRetried_ThenTransient_ARefusalIsNot()
    {
        var busy = new Handler(_ => Body("busy", HttpStatusCode.ServiceUnavailable));
        var refused = new Handler(_ => Body("bad", HttpStatusCode.BadRequest));
        using ProteomeCentralClient a = Central(busy);
        using ProteomeCentralClient b = Central(refused);

        var outage = Assert.ThrowsAsync<HttpRequestException>(async () => await a.SearchDatasetsAsync("x"));
        var refusal = Assert.ThrowsAsync<HttpRequestException>(async () => await b.SearchDatasetsAsync("x"));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(outage), Is.True);
            Assert.That(outage!.Message, Does.Contain("failed with status 503"));
            Assert.That(busy.Requests, Has.Count.EqualTo(4));
            Assert.That(ServiceFailure.IsTransient(refusal), Is.False);
            Assert.That(refused.Requests, Has.Count.EqualTo(1));
        });
    }

    [TestCase("not json")]
    [TestCase("""{"datasets":[]}""")]
    [TestCase("""{"result_set":{"datasets_title_list":["title"]},"datasets":[]}""")]
    [TestCase("""{"result_set":{"datasets_title_list":["dataset identifier"]},"datasets":[["   "]]}""")]
    [TestCase("""{"result_set":{"datasets_title_list":["dataset identifier"]},"datasets":["PXD000001"]}""")]
    public void Central_AnswerThatBreaksTheShape_IsAContractBreak(string body)
    {
        using ProteomeCentralClient client = Central(new Handler(_ => Body(body)));

        Assert.ThrowsAsync<MzLibException>(async () => await client.SearchDatasetsAsync("x"));
    }

    [Test]
    public void Central_ArgumentsAreCheckedBeforeAnyRequest()
    {
        var handler = new Handler(_ => ProteomeCentralPage(4));
        using ProteomeCentralClient client = Central(handler);

        Assert.Multiple(() =>
        {
            Assert.ThrowsAsync<ArgumentException>(async () => await client.SearchDatasetsAsync(null));
            Assert.ThrowsAsync<ArgumentException>(async () => await client.SearchDatasetsAsync(" ", new ProteomeCentralSearchFilter()));
            Assert.ThrowsAsync<ArgumentException>(async () => await client.SearchDatasetsAsync("x", new ProteomeCentralSearchFilter { Species = " " }));
            Assert.ThrowsAsync<ArgumentException>(async () => await client.SearchDatasetsAsync(new string('a', ProteomeCentralClient.MaxKeywordLength + 1)));
            Assert.ThrowsAsync<ArgumentOutOfRangeException>(async () => await client.SearchDatasetsAsync("x", pageSize: -1));
            Assert.That(handler.Requests, Is.Empty);
        });
    }

    [Test]
    public void Central_Constructors_AndDispose()
    {
        using (var owning = new ProteomeCentralClient())
            Assert.That(owning.MaxPages, Is.EqualTo(10000));
        Assert.Throws<ArgumentNullException>(() => _ = new ProteomeCentralClient(null!));

        ProteomeCentralClient client = Central(new Handler(_ => ProteomeCentralPage(4)));
        client.Dispose();
        client.Dispose();
        Assert.ThrowsAsync<ObjectDisposedException>(async () => await client.SearchDatasetsAsync("x"));
    }

    [Test]
    public void SearchResults_ToString_NameTheDataset()
    {
        Assert.That(new MassiveDatasetSearchResult { MassiveAccession = "MSV000000001", Title = "t" }.ToString(), Is.EqualTo("MSV000000001: t"));
        Assert.That(new ProteomeCentralDataset { Accession = "PXD000001", Repository = "PRIDE", Title = "t" }.ToString(), Is.EqualTo("PXD000001 (PRIDE): t"));
    }
}

/// <summary>
/// Live canaries for both searches. [Category("ExternalService")] keeps them out of the required run; a service
/// outage SKIPS, a format drift FAILS. One attempt each, a narrow keyword, and 6 s between tests (MassIVE rate-limits).
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("Massive")]
[ExcludeFromCodeCoverage]
public class MassiveSearchLiveTests
{
    [SetUp]
    public Task SetUp() => Task.Delay(TimeSpan.FromSeconds(6));

    [Test]
    public Task MassiveSearch_LiveZebrafish_FindsDatasetsProteomeCentralCannotSee() =>
        ExternalServiceTestHelper.RunAsync("MassIVE", async () =>
        {
            using var client = new MassiveClient { MaxRetries = 0 };

            IReadOnlyList<MassiveDatasetSearchResult> hits = await client.SearchDatasetsAsync("zebrafish", pageSize: 50);

            Assert.That(hits.Count, Is.GreaterThan(50), "more than one page");
            Assert.That(hits.Select(h => h.MassiveAccession).Distinct().Count(), Is.EqualTo(hits.Count));
            Assert.That(hits.Any(h => h.PxdAccession == null), Is.True, "MassIVE-only datasets (48 of 67 on 2026-11-10)");
        });

    [Test]
    public Task MassiveSearch_LiveNoMatch_IsEmpty() =>
        ExternalServiceTestHelper.RunAsync("MassIVE", async () =>
        {
            using var client = new MassiveClient { MaxRetries = 0 };

            Assert.That(await client.SearchDatasetsAsync("qqzzxxvvnomatch"), Is.Empty);
        });

    [Test]
    [Category("ProteomeCentral")]
    public Task CentralSearch_LiveRepositoryFilter_ReturnsOnlyThatRepository() =>
        ExternalServiceTestHelper.RunAsync("ProteomeCentral", async () =>
        {
            using var client = new ProteomeCentralClient { MaxRetries = 0 };

            IReadOnlyList<ProteomeCentralDataset> hits = await client.SearchDatasetsAsync("zebrafish",
                new ProteomeCentralSearchFilter { Repository = "MassIVE" }, pageSize: 5);

            Assert.That(hits.Count, Is.GreaterThan(5), "more than one page");
            Assert.That(hits.All(h => h.Repository == "MassIVE" && h.Accession.StartsWith("PXD")), Is.True);
        });
}
