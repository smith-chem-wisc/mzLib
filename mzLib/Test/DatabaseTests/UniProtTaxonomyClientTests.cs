using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using UsefulProteomicsDatabases.Taxonomy;

namespace Test.DatabaseTests;

/// <summary>
/// <see cref="UniProtTaxonomyClient"/> against recorded UniProt responses (DataFiles/Taxonomy, see its PROVENANCE.md),
/// served by a stub handler so no test here touches the network. The live canary is in
/// <see cref="UniProtTaxonomyClientLiveTests"/>.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class UniProtTaxonomyClientTests
{
    private const string Page2Url = "https://rest.uniprot.org/proteomes/search?query=reference%3Atrue&cursor=abc123&size=2";

    private sealed class StubHandler : HttpMessageHandler
    {
        private readonly Func<HttpRequestMessage, HttpResponseMessage> _responder;
        public List<string> RequestedUris { get; } = new();

        public StubHandler(Func<HttpRequestMessage, HttpResponseMessage> responder) => _responder = responder;

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            RequestedUris.Add(request.RequestUri.AbsoluteUri);
            return Task.FromResult(_responder(request));
        }
    }

    private static string Fixture(string name) =>
        File.ReadAllText(Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Taxonomy", name));

    private static HttpResponseMessage Json(string body, string nextLink = null, string release = "2026_04")
    {
        var response = new HttpResponseMessage(HttpStatusCode.OK) { Content = new StringContent(body, Encoding.UTF8, "application/json") };
        if (release != null)
            response.Headers.TryAddWithoutValidation("x-uniprot-release", release);
        if (nextLink != null)
            response.Headers.TryAddWithoutValidation("Link", $"<{nextLink}>; rel=\"next\"");
        return response;
    }

    /// <summary>The three recorded responses: two pages of proteomes joined by a Link header, then the organisms' ranks.</summary>
    private static HttpResponseMessage Recorded(HttpRequestMessage request)
    {
        string uri = request.RequestUri.AbsoluteUri;
        if (uri.Contains("/taxonomy/taxonIds/"))
            return Json(Fixture("uniprot-taxa.json"));
        if (uri.Contains("cursor="))
            return Json(Fixture("uniprot-proteomes-page2.json"));
        return Json(Fixture("uniprot-proteomes-page1.json"), nextLink: Page2Url);
    }

    private static UniProtTaxonomyClient Client(StubHandler handler)
    {
        var client = new UniProtTaxonomyClient(new HttpClient(handler));
        client.Engine.RetryDelay = (_, _) => Task.CompletedTask;
        return client;
    }

    [Test]
    public async Task BuildsTheTreeAboveEveryProteomeAcrossPages()
    {
        var handler = new StubHandler(Recorded);
        using UniProtTaxonomyClient client = Client(handler);

        UniProtTaxonomySnapshot snapshot = await client.FetchAsync(pageSize: 2);
        TaxonomyTree tree = snapshot.Tree;

        Assert.That(snapshot.Proteomes.Proteomes.Select(p => (p.ProteomeId, p.TaxonId)), Is.EqualTo(new[]
        {
            ("UP000000625", 83333L), ("UP000005640", 9606L), ("UP000464024", 2697049L),
        }));
        Assert.That(snapshot.Proteomes.Proteomes.All(p => tree.Contains(p.TaxonId)), Is.True,
            "every proteome's taxon is in the tree built from the same pages");
        Assert.That(snapshot.ProteomesWithoutLineage, Is.Empty);

        Assert.That(tree.AncestorAtRank(83333, TaxonomicRank.Family).Taxon!.ScientificName, Is.EqualTo("Enterobacteriaceae"));
        Assert.That(tree.AncestorAtRank(83333, TaxonomicRank.Species).Taxon!.ScientificName, Is.EqualTo("Escherichia coli"));
        Assert.That(tree.AncestorAtRank(83333, TaxonomicRank.Domain).Taxon!.ScientificName, Is.EqualTo("Bacteria"));
        Assert.That(tree.AncestorAtRank(9606, TaxonomicRank.Family).Taxon!.ScientificName, Is.EqualTo("Hominidae"));
        Assert.That(tree.AncestorAtRank(9606, TaxonomicRank.Class).Taxon!.ScientificName, Is.EqualTo("Mammalia"));
        Assert.That(tree.AncestorAtRank(2697049, TaxonomicRank.Order).Taxon!.ScientificName, Is.EqualTo("Nidovirales"));
        Assert.That(tree.AncestorAtRank(2697049, TaxonomicRank.Family).Taxon!.ScientificName, Is.EqualTo("Coronaviridae"));
        Assert.That(tree.AncestorAtRank(2697049, TaxonomicRank.Species).Taxon!.ScientificName, Is.EqualTo("Betacoronavirus pandemicum"));

        // The human and bacterial lineages share their first ancestor and nothing below it.
        Assert.That(tree.Ancestors(9606).Last().TaxonId, Is.EqualTo(tree.Ancestors(83333).Last().TaxonId));
        Assert.That(tree.Ancestors(9606).Select(n => n.TaxonId).Intersect(tree.Ancestors(83333).Select(n => n.TaxonId)).Count(), Is.EqualTo(1));
    }

    /// <summary>
    /// The proteome listing gives Homo sapiens' ancestors their ranks but not Homo sapiens its own, so without the
    /// second request "the species of the human proteome" would have no answer.
    /// </summary>
    [Test]
    public async Task ReadsEachOrganismsOwnRankFromTheTaxonomyEndpoint()
    {
        var handler = new StubHandler(Recorded);
        using UniProtTaxonomyClient client = Client(handler);

        TaxonomyTree tree = (await client.FetchAsync(pageSize: 2)).Tree;

        Assert.That(tree.TryGetNode(9606, out TaxonomyNode human), Is.True);
        Assert.That(human.Rank, Is.EqualTo("species"));
        Assert.That(tree.AncestorAtRank(9606, TaxonomicRank.Species).Taxon!.TaxonId, Is.EqualTo(9606));
        Assert.That(tree.AncestorAtRank(9606, TaxonomicRank.Genus).Taxon!.ScientificName, Is.EqualTo("Homo"));
        Assert.That(tree.TryGetNode(83333, out TaxonomyNode k12), Is.True);
        Assert.That(k12.Rank, Is.EqualTo("strain"));
        Assert.That(tree.TryGetNode(2697049, out TaxonomyNode sars2), Is.True);
        Assert.That(sars2.Rank, Is.EqualTo("no rank"));

        string taxaRequest = handler.RequestedUris.Single(u => u.Contains("/taxonomy/taxonIds/"));
        Assert.That(taxaRequest, Does.Contain("/taxonomy/taxonIds/9606,83333,2697049?"), "one request for all three, in identifier order");
    }

    [Test]
    public async Task AsksForTheQueryAndFollowsTheNextLinkExactly()
    {
        var handler = new StubHandler(Recorded);
        using UniProtTaxonomyClient client = Client(handler);

        UniProtTaxonomySnapshot snapshot = await client.FetchAsync("reference:true", pageSize: 2);

        Assert.That(handler.RequestedUris.Take(2), Is.EqualTo(new[]
        {
            "https://rest.uniprot.org/proteomes/search?query=reference%3Atrue&size=2",
            Page2Url,
        }));
        Assert.That(handler.RequestedUris, Has.Count.EqualTo(3));
        var expected = new TaxonomyProvenance("https://rest.uniprot.org/proteomes/search", "2026_04", "reference:true");
        Assert.That(snapshot.Tree.Provenance, Is.EqualTo(expected));
        Assert.That(snapshot.Proteomes.Provenance, Is.EqualTo(expected));
    }

    [Test]
    public async Task KeepsWhatUniProtSaysAboutEachProteome()
    {
        using UniProtTaxonomyClient client = Client(new StubHandler(Recorded));

        ProteomeTaxonMap proteomes = (await client.FetchAsync(pageSize: 2)).Proteomes;

        Assert.That(proteomes.TryGetProteome("UP000000625", out ProteomeTaxon k12), Is.True);
        Assert.That(k12, Is.EqualTo(new ProteomeTaxon("UP000000625", 83333, "Escherichia coli (strain K12)", 4403, "Reference proteome")));
    }

    [Test]
    public async Task TheSnapshotSurvivesBeingWrittenAndReadBack()
    {
        using UniProtTaxonomyClient client = Client(new StubHandler(Recorded));
        UniProtTaxonomySnapshot snapshot = await client.FetchAsync(pageSize: 2);
        string directory = Path.Combine(TestContext.CurrentContext.WorkDirectory, "UniProtTaxonomyClientTests_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(directory);
        try
        {
            string nodes = Path.Combine(directory, "taxonomy_nodes.tsv.gz");
            string proteomes = Path.Combine(directory, "proteome_taxa.tsv.gz");
            TaxonomySnapshotTsv.WriteNodes(nodes, snapshot.Tree);
            TaxonomySnapshotTsv.WriteProteomes(proteomes, snapshot.Proteomes);

            Assert.That(TaxonomyTree.Load(nodes).Nodes, Is.EqualTo(snapshot.Tree.Nodes));
            Assert.That(ProteomeTaxonMap.Load(proteomes).Proteomes, Is.EqualTo(snapshot.Proteomes.Proteomes));
        }
        finally
        {
            Directory.Delete(directory, recursive: true);
        }
    }

    // ---- answers the client cannot use: MzLibException ------------------------

    [Test]
    public void NoMatchingProteome_ThrowsMzLibException()
    {
        using UniProtTaxonomyClient client = Client(new StubHandler(_ => Json("{\"results\":[]}")));

        var ex = Assert.ThrowsAsync<MzLibException>(() => client.FetchAsync("upid:UP999999999"));
        Assert.That(ex.Message, Does.Contain("upid:UP999999999"));
    }

    [Test]
    [TestCase("<html>Service unavailable</html>", "not JSON")]
    [TestCase("{\"error\":\"nope\"}", "no 'results' list")]
    [TestCase("{\"results\":[{\"taxonomy\":{\"taxonId\":9606}}]}", "no identifier")]
    [TestCase("{\"results\":[{\"id\":\"UP1\"}]}", "UP1 with no taxon")]
    [TestCase("{\"results\":[{\"id\":\"UP1\",\"taxonomy\":{\"taxonId\":9606},\"taxonLineage\":[{\"scientificName\":\"x\"}]}]}", "lineage entry")]
    public void AnAnswerThatIsNotAProteomeListing_ThrowsMzLibException(string body, string expectedInMessage)
    {
        using UniProtTaxonomyClient client = Client(new StubHandler(_ => Json(body)));

        var ex = Assert.ThrowsAsync<MzLibException>(() => client.FetchAsync());
        Assert.That(ex.Message, Does.Contain(expectedInMessage));
    }

    /// <summary>Two proteomes that disagree about a taxon's parent cannot both be placed in one tree.</summary>
    [Test]
    public void ATaxonWithTwoDifferentParents_ThrowsMzLibException()
    {
        const string body = """
            {"results":[
             {"id":"UP1","taxonomy":{"taxonId":10,"scientificName":"A"},
              "taxonLineage":[{"taxonId":1,"scientificName":"Root","rank":"domain","hidden":false},{"taxonId":5,"scientificName":"G","rank":"genus","hidden":false}]},
             {"id":"UP2","taxonomy":{"taxonId":11,"scientificName":"B"},
              "taxonLineage":[{"taxonId":2,"scientificName":"Other","rank":"domain","hidden":false},{"taxonId":5,"scientificName":"G","rank":"genus","hidden":false}]}
            ]}
            """;
        using UniProtTaxonomyClient client = Client(new StubHandler(_ => Json(body)));

        var ex = Assert.ThrowsAsync<MzLibException>(() => client.FetchAsync());
        Assert.That(ex.Message, Does.Contain("taxon 5").And.Contain("UP2"));
    }

    [Test]
    public void AProteomeListedTwice_ThrowsMzLibException()
    {
        const string body = """
            {"results":[{"id":"UP1","taxonomy":{"taxonId":10,"scientificName":"A"}},{"id":"UP1","taxonomy":{"taxonId":10,"scientificName":"A"}}]}
            """;
        using UniProtTaxonomyClient client = Client(new StubHandler(_ => Json(body)));

        var ex = Assert.ThrowsAsync<MzLibException>(() => client.FetchAsync());
        Assert.That(ex.Message, Does.Contain("UP1 twice"));
    }

    /// <summary>A server that keeps handing back the same next link must not hold the fetch open for ever.</summary>
    [Test]
    public void ANextLinkThatRepeats_ThrowsMzLibException()
    {
        var handler = new StubHandler(request => request.RequestUri.AbsoluteUri.Contains("cursor=")
            ? Json(Fixture("uniprot-proteomes-page2.json"), nextLink: Page2Url)
            : Json(Fixture("uniprot-proteomes-page1.json"), nextLink: Page2Url));
        using UniProtTaxonomyClient client = Client(handler);

        var ex = Assert.ThrowsAsync<MzLibException>(() => client.FetchAsync(pageSize: 2));
        Assert.That(ex.Message, Does.Contain("did not end"));
        Assert.That(handler.RequestedUris, Has.Count.EqualTo(2));
    }

    [Test]
    public async Task AnEmptyPageEndsTheListingEvenWithANextLink()
    {
        var handler = new StubHandler(request =>
        {
            string uri = request.RequestUri.AbsoluteUri;
            if (uri.Contains("/taxonomy/taxonIds/"))
                return Json(Fixture("uniprot-taxa.json"));
            return uri.Contains("cursor=")
                ? Json("{\"results\":[]}", nextLink: Page2Url + "&again")
                : Json(Fixture("uniprot-proteomes-page1.json"), nextLink: Page2Url);
        });
        using UniProtTaxonomyClient client = Client(handler);

        UniProtTaxonomySnapshot snapshot = await client.FetchAsync(pageSize: 2);

        Assert.That(snapshot.Proteomes.Count, Is.EqualTo(2));
    }

    /// <summary>
    /// UniProt lists a few reference proteomes under a taxon NCBI has since merged or withdrawn: a number, with no
    /// name and no lineage (22 of 37,623 in release 2026_04; UP000031433 is one). Such a taxon gets no node, and the
    /// proteome is named so the caller knows it could not be placed.
    /// </summary>
    [Test]
    public async Task AProteomeWhoseTaxonHasNoNameOrLineage_IsRecordedButNotPlaced()
    {
        const string withdrawn = "{\"id\":\"UP000031433\",\"taxonomy\":{\"taxonId\":1510391},\"proteomeType\":\"Reference proteome\",\"proteinCount\":3151,\"taxonLineage\":null}";
        string page = Fixture("uniprot-proteomes-page2.json").Replace("\"results\": [", "\"results\": [" + withdrawn + ",");
        var handler = new StubHandler(request => request.RequestUri.AbsoluteUri.Contains("/taxonomy/taxonIds/")
            ? Json(Fixture("uniprot-taxa.json"))
            : Json(page));
        using UniProtTaxonomyClient client = Client(handler);

        UniProtTaxonomySnapshot snapshot = await client.FetchAsync();

        Assert.That(snapshot.ProteomesWithoutLineage, Is.EqualTo(new[] { "UP000031433" }));
        Assert.That(snapshot.Proteomes.TryGetProteome("UP000031433", out ProteomeTaxon proteome), Is.True);
        Assert.That(proteome, Is.EqualTo(new ProteomeTaxon("UP000031433", 1510391, "", 3151, "Reference proteome")));
        Assert.That(snapshot.Tree.Contains(1510391), Is.False);
        Assert.That(snapshot.Tree.AncestorAtRank(1510391, TaxonomicRank.Family).Status, Is.EqualTo(TaxonRankStatus.TaxonNotInTree));
        Assert.That(snapshot.Tree.Nodes.Select(n => n.ScientificName), Has.None.Empty);
        Assert.That(snapshot.Tree.Contains(2697049), Is.True, "the other proteome on the page is placed as usual");
    }

    /// <summary>A taxon the taxonomy endpoint does not return keeps the empty rank the listing gave it.</summary>
    [Test]
    public async Task AnOrganismTheTaxonomyEndpointOmits_KeepsNoRank()
    {
        var handler = new StubHandler(request => request.RequestUri.AbsoluteUri.Contains("/taxonomy/taxonIds/")
            ? Json("{\"results\":[]}")
            : Json(Fixture("uniprot-proteomes-page2.json")));
        using UniProtTaxonomyClient client = Client(handler);

        TaxonomyTree tree = (await client.FetchAsync()).Tree;

        Assert.That(tree.TryGetNode(2697049, out TaxonomyNode sars2), Is.True);
        Assert.That(sars2.Rank, Is.Empty);
        Assert.That(tree.AncestorAtRank(2697049, TaxonomicRank.Family).Taxon!.ScientificName, Is.EqualTo("Coronaviridae"));
    }

    [Test]
    public void TheTaxonomyEndpointDisagreeingAboutAParent_ThrowsMzLibException()
    {
        var handler = new StubHandler(request => request.RequestUri.AbsoluteUri.Contains("/taxonomy/taxonIds/")
            ? Json("{\"results\":[{\"taxonId\":2697049,\"rank\":\"no rank\",\"parent\":{\"taxonId\":42}}]}")
            : Json(Fixture("uniprot-proteomes-page2.json")));
        using UniProtTaxonomyClient client = Client(handler);

        var ex = Assert.ThrowsAsync<MzLibException>(() => client.FetchAsync());
        Assert.That(ex.Message, Does.Contain("2697049").And.Contain("42"));
    }

    // ---- outages: HttpRequestException ----------------------------------------

    [Test]
    [TestCase(HttpStatusCode.InternalServerError)]
    [TestCase(HttpStatusCode.ServiceUnavailable)]
    [TestCase(HttpStatusCode.BadRequest)]
    public void AFailedRequest_ThrowsHttpRequestException(HttpStatusCode status)
    {
        using UniProtTaxonomyClient client = Client(new StubHandler(_ => new HttpResponseMessage(status)));
        client.Engine.MaxRetries = 0;

        var ex = Assert.ThrowsAsync<HttpRequestException>(() => client.FetchAsync());
        Assert.That(ex.Message, Does.Contain(((int)status).ToString()));
    }

    [Test]
    public async Task ATransientFailureOnOnePageIsRetriedAndTheFetchCompletes()
    {
        int page2Attempts = 0;
        var handler = new StubHandler(request =>
        {
            if (request.RequestUri.AbsoluteUri.Contains("cursor=") && page2Attempts++ == 0)
                return new HttpResponseMessage(HttpStatusCode.ServiceUnavailable);
            return Recorded(request);
        });
        using UniProtTaxonomyClient client = Client(handler);

        UniProtTaxonomySnapshot snapshot = await client.FetchAsync(pageSize: 2);

        Assert.That(snapshot.Proteomes.Count, Is.EqualTo(3));
        Assert.That(page2Attempts, Is.EqualTo(2));
        Assert.That(handler.RequestedUris.Count(u => !u.Contains("cursor=") && u.Contains("/proteomes/")), Is.EqualTo(1),
            "only the failed page is asked for again");
    }

    // ---- arguments and cancellation -------------------------------------------

    [Test]
    [TestCase(null)]
    [TestCase(" ")]
    public void ABlankQuery_ThrowsArgumentException(string query)
    {
        using UniProtTaxonomyClient client = Client(new StubHandler(Recorded));

        Assert.ThrowsAsync<ArgumentException>(() => client.FetchAsync(query));
    }

    [Test]
    [TestCase(0)]
    [TestCase(501)]
    public void APageSizeOutOfRange_ThrowsArgumentOutOfRange(int pageSize)
    {
        using UniProtTaxonomyClient client = Client(new StubHandler(Recorded));

        Assert.ThrowsAsync<ArgumentOutOfRangeException>(() => client.FetchAsync(pageSize: pageSize));
    }

    [Test]
    public void Cancellation_StopsTheFetch()
    {
        var handler = new StubHandler(Recorded);
        using UniProtTaxonomyClient client = Client(handler);
        using var cancelled = new CancellationTokenSource();
        cancelled.Cancel();

        Assert.That(async () => await client.FetchAsync(cancellationToken: cancelled.Token),
            Throws.InstanceOf<OperationCanceledException>());
        Assert.That(handler.RequestedUris, Is.Empty);
    }

    [Test]
    public void ANullHttpClient_ThrowsArgumentNull()
    {
        Assert.Throws<ArgumentNullException>(() => new UniProtTaxonomyClient(null));
    }

    [Test]
    public void TheParameterlessClient_OwnsAndDisposesItsConnection()
    {
        var client = new UniProtTaxonomyClient();

        Assert.DoesNotThrow(client.Dispose);
    }
}

/// <summary>One request to the real service, so a change in UniProt's response shape is noticed. Skipped when UniProt is unreachable.</summary>
[TestFixture]
[Category("ExternalService")]
[Category("UniProt")]
[ExcludeFromCodeCoverage]
public class UniProtTaxonomyClientLiveTests
{
    [Test]
    public Task Live_TheHumanReferenceProteome_IsPlacedInItsFamily() =>
        ExternalServiceTestHelper.RunAsync("UniProt", async () =>
        {
            using var client = new UniProtTaxonomyClient();

            UniProtTaxonomySnapshot snapshot = await client.FetchAsync("upid:UP000005640");

            Assert.That(snapshot.Proteomes.TryGetTaxonId("UP000005640", out long human), Is.True);
            Assert.That(human, Is.EqualTo(9606));
            Assert.That(snapshot.Tree.AncestorAtRank(9606, TaxonomicRank.Species).Taxon!.TaxonId, Is.EqualTo(9606));
            Assert.That(snapshot.Tree.AncestorAtRank(9606, TaxonomicRank.Family).Taxon!.ScientificName, Is.EqualTo("Hominidae"));
            Assert.That(snapshot.Tree.Provenance.Release, Is.Not.Null.And.Not.Empty);
        });
}
