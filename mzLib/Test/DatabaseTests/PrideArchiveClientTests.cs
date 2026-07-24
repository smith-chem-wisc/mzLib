using NUnit.Framework;
using System;
using System.Diagnostics.CodeAnalysis;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using MzLibUtil;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

[TestFixture]
[ExcludeFromCodeCoverage]
public class PrideArchiveClientTests
{
    // ---- test doubles -------------------------------------------------------

    /// <summary>An HttpMessageHandler that returns a caller-supplied response and records request URIs.</summary>
    private sealed class StubHandler : HttpMessageHandler
    {
        private readonly Func<HttpRequestMessage, HttpResponseMessage> _responder;
        public List<string> RequestedUris { get; } = new();

        public StubHandler(Func<HttpRequestMessage, HttpResponseMessage> responder) => _responder = responder;

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            RequestedUris.Add(request.RequestUri.ToString());
            return Task.FromResult(_responder(request));
        }
    }

    private static HttpResponseMessage JsonResponse(string body, HttpStatusCode status = HttpStatusCode.OK, long? totalRecords = null)
    {
        var response = new HttpResponseMessage(status) { Content = new StringContent(body, Encoding.UTF8, "application/json") };
        if (totalRecords.HasValue)
            response.Headers.TryAddWithoutValidation("total_records", totalRecords.Value.ToString());
        return response;
    }

    /// <summary>A single PRIDE file object in the real v3 wire shape.</summary>
    private static string FileJson(string fileName, string category = "SEARCH", string categoryAccession = "PRIDE:0000408") =>
        $$"""
        {
          "projectAccessions": ["PXD012345"],
          "accession": "hashid",
          "fileCategory": { "@type": "CvParam", "cvLabel": "PRIDE", "accession": "{{categoryAccession}}", "name": "category", "value": "{{category}}" },
          "checksum": "",
          "publicFileLocations": [
            { "@type": "CvParam", "cvLabel": "PRIDE", "accession": "PRIDE:0000469", "name": "FTP Protocol", "value": "ftp://ftp.pride.ebi.ac.uk/{{fileName}}" },
            { "@type": "CvParam", "cvLabel": "PRIDE", "accession": "PRIDE:0000468", "name": "Aspera Protocol", "value": "prd_ascp@fasp.ebi.ac.uk:{{fileName}}" }
          ],
          "fileSizeBytes": 96358400,
          "fileName": "{{fileName}}",
          "compress": false,
          "submissionDate": "2019-01-15T09:42:57.000+00:00",
          "publicationDate": "2025-07-13T23:01:04.308+00:00",
          "updatedDate": "2019-01-15T09:47:55.000+00:00",
          "additionalAttributes": [],
          "totalDownloads": 19
        }
        """;

    private static string Array(params string[] fileJson) => "[" + string.Join(",", fileJson) + "]";

    // ---- deserialization ----------------------------------------------------

    [Test]
    public async Task GetProjectFilesAsync_DeserializesFields_AndResolvesFtpByAccession()
    {
        var handler = new StubHandler(_ => JsonResponse(Array(FileJson("run1.raw", "RAW", "PRIDE:0000404"))));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var files = await client.GetProjectFilesAsync("PXD012345");

        Assert.That(files.Count, Is.EqualTo(1));
        var file = files[0];
        Assert.That(file.FileName, Is.EqualTo("run1.raw"));
        Assert.That(file.FileSizeBytes, Is.EqualTo(96358400));
        Assert.That(file.FileCategory.Value, Is.EqualTo("RAW"));
        Assert.That(file.FileCategory.Accession, Is.EqualTo("PRIDE:0000404"));
        // resolve the FTP location by its stable accession, not by display name
        var ftp = file.PublicFileLocations.FirstOrDefault(l => l.Accession == "PRIDE:0000469");
        Assert.That(ftp, Is.Not.Null);
        Assert.That(ftp.Value, Does.StartWith("ftp://"));
        Assert.That(file.SubmissionDate.Year, Is.EqualTo(2019));
    }

    // ---- paging -------------------------------------------------------------

    [Test]
    public async Task GetProjectFilesAsync_PagesUntilTotalRecords_ConcatenatesAllFiles()
    {
        const long total = 3;
        var handler = new StubHandler(request =>
        {
            var uri = request.RequestUri.ToString();
            if (uri.Contains("page=0")) return JsonResponse(Array(FileJson("a"), FileJson("b")), totalRecords: total);
            if (uri.Contains("page=1")) return JsonResponse(Array(FileJson("c")), totalRecords: total);
            return JsonResponse("[]", totalRecords: total);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var files = await client.GetProjectFilesAsync("PXD012345", pageSize: 2);

        Assert.That(files.Select(f => f.FileName), Is.EqualTo(new[] { "a", "b", "c" }));
        Assert.That(handler.RequestedUris.Count, Is.EqualTo(2)); // stopped at total_records; no wasted 3rd request
    }

    [Test]
    public async Task GetProjectFilesAsync_NoTotalHeader_StopsOnShortPage()
    {
        var handler = new StubHandler(request =>
        {
            var uri = request.RequestUri.ToString();
            if (uri.Contains("page=0")) return JsonResponse(Array(FileJson("a"), FileJson("b"))); // full page, no header
            return JsonResponse(Array(FileJson("c")));                                            // short page -> last
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var files = await client.GetProjectFilesAsync("PXD012345", pageSize: 2);

        Assert.That(files.Count, Is.EqualTo(3));
        Assert.That(handler.RequestedUris.Count, Is.EqualTo(2));
    }

    [Test]
    public async Task GetProjectFilesAsync_UnparseableTotalHeader_StillTerminates()
    {
        var handler = new StubHandler(_ =>
        {
            var response = JsonResponse(Array(FileJson("a"))); // short page (1 < pageSize)
            response.Headers.TryAddWithoutValidation("total_records", "not-a-number");
            return response;
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var files = await client.GetProjectFilesAsync("PXD012345", pageSize: 2);

        Assert.That(files.Count, Is.EqualTo(1));
    }

    [Test]
    public async Task GetProjectFilesAsync_ServerCapsPageSize_StillReturnsFullManifest()
    {
        // Regression: PRIDE caps pageSize server-side (100 as of 2026-07-23) and then pages by the
        // capped size, so asking for 4 yields a 2-file page that is NOT the last page. Treating a
        // short page as terminal truncated the manifest silently -- it returned 2 of 3 files with
        // no error, while the doc promised the full manifest regardless of page size.
        // The stub serves at most 2 files per page whatever is asked for, standing in for that cap.
        const long total = 3;
        var handler = new StubHandler(request =>
        {
            var uri = request.RequestUri.ToString();
            if (uri.Contains("page=0")) return JsonResponse(Array(FileJson("a"), FileJson("b")), totalRecords: total);
            if (uri.Contains("page=1")) return JsonResponse(Array(FileJson("c")), totalRecords: total);
            return JsonResponse("[]", totalRecords: total);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        // Ask for more per page than the stub will ever hand back, exactly as a caller passing 500 does.
        var files = await client.GetProjectFilesAsync("PXD012345", pageSize: 4);

        // Parse each recorded request into its paging parameters, so the assertions pin the paging
        // semantics -- two requests, both at the CALLER's page size, pages 0 then 1 -- without coupling
        // to query-string order, spelling, or the presence of unrelated parameters. URL construction is
        // owned elsewhere and is explicitly not what this PR changes; a reordered or extended query
        // (an added sort/filter) must not read as a paging regression here.
        List<ILookup<string, string>> requestedQueries = handler.RequestedUris
            .Select(u => new Uri(u).Query.TrimStart('?')
                .Split('&', StringSplitOptions.RemoveEmptyEntries)
                .Select(pair => pair.Split('=', 2))
                .ToLookup(kv => kv[0], kv => kv.Length > 1 ? kv[1] : string.Empty))
            .ToList();

        // Read through the lookup rather than an indexer: a renamed or dropped paging parameter -- the
        // very drift this rewrite loosened the test against -- then shows up as a null in the compared
        // sequence and produces the intended message, instead of throwing KeyNotFoundException out of
        // the Assert.Multiple block. A repeated key would likewise break ToDictionary but not ToLookup.
        static string Param(ILookup<string, string> query, string name) => query[name].FirstOrDefault();

        // A precondition, kept outside the multiple block so the assertions below cannot run against
        // a request list of unexpected length.
        Assert.That(requestedQueries, Has.Count.EqualTo(2));

        Assert.Multiple(() =>
        {
            Assert.That(files.Select(f => f.FileName), Is.EqualTo(new[] { "a", "b", "c" }));
            // A SECOND page was requested, still at the CALLER's page size, after a first page came
            // back shorter than asked for -- that is the invariant, not the exact query text.
            Assert.That(requestedQueries.Select(q => Param(q, "pageSize")), Is.EqualTo(new[] { "4", "4" }));
            Assert.That(requestedQueries.Select(q => Param(q, "page")), Is.EqualTo(new[] { "0", "1" }));
        });
    }

    /// <summary>
    /// The header is authoritative for whether to keep paging, so when it overstates what the server
    /// will actually serve the empty-page check is the only thing left to stop the loop. That path
    /// became load-bearing with the capped-pageSize fix -- before it, the short-page break caught this
    /// case first -- so it is pinned here rather than left to the argument in the PR body.
    /// </summary>
    [Test]
    public async Task GetProjectFilesAsync_TotalRecordsOverstated_StopsOnEmptyPage()
    {
        const long overstatedTotal = 5; // the stub will only ever serve 2
        var handler = new StubHandler(request =>
            request.RequestUri.ToString().Contains("page=0")
                ? JsonResponse(Array(FileJson("a"), FileJson("b")), totalRecords: overstatedTotal)
                : JsonResponse("[]", totalRecords: overstatedTotal));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        // pageSize is deliberately ABOVE what the stub serves, and that is what makes this test
        // discriminate. Asking for 2 would make the old short-page check a no-op (2 < 2 is false),
        // so the pre-fix code produced an identical result and the test would pin nothing.
        var files = await client.GetProjectFilesAsync("PXD012345", pageSize: 4);

        Assert.Multiple(() =>
        {
            // The shortfall against the reported total is accepted, not thrown.
            Assert.That(files.Select(f => f.FileName), Is.EqualTo(new[] { "a", "b" }));
            // 2 requests, not 1: the overstated header kept us paging past the short page, and the
            // EMPTY page -- not the short one -- is what stopped us. This also pins the cost of the
            // fix, so an "optimisation" that restores the early break fails here.
            Assert.That(handler.RequestedUris.Count, Is.EqualTo(2));
        });
    }

    [Test]
    public async Task GetProjectFilesAsync_TotalRecordsExactlyMatched_DoesNotFetchExtraPage()
    {
        // The counterpart to the overstated case. Note the stub ignores `page` and serves the same
        // two files forever, never an empty page: that makes the total_records check the ONLY thing
        // that can stop this loop on the first request, so if it is ever weakened this test does not
        // fail by one wasted request -- the second request re-serves an identical page and trips the
        // identical-page guard instead. PagesUntilTotalRecords cannot pin that, because its fallback
        // returns [] and the empty-page check would rescue it.
        const long total = 2;
        var handler = new StubHandler(_ => JsonResponse(Array(FileJson("a"), FileJson("b")), totalRecords: total));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var files = await client.GetProjectFilesAsync("PXD012345", pageSize: 2);

        Assert.Multiple(() =>
        {
            Assert.That(files.Select(f => f.FileName), Is.EqualTo(new[] { "a", "b" }));
            Assert.That(handler.RequestedUris.Count, Is.EqualTo(1));
        });
    }

    [Test]
    public void GetProjectFilesAsync_ServerIgnoresPageButReportsTotal_ThrowsInsteadOfReturningDuplicates()
    {
        // A server that honours the capped pageSize but ignores `page` re-serves the same records on
        // every request while still reporting total_records. Left alone, files.Count accumulates the
        // duplicates and trips `files.Count >= total`, returning a manifest padded with duplicates and
        // silently missing the records past the cap as a "success".
        const long total = 3; // more than the stub will ever serve on its single, repeated page
        var handler = new StubHandler(_ => JsonResponse(Array(FileJson("a"), FileJson("b")), totalRecords: total));
        using var client = new PrideArchiveClient(new HttpClient(handler)); // MaxPages left at its default

        // MzLibException, not HttpRequestException: the server answered successfully and broke the
        // contract, which the class remarks reserve HttpRequestException from precisely so that
        // ExternalServiceTestHelper does not convert a real regression into a green skip.
        Assert.That(async () => await client.GetProjectFilesAsync("PXD012345", pageSize: 4),
            Throws.InstanceOf<MzLibException>());

        // The count is what makes this test discriminating rather than decorative. The repeated page is
        // detected on the SECOND request and aborts there. Delete the identical-page guard and the loop
        // instead runs all the way to the MaxPages backstop, so this assertion goes red -- asserting the
        // exception type alone would not, because MaxPages throws too.
        Assert.That(handler.RequestedUris.Count, Is.EqualTo(2));
    }

    /// <summary>
    /// Paging progress is judged per PAGE, never per file name, so no manifest entry is ever dropped.
    /// A manifest may legitimately repeat a leaf name (the same name under different
    /// publicFileLocations), and deciding membership by name uniqueness would silently lose the
    /// duplicate -- reintroducing the truncation this method exists to prevent.
    /// </summary>
    [Test]
    public async Task GetProjectFilesAsync_ManifestRepeatsAFileName_ReturnsEveryEntry()
    {
        const long total = 3;
        var handler = new StubHandler(request =>
            request.RequestUri.ToString().Contains("page=0")
                ? JsonResponse(Array(FileJson("a"), FileJson("dup"), FileJson("dup")), totalRecords: total)
                : JsonResponse("[]", totalRecords: total));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var files = await client.GetProjectFilesAsync("PXD012345", pageSize: 4);

        Assert.That(files.Select(f => f.FileName), Is.EqualTo(new[] { "a", "dup", "dup" }));
        Assert.That(files.Count, Is.EqualTo((int)total));
    }

    /// <summary>
    /// An entry whose fileName is absent from the JSON is a case DownloadFileAsync explicitly guards
    /// against, so the manifest must carry every such entry rather than collapsing them together.
    /// </summary>
    [Test]
    public async Task GetProjectFilesAsync_EntriesWithNoFileName_AreAllReturned()
    {
        const long total = 2;
        var handler = new StubHandler(request =>
            request.RequestUri.ToString().Contains("page=0")
                ? JsonResponse("""[{ "fileSizeBytes": 1 }, { "fileSizeBytes": 2 }]""", totalRecords: total)
                : JsonResponse("[]", totalRecords: total));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var files = await client.GetProjectFilesAsync("PXD012345", pageSize: 4);

        Assert.That(files.Count, Is.EqualTo(2), "neither unnamed entry may be dropped from the manifest");
    }

    /// <summary>
    /// The documented boundary of the narrowed guarantee: with no total_records header there is nothing
    /// authoritative to page against, so a server-side cap below the requested pageSize still truncates
    /// on the first short page. This pins the limitation the summary now admits to, in both directions.
    /// </summary>
    [Test]
    public async Task GetProjectFilesAsync_ServerCapsPageSizeWithoutTotalHeader_StopsOnShortPage()
    {
        // No total_records anywhere, and the stub caps every page at 2 however many are asked for.
        var handler = new StubHandler(request =>
            request.RequestUri.ToString().Contains("page=0")
                ? JsonResponse(Array(FileJson("a"), FileJson("b")))
                : JsonResponse(Array(FileJson("c"))));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var files = await client.GetProjectFilesAsync("PXD012345", pageSize: 4);

        Assert.Multiple(() =>
        {
            Assert.That(files.Select(f => f.FileName), Is.EqualTo(new[] { "a", "b" }));
            Assert.That(handler.RequestedUris.Count, Is.EqualTo(1), "the short first page terminates the fetch");
        });
    }

    [Test]
    public void GetProjectFilesAsync_ServerIgnoresPaging_ThrowsAfterMaxPages()
    {
        // Every page is full and carries no total_records header, so neither the total nor the
        // short-page stop ever fires; the MaxPages cap must terminate it.
        var handler = new StubHandler(_ => JsonResponse(Array(FileJson("a"), FileJson("b"))));
        using var client = new PrideArchiveClient(new HttpClient(handler)) { MaxPages = 3 };

        Assert.That(async () => await client.GetProjectFilesAsync("PXD012345", pageSize: 2),
            Throws.InstanceOf<HttpRequestException>());
        Assert.That(handler.RequestedUris.Count, Is.EqualTo(3));
    }

    // ---- empty / edge -------------------------------------------------------

    [Test]
    public async Task GetProjectFilesAsync_UnknownOrEmptyProject_ReturnsEmptyList()
    {
        var handler = new StubHandler(_ => JsonResponse("[]"));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var files = await client.GetProjectFilesAsync("PXDBOGUS999");

        Assert.That(files, Is.Empty);
    }

    [Test]
    public async Task GetProjectFilesAsync_JsonNullBody_ReturnsEmptyList()
    {
        var handler = new StubHandler(_ => JsonResponse("null"));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var files = await client.GetProjectFilesAsync("PXD012345");

        Assert.That(files, Is.Empty);
    }

    [Test]
    public async Task GetProjectFilesAsync_ExplicitJsonNullField_DoesNotOverwriteDefaults()
    {
        var handler = new StubHandler(_ => JsonResponse("""[{ "fileName": "x.raw", "checksum": null }]"""));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var files = await client.GetProjectFilesAsync("PXD012345");

        Assert.That(files[0].Checksum, Is.EqualTo(string.Empty)); // NullValueHandling.Ignore keeps the default
        Assert.That(files[0].FileCategory, Is.Not.Null);
        Assert.That(files[0].PublicFileLocations, Is.Not.Null);
    }

    // ---- error contract -----------------------------------------------------

    [Test]
    public void GetProjectFilesAsync_BlankAccession_ThrowsArgumentException()
    {
        using var client = new PrideArchiveClient(new HttpClient(new StubHandler(_ => JsonResponse("[]"))));
        Assert.That(async () => await client.GetProjectFilesAsync("   "), Throws.ArgumentException);
    }

    [TestCase(0)]
    [TestCase(-1)]
    public void GetProjectFilesAsync_NonPositivePageSize_ThrowsArgumentOutOfRange(int pageSize)
    {
        using var client = new PrideArchiveClient(new HttpClient(new StubHandler(_ => JsonResponse("[]"))));
        Assert.That(async () => await client.GetProjectFilesAsync("PXD012345", pageSize), Throws.InstanceOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void GetProjectFilesAsync_NonSuccessStatus_ThrowsHttpRequestException()
    {
        var handler = new StubHandler(_ => JsonResponse("server error", HttpStatusCode.InternalServerError));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        Assert.That(async () => await client.GetProjectFilesAsync("PXD012345"), Throws.InstanceOf<HttpRequestException>());
    }

    // ---- construction / lifetime -------------------------------------------

    [Test]
    public void Constructor_NullHttpClient_ThrowsArgumentNullException()
    {
        Assert.That(() => new PrideArchiveClient(null), Throws.ArgumentNullException);
    }

    [Test]
    public void DefaultConstructor_And_Dispose_DoNotThrow()
    {
        var client = new PrideArchiveClient();
        Assert.That(() => client.Dispose(), Throws.Nothing);
    }

    [Test]
    public void Dispose_CalledTwice_DoesNotThrow()
    {
        var client = new PrideArchiveClient(new HttpClient(new StubHandler(_ => JsonResponse("[]"))));
        client.Dispose();
        Assert.That(() => client.Dispose(), Throws.Nothing); // idempotent
    }

    [Test]
    public void Dispose_InjectedHttpClient_IsNotDisposed_AndRemainsUsable()
    {
        var handler = new StubHandler(_ => JsonResponse("[]"));
        var httpClient = new HttpClient(handler) { BaseAddress = new Uri(PrideArchiveClient.DefaultBaseAddress) };
        var client = new PrideArchiveClient(httpClient); // caller retains ownership of httpClient

        client.Dispose();

        // the injected HttpClient must survive the client's Dispose (a disposed client would throw here)
        Assert.That(async () => await httpClient.GetAsync("projects/x/files"), Throws.Nothing);
        httpClient.Dispose();
    }

    [Test]
    public void Constructor_InjectedHttpClientWithoutBaseAddress_SetsPrideDefault()
    {
        var httpClient = new HttpClient(new StubHandler(_ => JsonResponse("[]"))); // no BaseAddress
        using var client = new PrideArchiveClient(httpClient);
        Assert.That(httpClient.BaseAddress, Is.EqualTo(new Uri(PrideArchiveClient.DefaultBaseAddress)));
    }

    [TestCase(0)]
    [TestCase(1)]
    public void GetProjectFilesAsync_ServerIgnoresPaging_LowMaxPages_Throws(int maxPages)
    {
        // every page is full and carries no total_records header, so only the MaxPages cap can stop it
        var handler = new StubHandler(_ => JsonResponse(Array(FileJson("a"), FileJson("b"))));
        using var client = new PrideArchiveClient(new HttpClient(handler)) { MaxPages = maxPages };
        Assert.That(async () => await client.GetProjectFilesAsync("PXD012345", pageSize: 2),
            Throws.InstanceOf<HttpRequestException>());
    }
}

/// <summary>
/// Live canary against the real PRIDE Archive API. Carries [Category("ExternalService")] so CI runs
/// it in the dedicated, non-blocking external-service job (see .github/workflows/dotnet.yml) rather
/// than the required unit-test run; [Category("Pride")] allows selecting it on its own. The call is
/// routed through <see cref="ExternalServiceTestHelper.RunAsync"/> so a PRIDE outage (timeout / 5xx /
/// unreachable) SKIPS the test, while a genuine contract break (wrong URL, unparseable response,
/// missing value) still FAILS. Run on demand — or in the external-service job — to detect API drift.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideArchiveClientLiveTests
{
    [Test]
    public Task GetProjectFilesAsync_LivePxd012345_ReturnsFullManifest() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient();
            var files = await client.GetProjectFilesAsync("PXD012345");
            Assert.That(files.Count, Is.GreaterThan(1));
            Assert.That(files.All(f => !string.IsNullOrEmpty(f.FileName)));
        });

    /// <summary>
    /// Asking for more files per page than PRIDE will serve must not shorten the manifest. PRIDE caps
    /// pageSize server-side (100 as of 2026-07-23) and then pages by the capped size, so an
    /// over-large request comes back "short" on every page while more records remain. Counts are
    /// compared against the default-page-size call rather than a literal, because the project's file
    /// count drifts; the invariant — page size must not change the answer — does not.
    /// </summary>
    [Test]
    public Task GetProjectFilesAsync_LivePageSizeAboveServerCap_ReturnsSameManifestAsDefault() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            // The cap PRIDE was observed to apply, and a request deliberately above it. Named rather
            // than inline so a future cap change is a one-line edit.
            const int observedServerCap = 100;
            const int deliberatelyAboveCap = 500;

            using var client = new PrideArchiveClient();
            var atDefault = await client.GetProjectFilesAsync("PXD012345");
            var aboveCap = await client.GetProjectFilesAsync("PXD012345", pageSize: deliberatelyAboveCap);

            // Skip when the over-cap paging path cannot be exercised at all: if even the default call
            // (which pages at the cap and so returns the project's true count) comes back at or below
            // the cap, a single page already holds everything and pageSize cannot change the answer, so
            // there is nothing to distinguish. This deliberately treats an exactly-at-cap manifest as
            // inconclusive rather than a hard failure -- a genuinely 100-file project and a truncated
            // one are indistinguishable from the count alone, and resolving that ambiguity to red would
            // make the canary go permanently red on dataset drift. A real truncation regression cannot
            // hide here: it shows up below as aboveCap disagreeing with atDefault, which fails red.
            Assume.That(atDefault.Count, Is.GreaterThan(observedServerCap),
                "PXD012345 does not exceed the server cap, so an over-cap request cannot be exercised.");

            // The real invariant, and the loud-failure path: page size must not change the answer. A
            // manifest truncated at the cap for the over-large request fails here (red), not skips.
            Assert.That(aboveCap.Select(f => f.FileName), Is.EquivalentTo(atDefault.Select(f => f.FileName)));
        });
}
