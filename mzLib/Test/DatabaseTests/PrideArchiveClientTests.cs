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

    [Test]
    public void GetProjectFilesAsync_NonPositivePageSize_ThrowsArgumentOutOfRange()
    {
        using var client = new PrideArchiveClient(new HttpClient(new StubHandler(_ => JsonResponse("[]"))));
        Assert.That(async () => await client.GetProjectFilesAsync("PXD012345", 0), Throws.InstanceOf<ArgumentOutOfRangeException>());
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
}
