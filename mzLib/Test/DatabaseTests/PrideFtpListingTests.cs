using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Threading;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

/// <summary>
/// Tests for <see cref="PrideArchiveClient.GetProjectFilesFromFtpAsync"/> — the complete file list read
/// by walking a project's FTP directory tree, for the cases where PRIDE's REST manifest is incomplete.
/// Offline: a stub handler serves a recorded project payload and Apache autoindex pages, so no network.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class PrideFtpListingTests
{
    /// <summary>An HttpMessageHandler that dispatches by request URI and records what was requested.</summary>
    private sealed class StubHandler : HttpMessageHandler
    {
        private readonly Func<string, HttpResponseMessage> _byUri;
        public List<string> RequestedUris { get; } = new();

        public StubHandler(Func<string, HttpResponseMessage> byUri) => _byUri = byUri;

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            string uri = request.RequestUri!.ToString();
            RequestedUris.Add(uri);
            return Task.FromResult(_byUri(uri));
        }
    }

    private static HttpResponseMessage Ok(string body) =>
        new(HttpStatusCode.OK) { Content = new StringContent(body) };

    private const string ProjectJson =
        """{ "accession": "PXD000001", "publicationDate": "2012-03-13" }""";

    // A standard Apache autoindex table: a sort-header row (href "?C=..."), a "Parent Directory" row
    // (absolute href "/..."), three files, and a subdirectory. Each entry's size is the LAST
    // right-aligned cell; the first right-aligned cell is the last-modified date, as PRIDE serves it.
    private const string RootIndexHtml = """
        <html><head><title>Index of /pride/data/archive/2012/03/PXD000001</title></head><body>
        <table>
        <tr><th><a href="?C=N;O=D">Name</a></th><th><a href="?C=M;O=A">Last modified</a></th><th><a href="?C=S;O=A">Size</a></th></tr>
        <tr><th colspan="5"><hr></th></tr>
        <tr><td><a href="/pride/data/archive/2012/03/">Parent Directory</a></td><td>&nbsp;</td><td align="right">  - </td></tr>
        <tr><td><a href="README.txt">README.txt</a></td><td align="right">2021-10-20 04:56  </td><td align="right">1.6K</td></tr>
        <tr><td><a href="run1.raw">run1.raw</a></td><td align="right">2021-10-20 04:56  </td><td align="right">210M</td></tr>
        <tr><td><a href="hidden_from_rest.mzML">hidden_from_rest.mzML</a></td><td align="right">2021-10-20 04:56  </td><td align="right">429M</td></tr>
        <tr><td><a href="generated/">generated/</a></td><td align="right">2021-10-20 04:56  </td><td align="right">  - </td></tr>
        </table></body></html>
        """;

    private const string GeneratedIndexHtml = """
        <html><head><title>Index of /pride/data/archive/2012/03/PXD000001/generated</title></head><body>
        <table>
        <tr><th><a href="?C=N;O=D">Name</a></th></tr>
        <tr><td><a href="/pride/data/archive/2012/03/PXD000001/">Parent Directory</a></td><td align="right">  - </td></tr>
        <tr><td><a href="summary.mztab">summary.mztab</a></td><td align="right">2021-10-20 04:56  </td><td align="right">844K</td></tr>
        </table></body></html>
        """;

    private const string RootUrl = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/";
    private const string GeneratedUrl = RootUrl + "generated/";

    /// <summary>A client that serves the recorded project + autoindex pages, and 404s an unknown project.</summary>
    private static PrideArchiveClient StubbedClient(out StubHandler handler)
    {
        handler = new StubHandler(uri =>
        {
            if (uri.EndsWith("/projects/PXD000001", StringComparison.Ordinal)) return Ok(ProjectJson);
            if (uri.Contains("/projects/")) return new HttpResponseMessage(HttpStatusCode.NotFound);
            if (uri == RootUrl) return Ok(RootIndexHtml);
            if (uri == GeneratedUrl) return Ok(GeneratedIndexHtml);
            return new HttpResponseMessage(HttpStatusCode.NotFound);
        });
        return new PrideArchiveClient(new HttpClient(handler));
    }

    [Test]
    public async Task WalksTheWholeTreeAndReportsFilesTheRestManifestOmits()
    {
        using PrideArchiveClient client = StubbedClient(out _);

        List<PrideFtpFile> files = await client.GetProjectFilesFromFtpAsync("PXD000001");
        Dictionary<string, PrideFtpFile> byPath = files.ToDictionary(f => f.RelativePath);

        Assert.Multiple(() =>
        {
            // The complete set: three top-level files plus the one inside generated/. The "Parent
            // Directory" link, the column-sort links, and the subdirectory row itself are NOT files.
            Assert.That(byPath.Keys, Is.EquivalentTo(new[]
            {
                "README.txt", "run1.raw", "hidden_from_rest.mzML", "generated/summary.mztab",
            }));

            // The whole reason this method exists: a file the REST manifest hides is present here.
            Assert.That(byPath.ContainsKey("hidden_from_rest.mzML"), Is.True);

            // Recursion: a nested file carries its subdirectory in the relative path but a bare leaf name.
            PrideFtpFile nested = byPath["generated/summary.mztab"];
            Assert.That(nested.FileName, Is.EqualTo("summary.mztab"));
            Assert.That(nested.Url, Is.EqualTo(GeneratedUrl + "summary.mztab"));

            // The URL is the HTTPS location a caller can download from.
            Assert.That(byPath["run1.raw"].Url, Is.EqualTo(RootUrl + "run1.raw"));
        });
    }

    [Test]
    public async Task ApproximateSizesAreParsedFromTheIndexSuffixes()
    {
        using PrideArchiveClient client = StubbedClient(out _);

        Dictionary<string, PrideFtpFile> byPath =
            (await client.GetProjectFilesFromFtpAsync("PXD000001")).ToDictionary(f => f.RelativePath);

        Assert.Multiple(() =>
        {
            Assert.That(byPath["README.txt"].ApproximateSizeBytes, Is.EqualTo((long)Math.Round(1.6 * 1024)));
            Assert.That(byPath["run1.raw"].ApproximateSizeBytes, Is.EqualTo(210L * 1024 * 1024));
            Assert.That(byPath["hidden_from_rest.mzML"].ApproximateSizeBytes, Is.EqualTo(429L * 1024 * 1024));
        });
    }

    [Test]
    public void AnUnknownAccessionIsAnMzLibExceptionNotAnOutage()
    {
        using PrideArchiveClient client = StubbedClient(out _);

        // GetProjectAsync answers a 404 project with MzLibException (not HttpRequestException, which a
        // live test would skip). Walking the FTP tree must inherit that contract.
        Assert.ThrowsAsync<MzLibException>(async () =>
            await client.GetProjectFilesFromFtpAsync("PXD999999"));
    }

    [Test]
    public void ABlankAccessionIsAnArgumentExceptionBeforeAnyRequest()
    {
        using PrideArchiveClient client = StubbedClient(out StubHandler handler);

        Assert.ThrowsAsync<ArgumentException>(async () => await client.GetProjectFilesFromFtpAsync("   "));
        Assert.That(handler.RequestedUris, Is.Empty, "a blank accession must be rejected before any HTTP request");
    }

    [Test]
    public void AFailedDirectoryListingIsAnHttpRequestException()
    {
        // Project resolves, but its FTP directory returns non-success: that is an outage-shaped failure.
        var handler = new StubHandler(uri =>
            uri.EndsWith("/projects/PXD000001", StringComparison.Ordinal)
                ? Ok(ProjectJson)
                : new HttpResponseMessage(HttpStatusCode.ServiceUnavailable));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        Assert.ThrowsAsync<HttpRequestException>(async () =>
            await client.GetProjectFilesFromFtpAsync("PXD000001"));
    }
}

/// <summary>
/// Live canary against the real PRIDE FTP tree. [Category("ExternalService")] keeps it out of the
/// required unit run (see the sibling <see cref="PrideArchiveClientLiveTests"/>);
/// <see cref="ExternalServiceTestHelper.RunAsync"/> SKIPS it on a PRIDE outage but lets a genuine
/// format drift FAIL. This is what catches EBI changing its Apache autoindex layout — the one thing
/// the offline fixtures cannot.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideFtpListingLiveTests
{
    [Test]
    public Task GetProjectFilesFromFtp_LivePxd000001_IsMoreCompleteThanTheRestManifest() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient();

            List<PrideFtpFile> ftpFiles = await client.GetProjectFilesFromFtpAsync("PXD000001");
            List<PrideArchiveFile> restFiles = await client.GetProjectFilesAsync("PXD000001");

            // The whole point: the FTP tree is more complete than the REST manifest PRIDE serves.
            // Asserted as a comparison rather than a hard-coded 13, so it survives PRIDE re-curating the
            // dataset while still failing if the autoindex parse breaks (which would yield zero files).
            Assert.That(ftpFiles.Count, Is.GreaterThan(restFiles.Count),
                "the FTP tree should hold more files than the REST manifest");
            Assert.That(ftpFiles, Is.Not.Empty);
            Assert.That(ftpFiles.All(f => !string.IsNullOrEmpty(f.RelativePath) && !string.IsNullOrEmpty(f.Url)));
            Assert.That(ftpFiles.Any(f => f.ApproximateSizeBytes > 0), "at least one file should have a parsed size");
        });
}
