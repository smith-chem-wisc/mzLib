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

    /// <summary>A minimal Apache-autoindex page wrapping the given entry rows.</summary>
    private static string Index(params string[] rows) =>
        "<html><body><table>\n" +
        """<tr><th><a href="?C=N;O=D">Name</a></th></tr>""" + "\n" +
        string.Join("\n", rows) + "\n</table></body></html>";

    /// <summary>One autoindex row for an entry (a trailing "/" in <paramref name="href"/> makes it a directory).</summary>
    private static string Row(string href, string size) =>
        $"""<tr><td><a href="{href}">{href}</a></td><td align="right">2021-10-20 04:56  </td><td align="right">{size}</td></tr>""";

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

    [Test]
    public void ASelfReferentialSubdirectoryHitsTheDepthCapInsteadOfRecursingForever()
    {
        // A listing whose only entry links to a deeper copy of itself grows the URL each level, so the
        // visited-set cannot prune it — the depth cap is the backstop. It throws MzLibException (a broken
        // listing is a contract violation, NOT an HttpRequestException a live test would skip as an
        // outage), and the walk stays bounded to MaxFtpDirectoryDepth (64) levels plus the root = 65.
        int dirFetches = 0;
        var handler = new StubHandler(uri =>
        {
            if (uri.EndsWith("/projects/PXD000001", StringComparison.Ordinal)) return Ok(ProjectJson);
            if (uri.Contains("/projects/")) return new HttpResponseMessage(HttpStatusCode.NotFound);
            dirFetches++;
            return Ok(Index(Row("loop/", "-")));
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var ex = Assert.ThrowsAsync<MzLibException>(async () =>
            await client.GetProjectFilesFromFtpAsync("PXD000001"));
        Assert.Multiple(() =>
        {
            Assert.That(ex!.Message, Does.Contain("nesting exceeded"));
            Assert.That(dirFetches, Is.EqualTo(65), "the walk must stop at the depth cap, not run away");
        });
    }

    [Test]
    public async Task ARelativeParentLinkIsRefusedSoTheWalkStaysInsideTheProject()
    {
        // A "../" would resolve (via Uri normalisation) to the project's PARENT directory and pull in
        // unrelated files. It must be skipped; only plain child directories are followed.
        const string parentUrl = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/";
        var handler = new StubHandler(uri =>
        {
            if (uri.EndsWith("/projects/PXD000001", StringComparison.Ordinal)) return Ok(ProjectJson);
            if (uri == RootUrl) return Ok(Index(Row("../", "-"), Row("keep.raw", "1.0K")));
            if (uri == parentUrl) return Ok(Index(Row("escaped.raw", "9.9M")));
            return new HttpResponseMessage(HttpStatusCode.NotFound);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        List<PrideFtpFile> files = await client.GetProjectFilesFromFtpAsync("PXD000001");

        Assert.Multiple(() =>
        {
            Assert.That(handler.RequestedUris, Does.Not.Contain(parentUrl),
                "a '../' link must not be followed out of the project subtree");
            Assert.That(files.Select(f => f.RelativePath), Is.EquivalentTo(new[] { "keep.raw" }));
        });
    }

    [Test]
    public async Task PercentEncodedNamesAreDecodedWhileTheUrlStaysEncoded()
    {
        var handler = new StubHandler(uri =>
        {
            if (uri.EndsWith("/projects/PXD000001", StringComparison.Ordinal)) return Ok(ProjectJson);
            if (uri == RootUrl) return Ok(Index(Row("my%20run%20(1).raw", "1.0K")));
            return new HttpResponseMessage(HttpStatusCode.NotFound);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        PrideFtpFile file = (await client.GetProjectFilesFromFtpAsync("PXD000001")).Single();

        Assert.Multiple(() =>
        {
            // The human-readable name is percent-decoded; the URL keeps the encoding so it stays a
            // valid, downloadable path.
            Assert.That(file.FileName, Is.EqualTo("my run (1).raw"));
            Assert.That(file.RelativePath, Is.EqualTo("my run (1).raw"));
            Assert.That(file.Url, Is.EqualTo(RootUrl + "my%20run%20(1).raw"));
        });
    }

    [Test]
    public async Task TotalApproximateSizeBytesSumsTheWholeTree()
    {
        using PrideArchiveClient client = StubbedClient(out _);
        List<PrideFtpFile> files = await client.GetProjectFilesFromFtpAsync("PXD000001");

        // Pinned to the literal expected total (README 1.6K + run1 210M + hidden 429M + the nested
        // generated/summary 844K), so the assertion can fail on a real regression rather than
        // re-deriving the implementation's own Sum. Proves the helper sums the WHOLE tree, subdir included.
        const long expected = 1638L + 210L * 1024 * 1024 + 429L * 1024 * 1024 + 864256L;
        Assert.That(files.TotalApproximateSizeBytes(), Is.EqualTo(expected));
    }

    [Test]
    public async Task ADescriptionColumnDoesNotStealTheSizeCell()
    {
        // A stock Apache template also right-aligns the (empty) Description column, so the LAST
        // right-aligned cell is not the size. The parser picks the first cell that parses as a size.
        const string rowWithDescription =
            """<tr><td><a href="d.raw">d.raw</a></td><td align="right">2021-10-20 04:56  </td><td align="right">2.0M</td><td align="right">&nbsp;</td></tr>""";
        var handler = new StubHandler(uri =>
        {
            if (uri.EndsWith("/projects/PXD000001", StringComparison.Ordinal)) return Ok(ProjectJson);
            if (uri == RootUrl) return Ok(Index(rowWithDescription));
            return new HttpResponseMessage(HttpStatusCode.NotFound);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        PrideFtpFile file = (await client.GetProjectFilesFromFtpAsync("PXD000001")).Single();
        Assert.That(file.ApproximateSizeBytes, Is.EqualTo(2L * 1024 * 1024),
            "the empty Description column must not shadow the Size cell");
    }

    [Test]
    public async Task HtmlEntitiesInNamesAreDecoded()
    {
        var handler = new StubHandler(uri =>
        {
            if (uri.EndsWith("/projects/PXD000001", StringComparison.Ordinal)) return Ok(ProjectJson);
            if (uri == RootUrl) return Ok(Index(Row("a&amp;b.raw", "1.0K")));
            return new HttpResponseMessage(HttpStatusCode.NotFound);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        PrideFtpFile file = (await client.GetProjectFilesFromFtpAsync("PXD000001")).Single();
        Assert.That(file.FileName, Is.EqualTo("a&b.raw"), "an HTML entity in the href must decode for the name");
    }

    [Test]
    public async Task SizeSuffixesGAndBareByteCountsParse()
    {
        var handler = new StubHandler(uri =>
        {
            if (uri.EndsWith("/projects/PXD000001", StringComparison.Ordinal)) return Ok(ProjectJson);
            if (uri == RootUrl) return Ok(Index(Row("big.raw", "1.2G"), Row("tiny.txt", "512")));
            return new HttpResponseMessage(HttpStatusCode.NotFound);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        Dictionary<string, PrideFtpFile> byPath =
            (await client.GetProjectFilesFromFtpAsync("PXD000001")).ToDictionary(f => f.RelativePath);

        Assert.Multiple(() =>
        {
            Assert.That(byPath["big.raw"].ApproximateSizeBytes, Is.EqualTo((long)Math.Round(1.2 * 1024 * 1024 * 1024)));
            Assert.That(byPath["tiny.txt"].ApproximateSizeBytes, Is.EqualTo(512L), "a bare byte count carries no suffix");
        });
    }

    [Test]
    public void APreCancelledTokenStopsTheWalk()
    {
        using PrideArchiveClient client = StubbedClient(out _);
        using var cts = new CancellationTokenSource();
        cts.Cancel();

        // CatchAsync, not ThrowsAsync: the concrete type is TaskCanceledException, which derives from
        // OperationCanceledException but would fail an exact-type ThrowsAsync match.
        Assert.CatchAsync<OperationCanceledException>(async () =>
            await client.GetProjectFilesFromFtpAsync("PXD000001", cts.Token));
    }

    [Test]
    public void AProjectWithNoPublicationDateIsAnMzLibException()
    {
        // Without a publication date there is no year/month to locate the FTP directory. That is a
        // broken contract (MzLibException), reported before any FTP request is attempted.
        const string noDateJson = """{ "accession": "PXD000001" }""";
        var handler = new StubHandler(uri =>
            uri.EndsWith("/projects/PXD000001", StringComparison.Ordinal)
                ? Ok(noDateJson)
                : new HttpResponseMessage(HttpStatusCode.NotFound));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var ex = Assert.ThrowsAsync<MzLibException>(async () =>
            await client.GetProjectFilesFromFtpAsync("PXD000001"));
        Assert.That(ex!.Message, Does.Contain("no publication date"));
    }

    [Test]
    public void TotalApproximateSizeBytesRejectsANullSequence()
    {
        Assert.Throws<ArgumentNullException>(() => ((IEnumerable<PrideFtpFile>)null).TotalApproximateSizeBytes());
    }

    [Test]
    public async Task ATSuffixParsesAndAnUnrecognisedSizeCellIsZeroNotAFailure()
    {
        var handler = new StubHandler(uri =>
        {
            if (uri.EndsWith("/projects/PXD000001", StringComparison.Ordinal)) return Ok(ProjectJson);
            if (uri == RootUrl) return Ok(Index(Row("huge.raw", "1.5T"), Row("weird.dat", "N/A")));
            return new HttpResponseMessage(HttpStatusCode.NotFound);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        Dictionary<string, PrideFtpFile> byPath =
            (await client.GetProjectFilesFromFtpAsync("PXD000001")).ToDictionary(f => f.RelativePath);

        Assert.Multiple(() =>
        {
            Assert.That(byPath["huge.raw"].ApproximateSizeBytes,
                Is.EqualTo((long)Math.Round(1.5 * 1024 * 1024 * 1024 * 1024)));
            // A cell that is not a recognisable size leaves the file listed with an unknown (0) size
            // rather than dropping it or throwing — the list stays complete.
            Assert.That(byPath["weird.dat"].ApproximateSizeBytes, Is.EqualTo(0L));
        });
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

            HashSet<string> ftpNames = ftpFiles.Select(f => f.FileName).ToHashSet();

            // The whole point: the FTP tree is a strict SUPERSET of the REST manifest — it contains
            // every file REST reports AND more. Asserted structurally (not a hard-coded 13), so it
            // survives PRIDE re-curating the dataset while still failing if the parse breaks (zero files)
            // or ever drops a file REST does list.
            Assert.That(restFiles.Select(r => r.FileName), Is.SubsetOf(ftpNames),
                "the FTP tree should contain every file the REST manifest lists");
            Assert.That(ftpFiles.Count, Is.GreaterThan(restFiles.Count),
                $"FTP tree ({ftpFiles.Count}) should hold more files than the REST manifest ({restFiles.Count})");
            Assert.That(ftpFiles, Is.Not.Empty);
            Assert.That(ftpFiles.All(f => !string.IsNullOrEmpty(f.RelativePath) && !string.IsNullOrEmpty(f.Url)));
            Assert.That(ftpFiles.Any(f => f.ApproximateSizeBytes > 0), "at least one file should have a parsed size");
        });
}
