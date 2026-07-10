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
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

[TestFixture]
[ExcludeFromCodeCoverage]
public class PrideArchiveDownloadTests
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
            cancellationToken.ThrowIfCancellationRequested(); // honor cancellation like a real handler
            RequestedUris.Add(request.RequestUri.ToString());
            return Task.FromResult(_responder(request));
        }
    }

    /// <summary>A read-only stream that yields a few bytes and then throws, to simulate a mid-transfer failure.</summary>
    private sealed class ThrowingStream : Stream
    {
        private int _bytesBeforeThrow;
        public ThrowingStream(int bytesBeforeThrow) => _bytesBeforeThrow = bytesBeforeThrow;

        public override int Read(byte[] buffer, int offset, int count)
        {
            if (_bytesBeforeThrow <= 0)
                throw new IOException("simulated mid-stream failure");
            int n = Math.Min(count, _bytesBeforeThrow);
            for (int i = 0; i < n; i++) buffer[offset + i] = 0x41;
            _bytesBeforeThrow -= n;
            return n;
        }

        public override bool CanRead => true;
        public override bool CanSeek => false;
        public override bool CanWrite => false;
        public override long Length => throw new NotSupportedException();
        public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }
        public override void Flush() { }
        public override long Seek(long offset, SeekOrigin origin) => throw new NotSupportedException();
        public override void SetLength(long value) => throw new NotSupportedException();
        public override void Write(byte[] buffer, int offset, int count) => throw new NotSupportedException();
    }

    private static HttpResponseMessage Bytes(byte[] body, HttpStatusCode status = HttpStatusCode.OK) =>
        new(status) { Content = new ByteArrayContent(body) };

    private static HttpResponseMessage Json(string body, HttpStatusCode status = HttpStatusCode.OK) =>
        new(status) { Content = new StringContent(body, Encoding.UTF8, "application/json") };

    /// <summary>Builds a file DTO with the given category and public locations (accession, value) pairs.</summary>
    private static PrideArchiveFile MakeFile(string fileName, string category = "RAW",
        params (string accession, string value)[] locations) =>
        new()
        {
            FileName = fileName,
            FileSizeBytes = 10,
            FileCategory = new CvParam("PRIDE", "PRIDE:0000404", "category", category),
            PublicFileLocations = locations
                .Select(l => new CvParam("PRIDE", l.accession, "location", l.value))
                .ToList()
        };

    private static (string accession, string value) Ftp(string path) =>
        (PrideArchiveExtensions.FtpLocationAccession, "ftp://ftp.pride.ebi.ac.uk/" + path);

    private static (string accession, string value) Aspera(string path) =>
        (PrideArchiveExtensions.AsperaLocationAccession, "prd_ascp@fasp.ebi.ac.uk:" + path);

    private string _tempDir;

    [SetUp]
    public void SetUp()
    {
        _tempDir = Path.Combine(Path.GetTempPath(), "PrideArchiveDownloadTests", Guid.NewGuid().ToString("N"));
    }

    [TearDown]
    public void TearDown()
    {
        try { if (Directory.Exists(_tempDir)) Directory.Delete(_tempDir, recursive: true); }
        catch { /* best-effort cleanup */ }
    }

    // ---- URL resolution (pure) ---------------------------------------------

    [Test]
    public void TryGetHttpsDownloadUrl_UpgradesFtpLocationToHttps()
    {
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"), Aspera("pride/data/x/run1.raw"));

        Assert.That(file.TryGetHttpsDownloadUrl(out string url), Is.True);
        Assert.That(url, Is.EqualTo("https://ftp.pride.ebi.ac.uk/pride/data/x/run1.raw"));
    }

    [Test]
    public void TryGetHttpsDownloadUrl_PrefersExplicitHttpsOverFtp()
    {
        var file = MakeFile("run1.raw", "RAW",
            Ftp("pride/data/x/run1.raw"),
            ("PRIDE:0000000", "https://example.org/direct/run1.raw"));

        Assert.That(file.GetHttpsDownloadUrl(), Is.EqualTo("https://example.org/direct/run1.raw"));
    }

    [Test]
    public void GetHttpsDownloadUrl_AsperaOnly_ThrowsNotSupported()
    {
        var file = MakeFile("run1.raw", "RAW", Aspera("pride/data/x/run1.raw"));

        Assert.That(file.TryGetHttpsDownloadUrl(out _), Is.False);
        Assert.That(() => file.GetHttpsDownloadUrl(), Throws.InstanceOf<NotSupportedException>());
    }

    [Test]
    public void TryGetHttpsDownloadUrl_NoLocations_ReturnsFalse()
    {
        var file = MakeFile("run1.raw"); // no locations

        Assert.That(file.TryGetHttpsDownloadUrl(out string url), Is.False);
        Assert.That(url, Is.Null);
    }

    [Test]
    public void TryGetHttpsDownloadUrl_NullFile_Throws()
    {
        PrideArchiveFile file = null;
        Assert.That(() => file.TryGetHttpsDownloadUrl(out _), Throws.ArgumentNullException);
    }

    // ---- manifest filters (pure) -------------------------------------------

    [Test]
    public void WhereCategory_FiltersCaseInsensitively()
    {
        var files = new[]
        {
            MakeFile("a.raw", "RAW"),
            MakeFile("b.mzid", "SEARCH"),
            MakeFile("c.raw", "raw"),
        };

        var raws = files.WhereCategory("RAW").Select(f => f.FileName).ToArray();
        Assert.That(raws, Is.EqualTo(new[] { "a.raw", "c.raw" }));
    }

    [Test]
    public void WhereExtension_MatchesWithOrWithoutLeadingDot()
    {
        var files = new[] { MakeFile("a.raw"), MakeFile("b.RAW"), MakeFile("c.mzML") };

        Assert.That(files.WhereExtension("raw").Select(f => f.FileName), Is.EqualTo(new[] { "a.raw", "b.RAW" }));
        Assert.That(files.WhereExtension(".mzml").Select(f => f.FileName), Is.EqualTo(new[] { "c.mzML" }));
        Assert.That(files.WhereExtension("raw", "mzml").Count(), Is.EqualTo(3));
    }

    [Test]
    public void TotalSizeBytes_SumsAllFiles()
    {
        var files = new[] { MakeFile("a.raw"), MakeFile("b.raw"), MakeFile("c.raw") }; // 10 bytes each
        Assert.That(files.TotalSizeBytes(), Is.EqualTo(30));
    }

    [Test]
    public void Filters_NullOrBlankArguments_Throw()
    {
        IEnumerable<PrideArchiveFile> nil = null;
        Assert.That(() => nil.WhereCategory("RAW").ToArray(), Throws.ArgumentNullException);
        Assert.That(() => nil.WhereExtension("raw").ToArray(), Throws.ArgumentNullException);
        Assert.That(() => nil.TotalSizeBytes(), Throws.ArgumentNullException);
        Assert.That(() => new[] { MakeFile("a.raw") }.WhereCategory(" ").ToArray(), Throws.ArgumentException);
        Assert.That(() => new[] { MakeFile("a.raw") }.WhereExtension().ToArray(), Throws.ArgumentException);
    }

    // ---- DownloadFileAsync (offline) ---------------------------------------

    [Test]
    public async Task DownloadFileAsync_UpgradesFtp_WritesBytes_ReturnsPath()
    {
        byte[] payload = Encoding.UTF8.GetBytes("raw-file-contents");
        var handler = new StubHandler(_ => Bytes(payload));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        string path = await client.DownloadFileAsync(file, _tempDir);

        Assert.That(handler.RequestedUris.Single(), Is.EqualTo("https://ftp.pride.ebi.ac.uk/pride/data/x/run1.raw"));
        Assert.That(path, Is.EqualTo(Path.Combine(_tempDir, "run1.raw")));
        Assert.That(File.ReadAllBytes(path), Is.EqualTo(payload));
        Assert.That(File.Exists(path + ".partial"), Is.False); // temp file cleaned up
    }

    [Test]
    public async Task DownloadFileAsync_NotOverwrite_ExistingFile_SkipsRequest()
    {
        Directory.CreateDirectory(_tempDir);
        string existing = Path.Combine(_tempDir, "run1.raw");
        File.WriteAllText(existing, "old");
        var handler = new StubHandler(_ => Bytes(Encoding.UTF8.GetBytes("new")));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        string path = await client.DownloadFileAsync(file, _tempDir, overwrite: false);

        Assert.That(handler.RequestedUris, Is.Empty);              // no network call
        Assert.That(File.ReadAllText(path), Is.EqualTo("old"));    // untouched
    }

    [Test]
    public async Task DownloadFileAsync_Overwrite_ReplacesExistingFile()
    {
        Directory.CreateDirectory(_tempDir);
        string existing = Path.Combine(_tempDir, "run1.raw");
        File.WriteAllText(existing, "old");
        var handler = new StubHandler(_ => Bytes(Encoding.UTF8.GetBytes("new")));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        string path = await client.DownloadFileAsync(file, _tempDir); // overwrite defaults true

        Assert.That(File.ReadAllText(path), Is.EqualTo("new"));
    }

    [Test]
    public void DownloadFileAsync_NonSuccess_Throws_AndLeavesNoPartial()
    {
        var handler = new StubHandler(_ => Bytes(Encoding.UTF8.GetBytes("not found"), HttpStatusCode.NotFound));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        Assert.That(async () => await client.DownloadFileAsync(file, _tempDir), Throws.InstanceOf<HttpRequestException>());
        Assert.That(File.Exists(Path.Combine(_tempDir, "run1.raw.partial")), Is.False);
        Assert.That(File.Exists(Path.Combine(_tempDir, "run1.raw")), Is.False);
    }

    [Test]
    public void DownloadFileAsync_AsperaOnly_ThrowsNotSupported()
    {
        var handler = new StubHandler(_ => Bytes(Array.Empty<byte>()));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        var file = MakeFile("run1.raw", "RAW", Aspera("pride/data/x/run1.raw"));

        Assert.That(async () => await client.DownloadFileAsync(file, _tempDir), Throws.InstanceOf<NotSupportedException>());
        Assert.That(handler.RequestedUris, Is.Empty);
    }

    [Test]
    public void DownloadFileAsync_InvalidArguments_Throw()
    {
        using var client = new PrideArchiveClient(new HttpClient(new StubHandler(_ => Bytes(Array.Empty<byte>()))));
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        Assert.That(async () => await client.DownloadFileAsync(null, _tempDir), Throws.ArgumentNullException);
        Assert.That(async () => await client.DownloadFileAsync(file, "  "), Throws.ArgumentException);
        Assert.That(async () => await client.DownloadFileAsync(MakeFile("", "RAW", Ftp("x")), _tempDir), Throws.ArgumentException);
    }

    // ---- DownloadProjectFilesAsync (offline) -------------------------------

    [Test]
    public async Task DownloadProjectFilesAsync_FiltersByCategory_DownloadsSelectedOnly()
    {
        string manifest = "[" + string.Join(",", new[]
        {
            FileJson("keep.raw", "RAW", "PRIDE:0000404"),
            FileJson("drop.mzid", "SEARCH", "PRIDE:0000408"),
        }) + "]";

        var handler = new StubHandler(request =>
        {
            string uri = request.RequestUri.ToString();
            if (uri.Contains("/files")) return Json(manifest);           // manifest request
            return Bytes(Encoding.UTF8.GetBytes("bytes-of-" + uri));      // file download
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var paths = await client.DownloadProjectFilesAsync("PXD012345", _tempDir, f => f.FileCategory.Value == "RAW");

        Assert.That(paths.Count, Is.EqualTo(1));
        Assert.That(Path.GetFileName(paths[0]), Is.EqualTo("keep.raw"));
        Assert.That(File.Exists(Path.Combine(_tempDir, "keep.raw")), Is.True);
        Assert.That(File.Exists(Path.Combine(_tempDir, "drop.mzid")), Is.False);
        // the bytes written must be the payload served for the keep.raw download URL, not some other response
        string keepUrl = "https://ftp.pride.ebi.ac.uk/pride/data/x/keep.raw";
        Assert.That(File.ReadAllText(paths[0]), Is.EqualTo("bytes-of-" + keepUrl));
        // one manifest request + one download (the SEARCH file was filtered out before any download)
        Assert.That(handler.RequestedUris.Count(u => u.Contains("ftp.pride.ebi.ac.uk")), Is.EqualTo(1));
    }

    [Test]
    public void DownloadProjectFilesAsync_BlankDestination_Throws()
    {
        using var client = new PrideArchiveClient(new HttpClient(new StubHandler(_ => Json("[]"))));
        Assert.That(async () => await client.DownloadProjectFilesAsync("PXD012345", " "), Throws.ArgumentException);
    }

    // ---- URL resolution: host restriction and shadowing (regressions) ------

    [Test]
    public void TryGetHttpsDownloadUrl_FtpOnNonPrideHost_ReturnsFalse()
    {
        // REGRESSION: failed before the fix — any ftp:// host was upgraded to https://, producing a broken URL.
        var file = MakeFile("run1.raw", "RAW", (PrideArchiveExtensions.FtpLocationAccession, "ftp://ftp.example.org/x/run1.raw"));

        Assert.That(file.TryGetHttpsDownloadUrl(out string url), Is.False);
        Assert.That(url, Is.Null);
        Assert.That(() => file.GetHttpsDownloadUrl(), Throws.InstanceOf<NotSupportedException>());
    }

    [Test]
    public void TryGetHttpsDownloadUrl_MalformedFtpAccessionEntry_DoesNotShadowValidFtp()
    {
        // REGRESSION: failed before the fix — an FTP-accession entry with a null value shadowed the valid ftp:// URL.
        var file = new PrideArchiveFile
        {
            FileName = "run1.raw",
            FileSizeBytes = 10,
            FileCategory = new CvParam("PRIDE", "PRIDE:0000404", "category", "RAW"),
            PublicFileLocations = new List<CvParam>
            {
                null,                                                                                 // a null location entry must be skipped, not throw
                new("PRIDE", PrideArchiveExtensions.FtpLocationAccession, "location", null),           // malformed: no value
                new("PRIDE", PrideArchiveExtensions.FtpLocationAccession, "location", "ftp://ftp.pride.ebi.ac.uk/pride/data/x/run1.raw"),
            }
        };

        Assert.That(file.TryGetHttpsDownloadUrl(out string url), Is.True);
        Assert.That(url, Is.EqualTo("https://ftp.pride.ebi.ac.uk/pride/data/x/run1.raw"));
    }

    // ---- WhereExtension edge cases (regressions) ---------------------------

    [Test]
    public void WhereExtension_OnlyBlankExtensions_Throws()
    {
        // REGRESSION: failed before the fix — all-blank extensions silently matched zero files instead of throwing.
        var files = new[] { MakeFile("a.raw") };
        Assert.That(() => files.WhereExtension(" ", "").ToArray(), Throws.ArgumentException);
    }

    [Test]
    public void WhereExtension_MixOfBlankAndReal_FiltersOnRealOnly()
    {
        var files = new[] { MakeFile("a.raw"), MakeFile("b.mzML") };
        Assert.That(files.WhereExtension("", ".raw").Select(f => f.FileName), Is.EqualTo(new[] { "a.raw" }));
    }

    [Test]
    public void WhereExtension_ExtensionlessAndDoubleExtension_Handled()
    {
        var files = new PrideArchiveFile[] { null, MakeFile("README"), MakeFile("a.mzML.gz"), MakeFile("b.gz") };

        // a null file is skipped; an extensionless name is never a match; a double-extension name matches only its trailing extension
        Assert.That(files.WhereExtension(".gz").Select(f => f.FileName), Is.EqualTo(new[] { "a.mzML.gz", "b.gz" }));
        Assert.That(files.WhereExtension(".mzML").Select(f => f.FileName), Is.Empty);
    }

    [Test]
    public void WhereCategory_NullCategory_Skipped()
    {
        var withNullCategory = new PrideArchiveFile { FileName = "x.raw", FileCategory = null };
        var files = new PrideArchiveFile[] { null, withNullCategory, MakeFile("y.raw", "RAW") };

        // a null file, or one whose category is null (possible from the wire), is skipped, not an NRE
        Assert.That(() => files.WhereCategory("RAW").ToArray(), Throws.Nothing);
        Assert.That(files.WhereCategory("RAW").Select(f => f.FileName), Is.EqualTo(new[] { "y.raw" }));
    }

    // ---- DownloadFileAsync: path traversal, cancellation, mid-stream fault --

    [TestCase("../evil.raw")]
    [TestCase("sub/dir/file.raw")]
    public void DownloadFileAsync_FileNameNotBareLeaf_Throws_WritesNothingOutside(string fileName)
    {
        // REGRESSION: failed before the fix — a path-bearing FileName escaped destinationDirectory via Path.Combine.
        var handler = new StubHandler(_ => Bytes(Encoding.UTF8.GetBytes("x")));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        var file = MakeFile(fileName, "RAW", Ftp("pride/data/x/leaf.raw"));

        Assert.That(async () => await client.DownloadFileAsync(file, _tempDir), Throws.ArgumentException);
        Assert.That(handler.RequestedUris, Is.Empty);                        // rejected before any network call
        Assert.That(File.Exists(Path.Combine(_tempDir, "evil.raw")), Is.False);
    }

    [Test]
    public void DownloadFileAsync_PreCancelledToken_Throws_LeavesNoFile()
    {
        var handler = new StubHandler(_ => Bytes(Encoding.UTF8.GetBytes("data")));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));
        using var cts = new CancellationTokenSource();
        cts.Cancel();

        Assert.That(async () => await client.DownloadFileAsync(file, _tempDir, cancellationToken: cts.Token),
            Throws.InstanceOf<OperationCanceledException>());
        Assert.That(File.Exists(Path.Combine(_tempDir, "run1.raw")), Is.False);
        Assert.That(File.Exists(Path.Combine(_tempDir, "run1.raw.partial")), Is.False);
    }

    [Test]
    public void DownloadFileAsync_MidStreamFailure_LeavesNoPartial()
    {
        // the ".partial"-then-move design must clean up the sibling when the transfer faults mid-stream
        var handler = new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.OK)
        {
            Content = new StreamContent(new ThrowingStream(bytesBeforeThrow: 4))
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        Assert.That(async () => await client.DownloadFileAsync(file, _tempDir), Throws.InstanceOf<IOException>());
        Assert.That(File.Exists(Path.Combine(_tempDir, "run1.raw.partial")), Is.False);
        Assert.That(File.Exists(Path.Combine(_tempDir, "run1.raw")), Is.False);
    }

    // ---- DownloadProjectFilesAsync: null filter and empty selection --------

    [Test]
    public async Task DownloadProjectFilesAsync_NullFilter_DownloadsAllFiles()
    {
        string manifest = "[" + string.Join(",", new[]
        {
            FileJson("one.raw", "RAW", "PRIDE:0000404"),
            FileJson("two.mzid", "SEARCH", "PRIDE:0000408"),
        }) + "]";
        var handler = new StubHandler(request =>
        {
            string uri = request.RequestUri.ToString();
            if (uri.Contains("/files")) return Json(manifest);
            return Bytes(Encoding.UTF8.GetBytes("bytes-of-" + uri));
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var paths = await client.DownloadProjectFilesAsync("PXD012345", _tempDir); // no filter => download everything

        Assert.That(paths.Count, Is.EqualTo(2));
        Assert.That(File.Exists(Path.Combine(_tempDir, "one.raw")), Is.True);
        Assert.That(File.Exists(Path.Combine(_tempDir, "two.mzid")), Is.True);
        Assert.That(handler.RequestedUris.Count(u => u.Contains("ftp.pride.ebi.ac.uk")), Is.EqualTo(2));
    }

    [Test]
    public async Task DownloadProjectFilesAsync_FilterMatchesNothing_ReturnsEmpty_NoDownloads()
    {
        string manifest = "[" + FileJson("one.raw", "RAW", "PRIDE:0000404") + "]";
        var handler = new StubHandler(request =>
        {
            string uri = request.RequestUri.ToString();
            if (uri.Contains("/files")) return Json(manifest);
            return Bytes(Encoding.UTF8.GetBytes("bytes-of-" + uri));
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var paths = await client.DownloadProjectFilesAsync("PXD012345", _tempDir, _ => false); // matches nothing

        Assert.That(paths, Is.Empty);
        Assert.That(handler.RequestedUris.Count(u => u.Contains("ftp.pride.ebi.ac.uk")), Is.EqualTo(0)); // manifest only
    }

    /// <summary>A single PRIDE file object in the real v3 wire shape (FTP + Aspera locations).</summary>
    private static string FileJson(string fileName, string category, string categoryAccession) =>
        $$"""
        {
          "fileCategory": { "@type": "CvParam", "cvLabel": "PRIDE", "accession": "{{categoryAccession}}", "name": "category", "value": "{{category}}" },
          "publicFileLocations": [
            { "@type": "CvParam", "cvLabel": "PRIDE", "accession": "PRIDE:0000469", "name": "FTP Protocol", "value": "ftp://ftp.pride.ebi.ac.uk/pride/data/x/{{fileName}}" },
            { "@type": "CvParam", "cvLabel": "PRIDE", "accession": "PRIDE:0000468", "name": "Aspera Protocol", "value": "prd_ascp@fasp.ebi.ac.uk:pride/data/x/{{fileName}}" }
          ],
          "fileSizeBytes": 10,
          "fileName": "{{fileName}}"
        }
        """;
}

/// <summary>
/// Live canary for the download path against the real PRIDE Archive + FTP-over-HTTPS host. Carries
/// [Category("ExternalService")] so CI runs it in the non-blocking external-service job, and
/// [Category("Pride")] to select it alone. Routed through <see cref="ExternalServiceTestHelper.RunAsync"/>
/// so a PRIDE/EBI outage SKIPS, while a genuine contract break (no HTTPS-reachable location, wrong bytes)
/// FAILS. Downloads the smallest file in a small public project to a temp dir, then deletes it.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideArchiveDownloadLiveTests
{
    [Test]
    public Task DownloadFileAsync_LiveSmallestFile_WritesRealBytes() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient();
            var files = await client.GetProjectFilesAsync("PXD012345");
            var smallest = files.OrderBy(f => f.FileSizeBytes).First(f => f.TryGetHttpsDownloadUrl(out _));

            string dir = Path.Combine(Path.GetTempPath(), "PrideLiveDownload", Guid.NewGuid().ToString("N"));
            try
            {
                string path = await client.DownloadFileAsync(smallest, dir);
                Assert.That(File.Exists(path), Is.True);
                Assert.That(new FileInfo(path).Length, Is.GreaterThan(0));
            }
            finally
            {
                if (Directory.Exists(dir)) Directory.Delete(dir, recursive: true);
            }
        });
}
