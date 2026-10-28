using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Security.Cryptography;
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

    /// <summary>A handler that never answers, honouring only cancellation: HttpClient.Timeout is what ends it.</summary>
    private sealed class SilentHandler : HttpMessageHandler
    {
        protected override async Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            await Task.Delay(Timeout.Infinite, cancellationToken);
            throw new InvalidOperationException("unreachable");
        }
    }

    /// <summary>A read-only stream that yields a few bytes and then throws, to simulate a mid-transfer failure.</summary>
    private sealed class ThrowingStream : Stream
    {
        private int _bytesBeforeThrow;
        private readonly Func<Exception> _failure;
        public ThrowingStream(int bytesBeforeThrow, Func<Exception> failure = null)
        {
            _bytesBeforeThrow = bytesBeforeThrow;
            _failure = failure ?? (() => new IOException("simulated mid-stream failure"));
        }

        public override int Read(byte[] buffer, int offset, int count)
        {
            if (_bytesBeforeThrow <= 0)
                throw _failure();
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

    /// <summary>
    /// A stream that delivers <paramref name="bytesBeforeStall"/> bytes and then goes permanently silent,
    /// honouring only cancellation. This is what a stalled HTTP body looks like: the connection is alive
    /// and the read never completes, which is precisely the case no timeout covered before BodyStallTimeout.
    /// </summary>
    private sealed class StallingStream : Stream
    {
        private int _bytesBeforeStall;
        public StallingStream(int bytesBeforeStall) => _bytesBeforeStall = bytesBeforeStall;

        public override async ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
        {
            if (_bytesBeforeStall > 0)
            {
                int n = Math.Min(buffer.Length, _bytesBeforeStall);
                buffer.Span.Slice(0, n).Fill(0x41);
                _bytesBeforeStall -= n;
                return n;
            }

            // Never completes on its own — only the stall window (or the caller) can end this. The await
            // must PROPAGATE the cancellation: a ContinueWith here would swallow it and return 0, which the
            // copy loop reads as a clean end-of-stream and the download then "succeeds" on a truncated file.
            await Task.Delay(Timeout.Infinite, cancellationToken).ConfigureAwait(false);
            return 0; // unreachable
        }

        public override int Read(byte[] buffer, int offset, int count) =>
            ReadAsync(buffer.AsMemory(offset, count)).AsTask().GetAwaiter().GetResult();

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

    /// <summary>
    /// A stream that is SLOW but never silent: it delivers one byte per read, pausing
    /// <paramref name="pausePerRead"/> between them, for <paramref name="totalBytes"/> bytes. Used to prove
    /// the deadline bounds inactivity and not duration — the total run exceeds BodyStallTimeout while no
    /// single gap does.
    /// </summary>
    private sealed class SlowStream : Stream
    {
        private readonly TimeSpan _pausePerRead;
        private int _remaining;
        public SlowStream(int totalBytes, TimeSpan pausePerRead)
        {
            _remaining = totalBytes;
            _pausePerRead = pausePerRead;
        }

        public override async ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
        {
            if (_remaining <= 0)
                return 0;
            await Task.Delay(_pausePerRead, cancellationToken).ConfigureAwait(false);
            buffer.Span[0] = 0x41;
            _remaining--;
            return 1;
        }

        public override int Read(byte[] buffer, int offset, int count) =>
            ReadAsync(buffer.AsMemory(offset, count)).AsTask().GetAwaiter().GetResult();

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

        Assert.That(file.TryGetHttpsDownloadUrl(out string url), Is.True);
        Assert.That(url, Is.EqualTo("https://example.org/direct/run1.raw"));
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
    public async Task DownloadFileAsync_NotOverwrite_ExistingFile_NoHttpsLocation_SkipsWithoutThrowing()
    {
        // resume must skip an already-present file before resolving its URL, so an Aspera-only entry
        // (no HTTPS location) that is already on disk returns untouched instead of throwing NotSupportedException.
        Directory.CreateDirectory(_tempDir);
        string existing = Path.Combine(_tempDir, "run1.raw");
        File.WriteAllText(existing, "old");
        var handler = new StubHandler(_ => Bytes(Encoding.UTF8.GetBytes("new")));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        var file = MakeFile("run1.raw", "RAW", Aspera("pride/data/x/run1.raw")); // no HTTPS-reachable location

        string path = await client.DownloadFileAsync(file, _tempDir, overwrite: false);

        Assert.That(handler.RequestedUris, Is.Empty);            // skipped before any URL resolution or network call
        Assert.That(File.ReadAllText(path), Is.EqualTo("old"));  // left untouched
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
            if (uri.Contains("/files/checksum/")) return Text(ChecksumList("keep.raw", "drop.mzid")); // checksum list
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
    [TestCase(".")]   // bare "." survives Path.GetFileName unchanged but is a non-leaf the contract rejects
    [TestCase("..")]  // bare ".." likewise resolves to the parent directory; must be rejected
    public void DownloadFileAsync_FileNameNotBareLeaf_Throws_WritesNothingOutside(string fileName)
    {
        // REGRESSION: failed before the fix — a path-bearing FileName escaped destinationDirectory via Path.Combine.
        var handler = new StubHandler(_ => Bytes(Encoding.UTF8.GetBytes("x")));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        var file = MakeFile(fileName, "RAW", Ftp("pride/data/x/leaf.raw"));

        Assert.That(async () => await client.DownloadFileAsync(file, _tempDir), Throws.ArgumentException);
        Assert.That(handler.RequestedUris, Is.Empty);                        // rejected before any network call
        // the location an un-prevented Path.Combine would have written to (parent for "../evil.raw",
        // a nested subdir for "sub/dir/file.raw") — derived from the input so each case checks its own escape target
        string wouldBeTarget = Path.GetFullPath(Path.Combine(_tempDir, fileName));
        Assert.That(File.Exists(wouldBeTarget), Is.False);
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
    public void DownloadFileAsync_MidStreamFailure_ThrowsHttpRequestException_LeavesNoPartial()
    {
        // A dropped connection is an outage, not a local fault, so it must reach the caller as the contract's
        // transport type. Before this it escaped as a bare IOException, which neither the documented contract
        // nor ExternalServiceTestHelper covers -- an EBI reset reddened a live test instead of skipping it.
        // The ".partial"-then-move design must still clean up the sibling.
        var handler = new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.OK)
        {
            Content = new StreamContent(new ThrowingStream(bytesBeforeThrow: 4))
        });
        using var client = new PrideArchiveClient(new HttpClient(handler)) { MaxRetries = 0 };
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(file, _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(exception.InnerException, Is.InstanceOf<IOException>());
            Assert.That(exception.Message, Does.Contain("ended before the download completed"));
            Assert.That(File.Exists(Path.Combine(_tempDir, "run1.raw.partial")), Is.False);
            Assert.That(File.Exists(Path.Combine(_tempDir, "run1.raw")), Is.False);
        });
    }

    [Test]
    public void DownloadFileAsync_ResponseEndedPrematurely_ThrowsHttpRequestException()
    {
        // Seen live on a PXD015239 download: on .NET 10 a connection EBI drops mid-body surfaces as
        // HttpIOException(ResponseEnded) -- an IOException, not an HttpRequestException.
        var handler = new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.OK)
        {
            Content = new StreamContent(new ThrowingStream(bytesBeforeThrow: 4,
                () => new HttpIOException(HttpRequestError.ResponseEnded, "The response ended prematurely.")))
        });
        using var client = new PrideArchiveClient(new HttpClient(handler)) { MaxRetries = 0 };
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(file, _tempDir));

        Assert.That(exception.InnerException, Is.InstanceOf<HttpIOException>());
    }

    // ---- BodyStallTimeout: the response body gets an inactivity deadline ---

    [Test]
    public void DownloadFileAsync_BodyStalls_ThrowsHttpRequestExceptionAndLeavesNoPartial()
    {
        // Before BodyStallTimeout existed this hung until the socket gave up: the read is issued after
        // HttpClient.Timeout has already been spent on the headers (ResponseHeadersRead), so nothing
        // bounded it. Verified RED against the pre-fix CopyToAsync, where this test never returns.
        var handler = new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.OK)
        {
            Content = new StreamContent(new StallingStream(bytesBeforeStall: 4))
        });
        using var client = new PrideArchiveClient(new HttpClient(handler))
        {
            BodyStallTimeout = TimeSpan.FromMilliseconds(250),
            MaxRetries = 0
        };
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        var exception = Assert.ThrowsAsync<HttpRequestException>(
            async () => await client.DownloadFileAsync(file, _tempDir));

        Assert.Multiple(() =>
        {
            // HttpRequestException specifically, not MzLibException: silence from EBI is an outage, and
            // ExternalServiceTestHelper must be able to skip a live test on it rather than fail it.
            Assert.That(exception.Message, Does.Contain("delivered nothing"));
            Assert.That(File.Exists(Path.Combine(_tempDir, "run1.raw.partial")), Is.False);
            Assert.That(File.Exists(Path.Combine(_tempDir, "run1.raw")), Is.False);
        });
    }

    [Test]
    public void DownloadFileAsync_CallerCancelsDuringBody_StaysOperationCanceled()
    {
        // The stall window is LINKED to the caller's token, so both cancel the same read. They must still
        // be told apart: converting a deliberate cancellation into HttpRequestException would make a live
        // test skip (ExternalServiceTestHelper) when the test itself did the cancelling.
        var handler = new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.OK)
        {
            Content = new StreamContent(new StallingStream(bytesBeforeStall: 4))
        });
        using var client = new PrideArchiveClient(new HttpClient(handler))
        {
            BodyStallTimeout = TimeSpan.FromSeconds(30) // long: the CALLER must be what ends this, not the stall
        };
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));
        using var cts = new CancellationTokenSource(TimeSpan.FromMilliseconds(250));

        var exception = Assert.CatchAsync(async () => await client.DownloadFileAsync(file, _tempDir, cancellationToken: cts.Token));

        Assert.Multiple(() =>
        {
            Assert.That(exception, Is.InstanceOf<OperationCanceledException>());
            // Asserted explicitly, not left implied by the type above: this is the discrimination the
            // catch filter exists for, and it is what keeps a cancelled test from reading as an outage.
            Assert.That(exception, Is.Not.InstanceOf<HttpRequestException>());
        });
    }

    [Test]
    public async Task DownloadFileAsync_SlowButProgressingBody_Succeeds()
    {
        // The load-bearing one: this is an INACTIVITY deadline, not a total-duration cap. The transfer below
        // runs ~600ms against a 150ms deadline and must succeed, because no single gap reaches it. A total
        // cap would abort exactly the large healthy downloads PRIDE is full of.
        var handler = new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.OK)
        {
            Content = new StreamContent(new SlowStream(totalBytes: 12, pausePerRead: TimeSpan.FromMilliseconds(50)))
        });
        using var client = new PrideArchiveClient(new HttpClient(handler))
        {
            BodyStallTimeout = TimeSpan.FromMilliseconds(150)
        };
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        string path = await client.DownloadFileAsync(file, _tempDir);

        Assert.That(new FileInfo(path).Length, Is.EqualTo(12));
    }

    [Test]
    public void BodyStallTimeout_DefaultsToTwoMinutes()
    {
        // Pins the shipped default. Two minutes is ProteinDbRetriever's measured value (#1189), not a fresh
        // guess, and this is the whole guard against a future edit quietly making it generous again.
        using var client = new PrideArchiveClient(new HttpClient(new StubHandler(_ => Bytes(new byte[] { 1 }))));

        Assert.That(client.BodyStallTimeout, Is.EqualTo(TimeSpan.FromMinutes(2)));
    }

    // ---- Error contract: status codes, header timeout, no URL in a message --

    [Test]
    public void DownloadFileAsync_NonSuccess_CarriesStatusCode()
    {
        // No PRIDE throw used to set StatusCode, so a caller could only tell a 403 from a 503 by parsing the
        // message. Retrying, and a caller deciding whether to, both need the status as data.
        var handler = new StubHandler(_ => Bytes(Encoding.UTF8.GetBytes("forbidden"), HttpStatusCode.Forbidden));
        using var client = new PrideArchiveClient(new HttpClient(handler)) { MaxRetries = 0 };
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(file, _tempDir));

        Assert.That(exception.StatusCode, Is.EqualTo(HttpStatusCode.Forbidden));
    }

    [Test]
    public void GetProjectFilesAsync_NonSuccess_CarriesStatusCode()
    {
        // The REST path shares the rule: every status throw in the client carries the status.
        var handler = new StubHandler(_ => Json("{}", HttpStatusCode.ServiceUnavailable));
        using var client = new PrideArchiveClient(new HttpClient(handler)) { MaxRetries = 0 };

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetProjectFilesAsync("PXD000001"));

        Assert.That(exception.StatusCode, Is.EqualTo(HttpStatusCode.ServiceUnavailable));
    }

    [Test]
    public void DownloadFileAsync_ServerNeverAnswers_ThrowsHttpRequestException()
    {
        // HttpClient.Timeout expiring throws TaskCanceledException, indistinguishable by type from a caller's
        // cancellation, so it used to escape the contract as an OperationCanceledException.
        using var client = new PrideArchiveClient(new HttpClient(new SilentHandler()) { Timeout = TimeSpan.FromMilliseconds(200) }) { MaxRetries = 0 };
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(file, _tempDir));

        Assert.That(exception.Message, Does.Contain("did not respond"));
    }

    [Test]
    public void GetProjectAsync_ServerNeverAnswers_ThrowsHttpRequestException()
    {
        using var client = new PrideArchiveClient(new HttpClient(new SilentHandler()) { Timeout = TimeSpan.FromMilliseconds(200) }) { MaxRetries = 0 };

        Assert.That(async () => await client.GetProjectAsync("PXD000001"), Throws.InstanceOf<HttpRequestException>());
    }

    [Test]
    public void DownloadFileAsync_CallerCancelsWhileWaitingForHeaders_StaysOperationCanceled()
    {
        // The other half of the header-timeout filter: a caller's own cancellation is never converted.
        using var client = new PrideArchiveClient(new HttpClient(new SilentHandler()) { Timeout = TimeSpan.FromSeconds(30) });
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));
        using var cts = new CancellationTokenSource(TimeSpan.FromMilliseconds(200));

        var exception = Assert.CatchAsync(async () => await client.DownloadFileAsync(file, _tempDir, cancellationToken: cts.Token));

        Assert.Multiple(() =>
        {
            Assert.That(exception, Is.InstanceOf<OperationCanceledException>());
            Assert.That(exception, Is.Not.InstanceOf<HttpRequestException>());
        });
    }

    [Test]
    public void DownloadFileAsync_InjectedClientCancelsPendingRequests_StaysOperationCanceled()
    {
        // A cancellation the caller asks for through the HttpClient they injected, not through the token, is
        // still not a timeout: only HttpClient.Timeout expiring is reported as PRIDE not responding.
        var httpClient = new HttpClient(new SilentHandler()) { Timeout = TimeSpan.FromSeconds(30) };
        using var client = new PrideArchiveClient(httpClient);
        var file = MakeFile("run1.raw", "RAW", Ftp("pride/data/x/run1.raw"));

        var download = client.DownloadFileAsync(file, _tempDir);
        Thread.Sleep(200);
        httpClient.CancelPendingRequests();
        var exception = Assert.CatchAsync(async () => await download);

        Assert.Multiple(() =>
        {
            Assert.That(exception, Is.InstanceOf<OperationCanceledException>());
            Assert.That(exception, Is.Not.InstanceOf<HttpRequestException>());
        });
    }

    /// <summary>Every way a download can fail, each run against a URL whose query string is a credential.</summary>
    private static IEnumerable<TestCaseData> DownloadFailureModes()
    {
        yield return new TestCaseData((Func<HttpMessageHandler>)(() =>
            new StubHandler(_ => Bytes(Encoding.UTF8.GetBytes("no"), HttpStatusCode.Forbidden))), 30_000)
            .SetArgDisplayNames("status");
        yield return new TestCaseData((Func<HttpMessageHandler>)(() =>
            new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.OK) { Content = new StreamContent(new ThrowingStream(4)) })), 30_000)
            .SetArgDisplayNames("dropped");
        yield return new TestCaseData((Func<HttpMessageHandler>)(() =>
            new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.OK) { Content = new StreamContent(new StallingStream(4)) })), 30_000)
            .SetArgDisplayNames("stall");
        yield return new TestCaseData((Func<HttpMessageHandler>)(() => new SilentHandler()), 200)
            .SetArgDisplayNames("timeout");
    }

    [TestCaseSource(nameof(DownloadFailureModes))]
    public void DownloadFileAsync_Failure_NeverPutsUrlInMessage(Func<HttpMessageHandler> handler, int timeoutMs)
    {
        // PRIDE's reviewer-token route hands out download hrefs whose query string IS the token. An exception
        // message reaches logs, CI output and pasted issues, so no failure path
        // may carry the URL. Checked over ToString(), which includes every inner exception.
        const string secret = "S3CR3T-REVIEWER-TOKEN";
        using var client = new PrideArchiveClient(new HttpClient(handler()) { Timeout = TimeSpan.FromMilliseconds(timeoutMs) })
        {
            BodyStallTimeout = TimeSpan.FromMilliseconds(200),
            RetryDelay = (_, _) => Task.CompletedTask // retried for real, so the retry path is held to the rule too
        };
        var file = MakeFile("run1.raw", "RAW", ("PRIDE:0000000", $"https://private.example.org/files/run1.raw?token={secret}"));

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(file, _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(exception.ToString(), Does.Not.Contain(secret));
            // Still useful: the message names what failed and where.
            Assert.That(exception.Message, Does.Contain("'run1.raw' from private.example.org"));
        });
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
            if (uri.Contains("/files/checksum/")) return Text(ChecksumList("one.raw", "two.mzid"));
            if (uri.Contains("/files")) return Json(manifest);
            return Bytes(Encoding.UTF8.GetBytes("bytes-of-" + uri));
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var paths = await client.DownloadProjectFilesAsync("PXD012345", _tempDir); // no filter => download everything

        Assert.That(paths.Count, Is.EqualTo(2));
        Assert.That(File.Exists(Path.Combine(_tempDir, "one.raw")), Is.True);
        Assert.That(File.Exists(Path.Combine(_tempDir, "two.mzid")), Is.True);
        Assert.That(handler.RequestedUris.Count(u => u.Contains("ftp.pride.ebi.ac.uk")), Is.EqualTo(2));
        Assert.That(handler.RequestedUris.Count(u => u.Contains("/files/checksum/")), Is.EqualTo(1)); // once, not per file
    }

    [Test]
    public async Task DownloadProjectFilesAsync_FilterMatchesNothing_ReturnsEmpty_NoDownloads()
    {
        string manifest = "[" + FileJson("one.raw", "RAW", "PRIDE:0000404") + "]";
        var handler = new StubHandler(request =>
        {
            string uri = request.RequestUri.ToString();
            if (uri.Contains("/files/checksum/")) return Text(ChecksumList("one.raw"));
            if (uri.Contains("/files")) return Json(manifest);
            return Bytes(Encoding.UTF8.GetBytes("bytes-of-" + uri));
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var paths = await client.DownloadProjectFilesAsync("PXD012345", _tempDir, _ => false); // matches nothing

        Assert.That(paths, Is.Empty);
        Assert.That(handler.RequestedUris.Count(u => u.Contains("ftp.pride.ebi.ac.uk")), Is.EqualTo(0)); // manifest only
        Assert.That(handler.RequestedUris.Count(u => u.Contains("/files/checksum/")), Is.EqualTo(0)); // nothing to check
    }

    // ---- DownloadProjectFilesAsync: each file is checked against PRIDE's checksum list --------

    /// <summary>The bytes the project stubs below serve for a file: what <see cref="ChecksumList"/> records.</summary>
    private static byte[] Served(string fileName) =>
        Encoding.UTF8.GetBytes("bytes-of-https://ftp.pride.ebi.ac.uk/pride/data/x/" + fileName);

    /// <summary>PRIDE's checksum list, with the true MD5 and size of what is served for each named file.</summary>
    private static string ChecksumList(params string[] fileNames) =>
        PrideArchiveClient.ChecksumHeader + "\n" + string.Concat(fileNames.Select(name =>
            $"{name}\t{Convert.ToHexString(MD5.HashData(Served(name))).ToLowerInvariant()}\t{Served(name).Length}\n"));

    private static HttpResponseMessage Text(string body) =>
        new(HttpStatusCode.OK) { Content = new StringContent(body, Encoding.UTF8, "text/plain") };

    /// <summary>
    /// A project whose manifest lists <paramref name="fileNames"/>, whose checksum list is
    /// <paramref name="checksums"/>, and whose files are served as <paramref name="body"/> gives them
    /// (default <see cref="Served"/>). One attempt per download: a refusal must not wait out the backoffs.
    /// </summary>
    private static (PrideArchiveClient Client, StubHandler Handler) Project(string[] fileNames, string checksums,
        Func<string, byte[]> body = null)
    {
        string manifest = "[" + string.Join(",", fileNames.Select(n => FileJson(n, "RAW", "PRIDE:0000404"))) + "]";
        var handler = new StubHandler(request =>
        {
            string uri = request.RequestUri.ToString();
            if (uri.Contains("/files/checksum/")) return Text(checksums);
            if (uri.Contains("/files")) return Json(manifest);
            string name = uri[(uri.LastIndexOf('/') + 1)..];
            return Bytes((body ?? Served)(name));
        });
        return (new PrideArchiveClient(new HttpClient(handler)) { MaxRetries = 0 }, handler);
    }

    private static int Downloads(StubHandler handler) => handler.RequestedUris.Count(u => u.Contains("ftp.pride.ebi.ac.uk"));

    [TestCase("short", TestName = "DownloadProjectFilesAsync_AFileShorterThanItsChecksumRow_IsRefused")]
    [TestCase("ebi", TestName = "DownloadProjectFilesAsync_EbisErrorTextAtTheRecordedSize_IsRefused")]
    public void DownloadProjectFilesAsync_ADamagedFile_IsRefused_AndNotKept(string damage)
    {
        byte[] good = Served("one.raw");
        byte[] bad = damage == "short"
            ? good[..^1]
            : Encoding.UTF8.GetBytes("<Error><Code>ConnectionClosedException</Code>".PadRight(good.Length, ' '));
        Assert.That(damage == "short" || bad.Length == good.Length); // EBI's damage keeps the announced size
        var (client, handler) = Project(new[] { "one.raw" }, ChecksumList("one.raw"), _ => bad);
        using (client)
        {
            var e = Assert.ThrowsAsync<MzLibException>(async () => await client.DownloadProjectFilesAsync("PXD012345", _tempDir));
            Assert.That(e!.Message, Does.Contain("one.raw"));
        }
        Assert.That(File.Exists(Path.Combine(_tempDir, "one.raw")), Is.False);
        Assert.That(Downloads(handler), Is.EqualTo(1));
    }

    [Test]
    public async Task DownloadProjectFilesAsync_AFileWithNoChecksumRow_IsDownloadedUnchecked()
    {
        // generated/ files have no row: they cannot be checked, which is not a reason to refuse them.
        var (client, handler) = Project(new[] { "listed.raw", "generated.txt" }, ChecksumList("listed.raw"),
            name => name == "generated.txt" ? Encoding.UTF8.GetBytes("any length at all") : Served(name));
        using (client)
        {
            var paths = await client.DownloadProjectFilesAsync("PXD012345", _tempDir);
            Assert.That(paths.Select(Path.GetFileName), Is.EqualTo(new[] { "listed.raw", "generated.txt" }));
        }
        Assert.That(File.ReadAllText(Path.Combine(_tempDir, "generated.txt")), Is.EqualTo("any length at all"));
        Assert.That(Downloads(handler), Is.EqualTo(2));
    }

    [TestCase(false)]
    [TestCase(true)]
    public async Task DownloadProjectFilesAsync_AWrongMd5_IsRefusedOnlyWhenTheMd5IsAskedFor(bool verifyMd5)
    {
        // Right size, wrong bytes: only the MD5 can see it.
        byte[] good = Served("one.raw");
        byte[] wrong = (byte[])good.Clone();
        wrong[0] ^= 0x01;
        var (client, _) = Project(new[] { "one.raw" }, ChecksumList("one.raw"), _ => wrong);
        using (client)
        {
            if (verifyMd5)
            {
                var e = Assert.ThrowsAsync<MzLibException>(async () =>
                    await client.DownloadProjectFilesAsync("PXD012345", _tempDir, null, overwrite: true, verifyMd5: true));
                Assert.That(e!.Message, Does.Contain("MD5"));
                Assert.That(File.Exists(Path.Combine(_tempDir, "one.raw")), Is.False);
            }
            else
            {
                var paths = await client.DownloadProjectFilesAsync("PXD012345", _tempDir, null, overwrite: true, verifyMd5: false);
                Assert.That(File.ReadAllBytes(paths.Single()), Is.EqualTo(wrong));
            }
        }
    }

    // ---- DownloadProjectFilesByNameAsync -----------------------------------------------------

    [Test]
    public async Task DownloadProjectFilesByNameAsync_DownloadsTheNamedFiles_InManifestOrder_EachOnce()
    {
        var (client, handler) = Project(new[] { "a.raw", "b.raw", "c.raw" }, ChecksumList("a.raw", "b.raw", "c.raw"));
        using (client)
        {
            var paths = await client.DownloadProjectFilesByNameAsync("PXD012345", _tempDir, new[] { "c.raw", "a.raw", "c.raw" });
            Assert.That(paths.Select(Path.GetFileName), Is.EqualTo(new[] { "a.raw", "c.raw" }));
        }
        Assert.That(File.Exists(Path.Combine(_tempDir, "b.raw")), Is.False);
        Assert.That(Downloads(handler), Is.EqualTo(2));
        Assert.That(handler.RequestedUris.Count(u => u.Contains("/files/checksum/")), Is.EqualTo(1));
    }

    [Test]
    public void DownloadProjectFilesByNameAsync_NamesNotInTheProject_AreRefusedBeforeAnyDownload()
    {
        // "A.RAW" differs from "a.raw" only in case: PRIDE's names are matched ordinally, as everywhere in this client.
        var (client, handler) = Project(new[] { "a.raw", "b.raw" }, ChecksumList("a.raw", "b.raw"));
        using (client)
        {
            var e = Assert.ThrowsAsync<MzLibException>(async () =>
                await client.DownloadProjectFilesByNameAsync("PXD012345", _tempDir, new[] { "a.raw", "typo.raw", "A.RAW" }));
            Assert.That(e!.Message, Does.Contain("2 of the 3 file names"));
            Assert.That(e.Message, Does.Contain("'typo.raw'").And.Contain("'A.RAW'"));
            Assert.That(e.Message, Does.Not.Contain("'a.raw'"));
        }
        Assert.That(Downloads(handler), Is.EqualTo(0));
        Assert.That(handler.RequestedUris.Count(u => u.Contains("/files/checksum/")), Is.EqualTo(0));
        Assert.That(Directory.Exists(_tempDir) ? Directory.GetFiles(_tempDir) : Array.Empty<string>(), Is.Empty);
    }

    [TestCase(null, TestName = "DownloadProjectFilesByNameAsync_NullNames_ThrowsBeforeAnyRequest")]
    [TestCase("", TestName = "DownloadProjectFilesByNameAsync_NoNames_ThrowsBeforeAnyRequest")]
    [TestCase("a.raw| ", TestName = "DownloadProjectFilesByNameAsync_ABlankName_ThrowsBeforeAnyRequest")]
    public void DownloadProjectFilesByNameAsync_BadNames_ThrowArgumentException_BeforeAnyRequest(string joinedNames)
    {
        string[] names = joinedNames == null ? null
            : joinedNames.Length == 0 ? Array.Empty<string>()
            : joinedNames.Split('|');
        var (client, handler) = Project(new[] { "a.raw" }, ChecksumList("a.raw"));
        using (client)
        {
            var e = Assert.ThrowsAsync<ArgumentException>(async () =>
                await client.DownloadProjectFilesByNameAsync("PXD012345", _tempDir, names));
            Assert.That(e!.ParamName, Is.EqualTo("fileNames"));
        }
        Assert.That(handler.RequestedUris, Is.Empty);
    }

    [TestCase(" ", "dir", "accession")]
    [TestCase("PXD012345", " ", "destinationDirectory")]
    public void DownloadProjectFilesByNameAsync_BlankAccessionOrDirectory_ThrowsBeforeAnyRequest(string accession,
        string directory, string paramName)
    {
        var (client, handler) = Project(new[] { "a.raw" }, ChecksumList("a.raw"));
        using (client)
        {
            var e = Assert.ThrowsAsync<ArgumentException>(async () =>
                await client.DownloadProjectFilesByNameAsync(accession, directory == "dir" ? _tempDir : directory, new[] { "a.raw" }));
            Assert.That(e!.ParamName, Is.EqualTo(paramName));
        }
        Assert.That(handler.RequestedUris, Is.Empty);
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
            // one attempt: a PRIDE outage skips at once, not after the 5/20/60 s retry backoffs
            using var client = new PrideArchiveClient { MaxRetries = 0 };
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
