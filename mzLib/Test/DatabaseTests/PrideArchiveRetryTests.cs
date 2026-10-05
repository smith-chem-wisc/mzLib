using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Net.Http.Headers;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

/// <summary>
/// PR-B: bounded retry of transient failures, byte-range resume of a broken download, per-page retry of
/// the pager, and an empty FTP root listing treated as a transport failure. Offline: every response comes
/// from a scripted handler, and backoffs are recorded instead of slept through.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class PrideArchiveRetryTests
{
    // ---- test doubles -------------------------------------------------------

    /// <summary>
    /// Answers the Nth request with the Nth responder (the last one repeats) and records every request, with
    /// its Range and If-Range headers, so a test can assert both what was asked and in what order.
    /// </summary>
    private sealed class ScriptedHandler : HttpMessageHandler
    {
        private readonly Func<HttpRequestMessage, Task<HttpResponseMessage>>[] _script;
        public List<HttpRequestMessage> Requests { get; } = new();

        public ScriptedHandler(params Func<HttpRequestMessage, HttpResponseMessage>[] script) =>
            _script = script.Select(s => (Func<HttpRequestMessage, Task<HttpResponseMessage>>)(r => Task.FromResult(s(r)))).ToArray();

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            Requests.Add(request);
            return _script[Math.Min(Requests.Count - 1, _script.Length - 1)](request);
        }
    }

    /// <summary>Routes by URI so a multi-endpoint call (project, then FTP listing) can be scripted per endpoint.</summary>
    private sealed class RoutingHandler : HttpMessageHandler
    {
        private readonly Func<string, int, CancellationToken, Task<HttpResponseMessage>> _route;
        private readonly Dictionary<string, int> _counts = new();
        public List<string> RequestedUris { get; } = new();

        public RoutingHandler(Func<string, int, HttpResponseMessage> route) =>
            _route = (uri, n, _) => Task.FromResult(route(uri, n));

        public RoutingHandler(Func<string, int, CancellationToken, Task<HttpResponseMessage>> route) => _route = route;

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            string uri = request.RequestUri!.ToString();
            RequestedUris.Add(uri);
            _counts[uri] = _counts.TryGetValue(uri, out int n) ? n + 1 : 1;
            return _route(uri, _counts[uri], cancellationToken); // n is 1-based: the how-many-th request for this URI
        }
    }

    /// <summary>Yields the first <paramref name="count"/> bytes of <paramref name="data"/> and then throws, as a dropped connection does.</summary>
    private sealed class DroppingStream : Stream
    {
        private readonly byte[] _data;
        private readonly int _count;
        private int _position;

        public DroppingStream(byte[] data, int count)
        {
            _data = data;
            _count = count;
        }

        public override int Read(byte[] buffer, int offset, int count)
        {
            if (_position >= _count)
                throw new HttpIOException(HttpRequestError.ResponseEnded, "The response ended prematurely.");
            int n = Math.Min(count, _count - _position);
            Array.Copy(_data, _position, buffer, offset, n);
            _position += n;
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

    private static readonly byte[] Source = Enumerable.Range(0, 1000).Select(i => (byte)(i * 7 % 251)).ToArray();
    private const string ETag = "\"v1-abc\"";
    private const string FileUrl = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD000001/run1.raw";

    /// <summary>A whole-file 200 carrying <see cref="Source"/>'s announced length, a validator, and a body that drops after <paramref name="dropAfter"/> bytes (never, if null).</summary>
    private static HttpResponseMessage Full(int? dropAfter = null, string etag = ETag, DateTimeOffset? lastModified = null)
    {
        HttpContent content = dropAfter is int n
            ? new StreamContent(new DroppingStream(Source, n))
            : new ByteArrayContent(Source);
        content.Headers.ContentLength = Source.Length;
        if (lastModified.HasValue)
            content.Headers.LastModified = lastModified;
        var response = new HttpResponseMessage(HttpStatusCode.OK) { Content = content };
        if (etag != null)
            response.Headers.ETag = EntityTagHeaderValue.Parse(etag);
        return response;
    }

    /// <summary>A 206 serving <see cref="Source"/> from byte <paramref name="from"/> to the end, as Content-Range says.</summary>
    private static HttpResponseMessage Partial(long from, int? dropAfter = null)
    {
        byte[] rest = Source.Skip((int)from).ToArray();
        HttpContent content = dropAfter is int n ? new StreamContent(new DroppingStream(rest, n)) : new ByteArrayContent(rest);
        content.Headers.ContentLength = rest.Length;
        content.Headers.ContentRange = new ContentRangeHeaderValue(from, Source.Length - 1, Source.Length);
        var response = new HttpResponseMessage(HttpStatusCode.PartialContent) { Content = content };
        response.Headers.ETag = EntityTagHeaderValue.Parse(ETag);
        return response;
    }

    private static HttpResponseMessage Status(HttpStatusCode status, string body = "") =>
        new(status) { Content = new StringContent(body) };

    private static PrideArchiveFile MakeFile(string url = FileUrl) => new()
    {
        FileName = "run1.raw",
        FileSizeBytes = Source.Length,
        FileCategory = new CvParam("PRIDE", "PRIDE:0000404", "category", "RAW"),
        PublicFileLocations = { new CvParam("PRIDE", PrideArchiveExtensions.FtpLocationAccession, "FTP Protocol", url) },
    };

    private const string ProjectJson =
        """{ "accession": "PXD000001", "title": "t", "publicationDate": "2012-03-13" }""";

    private List<TimeSpan> _delays;
    private string _tempDir;

    /// <summary>A client over <paramref name="handler"/> whose backoffs are recorded, not slept through.</summary>
    private PrideArchiveClient ClientOver(HttpMessageHandler handler, int maxRetries = 3) =>
        new(new HttpClient(handler))
        {
            MaxRetries = maxRetries,
            RetryDelay = (delay, _) =>
            {
                _delays.Add(delay);
                return Task.CompletedTask;
            },
        };

    private string PartialPath => Path.Combine(_tempDir, "run1.raw.partial");
    private string ValidatorPath => PartialPath + PrideArchiveClient.ValidatorSuffix;
    private string DestinationPath => Path.Combine(_tempDir, "run1.raw");

    [SetUp]
    public void SetUp()
    {
        _delays = new List<TimeSpan>();
        _tempDir = Path.Combine(Path.GetTempPath(), "PrideArchiveRetryTests", Guid.NewGuid().ToString("N"));
    }

    [TearDown]
    public void TearDown()
    {
        try { if (Directory.Exists(_tempDir)) Directory.Delete(_tempDir, recursive: true); }
        catch { /* best-effort cleanup */ }
    }

    // ---- what is transient --------------------------------------------------

    [TestCase(HttpStatusCode.RequestTimeout)]
    [TestCase(HttpStatusCode.TooManyRequests)]
    [TestCase(HttpStatusCode.InternalServerError)]
    [TestCase(HttpStatusCode.BadGateway)]
    [TestCase(HttpStatusCode.ServiceUnavailable)]
    [TestCase(HttpStatusCode.GatewayTimeout)]
    public async Task RestRequest_TransientStatus_RecoversOnTheSecondAttempt(HttpStatusCode status)
    {
        var handler = new ScriptedHandler(_ => Status(status), _ => Status(HttpStatusCode.OK, ProjectJson));
        using var client = ClientOver(handler);

        PrideProject project = await client.GetProjectAsync("PXD000001");

        Assert.Multiple(() =>
        {
            Assert.That(project.Accession, Is.EqualTo("PXD000001"));
            Assert.That(handler.Requests, Has.Count.EqualTo(2));
            Assert.That(_delays, Is.EqualTo(new[] { TimeSpan.FromSeconds(5) }), "the first backoff is 5 s");
        });
    }

    [Test]
    public async Task Download_403FromTheFtpHost_IsRateLimitingAndIsRetried()
    {
        // EBI rate-limits its FTP host with 403, not 429 (field notes §1b), so there a 403 is transient.
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.Forbidden), _ => Full());
        using var client = ClientOver(handler);

        string path = await client.DownloadFileAsync(MakeFile(), _tempDir);

        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
            Assert.That(handler.Requests, Has.Count.EqualTo(2));
        });
    }

    [Test]
    public void RestRequest_403_IsARefusalAndIsNotRetried()
    {
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.Forbidden), _ => Status(HttpStatusCode.OK, ProjectJson));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetProjectAsync("PXD000001"));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.StatusCode, Is.EqualTo(HttpStatusCode.Forbidden));
            Assert.That(handler.Requests, Has.Count.EqualTo(1));
        });
    }

    [Test]
    public void Download_404_IsNotRetried()
    {
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.NotFound), _ => Full());
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
            Assert.That(handler.Requests, Has.Count.EqualTo(1));
            Assert.That(_delays, Is.Empty);
        });
    }

    [Test]
    public void BudgetExhausted_ThrowsTheLastStatusUnchanged_AfterTheFullBackoffSchedule()
    {
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.ServiceUnavailable));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetProjectAsync("PXD000001"));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.StatusCode, Is.EqualTo(HttpStatusCode.ServiceUnavailable));
            // pyMzLib's ClassifyError anchors a regex on this phrase (aging 013); it must survive retries.
            Assert.That(exception.Message, Does.Contain("failed with status 503"));
            Assert.That(handler.Requests, Has.Count.EqualTo(4), "one attempt and three retries");
            Assert.That(_delays, Is.EqualTo(new[] { TimeSpan.FromSeconds(5), TimeSpan.FromSeconds(20), TimeSpan.FromSeconds(60) }));
        });
    }

    [TestCase(HttpStatusCode.TooManyRequests, 7, 7)]
    [TestCase(HttpStatusCode.ServiceUnavailable, 600, 60)]
    [TestCase(HttpStatusCode.ServiceUnavailable, 0, 0)]
    public async Task RetryAfter_On429Or503_SetsTheWait_CappedAtSixtySeconds(HttpStatusCode status, int retryAfterSeconds, int expectedSeconds)
    {
        var handler = new ScriptedHandler(
            _ =>
            {
                HttpResponseMessage response = Status(status);
                response.Headers.RetryAfter = new RetryConditionHeaderValue(TimeSpan.FromSeconds(retryAfterSeconds));
                return response;
            },
            _ => Status(HttpStatusCode.OK, ProjectJson));
        using var client = ClientOver(handler);

        await client.GetProjectAsync("PXD000001");

        Assert.That(_delays, Is.EqualTo(new[] { TimeSpan.FromSeconds(expectedSeconds) }));
    }

    [Test]
    public async Task RetryAfter_OnAStatusOtherThan429Or503_IsIgnored()
    {
        var handler = new ScriptedHandler(
            _ =>
            {
                HttpResponseMessage response = Status(HttpStatusCode.BadGateway);
                response.Headers.RetryAfter = new RetryConditionHeaderValue(TimeSpan.FromSeconds(1));
                return response;
            },
            _ => Status(HttpStatusCode.OK, ProjectJson));
        using var client = ClientOver(handler);

        await client.GetProjectAsync("PXD000001");

        Assert.That(_delays, Is.EqualTo(new[] { TimeSpan.FromSeconds(5) }));
    }

    [Test]
    public void CancellationDuringABackoff_EndsItPromptly_AndSendsNothingMore()
    {
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.ServiceUnavailable));
        using var cts = new CancellationTokenSource(TimeSpan.FromMilliseconds(200));
        using var client = new PrideArchiveClient(new HttpClient(handler)); // the real 5 s backoff, really waited

        var started = DateTime.UtcNow;
        var exception = Assert.CatchAsync(async () => await client.GetProjectAsync("PXD000001", cts.Token));

        Assert.Multiple(() =>
        {
            Assert.That(exception, Is.InstanceOf<OperationCanceledException>());
            Assert.That(exception, Is.Not.InstanceOf<HttpRequestException>());
            Assert.That(DateTime.UtcNow - started, Is.LessThan(TimeSpan.FromSeconds(4)), "the 5 s backoff was cut short");
            Assert.That(handler.Requests, Has.Count.EqualTo(1));
        });
    }

    [Test]
    public void ABrokenContract_IsNeverRetried()
    {
        // 200 with a payload that is not a project: PRIDE answered, so a second ask will not mend it.
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.OK, "{}"));
        using var client = ClientOver(handler);

        Assert.ThrowsAsync<MzLibException>(async () => await client.GetProjectAsync("PXD000001"));
        Assert.That(handler.Requests, Has.Count.EqualTo(1));
    }

    [Test]
    public void AWriteSideIOException_IsNeverRetried()
    {
        // The partial is held open by someone else: that is the caller's disk, not an EBI outage.
        Directory.CreateDirectory(_tempDir);
        using var locker = new FileStream(PartialPath, FileMode.Create, FileAccess.Write, FileShare.None);
        var handler = new ScriptedHandler(_ => Full());
        using var client = ClientOver(handler);

        var exception = Assert.CatchAsync(async () => await client.DownloadFileAsync(MakeFile(), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(exception, Is.InstanceOf<IOException>());
            Assert.That(handler.Requests, Has.Count.EqualTo(1));
        });
    }

    // ---- "failed with status NNN" at every status-throw site (aging 013) ----

    private static IEnumerable<TestCaseData> StatusThrowSites()
    {
        yield return new TestCaseData((Func<PrideArchiveClient, Task>)(c => c.GetProjectAsync("PXD000001")))
            .SetArgDisplayNames("REST request");
        yield return new TestCaseData((Func<PrideArchiveClient, Task>)(c => c.GetProjectFilesAsync("PXD000001")))
            .SetArgDisplayNames("paged request");
        yield return new TestCaseData((Func<PrideArchiveClient, Task>)(c => c.GetProxiSpectrumAsync("mzspec:PXD000001:run:scan:1")))
            .SetArgDisplayNames("PROXI");
    }

    [TestCaseSource(nameof(StatusThrowSites))]
    public void StatusThrow_KeepsTheFailedWithStatusPhrase(Func<PrideArchiveClient, Task> call)
    {
        using var client = ClientOver(new ScriptedHandler(_ => Status(HttpStatusCode.BadRequest)));

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await call(client));

        Assert.That(exception!.Message, Does.Contain("failed with status 400"));
    }

    [Test]
    public void StatusThrow_Download_KeepsTheFailedWithStatusPhrase()
    {
        using var client = ClientOver(new ScriptedHandler(_ => Status(HttpStatusCode.Gone)));

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir));

        Assert.That(exception!.Message, Does.Contain("failed with status 410"));
    }

    [Test]
    public void StatusThrow_FtpDirectoryListing_KeepsTheFailedWithStatusPhrase()
    {
        var handler = new RoutingHandler((uri, _) =>
            uri.Contains("/projects/") ? Status(HttpStatusCode.OK, ProjectJson) : Status(HttpStatusCode.Gone));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetProjectFilesFromFtpAsync("PXD000001"));

        Assert.That(exception!.Message, Does.Contain("failed with status 410"));
    }

    // ---- resume ----------------------------------------------------------------

    [Test]
    public async Task ADroppedBody_IsResumedWithARangeUnderIfRange_AndTheFileIsByteIdentical()
    {
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400), _ => Partial(400));
        using var client = ClientOver(handler);

        string path = await client.DownloadFileAsync(MakeFile(), _tempDir);

        HttpRequestMessage resume = handler.Requests[1];
        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
            Assert.That(handler.Requests[0].Headers.Range, Is.Null, "the first request asks for the whole file");
            Assert.That(resume.Headers.Range?.ToString(), Is.EqualTo("bytes=400-"));
            Assert.That(resume.Headers.IfRange?.EntityTag?.Tag, Is.EqualTo(ETag));
            Assert.That(File.Exists(PartialPath), Is.False);
            Assert.That(File.Exists(ValidatorPath), Is.False, "the validator goes with the partial on success");
        });
    }

    [Test]
    public async Task NoETag_ResumesUnderLastModified()
    {
        var lastModified = new DateTimeOffset(2020, 1, 2, 3, 4, 5, TimeSpan.Zero);
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400, etag: null, lastModified: lastModified), _ => Partial(400));
        using var client = ClientOver(handler);

        await client.DownloadFileAsync(MakeFile(), _tempDir);

        Assert.That(handler.Requests[1].Headers.IfRange?.Date, Is.EqualTo(lastModified));
    }

    [Test]
    public async Task AWeakETag_IsNotUsed_LastModifiedIs()
    {
        // If-Range takes a strong validator only; a weak ETag there would be ignored or refused by the server.
        var lastModified = new DateTimeOffset(2020, 1, 2, 3, 4, 5, TimeSpan.Zero);
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400, etag: "W/\"weak\"", lastModified: lastModified), _ => Partial(400));
        using var client = ClientOver(handler);

        await client.DownloadFileAsync(MakeFile(), _tempDir);

        Assert.Multiple(() =>
        {
            Assert.That(handler.Requests[1].Headers.IfRange?.EntityTag, Is.Null);
            Assert.That(handler.Requests[1].Headers.IfRange?.Date, Is.EqualTo(lastModified));
        });
    }

    [Test]
    public async Task NoValidatorAtAll_RestartsFromZero_NeverSendsABareRange()
    {
        // A Range with no If-Range could splice two versions of a file. PRIDE's reviewer-token route sends no
        // validator (phred 006), so this is the path a private download takes.
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400, etag: null), _ => Full(etag: null));
        using var client = ClientOver(handler);

        string path = await client.DownloadFileAsync(MakeFile(), _tempDir);

        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
            Assert.That(handler.Requests[1].Headers.Range, Is.Null);
        });
    }

    [Test]
    public async Task A200ToARangedRequest_RestartsFromZero_AndTheFileIsByteIdentical()
    {
        // The file changed on the server, so If-Range made it answer with the whole new file.
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400), _ => Full());
        using var client = ClientOver(handler);

        string path = await client.DownloadFileAsync(MakeFile(), _tempDir);

        Assert.Multiple(() =>
        {
            Assert.That(handler.Requests[1].Headers.Range, Is.Not.Null);
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source), "the 200's body replaced the partial, not appended to it");
        });
    }

    [TestCase(300)]
    [TestCase(500)]
    public async Task A206AtTheWrongOffset_IsDiscarded_AndTheWholeFileIsFetched(int wrongStart)
    {
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400), _ => Partial(wrongStart), _ => Full());
        using var client = ClientOver(handler);

        string path = await client.DownloadFileAsync(MakeFile(), _tempDir);

        Assert.Multiple(() =>
        {
            Assert.That(handler.Requests, Has.Count.EqualTo(3));
            Assert.That(handler.Requests[2].Headers.Range, Is.Null);
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
            // One backoff, for the drop. The misplaced 206 is refused on sight, within the same attempt; it is
            // never spliced on and left for the length check to catch, which would cost a retry.
            Assert.That(_delays, Has.Count.EqualTo(1));
        });
    }

    [Test]
    public async Task A416_IsDiscarded_AndTheWholeFileIsFetched()
    {
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400), _ => Status(HttpStatusCode.RequestedRangeNotSatisfiable), _ => Full());
        using var client = ClientOver(handler);

        string path = await client.DownloadFileAsync(MakeFile(), _tempDir);

        Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
    }

    [Test]
    public void ALengthMismatch_Throws_AndLeavesNeitherPartialNorValidator()
    {
        // The body ends cleanly but short of the announced Content-Length, on every attempt.
        var handler = new ScriptedHandler(_ =>
        {
            HttpResponseMessage response = Full();
            response.Content = new ByteArrayContent(Source.Take(600).ToArray());
            response.Content.Headers.ContentLength = Source.Length;
            return response;
        });
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.Message, Does.Contain("came to 600 bytes where the server announced 1000"));
            Assert.That(handler.Requests.All(r => r.Headers.Range == null), Is.True, "mismatched bytes are never resumed");
            Assert.That(File.Exists(PartialPath), Is.False);
            Assert.That(File.Exists(ValidatorPath), Is.False);
            Assert.That(File.Exists(DestinationPath), Is.False);
        });
    }

    [Test]
    public void ATransportFailure_KeepsThePartialAndItsValidator_WhenRetriesRunOut()
    {
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400),
            r => Partial(r.Headers.Range!.Ranges.Single().From!.Value, dropAfter: 100));
        using var client = ClientOver(handler);

        Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(handler.Requests, Has.Count.EqualTo(4));
            Assert.That(new FileInfo(PartialPath).Length, Is.EqualTo(400 + 3 * 100), "every attempt's bytes were kept");
            Assert.That(File.ReadAllText(ValidatorPath), Is.EqualTo(ETag));
            Assert.That(File.Exists(DestinationPath), Is.False);
        });
    }

    [Test]
    public void ARefusalAfterAPartialTransfer_DeletesThePartialAndItsValidator()
    {
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400), _ => Status(HttpStatusCode.NotFound));
        using var client = ClientOver(handler);

        Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(File.Exists(PartialPath), Is.False);
            Assert.That(File.Exists(ValidatorPath), Is.False);
        });
    }

    [Test]
    public void CancellationAfterAPartialTransfer_DeletesThePartialAndItsValidator()
    {
        using var cts = new CancellationTokenSource();
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400));
        using var client = new PrideArchiveClient(new HttpClient(handler))
        {
            RetryDelay = (_, _) =>
            {
                cts.Cancel();
                return Task.FromCanceled(cts.Token);
            },
        };

        Assert.CatchAsync<OperationCanceledException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir, cancellationToken: cts.Token));

        Assert.Multiple(() =>
        {
            Assert.That(File.Exists(PartialPath), Is.False);
            Assert.That(File.Exists(ValidatorPath), Is.False);
        });
    }

    [Test]
    public async Task APartialLeftByAnEarlierCall_IsResumedByTheNextOne()
    {
        // What makes aging's overwrite:false re-run cheap after a crash: the second call pays only for the rest.
        var first = new ScriptedHandler(_ => Full(dropAfter: 400));
        using (var client = ClientOver(first, maxRetries: 0))
            Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir, overwrite: false));
        Assume.That(File.Exists(ValidatorPath), Is.True);

        var second = new ScriptedHandler(_ => Partial(400));
        using (var client = ClientOver(second))
        {
            string path = await client.DownloadFileAsync(MakeFile(), _tempDir, overwrite: false);

            Assert.Multiple(() =>
            {
                Assert.That(second.Requests, Has.Count.EqualTo(1));
                Assert.That(second.Requests[0].Headers.Range?.ToString(), Is.EqualTo("bytes=400-"));
                Assert.That(second.Requests[0].Headers.IfRange?.EntityTag?.Tag, Is.EqualTo(ETag));
                Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
                Assert.That(File.Exists(ValidatorPath), Is.False);
            });
        }
    }

    [Test]
    public async Task APartialWithNoValidator_IsNotResumed()
    {
        // Nothing proves orphan bytes are the start of the file the server holds now.
        Directory.CreateDirectory(_tempDir);
        File.WriteAllBytes(PartialPath, new byte[400]);
        var handler = new ScriptedHandler(_ => Full());
        using var client = ClientOver(handler);

        string path = await client.DownloadFileAsync(MakeFile(), _tempDir);

        Assert.Multiple(() =>
        {
            Assert.That(handler.Requests[0].Headers.Range, Is.Null);
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
        });
    }

    [Test]
    public void ATransportFailureWithNoValidator_LeavesNoPartial()
    {
        // A partial that can never be resumed is not worth keeping.
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400, etag: null));
        using var client = ClientOver(handler);

        Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir));

        Assert.That(File.Exists(PartialPath), Is.False);
    }

    [Test]
    public void RetryAndResumeTraces_NeverCarryTheUrl()
    {
        // phred PRIDE-Q3 / 006: a reviewer-token href is a complete credential on its own. Every attempt here
        // fails, through a fresh request and two resumes, and nothing thrown may carry the query string.
        const string secret = "S3CR3T-REVIEWER-TOKEN";
        var handler = new ScriptedHandler(_ => Full(dropAfter: 400), _ => Partial(400, dropAfter: 10),
            _ => Status(HttpStatusCode.ServiceUnavailable));
        using var client = ClientOver(handler);
        var file = MakeFile($"https://private.example.org/files/run1.raw?token={secret}");

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(file, _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(handler.Requests, Has.Count.EqualTo(4));
            Assert.That(handler.Requests[1].Headers.Range, Is.Not.Null, "a resume was attempted");
            Assert.That(exception!.ToString(), Does.Not.Contain(secret));
            Assert.That(exception.Message, Does.Contain("'run1.raw' from private.example.org"));
        });
    }

    // ---- the pager retries one page, not the fetch (aging REQ-PRIDE-1) -------

    private static string SearchPage(params string[] accessions) =>
        "[" + string.Join(",", accessions.Select(a => $$"""{ "accession": "{{a}}", "title": "{{a}}" }""")) + "]";

    private static HttpResponseMessage Page(string body, int total)
    {
        HttpResponseMessage response = Status(HttpStatusCode.OK, body);
        response.Headers.Add("total_records", total.ToString());
        return response;
    }

    [Test]
    public async Task ASearchPageThatTimesOutOnce_IsRetried_AndOnlyThatPageIsReRead()
    {
        // The real failure: about 6% of search requests stall past HttpClient.Timeout.
        var handler = new RoutingHandler(async (uri, n, token) =>
        {
            if (uri.Contains("page=1") && n == 1)
                await Task.Delay(Timeout.Infinite, token);
            return uri.Contains("page=0")
                ? Page(SearchPage("PXD000001", "PXD000002"), 3)
                : Page(SearchPage("PXD000003"), 3);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler) { Timeout = TimeSpan.FromMilliseconds(300) })
        {
            RetryDelay = (_, _) => Task.CompletedTask,
        };

        var hits = await client.SearchProjectsAsync("liver", pageSize: 2);

        Assert.Multiple(() =>
        {
            Assert.That(hits.Select(h => h.Accession), Is.EqualTo(new[] { "PXD000001", "PXD000002", "PXD000003" }));
            Assert.That(handler.RequestedUris.Count(u => u.Contains("page=0")), Is.EqualTo(1), "a page already held is never re-read");
            Assert.That(handler.RequestedUris.Count(u => u.Contains("page=1")), Is.EqualTo(2));
        });
    }

    [Test]
    public void ASearchPageThatFailsEveryTime_PropagatesAfterItsRetries()
    {
        var handler = new RoutingHandler((uri, _) => uri.Contains("page=0")
            ? Page(SearchPage("PXD000001", "PXD000002"), 3)
            : Status(HttpStatusCode.BadGateway));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.SearchProjectsAsync("liver", pageSize: 2));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.StatusCode, Is.EqualTo(HttpStatusCode.BadGateway));
            Assert.That(handler.RequestedUris.Count(u => u.Contains("page=1")), Is.EqualTo(4));
        });
    }

    // ---- an empty FTP root is a transport failure, an empty subdirectory is not (§4b) ----

    private static string Index(params string[] rows) =>
        "<html><body><table>\n" +
        """<tr><th><a href="?C=N;O=D">Name</a></th></tr>""" + "\n" +
        """<tr><td><a href="/pride/data/archive/2012/03/">Parent Directory</a></td><td align="right">  - </td></tr>""" + "\n" +
        string.Join("\n", rows) + "\n</table></body></html>";

    private static string Row(string href, string size) =>
        $"""<tr><td><a href="{href}">{href}</a></td><td align="right">2021-10-20 04:56  </td><td align="right">{size}</td></tr>""";

    private const string RootUrl = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/";

    [Test]
    public async Task AnEmptyRootListing_IsRetried_AndRecoversWhenTheListingReturns()
    {
        // PXD058248 (46 files) once listed nothing, and the [] that came back read as "no such project".
        var handler = new RoutingHandler((uri, n) =>
        {
            if (uri.Contains("/projects/")) return Status(HttpStatusCode.OK, ProjectJson);
            return Status(HttpStatusCode.OK, n == 1 ? Index() : Index(Row("run1.raw", "210M")));
        });
        using var client = ClientOver(handler);

        List<PrideFtpFile> files = await client.GetProjectFilesFromFtpAsync("PXD000001");

        Assert.Multiple(() =>
        {
            Assert.That(files.Select(f => f.RelativePath), Is.EqualTo(new[] { "run1.raw" }));
            Assert.That(handler.RequestedUris.Count(u => u == RootUrl), Is.EqualTo(2));
            Assert.That(handler.RequestedUris.Count(u => u.Contains("/projects/")), Is.EqualTo(1), "only the listing is retried");
        });
    }

    [Test]
    public void AnEmptyRootListingThatPersists_IsATransportFailure_NotAMissingProject()
    {
        var handler = new RoutingHandler((uri, _) =>
            uri.Contains("/projects/") ? Status(HttpStatusCode.OK, ProjectJson) : Status(HttpStatusCode.OK, Index()));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetProjectFilesFromFtpAsync("PXD000001"));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.StatusCode, Is.Null, "a transport failure carries no status");
            Assert.That(exception.Message, Does.Contain("PXD000001"));
            Assert.That(exception.Message, Does.Not.Contain("https://"), "the message names the project, not the URL");
            Assert.That(handler.RequestedUris.Count(u => u == RootUrl), Is.EqualTo(4));
        });
    }

    [Test]
    public async Task AnEmptySubdirectory_IsSilent_AndTheOtherFilesAreReturned()
    {
        // PXD001174's empty .wiff directory is real, and sdrf's phantom-file check depends on seeing it so.
        var handler = new RoutingHandler((uri, _) =>
        {
            if (uri.Contains("/projects/")) return Status(HttpStatusCode.OK, ProjectJson);
            if (uri == RootUrl) return Status(HttpStatusCode.OK, Index(Row("run1.raw", "210M"), Row("wiff/", "-")));
            return Status(HttpStatusCode.OK, Index());
        });
        using var client = ClientOver(handler);

        List<PrideFtpFile> files = await client.GetProjectFilesFromFtpAsync("PXD000001");

        Assert.Multiple(() =>
        {
            Assert.That(files.Select(f => f.RelativePath), Is.EqualTo(new[] { "run1.raw" }));
            Assert.That(handler.RequestedUris.Count(u => u == RootUrl + "wiff/"), Is.EqualTo(1), "an empty subdirectory is not retried");
        });
    }
}

/// <summary>
/// Live canary for byte-range resume against the real PRIDE FTP-over-HTTPS host. It skips on an EBI outage
/// (<see cref="ExternalServiceTestHelper.RunAsync"/>) and fails only if resume itself is broken.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideArchiveResumeLiveTests
{
    /// <summary>Records the Range header of every request it passes on.</summary>
    private sealed class RangeRecorder : DelegatingHandler
    {
        public List<RangeHeaderValue> Ranges { get; } = new();
        public RangeRecorder() : base(new HttpClientHandler()) { }

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            Ranges.Add(request.Headers.Range);
            return base.SendAsync(request, cancellationToken);
        }
    }

    [Test]
    public Task DownloadFileAsync_LiveHalfFile_IsResumedToTheSameBytes() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            string dir = Path.Combine(Path.GetTempPath(), "PrideLiveResume", Guid.NewGuid().ToString("N"));
            try
            {
                var recorder = new RangeRecorder();
                using var http = new HttpClient(recorder) { Timeout = TimeSpan.FromSeconds(100) };
                using var client = new PrideArchiveClient(http);

                var files = await client.GetProjectFilesAsync("PXD012345");
                var smallest = files.OrderBy(f => f.FileSizeBytes).First(f => f.TryGetHttpsDownloadUrl(out _) && f.FileSizeBytes > 1024);
                string path = await client.DownloadFileAsync(smallest, dir);
                byte[] whole = await File.ReadAllBytesAsync(path);

                // Cut the file at half and leave it as a .partial with the validator the host serves now.
                using var head = new HttpRequestMessage(HttpMethod.Head, smallest.GetHttpsDownloadUrl());
                using HttpResponseMessage headers = await http.SendAsync(head);
                headers.EnsureSuccessStatusCode();
                Assume.That(headers.Headers.ETag, Is.Not.Null, "the public host sent an ETag on 2026-09-23");
                string partial = path + ".partial";
                await File.WriteAllBytesAsync(partial, whole.Take(whole.Length / 2).ToArray());
                await File.WriteAllTextAsync(partial + PrideArchiveClient.ValidatorSuffix, headers.Headers.ETag!.ToString());
                File.Delete(path);
                recorder.Ranges.Clear();

                string resumed = await client.DownloadFileAsync(smallest, dir, overwrite: false);

                Assert.Multiple(() =>
                {
                    Assert.That(recorder.Ranges.Single()?.Ranges.Single().From, Is.EqualTo(whole.Length / 2));
                    Assert.That(File.ReadAllBytes(resumed), Is.EqualTo(whole));
                });
            }
            finally
            {
                if (Directory.Exists(dir)) Directory.Delete(dir, recursive: true);
            }
        });
}
