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
using System.Threading;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

/// <summary>
/// M0: <see cref="ServiceFailure.IsTransient"/> (bridge PRIDE-B3) says "ask again later" only for a failure the
/// retry loop judged transient and gave up on, and the shared <see cref="DownloadEngine"/> carries no PRIDE rule
/// of its own. Offline: every response comes from a scripted handler, and backoffs are not slept through.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ServiceFailureTests
{
    /// <summary>Answers the Nth request with the Nth responder (the last one repeats) and counts the requests.</summary>
    private sealed class ScriptedHandler : HttpMessageHandler
    {
        private readonly Func<HttpRequestMessage, CancellationToken, Task<HttpResponseMessage>>[] _script;
        public int Requests { get; private set; }

        public ScriptedHandler(params Func<HttpRequestMessage, HttpResponseMessage>[] script) =>
            _script = script.Select(s => (Func<HttpRequestMessage, CancellationToken, Task<HttpResponseMessage>>)((r, _) => Task.FromResult(s(r)))).ToArray();

        public ScriptedHandler(Func<HttpRequestMessage, CancellationToken, Task<HttpResponseMessage>> respond) => _script = new[] { respond };

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            Requests++;
            return _script[Math.Min(Requests - 1, _script.Length - 1)](request, cancellationToken);
        }
    }

    /// <summary>Yields the first <paramref name="count"/> bytes of <paramref name="data"/> and then throws, as a dropped connection does.</summary>
    private sealed class DroppingStream(byte[] data, int count) : MemoryStream(data, 0, count)
    {
        public override int Read(byte[] buffer, int offset, int length)
        {
            int read = base.Read(buffer, offset, length);
            return read > 0 ? read : throw new HttpIOException(HttpRequestError.ResponseEnded, "The response ended prematurely.");
        }

        public override async ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
        {
            int read = await base.ReadAsync(buffer, cancellationToken);
            return read > 0 ? read : throw new HttpIOException(HttpRequestError.ResponseEnded, "The response ended prematurely.");
        }
    }

    private static readonly byte[] Source = Enumerable.Range(0, 1000).Select(i => (byte)(i * 7 % 251)).ToArray();
    private const string FtpFileUrl = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD000001/run1.raw";
    private const string ProjectJson = """{ "accession": "PXD000001", "title": "t", "publicationDate": "2012-03-13" }""";

    private int _backoffs;
    private string _tempDir;

    private static HttpResponseMessage Status(HttpStatusCode status, string body = "") =>
        new(status) { Content = new StringContent(body) };

    private static HttpResponseMessage Full(byte[] body)
    {
        var content = new ByteArrayContent(body);
        content.Headers.ContentLength = body.Length;
        return new HttpResponseMessage(HttpStatusCode.OK) { Content = content };
    }

    private static HttpResponseMessage Dropping(int after)
    {
        var content = new StreamContent(new DroppingStream(Source, after));
        content.Headers.ContentLength = Source.Length;
        return new HttpResponseMessage(HttpStatusCode.OK) { Content = content };
    }

    private static PrideArchiveFile MakeFile() => new()
    {
        FileName = "run1.raw",
        FileSizeBytes = Source.Length,
        FileCategory = new CvParam("PRIDE", "PRIDE:0000404", "category", "RAW"),
        PublicFileLocations = { new CvParam("PRIDE", PrideArchiveExtensions.FtpLocationAccession, "FTP Protocol", FtpFileUrl) },
    };

    /// <summary>A client over <paramref name="handler"/> whose backoffs are counted, not slept through.</summary>
    private PrideArchiveClient ClientOver(HttpMessageHandler handler, int maxRetries = 3) =>
        new(new HttpClient(handler))
        {
            MaxRetries = maxRetries,
            RetryDelay = (_, _) =>
            {
                _backoffs++;
                return Task.CompletedTask;
            },
        };

    [SetUp]
    public void SetUp()
    {
        _backoffs = 0;
        _tempDir = Path.Combine(Path.GetTempPath(), "ServiceFailureTests", Guid.NewGuid().ToString("N"));
    }

    [TearDown]
    public void TearDown()
    {
        try { if (Directory.Exists(_tempDir)) Directory.Delete(_tempDir, recursive: true); }
        catch { /* best-effort cleanup */ }
    }

    // ---- true: a transient failure the retry loop gave up on -------------------------------

    [TestCase(HttpStatusCode.RequestTimeout)]
    [TestCase(HttpStatusCode.TooManyRequests)]
    [TestCase(HttpStatusCode.InternalServerError)]
    [TestCase(HttpStatusCode.ServiceUnavailable)]
    public void RestRequest_TransientStatusThatPersists_IsTransient(HttpStatusCode status)
    {
        var handler = new ScriptedHandler(_ => Status(status));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetProjectAsync("PXD000001"));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(exception), Is.True);
            Assert.That(exception!.StatusCode, Is.EqualTo(status), "the mark leaves the exception itself unchanged");
            Assert.That(handler.Requests, Is.EqualTo(4), "it was retried to the end first");
        });
    }

    [Test]
    public void Download_403FromTheFtpHostThatPersists_IsTransient()
    {
        // The case the mark exists for: the same 403 is a refusal on the REST host (next test), and nothing on the
        // exception says which host it came from.
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.Forbidden));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.StatusCode, Is.EqualTo(HttpStatusCode.Forbidden));
            Assert.That(ServiceFailure.IsTransient(exception), Is.True);
            Assert.That(handler.Requests, Is.EqualTo(4));
        });
    }

    [Test]
    public void Download_BodyThatKeepsDropping_IsTransient()
    {
        var handler = new ScriptedHandler(_ => Dropping(after: 400));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.StatusCode, Is.Null, "a transport failure carries no status");
            Assert.That(ServiceFailure.IsTransient(exception), Is.True);
            Assert.That(handler.Requests, Is.EqualTo(4));
        });
    }

    [Test]
    public void WithNoRetries_ATransientFailure_IsTransientAfterOneAttempt()
    {
        // MaxRetries = 0 is what every live test sets: the one attempt is also the last, so the loop gives up on it.
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.BadGateway));
        using var client = ClientOver(handler, maxRetries: 0);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetProjectAsync("PXD000001"));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(exception), Is.True);
            Assert.That(handler.Requests, Is.EqualTo(1));
            Assert.That(_backoffs, Is.Zero);
        });
    }

    // ---- false: everything else -------------------------------------------------------------

    [Test]
    public void RestRequest_403_IsNotTransient()
    {
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.Forbidden));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetProjectAsync("PXD000001"));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(exception), Is.False);
            Assert.That(handler.Requests, Is.EqualTo(1));
        });
    }

    [Test]
    public void Download_404_IsNotTransient()
    {
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.NotFound));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir));

        Assert.That(ServiceFailure.IsTransient(exception), Is.False);
    }

    [Test]
    public void ProxiResponseWithNoSpectrum_IsNotTransient_ThoughItIsAnHttpRequestExceptionWithNoStatus()
    {
        // Thrown after the retry loop, never retried: its shape is a transport failure's, so a verdict read from the
        // shape would call it transient. The mark is what tells them apart.
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.OK, "[]"));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () =>
            await client.GetProxiSpectrumAsync("mzspec:PXD000001:run1:scan:1"));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.StatusCode, Is.Null);
            Assert.That(ServiceFailure.IsTransient(exception), Is.False);
            Assert.That(handler.Requests, Is.EqualTo(1));
        });
    }

    [Test]
    public void DownloadThatFailsItsChecksOnEveryAttempt_IsNotTransient_ThoughItsLastAttemptWas()
    {
        // The user's W1 rule: a file PRIDE serves wrong every time is a data error, not an outage. The MzLibException
        // wraps the last attempt, which the loop gave up on as transient; only the outer exception is read.
        var handler = new ScriptedHandler(_ => Full(Source));
        using var client = ClientOver(handler);
        var expected = new PrideFileChecksum("run1.raw", Convert.ToHexString(MD5.HashData(Source)).ToLowerInvariant(), Source.Length + 1);

        var exception = Assert.ThrowsAsync<MzLibException>(async () => await client.DownloadFileAsync(MakeFile(), _tempDir, expected));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(exception), Is.False);
            Assert.That(ServiceFailure.IsTransient(exception!.InnerException), Is.True, "the attempt inside it was transient");
            Assert.That(handler.Requests, Is.EqualTo(4));
        });
    }

    [Test]
    public void CallerCancellation_IsNotTransient()
    {
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.OK, ProjectJson));
        using var client = ClientOver(handler);
        using var cancelled = new CancellationTokenSource();
        cancelled.Cancel();

        var exception = Assert.CatchAsync<OperationCanceledException>(async () => await client.GetProjectAsync("PXD000001", cancelled.Token));

        Assert.That(ServiceFailure.IsTransient(exception), Is.False);
    }

    [Test]
    public void AnHttpRequestExceptionNoClientThrew_AndNull_AreNotTransient()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(new HttpRequestException("not from a mzLib client")), Is.False);
            Assert.That(ServiceFailure.IsTransient(new HttpRequestException("503", null, HttpStatusCode.ServiceUnavailable)), Is.False);
            Assert.That(ServiceFailure.IsTransient(null), Is.False);
        });
    }

    // ---- the engine itself is repository-neutral ------------------------------------------------

    [Test]
    public void Engine_WithNoRuleOfItsOwn_Treats403FromPridesFtpHostAsARefusal()
    {
        // PRIDE's 403 rule is PRIDE's, passed in: an engine built for another service does not inherit it.
        var handler = new ScriptedHandler(_ => Status(HttpStatusCode.Forbidden));
        var engine = new DownloadEngine(new HttpClient(handler), "MassIVE") { RetryDelay = (_, _) => Task.CompletedTask };

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () =>
            await engine.GetSuccessAsync("https://ftp.pride.ebi.ac.uk/pride/data/archive/", r => $"failed with status {(int)r.StatusCode}", CancellationToken.None));

        Assert.Multiple(() =>
        {
            Assert.That(handler.Requests, Is.EqualTo(1), "not retried");
            Assert.That(ServiceFailure.IsTransient(exception), Is.False);
        });
    }

    [Test]
    public void Engine_NamesItsOwnServiceInATimeout()
    {
        var handler = new ScriptedHandler(async (_, token) =>
        {
            await Task.Delay(Timeout.Infinite, token);
            return Status(HttpStatusCode.OK);
        });
        var engine = new DownloadEngine(new HttpClient(handler) { Timeout = TimeSpan.FromMilliseconds(100) }, "MassIVE") { MaxRetries = 0 };

        var exception = Assert.ThrowsAsync<HttpRequestException>(async () =>
            await engine.GetSuccessAsync("https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/MSV000000001", r => "unused", CancellationToken.None));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.Message, Does.StartWith("MassIVE did not respond within"));
            Assert.That(exception.InnerException, Is.InstanceOf<TaskCanceledException>());
            Assert.That(ServiceFailure.IsTransient(exception), Is.True);
        });
    }
}
