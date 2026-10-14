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
using System.Security.Cryptography;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

[TestFixture]
[ExcludeFromCodeCoverage]
public class PrideChecksumTests
{
    // ---- test doubles -------------------------------------------------------

    /// <summary>Answers each request from a caller-supplied function and records what was asked.</summary>
    private sealed class StubHandler : HttpMessageHandler
    {
        private readonly Func<HttpRequestMessage, HttpResponseMessage> _responder;
        public List<HttpRequestMessage> Requests { get; } = new();
        public List<string> RequestedUris => Requests.Select(r => r.RequestUri!.ToString()).ToList();

        public StubHandler(Func<HttpRequestMessage, HttpResponseMessage> responder) => _responder = responder;

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            Requests.Add(request);
            return Task.FromResult(_responder(request));
        }
    }

    // ---- fixtures -------------------------------------------------------------

    /// <summary>
    /// PRIDE's checksum list for PXD010820, captured verbatim from
    /// <c>https://www.ebi.ac.uk/pride/ws/archive/v3/files/checksum/PXD010820</c> on 2026-10-06 (LF line ends,
    /// a trailing newline). The SDRF row's MD5 was checked against the downloaded file on 2026-09-24.
    /// </summary>
    private const string Pxd010820Checksums =
        "File-Name\tFile-MD5Checksum\tFile-Size\n" +
        "qe2_2017may21_09_dw19_04.raw\tcb1a1b4adccb075b175af9b07f6b2b6b\t1221458769\n" +
        "180418_SPROT_Human_Iso_53BP1_AUP000005640.fasta\t5216c718f08eaec311a8c4f45fbda04b\t47663988\n" +
        "README.txt\t7d4afc88dcffe035d881278712c3ea83\t1889\n" +
        "qe2_2017may21_11_dw19_05.raw\t8768c7bff19d8c525193e6744cf249f2\t1336526850\n" +
        "qe2_2017may21_21_dw19_10.raw\tbbc66f6ca267c9d3c1be0703ab61d21c\t1110735315\n" +
        "qe2_2017may21_03_dw19_01.raw\td5f697b9dce4b44c4e51d16aa308aa55\t1240908256\n" +
        "qe2_2017may21_05_dw19_02.raw\t016abc6af4e536083afa4237066384a7\t1392071168\n" +
        "Maxquant_Output.zip\tb685aaa69017a7894ec1d9cfc622b239\t8070973\n" +
        "qe2_2017may21_17_dw19_08.raw\taa34a5ee25bf67cbb160cafdb68db3a2\t1250816103\n" +
        "qe2_2017may21_07_dw19_03.raw\tec06afff160e7b412666403cb27a87ee\t1467060051\n" +
        "qe2_2017may21_15_dw19_07.raw\t1d8a2583e1098dc3a9ec8f652f6bfaed\t1221179766\n" +
        "qe2_2017may21_19_dw19_09.raw\t160d12f89d032a809ecf29a61a81acc3\t1238775209\n" +
        "qe2_2017may21_25_dw19_12.raw\t0ed27eb8f21dc686a74c3617221f885f\t1177650650\n" +
        "qe2_2017may21_13_dw19_06.raw\t099245bfe17c147044271622090ba83b\t1361400163\n" +
        "qe2_2017may21_23_dw19_11.raw\t97c1672ffb3d4d56b75f7a6c20359645\t1434266826\n" +
        "PXD010820_community_annotated.sdrf.tsv\t42a1b437add074fceea43b8bd3214908\t5842\n";

    private const string Header = "File-Name\tFile-MD5Checksum\tFile-Size\n";
    private const string ProjectJson = "{ \"accession\": \"PXD010820\", \"title\": \"t\" }";

    private static HttpResponseMessage Text(string body, HttpStatusCode status = HttpStatusCode.OK) =>
        new(status) { Content = new StringContent(body, Encoding.UTF8, "text/plain") };

    private static HttpResponseMessage Json(string body) =>
        new(HttpStatusCode.OK) { Content = new StringContent(body, Encoding.UTF8, "application/json") };

    /// <summary>A full download body, sent with an ETag so a validator sidecar is written beside the partial.</summary>
    private static HttpResponseMessage Download(byte[] body)
    {
        var response = new HttpResponseMessage(HttpStatusCode.OK) { Content = new ByteArrayContent(body) };
        response.Headers.ETag = EntityTagHeaderValue.Parse("\"v1\"");
        return response;
    }

    private static PrideArchiveClient ClientOver(StubHandler handler) =>
        new(new HttpClient(handler)) { MaxRetries = 0 };

    private static PrideArchiveFile MakeFile(string fileName, string url) =>
        new()
        {
            FileName = fileName,
            FileSizeBytes = 999_999, // deliberately NOT the download size, as PRIDE's REST value often is not
            FileCategory = new CvParam("PRIDE", "PRIDE:0000404", "category", "RAW"),
            PublicFileLocations = new List<CvParam> { new("PRIDE", "PRIDE:0000000", "location", url) }
        };

    private static readonly byte[] Body = Encoding.ASCII.GetBytes("ten bytes!");
    private static string Md5Of(byte[] bytes) => Convert.ToHexString(MD5.HashData(bytes)).ToLowerInvariant();
    private const string Url = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/10/PXD010820/run1.raw";

    private string _tempDir;
    private string Destination => Path.Combine(_tempDir, "run1.raw");

    /// <summary>How many responses a scripted handler has given; tests whose answers change per attempt count with it.</summary>
    private int _responses;

    [SetUp]
    public void SetUp()
    {
        _tempDir = Path.Combine(Path.GetTempPath(), "PrideChecksumTests", Guid.NewGuid().ToString("N"));
        _responses = 0;
    }

    [TearDown]
    public void TearDown()
    {
        try { if (Directory.Exists(_tempDir)) Directory.Delete(_tempDir, recursive: true); }
        catch (IOException) { /* best effort */ }
    }

    // ---- GetFileChecksumsAsync -----------------------------------------------

    [Test]
    public async Task GetFileChecksumsAsync_ParsesTheLiveList()
    {
        var handler = new StubHandler(_ => Text(Pxd010820Checksums));
        using var client = ClientOver(handler);

        var checksums = await client.GetFileChecksumsAsync("PXD010820");

        var sdrf = checksums["PXD010820_community_annotated.sdrf.tsv"];
        var raw = checksums["qe2_2017may21_07_dw19_03.raw"];
        Assert.Multiple(() =>
        {
            Assert.That(checksums, Has.Count.EqualTo(16));
            Assert.That(sdrf.FileName, Is.EqualTo("PXD010820_community_annotated.sdrf.tsv"));
            Assert.That(sdrf.Md5, Is.EqualTo("42a1b437add074fceea43b8bd3214908"));
            Assert.That(sdrf.SizeBytes, Is.EqualTo(5842));
            Assert.That(raw.SizeBytes, Is.EqualTo(1_467_060_051)); // larger than int.MaxValue / 2: a long, not an int
            Assert.That(handler.RequestedUris.Single(),
                Is.EqualTo("https://www.ebi.ac.uk/pride/ws/archive/v3/files/checksum/PXD010820"));
            // PRIDE answers "Accept: application/json" with 406 on this route.
            Assert.That(handler.Requests.Single().Headers.Accept, Is.Empty);
        });
    }

    [Test]
    public async Task GetFileChecksumsAsync_KeysAreCaseSensitive()
    {
        using var client = ClientOver(new StubHandler(_ => Text(Pxd010820Checksums)));

        var checksums = await client.GetFileChecksumsAsync("PXD010820");

        Assert.That(checksums.ContainsKey("readme.txt"), Is.False);
        Assert.That(checksums.ContainsKey("README.txt"), Is.True);
    }

    [Test]
    public async Task GetFileChecksumsAsync_ReadsCrLfLineEnds()
    {
        using var client = ClientOver(new StubHandler(_ => Text(Pxd010820Checksums.Replace("\n", "\r\n"))));

        var checksums = await client.GetFileChecksumsAsync("PXD010820");

        Assert.That(checksums, Has.Count.EqualTo(16));
        Assert.That(checksums["README.txt"].SizeBytes, Is.EqualTo(1889));
    }

    [Test]
    public void GetFileChecksumsAsync_EmptyBodyForAnUnknownAccession_Throws()
    {
        // PRIDE answers an accession it does not know with 200 and zero bytes (PXD999999999, measured live).
        var handler = new StubHandler(r => r.RequestUri!.AbsolutePath.Contains("/files/checksum/")
            ? Text("")
            : Text("", HttpStatusCode.NotFound));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<MzLibException>(async () => await client.GetFileChecksumsAsync("PXD999999999"));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.Message, Does.Contain("has no project with accession 'PXD999999999'"));
            Assert.That(handler.RequestedUris, Is.EqualTo(new[]
            {
                "https://www.ebi.ac.uk/pride/ws/archive/v3/files/checksum/PXD999999999",
                "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD999999999"
            }));
        });
    }

    [Test]
    public async Task GetFileChecksumsAsync_EmptyBodyForAnExistingProject_ReturnsEmpty()
    {
        var handler = new StubHandler(r => r.RequestUri!.AbsolutePath.Contains("/files/checksum/")
            ? Text("")
            : Json(ProjectJson));
        using var client = ClientOver(handler);

        var checksums = await client.GetFileChecksumsAsync("PXD010820");

        Assert.That(checksums, Is.Empty);
        Assert.That(handler.Requests, Has.Count.EqualTo(2));
    }

    [Test]
    public async Task GetFileChecksumsAsync_HeaderOnly_ReturnsEmptyWithoutAnExistenceCheck()
    {
        var handler = new StubHandler(_ => Text(Header));
        using var client = ClientOver(handler);

        var checksums = await client.GetFileChecksumsAsync("PXD010820");

        Assert.That(checksums, Is.Empty);
        Assert.That(handler.Requests, Has.Count.EqualTo(1));
    }

    [TestCase("File-Name\tMD5\tFile-Size\nrun1.raw\t" + "0123456789abcdef0123456789abcdef\t10\n", "starts with", TestName = "WrongHeader")]
    [TestCase(Header + "run1.raw\t0123456789abcdef0123456789abcdef\n", "malformed row at line 2", TestName = "TwoFields")]
    [TestCase(Header + "run1.raw\t0123456789abcdef0123456789abcdef\t10\textra\n", "malformed row at line 2", TestName = "FourFields")]
    [TestCase(Header + "\t0123456789abcdef0123456789abcdef\t10\n", "malformed row at line 2", TestName = "NoFileName")]
    [TestCase(Header + "run1.raw\t0123456789abcdef\t10\n", "not 32 hexadecimal characters", TestName = "ShortMd5")]
    [TestCase(Header + "run1.raw\t0123456789abcdef0123456789abcdeg\t10\n", "not 32 hexadecimal characters", TestName = "NonHexMd5")]
    [TestCase(Header + "run1.raw\t0123456789abcdef0123456789abcdef\t-10\n", "not a whole number of bytes", TestName = "NegativeSize")]
    [TestCase(Header + "run1.raw\t0123456789abcdef0123456789abcdef\t1.5\n", "not a whole number of bytes", TestName = "FractionalSize")]
    [TestCase(Header + "run1.raw\t0123456789abcdef0123456789abcdef\t10\nrun1.raw\t0123456789abcdef0123456789abcdef\t10\n", "lists 'run1.raw' twice", TestName = "DuplicateName")]
    [TestCase(Header + "run1.raw\t0123456789abcdef0123456789abcdef\t10\n\nrun2.raw\t0123456789abcdef0123456789abcdef\t10\n", "empty line at line 3", TestName = "BlankLineMidList")]
    public void GetFileChecksumsAsync_MalformedList_ThrowsAndIsNotRetried(string body, string expectedMessage)
    {
        var handler = new StubHandler(_ => Text(body));
        using var client = new PrideArchiveClient(new HttpClient(handler)) { RetryDelay = (_, _) => Task.CompletedTask };

        var exception = Assert.ThrowsAsync<MzLibException>(async () => await client.GetFileChecksumsAsync("PXD010820"));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.Message, Does.Contain(expectedMessage));
            Assert.That(exception.Message, Does.Contain("'PXD010820'"));
            Assert.That(handler.Requests, Has.Count.EqualTo(1)); // a broken list is not an outage
        });
    }

    [Test]
    public async Task GetFileChecksumsAsync_TransientFailure_IsRetried()
    {
        int calls = 0;
        var handler = new StubHandler(_ => ++calls == 1 ? Text("", HttpStatusCode.ServiceUnavailable) : Text(Pxd010820Checksums));
        using var client = new PrideArchiveClient(new HttpClient(handler)) { RetryDelay = (_, _) => Task.CompletedTask };

        var checksums = await client.GetFileChecksumsAsync("PXD010820");

        Assert.That(checksums, Has.Count.EqualTo(16));
        Assert.That(handler.Requests, Has.Count.EqualTo(2));
    }

    [TestCase(null)]
    [TestCase("")]
    [TestCase("  ")]
    public void GetFileChecksumsAsync_BlankAccession_Throws(string accession)
    {
        var handler = new StubHandler(_ => Text(Pxd010820Checksums));
        using var client = ClientOver(handler);

        Assert.ThrowsAsync<ArgumentException>(async () => await client.GetFileChecksumsAsync(accession));
        Assert.That(handler.Requests, Is.Empty);
    }

    // ---- DownloadFileAsync with a checksum ------------------------------------

    [Test]
    public void DownloadWithChecksum_NullChecksum_Throws()
    {
        var handler = new StubHandler(_ => Download(Body));
        using var client = ClientOver(handler);

        var exception = Assert.ThrowsAsync<ArgumentNullException>(async () =>
            await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir, expected: null));

        Assert.That(exception!.ParamName, Is.EqualTo("expected"));
        Assert.That(handler.Requests, Is.Empty);
    }

    [Test]
    public void DownloadWithChecksum_ChecksumForAnotherFile_Throws()
    {
        var handler = new StubHandler(_ => Download(Body));
        using var client = ClientOver(handler);
        var other = new PrideFileChecksum("run2.raw", Md5Of(Body), Body.Length);

        var exception = Assert.ThrowsAsync<ArgumentException>(async () =>
            await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir, other));

        Assert.That(exception!.ParamName, Is.EqualTo("expected"));
        Assert.That(handler.Requests, Is.Empty);
    }

    [TestCase(false)]
    [TestCase(true)]
    public async Task DownloadWithChecksum_Matching_WritesTheFile(bool verifyMd5)
    {
        using var client = ClientOver(new StubHandler(_ => Download(Body)));

        string path = await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(Body).ToUpperInvariant(), Body.Length), verifyMd5: verifyMd5);

        Assert.That(path, Is.EqualTo(Destination));
        Assert.That(File.ReadAllBytes(path), Is.EqualTo(Body));
    }

    [Test]
    public void DownloadWithChecksum_WrongSize_IsRetriedFromZeroThenThrowsAndLeavesNothing()
    {
        // The server's own Content-Length agrees with the body, so only the checksum list can catch this.
        var handler = new StubHandler(_ => Download(Body));
        using var client = new PrideArchiveClient(new HttpClient(handler)) { RetryDelay = (_, _) => Task.CompletedTask };

        var exception = Assert.ThrowsAsync<MzLibException>(async () =>
            await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
                new PrideFileChecksum("run1.raw", Md5Of(Body), Body.Length + 1)));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.Message, Does.Contain("came to 10 bytes where PRIDE's checksum list records 11"));
            Assert.That(handler.Requests, Has.Count.EqualTo(4)); // the first attempt and three retries
            Assert.That(handler.Requests.Select(r => r.Headers.Range), Is.All.Null); // each from byte zero, never resumed
            Assert.That(File.Exists(Destination), Is.False);
            Assert.That(File.Exists(Destination + ".partial"), Is.False);
            Assert.That(File.Exists(Destination + ".partial" + PrideArchiveClient.ValidatorSuffix), Is.False);
        });
    }

    [Test]
    public void DownloadWithChecksum_WrongMd5_ThrowsWhenAskedToVerify()
    {
        using var client = ClientOver(new StubHandler(_ => Download(Body)));
        var expected = new PrideFileChecksum("run1.raw", new string('0', 32), Body.Length);

        var exception = Assert.ThrowsAsync<MzLibException>(async () =>
            await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir, expected, verifyMd5: true));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.Message, Does.Contain($"has MD5 {Md5Of(Body)} where PRIDE's checksum list records {new string('0', 32)}"));
            Assert.That(File.Exists(Destination), Is.False);
            Assert.That(File.Exists(Destination + ".partial"), Is.False);
        });
    }

    [Test]
    public async Task DownloadWithChecksum_WrongMd5_IsNotCheckedByDefault()
    {
        using var client = ClientOver(new StubHandler(_ => Download(Body)));
        var expected = new PrideFileChecksum("run1.raw", new string('0', 32), Body.Length);

        string path = await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir, expected);

        Assert.That(File.ReadAllBytes(path), Is.EqualTo(Body));
    }

    [Test]
    public async Task DownloadWithChecksum_NoOverwrite_MatchingFile_IsKeptWithoutARequest()
    {
        Directory.CreateDirectory(_tempDir);
        File.WriteAllBytes(Destination, Body);
        var handler = new StubHandler(_ => Download(Body));
        using var client = ClientOver(handler);

        string path = await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(Body), Body.Length), overwrite: false, verifyMd5: true);

        Assert.That(path, Is.EqualTo(Destination));
        Assert.That(handler.Requests, Is.Empty);
    }

    [Test]
    public async Task DownloadWithChecksum_NoOverwrite_TruncatedFile_IsDownloadedAgain()
    {
        Directory.CreateDirectory(_tempDir);
        File.WriteAllBytes(Destination, Body.Take(4).ToArray());
        var handler = new StubHandler(_ => Download(Body));
        using var client = ClientOver(handler);

        await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(Body), Body.Length), overwrite: false);

        Assert.That(handler.Requests, Has.Count.EqualTo(1));
        Assert.That(File.ReadAllBytes(Destination), Is.EqualTo(Body));
    }

    [TestCase(true, 1)]
    [TestCase(false, 0)]
    public async Task DownloadWithChecksum_NoOverwrite_RightSizeWrongBytes_IsReplacedOnlyWhenVerifyingMd5(bool verifyMd5, int expectedRequests)
    {
        byte[] corrupt = Encoding.ASCII.GetBytes("TEN BYTES!");
        Directory.CreateDirectory(_tempDir);
        File.WriteAllBytes(Destination, corrupt);
        var handler = new StubHandler(_ => Download(Body));
        using var client = ClientOver(handler);

        await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(Body), Body.Length), overwrite: false, verifyMd5: verifyMd5);

        Assert.That(handler.Requests, Has.Count.EqualTo(expectedRequests));
        Assert.That(File.ReadAllBytes(Destination), Is.EqualTo(verifyMd5 ? Body : corrupt));
    }

    [Test]
    public async Task DownloadWithChecksum_NoOverwrite_MatchingFileWithNoHttpsLocation_IsKeptWithoutThrowing()
    {
        // The skip still runs before URL resolution, so an Aspera-only file already on disk is not an error.
        Directory.CreateDirectory(_tempDir);
        File.WriteAllBytes(Destination, Body);
        var file = MakeFile("run1.raw", "prd_ascp@fasp.ebi.ac.uk:pride/data/archive/run1.raw");
        file.PublicFileLocations[0] = new CvParam("PRIDE", PrideArchiveExtensions.AsperaLocationAccession, "location",
            "prd_ascp@fasp.ebi.ac.uk:pride/data/archive/run1.raw");
        using var client = ClientOver(new StubHandler(_ => Download(Body)));

        string path = await client.DownloadFileAsync(file, _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(Body), Body.Length), overwrite: false);

        Assert.That(path, Is.EqualTo(Destination));
    }

    [Test]
    public void DownloadWithChecksum_Mismatch_NeverPutsUrlInMessage()
    {
        // The same rule as DownloadFileAsync_Failure_NeverPutsUrlInMessage: a reviewer-token href's query string
        // is the credential, so the mismatch message names the file and host only.
        const string secret = "S3CR3T-REVIEWER-TOKEN";
        using var client = ClientOver(new StubHandler(_ => Download(Body)));
        var file = MakeFile("run1.raw", $"https://private.example.org/files/run1.raw?token={secret}");

        var exception = Assert.ThrowsAsync<MzLibException>(async () =>
            await client.DownloadFileAsync(file, _tempDir, new PrideFileChecksum("run1.raw", new string('0', 32), Body.Length), verifyMd5: true));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.ToString(), Does.Not.Contain(secret));
            Assert.That(exception.Message, Does.Contain("'run1.raw' from private.example.org"));
        });
    }

    /// <summary>A body stream that cancels <see cref="Source"/> once the last byte has been read.</summary>
    private sealed class CancelAtEndStream : MemoryStream
    {
        public CancellationTokenSource Source { get; }
        public CancelAtEndStream(byte[] bytes, CancellationTokenSource source) : base(bytes) => Source = source;

        public override async ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
        {
            int read = await base.ReadAsync(buffer, cancellationToken);
            if (read == 0)
                Source.Cancel();
            return read;
        }
    }

    [Test]
    public void DownloadWithChecksum_CancelledWhileVerifying_KeepsThePartialAndValidator()
    {
        // The whole file arrives, then the caller cancels during the MD5. A multi-gigabyte transfer is not thrown away.
        using var cts = new CancellationTokenSource();
        var handler = new StubHandler(_ =>
        {
            var response = new HttpResponseMessage(HttpStatusCode.OK) { Content = new StreamContent(new CancelAtEndStream(Body, cts)) };
            response.Headers.ETag = EntityTagHeaderValue.Parse("\"v1\"");
            return response;
        });
        using var client = ClientOver(handler);

        Assert.CatchAsync<OperationCanceledException>(async () =>
            await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
                new PrideFileChecksum("run1.raw", Md5Of(Body), Body.Length), verifyMd5: true, cancellationToken: cts.Token));

        Assert.Multiple(() =>
        {
            Assert.That(File.Exists(Destination), Is.False);
            Assert.That(File.ReadAllBytes(Destination + ".partial"), Is.EqualTo(Body));
            Assert.That(File.Exists(Destination + ".partial" + PrideArchiveClient.ValidatorSuffix), Is.True);
        });
    }

    [Test]
    public async Task DownloadWithChecksum_CompletePartialFromACancelledCheck_IsCheckedWithoutARequest()
    {
        Directory.CreateDirectory(_tempDir);
        File.WriteAllBytes(Destination + ".partial", Body);
        File.WriteAllText(Destination + ".partial" + PrideArchiveClient.ValidatorSuffix, "\"v1\"");
        var handler = new StubHandler(_ => Download(Body));
        using var client = ClientOver(handler);

        string path = await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(Body), Body.Length), verifyMd5: true);

        Assert.Multiple(() =>
        {
            Assert.That(handler.Requests, Is.Empty);
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Body));
            Assert.That(File.Exists(Destination + ".partial"), Is.False);
            Assert.That(File.Exists(Destination + ".partial" + PrideArchiveClient.ValidatorSuffix), Is.False);
        });
    }

    [Test]
    public void DownloadWithChecksum_CompletePartialWithTheWrongMd5_IsDeleted()
    {
        Directory.CreateDirectory(_tempDir);
        File.WriteAllBytes(Destination + ".partial", Body);
        File.WriteAllText(Destination + ".partial" + PrideArchiveClient.ValidatorSuffix, "\"v1\"");
        using var client = ClientOver(new StubHandler(_ => Download(Body)));

        Assert.ThrowsAsync<MzLibException>(async () =>
            await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
                new PrideFileChecksum("run1.raw", new string('0', 32), Body.Length), verifyMd5: true));

        Assert.Multiple(() =>
        {
            Assert.That(File.Exists(Destination), Is.False);
            Assert.That(File.Exists(Destination + ".partial"), Is.False);
            Assert.That(File.Exists(Destination + ".partial" + PrideArchiveClient.ValidatorSuffix), Is.False);
        });
    }

    // ---- EBI's error text inside a served file (2026-10-08) ---------------------

    /// <summary>
    /// The block EBI's storage wrote over bytes of PXD012307's FR1_young_11_2.raw (at offset 518,854,917 of
    /// 2,777,280,429), copied verbatim from the damaged download PXReprise kept. The file kept its full size.
    /// </summary>
    private const string EbiErrorBlock =
        "<Error><Code>ConnectionClosedException</Code><Message>Premature end of Content-Length delimited message body " +
        "(expected: 2,284,117,933; received: 25,692,421)</Message><ErrorMessage/><RequestId/></Error>";

    /// <summary>A deterministic "raw file" of <paramref name="length"/> bytes with no '&lt;' in it.</summary>
    private static byte[] RawBytes(int length) => Enumerable.Range(0, length).Select(i => (byte)(i % 59)).ToArray();

    /// <summary><paramref name="clean"/> with <see cref="EbiErrorBlock"/> written over it at <paramref name="offset"/>; same length.</summary>
    private static byte[] Damaged(byte[] clean, int offset)
    {
        byte[] damaged = (byte[])clean.Clone();
        Encoding.ASCII.GetBytes(EbiErrorBlock).CopyTo(damaged, offset);
        return damaged;
    }

    private static PrideArchiveClient ClientWithRetries(StubHandler handler) =>
        new(new HttpClient(handler)) { RetryDelay = (_, _) => Task.CompletedTask };

    [Test]
    public void DownloadWithChecksum_EbiErrorTextOnEveryAttempt_ThrowsAfterTheRetriesAndLeavesNothing()
    {
        byte[] clean = RawBytes(4096);
        byte[] damaged = Damaged(clean, 1000);
        var handler = new StubHandler(_ => Download(damaged));
        using var client = ClientWithRetries(handler);

        // The size matches, so before this check the damaged file was returned as complete.
        var exception = Assert.ThrowsAsync<MzLibException>(async () =>
            await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
                new PrideFileChecksum("run1.raw", Md5Of(clean), clean.Length)));

        Assert.Multiple(() =>
        {
            Assert.That(exception!.Message, Is.EqualTo(
                "The PRIDE download of 'run1.raw' from ftp.pride.ebi.ac.uk holds a server error response (\"<Error><Code>\") " +
                "at byte 1000, written into the file in place of its own bytes."));
            Assert.That(handler.Requests, Has.Count.EqualTo(4));
            Assert.That(handler.Requests.Select(r => r.Headers.Range), Is.All.Null);
            Assert.That(File.Exists(Destination), Is.False);
            Assert.That(File.Exists(Destination + ".partial"), Is.False);
            Assert.That(File.Exists(Destination + ".partial" + PrideArchiveClient.ValidatorSuffix), Is.False);
        });
    }

    [Test]
    public async Task DownloadWithChecksum_EbiErrorTextOnce_IsDownloadedAgainFromZero()
    {
        byte[] clean = RawBytes(4096);
        byte[] damaged = Damaged(clean, 1000);
        var handler = new StubHandler(_ => Download(_responses++ == 0 ? damaged : clean));
        using var client = ClientWithRetries(handler);

        string path = await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(clean), clean.Length));

        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(clean));
            Assert.That(handler.Requests, Has.Count.EqualTo(2));
            Assert.That(handler.Requests[1].Headers.Range, Is.Null);
            Assert.That(File.Exists(Destination + ".partial" + PrideArchiveClient.ValidatorSuffix), Is.False);
        });
    }

    [Test]
    public async Task DownloadWithChecksum_WrongMd5Once_IsDownloadedAgainFromZero()
    {
        byte[] clean = RawBytes(4096);
        byte[] damaged = Damaged(clean, 1000);
        var handler = new StubHandler(_ => Download(_responses++ == 0 ? damaged : clean));
        using var client = ClientWithRetries(handler);

        string path = await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(clean), clean.Length), verifyMd5: true);

        Assert.That(File.ReadAllBytes(path), Is.EqualTo(clean));
        Assert.That(handler.Requests, Has.Count.EqualTo(2));
    }

    [Test]
    public async Task DownloadWithChecksum_MatchingMd5_OutranksTheErrorTextScan()
    {
        // PRIDE's own MD5 proves the bytes are the depositor's, even if they happen to hold the marker.
        byte[] legitimate = Damaged(RawBytes(4096), 1000);
        var handler = new StubHandler(_ => Download(legitimate));
        using var client = ClientWithRetries(handler);

        string path = await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(legitimate), legitimate.Length), verifyMd5: true);

        Assert.That(File.ReadAllBytes(path), Is.EqualTo(legitimate));
        Assert.That(handler.Requests, Has.Count.EqualTo(1));
    }

    [Test]
    public void DownloadWithChecksum_EbiErrorTextAcrossTwoReads_IsFound()
    {
        // The copy loop reads 81,920 bytes at a time; this block starts 5 bytes before the end of the first read.
        byte[] damaged = Damaged(RawBytes(200_000), 81_915);
        using var client = ClientOver(new StubHandler(_ => Download(damaged)));

        var exception = Assert.ThrowsAsync<MzLibException>(async () =>
            await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
                new PrideFileChecksum("run1.raw", Md5Of(damaged), damaged.Length)));

        Assert.That(exception!.Message, Does.Contain("at byte 81915,"));
    }

    [Test]
    public async Task DownloadWithoutChecksum_IsNotScanned()
    {
        // The overload without a checksum is unchanged: no check, no retry, the bytes as served.
        byte[] damaged = Damaged(RawBytes(4096), 1000);
        var handler = new StubHandler(_ => Download(damaged));
        using var client = ClientWithRetries(handler);

        string path = await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir);

        Assert.That(File.ReadAllBytes(path), Is.EqualTo(damaged));
        Assert.That(handler.Requests, Has.Count.EqualTo(1));
    }

    /// <summary>A body that delivers its first <paramref name="count"/> bytes and then drops, as EBI's connections do.</summary>
    private sealed class DroppingStream(byte[] bytes, int count) : MemoryStream(bytes, 0, count)
    {
        public override async ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
        {
            int read = await base.ReadAsync(buffer, cancellationToken);
            return read > 0 ? read : throw new IOException("The response ended prematurely.");
        }
    }

    [Test]
    public async Task DownloadWithChecksum_EbiErrorTextInTheResumedStart_IsFoundAndDownloadedAgainFromZero()
    {
        // Attempt 1 delivers the damaged start and drops; attempt 2 resumes the clean rest. The scan saw only the
        // rest, so the whole file is read again, the error text found, and attempt 3 fetches it all afresh.
        byte[] clean = RawBytes(4096);
        byte[] damaged = Damaged(clean, 1000);
        var handler = new StubHandler(request =>
        {
            switch (_responses++)
            {
                case 0:
                    var dropped = new HttpResponseMessage(HttpStatusCode.OK) { Content = new StreamContent(new DroppingStream(damaged, 2048)) };
                    dropped.Content.Headers.ContentLength = clean.Length;
                    dropped.Headers.ETag = EntityTagHeaderValue.Parse("\"v1\"");
                    return dropped;
                case 1:
                    Assert.That(request.Headers.Range!.Ranges.Single().From, Is.EqualTo(2048));
                    var rest = new HttpResponseMessage(HttpStatusCode.PartialContent) { Content = new ByteArrayContent(clean[2048..]) };
                    rest.Content.Headers.ContentRange = new ContentRangeHeaderValue(2048, clean.Length - 1, clean.Length);
                    rest.Headers.ETag = EntityTagHeaderValue.Parse("\"v1\"");
                    return rest;
                default:
                    return Download(clean);
            }
        });
        using var client = ClientWithRetries(handler);

        string path = await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(clean), clean.Length));

        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(clean));
            Assert.That(handler.Requests, Has.Count.EqualTo(3));
            Assert.That(handler.Requests[2].Headers.Range, Is.Null);
        });
    }

    [Test]
    public async Task DownloadWithChecksum_EbiErrorTextInAPartialFromAnEarlierCall_IsFoundAndDownloadedAgainFromZero()
    {
        // An earlier call left the damaged start on disk with its validator. This call resumes it, so the bytes it
        // streams are only the rest: the file must be read again to be checked.
        byte[] clean = RawBytes(4096);
        Directory.CreateDirectory(_tempDir);
        File.WriteAllBytes(Destination + ".partial", Damaged(clean, 1000)[..2048]);
        File.WriteAllText(Destination + ".partial" + PrideArchiveClient.ValidatorSuffix, "\"v1\"");
        var handler = new StubHandler(_ =>
        {
            if (_responses++ > 0)
                return Download(clean);
            var rest = new HttpResponseMessage(HttpStatusCode.PartialContent) { Content = new ByteArrayContent(clean[2048..]) };
            rest.Content.Headers.ContentRange = new ContentRangeHeaderValue(2048, clean.Length - 1, clean.Length);
            rest.Headers.ETag = EntityTagHeaderValue.Parse("\"v1\"");
            return rest;
        });
        using var client = ClientWithRetries(handler);

        string path = await client.DownloadFileAsync(MakeFile("run1.raw", Url), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(clean), clean.Length));

        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(clean));
            Assert.That(handler.Requests, Has.Count.EqualTo(2));
            Assert.That(handler.Requests[0].Headers.Range!.Ranges.Single().From, Is.EqualTo(2048));
            Assert.That(handler.Requests[1].Headers.Range, Is.Null);
        });
    }

    [Test]
    public void DownloadDigest_FindsTheMarkerWhereverTheReadsSplitIt()
    {
        byte[] marker = PrideArchiveClient.DownloadDigest.ServerErrorMarker;
        byte[] damaged = Damaged(RawBytes(240), 20);

        // Every split of the file into three reads, including reads shorter than the marker and empty ones.
        for (int first = 0; first <= damaged.Length; first++)
        for (int second = first; second <= damaged.Length; second++)
        {
            using var digest = new PrideArchiveClient.DownloadDigest(verifyMd5: false);
            digest.Append(damaged.AsSpan(0, first));
            digest.Append(damaged.AsSpan(first, second - first));
            digest.Append(damaged.AsSpan(second));
            Assert.That(digest.ServerErrorAt, Is.EqualTo(20), $"reads split at {first} and {second}");
        }

        using var cleanDigest = new PrideArchiveClient.DownloadDigest(verifyMd5: false);
        foreach (byte b in RawBytes(64).Concat(marker[..^1]))
            cleanDigest.Append(new[] { b });
        Assert.That(cleanDigest.ServerErrorAt, Is.EqualTo(-1), "an unfinished marker at the end is not a hit");
    }

    [Test]
    public void PrideFileChecksum_NullFileNameOrMd5_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => new PrideFileChecksum(null!, Md5Of(Body), 10), Throws.TypeOf<ArgumentNullException>()
                .With.Property(nameof(ArgumentNullException.ParamName)).EqualTo("fileName"));
            Assert.That(() => new PrideFileChecksum("run1.raw", null!, 10), Throws.TypeOf<ArgumentNullException>()
                .With.Property(nameof(ArgumentNullException.ParamName)).EqualTo("md5"));
        });
    }
}

/// <summary>
/// Live canaries for the checksum list against the real PRIDE Archive. An outage SKIPS through
/// <see cref="ExternalServiceTestHelper.RunAsync"/>; a changed list format or a wrong MD5 FAILS.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideChecksumLiveTests
{
    private const string SdrfName = "PXD010820_community_annotated.sdrf.tsv";

    [Test]
    public Task GetFileChecksumsAsync_Live_ReadsTheList() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            // one attempt: a PRIDE outage skips at once, not after the 5/20/60 s retry backoffs
            using var client = new PrideArchiveClient { MaxRetries = 0 };

            var checksums = await client.GetFileChecksumsAsync("PXD010820");

            Assert.That(checksums, Has.Count.GreaterThanOrEqualTo(16));
            Assert.That(checksums[SdrfName].Md5, Is.EqualTo("42a1b437add074fceea43b8bd3214908"));
            Assert.That(checksums[SdrfName].SizeBytes, Is.EqualTo(5842));
        });

    [Test]
    public Task GetFileChecksumsAsync_Live_UnknownAccessionThrows() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient { MaxRetries = 0 };

            // Not Assert.ThrowsAsync: it would turn an outage's HttpRequestException into a failure instead of
            // letting it reach RunAsync, which skips.
            try
            {
                await client.GetFileChecksumsAsync("PXD999999999");
                Assert.Fail("An unknown accession must throw, not read as a project with no checksums.");
            }
            catch (MzLibException e)
            {
                Assert.That(e.Message, Does.Contain("has no project"));
            }
        });

    [Test]
    public Task DownloadWithChecksum_Live_VerifiesTheMd5() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient { MaxRetries = 0 };
            var checksums = await client.GetFileChecksumsAsync("PXD010820");
            var file = (await client.GetProjectFilesAsync("PXD010820")).Single(f => f.FileName == SdrfName);

            string dir = Path.Combine(Path.GetTempPath(), "PrideLiveChecksum", Guid.NewGuid().ToString("N"));
            try
            {
                string path = await client.DownloadFileAsync(file, dir, checksums[SdrfName], verifyMd5: true);
                Assert.That(new FileInfo(path).Length, Is.EqualTo(5842));
            }
            finally
            {
                if (Directory.Exists(dir)) Directory.Delete(dir, recursive: true);
            }
        });
}
