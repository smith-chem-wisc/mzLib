using FluentFTP.Exceptions;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Net.Http.Headers;
using System.Net.Security;
using System.Net.Sockets;
using System.Security.Authentication;
using System.Security.Cryptography;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

/// <summary>
/// M1a: <see cref="MassiveClient"/> offline. PROXI answers come from a stub handler (the record fixtures are
/// MassIVE's own bytes, contacts emptied: MassiveFixtures/PROVENANCE.md), the FTP server is a scripted
/// <see cref="IMassiveFtpConnection"/>, and backoffs are recorded, not slept through.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class MassiveClientTests
{
    // ---- doubles ------------------------------------------------------------------------------------------------

    /// <summary>Answers every request with <c>respond</c> and records the URIs asked for.</summary>
    private sealed class ProxiHandler(Func<HttpRequestMessage, HttpResponseMessage> respond) : HttpMessageHandler
    {
        public List<string> Requests { get; } = new();

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            Requests.Add(request.RequestUri!.AbsoluteUri);
            return Task.FromResult(respond(request));
        }
    }

    /// <summary>
    /// A scripted FTP server: folders answer LIST with their lines, files answer SIZE, MDTM and RETR, and a failure
    /// or a reply can be queued for a verb (CONNECT, LIST, CWD, SIZE, MDTM, RETR).
    /// </summary>
    private sealed class FakeFtp
    {
        public readonly Dictionary<string, string[]> Folders = new(StringComparer.Ordinal);
        public readonly Dictionary<string, byte[]> Files = new(StringComparer.Ordinal);
        public readonly Dictionary<string, string> Stamps = new(StringComparer.Ordinal);
        public readonly List<string> Log = new();
        public Func<string, long, Stream> OpenOverride;
        public Func<string, string[]> ListFallback;
        public bool ThrowOnDispose;
        public int Connections;
        public int Disposals;

        private readonly Dictionary<string, Queue<Func<Exception>>> _failures = new(StringComparer.Ordinal);
        private readonly Dictionary<string, Queue<MassiveFtpReply>> _replies = new(StringComparer.Ordinal);

        public void Fail(string verb, int times, Func<Exception> failure)
        {
            if (!_failures.TryGetValue(verb, out Queue<Func<Exception>> queue))
                _failures[verb] = queue = new Queue<Func<Exception>>();
            for (int i = 0; i < times; i++)
                queue.Enqueue(failure);
        }

        public void Answer(string verb, int times, MassiveFtpReply reply)
        {
            if (!_replies.TryGetValue(verb, out Queue<MassiveFtpReply> queue))
                _replies[verb] = queue = new Queue<MassiveFtpReply>();
            for (int i = 0; i < times; i++)
                queue.Enqueue(reply);
        }

        private void ThrowIfScripted(string verb)
        {
            if (_failures.TryGetValue(verb, out Queue<Func<Exception>> queue) && queue.Count > 0)
                throw queue.Dequeue()();
        }

        public Task<IMassiveFtpConnection> ConnectAsync(CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            Connections++;
            Log.Add("CONNECT");
            ThrowIfScripted("CONNECT");
            return Task.FromResult<IMassiveFtpConnection>(new Connection(this));
        }

        private sealed class Connection(FakeFtp server) : IMassiveFtpConnection
        {
            public Task<IReadOnlyList<string>> ListAsync(string directory, CancellationToken cancellationToken)
            {
                server.Log.Add("LIST " + directory);
                server.ThrowIfScripted("LIST");
                if (!server.Folders.TryGetValue(directory, out string[] lines) && (lines = server.ListFallback?.Invoke(directory)) == null)
                    throw new MassiveFtpReplyException(450, "Could not list: file access error: couldn't stat: file does not exist");
                return Task.FromResult<IReadOnlyList<string>>(lines);
            }

            public Task<MassiveFtpReply> ExecuteAsync(string command, CancellationToken cancellationToken)
            {
                server.Log.Add(command);
                string verb = command.Split(' ')[0];
                string argument = command.Substring(verb.Length + 1);
                server.ThrowIfScripted(verb);
                if (server._replies.TryGetValue(verb, out Queue<MassiveFtpReply> queue) && queue.Count > 0)
                    return Task.FromResult(queue.Dequeue());

                MassiveFtpReply reply = verb switch
                {
                    "CWD" => server.Folders.ContainsKey(argument) ? new(250, "CD worked on " + argument) : new(550, "CD issue: file does not exist"),
                    "SIZE" => server.Files.TryGetValue(argument, out byte[] bytes)
                        ? new(213, bytes.Length.ToString(CultureInfo.InvariantCulture))
                        : new(550, $"Couldn't access {argument}: file does not exist"),
                    "MDTM" => server.Files.ContainsKey(argument)
                        ? new(213, server.Stamps.TryGetValue(argument, out string stamp) ? stamp : "20150212140541")
                        : new(550, $"Couldn't access {argument}: file does not exist"),
                    _ => new(502, "Command not implemented"),
                };
                return Task.FromResult(reply);
            }

            public Task<Stream> OpenReadAsync(string path, long offset, CancellationToken cancellationToken)
            {
                server.Log.Add($"RETR {path} @{offset}");
                server.ThrowIfScripted("RETR");
                if (server.OpenOverride != null)
                    return Task.FromResult(server.OpenOverride(path, offset));
                if (!server.Files.TryGetValue(path, out byte[] bytes))
                    throw new MassiveFtpReplyException(550, "Could not access file: driver error: calling GetHandle: file does not exist");
                return Task.FromResult<Stream>(new MemoryStream(bytes, (int)offset, bytes.Length - (int)offset, writable: false));
            }

            public void Dispose()
            {
                server.Disposals++;
                if (server.ThrowOnDispose)
                    throw new MassiveFtpTransportException("The FTP connection failed: QUIT timed out", new TimeoutException());
            }
        }
    }

    /// <summary>Yields the first <paramref name="count"/> bytes and then fails as a dropped FTP data connection does.</summary>
    private sealed class DroppingStream(byte[] data, int start, int count) : MemoryStream(data, start, count, writable: false)
    {
        public override async ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
        {
            int read = await base.ReadAsync(buffer, cancellationToken);
            return read > 0 ? read : throw new MassiveFtpTransportException("The FTP data connection failed: connection reset", new SocketException(10054));
        }
    }

    /// <summary>Never delivers a byte: a stalled transfer. The await must propagate the cancellation it is given.</summary>
    private sealed class StallingStream : MemoryStream
    {
        public override async ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
        {
            await Task.Delay(Timeout.Infinite, cancellationToken);
            return 0;
        }
    }

    // ---- shared data -------------------------------------------------------------------------------------------

    private static readonly byte[] Source = Enumerable.Range(0, 1000).Select(i => (byte)(i * 7 % 251)).ToArray();
    private static readonly DateTime Now = new(2026, 11, 10, 13, 0, 0, DateTimeKind.Utc);
    private const string Root = "/v04/MSV000090383";
    private const string RawFile = Root + "/raw/wGH0001.raw";

    // Synthetic, in the shape of MassIVE's record for MSV000090383 (its identifier and FTP location terms only).
    private const string Msv000090383Record = """
        {"accession":[{"cvLabel":"MS","name":"MassIVE dataset identifier","accession":"MS:1002487","value":"MSV000090383"}],
         "datasetLink":[{"cvLabel":"MS","name":"Dataset FTP location","accession":"MS:1002852","value":"ftp://massive.ucsd.edu/v04/MSV000090383"}]}
        """;

    private List<TimeSpan> _backoffs;
    private string _tempDir;

    [SetUp]
    public void SetUp()
    {
        _backoffs = new List<TimeSpan>();
        _tempDir = Path.Combine(Path.GetTempPath(), "MassiveClientTests", Guid.NewGuid().ToString("N"));
    }

    [TearDown]
    public void TearDown()
    {
        try { if (Directory.Exists(_tempDir)) Directory.Delete(_tempDir, recursive: true); }
        catch { /* best-effort cleanup */ }
    }

    private static byte[] Fixture(string name) =>
        File.ReadAllBytes(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "MassiveFixtures", name));

    /// <summary>A response carrying <paramref name="body"/> under MassIVE's own label, which is wrong for half its records.</summary>
    private static HttpResponseMessage Proxi(byte[] body, HttpStatusCode status = HttpStatusCode.OK)
    {
        var content = new ByteArrayContent(body);
        content.Headers.ContentType = MediaTypeHeaderValue.Parse("application/json;charset=ISO-8859-1");
        return new HttpResponseMessage(status) { Content = content };
    }

    private static HttpResponseMessage Status(HttpStatusCode status, string body = "") =>
        new(status) { Content = new StringContent(body) };

    private static readonly byte[] NotFoundBody = Encoding.ASCII.GetBytes("""{"code":404,"title":"Not Found","message":"No data found for the specified parameters."}""");

    private MassiveClient ClientOver(FakeFtp ftp, HttpMessageHandler proxi = null, int maxRetries = 3, TimeSpan? stall = null) =>
        new(new HttpClient(proxi ?? new ProxiHandler(_ => Proxi(Encoding.UTF8.GetBytes(Msv000090383Record)))), ftp.ConnectAsync)
        {
            MaxRetries = maxRetries,
            RetryDelay = (delay, _) =>
            {
                _backoffs.Add(delay);
                return Task.CompletedTask;
            },
            Clock = () => Now,
            BodyStallTimeout = stall ?? TimeSpan.FromMinutes(2),
        };

    private static MassiveFile RawEntry(long size = 1000) => new("MSV000090383", Root, "raw/wGH0001.raw", size);

    private static FakeFtp ServerWithRawFile()
    {
        var ftp = new FakeFtp();
        ftp.Files[RawFile] = Source;
        ftp.Folders[Root + "/raw"] = new[] { "-rwxr-xr-x 1 ftp ftp         1000 Sep 25  2022 wGH0001.raw" };
        return ftp;
    }

    /// <summary>
    /// MSV000090383's tree. The root's lines are verbatim from MassIVE's FTP server (2026-11-10), FluentFTP's
    /// "Lines captured:" header included; the folders below are synthetic, in the same format.
    /// </summary>
    private static FakeFtp ServerWithTree()
    {
        var ftp = new FakeFtp();
        ftp.Folders[Root] = new[]
        {
            "Lines captured:",
            "Lrwxrwxrwx 1 ftp ftp           40 Sep 24  2024 ccms_peak",
            "drwxr-xr-x 1 ftp ftp            3 Sep 25  2022 ccms_parameters",
            "drwxr-xr-x 1 ftp ftp            3 Sep 25  2022 raw",
            "drwxr-xr-x 1 ftp ftp            5 Sep 25  2022 library",
            "drwxr-xr-x 1 ftp ftp           15 Oct 13  2022 updates",
        };
        ftp.Folders[Root + "/ccms_parameters"] = new[] { "-rwxr-xr-x 1 ftp ftp         1794 Sep 25  2022 params.xml" };
        ftp.Folders[Root + "/raw"] = new[]
        {
            "drwxr-xr-x 1 ftp ftp            3 Sep 25  2022 .",
            "drwxr-xr-x 1 ftp ftp            5 Sep 25  2022 ..",
            "-rwxr-xr-x 1 ftp ftp         1000 Sep 25  2022 wGH0001.raw",
            "drwxr-xr-x 1 ftp ftp            1 Sep 25  2022 License File",
        };
        ftp.Folders[Root + "/raw/License File"] = new[] { "-rwxr-xr-x 1 ftp ftp          456 Feb 12  2015 license.txt" };
        ftp.Folders[Root + "/library"] = Array.Empty<string>();
        ftp.Folders[Root + "/updates"] = new[] { "drwxr-xr-x 1 ftp ftp            4 Sep 27  2022 2022-09-27_RGHuffman_25d0cd1a" };
        ftp.Folders[Root + "/updates/2022-09-27_RGHuffman_25d0cd1a"] = new[] { "drwxr-xr-x 1 ftp ftp           23 Sep 27  2022 raw" };
        ftp.Folders[Root + "/updates/2022-09-27_RGHuffman_25d0cd1a/raw"] = new[] { "-rwxr-xr-x 1 ftp ftp         1000 Nov 10 12:30 wGH0697.raw" };
        return ftp;
    }

    // ---- the record --------------------------------------------------------------------------------------------

    [Test]
    public void Record_OrdinaryRecord_ReadsEveryField()
    {
        MassiveDataset d = MassiveDataset.FromProxi(Fixture("MSV000096226.proxi.json"), "MSV000096226");

        Assert.Multiple(() =>
        {
            Assert.That(d.MassiveAccession, Is.EqualTo("MSV000096226"));
            Assert.That(d.PxdAccession, Is.EqualTo("PXD057243"));
            Assert.That(d.Title, Is.EqualTo("Remodeling of the protein ubiquitylation landscape in the aging vertebrate brain"));
            Assert.That(d.TitleHadPrefix, Is.False);
            Assert.That(d.Summary, Does.StartWith("Whole proteome data from brains of old mice"));
            Assert.That(d.Keywords, Has.Count.EqualTo(7));
            Assert.That(d.Keywords[0], Is.EqualTo("DR"));
            Assert.That(d.Species, Is.EqualTo(new[] { new MassiveSpecies("10090", "Mus musculus", null) }));
            Assert.That(d.Instruments.Select(i => i.Accession), Is.EqualTo(new[] { "MS:1003028" }));
            Assert.That(d.Modifications.Select(m => m.Accession), Is.EqualTo(new[] { "UNIMOD:35", "UNIMOD:4" }));
            Assert.That(d.Publications, Is.Not.Null.And.Empty, "'no associated published manuscript' is stated as nothing, not unstated");
            Assert.That(d.DatasetUri, Is.EqualTo("https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000096226"));
            Assert.That(d.FtpLocation, Is.EqualTo("v07/MSV000096226"));
        });
    }

    [Test]
    public void Record_StringNullAndMissingFields_ReadAsUnstated_AnAllNullGroupAsEmpty()
    {
        MassiveDataset d = MassiveDataset.FromProxi(Fixture("MSV000078494.proxi.json"), "MSV000078494");

        Assert.Multiple(() =>
        {
            Assert.That(d.PxdAccession, Is.Null);
            Assert.That(d.Keywords, Is.Null, "the record has no keywords field");
            Assert.That(d.Modifications, Is.Null, "the record has no modifications field");
            Assert.That(d.Species, Is.Not.Null.And.Empty, "its only species group says \"null\", so the field is stated and the group dropped");
            Assert.That(d.Instruments, Is.Not.Null.And.Empty, "its only instrument is MS:1000031 = \"null\"");
            Assert.That(d.Publications, Has.Count.EqualTo(1));
            Assert.That(d.Publications[0].PubMedId, Is.EqualTo("22159717"));
            Assert.That(d.Publications[0].Doi, Is.Null);
            Assert.That(d.Publications[0].Reference, Does.StartWith("Naba A, Clauser KR"));
        });
    }

    [Test]
    public void Record_Latin1Bytes_UnderTheIsoLabel_DecodeAsWindows1252()
    {
        MassiveDataset d = MassiveDataset.FromProxi(Fixture("MSV000088323.proxi.json"), "MSV000088323");

        Assert.That(d.Summary, Does.Contain("Rappé, M. S."), "the one non-ASCII byte, 0xE9, is in the summary's citation list");
    }

    [Test]
    public void Record_Utf8Bytes_UnderTheIsoLabel_DecodeAsUtf8_NotByTheLabel()
    {
        MassiveDataset d = MassiveDataset.FromProxi(Fixture("MSV000079927.proxi.json"), "MSV000079927");

        Assert.Multiple(() =>
        {
            Assert.That(d.Summary, Does.Contain("CAUTI´s"));
            Assert.That(d.Summary, Does.Not.Contain("Â"), "decoded by its ISO-8859-1 label, the accent would become two characters");
            Assert.That(d.Species, Has.Count.EqualTo(4));
        });
    }

    [Test]
    public void Record_TitlePrefix_IsCutAndFlagged()
    {
        MassiveDataset d = MassiveDataset.FromProxi(Fixture("MSV000097834.proxi.json"), "MSV000097834");

        Assert.Multiple(() =>
        {
            Assert.That(d.TitleHadPrefix, Is.True);
            Assert.That(d.Title, Does.Not.StartWith("Title:"));
            Assert.That(d.Title, Is.Not.Empty);
        });
    }

    [Test]
    public void Record_InstrumentModelTerm_IsKeptWhenItNamesOne()
    {
        // Synthetic, in PROXI's shape: 7 of the 8 MS:1000031 terms in 121 records carry a real name like this one.
        const string json = """
            {"accession":[{"cvLabel":"MS","name":"MassIVE dataset identifier","accession":"MS:1002487","value":"MSV000000001"}],
             "instruments":[{"cvLabel":"MS","name":"instrument model","accession":"MS:1000031","value":"Bruker timsTOF Pro"},
                            {"cvLabel":"MS","name":"instrument model","accession":"MS:1000031","value":"null"}]}
            """;

        MassiveDataset d = MassiveDataset.FromProxi(Encoding.UTF8.GetBytes(json), "MSV000000001");

        Assert.That(d.Instruments.Select(i => i.Value), Is.EqualTo(new[] { "Bruker timsTOF Pro" }));
    }

    [Test]
    public void Record_PresentEmptyLists_AreEmpty_MissingOnesAreNull()
    {
        const string json = """
            {"accession":[{"cvLabel":"MS","name":"MassIVE dataset identifier","accession":"MS:1002487","value":"MSV000000001"}],
             "keywords":[], "species":[], "title":"", "summary":null}
            """;

        MassiveDataset d = MassiveDataset.FromProxi(Encoding.UTF8.GetBytes(json), "MSV000000001");

        Assert.Multiple(() =>
        {
            Assert.That(d.Keywords, Is.Not.Null.And.Empty);
            Assert.That(d.Species, Is.Not.Null.And.Empty);
            Assert.That(d.Title, Is.EqualTo(string.Empty), "stated as empty");
            Assert.That(d.Summary, Is.Null, "JSON null is unstated");
            Assert.That(d.Instruments, Is.Null);
            Assert.That(d.Publications, Is.Null);
            Assert.That(d.FtpLocation, Is.Null);
        });
    }

    [TestCase("[]")]
    [TestCase("{}")]
    [TestCase("{\"accession\":[{\"accession\":\"MS:1001919\",\"value\":\"PXD000001\"}]}")]
    [TestCase("not json")]
    public void Record_NotADatasetRecord_IsAContractBreak(string body)
    {
        Assert.Throws<MzLibException>(() => MassiveDataset.FromProxi(Encoding.UTF8.GetBytes(body), "MSV000000001"));
    }

    // ---- PROXI through the client ------------------------------------------------------------------------------

    [Test]
    public async Task TryGetDataset_UnknownAccession_IsNotFound_WithOneRequest()
    {
        var handler = new ProxiHandler(_ => Proxi(NotFoundBody, HttpStatusCode.NotFound));
        using var client = ClientOver(new FakeFtp(), handler);

        (bool found, MassiveDataset dataset) = await client.TryGetDatasetAsync("MSV999999999");

        Assert.Multiple(() =>
        {
            Assert.That(found, Is.False);
            Assert.That(dataset, Is.Null);
            Assert.That(handler.Requests, Has.Count.EqualTo(1), "a 404 is an answer, not retried, and ProteomeCentral is not asked");
        });
    }

    [Test]
    public async Task GetDataset_ThroughTheClient_IgnoresTheCharsetLabel()
    {
        var handler = new ProxiHandler(_ => Proxi(Fixture("MSV000079927.proxi.json")));
        using var client = ClientOver(new FakeFtp(), handler);

        MassiveDataset d = await client.GetDatasetAsync("PXD000164");

        Assert.Multiple(() =>
        {
            Assert.That(d.MassiveAccession, Is.EqualTo("MSV000079927"), "a PXD resolves to its MSV");
            Assert.That(d.Summary, Does.Contain("CAUTI´s"));
            Assert.That(handler.Requests.Single(), Is.EqualTo(MassiveClient.DefaultProxiBaseAddress + "datasets/PXD000164"));
        });
    }

    [Test]
    public void GetDataset_UnknownMsv_IsAnMzLibException()
    {
        using var client = ClientOver(new FakeFtp(), new ProxiHandler(_ => Proxi(NotFoundBody, HttpStatusCode.NotFound)));

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetAsync("MSV999999999"));

        Assert.That(e!.Message, Does.Contain("MassIVE has no dataset with accession 'MSV999999999'"));
    }

    [Test]
    public void GetDataset_PxdHostedElsewhere_NamesTheHost_FromProteomeCentral()
    {
        // datasetSummary as ProteomeCentral's PROXI sent it for PXD000001 (2026-11-10); the rest of its record is omitted.
        const string proteomeCentral = """{"datasetSummary":{"announceDate":"2024-10-07","hostingRepository":"PRIDE"}}""";
        var handler = new ProxiHandler(r => r.RequestUri!.Host == "massive.ucsd.edu"
            ? Proxi(NotFoundBody, HttpStatusCode.NotFound)
            : new HttpResponseMessage(HttpStatusCode.OK) { Content = new StringContent(proteomeCentral) });
        using var client = ClientOver(new FakeFtp(), handler);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetAsync("PXD000001"));

        Assert.Multiple(() =>
        {
            Assert.That(e!.Message, Does.Contain("'PXD000001' is hosted at PRIDE, not MassIVE"));
            Assert.That(handler.Requests[1], Is.EqualTo(MassiveClient.ProteomeCentralProxiBaseAddress + "datasets/PXD000001"));
        });
    }

    [Test]
    public void GetDataset_PxdNowhere_SaysSo()
    {
        var handler = new ProxiHandler(_ => Proxi(NotFoundBody, HttpStatusCode.NotFound));
        using var client = ClientOver(new FakeFtp(), handler);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetAsync("PXD999999"));

        Assert.That(e!.Message, Does.Contain("Neither MassIVE nor ProteomeCentral has a dataset with accession 'PXD999999'"));
    }

    [Test]
    public void GetDataset_ProteomeCentralDown_IsAskedOnce_AndTheAnswerIsStillMassIVEs()
    {
        var handler = new ProxiHandler(r => r.RequestUri!.Host == "massive.ucsd.edu"
            ? Proxi(NotFoundBody, HttpStatusCode.NotFound)
            : Status(HttpStatusCode.ServiceUnavailable));
        using var client = ClientOver(new FakeFtp(), handler);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetAsync("PXD000001"));

        Assert.Multiple(() =>
        {
            Assert.That(e!.Message, Does.Contain("ProteomeCentral could not be asked"));
            Assert.That(handler.Requests, Has.Count.EqualTo(2), "ProteomeCentral only names the host, so it is never retried");
            Assert.That(_backoffs, Is.Empty);
        });
    }

    [Test]
    public void GetDataset_ProxiOutage_IsRetried_ThenTransient_AndKeepsThePhrase()
    {
        var handler = new ProxiHandler(_ => Status(HttpStatusCode.ServiceUnavailable));
        using var client = ClientOver(new FakeFtp(), handler);

        var e = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetDatasetAsync("MSV000096226"));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(e), Is.True);
            Assert.That(e!.StatusCode, Is.EqualTo(HttpStatusCode.ServiceUnavailable));
            Assert.That(e.Message, Does.Contain("failed with status 503"));
            Assert.That(handler.Requests, Has.Count.EqualTo(4));
        });
    }

    [Test]
    public void GetDataset_ProxiRefusal_IsNotRetried_NorTransient()
    {
        var handler = new ProxiHandler(_ => Status(HttpStatusCode.BadRequest));
        using var client = ClientOver(new FakeFtp(), handler);

        var e = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetDatasetAsync("MSV000096226"));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(e), Is.False);
            Assert.That(e!.StatusCode, Is.EqualTo(HttpStatusCode.BadRequest));
            Assert.That(handler.Requests, Has.Count.EqualTo(1));
        });
    }

    [TestCase("")]
    [TestCase("   ")]
    [TestCase("MSV 000096226")]
    [TestCase("../MSV000096226")]
    [TestCase("000096226")]
    public void GetDataset_NotAnAccession_IsRefusedBeforeAnyRequest(string accession)
    {
        var handler = new ProxiHandler(_ => Status(HttpStatusCode.OK));
        using var client = ClientOver(new FakeFtp(), handler);

        Assert.ThrowsAsync<ArgumentException>(async () => await client.GetDatasetAsync(accession));
        Assert.That(handler.Requests, Is.Empty);
    }

    // ---- the listing -------------------------------------------------------------------------------------------

    [Test]
    public async Task GetDatasetFiles_WalksTheTree_ListsTheLink_AndReadsUpdatesFromThePath()
    {
        FakeFtp ftp = ServerWithTree();
        using var client = ClientOver(ftp);

        IReadOnlyList<MassiveFile> files = await client.GetDatasetFilesAsync("MSV000090383");

        Assert.That(files.Select(f => f.Path), Is.EqualTo(new[]
        {
            "ccms_peak",
            "ccms_parameters/params.xml",
            "raw/wGH0001.raw",
            "raw/License File/license.txt",
            "updates/2022-09-27_RGHuffman_25d0cd1a/raw/wGH0697.raw",
        }));

        MassiveFile link = files[0];
        MassiveFile update = files[4];
        Assert.Multiple(() =>
        {
            Assert.That(link.IsLink, Is.True, "MLSD would have called it a 40-byte file; FluentFTP's LIST parser would have dropped it");
            Assert.That(link.SizeBytes, Is.Null);
            Assert.That(link.Collection, Is.EqualTo("ccms_peak"));
            Assert.That(files[2].SizeBytes, Is.EqualTo(1000));
            Assert.That(files[2].Collection, Is.EqualTo("raw"));
            Assert.That(files[2].UpdateName, Is.Null);
            Assert.That(files[2].Modified, Is.EqualTo(new DateTime(2022, 9, 25)));
            Assert.That(files[3].FileName, Is.EqualTo("license.txt"));
            Assert.That(update.UpdateName, Is.EqualTo("2022-09-27_RGHuffman_25d0cd1a"));
            Assert.That(update.Collection, Is.EqualTo("raw"));
            Assert.That(update.Modified, Is.EqualTo(new DateTime(2026, 11, 10, 12, 30, 0)), "a time with no year is this year when that is not in the future");
            Assert.That(update.ServerPath, Is.EqualTo(Root + "/updates/2022-09-27_RGHuffman_25d0cd1a/raw/wGH0697.raw"));
            Assert.That(ftp.Log, Does.Not.Contain("LIST " + Root + "/ccms_peak"), "a link is never followed");
            Assert.That(ftp.Connections, Is.EqualTo(1), "one connection serves the whole listing");
        });
    }

    [Test]
    public void GetDatasetFiles_RootTheServerLacks_IsAnMzLibException_NotRetried()
    {
        var ftp = new FakeFtp(); // no folders at all: CWD answers 550, a LIST would answer 450
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetFilesAsync("MSV000090383"));

        Assert.Multiple(() =>
        {
            Assert.That(e!.Message, Does.Contain("has no folder '/v04/MSV000090383'"));
            Assert.That(ftp.Log.Count(l => l.StartsWith("CWD")), Is.EqualTo(1));
            Assert.That(ftp.Log.Any(l => l.StartsWith("LIST")), Is.False, "LIST would have said 450, FTP's 'try again later'");
            Assert.That(_backoffs, Is.Empty);
        });
    }

    [Test]
    public void GetDatasetFiles_EmptyRoot_IsRetried_ThenTransient()
    {
        var ftp = new FakeFtp();
        ftp.Folders[Root] = new[] { "Lines captured:" };
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetDatasetFilesAsync("MSV000090383"));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(e), Is.True);
            Assert.That(e!.Message, Does.Contain("came back empty"));
            Assert.That(ftp.Log.Count(l => l.StartsWith("LIST")), Is.EqualTo(4));
        });
    }

    [Test]
    public async Task GetDatasetFiles_TransientReplyOnce_IsRetried_AndTheListingCompletes()
    {
        FakeFtp ftp = ServerWithTree();
        ftp.Fail("LIST", 1, () => new MassiveFtpReplyException(450, "Could not list: file access error"));
        using var client = ClientOver(ftp);

        IReadOnlyList<MassiveFile> files = await client.GetDatasetFilesAsync("MSV000090383");

        Assert.Multiple(() =>
        {
            Assert.That(files, Has.Count.EqualTo(5));
            Assert.That(_backoffs, Is.EqualTo(new[] { TimeSpan.FromSeconds(5) }));
            Assert.That(ftp.Connections, Is.EqualTo(2), "a failed command drops the connection, and the retry opens a new one");
        });
    }

    [Test]
    public void GetDatasetFiles_TransientReplyThatPersists_IsTransient()
    {
        FakeFtp ftp = ServerWithTree();
        ftp.Fail("LIST", 100, () => new MassiveFtpReplyException(421, "Service not available, closing control connection"));
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetDatasetFilesAsync("MSV000090383"));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(e), Is.True);
            Assert.That(e!.StatusCode, Is.Null, "an FTP code is never passed off as an HTTP status");
            Assert.That(e.Message, Does.Contain("answered 421"));
            Assert.That(_backoffs, Has.Count.EqualTo(3));
        });
    }

    [Test]
    public void GetDatasetFiles_TransportFailureThatPersists_IsTransient()
    {
        FakeFtp ftp = ServerWithTree();
        ftp.Fail("CONNECT", 100, () => new MassiveFtpTransportException("The FTP connection failed: timed out", new TimeoutException()));
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<HttpRequestException>(async () => await client.GetDatasetFilesAsync("MSV000090383"));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(e), Is.True);
            Assert.That(e!.InnerException, Is.InstanceOf<MassiveFtpTransportException>());
            Assert.That(ftp.Connections, Is.EqualTo(4));
        });
    }

    [Test]
    public void GetDatasetFiles_PermanentReply_IsAnMzLibException_NotRetried()
    {
        FakeFtp ftp = ServerWithTree();
        ftp.Fail("LIST", 100, () => new MassiveFtpReplyException(530, "Not logged in"));
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetFilesAsync("MSV000090383"));

        Assert.Multiple(() =>
        {
            Assert.That(e!.Message, Does.Contain("refused").And.Contain("530"));
            Assert.That(_backoffs, Is.Empty);
        });
    }

    [Test]
    public void GetDatasetFiles_CertificateRefused_IsAnMzLibException_NotRetried()
    {
        var ftp = new FakeFtp();
        ftp.Fail("CONNECT", 100, () => new AuthenticationException("The remote certificate is invalid because of errors in the certificate chain"));
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetFilesAsync("MSV000090383"));

        Assert.Multiple(() =>
        {
            Assert.That(e!.Message, Does.Contain("TLS certificate check"));
            Assert.That(ftp.Connections, Is.EqualTo(1));
            Assert.That(ServiceFailure.IsTransient(e), Is.False);
        });
    }

    [Test]
    public void GetDatasetFiles_UnreadableLine_IsRefused_NotSkipped()
    {
        FakeFtp ftp = ServerWithTree();
        ftp.Folders[Root + "/raw"] = new[] { "type=file;size=1000; wGH0001.raw" };
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetFilesAsync("MSV000090383"));

        Assert.That(e!.Message, Does.Contain("cannot be read").And.Contain("type=file;size=1000; wGH0001.raw"));
    }

    [TestCase(null)]
    [TestCase("ftp://massive.ucsd.edu/v04/../etc")]
    public void GetDatasetFiles_RecordWithNoUsableFtpLocation_IsAnMzLibException(string location)
    {
        string record = location == null
            ? """{"accession":[{"accession":"MS:1002487","value":"MSV000090383"}]}"""
            : $$"""{"accession":[{"accession":"MS:1002487","value":"MSV000090383"}],"datasetLink":[{"accession":"MS:1002852","value":"{{location}}"}]}""";
        var ftp = new FakeFtp();
        using var client = ClientOver(ftp, new ProxiHandler(_ => Proxi(Encoding.UTF8.GetBytes(record))));

        Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetFilesAsync("MSV000090383"));
        Assert.That(ftp.Connections, Is.Zero);
    }

    [Test]
    public void GetDatasetFiles_OnlyEmptyFolders_IsAnMzLibException()
    {
        var ftp = new FakeFtp();
        ftp.Folders[Root] = new[] { "drwxr-xr-x 1 ftp ftp            3 Sep 25  2022 raw" };
        ftp.Folders[Root + "/raw"] = Array.Empty<string>();
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetFilesAsync("MSV000090383"));

        Assert.That(e!.Message, Does.Contain("only empty folders"));
    }

    [Test]
    public void Listing_ReadsEveryKind_AndAStandardLinkArrow()
    {
        List<MassiveFtpEntry> entries = MassiveFtpListing.Parse(new[]
        {
            "Lines captured:",
            "total 3",
            "Lrwxrwxrwx 1 ftp ftp           40 Sep 24  2024 ccms_peak",
            "lrwxrwxrwx 1 ftp ftp           12 Jan  5 09:15 latest -> raw/a.raw",
            "-rwxr-xr-x 1 ftp ftp          456 Feb 12  2015 license file.txt",
            "drwxr-xr-x 1 ftp ftp            3 Sep 25  2022 raw",
        }, "/v04/MSV000090383", Now);

        Assert.That(entries, Is.EqualTo(new[]
        {
            new MassiveFtpEntry("ccms_peak", MassiveFtpEntryKind.Link, 40, new DateTime(2024, 9, 24)),
            new MassiveFtpEntry("latest", MassiveFtpEntryKind.Link, 12, new DateTime(2026, 1, 5, 9, 15, 0)),
            new MassiveFtpEntry("license file.txt", MassiveFtpEntryKind.File, 456, new DateTime(2015, 2, 12)),
            new MassiveFtpEntry("raw", MassiveFtpEntryKind.Directory, 3, new DateTime(2022, 9, 25)),
        }));
    }

    [Test]
    public void Listing_TimeLaterThanNow_IsLastYear()
    {
        MassiveFtpEntry entry = MassiveFtpListing.Parse(new[] { "-rw-r--r-- 1 ftp ftp 1 Dec 30 23:00 a.raw" }, "/x", Now).Single();

        Assert.That(entry.Modified, Is.EqualTo(new DateTime(2025, 12, 30, 23, 0, 0)));
    }

    // ---- downloads ---------------------------------------------------------------------------------------------

    [Test]
    public async Task Download_KeepsThePathBelowTheDestination_AndLeavesNoPartial()
    {
        FakeFtp ftp = ServerWithRawFile();
        using var client = ClientOver(ftp);

        string path = await client.DownloadFileAsync(RawEntry(), _tempDir);

        Assert.Multiple(() =>
        {
            Assert.That(path, Is.EqualTo(Path.Combine(_tempDir, "raw", "wGH0001.raw")));
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
            Assert.That(File.Exists(path + ".partial"), Is.False);
            Assert.That(File.Exists(path + ".partial.validator"), Is.False);
            Assert.That(ftp.Log, Is.EqualTo(new[] { "CONNECT", "SIZE " + RawFile, "MDTM " + RawFile, $"RETR {RawFile} @0" }));
        });
    }

    [Test]
    public async Task Download_SameLeafInTwoFolders_DoNotOverwriteEachOther()
    {
        var ftp = new FakeFtp();
        byte[] other = Source.Reverse().ToArray();
        ftp.Files[Root + "/raw/a/x.raw"] = Source;
        ftp.Files[Root + "/raw/b/x.raw"] = other;
        using var client = ClientOver(ftp);

        string a = await client.DownloadFileAsync(new MassiveFile("MSV000090383", Root, "raw/a/x.raw", 1000), _tempDir);
        string b = await client.DownloadFileAsync(new MassiveFile("MSV000090383", Root, "raw/b/x.raw", 1000), _tempDir);

        Assert.That(a, Is.Not.EqualTo(b));
        Assert.That(File.ReadAllBytes(a), Is.EqualTo(Source));
        Assert.That(File.ReadAllBytes(b), Is.EqualTo(other));
    }

    [Test]
    public async Task Download_NoOverwrite_KeepsAFileOfTheListedSize_WithNoRequest_AndReplacesAnyOther()
    {
        FakeFtp ftp = ServerWithRawFile();
        using var client = ClientOver(ftp);
        string path = Path.Combine(_tempDir, "raw", "wGH0001.raw");
        Directory.CreateDirectory(Path.GetDirectoryName(path)!);

        File.WriteAllBytes(path, new byte[1000]);
        await client.DownloadFileAsync(RawEntry(), _tempDir, overwrite: false);
        Assert.That(ftp.Connections, Is.Zero, "a file of the listed size is kept with no request");

        File.WriteAllBytes(path, new byte[10]);
        await client.DownloadFileAsync(RawEntry(), _tempDir, overwrite: false);
        Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source), "a file of another size is downloaded again");
    }

    [Test]
    public void Download_Link_IsRefused_WithNoRequest()
    {
        var ftp = new FakeFtp();
        using var client = ClientOver(ftp);

        Assert.ThrowsAsync<NotSupportedException>(async () =>
            await client.DownloadFileAsync(new MassiveFile("MSV000090383", Root, "ccms_peak", null, isLink: true), _tempDir));
        Assert.That(ftp.Connections, Is.Zero);
    }

    [TestCase("raw/../../escape.raw")]
    [TestCase("raw/./x.raw")]
    [TestCase("raw//x.raw")]
    [TestCase("raw/a\\..\\..\\x.raw")]
    public void Download_PathThatCouldEscape_IsRefusedBeforeAnyRequest(string path)
    {
        var ftp = new FakeFtp();
        using var client = ClientOver(ftp);

        Assert.ThrowsAsync<ArgumentException>(async () =>
            await client.DownloadFileAsync(new MassiveFile("MSV000090383", Root, path, 1000), _tempDir));
        Assert.That(ftp.Connections, Is.Zero);
    }

    [Test]
    public async Task Download_DroppedTransfer_ResumesWithRest_InTheSameCall()
    {
        FakeFtp ftp = ServerWithRawFile();
        int opens = 0;
        ftp.OpenOverride = (_, offset) => opens++ == 0
            ? new DroppingStream(Source, 0, 400)
            : new MemoryStream(Source, (int)offset, Source.Length - (int)offset, writable: false);
        using var client = ClientOver(ftp);

        string path = await client.DownloadFileAsync(RawEntry(), _tempDir);

        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
            Assert.That(ftp.Log.Where(l => l.StartsWith("RETR")), Is.EqualTo(new[] { $"RETR {RawFile} @0", $"RETR {RawFile} @400" }));
            Assert.That(ftp.Connections, Is.EqualTo(2), "the dropped transfer's connection was replaced");
        });
    }

    [Test]
    public async Task Download_DroppedTransfer_KeepsThePartial_AndTheNextCallResumes()
    {
        FakeFtp ftp = ServerWithRawFile();
        ftp.OpenOverride = (_, _) => new DroppingStream(Source, 0, 400);
        using var client = ClientOver(ftp, maxRetries: 0);
        string path = Path.Combine(_tempDir, "raw", "wGH0001.raw");

        var e = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir));
        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(e), Is.True);
            Assert.That(e!.Message, Does.Contain("ended before the download completed"));
            Assert.That(new FileInfo(path + ".partial").Length, Is.EqualTo(400));
            Assert.That(File.ReadAllText(path + ".partial.validator"), Is.EqualTo("ftp size=1000 mdtm=20150212140541"));
        });

        ftp.OpenOverride = null;
        await client.DownloadFileAsync(RawEntry(), _tempDir);

        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
            Assert.That(ftp.Log.Last(), Is.EqualTo($"RETR {RawFile} @400"));
        });
    }

    [Test]
    public async Task Download_FileChangedSinceThePartial_StartsAgainFromZero()
    {
        FakeFtp ftp = ServerWithRawFile();
        string path = Path.Combine(_tempDir, "raw", "wGH0001.raw");
        Directory.CreateDirectory(Path.GetDirectoryName(path)!);
        File.WriteAllBytes(path + ".partial", Source.Take(400).ToArray());
        File.WriteAllText(path + ".partial.validator", "ftp size=1000 mdtm=20140101000000");
        using var client = ClientOver(ftp);

        await client.DownloadFileAsync(RawEntry(), _tempDir);

        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
            Assert.That(ftp.Log.Last(), Is.EqualTo($"RETR {RawFile} @0"), "the MDTM no longer matches");
        });
    }

    [Test]
    public async Task Download_PartialWithNoValidator_IsNotTrusted()
    {
        FakeFtp ftp = ServerWithRawFile();
        string path = Path.Combine(_tempDir, "raw", "wGH0001.raw");
        Directory.CreateDirectory(Path.GetDirectoryName(path)!);
        File.WriteAllBytes(path + ".partial", new byte[400]);
        using var client = ClientOver(ftp);

        await client.DownloadFileAsync(RawEntry(), _tempDir);

        Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
        Assert.That(ftp.Log.Last(), Is.EqualTo($"RETR {RawFile} @0"));
    }

    [Test]
    public void Download_ShortTransfer_IsRetriedFromZero_ThenTransient()
    {
        FakeFtp ftp = ServerWithRawFile();
        ftp.OpenOverride = (_, _) => new MemoryStream(Source, 0, 900, writable: false);
        using var client = ClientOver(ftp);
        string path = Path.Combine(_tempDir, "raw", "wGH0001.raw");

        var e = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(e), Is.True);
            Assert.That(e!.Message, Does.Contain("came to 900 bytes where the server's SIZE gives 1000"));
            Assert.That(ftp.Log.Where(l => l.StartsWith("RETR")), Is.All.EndsWith("@0"), "short bytes are never resumed");
            Assert.That(ftp.Log.Count(l => l.StartsWith("RETR")), Is.EqualTo(4));
            Assert.That(File.Exists(path + ".partial"), Is.False);
        });
    }

    [Test]
    public async Task Download_SizeDisagrees_ListingAgainAgrees_DownloadsTheServersFile()
    {
        FakeFtp ftp = ServerWithRawFile();
        using var client = ClientOver(ftp);

        string path = await client.DownloadFileAsync(RawEntry(size: 999), _tempDir);

        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(Source));
            Assert.That(ftp.Log, Does.Contain("LIST " + Root + "/raw"), "the folder was listed again once");
        });
    }

    [Test]
    public void Download_SizeDisagrees_AndListingAgainDisagrees_IsStale_NoTransfer()
    {
        FakeFtp ftp = ServerWithRawFile();
        ftp.Folders[Root + "/raw"] = new[] { "-rwxr-xr-x 1 ftp ftp          999 Sep 25  2022 wGH0001.raw" };
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.DownloadFileAsync(RawEntry(size: 999), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(e!.Message, Does.Contain("is stale"));
            Assert.That(ftp.Log.Any(l => l.StartsWith("RETR")), Is.False);
            Assert.That(_backoffs, Is.Empty);
        });
    }

    [Test]
    public void Download_FileTheServerNoLongerHas_IsAnMzLibException_NoTransfer()
    {
        var ftp = new FakeFtp();
        ftp.Folders[Root + "/raw"] = new[] { "-rwxr-xr-x 1 ftp ftp         1000 Sep 25  2022 other.raw" };
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(e!.Message, Does.Contain("no longer has"));
            Assert.That(ftp.Log.Count(l => l.StartsWith("LIST")), Is.EqualTo(1), "listed again once");
            Assert.That(ftp.Log.Any(l => l.StartsWith("RETR")), Is.False);
        });
    }

    [Test]
    public void Download_TransientReplyThatPersists_IsTransient()
    {
        FakeFtp ftp = ServerWithRawFile();
        ftp.Fail("RETR", 100, () => new MassiveFtpReplyException(425, "Can't open data connection"));
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(e), Is.True);
            Assert.That(ftp.Log.Count(l => l.StartsWith("RETR")), Is.EqualTo(4));
        });
    }

    [Test]
    public void Download_TransientReplyToSize_IsTransient()
    {
        FakeFtp ftp = ServerWithRawFile();
        ftp.Answer("SIZE", 100, new MassiveFtpReply(450, "Requested file action not taken"));
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir));

        Assert.That(ServiceFailure.IsTransient(e), Is.True);
        Assert.That(e!.Message, Does.Contain("answered 450"));
    }

    [Test]
    public void Download_PermanentReplyToRetr_IsAnMzLibException_NotRetried()
    {
        FakeFtp ftp = ServerWithRawFile();
        ftp.Fail("RETR", 100, () => new MassiveFtpReplyException(550, "Could not access file: permission denied"));
        using var client = ClientOver(ftp);

        Assert.ThrowsAsync<MzLibException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir));
        Assert.That(ftp.Log.Count(l => l.StartsWith("RETR")), Is.EqualTo(1));
    }

    [Test]
    public void Download_Stall_IsTransient_AndKeepsThePartial()
    {
        FakeFtp ftp = ServerWithRawFile();
        ftp.OpenOverride = (_, _) => new StallingStream();
        using var client = ClientOver(ftp, maxRetries: 0, stall: TimeSpan.FromMilliseconds(200));
        string path = Path.Combine(_tempDir, "raw", "wGH0001.raw");

        var e = Assert.ThrowsAsync<HttpRequestException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir));

        Assert.Multiple(() =>
        {
            Assert.That(ServiceFailure.IsTransient(e), Is.True);
            Assert.That(e!.Message, Does.Contain("delivered nothing"));
            Assert.That(File.Exists(path + ".partial.validator"), Is.True, "a transient failure keeps what the next call can resume");
        });
    }

    [Test]
    public void Download_Cancelled_KeepsThePartial_AndIsACancellation()
    {
        FakeFtp ftp = ServerWithRawFile();
        ftp.OpenOverride = (_, _) => new StallingStream();
        using var client = ClientOver(ftp);
        using var cancel = new CancellationTokenSource(TimeSpan.FromMilliseconds(200));
        string path = Path.Combine(_tempDir, "raw", "wGH0001.raw");

        Assert.CatchAsync<OperationCanceledException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir, cancellationToken: cancel.Token));

        Assert.That(File.Exists(path + ".partial"), Is.True);
        Assert.That(File.Exists(path + ".partial.validator"), Is.True);
    }

    [Test]
    public void Download_EveryFailureMessage_NamesNoUrlAndNoLogin()
    {
        var failures = new List<Exception>();
        foreach (Action<FakeFtp> script in new Action<FakeFtp>[]
                 {
                     f => f.Fail("CONNECT", 100, () => new MassiveFtpTransportException("The FTP connection failed: refused", new SocketException(10061))),
                     f => f.Fail("CONNECT", 100, () => new AuthenticationException("certificate chain")),
                     f => f.Fail("RETR", 100, () => new MassiveFtpReplyException(425, "Can't open data connection")),
                     f => f.Fail("RETR", 100, () => new MassiveFtpReplyException(550, "Could not access file")),
                     f => f.OpenOverride = (_, _) => new DroppingStream(Source, 0, 10),
                     f => f.Files.Remove(RawFile),
                 })
        {
            FakeFtp ftp = ServerWithRawFile();
            script(ftp);
            using var client = ClientOver(ftp, maxRetries: 0);
            failures.Add(Assert.CatchAsync<Exception>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir)));
        }

        Assert.Multiple(() =>
        {
            foreach (Exception failure in failures)
            {
                Assert.That(failure.ToString(), Does.Not.Contain("ftp://").And.Not.Contain("anonymous"));
                Assert.That(failure.Message, Does.Contain("raw/wGH0001.raw"), "the file is named by its path below the dataset");
            }
        });
    }

    [Test]
    public async Task DownloadDatasetFiles_SkipsLinks_AppliesTheFilter_InListingOrder()
    {
        FakeFtp ftp = ServerWithTree();
        ftp.Files[RawFile] = Source;
        ftp.Files[Root + "/raw/License File/license.txt"] = new byte[456];
        ftp.Files[Root + "/updates/2022-09-27_RGHuffman_25d0cd1a/raw/wGH0697.raw"] = Source;
        using var client = ClientOver(ftp);

        IReadOnlyList<string> paths = await client.DownloadDatasetFilesAsync("MSV000090383", _tempDir, f => f.Collection == "raw" || f.IsLink);
        IReadOnlyList<string> none = await client.DownloadDatasetFilesAsync("MSV000090383", _tempDir, f => f.FileName == "nothing.raw");

        Assert.Multiple(() =>
        {
            Assert.That(paths.Select(p => Path.GetRelativePath(_tempDir, p).Replace('\\', '/')), Is.EqualTo(new[]
            {
                "raw/wGH0001.raw",
                "raw/License File/license.txt",
                "updates/2022-09-27_RGHuffman_25d0cd1a/raw/wGH0697.raw",
            }), "the filter chose the link too, and a link is never downloaded");
            Assert.That(none, Is.Empty);
            Assert.That(ftp.Log.Any(l => l.Contains("ccms_peak") && l.StartsWith("RETR")), Is.False);
        });
    }

    // ---- the remaining branches --------------------------------------------------------------------------------

    [Test]
    public void GetDataset_EmptyBody_IsAContractBreak()
    {
        using var client = ClientOver(new FakeFtp(), new ProxiHandler(_ => Proxi(Array.Empty<byte>())));

        Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetAsync("MSV000096226"));
    }

    [Test]
    public void GetDataset_ProteomeCentralNamesMassive_SaysTheTwoDisagree()
    {
        var handler = new ProxiHandler(r => r.RequestUri!.Host == "massive.ucsd.edu"
            ? Proxi(NotFoundBody, HttpStatusCode.NotFound)
            : new HttpResponseMessage(HttpStatusCode.OK) { Content = new StringContent("""{"datasetSummary":{"hostingRepository":"MassIVE"}}""") });
        using var client = ClientOver(new FakeFtp(), handler);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetAsync("PXD057243"));

        Assert.That(e!.Message, Does.Contain("ProteomeCentral says MassIVE hosts 'PXD057243'"));
    }

    [Test]
    public void GetDataset_ProteomeCentralUnreachable_SaysItCouldNotBeAsked()
    {
        var handler = new ProxiHandler(r => r.RequestUri!.Host == "massive.ucsd.edu"
            ? Proxi(NotFoundBody, HttpStatusCode.NotFound)
            : throw new HttpRequestException("No such host is known."));
        using var client = ClientOver(new FakeFtp(), handler);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetAsync("PXD000001"));

        Assert.That(e!.Message, Does.Contain("ProteomeCentral could not be asked"));
    }

    [Test]
    public void GetDatasetFiles_CyclicListing_StopsAtTheDepthCap()
    {
        var ftp = new FakeFtp { ListFallback = _ => new[] { "drwxr-xr-x 1 ftp ftp            1 Sep 25  2022 a" } };
        ftp.Folders[Root] = new[] { "drwxr-xr-x 1 ftp ftp            1 Sep 25  2022 a" };
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetFilesAsync("MSV000090383"));

        Assert.That(e!.Message, Does.Contain("nests deeper than 64 folders"));
    }

    [Test]
    public async Task GetDatasetFiles_ConnectionThatFailsToClose_DoesNotHideTheRetry()
    {
        FakeFtp ftp = ServerWithTree();
        ftp.ThrowOnDispose = true;
        ftp.Fail("LIST", 1, () => new MassiveFtpTransportException("The FTP connection failed: reset", new SocketException(10054)));
        using var client = ClientOver(ftp);

        IReadOnlyList<MassiveFile> files = await client.GetDatasetFilesAsync("MSV000090383");

        Assert.That(files, Has.Count.EqualTo(5));
        Assert.That(ftp.Disposals, Is.EqualTo(1), "the failed connection was closed, and its failure to close was swallowed");
    }

    [TestCase("SIZE", "213", "not a size")]
    [TestCase("MDTM", "213", "2015 0212140541")]
    public void Download_ReplyThatIsNotAValue_IsAContractBreak(string verb, string code, string message)
    {
        FakeFtp ftp = ServerWithRawFile();
        ftp.Answer(verb, 1, new MassiveFtpReply(int.Parse(code), message));
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir));

        Assert.That(e!.Message, Does.Contain(verb).And.Contain(message));
    }

    [Test]
    public void Download_ModifiedTimeMissing_IsAnMzLibException()
    {
        FakeFtp ftp = ServerWithRawFile();
        ftp.Answer("MDTM", 1, new MassiveFtpReply(550, "file does not exist"));
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir));

        Assert.That(e!.Message, Does.Contain("no longer has"));
    }

    [Test]
    public void Download_FileChangingAgainAfterTheRelist_IsAnMzLibException()
    {
        // SIZE first says 999 and the folder, listed again, agrees; the transfer then brings the real 1000 bytes, and
        // the retry's SIZE says 1000: the file has changed twice within one call.
        FakeFtp ftp = ServerWithRawFile();
        ftp.Answer("SIZE", 1, new MassiveFtpReply(213, "999"));
        ftp.Folders[Root + "/raw"] = new[] { "-rwxr-xr-x 1 ftp ftp          999 Sep 25  2022 wGH0001.raw" };
        using var client = ClientOver(ftp);

        var e = Assert.ThrowsAsync<MzLibException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir));

        Assert.That(e!.Message, Does.Contain("changed while it was being downloaded"));
    }

    [Test]
    public void Download_ArgumentsAreCheckedBeforeAnyRequest()
    {
        var ftp = new FakeFtp();
        using var client = ClientOver(ftp);

        Assert.Multiple(() =>
        {
            Assert.ThrowsAsync<ArgumentNullException>(async () => await client.DownloadFileAsync(null, _tempDir));
            Assert.ThrowsAsync<ArgumentException>(async () => await client.DownloadFileAsync(RawEntry(), " "));
            Assert.ThrowsAsync<ArgumentException>(async () => await client.DownloadDatasetFilesAsync("MSV000090383", ""));
            Assert.That(ftp.Connections, Is.Zero);
        });
    }

    [Test]
    public void Client_Constructors_AndDispose()
    {
        using (var owning = new MassiveClient())
            Assert.That(owning.BodyStallTimeout, Is.EqualTo(TimeSpan.FromMinutes(2)));

        Assert.Throws<ArgumentNullException>(() => _ = new MassiveClient((HttpClient)null));

        var client = ClientOver(ServerWithRawFile());
        client.Dispose();
        client.Dispose();
        Assert.ThrowsAsync<ObjectDisposedException>(async () => await client.DownloadFileAsync(RawEntry(), _tempDir));
    }

    [TestCase("v07/MSV000096226", "v07/MSV000096226")]
    [TestCase("ftp://massive.ucsd.edu", "")]
    [TestCase("ftp://massive.ucsd.edu/v07/MSV%30%30%30096226/", "v07/MSV000096226")]
    public void Record_FtpLocation_IsTheRecordsPath(string location, string expected)
    {
        string json = $$"""{"accession":[{"accession":"MS:1002487","value":"MSV000096226"}],"datasetLink":[{"accession":"MS:1002852","value":"{{location}}"}]}""";

        Assert.That(MassiveDataset.FromProxi(Encoding.UTF8.GetBytes(json), "MSV000096226").FtpLocation, Is.EqualTo(expected));
    }

    // ---- MassiveFile -------------------------------------------------------------------------------------------

    [TestCase("", "/v07/MSV000096226", "raw/x.raw", 1L)]
    [TestCase("MSV000096226", "v07/MSV000096226", "raw/x.raw", 1L)]
    [TestCase("MSV000096226", "/v07/MSV000096226/", "raw/x.raw", 1L)]
    [TestCase("MSV000096226", "/v07/MSV000096226", "/raw/x.raw", 1L)]
    [TestCase("MSV000096226", "/v07/MSV000096226", "raw/", 1L)]
    [TestCase("MSV000096226", "/v07/MSV000096226", "raw/x.raw", null)]
    [TestCase("MSV000096226", "/v07/MSV000096226", "raw/x.raw", -1L)]
    public void MassiveFile_BadArguments_AreRefused(string accession, string root, string path, long? size)
    {
        Assert.Throws<ArgumentException>(() => _ = new MassiveFile(accession, root, path, size));
    }

    [Test]
    public void MassiveFile_ToString_NamesDatasetAndPath()
    {
        Assert.That(new MassiveFile("MSV000096226", "/v07/MSV000096226", "raw/x.raw", 1).ToString(), Is.EqualTo("MSV000096226:raw/x.raw"));
    }

    [TestCase("raw/x.raw", false, "raw", null)]
    [TestCase("README.txt", false, null, null)]
    [TestCase("ccms_peak", true, "ccms_peak", null)]
    [TestCase("updates/2023-01-30_jburton1_42cc6275/ccms_parameters/params.xml", false, "ccms_parameters", "2023-01-30_jburton1_42cc6275")]
    [TestCase("updates/2023-01-30_jburton1_42cc6275/notes.txt", false, null, "2023-01-30_jburton1_42cc6275")]
    public void MassiveFile_CollectionAndUpdate_ComeFromThePath(string path, bool isLink, string collection, string update)
    {
        var file = new MassiveFile("MSV000090777", "/v05/MSV000090777", path, isLink ? null : 1, isLink: isLink);

        Assert.That(file.Collection, Is.EqualTo(collection));
        Assert.That(file.UpdateName, Is.EqualTo(update));
    }

    // ---- the FluentFTP connection ------------------------------------------------------------------------------

    [Test]
    public void FluentFtp_Config_ValidatesCertificates_OverExplicitTls_InBinary()
    {
        var config = FluentFtpConnection.CreateConfig(TimeSpan.FromSeconds(100), TimeSpan.FromMinutes(2.5));

        Assert.Multiple(() =>
        {
            Assert.That(config.ValidateAnyCertificate, Is.False);
            Assert.That(config.EncryptionMode, Is.EqualTo(FluentFTP.FtpEncryptionMode.Explicit));
            Assert.That(config.DataConnectionEncryption, Is.True);
            Assert.That(config.DownloadDataType, Is.EqualTo(FluentFTP.FtpDataType.Binary));
            Assert.That(config.RetryAttempts, Is.EqualTo(1), "the client's retry loop is the only one");
            Assert.That(config.DataConnectionReadTimeout, Is.EqualTo(150000));
        });
    }

    [Test]
    public void FluentFtp_Client_HasACertificateHandler_ThatAcceptsOnlyACleanCertificate()
    {
        using var client = FluentFtpConnection.CreateClient(MassiveClient.FtpHost, MassiveClient.FtpPort, TimeSpan.FromSeconds(100), TimeSpan.FromMinutes(2));

        Assert.Multiple(() =>
        {
            Assert.That(client.ValidateCertificateHandlerExists, Is.True);
            Assert.That(FluentFtpConnection.AcceptCertificate(SslPolicyErrors.None), Is.True);
            Assert.That(FluentFtpConnection.AcceptCertificate(SslPolicyErrors.RemoteCertificateNameMismatch), Is.False);
            Assert.That(FluentFtpConnection.AcceptCertificate(SslPolicyErrors.RemoteCertificateChainErrors), Is.False);
            Assert.That(FluentFtpConnection.AcceptCertificate(SslPolicyErrors.RemoteCertificateNotAvailable), Is.False);
        });
    }

    [Test]
    public void FluentFtp_Translate_MapsEachFailureOntoTheSeam()
    {
        Assert.Multiple(() =>
        {
            Assert.That(FluentFtpConnection.Translate(new FtpCommandException("450", "Could not list")),
                Is.InstanceOf<MassiveFtpReplyException>().With.Property(nameof(MassiveFtpReplyException.Code)).EqualTo(450));
            Assert.That(FluentFtpConnection.Translate(new FtpCommandException("550", "Could not access file")),
                Is.InstanceOf<MassiveFtpReplyException>().With.Property(nameof(MassiveFtpReplyException.IsTransient)).False);
            Assert.That(FluentFtpConnection.Translate(new FtpInvalidCertificateException("RemoteCertificateNameMismatch")), Is.InstanceOf<AuthenticationException>());
            Assert.That(FluentFtpConnection.Translate(new FtpSecurityNotAvailableException()), Is.InstanceOf<AuthenticationException>());
            Assert.That(FluentFtpConnection.Translate(new AuthenticationException("x")), Is.InstanceOf<AuthenticationException>());
            Assert.That(FluentFtpConnection.Translate(new TimeoutException()), Is.InstanceOf<MassiveFtpTransportException>());
            Assert.That(FluentFtpConnection.Translate(new SocketException(10054)), Is.InstanceOf<MassiveFtpTransportException>());
            Assert.That(FluentFtpConnection.Translate(new FtpException("x")), Is.InstanceOf<MassiveFtpTransportException>());
            Assert.That(FluentFtpConnection.Translate(new IOException("x")), Is.InstanceOf<MassiveFtpTransportException>());
            Assert.That(FluentFtpConnection.Translate(new OperationCanceledException()), Is.Null, "a cancellation stays a cancellation");
            Assert.That(FluentFtpConnection.Translate(new InvalidOperationException()), Is.Null, "a bug is not an outage");
        });
    }
}

/// <summary>
/// The FluentFTP side of the seam, offline: a loopback server plays the control connection up to and through
/// <c>AUTH TLS</c>, so the connection's failure kinds, and the certificate check against a real TLS handshake, are
/// pinned without the network.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class MassiveFtpConnectionTests
{
    private static readonly TimeSpan Timeout = TimeSpan.FromSeconds(10);

    /// <summary>Accepts one control connection on a loopback port and runs <c>script</c> over it.</summary>
    private sealed class LoopbackServer : IDisposable
    {
        private readonly TcpListener _listener = new(IPAddress.Loopback, 0);

        public LoopbackServer(Func<NetworkStream, Task> script)
        {
            _listener.Start();
            Served = Task.Run(async () =>
            {
                using TcpClient client = await _listener.AcceptTcpClientAsync();
                using NetworkStream stream = client.GetStream();
                try { await script(stream); }
                catch (Exception) { /* the client hanging up mid-script is what several tests expect */ }
            });
        }

        public int Port => ((IPEndPoint)_listener.LocalEndpoint).Port;
        public Task Served { get; }
        public void Dispose() => _listener.Stop();
    }

    private static async Task WriteLineAsync(Stream stream, string line)
    {
        byte[] bytes = Encoding.ASCII.GetBytes(line + "\r\n");
        await stream.WriteAsync(bytes);
        await stream.FlushAsync();
    }

    /// <summary>Reads one line byte by byte, so nothing after it (a TLS handshake) is buffered away.</summary>
    private static async Task<string> ReadLineAsync(Stream stream)
    {
        var line = new StringBuilder();
        byte[] one = new byte[1];
        while (await stream.ReadAsync(one) == 1 && one[0] != '\n')
        {
            if (one[0] != '\r')
                line.Append((char)one[0]);
        }
        return line.ToString();
    }

    private static Task<IMassiveFtpConnection> Connect(LoopbackServer server) =>
        FluentFtpConnection.ConnectAsync("127.0.0.1", server.Port, Timeout, Timeout, CancellationToken.None);

    [Test]
    public void Connect_ServerWithoutTls_IsAnAuthenticationFailure()
    {
        using var server = new LoopbackServer(async s =>
        {
            await WriteLineAsync(s, "220 plain FTP only");
            while ((await ReadLineAsync(s)).Length > 0)
                await WriteLineAsync(s, "500 Command not understood");
        });

        Exception e = Assert.CatchAsync(async () => await Connect(server));

        Assert.That(e, Is.InstanceOf<AuthenticationException>(), e!.ToString());
    }

    [Test]
    public void Connect_CertificateThatFailsValidation_IsRefused()
    {
        // Self-signed, for "localhost", while the client asked for 127.0.0.1: two reasons for the check to refuse it.
        using RSA rsa = RSA.Create(2048);
        var request = new CertificateRequest("CN=localhost", rsa, HashAlgorithmName.SHA256, RSASignaturePadding.Pkcs1);
        using X509Certificate2 ephemeral = request.CreateSelfSigned(DateTimeOffset.UtcNow.AddDays(-1), DateTimeOffset.UtcNow.AddDays(1));
        using X509Certificate2 certificate = X509CertificateLoader.LoadPkcs12(ephemeral.Export(X509ContentType.Pfx), null);
        bool handshakeStarted = false;
        using var server = new LoopbackServer(async s =>
        {
            await WriteLineAsync(s, "220 SFTPGo_2.7.1");
            await ReadLineAsync(s); // AUTH TLS
            await WriteLineAsync(s, "234 AUTH command ok. Expecting TLS Negotiation.");
            handshakeStarted = true;
            using var tls = new SslStream(s);
            await tls.AuthenticateAsServerAsync(certificate);
        });

        Exception e = Assert.CatchAsync(async () => await Connect(server));

        Assert.Multiple(() =>
        {
            Assert.That(handshakeStarted, Is.True, "the refusal must come from the TLS handshake, not before it");
            Assert.That(e, Is.InstanceOf<AuthenticationException>(), e!.ToString());
        });
    }

    [Test]
    public void Connect_ServerThatHangsUpAtOnce_IsATransportFailure()
    {
        using var server = new LoopbackServer(_ => Task.CompletedTask);

        Exception e = Assert.CatchAsync(async () => await Connect(server));

        Assert.That(e, Is.InstanceOf<MassiveFtpTransportException>(), e!.ToString());
    }

    /// <summary>A stream whose reads and whose closing fail as told.</summary>
    private sealed class FailingStream(Exception onRead, Exception onDispose) : MemoryStream(new byte[] { 1, 2, 3 })
    {
        public override int Read(byte[] buffer, int offset, int count) => onRead != null ? throw onRead : base.Read(buffer, offset, count);

        public override ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default) =>
            onRead != null ? throw onRead : base.ReadAsync(buffer, cancellationToken);

        protected override void Dispose(bool disposing)
        {
            base.Dispose(disposing);
            if (onDispose != null)
                throw onDispose;
        }

        public override ValueTask DisposeAsync() => onDispose != null ? throw onDispose : base.DisposeAsync();
    }

    [Test]
    public async Task ReadStream_PassesBytesThrough_EveryWayItCanBeRead()
    {
        var stream = new MassiveFtpReadStream(new MemoryStream(new byte[] { 1, 2, 3, 4, 5, 6 }));
        byte[] buffer = new byte[2];

        int a = stream.Read(buffer, 0, 2);
        int b = await stream.ReadAsync(buffer, 0, 2, CancellationToken.None);
        int c = await stream.ReadAsync(buffer.AsMemory());
        await stream.DisposeAsync();

        Assert.Multiple(() =>
        {
            Assert.That(new[] { a, b, c }, Is.EqualTo(new[] { 2, 2, 2 }));
            Assert.That(buffer, Is.EqualTo(new byte[] { 5, 6 }));
            Assert.That(stream.CanRead && !stream.CanSeek && !stream.CanWrite, Is.True);
            Assert.Throws<NotSupportedException>(() => _ = stream.Length);
            Assert.Throws<NotSupportedException>(() => _ = stream.Position);
            Assert.Throws<NotSupportedException>(() => stream.Position = 1);
            Assert.Throws<NotSupportedException>(() => stream.Seek(0, SeekOrigin.Begin));
            Assert.Throws<NotSupportedException>(() => stream.SetLength(1));
            Assert.Throws<NotSupportedException>(() => stream.Write(buffer, 0, 1));
            Assert.DoesNotThrow(() => stream.Flush());
        });
    }

    [Test]
    public void ReadStream_ReadFailure_IsATransportFailure_ACancellationStaysOne()
    {
        var timedOut = new MassiveFtpReadStream(new FailingStream(new TimeoutException("Timed out trying to read data from the socket stream!"), null));
        var reset = new MassiveFtpReadStream(new FailingStream(new SocketException(10054), null));
        var replied = new MassiveFtpReadStream(new FailingStream(new FtpCommandException("426", "Connection closed; transfer aborted"), null));
        var cancelled = new MassiveFtpReadStream(new FailingStream(new OperationCanceledException(), null));

        Assert.Multiple(() =>
        {
            Assert.ThrowsAsync<MassiveFtpTransportException>(async () => await timedOut.ReadAsync(new byte[4].AsMemory()));
            Assert.Throws<MassiveFtpTransportException>(() => reset.Read(new byte[4], 0, 4));
            Assert.ThrowsAsync<MassiveFtpTransportException>(async () => await replied.ReadAsync(new byte[4].AsMemory()));
            Assert.ThrowsAsync<OperationCanceledException>(async () => await cancelled.ReadAsync(new byte[4].AsMemory()));
        });
    }

    [Test]
    public void ReadStream_FailureToClose_IsIgnored()
    {
        Assert.Multiple(() =>
        {
            Assert.DoesNotThrow(() => new MassiveFtpReadStream(new FailingStream(null, new IOException("226 never came"))).Dispose());
            Assert.DoesNotThrowAsync(async () => await new MassiveFtpReadStream(new FailingStream(null, new TimeoutException())).DisposeAsync());
        });
    }
}

/// <summary>
/// Live canaries against MassIVE itself. [Category("ExternalService")] keeps them out of the required run;
/// <see cref="ExternalServiceTestHelper.RunAsync"/> SKIPS them on an outage and lets a format drift FAIL. One attempt
/// each (<c>MaxRetries = 0</c>), few calls, and at least 6 s between PROXI calls.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("Massive")]
[ExcludeFromCodeCoverage]
public class MassiveClientLiveTests
{
    private string _tempDir;

    [SetUp]
    public async Task SetUp()
    {
        _tempDir = Path.Combine(Path.GetTempPath(), "MassiveClientLiveTests", Guid.NewGuid().ToString("N"));
        await Task.Delay(TimeSpan.FromSeconds(6)); // MassIVE rate-limits; space the PROXI calls
    }

    [TearDown]
    public void TearDown()
    {
        try { if (Directory.Exists(_tempDir)) Directory.Delete(_tempDir, recursive: true); }
        catch { /* best-effort cleanup */ }
    }

    [Test]
    public Task GetDataset_LivePxd_ResolvesToItsMsv() =>
        ExternalServiceTestHelper.RunAsync("MassIVE", async () =>
        {
            using var client = new MassiveClient { MaxRetries = 0 };

            MassiveDataset d = await client.GetDatasetAsync("PXD057243");

            Assert.That(d.MassiveAccession, Is.EqualTo("MSV000096226"));
            Assert.That(d.PxdAccession, Is.EqualTo("PXD057243"));
            Assert.That(d.FtpLocation, Is.EqualTo("v07/MSV000096226"));
            Assert.That(d.Title, Is.Not.Empty);
        });

    [Test]
    public Task TryGetDataset_LiveUnknownMsv_IsNotFound() =>
        ExternalServiceTestHelper.RunAsync("MassIVE", async () =>
        {
            using var client = new MassiveClient { MaxRetries = 0 };

            (bool found, MassiveDataset dataset) = await client.TryGetDatasetAsync("MSV999999999");

            Assert.That(found, Is.False);
            Assert.That(dataset, Is.Null);
        });

    [Test]
    public Task GetDataset_LivePrideHostedPxd_NamesPride() =>
        ExternalServiceTestHelper.RunAsync("MassIVE", async () =>
        {
            using var client = new MassiveClient { MaxRetries = 0 };

            var e = Assert.ThrowsAsync<MzLibException>(async () => await client.GetDatasetAsync("PXD000001"));

            if (e!.Message.Contains("could not be asked"))
                Assert.Ignore("ProteomeCentral was unavailable, so the host could not be named: " + e.Message);
            Assert.That(e.Message, Does.Contain("hosted at PRIDE"));
        });

    [Test]
    public Task GetDatasetFiles_LiveDatasetWithALink_ListsTheLinkAndExactSizes() =>
        ExternalServiceTestHelper.RunAsync("MassIVE", async () =>
        {
            using var client = new MassiveClient { MaxRetries = 0 };

            IReadOnlyList<MassiveFile> files = await client.GetDatasetFilesAsync("MSV000080763");

            MassiveFile link = files.Single(f => f.Path == "ccms_peak");
            MassiveFile parameters = files.Single(f => f.Path == "ccms_parameters/params.xml");
            Assert.That(link.IsLink, Is.True, "ccms_peak is a link MassIVE's FTP server will not open (2026-11-10)");
            Assert.That(parameters.SizeBytes, Is.EqualTo(7339));
            Assert.That(files.Single(f => f.FileName == "Saliva_txt.rar").SizeBytes, Is.EqualTo(1157656623));
        });

    [Test]
    public Task DownloadFile_LiveSmallFile_ArrivesWhole() =>
        ExternalServiceTestHelper.RunAsync("MassIVE", async () =>
        {
            using var client = new MassiveClient { MaxRetries = 0 };
            var file = new MassiveFile("MSV000080763", "/x01/MSV000080763", "ccms_parameters/params.xml", 7339);

            string path = await client.DownloadFileAsync(file, _tempDir);

            byte[] bytes = File.ReadAllBytes(path);
            Assert.That(bytes, Has.Length.EqualTo(7339));
            Assert.That(Encoding.ASCII.GetString(bytes, 0, 5), Is.EqualTo("<?xml"));
        });

    [Test]
    public Task DownloadFile_LiveUtf16TextFile_ArrivesByteForByte_InBinaryMode() =>
        ExternalServiceTestHelper.RunAsync("MassIVE", async () =>
        {
            using var client = new MassiveClient { MaxRetries = 0 };
            MassiveDataset d = await client.GetDatasetAsync("MSV000086944");
            var file = new MassiveFile(d.MassiveAccession, "/" + d.FtpLocation, "metadata/experimental.design.txt", 311842);

            string path = await client.DownloadFileAsync(file, _tempDir);

            // FTP's ASCII mode, like MassIVE's HTTPS endpoint (sdrf 028), would rewrite the line ends of a text file.
            byte[] bytes = File.ReadAllBytes(path);
            Assert.That(bytes, Has.Length.EqualTo(311842));
            Assert.That(bytes[0], Is.EqualTo(0xFF));
            Assert.That(bytes[1], Is.EqualTo(0xFE));
        });
}
